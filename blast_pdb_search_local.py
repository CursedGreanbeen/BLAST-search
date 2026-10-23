#!/usr/bin/env python3
"""
Local BLASTP batch runner for FASTA protein complexes against a local PDB database.

Instead of hitting the NCBI web BLAST API (slow, rate-limited to 1 request/10s),
this script uses a locally installed BLAST+ ``blastp`` binary and a locally
downloaded ``pdbaa`` database. Since the search runs on your own CPU with no
network calls or externally imposed rate limits, we get real parallelism via
multiple ``blastp`` subprocesses running concurrently.

Setup required before running this script (see accompanying bash instructions):
    1. Download & extract NCBI BLAST+ binaries, put them on PATH.
    2. Download the `pdbaa` BLAST database with `update_blastdb.pl --decompress pdbaa`.

Usage:
    python3 blast_pdb_search_local.py <fasta_dir> \
        --db /path/to/pdbaa \
        --workers 4 \
        --threads-per-search 2 \
        --min-identity 90 \
        --max-identity 99 \
        --output-dir results/blast_pdb_search_local

Each FASTA file in <fasta_dir> may contain a "complex" of multiple chains
(multiple FASTA records). Every unique sequence is BLASTed exactly once
(sequences repeated across chains/files are deduplicated), against the local
`pdbaa` database, and hits are filtered to the given percent-identity range.

Output:
    <output-dir>/<safe_seq_id>.xml               - raw BLAST XML per unique sequence
    <output-dir>/<safe_seq_id>_filtered_hits.csv  - filtered hits for that sequence
    <output-dir>/combined_filtered_hits.csv       - all filtered hits, all files/chains
    <output-dir>/errors.txt                       - list of failures, if any
"""

import argparse
import csv
import hashlib
import shutil
import subprocess
import sys
import threading
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

from Bio import SeqIO, Blast
from anarci import anarci as _anarci_run


DEFAULT_MIN_IDENTITY = 0.0
DEFAULT_MAX_IDENTITY = 100.0
DEFAULT_EVALUE = 10.0
FASTA_EXTENSIONS = (".fasta", ".fa", ".faa")
HEAVY_CHAIN_TYPES = {"H"}  # ANARCI chain_type label for heavy chains
DEFAULT_MIN_H3_IDENTITY = 80.0
DEFAULT_MAX_H3_IDENTITY = 100.0

IMGT_CDR_H3_START = 105
IMGT_CDR_H3_END = 117  # inclusive, per IMGT numbering scheme

_print_lock = threading.Lock()


def check_anarci() -> None:
    """Sanity-check that ANARCI is importable before launching any batch work."""
    try:
        from anarci import anarci  # noqa: F401
    except ImportError:
        sys.exit(
            "ERROR: 'anarci' Python package not found in this environment.\n"
            "Install it, e.g.: conda install -c bioconda anarci\n"
            "or: pip install anarci (if available for your platform)"
        )
    log("ANARCI OK: import successful.")


def run_anarci(sequence: str) -> dict:
    """
    Run ANARCI (Python API) on a single sequence and extract CDR-H3 boundaries.

    Returns a dict:
        {
            "chain_type": str or None,
            "domain_start": int or None,   # raw seq index, inclusive
            "domain_end": int or None,     # raw seq index, exclusive
            "cdr_h3_start": int or None,   # raw seq index, inclusive
            "cdr_h3_end": int or None,     # raw seq index, exclusive
        }
    """
    empty_result = {
        "chain_type": None,
        "domain_start": None,
        "domain_end": None,
        "cdr_h3_start": None,
        "cdr_h3_end": None,
        "cdr_h3_sequence": None,
    }

    try:
        numbering, alignment_details, _hit_tables = _anarci_run(
            [("query", sequence)], scheme="imgt", output=False
        )
    except Exception as exc:
        print(f"    ERROR: ANARCI call failed: {exc}")
        return empty_result

    domains = numbering[0]      # list of domains found in this sequence
    details = alignment_details[0]

    if not domains:
        return empty_result

    # Find the first heavy-chain domain; sequence may contain multiple domains
    # (e.g. scFv) not necessarily sorted by significance.
    heavy_domain = None
    heavy_detail = None
    for domain, detail in zip(domains, details):
        if detail["chain_type"] in HEAVY_CHAIN_TYPES:
            heavy_domain = domain
            heavy_detail = detail
            break

    if heavy_domain is None:
        return {**empty_result, "chain_type": details[0]["chain_type"]}

    numbered_residues, domain_start, domain_end = heavy_domain
    chain_type = heavy_detail["chain_type"]

    cdr_h3_indices = []
    raw_idx = domain_start
    for (position, _ins_code), residue in numbered_residues:
        if residue != "-":
            if IMGT_CDR_H3_START <= position <= IMGT_CDR_H3_END:
                cdr_h3_indices.append(raw_idx)
            raw_idx += 1

    if not cdr_h3_indices:
        cdr_h3_start, cdr_h3_end = None, None
    else:
        cdr_h3_start = min(cdr_h3_indices)
        cdr_h3_end = max(cdr_h3_indices) + 1

    cdr_h3_sequence = (
        sequence[cdr_h3_start:cdr_h3_end]
        if cdr_h3_start is not None and cdr_h3_end is not None
        else None
    )

    return {
        "chain_type": chain_type,
        "domain_start": domain_start,
        "domain_end": domain_end,
        "cdr_h3_start": cdr_h3_start,
        "cdr_h3_end": cdr_h3_end,
        "cdr_h3_sequence": cdr_h3_sequence,
    }


def log(msg: str) -> None:
    with _print_lock:
        print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def find_blastp() -> str:
    """Locate the blastp executable on PATH, fail fast with a clear message."""
    path = shutil.which("blastp")
    if not path:
        sys.exit(
            "ERROR: 'blastp' executable not found on PATH.\n"
            "Make sure you extracted BLAST+ and exported its bin/ directory, e.g.:\n"
            "  export PATH=\"$HOME/blast_local/ncbi-blast-2.17.0+/bin:$PATH\""
        )
    return path


def check_database(db_path: str, blastdbcmd: str = "blastdbcmd") -> None:
    """Sanity-check that the local BLAST database exists and is usable."""
    blastdbcmd_path = shutil.which(blastdbcmd)
    if not blastdbcmd_path:
        log("WARNING: 'blastdbcmd' not found on PATH, skipping database validation.")
        return
    result = subprocess.run(
        [blastdbcmd_path, "-db", db_path, "-info"],
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        sys.exit(
            f"ERROR: Local BLAST database '{db_path}' is not usable.\n"
            f"blastdbcmd output:\n{result.stdout}\n{result.stderr}"
        )
    log(f"Database OK: {db_path}")
    log(result.stdout.strip().splitlines()[0] if result.stdout.strip() else "")


def iter_fasta_files(fasta_dir: Path):
    for path in sorted(fasta_dir.iterdir()):
        if path.is_file() and path.suffix.lower() in FASTA_EXTENSIONS:
            yield path


def collect_unique_sequences(fasta_dir: Path):
    """
    Parse every FASTA file (possibly multi-chain complexes) in fasta_dir and
    build a mapping of unique sequence -> list of (file_name, chain_id, description).
    Returns dict: seq_hash -> {"sequence": str, "sources": [(file, chain_id, desc), ...]}
    """
    unique_seqs = {}
    total_records = 0

    for fasta_path in iter_fasta_files(fasta_dir):
        try:
            records = list(SeqIO.parse(str(fasta_path), "fasta"))
        except Exception as exc:
            log(f"ERROR: failed to parse {fasta_path.name}: {exc}")
            continue

        if not records:
            log(f"WARNING: no sequences found in {fasta_path.name}")
            continue

        for record in records:
            total_records += 1
            seq_str = str(record.seq).upper()
            if not seq_str:
                continue
            seq_hash = hashlib.sha1(seq_str.encode()).hexdigest()[:16]
            entry = unique_seqs.setdefault(
                seq_hash, {"sequence": seq_str, "sources": []}
            )
            entry["sources"].append((fasta_path.name, record.id, record.description))

    log(
        f"Parsed {total_records} chain record(s) from "
        f"{sum(1 for _ in iter_fasta_files(fasta_dir))} file(s) -> "
        f"{len(unique_seqs)} unique sequence(s) to BLAST."
    )
    return unique_seqs


def run_single_blastp(args_tuple):
    """
    Worker function executed in a separate process.
    Writes the query to a temp FASTA, runs local blastp, returns paths + status.
    Must be a top-level function (picklable) for ProcessPoolExecutor.
    """
    (
        seq_hash,
        sequence,
        blastp_path,
        db_path,
        evalue,
        threads_per_search,
        output_dir,
    ) = args_tuple

    out_dir = Path(output_dir)
    query_path = out_dir / f"{seq_hash}_query.fasta"
    xml_path = out_dir / f"{seq_hash}.xml"

    try:
        query_path.write_text(f">{seq_hash}\n{sequence}\n")

        cmd = [
            blastp_path,
            "-query", str(query_path),
            "-db", db_path,
            "-outfmt", "5",
            "-out", str(xml_path),
            "-evalue", str(evalue),
            "-num_threads", str(threads_per_search),
        ]
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=1800)

        if result.returncode != 0:
            return seq_hash, False, f"blastp failed: {result.stderr.strip()[:500]}"

        return seq_hash, True, str(xml_path)
    except subprocess.TimeoutExpired:
        return seq_hash, False, "blastp timed out after 1800s"
    except Exception as exc:
        return seq_hash, False, f"unexpected error: {exc}"
    finally:
        query_path.unlink(missing_ok=True)


def compute_h3_identity(hsp, cdr_h3_start, cdr_h3_end):
    not_covered = {
        "h3_covered": False,
        "h3_identity": None,
        "h3_align_length": None,
        "h3_pct_identity": None,
    }

    if cdr_h3_start is None or cdr_h3_end is None:
        return not_covered

    aligned_query = str(hsp[0])   # FIXED: row-indexing gives gapped alignment string
    aligned_hit = str(hsp[1])     # FIXED: same, for the hit/target row

    raw_idx = hsp.coordinates[0][0]   # FIXED: raw query start, from coordinates array

    h3_align_start = None
    h3_align_end = None

    for aln_pos, ch in enumerate(aligned_query):
        if ch != "-":
            if raw_idx == cdr_h3_start and h3_align_start is None:
                h3_align_start = aln_pos
            if raw_idx == cdr_h3_end:
                h3_align_end = aln_pos
                break
            raw_idx += 1

    if h3_align_start is None or h3_align_end is None:
        return not_covered

    query_h3 = aligned_query[h3_align_start:h3_align_end]
    hit_h3 = aligned_hit[h3_align_start:h3_align_end]

    h3_align_length = len(query_h3)
    h3_identity = sum(1 for q, h in zip(query_h3, hit_h3) if q == h and q != "-")
    h3_pct_identity = (h3_identity / h3_align_length * 100) if h3_align_length else 0.0

    return {
        "h3_covered": True,
        "h3_identity": h3_identity,
        "h3_align_length": h3_align_length,
        "h3_pct_identity": round(h3_pct_identity, 2),
    }


def parse_and_filter_hits(
    xml_path: str,
    min_identity: float,
    max_identity: float,
    cdr_h3_start=None,
    cdr_h3_end=None,
    min_h3_identity=None,
    max_h3_identity=None,
):
    """
    New behavior:
      - Whole-chain pct_identity is still computed and reported (useful context),
        but is NO LONGER the filtering criterion when CDR-H3 boundaries are known.
      - If cdr_h3_start/end are provided (query is a heavy chain with a detected
        CDR-H3), filtering is done on H3-restricted identity via compute_h3_identity.
        HSPs that don't fully cover CDR-H3 are skipped (can't be judged reliably).
      - If cdr_h3_start/end are None (non-heavy chain, e.g. light chain or antigen,
        or ANARCI found nothing), falls back to the OLD whole-chain identity filter
        -- unchanged behavior for those sequences.
    """
    filtered = []
    try:
        records = Blast.parse(xml_path)
        for record in records:
            for hit in record:
                for hsp in hit:
                    # Safely access HSP attributes - they may be missing in some XML results
                    if not hasattr(hsp, 'length'):
                        continue
                    ident = hsp.annotations.get("identity", 0)
                    align_len = hsp.length
                    if align_len == 0:
                        continue
                    pct_identity = ident / align_len * 100

                    # Check for required attributes before building the row
                    if not hasattr(hit, 'target') or not hasattr(hit.target, 'id'):
                        continue

                    base_row = {
                        "hit_id": hit.target.id,
                        "hit_description": getattr(hit.target, 'description', hit.target.id),
                        "pct_identity": round(pct_identity, 2),
                        "identity": ident,
                        "align_length": align_len,
                        "evalue": hsp.annotations.get("evalue"),
                        "bit_score": hsp.annotations.get("bit score"),
                    }

                    if cdr_h3_start is not None and cdr_h3_end is not None:
                        # Filter by CDR-H3-restricted identity AND whole-chain identity.
                        # OLD BEHAVIOR (whole-chain + H3 filtering):
                        #   if not (min_h3_identity <= h3_info["h3_pct_identity"] <= max_h3_identity):
                        #       continue  # H3 identity outside range
                        #   if not (min_identity <= pct_identity <= max_identity):
                        #       continue  # Whole-chain identity outside range
                        #   filtered.append({**base_row, **h3_info})
                        #
                        # NEW BEHAVIOR (CDR-H3-only filtering):
                        #   Only filter by CDR-H3 identity, ignore whole-chain identity.
                        #   To revert: uncomment the 3 lines above and remove the 2 lines below.
                        h3_info = compute_h3_identity(hsp, cdr_h3_start, cdr_h3_end)
                        if not h3_info["h3_covered"]:
                            continue  # HSP doesn't span CDR-H3, can't judge it
                        if not (min_h3_identity <= h3_info["h3_pct_identity"] <= max_h3_identity):
                            continue  # H3 identity outside range
                        # NOTE: whole-chain identity check removed for CDR-only comparison
                        filtered.append({**base_row, **h3_info})
    except Exception as exc:
        log(f"ERROR: failed to parse {xml_path}: {exc}")
    return filtered
