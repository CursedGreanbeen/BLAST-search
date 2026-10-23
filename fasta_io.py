"""FASTA parsing / I-O helpers.

Pure file-reading concerns: discovering FASTA files in a directory and
deduplicating sequences across files/chains before they get sent to BLAST
or ANARCI. No BLAST, no ANARCI, no subprocess calls in here.
"""

import hashlib
from pathlib import Path

from Bio import SeqIO

from logging_utils import log

FASTA_EXTENSIONS = (".fasta", ".fa", ".faa")


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
