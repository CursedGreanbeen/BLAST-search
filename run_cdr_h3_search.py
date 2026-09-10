#!/usr/bin/env python3
"""
CLI entry point for local BLASTP + ANARCI CDR-H3 batch search.
See blast_pdb_search_local.py for the underlying implementation.
"""

import argparse
import csv
import sys
import time
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

from blast_pdb_search_local import (
    log,
    find_blastp,
    check_database,
    check_anarci,
    collect_unique_sequences,
    run_single_blastp,
    run_anarci,
    parse_and_filter_hits,
    DEFAULT_MIN_IDENTITY,
    DEFAULT_MAX_IDENTITY,
    DEFAULT_EVALUE,
    DEFAULT_MIN_H3_IDENTITY,
    DEFAULT_MAX_H3_IDENTITY,
)

CSV_FIELDNAMES = [
    "source_file", "chain_id", "hit_id", "hit_description",
    "pct_identity", "identity", "align_length", "evalue", "bit_score",
    "h3_covered", "h3_identity", "h3_align_length", "h3_pct_identity",
]


def main():
    parser = argparse.ArgumentParser(
        description="Batch local BLASTP search of FASTA complexes against a local PDB database."
    )
    parser.add_argument("fasta_dir", type=str, help="Directory containing .fasta files")
    parser.add_argument(
        "--db", type=str, required=True,
        help="Path to local BLAST database (e.g. ~/blast_local/db/pdbaa)"
    )
    parser.add_argument(
        "--workers", type=int, default=4,
        help="Number of blastp processes to run concurrently (default: 4)"
    )
    parser.add_argument(
        "--threads-per-search", type=int, default=2,
        help="Threads passed to each blastp process via -num_threads (default: 2). "
             "workers * threads_per_search should not exceed your CPU core count."
    )
    parser.add_argument(
        "--min-identity", type=float, default=DEFAULT_MIN_IDENTITY,
        help=f"Minimum percent identity to keep a hit (default: {DEFAULT_MIN_IDENTITY})"
    )
    parser.add_argument(
        "--max-identity", type=float, default=DEFAULT_MAX_IDENTITY,
        help=f"Maximum percent identity to keep a hit (default: {DEFAULT_MAX_IDENTITY})"
    )
    parser.add_argument(
        "--min-h3-identity", type=float, default=DEFAULT_MIN_H3_IDENTITY,
        help=f"Minimum CDR-H3 percent identity to keep a hit (default: {DEFAULT_MIN_H3_IDENTITY})"
    )
    parser.add_argument(
        "--max-h3-identity", type=float, default=DEFAULT_MAX_H3_IDENTITY,
        help=f"Maximum CDR-H3 percent identity to keep a hit (default: {DEFAULT_MAX_H3_IDENTITY})"
    )
    parser.add_argument(
        "--require-h3", action="store_true",
        help="Skip sequences where ANARCI found no heavy-chain CDR-H3 (no whole-chain fallback)."
    )
    parser.add_argument(
        "--evalue", type=float, default=DEFAULT_EVALUE,
        help=f"E-value threshold passed to blastp (default: {DEFAULT_EVALUE})"
    )
    parser.add_argument(
        "--output-dir", type=str, default="results/blast_pdb_search_local",
        help="Directory to write results into"
    )
    args = parser.parse_args()

    fasta_dir = Path(args.fasta_dir)
    if not fasta_dir.is_dir():
        sys.exit(f"ERROR: '{fasta_dir}' is not a directory.")

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    blastp_path = find_blastp()
    check_database(args.db)
    check_anarci()

    unique_seqs = collect_unique_sequences(fasta_dir)
    if not unique_seqs:
        sys.exit("No sequences found to BLAST. Exiting.")

    log(f"Running ANARCI on {len(unique_seqs)} unique sequence(s)...")
    anarci_results = {}
    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        future_to_hash = {
            executor.submit(run_anarci, data["sequence"]): seq_hash
            for seq_hash, data in unique_seqs.items()
        }
        for future in as_completed(future_to_hash):
            seq_hash = future_to_hash[future]
            try:
                anarci_results[seq_hash] = future.result()
            except Exception as exc:
                log(f"WARNING: ANARCI worker crashed for {seq_hash}: {exc}")
                anarci_results[seq_hash] = {
                    "chain_type": None, "domain_start": None, "domain_end": None,
                    "cdr_h3_start": None, "cdr_h3_end": None, "cdr_h3_sequence": None,
                }

    n_heavy = sum(1 for r in anarci_results.values() if r["cdr_h3_start"] is not None)
    log(f"ANARCI found CDR-H3 in {n_heavy}/{len(unique_seqs)} sequence(s).")

    h3_groups = defaultdict(list)  # cdr_h3_sequence -> list of seq_hash
    h3_to_pdbs = defaultdict(set)  # cdr_h3_sequence -> set of PDB IDs (from filtered hits)
    for seq_hash, result in anarci_results.items():
        h3_seq = result["cdr_h3_sequence"]
        if h3_seq:
            h3_groups[h3_seq].append(seq_hash)

    n_unique_h3 = len(h3_groups)
    n_shared = sum(1 for hashes in h3_groups.values() if len(hashes) > 1)
    log(f"{n_unique_h3} unique CDR-H3 sequence(s) found; {n_shared} shared by multiple chains.")

    if args.require_h3:
        skipped = [h for h, r in anarci_results.items() if r["cdr_h3_start"] is None]
        for seq_hash in skipped:
            del unique_seqs[seq_hash]
        if skipped:
            log(f"--require-h3 set: skipping {len(skipped)} sequence(s) with no detected CDR-H3.")

    tasks = [
        (
            seq_hash,
            data["sequence"],
            blastp_path,
            args.db,
            args.evalue,
            args.threads_per_search,
            str(output_dir),
        )
        for seq_hash, data in unique_seqs.items() if anarci_results[seq_hash]["cdr_h3_start"] is not None
    ]

    start_time = time.time()
    errors = []
    combined_rows = []

    log(
        f"Launching {len(tasks)} local blastp job(s) with "
        f"{args.workers} worker process(es) "
        f"({args.threads_per_search} thread(s) each)..."
    )

    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        future_to_hash = {
            executor.submit(run_single_blastp, task): task[0] for task in tasks
        }
        completed = 0
        for future in as_completed(future_to_hash):
            seq_hash = future_to_hash[future]
            completed += 1
            try:
                seq_hash, ok, info = future.result()
            except Exception as exc:
                ok, info = False, f"worker crashed: {exc}"

            sources = unique_seqs[seq_hash]["sources"]
            source_desc = "; ".join(f"{f}:{c}" for f, c, _ in sources)

            if not ok:
                log(f"[{completed}/{len(tasks)}] FAILED {seq_hash} ({source_desc}): {info}")
                errors.append(f"{seq_hash} ({source_desc}): {info}")
                continue

            xml_path = info
            log(f"[{completed}/{len(tasks)}] done {seq_hash} ({source_desc}) -> {xml_path}")

            h3_info = anarci_results[seq_hash]
            filtered_hits = parse_and_filter_hits(
                xml_path,
                args.min_identity, args.max_identity,
                cdr_h3_start=h3_info["cdr_h3_start"],
                cdr_h3_end=h3_info["cdr_h3_end"],
                min_h3_identity=args.min_h3_identity,
                max_h3_identity=args.max_h3_identity,
            )

            # Collect PDB IDs for this CDR-H3 sequence
            h3_seq = anarci_results[seq_hash]["cdr_h3_sequence"]
            if h3_seq:
                for hit in filtered_hits:
                    h3_to_pdbs[h3_seq].add(hit["hit_id"])

            per_seq_csv = output_dir / f"{seq_hash}_filtered_hits.csv"
            with open(per_seq_csv, "w", newline="") as f:
                writer = csv.DictWriter(f, fieldnames=CSV_FIELDNAMES)
                writer.writeheader()
                for hit in filtered_hits:
                    for file_name, chain_id, _desc in sources:
                        row = {"source_file": file_name, "chain_id": chain_id, **hit}
                        writer.writerow(row)
                        combined_rows.append(row)

    combined_csv = output_dir / "combined_filtered_hits.csv"
    with open(combined_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=CSV_FIELDNAMES)
        writer.writeheader()
        writer.writerows(combined_rows)

    h3_summary_csv = output_dir / "unique_cdr_h3_summary.csv"
    with open(h3_summary_csv, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["cdr_h3_sequence", "n_chains_sharing", "pdb_ids", "seq_hashes", "example_sources"])
        for h3_seq, hashes in h3_groups.items():
            pdb_ids = ";".join(sorted(h3_to_pdbs.get(h3_seq, set())))
            example_sources = "; ".join(
                f"{fn}:{cid}" for fn, cid, _ in unique_seqs[hashes[0]]["sources"]
            )
            writer.writerow([h3_seq, len(hashes), pdb_ids, ";".join(hashes), example_sources])

    if errors:
        with open(output_dir / "errors.txt", "w") as f:
            f.write("\n".join(errors))

    elapsed = time.time() - start_time
    log(
        f"Done in {elapsed:.1f}s. "
        f"{len(tasks) - len(errors)}/{len(tasks)} sequences succeeded. "
        f"{len(combined_rows)} filtered hit row(s) written to {combined_csv}."
    )
    if errors:
        log(f"{len(errors)} sequence(s) failed, see {output_dir / 'errors.txt'}")


if __name__ == "__main__":
    main()