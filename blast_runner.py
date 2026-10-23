"""Local blastp execution: finding the binary, validating the database,
and running individual searches as subprocesses.

Pure "talk to the blastp/blastdbcmd command-line tools" concerns. No ANARCI,
no FASTA parsing, no CDR-H3 logic in here.
"""

import shutil
import subprocess
import sys
from pathlib import Path

from logging_utils import log

DEFAULT_MIN_IDENTITY = 90.0
DEFAULT_MAX_IDENTITY = 99.0
DEFAULT_EVALUE = 10.0


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
