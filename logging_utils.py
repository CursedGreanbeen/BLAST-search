"""Shared logging helper used across all modules in this pipeline.

Kept in its own module (rather than living inside any single module) so that
fasta_io, anarci_utils, blast_runner and h3_matching can all import the same
log() / _print_lock without creating circular imports between each other.
"""

import threading
import time

_print_lock = threading.Lock()


def log(msg: str) -> None:
    with _print_lock:
        print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)
