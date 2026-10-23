"""ANARCI integration: availability check + CDR-H3 extraction.

Everything here concerns running ANARCI on a single sequence and mapping its
IMGT-numbered output back onto raw sequence indices to find CDR-H3
boundaries. No BLAST or file I/O concerns live in this module.
"""

import sys

from anarci import anarci as _anarci_run

from logging_utils import log

HEAVY_CHAIN_TYPES = {"H"}  # ANARCI chain_type label for heavy chains

IMGT_CDR_H3_START = 105
IMGT_CDR_H3_END = 117  # inclusive, per IMGT numbering scheme


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
