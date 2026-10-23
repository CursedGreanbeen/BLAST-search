"""CDR-H3-aware BLAST hit filtering.

Extracts the aligned CDR-H3 substrings from a BLAST HSP (using the query's
known ANARCI-derived CDR-H3 boundaries) and decides whether a hit counts as
"the same" H3 using the length-bucketed edit-distance threshold from
normalize_cdr_h3_length, instead of a flat percent-identity cutoff.
"""

from Bio import Blast

from logging_utils import log
from normalize_cdr_h3_length import is_same_h3


def compute_h3_identity(hsp, cdr_h3_start, cdr_h3_end):
    """Slice out the query/hit CDR-H3 substrings from a BLAST HSP's gapped
    alignment, using the query's raw-sequence CDR-H3 boundaries (from ANARCI)
    mapped onto alignment-column coordinates.

    Still reports whole-alignment-style bookkeeping (h3_align_length,
    h3_identity, h3_pct_identity) for context/debugging, but the actual
    match/no-match decision now happens in parse_and_filter_hits via
    is_same_h3(), not via h3_pct_identity.
    """
    not_covered = {
        "h3_covered": False,
        "h3_identity": None,
        "h3_align_length": None,
        "h3_pct_identity": None,
        "query_h3": None,
        "hit_h3": None,
    }

    if cdr_h3_start is None or cdr_h3_end is None:
        return not_covered

    aligned_query = str(hsp[0])   # gapped alignment string, query row
    aligned_hit = str(hsp[1])     # gapped alignment string, hit/target row

    raw_idx = hsp.coordinates[0][0]   # raw query start, from coordinates array

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
        # Gap characters stripped before handing off to is_same_h3, since
        # edit distance should operate on raw residues, not alignment columns.
        "query_h3": query_h3.replace("-", ""),
        "hit_h3": hit_h3.replace("-", ""),
    }


def parse_and_filter_hits(
    xml_path: str,
    min_identity: float,
    max_identity: float,
    cdr_h3_start=None,
    cdr_h3_end=None,
):
    """
    Behavior:
      - Whole-chain pct_identity is still computed and reported (useful context),
        but is NOT the filtering criterion when CDR-H3 boundaries are known.
      - If cdr_h3_start/end are provided (query is a heavy chain with a detected
        CDR-H3), filtering is done via is_same_h3() on the raw (gap-stripped)
        query/hit CDR-H3 substrings extracted by compute_h3_identity. HSPs that
        don't fully cover CDR-H3 are skipped (can't be judged reliably).
      - If cdr_h3_start/end are None (non-heavy chain, e.g. light chain or
        antigen, or ANARCI found nothing), no CDR-H3 filtering is applied here
        -- caller is responsible for deciding what to do with such sequences
        (see --require-h3 in run_cdr_h3_search.py).
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
                        # Filter by CDR-H3 edit-distance match instead of a
                        # flat whole-chain or H3 percent-identity cutoff.
                        h3_info = compute_h3_identity(hsp, cdr_h3_start, cdr_h3_end)
                        if not h3_info["h3_covered"]:
                            continue  # HSP doesn't span CDR-H3, can't judge it
                        if not is_same_h3(h3_info["query_h3"], h3_info["hit_h3"]):
                            continue
                        filtered.append({**base_row, **h3_info})

    except Exception as exc:
        log(f"ERROR: failed to parse {xml_path}: {exc}")
    return filtered
