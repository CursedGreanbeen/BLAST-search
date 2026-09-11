from dataclasses import dataclass
from rapidfuzz.distance import Levenshtein


@dataclass(frozen=True)
class LengthBucket:
    min_len: int
    max_len: int  # inclusive
    max_edits: int


# Default bucket table encoding "short H3 => near-zero tolerance".
# Bounds are inclusive on both ends. Tune freely — this is just a first cut.
DEFAULT_BUCKETS = (
    LengthBucket(1, 4, 0),
    LengthBucket(5, 8, 1),
    LengthBucket(9, 14, 2),
    LengthBucket(15, 999, 3),
)


def max_allowed_edits(length: int, buckets=DEFAULT_BUCKETS) -> int:
    """Return the max allowed edit distance for a given H3 length,
    based on which bucket it falls into."""
    for bucket in buckets:
        if bucket.min_len <= length <= bucket.max_len:
            return bucket.max_edits
    # Length didn't fall into any bucket (gap in table or negative/zero length) —
    # fail loudly rather than silently allowing/disallowing a match.
    raise ValueError(f"No bucket defined for length={length}")


def edit_distance(a: str, b: str) -> int:
    """Raw Levenshtein distance between two (unaligned) H3 substrings."""
    return Levenshtein.distance(a, b)


def is_same_h3(query_h3: str, hit_h3: str, buckets=DEFAULT_BUCKETS) -> bool:
    """Decide whether query_h3 and hit_h3 count as 'the same' CDR-H3,
    using a length-bucketed absolute-edit-count threshold.

    Bucket is chosen from the shorter of the two sequences, so a short
    H3 can't inherit a more lenient tolerance just because it happens
    to be compared against a longer one.
    """
    if not query_h3 or not hit_h3:
        raise ValueError("H3 sequences must be non-empty")

    bucket_length = min(len(query_h3), len(hit_h3))
    threshold = max_allowed_edits(bucket_length, buckets)
    distance = edit_distance(query_h3, hit_h3)

    return distance <= threshold