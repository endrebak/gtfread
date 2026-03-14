"""Public package API for gtfread."""

try:
    from ._parser import parse_chunk_columns
except ImportError:
    from ._fallback import parse_chunk_columns
from .readers import (
    find_first_data_line_index,
    parse_kv_fields,
    read_gtf,
    read_gtf_full,
    read_gtf_full_parallel,
    read_gtf_restricted,
    to_rows,
    to_rows_keep_duplicates,
)

__all__ = [
    "find_first_data_line_index",
    "parse_chunk_columns",
    "parse_kv_fields",
    "read_gtf",
    "read_gtf_full",
    "read_gtf_full_parallel",
    "read_gtf_restricted",
    "to_rows",
    "to_rows_keep_duplicates",
]
