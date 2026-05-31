#!/usr/bin/env python3
"""
tsv_to_parquet.py
-----------------
Convert a delimited flat file (TSV, CSV, or .gz variants) to a Parquet file
or Hive-partitioned Parquet directory using DuckDB.

Examples
--------
# Simple conversion, no partitioning
python tsv_to_parquet.py input.tsv.gz output.parquet

# Partition by one column
python tsv_to_parquet.py input.tsv.gz output_dir/ --partition-by chrom

# Partition by multiple columns
python tsv_to_parquet.py input.tsv.gz output_dir/ --partition-by chrom strand

# Override delimiter and force specific column types
python tsv_to_parquet.py input.csv.gz output.parquet \\
    --delimiter , \\
    --dtypes chrom:VARCHAR pos:INTEGER score:DOUBLE

# Dry run — print inferred schema and row count, write nothing
python tsv_to_parquet.py input.tsv.gz output.parquet --dry-run

# Tune Parquet row group size and compression codec
python tsv_to_parquet.py input.tsv.gz output.parquet \\
    --row-group-size 500000 \\
    --compression zstd
"""

import argparse
import sys
from pathlib import Path

try:
    import duckdb
except ImportError:
    sys.exit("duckdb is required: pip install duckdb")


# ---------------------------------------------------------------------------
# Delimiter inference
# ---------------------------------------------------------------------------

def infer_delimiter(path: str, explicit: str | None) -> str:
    if explicit is not None:
        # Allow the user to write 'tab' or '\\t' on the CLI
        if explicit.lower() == "tab":
            return "\t"
        return explicit.replace("\\t", "\t")
    stem = Path(path).stem  # strips one extension; handles .tsv.gz → .tsv
    if Path(stem).suffix.lower() in (".tsv", ".txt"):
        return "\t"
    return ","  # default to comma for .csv and unknowns


# ---------------------------------------------------------------------------
# dtype override parsing  (e.g. "chrom:VARCHAR" "pos:INTEGER")
# ---------------------------------------------------------------------------

def parse_dtypes(dtype_args: list[str] | None) -> dict[str, str]:
    if not dtype_args:
        return {}
    result = {}
    for token in dtype_args:
        if ":" not in token:
            sys.exit(f"--dtypes entries must be 'column:TYPE', got: {token!r}")
        col, dtype = token.split(":", 1)
        result[col.strip()] = dtype.strip().upper()
    return result


def dtype_map_to_sql(dtype_map: dict[str, str]) -> str:
    """Convert {col: TYPE} dict to DuckDB columns= syntax fragment."""
    if not dtype_map:
        return ""
    entries = ", ".join(f"'{col}': '{dtype}'" for col, dtype in dtype_map.items())
    return f", columns={{{entries}}}"


# ---------------------------------------------------------------------------
# Validation helpers
# ---------------------------------------------------------------------------

def validate_output_path(output: str, partition_by: list[str] | None) -> None:
    """Exit early with a clear message if the output path is inconsistent."""
    path = Path(output)
    if not partition_by and path.suffix.lower() != ".parquet":
        sys.exit(
            f"Error: output path {output!r} does not end in .parquet.\n"
            "When --partition-by is not supplied the output must be a single\n"
            "Parquet file. Either add a .parquet extension or supply --partition-by."
        )


def validate_partition_columns(
    con: duckdb.DuckDBPyConnection,
    read_clause: str,
    partition_by: list[str],
) -> None:
    """Check that every requested partition column exists in the input schema."""
    rows = con.execute(f"DESCRIBE SELECT * FROM {read_clause} LIMIT 0").fetchall()
    available = {r[0] for r in rows}
    missing = [c for c in partition_by if c not in available]
    if missing:
        sys.exit(
            f"Error: partition column(s) not found in input: {', '.join(missing)}\n"
            f"Available columns: {', '.join(sorted(available))}"
        )


# ---------------------------------------------------------------------------
# Core logic
# ---------------------------------------------------------------------------

def build_read_clause(input_path: str, delimiter: str, dtype_sql: str) -> str:
    escaped_delim = delimiter.replace("'", "\\'").replace("\t", "\\t")
    return (
        f"read_csv('{input_path}', "
        f"delim='{escaped_delim}', "
        f"header=true, "
        f"auto_detect=true"
        f"{dtype_sql})"
    )


def run_dry_run(con: duckdb.DuckDBPyConnection, read_clause: str) -> None:
    print("\n── Inferred Schema ──────────────────────────────────────")
    schema = con.execute(f"DESCRIBE SELECT * FROM {read_clause} LIMIT 0").fetchall()
    col_w = max(len(r[0]) for r in schema) + 2
    for row in schema:
        col, dtype = row[0], row[1]
        print(f"  {col:<{col_w}} {dtype}")

    print("\n── Row Count ────────────────────────────────────────────")
    count = con.execute(f"SELECT COUNT(*) FROM {read_clause}").fetchone()[0]
    print(f"  {count:,} rows")
    print("\n(dry run — nothing written)\n")


def run_conversion(
    con: duckdb.DuckDBPyConnection,
    read_clause: str,
    output: str,
    partition_by: list[str] | None,
    row_group_size: int,
    compression: str,
) -> None:
    output_path = Path(output)

    if partition_by:
        # Partitioned output → directory
        output_path.mkdir(parents=True, exist_ok=True)
        partition_cols = ", ".join(partition_by)
        sql = f"""
            COPY (SELECT * FROM {read_clause})
            TO '{output_path}'
            (
                FORMAT PARQUET,
                PARTITION_BY ({partition_cols}),
                ROW_GROUP_SIZE {row_group_size},
                COMPRESSION {compression},
                OVERWRITE_OR_IGNORE true
            )
        """
    else:
        # Flat single-file output
        output_path.parent.mkdir(parents=True, exist_ok=True)
        sql = f"""
            COPY (SELECT * FROM {read_clause})
            TO '{output_path}'
            (
                FORMAT PARQUET,
                ROW_GROUP_SIZE {row_group_size},
                COMPRESSION {compression}
            )
        """

    print(f"Writing to: {output_path}")
    if partition_by:
        print(f"Partitioning by: {', '.join(partition_by)}")
    print(f"Compression: {compression}  |  Row group size: {row_group_size:,}")
    print("Running… ", end="", flush=True)

    con.execute(sql)
    print("done.")

    # Report output size
    if output_path.is_dir():
        total = sum(f.stat().st_size for f in output_path.rglob("*.parquet"))
        n_files = sum(1 for _ in output_path.rglob("*.parquet"))
        print(f"Output: {n_files} Parquet file(s), {total / 1_048_576:.1f} MB total")
    else:
        size = output_path.stat().st_size
        print(f"Output: {size / 1_048_576:.1f} MB")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="tsv_to_parquet.py",
        description="Convert a TSV/CSV (optionally .gz) to Parquet using DuckDB.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__.split("Examples")[1] if "Examples" in __doc__ else "",
    )

    parser.add_argument(
        "input",
        help="Input file: .tsv, .csv, .tsv.gz, .csv.gz, etc.",
    )
    parser.add_argument(
        "output",
        help=(
            "Output path. Use a .parquet extension for a single file, "
            "or a directory path when --partition-by is supplied."
        ),
    )
    parser.add_argument(
        "--partition-by",
        nargs="+",
        metavar="COLUMN",
        default=None,
        help="One or more column names to partition the output by (Hive-style).",
    )
    parser.add_argument(
        "--delimiter",
        default=None,
        metavar="CHAR",
        help=(
            "Field delimiter. Inferred from file extension if omitted "
            "(.tsv/.txt → TAB, .csv → comma). Use 'tab' or '\\t' for tab."
        ),
    )
    parser.add_argument(
        "--dtypes",
        nargs="+",
        metavar="COL:TYPE",
        default=None,
        help=(
            "Override inferred column types. Format: COLUMN:DUCKDB_TYPE. "
            "Example: --dtypes chrom:VARCHAR pos:INTEGER score:DOUBLE"
        ),
    )
    parser.add_argument(
        "--row-group-size",
        type=int,
        default=122_880,
        metavar="N",
        help="Parquet row group size (default: 122880 ≈ DuckDB default).",
    )
    parser.add_argument(
        "--compression",
        default="snappy",
        choices=["snappy", "zstd", "gzip", "lz4", "none"],
        help="Parquet compression codec (default: snappy).",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print inferred schema and row count without writing any output.",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=None,
        metavar="N",
        help="Number of DuckDB threads (default: all available cores).",
    )

    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    delimiter = infer_delimiter(args.input, args.delimiter)
    dtype_map = parse_dtypes(args.dtypes)
    dtype_sql = dtype_map_to_sql(dtype_map)
    read_clause = build_read_clause(args.input, delimiter, dtype_sql)

    # --- output path guard (no I/O needed) ---
    validate_output_path(args.output, args.partition_by)

    con = duckdb.connect()
    if args.threads is not None:
        con.execute(f"SET threads = {args.threads}")

    # --- partition column guard (requires reading input schema) ---
    if args.partition_by:
        validate_partition_columns(con, read_clause, args.partition_by)

    if args.dry_run:
        run_dry_run(con, read_clause)
    else:
        run_conversion(
            con,
            read_clause,
            args.output,
            args.partition_by,
            args.row_group_size,
            args.compression,
        )

    con.close()


if __name__ == "__main__":
    main()
