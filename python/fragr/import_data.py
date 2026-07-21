"""Import PeakScanner *combined table* .csv files (file or directory)."""

import io
import re
from pathlib import Path

import polars as pl

_LEGACY_COLUMNS = {
    "Status": pl.Utf8,
    "Sample Name": pl.Utf8,
    "Sample Type": pl.Utf8,
    "Size Standard": pl.Utf8,
    "Analysis Method": pl.Utf8,
    "Offscale": pl.Utf8,
    "Quality": pl.Utf8,
    "Dye/Sample Peak": pl.Utf8,
    "Sample File Name": pl.Utf8,
    "Size": pl.Float64,
    "Height": pl.Float64,
    "Area in Point": pl.Float64,
    "Area in BP": pl.Float64,
    "Data Point": pl.Float64,
    "Begin Point": pl.Float64,
    "Begin BP": pl.Float64,
    "End Point": pl.Float64,
    "End BP": pl.Float64,
    "Width in Point": pl.Float64,
    "Width in BP": pl.Float64,
    "User Comments": pl.Utf8,
    "User Edit": pl.Utf8,
}

_LEGACY_RENAME = {
    "Status": "status",
    "Sample Name": "sample_id",
    "Sample Type": "type",
    "Size Standard": "standard",
    "Analysis Method": "method",
    "Offscale": "offscale",
    "Quality": "quality",
    "Dye/Sample Peak": "dye_peak",
    "Sample File Name": "ABIF_name",
    "Size": "bp",
    "Height": "height",
    "Area in Point": "scan_area",
    "Area in BP": "bp_area",
    "Data Point": "scan",
    "Begin Point": "scan_begin",
    "Begin BP": "bp_begin",
    "End Point": "scan_end",
    "End BP": "bp_end",
    "Width in Point": "scan_width",
    "Width in BP": "bp_width",
    "User Comments": "comment",
    "User Edit": "edit",
}

_LEGACY_ORDER = [
    "ABIF_name",
    "sample_id",
    "standard",
    "method",
    "quality",
    "offscale",
    "dye_peak",
    "height",
    "bp",
    "bp_begin",
    "bp_end",
    "bp_width",
    "bp_area",
    "scan",
    "scan_begin",
    "scan_end",
    "scan_width",
    "scan_area",
]

_CLOUD_COLUMNS = {
    "Sample Filename": pl.Utf8,
    "Dye Color": pl.Utf8,
    "Dye, Sample Peak": pl.Utf8,
    "Size": pl.Float64,
    "Height": pl.Float64,
    "Area (Data Point)": pl.Float64,
    "Area (Base Pairs)": pl.Float64,
    "Data Point": pl.Float64,
    "Begin Point (Data Point)": pl.Float64,
    "End Point (Data Point)": pl.Float64,
    "Begin Point (Base Pairs)": pl.Float64,
    "End Point (Base Pairs)": pl.Float64,
}

_CLOUD_RENAME = {
    "Sample Filename": "ABIF_name",
    "Dye Color": "dye_color",
    "Dye, Sample Peak": "dye_peak",
    "Size": "bp",
    "Height": "height",
    "Area (Data Point)": "scan_area",
    "Area (Base Pairs)": "bp_area",
    "Data Point": "scan",
    "Begin Point (Data Point)": "scan_begin",
    "End Point (Data Point)": "scan_end",
    "Begin Point (Base Pairs)": "bp_begin",
    "End Point (Base Pairs)": "bp_end",
}

_CLOUD_ORDER = [
    "ABIF_name",
    "dye_peak",
    "height",
    "bp",
    "bp_begin",
    "bp_end",
    "bp_area",
    "scan",
    "scan_begin",
    "scan_end",
    "scan_area",
]

_COLUMN_SPECS = {
    "legacy": {"columns": _LEGACY_COLUMNS, "rename": _LEGACY_RENAME, "order": _LEGACY_ORDER},
    "cloud": {"columns": _CLOUD_COLUMNS, "rename": _CLOUD_RENAME, "order": _CLOUD_ORDER},
}


def _final_order(order):
    result = ["csv_name"]
    for col in order:
        result += ["dye", "peak"] if col == "dye_peak" else [col]
    return result


def _read_one(file_path, spec, delim):
    columns = spec["columns"]
    with open(file_path, "rb") as fh:
        raw = fh.read()
    # Normalize line endings: PeakScanner desktop exports sometimes use bare
    # CR (classic Mac) endings, which polars' CSV reader treats as data
    # rather than a row separator.
    normalized = raw.replace(b"\r\n", b"\n").replace(b"\r", b"\n")
    # Read every target column as text: PeakScanner Cloud prefixes an
    # asterisk on Size for peaks it used as size-standard calibration
    # points (e.g. "*35.00"), which numeric parsing at read time can't handle.
    df = pl.read_csv(
        io.BytesIO(normalized),
        separator=delim,
        columns=list(columns.keys()),
        schema_overrides={name: pl.Utf8 for name in columns},
    )
    float_cols = [name for name, dtype in columns.items() if dtype == pl.Float64]
    df = df.with_columns(
        pl.col(name).str.strip_chars().str.strip_chars_start("*").cast(pl.Float64, strict=False)
        for name in float_cols
    )
    df = (
        df.rename(spec["rename"])
        .select(spec["order"])
        .with_columns(
            pl.col("dye_peak")
            .str.splitn(",", 2)
            .struct.rename_fields(["dye", "peak"])
        )
        .unnest("dye_peak")
        .with_columns(
            pl.col("peak").str.strip_chars().cast(pl.Int64, strict=False),
            pl.col("bp").fill_null(0),
            pl.lit(Path(file_path).name).alias("csv_name"),
        )
    )
    return df.select(_final_order(spec["order"]))


def read_peakscanner(path, pattern=None, type="legacy", delim=","):
    """Read one PeakScanner *combined table* .csv, or a whole directory of them.

    Args:
        path: path to a single .csv file, or a directory containing them.
        pattern: regex used to select files when `path` is a directory;
            defaults to all .csv files.
        type: 'legacy' (desktop app) or 'cloud' (online app) export format.
        delim: field delimiter, only used when `path` is a single file.

    Returns:
        A polars DataFrame if `path` is a file, or a dict mapping each
        matched file's path to its polars DataFrame if `path` is a directory.
    """
    if type not in ("legacy", "cloud"):
        raise ValueError("Arg `type` must be either 'legacy' (for desktop app) or 'cloud'.")

    spec = _COLUMN_SPECS[type]
    path = Path(path)

    if path.is_dir():
        regex = re.compile(pattern if pattern else r"\.csv")
        files = sorted(p for p in path.iterdir() if p.is_file() and regex.search(p.name))
        return {str(p): _read_one(p, spec, delim) for p in files}

    if path.is_file():
        return _read_one(path, spec, delim)

    raise FileNotFoundError(f"No such file or directory: {path}")
