from pathlib import Path

import polars as pl
import pytest
from fragr.import_data import read_peakscanner

FIXTURES = Path(__file__).resolve().parents[2] / "tests" / "fixtures"


def test_reads_legacy_csv_with_cr_only_line_endings():
    df = read_peakscanner(FIXTURES / "PeakScanner_legacy.csv", type="legacy")
    assert df.shape == (23305, 20)
    assert df["ABIF_name"].null_count() == 0


def test_reads_legacy_tsv_matches_csv_counterpart():
    csv_df = read_peakscanner(FIXTURES / "PeakScanner_legacy.csv", type="legacy")
    tsv_df = read_peakscanner(
        FIXTURES / "PeakScanner_legacy.tsv", type="legacy", delim="\t"
    )
    assert tsv_df.shape == csv_df.shape
    assert tsv_df.drop("csv_name").equals(csv_df.drop("csv_name"))


def test_reads_agoactivity_csv():
    df = read_peakscanner(FIXTURES / "PeakScanner_agoactivity.csv", type="legacy")
    assert df.shape == (3791, 20)


def test_reads_real_legacy_export():
    df = read_peakscanner(
        FIXTURES / "PeakScanner1.0_legacy.tsv", type="legacy", delim="\t"
    )
    assert df.shape == (7293, 20)
    assert df["ABIF_name"].n_unique() == 16


def test_reads_real_cloud_export_and_preserves_ladder_peak_sizes():
    df = read_peakscanner(FIXTURES / "PeakScanner3.0_cloud.csv", type="cloud")
    assert df.shape == (2073, 13)
    assert df["ABIF_name"].n_unique() == 8

    # PeakScanner Cloud prefixes an asterisk on Size for peaks it used as
    # size-standard calibration points (e.g. "*35.00"); these must be kept
    # as their real numeric value, not dropped to null/0.
    ladder_35bp = df.filter((pl.col("dye") == "O") & (pl.col("bp") == 35.0))
    assert ladder_35bp.height == 8


def test_mismatched_export_raises():
    # PeakScanner_cloud_fail.csv is actually a sample-list export, not a
    # combined table: it has none of the peak-level columns either format
    # requires, so it must fail rather than silently returning junk.
    with pytest.raises(pl.exceptions.ColumnNotFoundError):
        read_peakscanner(FIXTURES / "PeakScanner_cloud_fail.csv", type="cloud")
    with pytest.raises(pl.exceptions.ColumnNotFoundError):
        read_peakscanner(FIXTURES / "PeakScanner_cloud_fail.csv", type="legacy")


def test_invalid_type_raises():
    with pytest.raises(ValueError):
        read_peakscanner(FIXTURES / "PeakScanner_legacy.csv", type="bogus")


def test_missing_path_raises():
    with pytest.raises(FileNotFoundError):
        read_peakscanner(FIXTURES / "does_not_exist.csv", type="legacy")


def test_directory_import_returns_dict_keyed_by_path():
    result = read_peakscanner(
        FIXTURES, pattern=r"PeakScanner_(legacy|agoactivity)\.csv$", type="legacy"
    )
    assert set(result.keys()) == {
        str(FIXTURES / "PeakScanner_agoactivity.csv"),
        str(FIXTURES / "PeakScanner_legacy.csv"),
    }
    for df in result.values():
        assert df.shape[1] == 20
