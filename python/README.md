# fragr (Python)

## Overview

`fragr` is a Python package for wrangling capillary electrophoresis data
that has been called in and exported from PeakScanner. It's a Python port
of the [`fragr` R package](../README.md), built incrementally alongside it.

Currently implemented:

  - `read_peakscanner()` reads a single PeakScanner *combined table* file,
    or a whole directory of them, into [polars](https://pola.rs) DataFrames.
    It auto-detects file vs. directory input and supports both the
    'legacy' (desktop app) and 'cloud' (online app) export column layouts.

Not yet ported from the R package: `subfix()`, `prodify()`,
`act_specific()`, `act_offsense()`, `act_offanti()`.

## Installation

This package uses [uv](https://docs.astral.sh/uv/) for dependency
management. From the `python/` directory:

```sh
uv sync
```

This installs `fragr` in editable mode into `.venv`, along with its
dependencies (`polars`) and dev dependencies (`pytest`).

## Usage

```python
from fragr.import_data import read_peakscanner

# a single combined-table export
df = read_peakscanner("path/to/export.csv", type="legacy")

# a directory of exports, matched by regex (defaults to all *.csv files)
tables = read_peakscanner("path/to/exports/", type="cloud", pattern=r"\.csv$")

# tab-delimited exports
df = read_peakscanner("path/to/export.tsv", type="legacy", delim="\t")
```

## Testing

```sh
cd python
uv run pytest -v
```

Tests live in `tests/` and read fixture PeakScanner exports from
`tests/fixtures/`.
