
# fragr <img src='man/figures/logo.png' align="right" height="139" />

## Overview

`fragr` is an R package explicitly for wrangling capillary
electrophoresis data that has been called in and exported from
PeakScanner.

  - `read_PeakScanner()` imports csv-format PeakScanner *combined table*
    files into a list
  - `subfix()` repairs instances where the substrate peak is called
    twice (or more) from shoulders
  - `prodify()` frames relative peak area in terms of product peak area
    only
  - `act_specific()` defines specific activity
  - `act_offsense()` defines all other activity on the same strand
  - `act_offanti()` for two-dye experiments, defines activity on the
    opposite strand
