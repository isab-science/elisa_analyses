# ELISA Analyses

This repository contains a generalized ELISA analysis script for comparing an antigen-coated plate with a matched blocked control plate.

## What It Expects

- Two `.xlsx` workbooks: one for the antigen of interest and one for the blocked control.
- The first worksheet in each workbook is used.
- The `patient_id` column identifies samples.
- A column containing `LOG10(EC50)` supplies the fitted EC50 value.
- Measurement columns default to `K:N`.
- The row-1 headers of `K:N` are interpreted as dilution fractions.
  `0.02` is treated as a 1:50 starting dilution, and the script converts those header values into dilution factors.

## What It Does

- Merges antigen and blocked measurements by sample ID.
- Converts the exported EC50 term back into dilution-factor units.
- Falls back to an EC50 estimate from the selected measurements if an exported EC50 is missing or clearly failed.
- Averages semicolon-separated duplicate values inside cells.
- Averages duplicate sample rows after grouping by sample ID.
- Produces EC50 tables, scatter plots, histograms, exploratory plots, and figure legends.

## Usage

```powershell
python fraa_reproducible_analysis.py --antigen-input "C:\path\to\curves_APOER2.xlsx" --blocked-input "C:\path\to\curves_Blocked.xlsx" --antigen-name "ApoER2"
```

You can also run the script without input arguments and pick both files from file dialogs:

```powershell
python fraa_reproducible_analysis.py
```

Optional arguments:

- `--measurement-columns K:N`
- `--outdir path\to\output_folder`

## Outputs

The script writes an output folder containing:

- `parsed_ec50_data.csv`
- `patient_ec50_table.csv`
- `patient_ec50_table.xlsx`
- `figure_1_histogram_<antigen>_ec50.svg`
- `figure_1_histogram_<antigen>_ec50_with_gaussian.svg`
- `figure_2_scatter_grey_subset.svg`
- `figure_2_scatter_with_uncertainty.svg`
- `figure_2_scatter_with_uncertainty_source_data.csv`
- `exploratory_*.svg`
- `exploratory_cohort_medians.csv`
- `detailed_legends.txt`

## Validated Example

The script was tested locally with:

- `C:\Users\aag\Downloads\curves_APOER2.xlsx`
- `C:\Users\aag\Downloads\curves_Blocked.xlsx`

using `--antigen-name ApoER2`.
