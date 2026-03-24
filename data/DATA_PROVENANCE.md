# Data Provenance

## Source

CDC PRAMStat (Pregnancy Risk Assessment Monitoring System Statistics)
Years: 2000–2011
Unit of analysis: State-year aggregate prevalence estimates
Geographic coverage: 40 US locations (states + DC)

## Original Download

Data was downloaded from catalog.data.gov CKAN API endpoints for each year
(2000–2011). Each download includes CSV and JSON formats. SHA-256 checksums
were computed at download time and stored in `download_manifest.csv`.

## Pipeline

1. **Ingestion**: Raw files copied to project tree with integrity verification
2. **Schema validation**: Column headers harmonized across years, CSV/JSON
   cross-checked for consistency, problematic years excluded
3. **Master dataset construction**: All yearly CSVs merged into long-format
   `analysis_master.csv` with quality flags and domain flags
4. **Panel construction**: Unstratified (overall) YES-response rows pivoted
   to state-year × question_id wide format

## Master Dataset Location

Parent project: `../analysis_output/pramstat_2000_2011/run_20260228T185246Z/analysis_master.csv`
Copied to this repository: `data/analysis_master.csv`

## Filtering Criteria for Panel

- `break_out` is NaN AND `break_out_category` is NaN (unstratified rows)
- `binary_response_class` == 'yes'
- Pivoted: rows = (location_abbr, year), columns = question_id, values = estimate
- Result: 329 state-year rows × 229 question_id columns

## Integrity

The master dataset hash and row/column counts are logged in every analysis
run's metadata JSON for verification.
