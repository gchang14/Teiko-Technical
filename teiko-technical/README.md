# Teiko Clinical Immune Profiling – Technical Assignment

## Overview

This project analyses immune cell population data from a clinical trial to understand how the drug candidate **miraclib** affects immune profiles.  The pipeline covers four parts:

| Part | Description | Output |
|------|-------------|--------|
| 1 | SQLite database initialisation & data loading | `teiko.db` |
| 2 | Relative cell-population frequencies per sample | `outputs/part2_frequencies.csv` |
| 3 | Statistical comparison of responders vs non-responders (melanoma · PBMC · miraclib) | `outputs/part3_boxplot.png`, `outputs/part3_stats.csv` |
| 4 | Baseline subset analysis + key numeric answer | `outputs/part4_subset.csv` |

---

## Quick Start (GitHub Codespaces)

```bash
make setup      # install dependencies
make pipeline   # run full analysis (Parts 1–4)
make dashboard  # start the interactive dashboard on port 8050
```

---

## Database Schema

```
projects    (project_id PK)
    │
subjects    (subject_id PK, project_id FK,
    │        condition, age, sex, treatment, response)
    │
samples     (sample_id PK, subject_id FK,
    │        sample_type, time_from_treatment_start)
    │
cell_counts (id PK AUTOINCREMENT, sample_id FK,
             population, count)
```

### Rationale

**Normalisation into four tables** avoids repeating subject-level attributes (age, sex, condition, treatment, response) across every sample row.  With thousands of subjects each having multiple timepoints, this keeps the data consistent and compact.

**Long-format `cell_counts`** stores one row per (sample, population) pair rather than one column per population.  Adding new cell types (e.g. `regulatory_t_cell`) never requires an `ALTER TABLE`; it is just new rows.  This also makes aggregate queries (`GROUP BY population`) trivial.

**Scalability to hundreds of projects / thousands of samples:**
- Indexes on `condition`, `treatment`, `response`, `sample_type`, and `time_from_treatment_start` keep analytic `WHERE` and `JOIN` clauses fast.
- The long-format design keeps `cell_counts` append-only; new panels or populations integrate without schema changes.
- If read throughput becomes a bottleneck, materialised views or a columnar store (DuckDB, Parquet) can be layered on top of the same schema.

---

## Code Structure

```
.
├── load_data.py       # Part 1 – schema + CSV ingestion
├── analysis.py        # Parts 2, 3, 4 – analysis & output files
├── dashboard.py       # Interactive Dash dashboard
├── cell-count.csv     # Input data
├── requirements.txt
├── Makefile
└── outputs/
    ├── part2_frequencies.csv
    ├── part3_boxplot.png
    ├── part3_stats.csv
    └── part4_subset.csv
```

**`load_data.py`** – single responsibility: schema creation and bulk load.  Runs in seconds even for large datasets because it batches all `cell_counts` inserts via `executemany`.

**`analysis.py`** – reads from the database, keeps logic for each Part in its own function (`part2_frequency_table`, `part3_stats`, `part4_subset`).  All outputs are written to `outputs/`.

**`dashboard.py`** – Dash application with three tabs (one per analytical part).  Reads from the same SQLite database; all filtering is done in-process with pandas.

---

## Key Results

### Part 3 – Statistical Analysis
Mann-Whitney U test (two-sided, α = 0.05), melanoma PBMC miraclib samples only:

| Population | Mean Resp | Mean Non-Resp | p-value | Significant? |
|------------|-----------|---------------|---------|--------------|
| B Cell     | 9.80%     | 10.00%        | 0.0557  | No |
| CD8 T Cell | 24.88%    | 24.94%        | 0.6391  | No |
| **CD4 T Cell** | **30.54%** | **29.90%** | **0.0133** | **Yes ✓** |
| NK Cell    | 14.84%    | 15.07%        | 0.1211  | No |
| Monocyte   | 19.94%    | 20.08%        | 0.1631  | No |

Only **CD4 T cells** show a statistically significant difference between responders and non-responders (p = 0.013).  Responders have a modestly higher CD4 T cell frequency, which is consistent with the known role of CD4 helper T cells in supporting anti-tumour immunity.

### Part 4 – Key Numeric Answer
**Average B cell count for melanoma male responders at t = 0: `10401.28`**  
(n = 184 samples)
