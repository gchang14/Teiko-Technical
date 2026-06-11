"""
load_data.py  –  Part 1: Data Management
=========================================
Initialises teiko.db (SQLite) and loads cell-count.csv.

Schema
------
  projects    one row per project
  subjects    one row per subject; FK → projects
  samples     one row per biological sample; FK → subjects
  cell_counts long-format cell counts: one row per (sample, population)

Normalising subjects out of the sample rows avoids repeating age/sex/
condition/treatment/response for every timepoint.  Long-format cell_counts
means adding new populations never requires a schema change.  All columns
used in WHERE / JOIN / GROUP BY carry an index.
"""

import csv, os, sqlite3

DB_PATH  = "teiko.db"
CSV_PATH = "cell-count.csv"
POPULATIONS = ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]

DDL = """
PRAGMA foreign_keys = ON;

CREATE TABLE IF NOT EXISTS projects (
    project_id  TEXT PRIMARY KEY
);

CREATE TABLE IF NOT EXISTS subjects (
    subject_id  TEXT PRIMARY KEY,
    project_id  TEXT NOT NULL REFERENCES projects(project_id),
    condition   TEXT NOT NULL,
    age         INTEGER,
    sex         TEXT,
    treatment   TEXT,
    response    TEXT        -- 'yes' | 'no' | NULL (healthy/untreated)
);

CREATE TABLE IF NOT EXISTS samples (
    sample_id                 TEXT PRIMARY KEY,
    subject_id                TEXT NOT NULL REFERENCES subjects(subject_id),
    sample_type               TEXT NOT NULL,
    time_from_treatment_start INTEGER NOT NULL
);

CREATE TABLE IF NOT EXISTS cell_counts (
    id          INTEGER PRIMARY KEY AUTOINCREMENT,
    sample_id   TEXT    NOT NULL REFERENCES samples(sample_id),
    population  TEXT    NOT NULL,
    count       INTEGER NOT NULL
);

CREATE INDEX IF NOT EXISTS idx_sub_project   ON subjects(project_id);
CREATE INDEX IF NOT EXISTS idx_sub_condition ON subjects(condition);
CREATE INDEX IF NOT EXISTS idx_sub_treatment ON subjects(treatment);
CREATE INDEX IF NOT EXISTS idx_sub_response  ON subjects(response);
CREATE INDEX IF NOT EXISTS idx_sub_sex       ON subjects(sex);

CREATE INDEX IF NOT EXISTS idx_smp_subject   ON samples(subject_id);
CREATE INDEX IF NOT EXISTS idx_smp_type      ON samples(sample_type);
CREATE INDEX IF NOT EXISTS idx_smp_time      ON samples(time_from_treatment_start);

CREATE INDEX IF NOT EXISTS idx_cc_sample     ON cell_counts(sample_id);
CREATE INDEX IF NOT EXISTS idx_cc_population ON cell_counts(population);
"""

def load(db_path=DB_PATH, csv_path=CSV_PATH):
    if os.path.exists(db_path):
        os.remove(db_path)

    conn = sqlite3.connect(db_path)
    conn.executescript(DDL)

    projects_seen, subjects_seen = set(), set()
    cc_rows = []

    with open(csv_path, newline="") as fh:
        for row in csv.DictReader(fh):
            pid, sid, smid = row["project"], row["subject"], row["sample"]

            if pid not in projects_seen:
                conn.execute("INSERT OR IGNORE INTO projects VALUES (?)", (pid,))
                projects_seen.add(pid)

            if sid not in subjects_seen:
                conn.execute(
                    "INSERT OR IGNORE INTO subjects VALUES (?,?,?,?,?,?,?)",
                    (sid, pid, row["condition"],
                     int(row["age"]) if row["age"] else None,
                     row["sex"], row["treatment"],
                     row["response"] if row["response"] else None),
                )
                subjects_seen.add(sid)

            conn.execute(
                "INSERT OR IGNORE INTO samples VALUES (?,?,?,?)",
                (smid, sid, row["sample_type"],
                 int(row["time_from_treatment_start"])),
            )

            for pop in POPULATIONS:
                cc_rows.append((smid, pop, int(row[pop])))

    conn.executemany(
        "INSERT INTO cell_counts (sample_id, population, count) VALUES (?,?,?)",
        cc_rows,
    )
    conn.commit()

    cur = conn.cursor()
    stats = {t: cur.execute(f"SELECT COUNT(*) FROM {t}").fetchone()[0]
             for t in ("projects","subjects","samples","cell_counts")}
    conn.close()

    print(f"✓ Database ready: {db_path}")
    for t, n in stats.items():
        print(f"  {t:12s}: {n:,}")

if __name__ == "__main__":
    load()
