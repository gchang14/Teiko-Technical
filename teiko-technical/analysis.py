"""
analysis.py  –  Parts 2, 3 & 4
================================
Reads from teiko.db (created by load_data.py) and produces:

  Part 2 – summary table of relative cell-population frequencies per sample
            saved to  outputs/part2_frequencies.csv

  Part 3 – Mann-Whitney U tests comparing responders vs non-responders
            (melanoma, PBMC, miraclib only); boxplot saved to
            outputs/part3_boxplot.png  and stats to
            outputs/part3_stats.csv

  Part 4 – subset analysis: melanoma PBMC baseline (t=0) miraclib samples
            saved to  outputs/part4_subset.csv
            plus the specific numeric answer printed to stdout.
"""

import os, sqlite3
import pandas as pd
from scipy import stats as scipy_stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

DB_PATH = "teiko.db"
OUT_DIR = "outputs"
POPULATIONS = ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]

os.makedirs(OUT_DIR, exist_ok=True)


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def get_conn():
    conn = sqlite3.connect(DB_PATH)
    conn.row_factory = sqlite3.Row
    return conn


# ---------------------------------------------------------------------------
# Part 2 – frequency table
# ---------------------------------------------------------------------------

def part2_frequency_table(conn) -> pd.DataFrame:
    """
    For every sample compute total cells, then per-population count & %.
    Returns a long-format DataFrame with columns:
        sample, total_count, population, count, percentage
    """
    sql = """
        SELECT
            s.sample_id   AS sample,
            cc.population AS population,
            cc.count      AS count
        FROM cell_counts cc
        JOIN samples s ON s.sample_id = cc.sample_id
        ORDER BY s.sample_id, cc.population
    """
    df = pd.read_sql_query(sql, conn)

    totals = df.groupby("sample")["count"].sum().rename("total_count")
    df = df.join(totals, on="sample")
    df["percentage"] = (df["count"] / df["total_count"] * 100).round(4)

    df = df[["sample", "total_count", "population", "count", "percentage"]]

    out = os.path.join(OUT_DIR, "part2_frequencies.csv")
    df.to_csv(out, index=False)
    print(f"✓ Part 2: frequency table saved → {out}  ({len(df):,} rows)")
    return df


# ---------------------------------------------------------------------------
# Part 3 – statistical analysis
# ---------------------------------------------------------------------------

def part3_stats(conn, freq_df: pd.DataFrame):
    """
    Filter: melanoma, PBMC, miraclib, response in (yes/no).
    Test each population with Mann-Whitney U (two-sided).
    Produce boxplot + CSV of stats.
    """
    # Pull subject/sample metadata
    meta_sql = """
        SELECT
            s.sample_id,
            sub.condition,
            sub.treatment,
            sub.response,
            s.sample_type
        FROM samples s
        JOIN subjects sub ON sub.subject_id = s.subject_id
        WHERE sub.condition  = 'melanoma'
          AND sub.treatment  = 'miraclib'
          AND sub.response   IN ('yes','no')
          AND s.sample_type  = 'PBMC'
    """
    meta = pd.read_sql_query(meta_sql, conn)

    merged = freq_df.merge(meta, on="sample_id" if "sample_id" in freq_df.columns else None,
                           left_on="sample", right_on="sample_id")

    # ---- boxplot ----
    fig, axes = plt.subplots(1, 5, figsize=(18, 6))
    fig.suptitle(
        "Cell Population Frequencies: Responders vs Non-Responders\n"
        "(Melanoma · PBMC · Miraclib)",
        fontsize=13, fontweight="bold"
    )

    stat_rows = []
    resp_color    = "#4C9BE8"
    nonresp_color = "#E8724C"

    for ax, pop in zip(axes, POPULATIONS):
        pop_data = merged[merged["population"] == pop]
        yes_vals = pop_data[pop_data["response"] == "yes"]["percentage"].values
        no_vals  = pop_data[pop_data["response"] == "no" ]["percentage"].values

        stat, pval = scipy_stats.mannwhitneyu(yes_vals, no_vals, alternative="two-sided")

        # significance label
        if pval < 0.001:
            sig = "***"
        elif pval < 0.01:
            sig = "**"
        elif pval < 0.05:
            sig = "*"
        else:
            sig = "ns"

        bp = ax.boxplot(
            [yes_vals, no_vals],
            patch_artist=True,
            medianprops=dict(color="black", linewidth=2),
            whiskerprops=dict(linewidth=1.2),
            capprops=dict(linewidth=1.2),
            flierprops=dict(marker="o", markersize=2, alpha=0.4),
        )
        bp["boxes"][0].set_facecolor(resp_color)
        bp["boxes"][1].set_facecolor(nonresp_color)

        ax.set_title(pop.replace("_", " ").title(), fontsize=10, fontweight="bold")
        ax.set_xticks([1, 2])
        ax.set_xticklabels(["Resp\n(n={})".format(len(yes_vals)),
                             "Non-resp\n(n={})".format(len(no_vals))], fontsize=8)
        ax.set_ylabel("Frequency (%)" if pop == POPULATIONS[0] else "")
        ax.set_xlabel(f"p={pval:.4f}  {sig}", fontsize=8)

        stat_rows.append({
            "population":       pop,
            "n_responders":     len(yes_vals),
            "n_nonresponders":  len(no_vals),
            "mean_resp":        round(yes_vals.mean(), 4),
            "mean_nonresp":     round(no_vals.mean(),  4),
            "median_resp":      round(np.median(yes_vals), 4),
            "median_nonresp":   round(np.median(no_vals),  4),
            "mannwhitney_stat": round(stat, 2),
            "p_value":          round(pval, 6),
            "significant_0.05": pval < 0.05,
        })

    legend_handles = [
        mpatches.Patch(facecolor=resp_color,    label="Responders"),
        mpatches.Patch(facecolor=nonresp_color, label="Non-responders"),
    ]
    fig.legend(handles=legend_handles, loc="lower center", ncol=2,
               bbox_to_anchor=(0.5, -0.02), fontsize=9)
    plt.tight_layout(rect=[0, 0.04, 1, 1])

    plot_out = os.path.join(OUT_DIR, "part3_boxplot.png")
    fig.savefig(plot_out, dpi=150, bbox_inches="tight")
    plt.close()

    stats_df = pd.DataFrame(stat_rows)
    stats_out = os.path.join(OUT_DIR, "part3_stats.csv")
    stats_df.to_csv(stats_out, index=False)

    print(f"✓ Part 3: boxplot saved      → {plot_out}")
    print(f"✓ Part 3: stats table saved  → {stats_out}")
    print()
    print("  Statistical results (Mann-Whitney U, two-sided):")
    print(f"  {'Population':<15} {'Mean Resp':>10} {'Mean Non-R':>11} {'p-value':>9}  {'Sig?':>5}")
    print("  " + "-"*55)
    for r in stat_rows:
        sig_str = "YES *" if r["significant_0.05"] else "no"
        print(f"  {r['population']:<15} {r['mean_resp']:>10.2f}% {r['mean_nonresp']:>10.2f}%"
              f" {r['p_value']:>9.4f}  {sig_str:>5}")

    return stats_df


# ---------------------------------------------------------------------------
# Part 4 – data subset analysis
# ---------------------------------------------------------------------------

def part4_subset(conn) -> pd.DataFrame:
    """
    Melanoma PBMC samples at baseline (time=0) treated with miraclib.
    Breaks down by project, response, and sex.
    Also computes avg b_cell for male responders.
    """
    sql = """
        SELECT
            s.sample_id,
            sub.subject_id,
            sub.project_id,
            sub.condition,
            sub.sex,
            sub.response,
            sub.treatment,
            s.sample_type,
            s.time_from_treatment_start
        FROM samples s
        JOIN subjects sub ON sub.subject_id = s.subject_id
        WHERE sub.condition               = 'melanoma'
          AND s.sample_type               = 'PBMC'
          AND s.time_from_treatment_start = 0
          AND sub.treatment               = 'miraclib'
        ORDER BY s.sample_id
    """
    subset = pd.read_sql_query(sql, conn)

    # pull b_cell counts for this subset
    placeholders = ",".join("?" * len(subset))
    sample_ids   = subset["sample_id"].tolist()
    cc_sql = f"""
        SELECT sample_id, count AS b_cell_count
        FROM cell_counts
        WHERE sample_id IN ({placeholders})
          AND population = 'b_cell'
    """
    bcells = pd.read_sql_query(cc_sql, conn, params=sample_ids)
    subset = subset.merge(bcells, on="sample_id")

    out = os.path.join(OUT_DIR, "part4_subset.csv")
    subset.to_csv(out, index=False)
    print(f"✓ Part 4: subset saved        → {out}  ({len(subset):,} samples)")

    # ---- breakdown ----
    print()
    print("  Part 4 – Baseline melanoma PBMC miraclib samples")
    print(f"  Total samples: {len(subset)}")

    print()
    print("  Samples per project:")
    for proj, cnt in subset.groupby("project_id").size().items():
        print(f"    {proj}: {cnt}")

    print()
    print("  Subjects by response:")
    for resp, cnt in subset.groupby("response").size().items():
        label = "responders" if resp == "yes" else "non-responders"
        print(f"    {label}: {cnt}")

    print()
    print("  Subjects by sex:")
    for sex, cnt in subset.groupby("sex").size().items():
        label = "male" if sex == "M" else "female"
        print(f"    {label}: {cnt}")

    # ---- key numeric answer ----
    male_resp = subset[(subset["sex"] == "M") & (subset["response"] == "yes")]
    avg_bcell = male_resp["b_cell_count"].mean()
    print()
    print(f"  ★ Avg B cells – melanoma male responders at t=0: {avg_bcell:.2f}")
    print(f"    (n = {len(male_resp)} samples)")

    return subset


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def run_all():
    conn = get_conn()
    print("=" * 60)
    print("PART 2 – Frequency Table")
    print("=" * 60)
    freq_df = part2_frequency_table(conn)

    print()
    print("=" * 60)
    print("PART 3 – Statistical Analysis")
    print("=" * 60)
    part3_stats(conn, freq_df)

    print()
    print("=" * 60)
    print("PART 4 – Data Subset Analysis")
    print("=" * 60)
    part4_subset(conn)

    conn.close()

if __name__ == "__main__":
    run_all()
