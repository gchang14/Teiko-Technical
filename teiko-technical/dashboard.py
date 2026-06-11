"""
dashboard.py  –  Interactive Dash dashboard for all four parts
"""

import sqlite3, os
import pandas as pd
import numpy as np
from scipy import stats as scipy_stats

import dash
from dash import dcc, html, dash_table, Input, Output, callback
import plotly.graph_objects as go
import plotly.express as px

DB_PATH = "teiko.db"
POPULATIONS = ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]
POP_LABELS  = {
    "b_cell":      "B Cell",
    "cd8_t_cell":  "CD8 T Cell",
    "cd4_t_cell":  "CD4 T Cell",
    "nk_cell":     "NK Cell",
    "monocyte":    "Monocyte",
}

# ── colours ────────────────────────────────────────────────────────────────
RESP_COLOR    = "#4C9BE8"
NONRESP_COLOR = "#E8724C"
BG            = "#F7F9FC"
CARD_BG       = "#FFFFFF"
ACCENT        = "#1A3C6E"

# ── load data once ─────────────────────────────────────────────────────────

def load_freq_df():
    conn = sqlite3.connect(DB_PATH)
    sql = """
        SELECT
            s.sample_id   AS sample,
            sub.subject_id,
            sub.project_id,
            sub.condition,
            sub.sex,
            sub.treatment,
            sub.response,
            s.sample_type,
            s.time_from_treatment_start,
            cc.population,
            cc.count
        FROM cell_counts cc
        JOIN samples  s   ON s.sample_id  = cc.sample_id
        JOIN subjects sub ON sub.subject_id = s.subject_id
    """
    df = pd.read_sql_query(sql, conn)
    conn.close()
    totals = df.groupby("sample")["count"].sum().rename("total_count")
    df = df.join(totals, on="sample")
    df["percentage"] = (df["count"] / df["total_count"] * 100).round(4)
    return df

FREQ_DF = load_freq_df()

# ── app ─────────────────────────────────────────────────────────────────────

app = dash.Dash(__name__, title="Teiko – Immune Profiling Dashboard")

def card(children, style=None):
    base = {
        "background": CARD_BG,
        "borderRadius": "10px",
        "boxShadow": "0 2px 8px rgba(0,0,0,0.08)",
        "padding": "20px 24px",
        "marginBottom": "20px",
    }
    if style:
        base.update(style)
    return html.Div(children, style=base)

def section_title(text):
    return html.H3(text, style={
        "color": ACCENT, "marginTop": 0, "marginBottom": "14px",
        "fontSize": "16px", "fontWeight": "700", "letterSpacing": "0.3px",
    })

app.layout = html.Div(style={"background": BG, "minHeight": "100vh",
                              "fontFamily": "'Segoe UI', Arial, sans-serif",
                              "padding": "0"}, children=[

    # ── header ──────────────────────────────────────────────────────────────
    html.Div(style={
        "background": ACCENT, "color": "white",
        "padding": "18px 36px", "marginBottom": "28px",
        "display": "flex", "alignItems": "center", "gap": "16px",
    }, children=[
        html.H1("🔬 Teiko Clinical Immune Profiling",
                style={"margin": 0, "fontSize": "22px", "fontWeight": "700"}),
        html.Span("Interactive analysis dashboard",
                  style={"opacity": "0.75", "fontSize": "13px"}),
    ]),

    html.Div(style={"maxWidth": "1400px", "margin": "0 auto", "padding": "0 24px"}, children=[

        # ── tabs ────────────────────────────────────────────────────────────
        dcc.Tabs(id="tabs", value="tab-2", style={"marginBottom": "20px"},
                 colors={"border": "#ddd", "primary": ACCENT, "background": BG},
        children=[
            dcc.Tab(label="Part 2 – Frequency Table", value="tab-2"),
            dcc.Tab(label="Part 3 – Responder Analysis", value="tab-3"),
            dcc.Tab(label="Part 4 – Subset Analysis", value="tab-4"),
        ]),

        html.Div(id="tab-content"),
    ]),
])


# ── Tab layouts ─────────────────────────────────────────────────────────────

def layout_part2():
    conditions  = sorted(FREQ_DF["condition"].unique())
    projects    = sorted(FREQ_DF["project_id"].unique())
    sample_types= sorted(FREQ_DF["sample_type"].unique())

    return html.Div([
        card([
            section_title("Filters"),
            html.Div(style={"display": "flex", "gap": "24px", "flexWrap": "wrap"}, children=[
                html.Div([
                    html.Label("Condition", style={"fontWeight": "600", "fontSize": "13px"}),
                    dcc.Dropdown(id="p2-condition",
                                 options=[{"label": "All", "value": "All"}] +
                                         [{"label": c, "value": c} for c in conditions],
                                 value="All", clearable=False, style={"width": "160px"}),
                ]),
                html.Div([
                    html.Label("Project", style={"fontWeight": "600", "fontSize": "13px"}),
                    dcc.Dropdown(id="p2-project",
                                 options=[{"label": "All", "value": "All"}] +
                                         [{"label": p, "value": p} for p in projects],
                                 value="All", clearable=False, style={"width": "160px"}),
                ]),
                html.Div([
                    html.Label("Sample Type", style={"fontWeight": "600", "fontSize": "13px"}),
                    dcc.Dropdown(id="p2-stype",
                                 options=[{"label": "All", "value": "All"}] +
                                         [{"label": s, "value": s} for s in sample_types],
                                 value="All", clearable=False, style={"width": "160px"}),
                ]),
            ]),
        ]),

        card([
            section_title("Relative Frequency by Population"),
            dcc.Graph(id="p2-bar", style={"height": "420px"}),
        ]),

        card([
            section_title("Frequency Table (sample × population)"),
            html.Div(id="p2-table"),
        ]),
    ])


def layout_part3():
    return html.Div([
        card([
            section_title("Analysis: Melanoma · PBMC · Miraclib – Responders vs Non-Responders"),
            html.P("Comparing relative cell-population frequencies between patients who responded "
                   "to miraclib and those who did not.  Statistical test: Mann-Whitney U (two-sided).",
                   style={"color": "#555", "fontSize": "13px", "margin": "0 0 16px 0"}),
            dcc.Graph(id="p3-boxplot", style={"height": "500px"}),
        ]),
        card([
            section_title("Statistical Summary"),
            html.Div(id="p3-stats-table"),
        ]),
    ])


def layout_part4():
    return html.Div([
        card([
            section_title("Baseline Subset: Melanoma · PBMC · t = 0 · Miraclib"),
            html.P("Filters applied: condition = melanoma, sample_type = PBMC, "
                   "time_from_treatment_start = 0, treatment = miraclib.",
                   style={"color": "#555", "fontSize": "13px", "margin": "0 0 20px 0"}),
            html.Div(style={"display": "flex", "gap": "16px", "flexWrap": "wrap"}, children=[
                html.Div(id="p4-summary-cards"),
            ]),
        ]),
        html.Div(style={"display": "grid", "gridTemplateColumns": "1fr 1fr 1fr",
                        "gap": "16px"}, children=[
            card([section_title("Samples per Project"), dcc.Graph(id="p4-project-pie", style={"height": "280px"})]),
            card([section_title("Responders vs Non-Responders"), dcc.Graph(id="p4-response-pie", style={"height": "280px"})]),
            card([section_title("Sex Breakdown"), dcc.Graph(id="p4-sex-pie", style={"height": "280px"})]),
        ]),
        card([
            section_title("B Cell Counts at Baseline – Responders vs Non-Responders"),
            dcc.Graph(id="p4-bcell-box", style={"height": "380px"}),
            html.Div(id="p4-bcell-answer", style={
                "background": "#EAF3FB", "borderRadius": "8px",
                "padding": "14px 18px", "marginTop": "14px",
                "fontWeight": "600", "color": ACCENT, "fontSize": "15px",
            }),
        ]),
    ])


# ── Tab routing ─────────────────────────────────────────────────────────────

@app.callback(Output("tab-content", "children"), Input("tabs", "value"))
def render_tab(tab):
    if tab == "tab-2":
        return layout_part2()
    if tab == "tab-3":
        return layout_part3()
    return layout_part4()


# ── Part 2 callbacks ────────────────────────────────────────────────────────

@app.callback(
    Output("p2-bar",   "figure"),
    Output("p2-table", "children"),
    Input("p2-condition", "value"),
    Input("p2-project",   "value"),
    Input("p2-stype",     "value"),
)
def update_part2(condition, project, stype):
    df = FREQ_DF.copy()
    if condition != "All": df = df[df["condition"]   == condition]
    if project   != "All": df = df[df["project_id"]  == project]
    if stype     != "All": df = df[df["sample_type"] == stype]

    # mean % per population across filtered samples
    agg = df.groupby("population")["percentage"].mean().reset_index()
    agg["pop_label"] = agg["population"].map(POP_LABELS)
    agg = agg.sort_values("percentage", ascending=False)

    fig = go.Figure(go.Bar(
        x=agg["pop_label"], y=agg["percentage"],
        marker_color=ACCENT, text=agg["percentage"].round(2).astype(str) + "%",
        textposition="outside",
    ))
    fig.update_layout(
        yaxis_title="Mean Relative Frequency (%)",
        xaxis_title="Cell Population",
        plot_bgcolor="white", paper_bgcolor="white",
        margin=dict(t=20, b=40),
        yaxis=dict(gridcolor="#eee"),
    )

    # table (limit to 200 rows for performance)
    tbl_df = df[["sample","total_count","population","count","percentage"]].head(200)
    tbl_df["population"] = tbl_df["population"].map(POP_LABELS)
    tbl_df = tbl_df.rename(columns={
        "sample": "Sample", "total_count": "Total Count",
        "population": "Population", "count": "Count", "percentage": "% Frequency",
    })
    table = dash_table.DataTable(
        data=tbl_df.to_dict("records"),
        columns=[{"name": c, "id": c} for c in tbl_df.columns],
        page_size=15,
        style_table={"overflowX": "auto"},
        style_cell={"fontSize": "13px", "padding": "6px 10px", "fontFamily": "inherit"},
        style_header={"fontWeight": "700", "backgroundColor": "#f0f4fa",
                      "color": ACCENT, "border": "1px solid #ddd"},
        style_data_conditional=[
            {"if": {"row_index": "odd"}, "backgroundColor": "#fafbfd"}
        ],
    )
    return fig, table


# ── Part 3 callbacks ────────────────────────────────────────────────────────

@app.callback(
    Output("p3-boxplot",    "figure"),
    Output("p3-stats-table","children"),
    Input("tabs", "value"),
)
def update_part3(tab):
    if tab != "tab-3":
        return go.Figure(), html.Div()

    df = FREQ_DF[
        (FREQ_DF["condition"]   == "melanoma") &
        (FREQ_DF["treatment"]   == "miraclib") &
        (FREQ_DF["response"].isin(["yes","no"])) &
        (FREQ_DF["sample_type"] == "PBMC")
    ].copy()

    fig = go.Figure()
    stat_rows = []
    positions = {p: i+1 for i, p in enumerate(POPULATIONS)}

    for pop in POPULATIONS:
        pop_data  = df[df["population"] == pop]
        yes_vals  = pop_data[pop_data["response"] == "yes"]["percentage"].values
        no_vals   = pop_data[pop_data["response"] == "no" ]["percentage"].values
        stat, pval = scipy_stats.mannwhitneyu(yes_vals, no_vals, alternative="two-sided")

        label = POP_LABELS[pop]
        fig.add_trace(go.Box(
            y=yes_vals, name=f"{label}<br>Resp", x0=label,
            offsetgroup="Responders",
            marker_color=RESP_COLOR, line_color=RESP_COLOR,
            boxmean=False, showlegend=(pop == POPULATIONS[0]),
            legendgroup="resp", legendgrouptitle_text="",
        ))
        fig.add_trace(go.Box(
            y=no_vals, name=f"{label}<br>Non-Resp", x0=label,
            offsetgroup="Non-Responders",
            marker_color=NONRESP_COLOR, line_color=NONRESP_COLOR,
            boxmean=False, showlegend=(pop == POPULATIONS[0]),
            legendgroup="nonresp",
        ))

        sig = ("***" if pval < 0.001 else "**" if pval < 0.01
               else "*" if pval < 0.05 else "ns")
        stat_rows.append({
            "Population":      POP_LABELS[pop],
            "N Resp":          len(yes_vals),
            "N Non-Resp":      len(no_vals),
            "Mean Resp (%)":   round(float(yes_vals.mean()), 2),
            "Mean Non-Resp (%)": round(float(no_vals.mean()), 2),
            "p-value":         round(pval, 6),
            "Significant (α=0.05)": "✓ YES" if pval < 0.05 else "—",
            "Annotation":      sig,
        })

    fig.update_layout(
        boxmode="group",
        yaxis_title="Relative Frequency (%)",
        xaxis_title="Cell Population",
        plot_bgcolor="white", paper_bgcolor="white",
        legend=dict(title="Group"),
        margin=dict(t=30, b=40),
        yaxis=dict(gridcolor="#eee"),
    )

    stats_df = pd.DataFrame(stat_rows)
    table = dash_table.DataTable(
        data=stats_df.to_dict("records"),
        columns=[{"name": c, "id": c} for c in stats_df.columns],
        style_table={"overflowX": "auto"},
        style_cell={"fontSize": "13px", "padding": "6px 12px", "fontFamily": "inherit"},
        style_header={"fontWeight": "700", "backgroundColor": "#f0f4fa",
                      "color": ACCENT, "border": "1px solid #ddd"},
        style_data_conditional=[
            {"if": {"filter_query": '{Significant (α=0.05)} = "✓ YES"'},
             "backgroundColor": "#EAF6EA", "color": "#1A6E2E", "fontWeight": "600"},
            {"if": {"row_index": "odd"}, "backgroundColor": "#fafbfd"},
        ],
    )
    return fig, table


# ── Part 4 callbacks ────────────────────────────────────────────────────────

@app.callback(
    Output("p4-summary-cards", "children"),
    Output("p4-project-pie",   "figure"),
    Output("p4-response-pie",  "figure"),
    Output("p4-sex-pie",       "figure"),
    Output("p4-bcell-box",     "figure"),
    Output("p4-bcell-answer",  "children"),
    Input("tabs", "value"),
)
def update_part4(tab):
    if tab != "tab-4":
        empty = go.Figure()
        return html.Div(), empty, empty, empty, empty, ""

    df = FREQ_DF[
        (FREQ_DF["condition"]               == "melanoma") &
        (FREQ_DF["sample_type"]             == "PBMC") &
        (FREQ_DF["time_from_treatment_start"]== 0) &
        (FREQ_DF["treatment"]               == "miraclib") &
        (FREQ_DF["population"]              == "b_cell")
    ].copy()

    total = len(df)

    def metric_card(label, value):
        return html.Div(style={
            "background": "#EAF3FB", "borderRadius": "8px",
            "padding": "14px 20px", "minWidth": "160px",
        }, children=[
            html.Div(str(value), style={"fontSize": "28px", "fontWeight": "800", "color": ACCENT}),
            html.Div(label, style={"fontSize": "12px", "color": "#555", "marginTop": "2px"}),
        ])

    summary = html.Div(style={"display": "flex", "gap": "12px", "flexWrap": "wrap"}, children=[
        metric_card("Total Samples",       total),
        metric_card("Projects",            df["project_id"].nunique()),
        metric_card("Responders",          (df["response"] == "yes").sum()),
        metric_card("Non-Responders",      (df["response"] == "no").sum()),
        metric_card("Male",                (df["sex"] == "M").sum()),
        metric_card("Female",              (df["sex"] == "F").sum()),
    ])

    pie_colors = [ACCENT, RESP_COLOR, NONRESP_COLOR, "#6BBF59", "#F5A623"]

    def pie(series, title):
        vc = series.value_counts()
        fig = go.Figure(go.Pie(
            labels=vc.index, values=vc.values,
            marker_colors=pie_colors, hole=0.4,
            textinfo="label+percent",
        ))
        fig.update_layout(
            showlegend=False, margin=dict(t=10, b=10, l=10, r=10),
            paper_bgcolor="white",
        )
        return fig

    proj_fig = pie(df["project_id"], "Project")
    resp_fig = pie(df["response"].map({"yes": "Responders", "no": "Non-Resp"}), "Response")
    sex_fig  = pie(df["sex"].map({"M": "Male", "F": "Female"}), "Sex")

    # B cell boxplot by response
    bcell_fig = go.Figure()
    for resp, color, label in [("yes", RESP_COLOR, "Responders"), ("no", NONRESP_COLOR, "Non-Responders")]:
        vals = df[df["response"] == resp]["count"].values
        bcell_fig.add_trace(go.Box(
            y=vals, name=label, marker_color=color, line_color=color,
            boxpoints="outliers",
        ))
    bcell_fig.update_layout(
        yaxis_title="B Cell Count", plot_bgcolor="white", paper_bgcolor="white",
        margin=dict(t=20, b=20), yaxis=dict(gridcolor="#eee"),
    )

    # key answer
    male_resp = df[(df["sex"] == "M") & (df["response"] == "yes")]
    avg = male_resp["count"].mean()
    answer = f"★  Average B cell count – Melanoma male responders at t = 0:  {avg:.2f}  (n = {len(male_resp)} samples)"

    return summary, proj_fig, resp_fig, sex_fig, bcell_fig, answer


if __name__ == "__main__":
    app.run(debug=False, host="0.0.0.0", port=8050)
