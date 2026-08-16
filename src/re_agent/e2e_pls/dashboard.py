"""Judge-facing dashboard: score a sequence's 9-mers, then "Steer to Safety".

    uv run streamlit run src/re_agent/e2e_pls/dashboard.py

All chart data comes from real calibrated-head scores and a real steering
trace -- there is no hard-coded animation or fabricated safe cluster.
Plotting/scoring helpers below are plain functions (no Streamlit calls) so
they're importable and unit-testable; the Streamlit UI lives in `render()`,
run only when this file is executed directly by `streamlit run`.
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go

from re_agent.e2e_pls import fixtures, schema
from re_agent.e2e_pls import score as score_mod
from re_agent.e2e_pls.encoder import EmbeddingCache, ProteinEncoder
from re_agent.e2e_pls.model import ThreeHeadModel
from re_agent.e2e_pls.steer import CLAIM_DISCLAIMER, SteeringConfig, SteeringTrace, steer_to_safety
from re_agent.e2e_pls.train import DEFAULT_CACHE_PATH, DEFAULT_OUTPUT_DIR

DEFAULT_HLA = "HLA-A*02:01"


def load_heads(checkpoint_dir: str | Path) -> ThreeHeadModel | None:
    try:
        return ThreeHeadModel.load(checkpoint_dir)
    except FileNotFoundError:
        return None


def fit_reference_latent_map(
    heads: ThreeHeadModel, reference_df: pd.DataFrame, client: ProteinEncoder
) -> tuple[object, np.ndarray, np.ndarray]:
    """Fixed 2D PCA fit once on the reference dataset's MHC-latent embeddings.

    Returns (fitted PCA, reference coords (N,2), reference risk proxy (N,)) --
    the "safe/risky reference zones" backdrop the demo point is plotted against.
    """
    from sklearn.decomposition import PCA

    mer_vecs = np.stack(
        [
            client.embed_pooled(r.peptide, r.n_flank or "", r.c_flank or "", "mean_9mer")
            for r in reference_df.itertuples()
        ]
    )
    z = heads.mhc.project(mer_vecs)
    pca2 = PCA(n_components=2, random_state=0).fit(z)
    coords = pca2.transform(z)
    risk_proxy = np.clip(1 - reference_df["mhc_percentile"].values / 100, 0, 1)
    return pca2, coords, risk_proxy


def steering_path_points(
    trace: SteeringTrace, heads: ThreeHeadModel, client: ProteinEncoder, pca2
) -> list[dict]:
    """Latent-map (x, y) for the target window before and after each accepted mutation."""
    points = []
    for label, seq in [("start", trace.input_sequence)] + [
        (f"step {i}", _apply_mutations(trace.input_sequence, trace.mutations[: i + 1]))
        for i in range(len(trace.mutations))
    ]:
        window = seq[trace.target_window_start : trace.target_window_end]
        n_flank = seq[
            max(
                0, trace.target_window_start - score_mod.DEFAULT_FLANK_LEN
            ) : trace.target_window_start
        ]
        c_flank = seq[
            trace.target_window_end : trace.target_window_end + score_mod.DEFAULT_FLANK_LEN
        ]
        mer_vec = client.embed_pooled(window, n_flank, c_flank, "mean_9mer")
        z = heads.mhc.project(mer_vec[None, :])
        xy = pca2.transform(z)[0]
        points.append({"label": label, "x": float(xy[0]), "y": float(xy[1])})
    return points


def _apply_mutations(sequence: str, mutations: list[dict]) -> str:
    chars = list(sequence)
    for m in mutations:
        chars[m["position"]] = m["to"]
    return "".join(chars)


def build_latent_map_figure(
    coords: np.ndarray, risk_proxy: np.ndarray, path_points: list[dict] | None = None
) -> go.Figure:
    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=coords[:, 0],
            y=coords[:, 1],
            mode="markers",
            marker=dict(
                color=risk_proxy,
                colorscale="RdYlGn_r",
                size=6,
                opacity=0.5,
                colorbar=dict(title="reference risk"),
            ),
            name="reference set (fixed PCA)",
        )
    )
    if path_points:
        fig.add_trace(
            go.Scatter(
                x=[p["x"] for p in path_points],
                y=[p["y"] for p in path_points],
                mode="markers+lines+text",
                text=[p["label"] for p in path_points],
                textposition="top center",
                marker=dict(color="blue", size=11, symbol="diamond"),
                line=dict(color="blue"),
                name="steering path",
            )
        )
    fig.update_layout(
        xaxis_title="PC1", yaxis_title="PC2", title="MHC latent map (fixed PCA on reference set)"
    )
    return fig


def build_animated_path_figure(
    coords: np.ndarray, risk_proxy: np.ndarray, path_points: list[dict]
) -> go.Figure:
    """Real recorded latent coordinates, revealed frame by frame."""
    backdrop = go.Scatter(
        x=coords[:, 0],
        y=coords[:, 1],
        mode="markers",
        marker=dict(color=risk_proxy, colorscale="RdYlGn_r", size=6, opacity=0.35),
        name="reference",
    )
    frames = []
    for i in range(1, len(path_points) + 1):
        path_trace = go.Scatter(
            x=[p["x"] for p in path_points[:i]],
            y=[p["y"] for p in path_points[:i]],
            mode="markers+lines+text",
            text=[p["label"] for p in path_points[:i]],
            textposition="top center",
            marker=dict(color="blue", size=12, symbol="diamond"),
            name="path",
        )
        frames.append(go.Frame(data=[backdrop, path_trace], name=str(i)))

    fig = go.Figure(
        data=frames[0].data if frames else [backdrop],
        frames=frames,
        layout=go.Layout(
            xaxis_title="PC1",
            yaxis_title="PC2",
            title="Original -> variant steering path",
            updatemenus=[
                dict(
                    type="buttons",
                    buttons=[
                        dict(
                            label="Play",
                            method="animate",
                            args=[None, {"frame": {"duration": 700, "redraw": True}}],
                        )
                    ],
                )
            ],
            sliders=[
                dict(
                    steps=[
                        dict(
                            method="animate",
                            args=[[f.name], {"frame": {"duration": 0, "redraw": True}}],
                            label=f.name,
                        )
                        for f in frames
                    ]
                )
            ],
        ),
    )
    return fig


def build_convergence_figure(path_rows: list[dict]) -> go.Figure:
    accepted = [r for r in path_rows if r["accepted"]]
    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=[r["step"] for r in accepted],
            y=[r["composite_risk"] for r in accepted],
            name="composite risk",
            mode="lines+markers",
        )
    )
    fig.add_trace(
        go.Scatter(
            x=[r["step"] for r in accepted],
            y=[r["latent_distance"] for r in accepted],
            name="latent distance",
            mode="lines+markers",
            yaxis="y2",
        )
    )
    fig.update_layout(
        xaxis_title="step",
        yaxis=dict(title="composite risk"),
        yaxis2=dict(title="latent distance", overlaying="y", side="right"),
        title="Steering convergence (accepted steps)",
    )
    return fig


def sequence_diff(original: str, variant: str) -> list[dict]:
    return [
        {"position": i, "original": o, "variant": v}
        for i, (o, v) in enumerate(zip(original, variant, strict=True))
        if o != v
    ]


def benchmark_surrogate(
    sequence: str,
    hla_allele: str,
    client: ProteinEncoder,
    heads: ThreeHeadModel,
    teacher_seconds_per_window: float | None = None,
) -> dict:
    """Measured surrogate throughput; only computes a speedup ratio if the
    caller supplies an actual measured teacher-per-window time -- never fabricated.
    """
    t0 = time.monotonic()
    scores = score_mod.score_sequence(sequence, hla_allele, client, heads)
    elapsed = time.monotonic() - t0
    result = {
        "n_windows": len(scores),
        "surrogate_seconds": elapsed,
        "surrogate_windows_per_second": (len(scores) / elapsed) if elapsed > 0 else None,
    }
    if teacher_seconds_per_window is not None:
        teacher_total = teacher_seconds_per_window * len(scores)
        result["teacher_seconds_observed"] = teacher_total
        result["speedup_ratio"] = (teacher_total / elapsed) if elapsed > 0 else None
    return result


def render() -> None:
    import streamlit as st

    st.set_page_config(page_title="E2E-PLS: MHC-I processing surrogate", layout="wide")
    st.title("E2E-PLS — MHC-I processing surrogate + Steer to Safety")
    st.caption(CLAIM_DISCLAIMER)

    with st.sidebar:
        st.header("Runtime")
        checkpoint_dir = st.text_input("Checkpoint dir", str(DEFAULT_OUTPUT_DIR))
        encoder_mode = st.selectbox("Encoder mode", ["esm2", "mock", "modal"], index=0)
        st.caption(
            "esm2 = real ESM-2 (t33, 650M) loaded locally, ~2.5 GB download on first use. "
            "mock = deterministic offline stand-in for tests. "
            "modal = deployed ESM3-open, requires the user's own HF_TOKEN + Modal auth."
        )

    heads = load_heads(checkpoint_dir)
    if heads is None:
        st.error(
            f"No checkpoint at `{checkpoint_dir}`. "
            "Train one first: `uv run python -m re_agent.e2e_pls.train`"
        )
        st.stop()

    cache = EmbeddingCache(DEFAULT_CACHE_PATH)
    client = ProteinEncoder(mode=encoder_mode, cache=cache)

    st.subheader("Checkpoint")
    cols = st.columns(4)
    cols[0].metric("Model version", heads.model_version)
    cols[1].metric("Encoder", heads.encoder_model_id)
    cols[2].metric("Dataset hash", heads.dataset_version_hash or "n/a")
    cols[3].metric("HLA alleles", ", ".join(heads.mhc.centroids) or "none")

    st.subheader("Input sequence")
    input_mode = st.radio("Source", ["Demo seed", "Paste sequence"], horizontal=True)
    if input_mode == "Demo seed":
        demo_key = st.selectbox(
            "Demo seed",
            list(fixtures.DEMO_SEQUENCES),
            format_func=lambda k: f"{k} ({fixtures.DEMO_SEQUENCES[k]['kind']})",
        )
        demo = fixtures.DEMO_SEQUENCES[demo_key]
        sequence = demo["sequence"]
        st.caption(demo["description"])
    else:
        sequence = st.text_area("Sequence (single-letter amino acids)", "").strip().upper()

    available_alleles = list(heads.mhc.centroids) or [DEFAULT_HLA]
    hla_allele = st.selectbox("HLA allele", available_alleles)

    if not sequence or len(sequence) < 9:
        st.info("Enter or select a sequence of at least 9 residues.")
        st.stop()
    if not set(sequence) <= schema.CANONICAL_RESIDUES:
        st.error("Sequence contains non-canonical residues.")
        st.stop()

    bench = benchmark_surrogate(sequence, hla_allele, client, heads)
    scores = score_mod.score_sequence(sequence, hla_allele, client, heads)
    protein_risk = score_mod.aggregate_protein_risk(scores)

    scoring_ms = bench["surrogate_seconds"] * 1000
    st.subheader(f"Checkpoint scores — {bench['n_windows']} windows in {scoring_ms:.0f} ms")
    score_df = pd.DataFrame([s.to_dict() for s in scores])
    max_idx = score_df["composite_risk"].idxmax()
    st.dataframe(
        score_df.style.apply(
            lambda row: ["background-color: #ffd6d6" if row.name == max_idx else "" for _ in row],
            axis=1,
        ),
        width="stretch",
    )
    w = protein_risk.max_risk_window
    st.metric(
        "Highest-risk window",
        f"{w.peptide} @ {w.start}-{w.end}",
        f"risk={protein_risk.max_risk:.3f}",
    )

    st.subheader("Latent map")
    reference_df = fixtures.load_dev_fixture()
    pca2, coords, risk_proxy = fit_reference_latent_map(heads, reference_df, client)
    st.plotly_chart(build_latent_map_figure(coords, risk_proxy), width="stretch")

    st.subheader("Steer to Safety")
    with st.form("steer_form"):
        max_mutations = st.slider("Max mutations", 1, 3, 3)
        preservation_threshold = st.slider(
            "Preservation threshold (max per-residue log-lik drop)", 0.1, 2.0, 0.5
        )
        protected_str = st.text_input(
            "Protected positions (comma-separated, 0-indexed in the full sequence)", ""
        )
        submitted = st.form_submit_button("Steer to Safety")

    if submitted:
        protected = frozenset(int(p) for p in protected_str.split(",") if p.strip())
        config = SteeringConfig(
            hla_allele=hla_allele,
            max_mutations=max_mutations,
            preservation_threshold=preservation_threshold,
            protected_positions=protected,
        )
        with st.spinner("Steering..."):
            trace = steer_to_safety(sequence, w.start, w.end, client, heads, config)
        st.session_state["trace"] = trace

    trace: SteeringTrace | None = st.session_state.get("trace")
    if trace is None:
        return

    st.markdown(f"**Mutations applied:** {len(trace.mutations)} — {trace.mutations}")
    before, after = st.columns(2)
    before.metric("Composite risk before", f"{trace.initial_score['composite_risk']:.3f}")
    after.metric(
        "Composite risk after",
        f"{trace.final_score['composite_risk']:.3f}",
        f"{trace.final_score['composite_risk'] - trace.initial_score['composite_risk']:.3f}",
    )

    path_points = steering_path_points(trace, heads, client, pca2)
    st.plotly_chart(build_animated_path_figure(coords, risk_proxy, path_points), width="stretch")
    st.plotly_chart(build_convergence_figure(trace.optimization_path()), width="stretch")

    st.subheader("Sequence diff")
    diff = sequence_diff(trace.input_sequence, trace.output_sequence)
    st.dataframe(pd.DataFrame(diff), width="stretch")

    st.subheader("Preservation checks")
    protected_ok = all(
        m["position"] not in trace.constraints["protected_positions"] for m in trace.mutations
    )
    st.write(
        {
            "mutations_within_budget": len(trace.mutations) <= trace.constraints["max_mutations"],
            "protected_positions_untouched": protected_ok,
            "preservation_threshold": trace.constraints["preservation_threshold"],
        }
    )

    st.download_button(
        "Download trace (JSON)",
        data=json.dumps(trace.to_dict(), indent=2),
        file_name="steering_trace.json",
        mime="application/json",
    )


if __name__ == "__main__":
    render()
