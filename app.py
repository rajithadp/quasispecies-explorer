"""Quasispecies Explorer (v2).  Run with:  streamlit run app.py"""
import json

import pandas as pd
import streamlit as st

from quasigraph import (
    analyze,
    frequency_spectrum,
    mutations_vs_consensus,
    parse_haplotypes,
    read_ivar_variants,
    site_summary,
    site_table,
    to_gfa,
)
from quasigraph.plot import draw_graph, plot_site_positions, plot_site_spectrum

st.set_page_config(page_title="Quasispecies Explorer", layout="wide")
st.title("🧬 Quasispecies Explorer")
st.markdown("Upload aligned haplotypes (FASTA, with frequencies in the headers if available) and optional metadata.")

# --- SIDEBAR ---
st.sidebar.header("1. Data input")
uploaded_fasta = st.sidebar.file_uploader("Haplotype FASTA, or iVar variants table (.tsv)", type=["fasta", "fa", "fna", "tsv"])
uploaded_meta = st.sidebar.file_uploader("Metadata CSV (optional, FASTA mode only)", type=["csv"])
site_mode = uploaded_fasta is not None and uploaded_fasta.name.lower().endswith(".tsv")

st.sidebar.header("2. Parameters")
if site_mode:
    site_min_freq = st.sidebar.slider("Minimum minor-allele frequency", 0.005, 0.10, 0.02, step=0.005, format="%.3f",
                                      help="Variants near iVar's 1% threshold are partly sequencing noise; 2-3% is common practice.")
    site_min_depth = st.sidebar.number_input("Minimum read depth at a site", min_value=0, value=100, step=10)
    genome_length = st.sidebar.number_input("Genome length (bp)", min_value=1, value=29903, step=1,
                                            help="Used to turn per-site heterozygosity into a per-base diversity estimate.")
    snv_only = st.sidebar.checkbox("Single-nucleotide variants only", value=True)
else:
    rhl_threshold = st.sidebar.slider("Rare haplotype threshold", 0.001, 0.20, 0.05, step=0.001, format="%.3f")
    min_minor = st.sidebar.slider(
        "Ignore sites whose minor alleles total less than", 0.0, 0.10, 0.0, step=0.005, format="%.3f",
        help="0 = no filter. This is the noise filter that used to be hard-coded at 5%.")

if uploaded_fasta is None:
    st.info("👈 Upload a haplotype FASTA (from CliqueSNV) or an iVar variants table (.tsv) in the sidebar to begin.")
    st.stop()

if site_mode:
    try:
        variants = read_ivar_variants(uploaded_fasta.getvalue().decode("utf-8"))
        sites = site_table(variants, min_depth=int(site_min_depth), snv_only=snv_only)
        summ = site_summary(sites, int(genome_length), site_min_freq)
    except ValueError as e:
        st.error(str(e))
        st.stop()

    st.info("**Site-level mode.** Frequencies come from iVar, one genome position at a time. Linkage between sites "
            "(haplotypes) cannot be seen in this view. Variants close to iVar's 1% threshold are partly sequencing noise.")
    st.header("📊 Site-level metrics")
    c = st.columns(5)
    c[0].metric(f"Variable sites (≥{site_min_freq:.1%})", summ["variable_sites"])
    c[1].metric("Fixed differences vs reference", summ["fixed_differences"],
                help="Positions where the most common allele is not the reference base (lineage-defining, not diversity).")
    c[2].metric("Mean minor-allele frequency", f"{summ['mean_minor_freq']:.2%}")
    c[3].metric("π estimate", f"{summ['nucleotide_diversity_estimate']:.2e}")
    c[4].metric("Median depth", f"{summ['median_depth']:,}")
    if sites.empty:
        st.warning("No variant positions remain after the depth and SNV filters.")
        st.stop()

    st.divider()
    t1, t2 = st.tabs(["📍 Along the genome", "📈 Frequency spectrum"])
    with t1:
        st.pyplot(plot_site_positions(sites, int(genome_length), site_min_freq))
        variable = sites[(~sites["fixed_difference"]) & (sites["minor_freq"] >= site_min_freq)].sort_values(
            "minor_freq", ascending=False)
        st.subheader("Variable sites")
        st.dataframe(variable, use_container_width=True)
        with st.expander("Fixed differences from the reference"):
            st.dataframe(sites[sites["fixed_difference"]], use_container_width=True)
    with t2:
        spectrum = frequency_spectrum(sites)
        st.pyplot(plot_site_spectrum(spectrum))
        st.caption("All polymorphic positions after the depth and SNV filters, regardless of the frequency threshold above.")
    st.download_button("Download site table (CSV)", sites.to_csv(index=False), "site_table.csv")
    st.sidebar.download_button("Download summary (JSON)", json.dumps(summ, indent=2), "site_summary.json")
    st.stop()

try:
    hs = parse_haplotypes(uploaded_fasta.getvalue().decode("utf-8"))
    res = analyze(hs.seqs, hs.freqs, rhl_threshold=rhl_threshold, min_minor_freq=min_minor)
except ValueError as e:  # unaligned input, bad or mixed frequencies, empty file
    st.error(str(e))
    st.stop()

s = res["summary"]

# --- PROVENANCE ---
if hs.synthetic:
    st.warning(
        "These haplotypes are **synthetic**: they were sampled from allele frequencies with no linkage between sites. "
        "Haplotype-level numbers (counts, rare haplotype load, effective haplotypes) depend on how many were generated "
        "and are hidden. Per-site entropy and the variable-site table reflect the real allele frequencies."
    )
elif hs.weights_inferred:
    st.info("No frequencies found in the FASTA headers, so every record is weighted equally. "
            "If each record is one read or one clone this is fine; for haplotype callers such as CliqueSNV, "
            "use output that keeps the frequency (e.g. `>1_fr_0.58` or `freq=0.58`).")

# --- METRICS ---
st.header("📊 Population metrics")
c = st.columns(5)
if hs.synthetic:
    c[0].metric("Variable sites", s.segregating_sites)
    c[1].metric("Mean site entropy", f"{s.mean_site_entropy_bits:.4f} bits")
    c[2].metric("Nucleotide diversity (π)", f"{s.nucleotide_diversity:.5f}")
else:
    c[0].metric("Distinct haplotypes", s.n_haplotypes)
    c[1].metric("Effective haplotypes", f"{s.effective_haplotypes:.2f}")
    c[2].metric("Nucleotide diversity (π)", f"{s.nucleotide_diversity:.5f}")
    c[3].metric(f"Rare haplotype load (<{rhl_threshold:.1%})", f"{s.rhl_mass:.1%}",
                help="Share of the population carried by haplotypes below the threshold.")
    c[4].metric("Rare haplotypes", f"{s.n_rare_haplotypes} of {s.n_haplotypes}")
st.caption(f"Mean per-site entropy: {s.mean_site_entropy_bits:.4f} bits · variable sites: {s.segregating_sites} of {s.alignment_length}")

st.divider()
tab1, tab2, tab3 = st.tabs(["🔥 Entropy & hotspots", "🕸️ Variation graph", "🌡️ Metadata"])

with tab1:
    st.line_chart(res["site_entropy"].set_index("position")["entropy_bits"])
    st.subheader("Most variable sites")
    st.dataframe(res["top_sites"], use_container_width=True)
    st.download_button("Download site entropy (CSV)", res["site_entropy"].to_csv(index=False), "site_entropy.csv")

with tab2:
    vg = res["graph"]
    st.caption("Each distinct haplotype is a path from left to right. Grey = positions shared by all haplotypes, "
               "red = variable positions (bubbles). Line width and node size show frequency.")
    if vg.n_bubbles > 40:
        st.warning(f"{vg.n_bubbles} variable blocks: the drawing will be crowded. Use a shorter region or fewer haplotypes.")
    st.pyplot(draw_graph(vg))
    st.json(vg.stats())
    st.download_button("Download graph (GFA)", to_gfa(vg), "quasispecies.gfa")

with tab3:
    if uploaded_meta is None:
        st.info("Upload a metadata CSV in the sidebar. Columns: HaplotypeID, Temperature_C, Host.")
    else:
        meta_df = pd.read_csv(uploaded_meta)
        missing = {"HaplotypeID", "Temperature_C", "Host"} - set(meta_df.columns)
        if missing:
            st.error(f"Metadata is missing column(s): {', '.join(sorted(missing))}")
        else:
            dist = mutations_vs_consensus(hs.seqs, hs.freqs)
            merged = pd.merge(pd.DataFrame({"HaplotypeID": hs.names, "Mutations_vs_Consensus": dist}),
                              meta_df, on="HaplotypeID")
            if merged.empty:
                st.error("No HaplotypeID in the metadata matches the FASTA names.")
            else:
                try:
                    import plotly.express as px
                    st.plotly_chart(px.box(merged, x="Temperature_C", y="Mutations_vs_Consensus", color="Host",
                                           points="all", title="Mutations vs consensus, by temperature and host"),
                                    use_container_width=True)
                except ImportError:
                    st.caption("Install plotly to see the box plot.")
                st.dataframe(merged, use_container_width=True)

st.sidebar.download_button("Download summary (JSON)", json.dumps(s.as_dict(), indent=2), "summary.json")
