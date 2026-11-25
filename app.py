import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from Bio import SeqIO
from scipy.stats import entropy
import io

st.set_page_config(page_title="Quasispecies Explorer", layout="wide")

st.title("🧬 Quasispecies Explorer")
st.markdown("Upload aligned FASTA haplotypes and optional metadata to quantify and visualize viral diversity.")

# --- SIDEBAR ---
st.sidebar.header("1. Data Input")
uploaded_fasta = st.sidebar.file_uploader("Upload Haplotype FASTA", type=["fasta", "fa", "fna"])
uploaded_meta = st.sidebar.file_uploader("Upload Metadata (CSV) - Optional", type=["csv"])

if uploaded_fasta is not None:
    # 1. Parse FASTA
    fasta_io = io.StringIO(uploaded_fasta.getvalue().decode("utf-8"))
    records = list(SeqIO.parse(fasta_io, "fasta"))
    sequences = [str(rec.seq) for rec in records]
    hap_ids = [rec.id for rec in records]
    seq_length = len(sequences[0])
    
    # 2. Parse Metadata (if provided)
    meta_df = None
    if uploaded_meta is not None:
        meta_df = pd.read_csv(uploaded_meta)

    # 3. Calculate Metrics
    num_haplotypes = len(set(sequences))
    total_reads = len(sequences)
    
    # Haplotype Frequencies
    hap_counts = pd.Series(sequences).value_counts()
    hap_freqs = hap_counts / total_reads
    
    # Rare Haplotype Load (RHL) - Fraction of haplotypes at < 5% frequency
    rht_threshold = 0.05
    rare_haplotypes = hap_freqs[hap_freqs < rht_threshold]
    rhl = rare_haplotypes.sum()
    
    # Per-site Entropy & Nucleotide Diversity
    entropies, pi_values = [], []
    for i in range(seq_length):
        column = [seq[i] for seq in sequences]
        counts = pd.Series(column).value_counts(normalize=True)
        entropies.append(entropy(counts, base=2))
        
        sum_sq = sum(p**2 for p in counts)
        pi = (total_reads / (total_reads - 1)) * (1 - sum_sq) if total_reads > 1 else 0
        pi_values.append(pi)

    mean_entropy = np.mean(entropies)
    mean_pi = np.mean(pi_values)

    # --- DISPLAY METRICS ---
    st.header("📊 Population Metrics")
    c1, c2, c3, c4, c5 = st.columns(5)
    c1.metric("Total Haplotypes", num_haplotypes)
    c2.metric("Total Sequences", total_reads)
    c3.metric("Mean Shannon Entropy", f"{mean_entropy:.3f} bits")
    c4.metric("Nucleotide Diversity (π)", f"{mean_pi:.5f}")
    c5.metric(f"Rare Haplotype Load (<{rht_threshold*100}%)", f"{rhl:.2%}")

    st.divider()

    # --- VISUALIZATIONS ---
    tab1, tab2, tab3 = st.tabs(["🔥 Entropy & Hotspots", "🕸️ Variation Graph", "🌡️ Metadata Correlation"])
    
    with tab1:
        df_entropy = pd.DataFrame({"Position": range(1, seq_length + 1), "Entropy": entropies})
        fig_entropy = px.line(df_entropy, x="Position", y="Entropy", 
                              title="Shannon Entropy Across Genome",
                              labels={"Entropy": "Entropy (bits)", "Position": "Genomic Position (bp)"})
        fig_entropy.update_traces(line_color='#e74c3c', fill='tozeroy')
        st.plotly_chart(fig_entropy, use_container_width=True)

    with tab2:
        st.subheader("Haplotype Paths (Bubbles = Variants)")
        fig_graph = go.Figure()
        
        # We will draw the backbone
        fig_graph.add_trace(go.Scatter(x=[1, seq_length], y=[0, 0], mode='lines', 
                                       name='Consensus Backbone', line=dict(color='black', width=3)))
        
        # Draw variant paths (colored by frequency)
        colors = px.colors.sequential.Viridis
        variant_positions = [i for i, e in enumerate(entropies) if e > 0.1] # Only major hotspots
        
        for i, pos in enumerate(variant_positions):
            # Find the variant allele
            column = [seq[pos] for seq in sequences]
            ref_allele = column[0]
            var_allele = pd.Series(column).value_counts().index[1]
            
            # Get frequency
            freq = pd.Series(column).value_counts(normalize=True)[var_allele]
            color = colors[int(freq * 10) if int(freq*10) < len(colors) else -1]
            
            # Draw the bubble
            fig_graph.add_trace(go.Scatter(
                x=[pos, pos, pos], 
                y=[0, 1, 0], 
                mode='lines+markers', 
                name=f"{ref_allele}->{var_allele} at {pos} ({freq:.1%})",
                line=dict(color=color, width=2),
                marker=dict(size=10, color=color),
                hoverinfo='name'
            ))
            
        fig_graph.update_layout(title="Variation Graph (Colored by Frequency)",
                                xaxis_title="Genomic Position", 
                                yaxis_title="Bubble Level",
                                height=400)
        st.plotly_chart(fig_graph, use_container_width=True)

    with tab3:
        st.subheader("Quasispecies vs. Environmental Metadata")
        if meta_df is not None:
            # Merge metrics with metadata
            # We need per-haplotype diversity, so we'll calculate a simple proxy: sequence divergence from consensus
            consensus = "".join([pd.Series([seq[i] for seq in sequences]).mode()[0] for i in range(seq_length)])
            
            hap_metrics = []
            for seq, hap_id in zip(sequences, hap_ids):
                diff = sum(1 for a, b in zip(seq, consensus) if a != b)
                hap_metrics.append({"HaplotypeID": hap_id, "Mutations_vs_Consensus": diff})
            
            df_hap_metrics = pd.DataFrame(hap_metrics)
            merged_df = pd.merge(df_hap_metrics, meta_df, on="HaplotypeID")
            
            # Plot: Temperature vs Mutations
            fig_scatter = px.box(merged_df, x="Temperature_C", y="Mutations_vs_Consensus", 
                                 color="Host", points="all",
                                 title="Mutations per Haplotype vs. Temperature & Host",
                                 labels={"Temperature_C": "Temperature (°C)", "Mutations_vs_Consensus": "Mutations vs Consensus"})
            st.plotly_chart(fig_scatter, use_container_width=True)
            st.dataframe(merged_df.head(10), use_container_width=True)
        else:
            st.info("Upload a metadata CSV in the sidebar to see correlations.")
