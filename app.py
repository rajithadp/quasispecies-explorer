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
st.markdown("Upload an aligned FASTA file of viral haplotypes to calculate diversity metrics and visualize the variation graph.")

# --- SIDEBAR ---
st.sidebar.header("Data Input")
uploaded_file = st.sidebar.file_uploader("Upload FASTA file", type=["fasta", "fa", "fna"])

if uploaded_file is not None:
    # 1. Parse the FASTA
    fasta_io = io.StringIO(uploaded_file.getvalue().decode("utf-8"))
    records = list(SeqIO.parse(fasta_io, "fasta"))
    
    if len(records) == 0:
        st.error("No sequences found in the FASTA file.")
        st.stop()
        
    sequences = [str(rec.seq) for rec in records]
    seq_length = len(sequences[0])
    
    # Ensure all sequences are aligned (same length)
    if not all(len(s) == seq_length for s in sequences):
        st.error("Sequences are not aligned (they have different lengths). Please provide an MSA.")
        st.stop()

    # 2. Calculate Metrics
    num_haplotypes = len(set(sequences))
    total_reads = len(sequences)
    
    # Per-site Shannon Entropy & Nucleotide Diversity
    entropies = []
    pi_values = []
    
    for i in range(seq_length):
        column = [seq[i] for seq in sequences]
        # Entropy
        counts = pd.Series(column).value_counts(normalize=True)
        entropies.append(entropy(counts, base=2))
        
        # Nucleotide Diversity (Pi) at this site
        # Pi = (n / (n-1)) * (1 - sum(p_i^2))
        if total_reads > 1:
            sum_sq = sum(p**2 for p in counts)
            pi = (total_reads / (total_reads - 1)) * (1 - sum_sq)
            pi_values.append(pi)
        else:
            pi_values.append(0)

    mean_entropy = np.mean(entropies)
    mean_pi = np.mean(pi_values)

    # --- DISPLAY METRICS ---
    st.header("📊 Population Metrics")
    col1, col2, col3, col4 = st.columns(4)
    col1.metric("Total Haplotypes", num_haplotypes)
    col2.metric("Total Sequences", total_reads)
    col3.metric("Mean Shannon Entropy", f"{mean_entropy:.3f} bits")
    col4.metric("Nucleotide Diversity (π)", f"{mean_pi:.5f}")

    st.divider()

    # --- VISUALIZATIONS ---
    col_left, col_right = st.columns(2)

    with col_left:
        st.subheader("🔥 Mutational Hotspots (Entropy)")
        df_entropy = pd.DataFrame({"Position": range(1, seq_length + 1), "Entropy": entropies})
        fig_entropy = px.line(df_entropy, x="Position", y="Entropy", 
                              title="Shannon Entropy Across Genome",
                              labels={"Entropy": "Entropy (bits)", "Position": "Genomic Position (bp)"})
        fig_entropy.update_traces(line_color='#e74c3c', fill='tozeroy')
        st.plotly_chart(fig_entropy, use_container_width=True)

    with col_right:
        st.subheader("🕸️ Variation Graph")
        # A simplified variation graph builder
        # We look at columns with entropy > 0, these are our "bubbles"
        variant_positions = [i for i, e in enumerate(entropies) if e > 0]
        
        fig_graph = go.Figure()
        
        # Draw paths for each unique haplotype
        unique_seqs = list(set(sequences))
        colors = px.colors.qualitative.Plotly
        
        for idx, seq in enumerate(unique_seqs):
            color = colors[idx % len(colors)]
            # Draw the "backbone" (y=0 for reference base, y=1 for variant base)
            # To make it look like a graph, we draw horizontal lines and vertical jumps
            x_vals = []
            y_vals = []
            
            ref_seq = sequences[0]
            
            for i in range(seq_length):
                x_vals.append(i + 1)
                if seq[i] == ref_seq[i]:
                    y_vals.append(0 + (idx * 0.1)) # Baseline for backbone
                else:
                    y_vals.append(1 + (idx * 0.1)) # Branch up for variant
                    
            fig_graph.add_trace(go.Scatter(x=x_vals, y=y_vals, mode='lines', 
                                           name=f"Haplotype {idx+1}",
                                           line=dict(color=color, width=2, shape='hv')))
            
        fig_graph.update_layout(title="Haplotype Paths (Bubbles = Variants)",
                                xaxis_title="Genomic Position", 
                                yaxis_title="Graph Node Level",
                                showlegend=True,
                                height=400)
        st.plotly_chart(fig_graph, use_container_width=True)

else:
    st.info("👈 Please upload a FASTA file in the sidebar to begin. You can use the `simulated_haplotypes.fasta` generated in the previous step.")
