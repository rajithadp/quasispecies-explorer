# 🧬 Quasispecies Explorer

An end-to-end bioinformatics pipeline and interactive web application for quantifying and visualizing intra-host viral quasispecies from raw FASTQ reads.

This tool addresses the challenge of viral diversity by moving beyond consensus sequences. It maps raw reads, trims amplicon primers, calls low-frequency variants, reconstructs haplotypes, and provides an intuitive web dashboard to calculate diversity metrics and correlate them with environmental metadata.

## 🚀 Features

*   **Raw Data Processing:** Automated mapping (`minimap2`), ARTIC primer trimming (`iVar`), and variant calling (`iVar variants`).
*   **Smart Haplotype Reconstruction:** Utilizes `CliqueSNV` to find physically linked mutations. If the sample is too diverse for strict linkage analysis, it automatically falls back to a frequency-based reconstruction algorithm to simulate the quasispecies cloud.
*   **Population Metrics:** Calculates Total Haplotypes, Mean Shannon Entropy, Nucleotide Diversity (π), and Rare Haplotype Load (RHL).
*   **Interactive Visualizations:**
    *   **Mutational Hotspots:** Per-site Shannon entropy plot to identify hypervariable regions.
    *   **Variation Graph:** A dynamic graph showing the consensus backbone and the top 10 mutational bubbles (colored by frequency).
    *   **Metadata Correlation:** Imports environmental metadata (e.g., temperature, host) to visualize correlations with sequence diversity.

## 🛠️ Pipeline Architecture

1.  **Map:** Raw paired-end FASTQ reads are mapped to a reference genome using `minimap2`.
2.  **Trim:** Synthetic ARTIC amplicon primers are removed using `ivar trim` to prevent false-positive variant calls.
3.  **Call:** Low-frequency variants (down to 1%) are identified using `ivar variants`.
4.  **Reconstruct:** `CliqueSNV` attempts to link mutations that co-occur on the same reads. If it fails (due to high diversity), the pipeline filters variants >2% and generates a simulated haplotype cloud based on allele frequencies.
5.  **Visualize:** A Streamlit dashboard processes the output FASTA and metadata.

## ⚙️ Installation

We recommend using a Conda environment to manage the bioinformatics dependencies.

```bash
conda create -n quasispecies -c bioconda -c conda-forge python=3.10 minimap2 samtools ivar cliquesnv biopython pandas numpy scipy streamlit plotly
conda activate quasispecies
