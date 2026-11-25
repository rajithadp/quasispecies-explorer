# Quasispecies Explorer

A lightweight, interactive tool for quantifying and visualizing viral quasispecies diversity from aligned haplotype FASTA files.This tool takes raw FASTQ reads, maps them to a reference genome, trims primers, calls variants using iVar, and reconstructs haplotypes. It includes a fallback algorithm to handle highly diverse quasispecies samples that standard haplotype reconstructors struggle with.


## Features
*   **Population Metrics:** Calculates the total number of unique haplotypes, mean Shannon entropy, and nucleotide diversity (π).
*   **Mutational Hotspots:** Interactive Plotly graph showing per-site Shannon entropy across the viral genome to identify regions of high variation.
*   **Variation Graph:** A dynamic graph representation where shared genomic segments are baseline nodes and variants branch upwards, visually demonstrating the quasispecies structure.

## Usage
1. Install requirements: `pip install streamlit pandas numpy biopython scipy plotly`
2. Generate simulated data: `python simulate_data.py`
3. Run the app: `streamlit run app.py`
4. Upload `simulated_haplotypes.fasta` to the web interface.
