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
```

## ▶️ Usage

1. **Run the Pipeline**
To process raw FASTQ reads, you need the forward (R1) and reverse (R2) reads, and a reference genome.

```bash
python fastq_to_fasta.py \
  --fastq1 raw_reads/SRR12447392_1.fastq \
  --fastq2 raw_reads/SRR12447392_2.fastq \
  --ref MN908947.3.fasta \
  --out real_haplotypes.fasta
```

2. **Launch the Dashboard**
Once the pipeline finishes and generates the real_haplotypes.fasta file, launch the web interface:

```bash
streamlit run app.py
```

Upload the generated FASTA file, and optionally upload a CSV containing metadata (columns must include HaplotypeID, Temperature_C, and Host) to see the correlation plots.

## 🧪 Testing with Simulated Data

If you want to test the tool without downloading large public datasets, use the built-in simulation scripts:

```bash
python simulate_data.py
python app.py # or streamlit run app.py
```

Upload the generated simulated_haplotypes.fasta to see a controlled, simple quasispecies cloud.

## 📊 Expected Outputs

* **Haplotype Count:** The number of unique viral genomes in the population.
* **Shannon Entropy:** The uncertainty (variation) at each genomic position. High spikes indicate mutational hotspots.
* **Nucleotide Diversity (π):** The average number of pairwise differences per site across the population.
* **Rare Haplotype Load (RHL):** The fraction of haplotypes below a 5% threshold, which is a strong marker of viral adaptability and mutagenic pressure.

## 👤 Author

Rajitha - [GitHub Profile](https://github.com/rajithadp)

## 📄 License

This project is licensed under the MIT License.
