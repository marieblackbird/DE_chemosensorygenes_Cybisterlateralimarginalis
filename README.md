# Expression of chemosensory genes in larval and adult cephalic appendages in a diving beetle _Cybister lateralimarginalis_

This repository contains scripts and resources for the assembly, annotation, quantification, and differential expression analysis of the transcriptome of _Cybister lateralimarginalis_ using Illumina sequencing data from larvae and adults, categorized by organ.

## Contents

1. **Transcriptome Assembly**
   - **`1_Assembly.sh`**: A shell script for assembling the transcriptome of _Cybister lateralimarginalis_ using Illumina sequencing data.

2. **Annotation of Chemosensory Genes**
   - **`2_Annotation.sh`**: A shell script for the annotation of chemosensory genes in the assembled transcriptome.

3. **Quantification**
   - **`3_Quantification.sh`**: A shell script for quantifying gene expression levels across different samples.

4. **Differential Expression Analysis**
   - **`4_DE_Analysis.R`**: An R script for performing differential expression analysis between different organs and developmental stages.

   **Utility Scripts**
   - **`utils/`**: Small utility scripts used during the annotation and analysis process.

## Usage

Each script in this repository is intended to be run sequentially to perform the full analysis pipeline. Start with `1_Assembly.sh` to assemble the transcriptome, then proceed with `2_Annotation.sh` to annotate the genes, followed by `3_Quantification.sh` for quantification, and finally, use `4_DE_Analysis.R` for differential expression analysis.

The `utils` directory contains small helper scripts like `parse_blast_output.py`, which are called during the main analysis steps to handle specific tasks such as file parsing.

## Getting Started

To get started, clone the repository and follow the steps outlined in the scripts. Each script is designed to be self-contained and can be run independently, assuming the necessary dependencies are installed and the input files are correctly formatted.

## Dependencies

- **Python**: Required for `parse_blast_output.py`.
- **R**: Required for `4_DE_Analysis.R`.
- **Shell**: Required for the assembly, annotation, and quantification scripts (`1_Assembly.sh`, `2_Annotation.sh`, `3_Quantification.sh`).

# About this Repository

This repository contains the research data associated with the project Evo‐AquaSense, funded by the National Agency for Research (ANR) ANR-19-CE02-003-01. The Evo‐AquaSense project was conducted by the Institut de Systématique, Evolution, Biodiversité (ISYEB UMR 7205) and the Institute of Ecology and Environmental Sciences of Paris (iEES-Paris).
