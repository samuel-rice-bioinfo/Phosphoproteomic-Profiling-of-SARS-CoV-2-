# Phosphoproteomics Analysis of SARS-CoV-2 Infection in Human Cells

This repository contains the Comet MS/MS configuration file and associated analysis scripts and summaries for a phosphoproteomics project focused on identifying phosphorylation events in human H1299 cells following SARS-CoV-2 infection.

---

## Project Overview

Mass spectrometry (MS)-based proteomics was used to identify and quantify phosphorylated peptides in SARS-CoV-2-infected cells compared to a control. Tandem Mass Tag (TMT) labeling enabled relative quantification, and the Comet search engine was employed for peptide-spectrum matching against the human proteome. The goal was to localize phosphorylation sites and explore infection-induced regulatory changes in cellular signaling.

---

## Comet Search Setup

The Comet parameter file (`comet.params`) includes the following key specifications:

* **Database:** Human proteome (`UP000005640_9606.fasta`)
* **TMTpro Modifications:** N-term (+304.207146), Lysine (+304.207146)
* **Phosphorylation Variable Mod:** STY residues (+79.966331)
* **Enzyme:** Trypsin (fully digested, max 2 missed cleavages)
* **Mass Tolerance:** ±10 ppm
* **Fragment Tolerance:** 0.02 Da, with b- and y-ions
* **FDR Estimation:** Decoy search enabled (target-decoy strategy)

---

## Analysis Workflow

1. **Peptide Identification:**
   Comet was run on `.mgf` files using the configured parameter file, identifying peptide-spectrum matches (PSMs).

2. **FDR Filtering:**
   PSMs were filtered at 1%, 5%, and 10% global FDR thresholds using decoy counts to control confidence in identifications.

3. **Phosphosite Localisation:**
   Lorikeet was used to manually inspect and validate phosphorylation sites on high-confidence peptides by visualizing b/y-ion fragmentation patterns.

4. **Quantitative Analysis (Part B):**

   * Log-transformed and median-normalized TMT intensity data
   * Differential analysis via unpaired t-tests
   * Volcano plots comparing SARS-CoV-1 vs control and SARS-CoV-2 vs control
   * Functional enrichment using DAVID with pathway clustering and heatmap visualization

---

## Key Results

* Over 10,000 PSMs identified, with >3,300 phosphopeptides detected
* High-confidence (1% FDR) peptides used for site-specific phosphorylation analysis
* Lorikeet confirmed multiple plausible phosphorylation sites on ZN335\_HUMAN, with Threonine 15 showing the strongest ion current support
* Volcano plots revealed increased phosphorylation activity in SARS-CoV-2 infection vs SARS-CoV-1
* Functional enrichment highlighted nuclear regulatory pathways (SARS-CoV-1) and cytoskeletal/signaling pathways (SARS-CoV-2)

---

## Tools and Technologies Used

* Comet (v2024.02 rev. 0)
* Lorikeet Viewer
* RStudio (ggplot2, heatmaps, data cleaning)
* DAVID Bioinformatics Resources
* UniProt, STRING, NCBI (sequence resources)

---

## Files Included

* `comet.params`: MS/MS search parameters
* `Phosphoproteomics_Part_A`: Identification & site localisation (Lorikeet)
* `Phosphoproteomics_Part_B_task1`: Differential expression & volcano plots
* `Phosphoproteomics_Part_B_task2_heatmap`: DAVID enrichment heatmap
* `Phosphoproteomics_Part_B_task2_david`: DAVID functional annotation summary

---

## Reference Highlights

* Gordon et al., 2020 – Protein interaction mapping of SARS-CoV-2
* Bouhaddou et al., 2020 – Global phosphorylation landscape of infection
* Sharma et al., 2012 – Lorikeet and phosphosite localisation
* Elias & Gygi, 2007 – Target-decoy FDR estimation
* Thompson et al., 2003 – TMT quantification method

