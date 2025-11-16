# SNCM Genus-Level Microbiome Analysis

This repository contains the Sloan Neutral Community Model (SNCM) analysis pipeline for microbiome data at the genus level. The pipeline supports per-treatment and per-niche analyses with plots and summary statistics.

---

## 📂 Repository Structure

├── data/ # Input data
│ ├── phyloseq_object.rds
│ └── metadata.csv
├── out/ # Output folder for figures and tables
├── scripts/
│ └── SNCM_genus_analysis.R
├── README.md
└── .gitignore

yaml
Copy code

---

## 🛠️ Requirements

R (≥ 4.2) with the following packages:

```r
install.packages(c("dplyr", "ggplot2", "phyloseq", "MicEco"))
⚡ Usage
Place your input files in the data/ folder:

phyloseq_object.rds — phyloseq object containing ASV counts and taxonomy.

metadata.csv — sample metadata with columns like Treatment and Niche.

Run the SNCM analysis script:

r
Copy code
source("scripts/SNCM_genus_analysis.R")
Output will be saved in the out/ folder:

SNCM_AllNiches_Genus_fitclass.csv — SNCM results per genus.

SNCM_AllNiches_Genus_m_R2_summary.csv — immigration rate (m) and R² summary.

Faceted plots per treatment and 8-panel combined figure.

sessionInfo.txt for reproducibility.

🖼️ Outputs
Plots for individual treatments and combined dataset.

8-panel faceted plots showing m and R² per treatment.

CSV summaries for further analysis or reporting.

🧪 Notes
Low-prevalence taxa (present in <3 samples) are filtered.

Genera are classified as Above, Neutral, or Below the neutral model.

Manual pseudo-R² is computed for model quality assessment.