# An Integrative Systems Biology Framework for Prioritizing Candidate Molecular Mediators of Photobiomodulation in Meningeal Lymphatic Vessels

[![Repository](https://img.shields.io/badge/GitHub-TeslaHelper%2FBioinformatic-181717?logo=github)](https://github.com/TeslaHelper/Bioinformatic)
[![GEO GSE245658](https://img.shields.io/badge/NCBI_GEO-GSE245658-2ea44f)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE245658)
[![GEO GSE217276](https://img.shields.io/badge/NCBI_GEO-GSE217276-2ea44f)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE217276)

This repository contains the complete, modularized source datasets, split CSV tables, supplementary tables, and executable R/Python scripts for the manuscript:

> **An Integrative Systems Biology Framework for Prioritizing Candidate Molecular Mediators of Photobiomodulation in Meningeal Lymphatic Vessels**

Every sheet from the master workbooks and supplementary documents has been cleaned, standardized (with Row 1 as a machine-readable header), and organized into **9 numbered thematic modules** that follow the exact progression of the manuscript (`Sections 2.1–2.12` and `Sections 3.1–3.11`).

---

## Repository Architecture

```text
Bioinformatic/
├── README.md                                                     # Master documentation & manuscript crosswalk
├── 01_Interactome_Reconstruction_and_Topology/                   # Sections 2.1–2.4, 3.1–3.2 | Figs 1–3 | Table 1 | Tables S1, S5, S13, S14
├── 02_Candidate_Prioritization_and_AlphaMissense/                # Sections 2.5–2.6, 3.3 | Figs 3–4 | Table 2 | Tables S5–S6
├── 03_Functional_Enrichment_and_Multifunctionality/              # Sections 2.7, 2.9, 3.4–3.6 | Figs 5–7 | Tables S2–S4
├── 04_In_Silico_Perturbation_and_Path_Redundancy/                # Sections 2.8, 3.7 | Fig 8 | Table S7
├── 05_Boolean_Network_Dynamics/                                  # Sections 2.10, 3.8 | Figs 9–10 | Table 3 | Tables S10–S11
├── 06_Transcriptomic_Benchmarking_GSE245658_GSE217276/           # Sections 2.11, 3.9 | Figs 11–12 | Tables 4–5 | Table S8 | Figs S1–S3
├── 07_Continuous_ODE_Modeling_and_Optuna_Calibration/            # Sections 2.12, 3.10–3.11 | Figs 13–14 | Table 6 | Tables S9, S12
├── 08_Executable_Analysis_Scripts/                               # R (BoolNet, DESeq2, clusterProfiler) & Python (Optuna ODE) scripts
```

---

## Complete File-to-Manuscript Crosswalk

| Module Folder | File Name | Manuscript Figure / Table | Dimensions | Description |
| :--- | :--- | :--- | :--- | :--- |
| `01_Interactome_Reconstruction_and_Topology` | [`TableS01_15_Mechanistic_Candidate_Genes.csv`](01_Interactome_Reconstruction_and_Topology/TableS01_15_Mechanistic_Candidate_Genes.csv) | **Supplementary Table S1 / Section 2.1** | 14 rows × 3 cols | 15-gene mechanistic candidate subset with functional categories and biological roles |
| `01_Interactome_Reconstruction_and_Topology` | [`Table01_Fig01BC_Bridge_Interactome_Topology_vs_ER_Null.csv`](01_Interactome_Reconstruction_and_Topology/Table01_Fig01BC_Bridge_Interactome_Topology_vs_ER_Null.csv) | **Table 1 / Figure 1B–C** | 8 rows × 3 cols | 205-node, 4,074-edge bridge interactome global topology vs. 1,000 Erdos-Renyi random graphs |
| `01_Interactome_Reconstruction_and_Topology` | [`Fig01A_205Node_Degree_Distribution.csv`](01_Interactome_Reconstruction_and_Topology/Fig01A_205Node_Degree_Distribution.csv) | **Figure 1A** | 89 rows × 2 cols | Degree frequency distribution across all 205 bridge-interactome nodes |
| `01_Interactome_Reconstruction_and_Topology` | [`Fig02_TableS05_15Candidate_105_Pairwise_STRING_Edges.csv`](01_Interactome_Reconstruction_and_Topology/Fig02_TableS05_15Candidate_105_Pairwise_STRING_Edges.csv) | **Figure 2 / Supplementary Table S5** | 105 rows × 5 cols | All 105 undirected STRING v11.5 pairwise interaction scores among the 15 mechanistic candidates |
| `01_Interactome_Reconstruction_and_Topology` | [`TableS13_Full_205Node_Topology_and_NetworkScores.csv`](01_Interactome_Reconstruction_and_Topology/TableS13_Full_205Node_Topology_and_NetworkScores.csv) | **Figure 3 / Supplementary Table S13** | 205 rows × 11 cols | Complete 205-node bridge interactome degree, betweenness, PageRank, N_PBM, N_Lymph, Bridge Score, and NetworkScore* |
| `01_Interactome_Reconstruction_and_Topology` | [`TableS14_Full_205Node_4074_STRING_Edges.csv`](01_Interactome_Reconstruction_and_Topology/TableS14_Full_205Node_4074_STRING_Edges.csv) | **Supplementary Table S14 / Section 2.3** | 4074 rows × 3 cols | Complete 4,074 interaction edges (combined_score >= 0.400) of the 205-node human PBM-MLV bridge interactome |
| `02_Candidate_Prioritization_and_AlphaMissense` | [`Table02_Top20_GoldenScore_Prioritization_Harmonized.csv`](02_Candidate_Prioritization_and_AlphaMissense/Table02_Top20_GoldenScore_Prioritization_Harmonized.csv) | **Table 2 / Figure 4 / Supplementary Table S6** | 20 rows × 16 cols | Integrated Golden Score prioritization across all 20 candidates with N_PBM, N_Lymph, scaled 70%/30% components, and rank ranges |
| `02_Candidate_Prioritization_and_AlphaMissense` | [`Table02_Source_20Gene_Topology_and_AlphaMissense.csv`](02_Candidate_Prioritization_and_AlphaMissense/Table02_Source_20Gene_Topology_and_AlphaMissense.csv) | **Table 2 (Source Rederivation)** | 20 rows × 11 cols | Full-205 network rederivation of topology and verified AlphaMissense EVS across the 20 candidates |
| `02_Candidate_Prioritization_and_AlphaMissense` | [`Fig02_15Candidate_Node_Properties.csv`](02_Candidate_Prioritization_and_AlphaMissense/Fig02_15Candidate_Node_Properties.csv) | **Figure 2** | 15 rows × 12 cols | Node centrality and score attributes for the 15-gene mechanistic subset |
| `02_Candidate_Prioritization_and_AlphaMissense` | [`Fig04A_AlphaMissense_vs_Topology_Scatter.csv`](02_Candidate_Prioritization_and_AlphaMissense/Fig04A_AlphaMissense_vs_Topology_Scatter.csv) | **Figure 4** | 15 rows × 7 cols | AlphaMissense mean pathogenicity score vs. network topology scatter data |
| `02_Candidate_Prioritization_and_AlphaMissense` | [`Fig04B_GoldenScore_Stacked_Components.csv`](02_Candidate_Prioritization_and_AlphaMissense/Fig04B_GoldenScore_Stacked_Components.csv) | **Figure 4** | 15 rows × 9 cols | Stacked component decomposition of Golden Score across candidates |
| `02_Candidate_Prioritization_and_AlphaMissense` | [`TableS05_GoldenScore_Weight_Sensitivity_Rank_Matrix.csv`](02_Candidate_Prioritization_and_AlphaMissense/TableS05_GoldenScore_Weight_Sensitivity_Rank_Matrix.csv) | **Supplementary Table S5 / Section 3.3** | 15 rows × 6 cols | Candidate rank matrix across 5 alpha:beta weighting splits (70:30, 60:40, 80:20, 50:50, 90:10) |
| `02_Candidate_Prioritization_and_AlphaMissense` | [`TableS06_GoldenScore_Weight_Sensitivity_LongForm.csv`](02_Candidate_Prioritization_and_AlphaMissense/TableS06_GoldenScore_Weight_Sensitivity_LongForm.csv) | **Supplementary Table S6 / Section 3.3** | 75 rows × 6 cols | Long-form sensitivity table across all 5 alpha:beta weighting scenarios |
| `03_Functional_Enrichment_and_Multifunctionality` | [`Fig05_GO_BP_Highlighted_Terms.csv`](03_Functional_Enrichment_and_Multifunctionality/Fig05_GO_BP_Highlighted_Terms.csv) | **Figure 5 / Section 3.4** | 12 rows × 10 cols | Principal Gene Ontology Biological Process terms displayed in Figure 5 |
| `03_Functional_Enrichment_and_Multifunctionality` | [`Fig06_KEGG_Highlighted_Pathways.csv`](03_Functional_Enrichment_and_Multifunctionality/Fig06_KEGG_Highlighted_Pathways.csv) | **Figure 6 / Section 3.4** | 12 rows × 9 cols | Principal KEGG pathways displayed in Figure 6 dot plot |
| `03_Functional_Enrichment_and_Multifunctionality` | [`TableS02_GO_Enrichment_All_1857_Significant_Terms.csv`](03_Functional_Enrichment_and_Multifunctionality/TableS02_GO_Enrichment_All_1857_Significant_Terms.csv) | **Supplementary Table S2 / Section 3.4** | 1857 rows × 13 cols | Complete whole-genome GO enrichment output (1,857 significant terms: 1,732 BP + 64 MF + 61 CC) |
| `03_Functional_Enrichment_and_Multifunctionality` | [`TableS03b_KEGG_Enrichment_All_119_Significant_Pathways.csv`](03_Functional_Enrichment_and_Multifunctionality/TableS03b_KEGG_Enrichment_All_119_Significant_Pathways.csv) | **Supplementary Table S3b / Section 3.4** | 119 rows × 14 cols | Complete whole-genome KEGG pathway enrichment output (119 significant pathways) |
| `03_Functional_Enrichment_and_Multifunctionality` | [`Fig07_TableS03a_Restricted_205Node_Universe_Comparison.csv`](03_Functional_Enrichment_and_Multifunctionality/Fig07_TableS03a_Restricted_205Node_Universe_Comparison.csv) | **Figure 7 / Supplementary Table S3a / Section 3.5** | 4 rows × 4 cols | Comparison of significant GO BP and KEGG terms under whole-genome vs. 205-node bridge-network background |
| `03_Functional_Enrichment_and_Multifunctionality` | [`TableS04_GenomeWide_GO_Multifunctionality_Index_18870_Genes.csv`](03_Functional_Enrichment_and_Multifunctionality/TableS04_GenomeWide_GO_Multifunctionality_Index_18870_Genes.csv) | **Supplementary Table S4 / Section 3.6** | 18870 rows × 4 cols | Genome-wide Gene Ontology BP annotation counts and multifunctionality percentiles across 18,870 human genes |
| `04_In_Silico_Perturbation_and_Path_Redundancy` | [`Fig08_TableS07_COX4I1_FLT4_Path_Redundancy_Loss.csv`](04_In_Silico_Perturbation_and_Path_Redundancy/Fig08_TableS07_COX4I1_FLT4_Path_Redundancy_Loss.csv) | **Figure 8 / Supplementary Table S7 / Section 3.7** | 19 rows × 6 cols | Single-node deletion path-redundancy loss across COX4I1-FLT4 simple paths (cutoff <= 5 edges; baseline 573 paths) |
| `05_Boolean_Network_Dynamics` | [`Table03_Boolean_Network_Architecture_and_Rules.csv`](05_Boolean_Network_Dynamics/Table03_Boolean_Network_Architecture_and_Rules.csv) | **Table 3 / Equations (6)–(18) / Section 2.10** | 14 rows × 5 cols | Complete 13-node + auxiliary memory M synchronous Boolean transition rules and post-translational state definitions |
| `05_Boolean_Network_Dynamics` | [`Fig09A_TableS10a_Boolean_Trajectory_WildType.csv`](05_Boolean_Network_Dynamics/Fig09A_TableS10a_Boolean_Trajectory_WildType.csv) | **Figure 9A / Supplementary Table S10a** | 15 rows × 14 cols | Wild-type memory-augmented synchronous Boolean state trajectory |
| `05_Boolean_Network_Dynamics` | [`Fig09B_TableS10b_Boolean_Trajectory_ACTB_Clamped.csv`](05_Boolean_Network_Dynamics/Fig09B_TableS10b_Boolean_Trajectory_ACTB_Clamped.csv) | **Figure 9B / Supplementary Table S10b** | 15 rows × 14 cols | ACTB-clamped (ACTB_KO) Boolean trajectory showing complete downstream blockade |
| `05_Boolean_Network_Dynamics` | [`Fig09C_TableS10c_Boolean_Trajectory_IL6_Clamped.csv`](05_Boolean_Network_Dynamics/Fig09C_TableS10c_Boolean_Trajectory_IL6_Clamped.csv) | **Figure 9C / Supplementary Table S10c** | 15 rows × 14 cols | IL6-clamped (IL6_KO) Boolean trajectory demonstrating compensation by the parallel TNF branch |
| `05_Boolean_Network_Dynamics` | [`Fig09D_Fig10B_TableS11a_Boolean_Trajectory_MemoryFree_Naive.csv`](05_Boolean_Network_Dynamics/Fig09D_Fig10B_TableS11a_Boolean_Trajectory_MemoryFree_Naive.csv) | **Figures 9D & 10A–B / Supplementary Table S11a** | 15 rows × 5 cols | Memory-free (Hysteresis_Naive) Boolean trajectory showing full collapse at t = 13 |
| `05_Boolean_Network_Dynamics` | [`Fig10C_TableS11b_Boolean_Trajectory_Memory_Augmented.csv`](05_Boolean_Network_Dynamics/Fig10C_TableS11b_Boolean_Trajectory_Memory_Augmented.csv) | **Figure 10C / Supplementary Table S11b** | 15 rows × 5 cols | Memory-augmented hysteresis trajectory maintaining downstream drainage after pulse withdrawal |
| `06_Transcriptomic_Benchmarking_GSE245658_GSE217276` | [`Table04_Transcriptomic_Experimental_Design_and_QC.csv`](06_Transcriptomic_Benchmarking_GSE245658_GSE217276/Table04_Transcriptomic_Experimental_Design_and_QC.csv) | **Table 4 / Sections 2.11 & 3.9** | 7 rows × 3 cols | Experimental design, sample sizes, PCA outlier handling, and DESeq2 modeling specifications for GSE245658 and GSE217276 |
| `06_Transcriptomic_Benchmarking_GSE245658_GSE217276` | [`Table05_Dual_Dataset_Candidate_and_Junctional_DESeq2_Summary.csv`](06_Transcriptomic_Benchmarking_GSE245658_GSE217276/Table05_Dual_Dataset_Candidate_and_Junctional_DESeq2_Summary.csv) | **Table 5 / Figures 11–12 / Section 3.9** | 23 rows × 12 cols | Complete DESeq2 statistics across GSE245658 and GSE217276 for all 15 mechanistic candidates and 8 reference junctional/signaling genes |
| `06_Transcriptomic_Benchmarking_GSE245658_GSE217276` | [`Fig11_Fig12_Dual_Dataset_15Candidate_DESeq2.csv`](06_Transcriptomic_Benchmarking_GSE245658_GSE217276/Fig11_Fig12_Dual_Dataset_15Candidate_DESeq2.csv) | **Figures 11 & 12** | 27 rows × 11 cols | Dual-dataset differential expression summary plotted in Figures 11 and 12 |
| `06_Transcriptomic_Benchmarking_GSE245658_GSE217276` | [`TableS08_GSE245658_Full_Transcriptome_DESeq2_11804_Genes.csv`](06_Transcriptomic_Benchmarking_GSE245658_GSE217276/TableS08_GSE245658_Full_Transcriptome_DESeq2_11804_Genes.csv) | **Supplementary Table S8 / Figure 11A** | 11804 rows × 9 cols | Complete genome-wide DESeq2 differential expression output for GSE245658 (FAD_L vs. FAD) |
| `06_Transcriptomic_Benchmarking_GSE245658_GSE217276` | [`FigS01_GSE245658_VST_PCA_Coordinates.csv`](06_Transcriptomic_Benchmarking_GSE245658_GSE217276/FigS01_GSE245658_VST_PCA_Coordinates.csv) | **Supplementary Figure S1** | 10 rows × 6 cols | Sample-level PCA coordinates (PC1 and PC2) across WT, FAD, and FAD_L libraries in GSE245658 |
| `06_Transcriptomic_Benchmarking_GSE245658_GSE217276` | [`FigS02_GSE245658_Volcano_Candidates_and_Top_DEGs.csv`](06_Transcriptomic_Benchmarking_GSE245658_GSE217276/FigS02_GSE245658_Volcano_Candidates_and_Top_DEGs.csv) | **Supplementary Figure S2** | 15 rows × 5 cols | Volcano plot source table of candidate genes and top genome-wide differentially expressed transcripts |
| `06_Transcriptomic_Benchmarking_GSE245658_GSE217276` | [`FigS03_GSE245658_Candidate_VST_Zscore_Heatmap.csv`](06_Transcriptomic_Benchmarking_GSE245658_GSE217276/FigS03_GSE245658_Candidate_VST_Zscore_Heatmap.csv) | **Supplementary Figure S3** | 10 rows × 11 cols | Gene-wise z-scores of variance-stabilizing-transformed (VST) counts across all 10 GSE245658 libraries |
| `07_Continuous_ODE_Modeling_and_Optuna_Calibration` | [`Table06_Comparative_ODE_Regimes_Primary_vs_Secondary.csv`](07_Continuous_ODE_Modeling_and_Optuna_Calibration/Table06_Comparative_ODE_Regimes_Primary_vs_Secondary.csv) | **Table 6 / Equations (19)–(21) / Section 3.10** | 11 rows × 6 cols | Side-by-side comparative parameterization of Primary (Extended Remodeling) and Secondary (Physiological Turnover) ODE regimes |
| `07_Continuous_ODE_Modeling_and_Optuna_Calibration` | [`Table06A_Primary_Regime_Optuna_Parameters.csv`](07_Continuous_ODE_Modeling_and_Optuna_Calibration/Table06A_Primary_Regime_Optuna_Parameters.csv) | **Table 6 (Primary Regime Source)** | 5 rows × 6 cols | Primary extended-remodeling Optuna TPE parameter values |
| `07_Continuous_ODE_Modeling_and_Optuna_Calibration` | [`TableS09_Secondary_Regime_Physiological_Turnover_Parameters.csv`](07_Continuous_ODE_Modeling_and_Optuna_Calibration/TableS09_Secondary_Regime_Physiological_Turnover_Parameters.csv) | **Supplementary Table S9 / Table 6 (Secondary Regime)** | 5 rows × 6 cols | Secondary physiological-turnover Optuna parameter values (k_dep = 0.3726 h^-1, t1/2 ~ 112 min) |
| `07_Continuous_ODE_Modeling_and_Optuna_Calibration` | [`Fig13A_TableS12a_ODE_72h_Trajectory_1441_Points.csv`](07_Continuous_ODE_Modeling_and_Optuna_Calibration/Fig13A_TableS12a_ODE_72h_Trajectory_1441_Points.csv) | **Figure 13A / Supplementary Table S12a** | 1441 rows × 4 cols | Continuous 72-h ODE trajectory (1,441 time points, dt = 0.05 h) for S(t), A(t), and D(t) |
| `07_Continuous_ODE_Modeling_and_Optuna_Calibration` | [`Fig13B_TableS12b_ODE_Sensitivity_kdep_vs_Persistence.csv`](07_Continuous_ODE_Modeling_and_Optuna_Calibration/Fig13B_TableS12b_ODE_Sensitivity_kdep_vs_Persistence.csv) | **Figure 13B / Supplementary Table S12b** | 80 rows × 3 cols | Parameter sensitivity sweep of actin depolymerization rate k_dep vs. A(t) > 0.5 persistence duration |
| `07_Continuous_ODE_Modeling_and_Optuna_Calibration` | [`Fig14_Multiphase_Mechanism_Schematic_Nodes.csv`](07_Continuous_ODE_Modeling_and_Optuna_Calibration/Fig14_Multiphase_Mechanism_Schematic_Nodes.csv) | **Figure 14 / Section 3.11** | 14 rows × 4 cols | Multiphase systems-biology schematic node and phase table underlying Figure 14 |
| `08_Executable_Analysis_Scripts` | [`boolean_dynamics.R`](08_Executable_Analysis_Scripts/boolean_dynamics.R) | **Table 3 / Figures 9–10 / Section 2.10** | 199 rows × 1 cols | BoolNet v2.1.9 synchronous 13-node + auxiliary memory M Boolean network simulation script |
| `08_Executable_Analysis_Scripts` | [`run_optuna_actb.py`](08_Executable_Analysis_Scripts/run_optuna_actb.py) | **Table 6 / Figure 13 / Section 2.12** | 223 rows × 1 cols | Optuna v3.5 TPE calibration (400 trials, seed=42) and forward Euler ODE simulation script |
| `08_Executable_Analysis_Scripts` | [`deseq2_verification_GSE245658.R`](08_Executable_Analysis_Scripts/deseq2_verification_GSE245658.R) | **Tables 4–5 / Figure 11A / Section 2.11** | 303 rows × 1 cols | DESeq2 v1.42 reanalysis script for GSE245658 (5xFAD + 808 nm NIR light vs. untreated 5xFAD) |
| `08_Executable_Analysis_Scripts` | [`validate_hub_genes_217276.py`](08_Executable_Analysis_Scripts/validate_hub_genes_217276.py) | **Tables 4–5 / Figure 11B & 12 / Section 2.11** | 86 rows × 1 cols | Differential expression benchmark script for GSE217276 (intracranial AAV-IFNg vs. AAV-vehicle) |
| `08_Executable_Analysis_Scripts` | [`universe_correction.R`](08_Executable_Analysis_Scripts/universe_correction.R) | **Figure 7 / Supplementary Table S3 / Section 2.7** | 132 rows × 1 cols | Restricted 205-node reference universe over-representation and rank-based GSEA script |
---

## Reproducing the Computational Analyses

### 1. Synchronous Boolean Network Simulation (`08_Executable_Analysis_Scripts/boolean_dynamics.R`)
- **Environment**: R (>= 4.2) with `BoolNet` (v2.1.9)
- **Command**:
  ```bash
  Rscript 08_Executable_Analysis_Scripts/boolean_dynamics.R
  ```
- **Outputs**: Reproduces `Table 3` (`Equations 6–18`) and all five synchronous trajectories (`WT`, `ACTB_KO`, `IL6_KO`, `Hysteresis_Naive`, `Hysteresis_Memory`) in `05_Boolean_Network_Dynamics/`.

### 2. Dual-Regime Continuous ODE Calibration (`08_Executable_Analysis_Scripts/run_optuna_actb.py`)
- **Environment**: Python (>= 3.10) with `optuna` (v3.5), `numpy`, `pandas`, `scipy`
- **Command**:
  ```bash
  python3 08_Executable_Analysis_Scripts/run_optuna_actb.py
  ```
- **Outputs**: Reproduces `Table 6` (`Equations 19–21`), the 1,441-point 72-h forward Euler trajectory (`Fig13A_TableS12a_ODE_72h_Trajectory_1441_Points.csv`), and the $k_{\text{dep}}$ persistence sensitivity curve (`Fig13B_TableS12b_ODE_Sensitivity_kdep_vs_Persistence.csv`).

### 3. Transcriptomic Reanalysis & Positive Inflammatory Benchmark (`deseq2_verification_GSE245658.R` & `validate_hub_genes_217276.py`)
- **Environment**: R with `DESeq2` (v1.42) and `apeglm`; Python with `pydeseq2` / `pandas`
- **Command**:
  ```bash
  Rscript 08_Executable_Analysis_Scripts/deseq2_verification_GSE245658.R
  python3 08_Executable_Analysis_Scripts/validate_hub_genes_217276.py
  ```
- **Outputs**: Reproduces `Tables 4–5`, `Figures 11–12`, `Supplementary Table S8`, and `Supplementary Figures S1–S3` in `06_Transcriptomic_Benchmarking_GSE245658_GSE217276/`.

### 4. Restricted 205-Node Reference Universe Enrichment (`08_Executable_Analysis_Scripts/universe_correction.R`)
- **Environment**: R with `clusterProfiler` (v4.8.3) and `org.Hs.eg.db` (v3.18.0)
- **Command**:
  ```bash
  Rscript 08_Executable_Analysis_Scripts/universe_correction.R
  ```
- **Outputs**: Reproduces `Figure 7` and `Supplementary Table S3` in `03_Functional_Enrichment_and_Multifunctionality/`.
