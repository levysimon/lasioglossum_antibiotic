# Antibiotic impact on gut microbiota and social behavior in the wild bee *Lasioglossum malachurum*
---

## Summary

We combined behavioral and activity tracking, survival analysis, qPCR, and amplicon sequencing to evaluate the effects of antibiotic exposure on the gut microbiota and social dynamics of the wild bee *Lasioglossum malachurum*. Colonies were exposed to a 3 × 2 factorial design varying antibiotic dose (zero, low, high) and gut re-inoculation status (no, yes). All analyses are fully reproducible from the raw data using the R scripts provided in each subdirectory.

------------------------------------------------------------------------

## Repository Structure

-   **`activity-and-behavior/`**: Parses raw behavioral event logs (30 min arenas), loads inter-individual distance and locomotor activity data, and fits Negative Binomial GLMMs (`glmmTMB`) and Gaussian LMMs (`lmer`) under 3 × 2 and 2 × 2 factorial designs. Includes a multivariate PCA on 7 metrics with PERMANOVA, and produces a combined 6-panel figure.

-   **`amplicon-sequencing/`**: Runs a full DADA2 pipeline (primer trimming with Cutadapt, quality filtering, error learning, denoising, merging, chimera removal) and assigns taxonomy against the Silva v138.2 database. Constructs a phyloseq object, removes contaminants using decontam (frequency and prevalence methods), and normalizes ASV counts by qPCR 16S copy numbers to produce absolute abundances. Community composition is visualized via PCoA (Bray-Curtis) and stacked bar plots. Differences between treatment groups are tested using PERMANOVA (`adonis2`) and homogeneity of dispersion (`betadisper`). Differentially abundant ASVs are identified using permutation ANOVAs with Tukey post-hoc tests and FDR correction, and a Wolbachia-free sensitivity analysis is run in parallel.

-   **`qpcr/`**: Aggregates raw plate data, handles absolute quantification via standard curves, and fits Linear Mixed-Effects Models (`lme4::lmer`) to evaluate 16S rRNA gene copy numbers.

-   **`survival-analysis/`**: Cleans longitudinal monitoring data, runs Log-Rank tests, and generates Kaplan-Meier survival curves.

------------------------------------------------------------------------

## How to Run

Each subdirectory contains a self-contained R Markdown (`.Rmd`) script. To reproduce the analyses:

1.  Clone or download the repository.
2.  Open the project in RStudio and set the working directory to the repository root.
3.  All file paths are managed with the `here` package and resolve automatically relative to the project root.
4.  Run each `.Rmd` file in the following recommended order:
    1.  `survival-analysis/`
    2.  `activity-and-behavior/`
    3.  `qpcr/`
    4.  `amplicon-sequencing/`
5.  For `amplicon-sequencing/`, raw sequencing data must first be downloaded from the SRA database (accession: **SRP666708**) and placed in `amplicon-sequencing/data/RawFiles/`. The Silva v138.2 reference databases must be downloaded from <https://zenodo.org/records/14169026> and placed in `amplicon-sequencing/data/reference/`. The processed phyloseq object (`ps_dada2taxa_Lasio.rds`) is already available, the DADA2 steps can be skipped and the analysis started directly from the `# The analysis can be started from here` section.

------------------------------------------------------------------------

## Software and Package Versions

All analyses were run in **R v4.4.3**. Key packages:

| Package                  | Version     | Use                      |
|--------------------------|-------------|--------------------------|
| `dada2`                  | ≥1.26       | Amplicon denoising       |
| `phyloseq`               | ≥1.42       | Microbiome data handling |
| `decontam`               | ≥1.18       | Contaminant removal      |
| `DESeq2`                 | ≥1.38       | Differential abundance   |
| `vegan`                  | ≥2.6        | PERMANOVA, ordination    |
| `glmmTMB`                | ≥1.1        | Negative binomial GLMMs  |
| `lme4` / `lmerTest`      | ≥1.1        | Linear mixed models      |
| `emmeans`                | ≥1.8        | Post-hoc comparisons     |
| `car`                    | ≥3.1        | Wald Chi-square tests    |
| `DHARMa`                 | ≥0.4        | Residual diagnostics     |
| `survival` / `survminer` | ≥3.5 / ≥0.4 | Survival analysis        |
| `dimensio`               | ≥0.4        | PCA                      |
| `ggplot2`                | ≥3.4        | Visualizations           |
| `patchwork` / `cowplot`  | ≥1.1        | Figure assembly          |
| `here`                   | ≥1.0        | Path management          |

## Full session info for each module is printed at the end of each `.Rmd` file.

## Data Availability

Raw sequencing data are deposited at the SRA database under accession **SRP666708**.

All other raw data files (behavioral logs, qPCR plates, survival monitoring) are included in the respective `data/` subdirectories.

Video recordings of the bees behaviors from which the behavioral event logs and activity-tracking outputs in **`activity-and-behavior/`** were derived are archived on Zenodo: <https://zenodo.org/records/20828421>

------------------------------------------------------------------------

## Links

-   Silva v138.2 reference database: <https://zenodo.org/records/14169026>
-   SRA accession: SRP666708

------------------------------------------------------------------------

*Author information and publication links will be added upon acceptance.*
