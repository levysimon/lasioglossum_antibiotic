# Antibiotic impact on gut microbiota and social behavior in the wild bee *Lasioglassum malachurum*

We combined behavioral and activity tracking, survival analysis, qPCR, and amplicon sequencing to evaluate the effects of antibiotic exposure on the gut microbiota and social dynamics of the wild bee *Lasioglossum malachurum*.


## Pipeline Structure
* **`activity-and-behavior/`**: Parses raw behavioral event logs (30 min arenas), loads inter-individual distance and locomotor activity data, and fits Negative Binomial GLMMs (glmmTMB) and Gaussian LMMs (lmer) under 3 × 2 and 2 × 2 factorial designs. Includes a multivariate PCA on 7 metrics with PERMANOVA, and produces a combined 6-panel figure.

* **`amplicon-sequencing/`**: Runs a full DADA2 pipeline (primer trimming with Cutadapt, quality filtering, error learning, denoising, merging, chimera removal) and assigns taxonomy against the Silva v138.2 database. Constructs a phyloseq object, removes contaminants using decontam (frequency and prevalence methods), and normalizes ASV counts by qPCR 16S copy numbers to produce absolute abundances. Community composition is visualized via PCoA (Bray-Curtis) and stacked bar plots. Differences between treatment groups are tested using PERMANOVA (adonis2) and homogeneity of dispersion (betadisper). Differentially abundant ASVs are identified using permutation ANOVAs with Tukey post-hoc tests and FDR correction, and a Wolbachia-free sensitivity analysis is run in parallel.
  
* **`qpcr/`**: Aggregates raw plate data, handles absolute quantification via standard curves, and fits Linear Mixed-Effects Models (`lme4::lmer`) to evaluate 16S rRNA gene copy numbers.
* **`survival-analysis/`**: Cleans longitudinal monitoring data, runs Log-Rank tests, and generates Kaplan-Meier survival curves.


------------------------------------------------------------------------
