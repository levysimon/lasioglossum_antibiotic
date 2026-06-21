# Antibiotic impact on gut microbiota and social behavior in the wild bee *Lasioglassum malachurum*

We combined behavioral and activity tracking, survival analysis, qPCR, and amplicon sequencing to evaluate the effects of antibiotic exposure on the gut microbiota and social dynamics of the wild bee *Lasioglossum malachurum*.


## Pipeline Structure
* **`amplicon-sequencing/`**: Houses the DADA2 and Phyloseq workflows for quality filtering, taxonomic assignment against the Silva database, and community composition analysis.
* **`qpcr/`**: Aggregates raw plate data, handles absolute quantification via standard curves, and fits Linear Mixed-Effects Models (`lme4::lmer`) to evaluate 16S rRNA gene copy numbers.
* **`survival-analysis/`**: Cleans longitudinal monitoring data, runs Log-Rank tests, and generates Kaplan-Meier survival curves.


------------------------------------------------------------------------
