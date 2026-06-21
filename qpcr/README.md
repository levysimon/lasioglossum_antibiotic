# qPCR Quantification

This pipeline aggregates raw qPCR plate data from `qpcr/data/rawdata/` to calculate absolute bacterial 16S rRNA gene copy numbers. 
It automatically cleans raw $C_t$ signals, standardizes replicates, and accounts for false positives against negative controls. 
Using Linear Mixed Models (`lmer`), it tests the statistical impacts of antibiotic doses and gut reinoculations while exporting exploratory plots and ANOVA logs directly to `qpcr/result/`.

