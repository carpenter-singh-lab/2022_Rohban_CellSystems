# 2016_08_22_Gene_Compound

System Requirement:

R (3.6.0)
  dplyr (0.8.5)
  stringr (1.4.0)
  reshape2 (1.4.3)
  htmlTable (2.2.1)
  foreach (1.4.4)
  doMC (1.3.7)
  magrittr (1.5)
  readr (1.3.1)
  ggplot2 (3.3.0)
 
 All analyses are based on the R scripts in the archive folder. These are the main scripts/notebooks:
 
 CDRP_per_gene_analysis.Rmd : querying each gene for the compound matches and doing enrichment analysis based on known valid compound/gene pairs
 
 cdrp_rep_corr_against_density.R : Replicate Corr. vs Cell Count plot to verify that the high replicate correlation is not due to toxicity 
  
[Confluence project page](https://broadinstitute.atlassian.net/wiki/spaces/IP/pages/109838539/Gene-Compound+Association) - Broad internal lab notebook

The file in this repository named "corr_mat.csv.zip" contains all 69 genes in the experiment and their correlations to the compounds in the experiment (the 15,863 impactful compounds of the 30,616 set).

