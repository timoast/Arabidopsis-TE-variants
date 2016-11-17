## Code for running analyses on processed datasets

### TE insertion distribution  
**variants_per_accession.R**: Identify the number of true TE insertions and true TE deletions called for each accession in the population.  
**coverage_vs_variants.R**: Examine the relationship between sequencing depth of coverage and the number of TE variants found.  
**TE_insertion_distribution.R**: Find the minor allele frequency distribution, distribution of TE variant lengths, and the genomic distribution of TE insertions, TE deletions, rare variants and common variants.  
**TE_genomic_features.Rmd**: Find the frequency of TE variants occurring in different genomic features, such as genes, DHS sites, other TEs, pseudogenes and intergenic regions.  
**kee_regions.R**: Checks for TEPAV enrichments in the KNOT-engaged elements, calculate p-values and plot results as a histogram and a table.  
**quadrana_overlap.sh**: Find the overlap between our set of TE variants and those found by Quadrana et al.  
**coverage_vs_variants.R**: Plot the number of variants found in each accession vs the sequencing depth of coverage for that accession.  
**family_enrichment.R**: Calculate TE family and superfamily enrichments within the set of TE variants, for all TE variants as well as insertions and deletions separately. Plot the results as a barchart.  

### DNA methylation    
**burden_mC_DMRs.R**: Run burden analysis for TE insertions and deletions on C-DMR or CG-DMR methylation level and plot the results.  
**mC_heatmaps.R**: Plot the DNA methylation level in regions flanking TE variants as a heatmap.  
**trans_dmr_pval.R**: All-vs-all correlations for TE presence/absence and DMR methylation level. Calculate p-values for each pairwise correlation between TE and DMR, based on bootstrap. Save p-values and coordinates of significant correlations and strong correlations to specified file.  
**run_trans_analysis.sh**: Run all pairwise comparisons between chromosomes for the trans-DMR analysis.  
**plot_trans_scatter.R**: Plot significant TE-DMR associations as scatterplot.  
**TE_DMR_analysis.R**: Analyse the proportion of DMRs that are physically close to a TE variant and compare the the expected proportion. For DMRs within 1 kb of a TE variant, compare the DNA methylation levels in these DMRs in the presence/absence of the nearby TE. Plot the results as a cumulative density function, scatterplot, boxplot and density plot.  
**flanking_mC_correlations.R**: Correlate flanking +/- 200 bp mC levels with TE presence / absence, generate plots.  

### Gene expression  
**gene_expression_burden.R**: Run burden analysis for TE insertions and deletions on gene expression levels and plot the results.  
**gene_expression_volcano.R**: Run significance testing for effect of TE variants on nearby gene expression, plot results as a series of volcano plots.  
**RNA_heatmaps.R**: Plot heatmaps showing expression of selected genes with TE insertion effects on gene expression.
**resistance_phenotype.R**:  Download data from the plos genetics paper with pathogen resistance phenotype data, plot the percentage of accessions that were resistant to each pathogen and compare those accessions containing a TE insertion in AtRLP18 to those without the TE insertion.  