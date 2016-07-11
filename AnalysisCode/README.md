## Code for running analyses on processed datasets

### TE insertion distribution  
**TE_insertion_distribution.R**: Find the minor allele frequency distribution, distribution of TE variant lengths, and the genomic distribution of TE insertions, TE deletions, rare variants and common variants. Find the frequency of TE variants occurring in different genomic features, such as genes, DHS sites, other TEs or pseudogenes. Compare rare variants and common variants, and compare TE insertions to TE deletions.  
**kee_regions.R**: Checks for TEPAV enrichments in the KNOT-engaged elements, calculate p-values and plot results as a histogram and a table.  
**quadrana_overlap.sh**: Find the overlap between our set of TE variants and those found by Quadrana et al.  
**coverage_vs_variants.R**: Plot the number of variants found in each accession vs the sequencing depth of coverage for that accession.  
**family_enrichment.R**: Calculate TE family and superfamily enrichments within the set of TE variants, for all TE variants as well as insertions and deletions separately. Plot the results as a barchart.  

### DNA methylation    
**burden_mC_DMRs.R**: Run burden analysis for TE insertions and deletions on C-DMR or CG-DMR methylation level and plot the results.  
**mC_heatmaps.R**: Plot the DNA methylation level in regions flanking TE variants as a heatmap.  
**trans_DMR_analysis.R**: All-vs-all correlations for TE presence/absence and DMR methylation level. Plot results as heatmap.  *Note - this script cannot be run directly to get the figures as output*  
**TE_DMR_analysis.R**: Analyse the proportion of DMRs that are physically close to a TE variant and compare the the expected proportion. For DMRs within 1 kb of a TE variant, compare the DNA methylation levels in these DMRs in the presence/absence of the nearby TE. Plot the results as a cumulative density function, scatterplot, boxplot and density plot.  
**flanking_mC_correlations.R**: Correlate flanking +/- 200 bp mC levels with TE presence / absence, generate plots.  

### Gene expression  
**gene_expression_burden.R**: Run burden analysis for TE insertions and deletions on gene expression levels and plot the results.  
**gene_expression_volcano.R**: Run significance testing for effect of TE variants on nearby gene expression, plot results as a series of volcano plots.  
**RNA_heatmaps.R**: Plot heatmaps showing expression of selected genes with TE insertion effects on gene expression.
**resistance_phenotype.R**:  Download data from the plos genetics paper with pathogen resistance phenotype data, plot the percentage of accessions that were resistant to each pathogen and compare those accessions containing a TE insertion in AtRLP18 to those without the TE insertion.  