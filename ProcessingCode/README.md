## Code for processing raw data into processed datasets.

### TE variant processing  
**process_TEPAV_data.R**: Classify TE variants as being either due to the deletion of a common TE, or the insertion of a new TE. Add in SNP association data.   
**split_true_deletions.sh**:  Create separate files for TE insertions and TE deletions.  
**reformat_tepav_matrix.py**: Reformat the TE variants file into a binary matrix.  

### Gene feature intersections with TE variants  
**make_genomic_feature_files.sh**:  Create annotation files for gene upstream and downstream regions, exons, introns, UTRs, pseudogenes, DHS sites and TEs. Intersect each with the TE variants. Count the number of intersections for each feature.  
**add_feature_between.py**: Add intron coordinates to gff file, ie the coordinates between exons of the same gene.  
**gene_flanking.py**: Get the upstream or downstream 2 kb regions of a gene in a gff file.  

### DNA methylation    
**bootstrap_te_variants.sh**: Randomly select regions of the genome and intersect with DMR coordinates to determine the probability of the observed TE variant overlaps with DMRs.  
**choose_random_coords.py**: Randomly select genomic coordinates of a given length and number.  
**closest_DMRs.sh**: Get the closest DMR to each TE variant, and the distnance to that DMR.  
**get_mc_DMRs_context.py**: Get the DNA methylation level within given coordinates (eg DMRs) from a MySQL database for a given sequence context.  
**get_mc_DMRs.py**: Get the DNA methylation level within given coordinates (eg DMRs) from a MySQL database in all sequence contexts.  **get_mc_flanking.py**: Get DNA methylation levels flanking a given set of coordinates, excluding the coordinates themselves.  
**get_mc_te.py**: Get DNA methylation levels flanking a given set of coordinates, including the coordinates themselves.  
**run_flanking_mc.sh**: run get\_mc\_flanking.py for TE insertions and TE deletions in each DNA methylation context.  

### TE variant stats  
**get_counts.py**: Count values for each accession and print to stdout.  
**get_cov.py**:  Process the sequencing coverage for each accession.  
**get_coverage_all_accessions.sh**: Find the sequencing depth of coverage for all accessions.  

### KEE regions  
**intersect_KEE_regions.sh**:  Find the number of TE variants in each Knot-engaged element region, and randomly select regions from the genome to make a comparison.  