#! /bin/sh

# In AnalysisCode directory

# C-DMR
for x in $(seq 5); do
    for i in $(seq 5); do
        Rscript trans_dmr_pval.R ../ProcessedData/c_dmr_allC.tsv.gz chr$x chr$i \
        ../ProcessedData/trans_analysis_data/cdmr_chr$x_vs_chr$i_sig_pvals.tsv \
        ../ProcessedData/trans_analysis_data/cdmr_chr$x_vs_chr$i_high_cor.tsv \
        ../ProcessedData/trans_analysis_data/cdmr_chr$x_vs_chr$i_pvals_insertions.tsv \
        ../ProcessedData/trans_analysis_data/cdmr_chr$x_vs_chr$i_pvals_deletions.tsv &
    done
    wait
done

cd ../ProcessedData/trans_analysis_data
cat *_sig_pvals.tsv | grep -v "dmr_chr" > ../cdmr_sig_pvals.tsv
cat *_high_cor.tsv | grep -v "dmr_chr" > ../cdmr_high_cor.tsv
rm *_sig_pvals.tsv *_high_cor.tsv
cd -

# CG-DMR
for x in $(seq 5); do
    for i in $(seq 5); do
        Rscript trans_dmr_pval.R ../ProcessedData/cg_dmrs_allC.tsv.gz chr$x chr$i \
        ../ProcessedData/trans_analysis_data/cgdmr_chr$x_vs_chr$i_sig_pvals.tsv \
        ../ProcessedData/trans_analysis_data/cgdmr_chr$x_vs_chr$i_high_cor.tsv \
        ../ProcessedData/trans_analysis_data/cgdmr_chr$x_vs_chr$i_pvals_insertions.tsv \
        ../ProcessedData/trans_analysis_data/cgdmr_chr$x_vs_chr$i_pvals_deletions.tsv &
    done
    wait
done

cd ../ProcessedData/trans_analysis_data
cat *_sig_pvals.tsv | grep -v "dmr_chr" > ../cgdmr_sig_pvals.tsv
cat *_high_cor.tsv | grep -v "dmr_chr" > ../cgdmr_high_cor.tsv
rm *_sig_pvals.tsv *_high_cor.tsv
gzip *
cd -