# This script gets the DNA methylation levels 
# for C-DMRs and CG-DMRs in each accession out of the 
# database created using previous scripts, for each
# DNA methylation context, and for all contexts combined
# MySQL database parameters must be added before running

# C-DMR methylation levels
python get_mc_DMRs.py -x [host] \
-u [user] -p [password] - d [mC_database] -f ../RawData/c_dmrs.tsv.gz \
> ../ProcessedData/c_dmr_allC.tsv

python get_mc_DMRs_context.py -x [host] \
-u [user] -p [password] - d [mC_database] -f ../RawData/c_dmrs.tsv.gz -c CG \
> ../ProcessedData/c_dmr_CG.tsv

python get_mc_DMRs_context.py -x [host] \
-u [user] -p [password] - d [mC_database] -f ../RawData/c_dmrs.tsv.gz -c CHG \
> ../ProcessedData/c_dmr_CHG.tsv

python get_mc_DMRs_context.py -x [host] \
-u [user] -p [password] - d [mC_database] -f ../RawData/c_dmrs.tsv.gz -c CHH \
> ../ProcessedData/c_dmr_CHH.tsv

# CG-DMR methylation levels
python get_mc_DMRs.py -x [host] \
-u [user] -p [password] - d [mC_database] -f ../RawData/cg_dmrs.tsv.gz \
> ../ProcessedData/cg_dmr_allC.tsv

python get_mc_DMR_context.py -x [host] \
-u [user] -p [password] - d [mC_database] -f ../RawData/cg_dmrs.tsv.gz -c CG \
> ../ProcessedData/cg_dmr_CG.tsv

python get_mc_DMR_context.py -x [host] \
-u [user] -p [password] - d [mC_database] -f ../RawData/cg_dmrs.tsv.gz -c CHG \
> ../ProcessedData/cg_dmr_CHG.tsv

python get_mc_DMR_context.py -x [host] \
-u [user] -p [password] - d [mC_database] -f ../RawData/cg_dmrs.tsv.gz -c CHH \
> ../ProcessedData/cg_dmr_CHH.tsv