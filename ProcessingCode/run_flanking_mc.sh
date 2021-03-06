# No insertions
# allC
python get_mc_flanking.py \
-a [host] \
-p [password] \
-u [user] \
-d [database] \
-i ../ProcessedData/TEPAV_no_insertion.tsv.gz \
-o ../ProcessedData/flanking_allC_no_insertion.tsv \
-b 200 \
-n 20 \
-m ../RawData/all_mc_tables.txt

# CG
python get_mc_flanking.py \
-a [host] \
-p [password] \
-u [user] \
-d [database] \
-i ../ProcessedData/TEPAV_no_insertion.tsv.gz \
-o ../ProcessedData/flanking_CG_no_insertion.tsv \
-b 200 \
-n 20 \
-c CG \
-m ../RawData/all_mc_tables.txt

# CHG
python get_mc_flanking.py \
-a [host] \
-p [password] \
-u [user] \
-d [database] \
-i ../ProcessedData/TEPAV_no_insertion.tsv.gz \
-o ../ProcessedData/flanking_CHG_no_insertion.tsv \
-b 200 \
-n 20 \
-c CHG \
-m ../RawData/all_mc_tables.txt

# CHH
python get_mc_flanking.py \
-a [host] \
-p [password] \
-u [user] \
-d [database] \
-i ../ProcessedData/TEPAV_no_insertion.tsv.gz \
-o ../ProcessedData/flanking_CHH_no_insertion.tsv \
-b 200 \
-n 20 \
-c CHH \
-m ../RawData/all_mc_tables.txt

# True deletions
# allC
python get_mc_flanking.py \
-a [host] \
-p [password] \
-u [user] \
-d [database] \
-i ../ProcessedData/TEPAV_true_deletion.tsv.gz \
-o ../ProcessedData/flanking_allC_true_deletion.tsv \
-b 200 \
-n 20 \
-m ../RawData/all_mc_tables.txt

# CG
python get_mc_flanking.py \
-a [host] \
-p [password] \
-u [user] \
-d [database] \
-i ../ProcessedData/TEPAV_true_deletion.tsv.gz \
-o ../ProcessedData/flanking_CG_true_deletion.tsv \
-b 200 \
-n 20 \
-c CG \
-m ../RawData/all_mc_tables.txt

# CHG
python get_mc_flanking.py \
-a [host] \
-p [password] \
-u [user] \
-d [database] \
-i ../ProcessedData/TEPAV_true_deletion.tsv.gz \
-o ../ProcessedData/flanking_CHG_true_deletion.tsv \
-b 200 \
-n 20 \
-c CHG \
-m ../RawData/all_mc_tables.txt

# CHH
python get_mc_flanking.py \
-a [host] \
-p [password] \
-u [user] \
-d [database] \
-i ../ProcessedData/TEPAV_true_deletion.tsv.gz \
-o ../ProcessedData/flanking_CHH_true_deletion.tsv \
-b 200 \
-n 20 \
-c CHH \
-m ../RawData/all_mc_tables.txt

# Both types combined
# allC
python get_mc_flanking.py \
-a [host] \
-p [password] \
-u [user] \
-d [database] \
-i ../RawData/TEPID_TEPAV.tsv.gz \
-o ../ProcessedData/flanking_allC.tsv \
-b 200 \
-n 20 \
-m ../RawData/all_mc_tables.txt

#CG
python get_mc_flanking.py \
-a [host] \
-p [password] \
-u [user] \
-d [database] \
-i ../RawData/TEPID_TEPAV.tsv.gz \
-o ../ProcessedData/flanking_CG.tsv \
-b 200 \
-n 20 \
-c CG \
-m ../RawData/all_mc_tables.txt

#CHG
python get_mc_flanking.py \
-a [host] \
-p [password] \
-u [user] \
-d [database] \
-i ../RawData/TEPID_TEPAV.tsv.gz \
-o ../ProcessedData/flanking_CHG.tsv \
-b 200 \
-n 20 \
-c CHG \
-m ../RawData/all_mc_tables.txt

#CHH
python get_mc_flanking.py \
-a [host] \
-p [password] \
-u [user] \
-d [database] \
-i ../RawData/TEPID_TEPAV.tsv.gz \
-o ../ProcessedData/flanking_CHH.tsv \
-b 200 \
-n 20 \
-c CHH \
-m ../RawData/all_mc_tables.txt