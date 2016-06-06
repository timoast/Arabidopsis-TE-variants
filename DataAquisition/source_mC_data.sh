# this file will download DNA methylation data for accessions used in Schmitz et al. 2013
# and create a MySQL database to hold the data
# The script must be edited before running to include MySQL login details

mkdir mC
cd mC
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE43nnn/GSE43857/suppl/GSE43857_RAW.tar
tar -x -f GSE43857_RAW.tar
gunzip *.gz
rm GSE43857_RAW.tar
python ../build_tables.py -u [username] -p [mysql_password] -d population_epigenetics --host [host]
cd ..
