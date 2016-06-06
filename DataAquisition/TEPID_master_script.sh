# This script will run TEPID on all PE reads from Schmitz et al. 2013
# and identify TE insertion and TE absences with repect to the reference genome

cd ../RawData

# download reference genome
wget ftp://igenome:G3nom3s4u@ussd-ftp.illumina.com/Arabidopsis_thaliana/NCBI/TAIR10/Arabidopsis_thaliana_NCBI_TAIR10.tar.gz
tar -xvzf Arabidopsis_thaliana_NCBI_TAIR10.tar.gz

# build reference index
bowtie2-build TAIR10/Sequence/WholGenomeFasta/genome.fa tair10
yaha -g TAIR10/Sequence/WholGenomeFasta/genome.fa

# using one core for mapping
tepid-map -x tair10 -y ./TAIR10/Sequence/WholeGenomeFasta/genome.X15_01_65525S -p 1 -s 200 -r -z


discover () {
  for myfile in $(ls -d *.split.bam);do
    fname=(${myfile//.split.bam/ })
    tepy-discover -n $fname -c "${fname}.bam" -s $myfile -t ../../../TAIR9_TE.bed.gz
  done
}

# discover TEs
for directory in ./*; do
    if [ -d "$directory" ]; then
        cd $directory
        discover
        cd ..
    fi
done

# merge all calls
python merge_deletions.py -f deletions
python merge_insertions.py -f insertions
python create_poly_te.py

# run refine step
refine () {
  for myfile in $(ls -d *.split.bam);do
    fname=(${myfile//.split.bam/ })
    tepy-refine -a ../../../accession_names.tsv -i ../../insertions.bed -d ../../deletions.bed -n $fname -c "${fname}.bam" -s $myfile -t ../../../TAIR9_TE.bed.gz
  done
}

for directory in ./*; do
    if [ -d "$directory" ]; then
        cd $directory
        map
        cd ..
    fi
done

python genotype.py -i -a ambiguous -m merged_insertions.