# This file will download PE sequencing data for all Arabidopsis accessions
# in Schmitz et al. 2013, and extract the fastq files

cd ../RawData
while read acc; do
  echo $acc
  wget -r --no-parent --no-directories "ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR492/${acc}/*"
done < ./SraAccList.txt

for myfile in ./*.sra; do
  fname=(${myfile//.sra/ })
  mkdir $fname
  mv $myfile $fname
done

for directory in ./*; do
    if [ -d "$directory" ]; then
        cd $directory
        for myfile in $(ls -d *.sra);do
            fastq-dump --split-3 -v $myfile
            rm -f $myfile
        done
        cd ..
    fi
done
cd ../DataAquisition