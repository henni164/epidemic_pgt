#!/bin/bash -l
#SBATCH --job-name=process_vcfs
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=20GB
#SBATCH --time=1:00:00

rm -f all_varpos_index.txt

LIST=($(cut -f 1 allhaps.txt))

#gzip -d *.vcf.gz

for f in ${LIST[@]}; do sed -i 's/#CHROM/CHROM/g' ${f}.vcf; done

for f in ${LIST[@]}; do grep -v "#" ${f}.vcf | head -n 1 | cut -f 3,10-62 | sed 's/\t/\n/g' | awk 'FNR==NR{a[$1]=$2;next}{$1=a[$1]}1' haplookup.txt - | tr "\n" "\t" | sed 's/^\t/ID\t/g' > ${f}_setup.txt; done

for f in ${LIST[@]}; do grep -v "#" ${f}.vcf | grep -v "CHROM" | cut -f 3,10-62 >> ${f}_setup.txt; done

for f in ${LIST[@]}; do sed -i 's/\t>/\n>/g' ${f}_setup.txt; done

for f in ${LIST[@]}; do grep -v "#" ${f}.vcf | awk -v OFS="\t" -v HAP=$f '{print HAP,$1,$2,$3}' - >> all_varpos_index.txt; done
