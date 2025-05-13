module load samtools

module load minimap2

#run one Avr at a time
avr="AvrSr13"
#avr="AvrSr22"
#avr="AvrSr27"
#avr="AvrSr33"
#avr="AvrSr50_35"
#avr="AvrSr62"

LIST=($(cut -f 1 -d"_" ${avr}_intervals.txt))

for f in ${LIST[@]}; do

	interval=$(grep -E "${f}" ${avr}_intervals.txt | cut -f 2 -d"_")
	echo ${interval}
	samtools faidx ${f}.fasta ${interval} | sed "s/\(>chr[0-9]*\).*$/\1_${f}/g" > ${f}_${avr}_interval.fasta

done

cat hap*_${avr}_interval.fasta > ${avr}_intervals.fasta

minimap2 -X -N 1 -g 300 -m 20 -O7 -E7 ${avr}_intervals.fasta ${avr}_intervals.fasta > ${avr}_intervals.paf

sed 's/[_:-]/\t/g' ${avr}_intervals.txt | awk -v OFS="\t" '{print $2"_"$1,$3,$4}' > ${avr}_intervals_for_gggenomes.txt
