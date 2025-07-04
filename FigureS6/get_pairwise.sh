LIST=($( cut -f 1 haplist.txt))
haplengths=($(cut -f 2 haplengths.txt))

for i in ${!LIST[@]}; do
	for f in ${LIST[@]}; do 
		if [ ${f} == ${LIST[${i}]} ]; then
			:
		else
			awk -v OFS="\t" -v mylength="${haplengths[${i}]}" '{print (($3 - $2) / mylength) * 100}' ./forSfig/${f}_${LIST[${i}]}_shared.txt | awk -v OFS="\t" -v hap1="${f}" -v hap2="${LIST[${i}]}" '{s+=$1} END {print s,hap1,hap2}'
		fi
	done
done
