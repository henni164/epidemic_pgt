#!/bin/bash -l
#SBATCH --job-name=mummer
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem=40GB
#SBATCH --time=24:00:00

#sbatch -a 0-5 --export INPUTLIST="./pgthaps.txt" get_mummer.sh

module load mummer

if [ -s "$INPUTLIST" ]
then
    SAMPLES=( $(cut -f 1 $INPUTLIST) );

    if [ ! -z "$SLURM_ARRAY_TASK_ID" ]
    then
        i=$SLURM_ARRAY_TASK_ID

        let j=$i+1

        echo "first ref is ${SAMPLES[i]}"
        echo "ref1 is index ${i}"
        echo "ref2 is index ${j}"
        hap1=${SAMPLES[i]}
        hap2=${SAMPLES[j]}
        echo "doing ${hap1} versus ${hap2} now"
	   nucmer -t 16 --maxmatch -c 100 -b 500 -l 50 -p ${hap1}_vs_${hap2} ${hap1}.fasta ${hap2}.fasta
	   delta-filter -m -i 60 -l 1000 ${hap1}_vs_${hap2}.delta > ${hap1}_vs_${hap2}.filtered.delta
	   show-coords -THrd ${hap1}_vs_${hap2}.filtered.delta > ${hap1}_vs_${hap2}.filtered.coords

    else
        echo "Error: Missing array index as SLURM_ARRAY_TASK_ID"
    fi
else
    echo "Error: Missing inputs file list as --export env INPUTLIST or file empty"
fi
