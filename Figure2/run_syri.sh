#!/bin/bash -l
#SBATCH --job-name=syri
#SBATCH --nodes=1
#SBATCH --ntasks=18
#SBATCH --mem=20GB
#SBATCH --time=1:00:00

#to submit: sbatch -a 0-5 --export INPUTLIST="./coords_files.txt" run_syri.sh

module load syri
module load mummer

if [ -s "$INPUTLIST" ]
then
    SAMPLES=( $(cut -f 1 $INPUTLIST) );

    if [ ! -z "$SLURM_ARRAY_TASK_ID" ]
    then
        i=$SLURM_ARRAY_TASK_ID

        hap1=($(echo ${SAMPLES[i]} | cut -f 1 -d"v" | sed 's/_$//g'))
        hap2=($(echo ${SAMPLES[i]} | cut -f 2 -d"v" | sed 's/^s_//g' | cut -f 1 -d"."))

		syri -c ${SAMPLES[i]} -r ${hap1}.fasta -q ${hap2}.fasta -d ${hap1}_vs_${hap2}.filtered.delta -F T --prefix ${hap1}_vs_${hap2} --nc 18 -s /apps/mummer/4.0.0rc1/bin/show-snps

    else
        echo "Error: Missing array index as SLURM_ARRAY_TASK_ID"
    fi
else
    echo "Error: Missing inputs file list as --export env INPUTLIST or file empty"
fi
