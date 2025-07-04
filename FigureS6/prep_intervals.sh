#!/bin/bash -l
#SBATCH --job-name=bigplot
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=2:00:00
#SBATCH --mem=50GB

#to submit: sbatch -a 0-6 --export INPUTLIST="./allhaps.txt" prep_intervals.sh

module load R

if [ -s "$INPUTLIST" ]
then
    SAMPLES=( $(cut -f 1 $INPUTLIST) );

    if [ ! -z "$SLURM_ARRAY_TASK_ID" ]
    then
        i=$SLURM_ARRAY_TASK_ID

	   Rscript prep_intervals.R -i ${SAMPLES[i]}

    else
        echo "Error: Missing array index as SLURM_ARRAY_TASK_ID"
    fi
else
    echo "Error: Missing inputs file list as --export env INPUTLIST or file empty"
fi
