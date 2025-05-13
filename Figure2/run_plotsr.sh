#!/bin/bash -l
#SBATCH --job-name=plotsr
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=20GB
#SBATCH --time=1:00:00

module load plotsr

plotsr --itx --chrord chrord.txt --sr hap02_vs_hap01syri.out \
--sr hap01_vs_hap03syri.out \
--sr hap03_vs_hap04syri.out \
--sr hap04_vs_hap05syri.out \
--sr hap05_vs_hap06syri.out \
--sr hap06_vs_hap07syri.out \
--genomes genomes.txt --cfg config.cfg -s 50000 \
-o Figure2B_synteny.png -W 11 -H 4 -f 6 -d 600 -S 1.5
