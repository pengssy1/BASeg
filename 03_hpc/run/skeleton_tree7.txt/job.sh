#!/bin/bash
#PBS -N skeleton_tree7.txt
#PBS -l nodes=1:ppn=1
#PBS -l walltime=72:00:00
#PBS -l mem=32gb

ml purge; ml R-bundle-Bioconductor/3.15-foss-2021b-R-4.2.0
# ml purge; ml R-bundle-Bioconductor/3.15-foss-2021b-R-4.2.0
cd $PBS_O_WORKDIR
# Rscript run.r


Rscript /data/gent/vo/000/gvo00074/Xipeng/creat/run/skeleton_tree7.txt/run.r