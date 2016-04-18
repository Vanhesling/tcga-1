#!/bin/bash

#PBS -N czysz_tcga
#PBS -S /bin/bash
#PBS -l walltime=30:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=32gb

#PBS -o $HOME/tcga.out
#PBS -e $HOME/tcga.err

module load R/3.1.0

Rscript /home/t.cri.cczysz/tcga/mirna_tcga.R
