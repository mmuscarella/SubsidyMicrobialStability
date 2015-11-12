#!/bin/bash
#PBS -k o
#PBS -l nodes=1:ppn=8,vmem=100gb,walltime=24:00:00
#PBS -M mmuscare@indiana.edu,lennonj@indiana.edu
#PBS -m abe
#PBS -j oe
cd /N/dc2/projects/Lennon_Sequences/2015_PondsMicrobialStability
module load gcc
module load mothur/1.32.1 
mothur Ponds.Bacteria.Batch
