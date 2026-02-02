#!/bin/sh
# Grid Engine options (lines prefixed with #$)
#$ -N WHOLEIDR
#$ -cwd
#$ -l h_rt=72:00:00
#$ -l h_vmem=32G

# Initialise the environment modules
. /etc/profile.d/modules.sh

HD=/exports/igmm/eddie/marsh-lab/users/snouto
SD=/exports/eddie/scratch/s2273299
JD=$HD/wholeidr
OUTFOLDER=$JD/proteins
PROTEINS=$JD/proteins.txt
FASTA=$JD/uniprots.seq.2024.fasta
ALPHAFOLD=/exports/igmm/eddie/marsh-lab/data/human_alphafold/alphafold_human_residues.out
PEXEC=/home/s2273299/.conda/envs/features/bin/python

mkdir -p  $OUTFOLDER

$PEXEC $JD/annotate.location.py --fasta=$FASTA --alphafold=$ALPHAFOLD --proteins=$PROTEINS --position=$SGE_TASK_ID --output=$OUTFOLDER
