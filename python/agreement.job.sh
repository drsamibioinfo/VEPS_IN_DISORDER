#!/bin/sh
# Grid Engine options (lines prefixed with #$)
#$ -N agreement
#$ -cwd
#$ -l h_rt=48:00:00
#$ -l h_vmem=120G

# Initialise the environment modules
. /etc/profile.d/modules.sh

HD=/exports/igmm/eddie/marsh-lab/users/snouto
JD=$HD/secondpaper/agreement
OUT=$JD/outs
DATA=$HD/vectorized/all.variants.csv.gz
PEXEC=/home/s2273299/.conda/envs/features/bin/python

mkdir -p $OUT

$PEXEC $JD/agreement.py --mutations=$DATA --output=$OUT
