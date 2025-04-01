#!/bin/sh
# Grid Engine options (lines prefixed with #$)
#$ -N METRICS
#$ -cwd
#$ -l h_rt=12:00:00
#$ -l h_vmem=128G

# Initialise the environment modules
. /etc/profile.d/modules.sh

HD=/exports/igmm/eddie/marsh-lab/users/snouto
JD=$HD/genesis/metrics_latest
DATA=$HD/vectorized/all.variants.csv.gz
PEXEC=/home/s2273299/.conda/envs/features/bin/python


$PEXEC $JD/metrics.py --missense=$DATA --output=$JD
