#!/bin/bash
cd /SSD/reg259/gene_drive_sims/20251021_figure2C-D
MODEL=/SSD/reg259/gene_drive_sims/WF_multilocus_cluster.slim

for COST in 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9;
do
    for REPEAT in {1..20};
    do
        slim -d drive_coeff=${COST} $MODEL > stdout_COST${COST}_rep${REPEAT}.txt
    done
done