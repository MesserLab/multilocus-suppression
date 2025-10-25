#!/bin/bash
cd /SSD/reg259/gene_drive_sims/20251021_figure2A-B
MODEL=/SSD/reg259/gene_drive_sims/WF_multilocus_cluster.slim

for R1 in 1e-2 1e-3 1e-4 1e-5 1e-6;
do
    for REPEAT in {1..20};
    do
        slim -d resistance_rate=${R1} $MODEL > stdout_R1${R1}_rep${REPEAT}.txt
    done
done