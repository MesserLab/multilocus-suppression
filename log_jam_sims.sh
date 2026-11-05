#!/usr/bin/env bash
# Run N_RUNS copies of the same simulation at once. All runs save into
# DATA_PATH; the haplotype printout of each run goes to log-jam-sim1.txt,
# log-jam-sim2.txt, ... Numbers already in use are skipped (never overwritten).
set -euo pipefail
cd "$(dirname "$0")"

# define parameters
N_RUNS=9
DATA_PATH=data/log-jam-sims
RECOMB=0.0000001
LOCI=4
COST=0.1
MODEL=WF_multilocus_cluster.slim

# run the model
mkdir -p $DATA_PATH
k=0
started=0
while (( started < N_RUNS )); do
    k=$((k + 1))
    OUT=${DATA_PATH}/log-jam-sim${k}
    [[ -e ${OUT}.txt ]] && continue
    setsid nohup slim -d "path='${DATA_PATH}/'" \
        -d recomb_rate=${RECOMB} \
        -d num_target_loci=${LOCI} \
        -d broken_coeff=${COST} \
        -d PRINT_COMMON_HAPLOTYPES=T \
        $MODEL > ${OUT}.txt 2>&1 < /dev/null &
    echo "started ${OUT}.txt (pid $!)"
    started=$((started + 1))
    sleep 2   # SLiM names its output by date_time; keep the runs' names distinct
done
