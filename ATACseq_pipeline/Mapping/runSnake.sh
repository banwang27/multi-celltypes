#!/bin/bash
#SBATCH --job-name=snakemake
#SBATCH --time=7-00:00:00
##SBATCH --qos=long
#SBATCH -p hbfraser
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=astarr97@stanford.edu

snakemake --unlock

echo "---start running ----"

snakemake --keep-going --rerun-incomplete --configfile config.yaml -j 500 --latency-wait 60 --cluster "sbatch -p hbfraser --time {params.run_time} --cpus-per-task {params.cores} --mem {params.memory}000"

echo "---end running ----"