#!/bin/bash
#SBATCH --job-name=snakemake
#SBATCH --error=job_errors/snakemake.%j.err
#SBATCH --output=job_errors/snakemake.%j.out
#SBATCH --time=7-00:00:00
##SBATCH --qos=long
#SBATCH -p hbfraser
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=banwang3@stanford.edu

echo "---start running ----"
ml contribs
ml poldrack
ml anaconda
source activate humanzee
export PYTHONPATH="/home/users/banwang3/.conda/envs/humanzee/lib/python3.6/site-packages"

snakemake --keep-going --configfile config.yaml -j 500 --latency-wait 60 --cluster "sbatch -p hbfraser --job-name {params.job_name} -o {params.job_out_dir}/{params.job_out_file}.out -e {params.job_out_dir}/{params.job_out_file}.error --time {params.run_time} --cpus-per-task {params.cores} --mem {params.memory}000"

echo "---end running ----"
