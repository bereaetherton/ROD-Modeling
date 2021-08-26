#!/bin/sh
#SBATCH --account=epi                # what account you are with
#SBATCH --qos=epi-b                  # use which account 
#SBATCH --job-name=RODanaIVP            # Job name
#SBATCH --mail-type=ALL               # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=betherton@ufl.edu           # Where to send mail	
#SBATCH --nodes=1                   # Use one node
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=12             # Use 1 core
#SBATCH --mem=3GB                   # Memory limit
#SBATCH --time=2:00:00               # Time limit hrs:min:sec
#SBATCH --output=RODanaIVP%j.out   # Standard output and error log
#SBATCH --array=1-144%100

pwd;hostname;date

module load ufrc
module load R/3.6

ARG=$(sed -n ${SLURM_ARRAY_TASK_ID}p IVPparams.txt)
echo ${ARG} 
Rscript --vanilla ana.for.hpg.ivp.R ${ARG}
date
