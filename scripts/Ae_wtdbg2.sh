#!/bin/bash
## Job Name
#SBATCH --job-name=Ae_wtdbg2
## Resources
## Nodes
#SBATCH --nodes=1
## Walltime (360 hours)
#SBATCH --time=360:00:00
## Memory per node
#SBATCH --mem=128G
##turn on e-mail notification
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jldimond@uw.edu
## Specify the working directory for this job
#SBATCH --workdir=/gscratch/scrubbed/jldimond/Ae_wtdbg2_103118

/gscratch/srlab/programs/wtdbg-2.1_x64-linux/wtdbg2 -t 28 -S 3 -A -i /gscratch/srlab/data/Ae/merged2.fq -fo Ae
