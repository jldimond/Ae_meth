#!/bin/bash
## Job Name
#SBATCH --job-name=Ae_wtpoa
## Resources
## Nodes
#SBATCH --nodes=1
## Walltime (50 hours)
#SBATCH --time=50:00:00
## Memory per node
#SBATCH --mem=120G
##turn on e-mail notification
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jldimond@uw.edu
## Specify the working directory for this job
#SBATCH --workdir=/gscratch/scrubbed/jldimond/Ae_wtdbg2_103118

/gscratch/srlab/programs/wtdbg-2.1_x64-linux/wtpoa-cns -t 28 -i /gscratch/scrubbed/jldimond/Ae_wtdbg2_103118/Ae.ctg.lay -fo Ae.ctg.lay.fa