#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=2GB
#SBATCH --time=24:00:00
#SBATCH --account=knayak_127

module load matlab/2019a
matlab -nosplash -nodesktop -r demo_25