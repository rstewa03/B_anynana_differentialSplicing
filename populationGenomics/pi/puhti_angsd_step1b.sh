#!/bin/bash
#SBATCH --job-name=angsd_1b
#SBATCH --account=project_2002170
#SBATCH --time=36:00:00
#SBATCH --mem-per-cpu=8G   # per core
#SBATCH --partition=small
#SBATCH --cpus-per-task=20   # max 40
#SBATCH --ntasks=1
#SBATCH --output=slurm_%j_%x.out

module load biokit
module load angsd


echo $(date)

angsd_input=/scratch/project_2002170/vicencio/anynana/angsd_output/Bany_NCBIref_genomic_angsd_step1a.saf.idx
angsd_output=/scratch/project_2002170/vicencio/anynana/angsd_output/Bany_NCBIref_genomic_angsd_step1b
srun realSFS $angsd_input -P $SLURM_CPUS_PER_TASK  -fold 1 > $angsd_output      # command works

echo $(date)






