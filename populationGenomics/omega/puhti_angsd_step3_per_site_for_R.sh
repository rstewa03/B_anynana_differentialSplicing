#!/bin/bash
#SBATCH --job-name=angsd_3c
#SBATCH --account=project_2002170
#SBATCH --time=36:00:00
#SBATCH --mem-per-cpu=2G   # per core
#SBATCH --partition=small
#SBATCH --cpus-per-task=20   # max 40
#SBATCH --ntasks=1
#SBATCH --output=slurm_%j_%x.out

module load biokit
module load angsd


angsd_thetas_per_site=/scratch/project_2002170/vicencio/anynana/angsd_output/Bany_NCBIref_genomic_angsd_step2.thetas.idx
angsd_thetas_per_site_out=/scratch/project_2002170/vicencio/anynana/angsd_output/Bany_NCBIref_genomic_angsd_step2.thetas_per_site.gz

echo $(date)

thetaStat print $angsd_thetas_per_site | gzip > $angsd_thetas_per_site_out

echo $(date)



