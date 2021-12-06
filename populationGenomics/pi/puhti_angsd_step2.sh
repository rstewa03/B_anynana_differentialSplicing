#!/bin/bash
#SBATCH --job-name=angsd_2
#SBATCH --account=project_2002170
#SBATCH --time=36:00:00
#SBATCH --mem-per-cpu=2G   # per core
#SBATCH --partition=small
#SBATCH --cpus-per-task=20   # max 40
#SBATCH --ntasks=1
#SBATCH --output=slurm_%j_%x.out

module load biokit
module load angsd

angsd_saf=/scratch/project_2002170/vicencio/anynana/angsd_output/Bany_NCBIref_genomic_angsd_step1a.saf.idx
angsd_sfs=/scratch/project_2002170/vicencio/anynana/angsd_output/Bany_NCBIref_genomic_angsd_step1b.sfs
angsd_output=/scratch/project_2002170/vicencio/anynana/angsd_output/Bany_NCBIref_genomic_angsd_step2

echo $(date)

realSFS saf2theta $angsd_saf -outname out -sfs $angsd_sfs -fold 1

echo $(date)






