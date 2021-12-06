#!/bin/bash
#SBATCH --job-name=angsd_1a
#SBATCH --account=project_2002170
#SBATCH --time=36:00:00
#SBATCH --mem-per-cpu=2G   # per core
#SBATCH --partition=small
#SBATCH --cpus-per-task=20   # max 40
#SBATCH --ntasks=1
#SBATCH --output=slurm_%j_%x.out

module load biokit
module load angsd

bamlist=bamlist.txt
refgenome=GCF_900239965.1_Bicyclus_anynana_v1.2_genomic.fna
angsd_output=angsd_output/Bany_NCBIref_genomic_angsd_step1a

echo $(date)

srun angsd -P $SLURM_CPUS_PER_TASK -b $bamlist \
	-ref $refgenome \
	-anc $refgenome \
	-out $angsd_output \
	-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 -minQ 13 -minMapQ 20 \
	-minInd 5 \
	-GL 1 -doSaf 1

echo $(date)






