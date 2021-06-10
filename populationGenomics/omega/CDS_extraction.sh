# TOOLS=/path/to/software
# Bicyclus_anynana=/path/to/project

cd $Bicyclus_anynana/SCO_Parargae_aegeria

# use the OrthoVenn2 SCO table output
scp protein_IDs.Pararge_Bicyclus.tsv $Bicyclus_anynana/SCO_Parargae_aegeria


head protein_IDs.Pararge_Bicyclus.tsv
jg12524.t1	LOC112042788
jg13290.t1	LOC112042789
jg13299.t1	LOC112042790
jg26874.t1	LOC112042791

# data
mkdir SCO_per_species_sets
cd SCO_per_species_sets

# CDS sequence sets
# Paegeria
ln -s $Bicyclus_anynana/Parargae_aegeria/Pararge_aegeria_v2.fa.longest_isoform.cds.fa .
# Bany
ln -s $Bicyclus_anynana/AS_genomic_data/GCF_900239965.1_Bicyclus_anynana_v1.2_genomic.fna.longest_isoform.cds.fa .
sed 's/rna.*.gene=//g' GCF_900239965.1_Bicyclus_anynana_v1.2_genomic.fna.longest_isoform.cds.fa > GCF_900239965.1_Bicyclus_anynana_v1.2_genomic.fna.longest_isoform.cds.LOC.fa

# index CDS sets
cat Pararge_aegeria_v2.fa.longest_isoform.cds.fa GCF_900239965.1_Bicyclus_anynana_v1.2_genomic.fna.longest_isoform.cds.LOC.fa > CDS_combined.fa
samtools faidx CDS_combined.fa

# grab CDS from 2 species

# search and extract
# because echo changes tab to space
mkdir SCO_combined_fastas
while read p; do set=$(echo $p); Bid=$(echo $set| cut -f2 -d ' ' ); samtools faidx SCO_per_species_sets/CDS_combined.fa $set -o SCO_combined_fastas/$Bid.fasta ; done < protein_IDs.Pararge_Bicyclus.tsv

# assess
ls SCO_combined_fastas/ | wc -l
9362

grep '>' SCO_combined_fastas/*fasta | wc -l
18724

# get lengths
for p in *.fasta; do $TOOLS/exonerate-2.2.0-x86_64/bin/fastalength $p; done
