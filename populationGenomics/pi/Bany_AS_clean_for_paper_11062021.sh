##### Whole-genome sequencing data done at Centre for Genomic Research (CGR) University of Liverpool 2020
##### Nucleotide diversity analysis was performed by Vicencio Oostra (oostra@liverpool.ac.uk) 2021

### Data info from manuscript: 
## "We generated whole-genome DNA resequencing data from a wild B. anynana population from Zomba, Malawi (15°22’S, 35°19’E). The Zomba population was collected and brought to 
## the laboratory in March 2007 by Maaike de Jong and has been studied previously in phylogeographic analyses using candidate genes (Jong et al. 2011; de Jong et al. 2013).
## We selected five females for whole-genome sequencing. 
## Illumina fragment libraries were prepared by the Centre for Genomic Research (University of Liverpool) using the NEBNext Ultra II FS kit (New England BioLabs) on the Mosquito 
## platform with 1/10 volumes (SPT Labtech), and quantified using qPCR. Libraries were sequenced over ¼ of a lane on the Illumina NovaSeq using S4 chemistry, with 2x150 bp paired-end reads. 
## Sequencing depth after trimming ranged 53-96 million reads (mean: 77M)."


### software and versions
# gffread-0.9.12.Linux_x86_64
# bwa 0.7.17-r1198-dirty
# samtools 1.10-13-ga2916aa
# picard-tools-2.0.1.zip
# angsd version 0.933 build 15 Sep 2020


### Trimming of raw reads was done by bioinformatics support at CGR as follows
# "The raw Fastq files are trimmed for the presence of Illumina adapter sequences using Cutadapt version 1.2.1. The option -O 3 was used, so the 3' end of any reads which match the adapter sequence for 3 bp. or more are trimmed."
# "The reads are further trimmed using Sickle version 1.200 with a minimum window quality score of 20. Reads shorter than 15 bp. after trimming were removed."


### reference genome assembly and annotation: Bicyclus anynana v1.2 (abbr. as v1) from NCBI: https://www.ncbi.nlm.nih.gov/genome/?term=bicyclus+anynana
## genome file:
GCF_900239965.1_Bicyclus_anynana_v1.2_genomic.fna
## annotation to calculate per-gene polymorphisms for coding and full gene sequence
## "The gff (accessed from NCBI PRJNA434100 19 Dec 2019) was converted to a gtf using gffread. The gene IDs that were used throughout the analyses were the LOC IDs (e.g. LOC112042788). "
## filtered to only have longest isoform per locus
GCF_900239965.1_Bicyclus_anynana_v1.2_genomic.longest_isoform.gff

## produce coding sequence from gff and genome
# files
genome=GCF_900239965.1_Bicyclus_anynana_v1.2_genomic.fna
gff=GCF_900239965.1_Bicyclus_anynana_v1.2_genomic.longest_isoform.gff
output_cds=GCF_900239965.1_Bicyclus_anynana_v1.2_genomic.longest_isoform_cds.fa
# output coding sequence;
gffread -x $output_cds -g $genome $gff
grep ">" GCF_900239965.1_Bicyclus_anynana_v1.2_genomic.longest_isoform_cds.fa | wc -l
cut -f3 GCF_900239965.1_Bicyclus_anynana_v1.2_genomic.longest_isoform.gff | grep mRNA | wc -l

### Method: ANGSD
## map reads to v1 reference
## use angsd to calculate thetas (pi, watterson's theta) for coding sequence


## map reads to Bany v1 reference with bwa; 
# index reference for bwa
ncbigenome=GCF_900239965.1_Bicyclus_anynana_v1.2_genomic.fna
bwa index $ncbigenome
# initial mapping
ncbigenome=GCF_900239965.1_Bicyclus_anynana_v1.2_genomic.fna
ref_genome=$ncbigenome
readsdir=Trimmed/
bamdir=mapped_reads/NCBI_genome_v1.2_bwa/
while read lib ; do
    echo $lib
    read1=$readsdir$lib"_L002_R1_001.fastq.gz"
    read2=$readsdir$lib"_L002_R2_001.fastq.gz"
    outputbam=$bamdir$lib".bam"
    ls $read1
    ls $read2
    echo $outputbam
    bwa mem -t 7 $ref_genome $read1 $read2 | samtools view -bS - > $outputbam
done < libraries_anynana.txt

## sort, remove duplicates, index in one go
bamdir=mapped_reads/NCBI_genome_v1.2_bwa/
while read lib ; do
    echo $lib
    inputbam=$bamdir$lib".bam"
    sortedbam=$bamdir$lib"_sorted.bam"
    dedupbam=$bamdir$lib"_sorted_dedupped.bam"
    outputmatrixfile=$bamdir$lib"_sorted_dedupped_matrix.txt"
    ls $inputbam
    echo $(date)
    samtools sort -@ 7 $inputbam > $sortedbam
    java -Xms4g -Xmx30g -jar picard.jar MarkDuplicates INPUT=$sortedbam \
	    OUTPUT=$dedupbam \
	    METRICS_FILE=$outputmatrixfile \
	    REMOVE_DUPLICATES=true 
    samtools index -b $dedupbam
done < libraries_anynana.txt

# add read groups
bamdir=mapped_reads/NCBI_genome_v1.2_bwa/
while read lib ; do
    echo $lib
    dedupbam=$bamdir$lib"_sorted_dedupped.bam"
    outputbam=$bamdir$lib"_sorted_dedupped_addRG.bam"
    RG=$(echo $lib | cut -f1-3 -d'_')
    ls $dedupbam
    echo $RG
    echo $(date)
    java -Xms4g -Xmx30g -jar picard.jar AddOrReplaceReadGroups INPUT=$dedupbam OUTPUT=$outputbam \
         RGID=anynana RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=$RG
done < libraries_anynana.txt
# and index
bamdir=mapped_reads/NCBI_genome_v1.2_bwa/
while read lib ; do
    echo $lib
    outputbam=$bamdir$lib"_sorted_dedupped_addRG.bam"
    ls $outputbam
    samtools index $outputbam
done < libraries_anynana.txt




#### run angsd to get thetas
## 1) global estimate of SFS using priors across whole genome
## 2) calculate per-site thetas
## 3) estimate thetas in coding sequence of genes
## 4) combine thetas across different regions (exons) of same gene, using R script

### All angsd commands run on cluster using slurm scripts

## step 1: global estimate of SFS using whole genome

# 1a command:
bamlist=bamlist.txt
refgenome=GCF_900239965.1_Bicyclus_anynana_v1.2_genomic.fna
angsd_output=angsd_output/Bany_NCBIref_genomic_angsd_step1a
angsd -b $bamlist \
	-ref $refgenome \
	-anc $refgenome \
	-out $angsd_output \
	-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 -minQ 13 -minMapQ 20 \
	-minInd 5 \
	-GL 1 -doSaf 1

# 1b: Obtain the maximum likelihood estimate of the SFS using the realSFS program
angsd_input=angsd_output/Bany_NCBIref_genomic_angsd_step1a.saf.idx
angsd_output=angsd_output/Bany_NCBIref_genomic_angsd_step1b
realSFS $angsd_input  -fold 1 > $angsd_output

## step 2: Calculate per-site thetas.
angsd_saf=angsd_output/Bany_NCBIref_genomic_angsd_step1a.saf.idx
angsd_sfs=angsd_output/Bany_NCBIref_genomic_angsd_step1b.sfs
angsd_output=angsd_output/Bany_NCBIref_genomic_angsd_step2
realSFS saf2theta $angsd_saf -outname out -sfs $angsd_sfs -fold 1


## step 3: Calculate neutrality tests statistics
angsd_thetas_per_site=angsd_output/Bany_NCBIref_genomic_angsd_step2.thetas.idx
angsd_thetas_per_site_out=angsd_output/Bany_NCBIref_genomic_angsd_step2.thetas_per_site.gz
thetaStat print $angsd_thetas_per_site | gzip > $angsd_thetas_per_site_out


### step 4: calculate thetas per exon, then combine those across different regions (exons) of same gene
## first extract locations of exons, then use those in R script to calculate thetas per gene manually in R from per-site thetas

## make bed file per gene from gff of just the coding sequence, to use in R script
gff=GCF_900239965.1_Bicyclus_anynana_v1.2_genomic.longest_isoform.gff
cut -f1-5 $gff | grep -n "\bCDS\b" | cut -f1 -d':' > CDS_lines.txt
while read line  ; do
    head -n $line $gff | tail -1 >> CDS.gff
done < CDS_lines.txt
cut -f1,4,5 CDS.gff > CDS.bed
cut -f9 CDS.gff | awk -F 'gene=' '{print $2}' | awk -F ';' '{print $1}' > CDS_LOC_names_corrected.txt
paste CDS.bed CDS_LOC_names_corrected.txt > CDS_bed_and_names_for_R.txt

## make list of scaffolds form genome index
genome=GCF_900239965.1_Bicyclus_anynana_v1.2_genomic.fna.fai
cut -f1 $genome > GCF_900239965.1_Bicyclus_anynana_v1.2_genomic_fna_scaffolds_list.txt

## loop over scaffolds, and for each scaffold a) subset per-site thetas angsd file; b) call R script which calculates per-gene-thetas from per-site-thetas
scaffolds=Bany_genome_v1.2_from_NCBI/GCF_900239965.1_Bicyclus_anynana_v1.2_genomic_fna_scaffolds_list.txt
thetas=Bany_NCBIref_genomic_angsd_step2.thetas_per_site
angsd_temp_folder=angsd_per_scaffold_temp_folder/
CDS=CDS_bed_and_names_for_R.txt
while read chrom ; do
  echo scaffold $chrom
  # subset thetas for current scaffold into separate file (still per-site) which is used by the R script
  grep $chrom $thetas > $angsd_temp_folder$chrom"_angsd_thetas.txt"
  # calculate thetas for coding sequence in R
  Rscript --vanilla calculate_per-gene_thetas_from_angsd_per-site_thetas_headless_CDS.R $chrom
done < $scaffolds



## post-processing:
# 1. paste together all per-CDS thetas which are now in separate files per scaffold;
# 15642 whole genes; 103293 exons
cd angsd_per_scaffold_temp_folder/
files=angsd_thetas_CDS*
outfile=angsd_thetas_genes_CDS_exons.txt
head -n1 angsd_thetas_CDS_NW_019862752.1.txt > $outfile  # header
for f in $files; do
    tail -n +2 $f >> $outfile
done
# 2. calculate average theta and pi per CDS gene using exon-level data (use length-weighted average)
# do in R, see calculate_per-gene_thetas_from_angsd_per-exon_thetas.R. output file "angsd_thetas_genes_CDS.txt"
