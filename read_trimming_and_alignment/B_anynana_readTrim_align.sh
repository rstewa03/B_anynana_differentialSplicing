############################
#### B_anynana pipeline ####
############################

export BA_TOOLS=mypath/programs
export RNA_TOOLS=/cerberus/projects/racste/HISAT2_Tutorial/student_tools

export BA_HOME=mypath/B_anynana
export BA_DATA_DIR=$BA_HOME/RNA_data/raw
export BA_DATA_TRIM=$BA_HOME/RNA_data/trim
export BA_REFS_DIR=$BA_HOME/refs
export BA_REF_INDEX=$BA_REFS_DIR/GCF_900239965.1_Bicyclus_anynana_v1.2_genomic
export BA_REF_FASTA=$BA_REF_INDEX.fa
# /cerberus/projects/racste/B_anynana/GCF_900239965.1_Bicyclus_anynana_v1.2_genomic.fna
export BA_REF_GTF=$BA_REF_INDEX.gtf
# /cerberus/projects/racste/B_anynana/refs/GCF_900239965.1_Bicyclus_anynana_v1.2_genomic.gtf
export BA_ALIGN_DIR=$BA_HOME/alignments/star



########################################
#### Quality trim raw RNA-seq files ####
########################################
# clean and trim the reads within the folder containing the raw read files
cd $BA_DATA_DIR

# the script we will be using
# /data/programs/scripts/getlog_noclonefilter_Q20.py

# getting list of all the files to clean
ls *gz > fastq.gz_raw_all

# run the cleaning script on this subset of files.
python /data/programs/scripts/getlog_noclonefilter_Q20.py fastq.gz_raw_all
# log.txt file example output:
# Phred-encoding
# ==============
#   sanger  fastq   gz      single-ended    100bp
#
# Adapter filter
# ==============
#   java -Djava.library.path=/data/programs/bbmap_35.69/jni/ -ea -Xmx182201m -Xms182201m -cp /data/programs/bbmap_35.69/current/ jgi.BBDukF in=SRR5430602_1.fastq.gz in2=SRR5430602_2.fastq.gz out=SRR5430602_1.fastq.gztemp_1.fq.gz out2=SRR5430602_2.fastq.gztemp_2.fq.gz ref=/data/programs/bbmap_35.69/resources/truseq.fa.gz,/data/programs/bbmap_35.69/resources/nextera.fa.gz ktrim=r k=23 mink=11 hdist=1 overwrite=t
# Executing jgi.BBDukF [in=SRR5430602_1.fastq.gz, in2=SRR5430602_2.fastq.gz, out=SRR5430602_1.fastq.gztemp_1.fq.gz, out2=SRR5430602_2.fastq.gztemp_2.fq.gz, ref=/data/programs/bbmap_35.69/resources/truseq.fa.gz,/data/programs/bbmap_35.69/resources/nextera.fa.gz, ktrim=r, k=23, mink=11, hdist=1, overwrite=t]
#
# BBDuk version 35.69
# maskMiddle was disabled because useShortKmers=true
# Initial:
#   Memory: max=191058m, free=190975m, used=83m
#
# Added 140059 kmers; time:       0.185 seconds.
# Memory: max=191058m, free=190807m, used=251m
#
# Input is being processed as paired
# Started output streams: 0.046 seconds.
# Processing time:                64.982 seconds.
#
# Input:                          30487822 reads          3048782200 bases.
# KTrimmed:                       36759 reads (0.12%)     409898 bases (0.01%)
# Result:                         30487822 reads (100.00%)        3048372302 bases (99.99%)
#
# Time:                           65.239 seconds.
# Reads Processed:      30487k    467.33k reads/sec
# Bases Processed:       3048m    46.73m bases/sec
#
# Low quality filter
# ==================
#   java -Djava.library.path=/data/programs/bbmap_35.69/jni/ -ea -Xmx182150m -Xms182150m -cp /data/programs/bbmap_35.69/current/ jgi.BBDuk2 in=SRR5430602_1.fastq.gztemp_1.fq.gz in2=SRR5430602_2.fastq.gztemp_2.fq.gz out=SRR5430602_1.fastq.filt.fq.gz out2=SRR5430602_2.fastq.filt.fq.gz ref=/data/programs/bbmap_35.69//resources/phix174_ill.ref.fa.gz k=27 hdist=1 qtrim=rl trimq=20 minlen=40 qout=33 overwrite=t
# Executing jgi.BBDuk2 [in=SRR5430602_1.fastq.gztemp_1.fq.gz, in2=SRR5430602_2.fastq.gztemp_2.fq.gz, out=SRR5430602_1.fastq.filt.fq.gz, out2=SRR5430602_2.fastq.filt.fq.gz, ref=/data/programs/bbmap_35.69//resources/phix174_ill.ref.fa.gz, k=27, hdist=1, qtrim=rl, trimq=20, minlen=40, qout=33, overwrite=t]
#
# BBDuk2 version 35.69
# k=27
# maskMiddle=true
# hamming distance=1
# kfiltering using 1 reference.
# quality-trimming both ends to Q20
#
# Initial:
#   Memory: max=191025m, free=190941m, used=84m
#
# Added 423440 kmers; time:       0.218 seconds.
# Memory: max=191025m, free=190790m, used=235m
#
# Input is being processed as paired
# Started output streams: 1.564 seconds.
# Processing time:                46.121 seconds.
#
# Input:                          30487822 reads          3048372302 bases.
# Contaminants:                   0 reads (0.00%)         0 bases (0.00%)
# QTrimmed:                       8541561 reads (28.02%)  319131323 bases (10.47%)
# Result:                         28604510 reads (93.82%)         2729240979 bases (89.53%)
#
# Time:                           47.944 seconds.
# Reads Processed:      30487k    635.90k reads/sec
# Bases Processed:       3048m    63.58m bases/sec


# Convert the log.txt file to a table (modified from Wheat via Dort)
grep -E '_1.fastq.gztemp_1.fq.gz out2=|Input:|Result:' log.txt | cut -f3 -d '=' |  cut -f1 -d "_" | tr -d '\r\n'| tr -s " " |\
#remove tabs, make all lines begin with SRR files, remove blank lines
awk '{gsub("\t","",$0); print;}' | sed 's/SRR/\nSRR/g' | grep "\S" |\
#the following line replaces filler text with commas so we have a usable csv. Spaces matter.
sed -e 's/ reads /,/g' -e 's/Input: /,/g' -e 's/ reads /,/g' -e 's/ bases.Result: /,/g' -e 's/ reads /,/g' -e 's/ bases /,/g' -e 's/Input: /,/g'  -e 's/ reads /,/g' -e 's/ bases.Result: /,/g' -e 's/) /,/g' -e 's/ reads /,/g' -e 's/) bases //g' -e 's/(//g' -e 's/)//g' |\
sed '1 i file,adapter_in_reads,adapter_in_bases,adapter_out_reads,adapter_out_pct_reads,adapter_out_bases,adapter_out_bases_pct,lowq_in_reads,lowq_in_bases,lowq_out_reads,lowq_out_reads_pct,lowq_out_bases,lowq_out_bases_pct' > filt_qtrim_log.csv

##################################################
#### STAR Index of the genome and gtf files ####
##################################################
cd $BA_ALIGN_DIR

mkdir -p $BA_HOME/alignments/star/genome_index_STAR_r1
$BA_TOOLS/STAR/bin/Linux_x86_64/STAR --runThreadN 4 --runMode genomeGenerate --genomeDir genome_index_STAR_r1 --genomeFastaFiles /cerberus/projects/racste/B_anynana/refs/hisat_index/GCF_900239965.1_Bicyclus_anynana_v1.2_genomic.fna --sjdbGTFfile /cerberus/projects/racste/B_anynana/refs/hisat_index/GCF_900239965.1_Bicyclus_anynana_v1.2_genomic.gtf --sjdbOverhang 99

#########################################
#### Alignment with STAR: First Pass ####
#########################################
# example information# /data/programs/STAR/bin/Linux_x86_64/STAR --runThreadN 8 --genomeDir /cerberus/projects/racste/B_anynana/alignments/star/genome_index_STAR_r1 --outFileNamePrefix SRR5430602 --readFilesIn /cerberus/projects/racste/B_anynana/RNA_data/raw/SRR5430602_1.fastq.filt.fq.gz /cerberus/projects/racste/B_anynana/RNA_data/raw/SRR5430602_2.fastq.gz --outSJfilterReads Unique

# move to the raw data directory to
cd $BA_DATA_TRIM
rm $BA_ALIGN_DIR/jobs.txt
for I in `ls -1 *_1.fastq.filt.fq.gz | sed 's/_1.fastq.filt.fq.gz//'`;
do echo "$BA_TOOLS/STAR/bin/Linux_x86_64/STAR --runThreadN 8 --genomeDir $BA_HOME/alignments/star/genome_index_STAR_r1 --outFileNamePrefix $I --readFilesCommand zcat --readFilesIn $BA_DATA_TRIM/$I\_1.fastq.filt.fq.gz $BA_DATA_TRIM/$I\_2.fastq.filt.fq.gz --outSJfilterReads Unique" >> $BA_ALIGN_DIR/jobs.txt; done;
sed -i -e 's/\\//g' $BA_ALIGN_DIR/jobs.txt
head -1 $BA_ALIGN_DIR/jobs.txt

mkdir -p $BA_ALIGN_DIR/1st_SJ
cd $BA_ALIGN_DIR/1st_SJ
parallel -j20 < ../jobs.txt

rm *Aligned.out.sam

##########################################
#### Alignment with STAR: Second Pass ####
##########################################

ls 1st_SJ/*SJ.out.tab > SJ_out_all
SJ_out_list=$BA_ALIGN_DIR/SJ_out_all

cd $BA_DATA_TRIM
rm $BA_ALIGN_DIR/jobs*

for I in `ls -1 *_1.fastq.filt.fq.gz | sed 's/_1.fastq.filt.fq.gz//'`;
do echo "$BA_TOOLS/STAR/bin/Linux_x86_64/STAR --runThreadN 3 --genomeDir /cerberus/projects/racste/B_anynana/alignments/star/genome_index_STAR_r1 --outFileNamePrefix $I --readFilesCommand zcat --readFilesIn /cerberus/projects/racste/B_anynana/RNA_data/trim/$I\_1.fastq.filt.fq.gz /cerberus/projects/racste/B_anynana/RNA_data/trim/$I\_2.fastq.filt.fq.gz --outSJfilterReads Unique --outSAMstrandField intronMotif --outFilterMultimapNmax 1 -sjdbFileChrStartEnd SJ_out_list" >> $BA_HOME/alignments/star/jobs_2.txt; done;
sed -i -e 's/\\//g' $BA_HOME/alignments/star/jobs_2.txt
head -1 $BA_HOME/alignments/star/jobs_2.txt

cd $BA_ALIGN_DIR
parallel -j20 < jobs_2.txt

### Rename Files
# e.g. rename s/SRR5430/TH_WET_29_/ SRR5430602*

## Parse log files
# Convert the log.txt file to a table (modified from Wheat via Dort)
echo Run,Input_reads,uniq_map_reads,uniq_map_perc,splices_annot,splices_noncanon > Log.final.out.csv

for I in `ls -1 *Log.final.out | sed 's/Log.final.out//'`;
do
grep -E 'Number of input |Uniquely mapped|sjdb|Non-canonical' $I\Log.final.out | cut -f2 -d '|' | tr -d '\r\n'| tr -s " " |\
awk '{gsub("\t",",",$0); print;}' | sed 's/^/'$I'/' > temp1.csv
cp Log.final.out.csv temp2.csv;
cat temp2.csv temp1.csv > Log.final.out.csv;
rm temp1.csv temp2.csv;
#star
#cat star_log.csv temp1.csv > star_log.csv;
done

############################
#### Sort Aligned Files ####
############################
cd $BA_ALIGN_DIR

ls *Aligned.out.sam | parallel -j20 "$BA_TOOLS/samtools-1.9/samtools view -bS {.}.sam > {.}.bam"
ls *Aligned.out.bam | parallel -j20 "$BA_TOOLS/samtools-1.9/samtools sort -o {.}_sorted.bam -T {.}_temp {.}.bam"
ls *Aligned.out_sorted.bam | parallel -j20 "$BA_TOOLS/samtools-1.9/samtools index {.}.bam"

rm *Aligned.out.bam

####################
#### Downstream ####
####################
# edgeR: continue to rStudio and r Scripts
# rMATS: continue to rMATS
