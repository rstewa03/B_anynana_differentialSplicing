# Version 2.0.2

# move files to working directory
# beware, JUM will delete extra files so do not run in same folder as log files, etc.
export BA_HOME=/cerberus/projects/racste/B_anynana
export BA_REFS_DIR=$BA_HOME/refs
export BA_TOOLS=/data/programs
export BA_ALIGN_DIR=$BA_HOME/alignments/star


# Step 1
# Run JUM_2-1.sh:

mkdir -p $BA_ALIGN_DIR/JUM_DRYvWET_multi
cd $BA_ALIGN_DIR/JUM_DRYvWET_multi
mv ../*sam ../*bam ../*out.tab .
$BA_TOOLS/JUM/JUM_2-1.sh

# Step 2
# Create subdirectories for each condition. For example, control, treatA, and treatB:

mkdir -p $BA_ALIGN_DIR/JUM_DRYvWET_multi/AB_DRY
mkdir -p $BA_ALIGN_DIR/JUM_DRYvWET_multi/AB_WET
mkdir -p $BA_ALIGN_DIR/JUM_DRYvWET_multi/TH_DRY
mkdir -p $BA_ALIGN_DIR/JUM_DRYvWET_multi/TH_WET

# Step 3
# Copy files with suffix SJ.out.tab_strand_symbol_scaled from the current directory to the corresponding condition subdirectories, respectively:

cp AB_DRY*.out.tab_strand_symbol_scaled AB_DRY/
cp AB_WET*.out.tab_strand_symbol_scaled AB_WET/
cp TH_DRY*.out.tab_strand_symbol_scaled TH_DRY/
cp TH_WET*.out.tab_strand_symbol_scaled TH_WET/

# Step 4
# Copy the file UNION_junc_coor_with_junction_ID.txt to each of the condition subdirectories:

cp UNION_junc_coor_with_junction_ID.txt AB_DRY/
cp UNION_junc_coor_with_junction_ID.txt AB_WET/
cp UNION_junc_coor_with_junction_ID.txt TH_DRY/
cp UNION_junc_coor_with_junction_ID.txt TH_WET/

# Step 5
# In each of the subdirectories, run JUM-2-2.sh as follows and then return to the current/main directory:
#--Folder: path of the downloaded JUM package
#--Threshold - for junction filtering: JUM will filter for splice junctions that have more than this # of unique reads mapped to it in at least #file_number samples out of all replicates under the condition as valid junctions for downstream analysis
#--Filenum - for junction filtering: JUM will filter for splice junctions that have more than #read_threshold of unique reads mapped to it in at least this # samples out of all replicates under the condition as valid junctions for downstream analysis
#--Condition: the name of the condition, for example, control

cd $BA_ALIGN_DIR/JUM_DRYvWET_multi/AB_DRY
$BA_TOOLS/JUM/JUM_2-2.sh --Folder $BA_TOOLS/JUM --Threshold 5 --Filenum 33 --Condition AB_DRY

cd $BA_ALIGN_DIR/JUM_DRYvWET_multi/AB_WET
$BA_TOOLS/JUM/JUM_2-2.sh --Folder $BA_TOOLS/JUM --Threshold 5 --Filenum 34 --Condition AB_WET

cd $BA_ALIGN_DIR/JUM_DRYvWET_multi/TH_DRY
$BA_TOOLS/JUM/JUM_2-2.sh --Folder $BA_TOOLS/JUM --Threshold 5 --Filenum 34 --Condition TH_DRY

cd $BA_ALIGN_DIR/JUM_DRYvWET_multi/TH_WET
$BA_TOOLS/JUM/JUM_2-2.sh --Folder $BA_TOOLS/JUM --Threshold 5 --Filenum 34 --Condition TH_WET

# Step 6
#Copy the files with suffix junction_counts.txt and formatted.txt from each of the condition subdirectories to the main/parent directory:
cd $BA_ALIGN_DIR/JUM_DRYvWET_multi/
cp AB_DRY/*junction_counts.txt .
cp AB_WET/*junction_counts.txt .
cp TH_DRY/*junction_counts.txt .
cp TH_WET/*junction_counts.txt .

cp AB_DRY/*formatted.txt .
cp AB_WET/*formatted.txt .
cp TH_DRY/*formatted.txt .
cp TH_WET/*formatted.txt .

# Step 7
# Process the bam files and the sorting first. For each *Aligned.out_sorted.bam file,
ls *Aligned.out_sorted.bam | parallel -j20 "bedtools bamtobed -i {.}.bam > {.}.bed"
ls *Aligned.out_sorted.bed | parallel -j20 "sort -k1,1 -k2,2n {.}.bed > {.}_resort.bed"
rm *_sorted.bed
rename s/_sorted_resort.bed/.sorted.bed/ *_resort.bed

# Step 8
# Run JUM_A_multi_1.sh under the main directory:
# /user/home/JUM_2.0.2/JUM_A_multi_1.sh --Folder "directory" --JuncThreshold "junction_read_threshold" --fileNum_threshold "file_number" --IRthreshold "IR_read_threshold" --Readlength "read_length" --Thread "thread_num"
  #--Folder: path of the downloaded JUM package.
  #--JuncThreshold: as in step 5
  #--fileNum_threshold: as fileNum in step 5.  In the situation when this parameter is different for different conditions, choose the biggest fileNum in step 5.
  #--IRthreshold - IR filter: JUM will filter for IR events that have more than this # of unique reads mapped to the upstream exon-intron and downstream intron-exon boundaries in as potential true IR events
  #--Readlength: the length of the RNA-seq reads
  #--Thread: number of threads for multi-threading processing of sam/bam files

cd $BA_ALIGN_DIR/JUM_DRYvWET_multi
$BA_TOOLS/JUM/JUM_A_multi_1.sh --Folder $BA_TOOLS/JUM --JuncThreshold 5 --fileNum_threshold 34 --IRthreshold 5 --Readlength 100 --Thread 20

# Step 9
# Create subdirectories for each condition for downstream intron retention analysis. For example, control_IR, treatA_IR, and treatB_IR:

mkdir -p $BA_ALIGN_DIR/JUM_DRYvWET_multi/AB_DRY_IR
mkdir -p $BA_ALIGN_DIR/JUM_DRYvWET_multi/AB_WET_IR
mkdir -p $BA_ALIGN_DIR/JUM_DRYvWET_multi/TH_DRY_IR
mkdir -p $BA_ALIGN_DIR/JUM_DRYvWET_multi/TH_WET_IR

# Step 10
# Copy files with suffix junction_counts_combined_intron_retention_event_list.txt from the current directory to the corresponding condition subdirectories, respectively:

cp AB_DRY*junction_counts_combined_intron_retention_event_list.txt AB_DRY_IR/
cp AB_WET*junction_counts_combined_intron_retention_event_list.txt AB_WET_IR/
cp TH_DRY*junction_counts_combined_intron_retention_event_list.txt TH_DRY_IR/
cp TH_WET*junction_counts_combined_intron_retention_event_list.txt TH_WET_IR/

# Step 11
# In each of the subdirectories, run the following perl script as follows and then return to the current/main directory:
# bash /user/home/JUM_2.0.2/Identify_intron_retention_event_exist_in_all_samples.pl *junction_counts_combined_intron_retention_event_list.txt "combined_file_name" "file_number"
#"combined_file_name": customized output name for that condition, see an example below
#"file_number": same as fileNum in step 5 for each condition, respectively.

cd $BA_ALIGN_DIR/JUM_DRYvWET_multi/AB_DRY_IR
perl $BA_TOOLS/JUM/Identify_intron_retention_event_exist_in_all_samples.pl *junction_counts_combined_intron_retention_event_list.txt AB_DRY_junction_counts_intron_retention_in_all_samples_list.txt 33

cd $BA_ALIGN_DIR/JUM_DRYvWET_multi/AB_WET_IR
perl $BA_TOOLS/JUM/Identify_intron_retention_event_exist_in_all_samples.pl *junction_counts_combined_intron_retention_event_list.txt AB_WET_junction_counts_intron_retention_in_all_samples_list.txt 34

cd $BA_ALIGN_DIR/JUM_DRYvWET_multi/TH_DRY_IR
perl $BA_TOOLS/JUM/Identify_intron_retention_event_exist_in_all_samples.pl *junction_counts_combined_intron_retention_event_list.txt TH_DRY_junction_counts_intron_retention_in_all_samples_list.txt 34

cd $BA_ALIGN_DIR/JUM_DRYvWET_multi/TH_WET_IR
perl $BA_TOOLS/JUM/Identify_intron_retention_event_exist_in_all_samples.pl *junction_counts_combined_intron_retention_event_list.txt TH_WET_junction_counts_intron_retention_in_all_samples_list.txt 34

# Step 12
# Combine output files from the previous step as follows:
cd $BA_ALIGN_DIR/JUM_DRYvWET_multi
cat AB_DRY_IR/AB_DRY_junction_counts_intron_retention_in_all_samples_list.txt AB_WET_IR/AB_WET_junction_counts_intron_retention_in_all_samples_list.txt TH_DRY_IR/TH_DRY_junction_counts_intron_retention_in_all_samples_list.txt TH_WET_IR/TH_WET_junction_counts_intron_retention_in_all_samples_list.txt | sort -u > All_junction_counts_intron_retention_in_all_samples_sorted_list.txt

# Step 13
# Run JUM_A_multi_2.sh under the main directory:
# $BA_TOOLS/JUM/JUM_A_multi_2.sh --Folder "directory" --JuncThreshold "junction_read_threshold" --fileNum_threshold "file_number" --IRthreshold "IR_read_threshold" --Readlength "read_length" --Thread "thread_num"
      #--Folder: path of the downloaded JUM package.
      #--JuncThreshold: as in step 5
      #--fileNum_threshold: as fileNum in step 5.  In the situation when this parameter is different for different conditions, choose the biggest fileNum in step 5.
      #--IRthreshold - IR filter: JUM will filter for IR events that have more than this # of unique reads mapped to the upstream exon-intron and downstream intron-exon boundaries in as potential true IR events
      #--Readlength: the length of the RNA-seq reads
      #--Thread: number of threads for multi-threading processing of sam/bam files

$BA_TOOLS/JUM/JUM_A_multi_2.sh --Folder $BA_TOOLS/JUM --JuncThreshold 5 --fileNum_threshold 34 --IRthreshold 5 --Readlength 100 --Thread 20 # note, this did not work all the way (11 JUne 2020). had to run the second half of the script by hand.

# after running JUM_A_multi_2.sh, make sure that the main folder contains all files:
# *Aligned.out.sam 
# *Aligned.out_sorted.bam* 
# *formatted.txt 
# *SJ.out.tab 

# and the JUM_diff directory contains:
# *combined_count.txt 
# combined_AS_JUM.gff 
# *Aligned.out_coverage.bed 
# *profiled_total_AS_event_junction_first_processing_for_JUM_reference_building.txt
# UNION_junc_coor_with_junction_ID_more_than_"$threshold"_read_in_at_least_"$file_num"_samples.txt

# Step 14
# Re-distribute the files in JUM_diff/ into subdirectories for comparing different conditions.
# Each subdirectory should contain input files for the specific condition and the control condition under comparison.
# Take the subdirectory ctrl_vs_treatA_JUM_diff/ for example, it should contain the following files for treatA and control conditions:
mkdir -p $BA_ALIGN_DIR/JUM_DRYvWET_multi/JUM_diff/AB_DRY_vs_AB_WET
mkdir -p $BA_ALIGN_DIR/JUM_DRYvWET_multi/JUM_diff/TH_DRY_vs_TH_WET
mkdir -p $BA_ALIGN_DIR/JUM_DRYvWET_multi/JUM_diff/AB_DRY_vs_TH_DRY
mkdir -p $BA_ALIGN_DIR/JUM_DRYvWET_multi/JUM_diff/AB_WET_vs_TH_WET

cd $BA_ALIGN_DIR/JUM_DRYvWET_multi/JUM_diff
cp UNION_junc_coor_with_junction_ID_more_than_5_read_in_at_least_34_samples.txt AB_DRY_vs_AB_WET/
cp UNION_junc_coor_with_junction_ID_more_than_5_read_in_at_least_34_samples.txt AB_DRY_vs_TH_DRY/
cp UNION_junc_coor_with_junction_ID_more_than_5_read_in_at_least_34_samples.txt TH_DRY_vs_TH_WET/
cp UNION_junc_coor_with_junction_ID_more_than_5_read_in_at_least_34_samples.txt AB_WET_vs_TH_WET/

cp more_than_5_profiled_total_AS_event_junction_first_processing_for_JUM_reference_building.txt AB_DRY_vs_AB_WET/
cp more_than_5_profiled_total_AS_event_junction_first_processing_for_JUM_reference_building.txt AB_DRY_vs_TH_DRY/
cp more_than_5_profiled_total_AS_event_junction_first_processing_for_JUM_reference_building.txt TH_DRY_vs_TH_WET/
cp more_than_5_profiled_total_AS_event_junction_first_processing_for_JUM_reference_building.txt AB_WET_vs_TH_WET/

cp combined_AS_JUM.gff AB_DRY_vs_AB_WET/
cp combined_AS_JUM.gff AB_DRY_vs_TH_DRY/
cp combined_AS_JUM.gff TH_DRY_vs_TH_WET/
cp combined_AS_JUM.gff AB_WET_vs_TH_WET/

cp AB_DRY*_combined_count.txt AB_DRY_vs_AB_WET/
cp AB_DRY*_combined_count.txt AB_DRY_vs_TH_DRY/

cp AB_WET*_combined_count.txt AB_DRY_vs_AB_WET/
cp AB_WET*_combined_count.txt AB_WET_vs_TH_WET/

cp TH_DRY*_combined_count.txt AB_DRY_vs_TH_DRY/
cp TH_DRY*_combined_count.txt TH_DRY_vs_TH_WET/

cp TH_WET*_combined_count.txt AB_WET_vs_TH_WET/
cp TH_WET*_combined_count.txt TH_DRY_vs_TH_WET/

cp AB_DRY*Aligned.out_coverage.bed AB_DRY_vs_AB_WET/
cp AB_DRY*ligned.out_coverage.bed AB_DRY_vs_TH_DRY/

cp AB_WET*Aligned.out_coverage.bed AB_DRY_vs_AB_WET/
cp AB_WET*Aligned.out_coverage.bed AB_WET_vs_TH_WET/

cp TH_DRY*Aligned.out_coverage.bed AB_DRY_vs_TH_DRY/
cp TH_DRY*Aligned.out_coverage.bed TH_DRY_vs_TH_WET/

cp TH_WET*Aligned.out_coverage.bed AB_WET_vs_TH_WET/
cp TH_WET*Aligned.out_coverage.bed TH_DRY_vs_TH_WET/

# Run R_script_JUM.R in each subdirectory to produce a file called AS_differential.txt
cd $BA_ALIGN_DIR/JUM_DRYvWET_multi/JUM_diff/AB_DRY_vs_AB_WET
ls *.bed | grep -E '.bed' | sed 's/Aligned.out_coverage.bed//'| sed 's/$/ /' | sed 's/.*/&&/' | cut -f1-5 -d "_" | sed '1 i \\t   condition' > experiment_design.txt
Rscript $BA_TOOLS/JUM/R_script_JUM.R experiment_design.txt > outputFile.Rout 2> errorFile.Rout

cd $BA_ALIGN_DIR/JUM_DRYvWET_multi/JUM_diff/TH_DRY_vs_TH_WET
ls *.bed | grep -E '.bed' | sed 's/Aligned.out_coverage.bed//'| sed 's/$/ /' | sed 's/.*/&&/' | cut -f1-5 -d "_" | sed '1 i \\t   condition' > experiment_design.txt
Rscript $BA_TOOLS/JUM/R_script_JUM.R experiment_design.txt > outputFile.Rout 2> errorFile.Rout

cd $BA_ALIGN_DIR/JUM_DRYvWET_multi/JUM_diff/AB_DRY_vs_TH_DRY
ls *.bed | grep -E '.bed' | sed 's/Aligned.out_coverage.bed//'| sed 's/$/ /' | sed 's/.*/&&/' | cut -f1-4 -d "_" | sed '1 i \\t   condition' > experiment_design.txt
Rscript $BA_TOOLS/JUM/R_script_JUM.R experiment_design.txt > outputFile.Rout 2> errorFile.Rout

cd $BA_ALIGN_DIR/JUM_DRYvWET_multi/JUM_diff/AB_WET_vs_TH_WET
ls *.bed | grep -E '.bed' | sed 's/Aligned.out_coverage.bed//'| sed 's/$/ /' | sed 's/.*/&&/' | cut -f1-4 -d "_" | sed '1 i \\t   condition' > experiment_design.txt
Rscript $BA_TOOLS/JUM/R_script_JUM.R experiment_design.txt > outputFile.Rout 2> errorFile.Rout

# Steps 15 & 16
# Prior to completing the analysis (JUM_C.sh), prepare a the refFlat.txt 
# cd $RNA_TOOLS
# rsync -aP rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/linux.x86_64/gff3ToGenePred ./
# cd $BA_REFS_DIR 
# $RNA_TOOLS/gff3ToGenePred  -genePredExt -geneNameAsName2 GCF_900239965.1_Bicyclus_anynana_v1.2_genomic.gff refFlat.tmp.txt
# paste <(cut -f 12 refFlat.tmp.txt) <(cut -f 1-10 refFlat.tmp.txt) > refFlat.txt
# rm refFlat.tmp.txt

# Now run JUM_B.sh in each JUM_diff subdirectory
# We use a pvalue cutoff of 1 to quantify all detected splice events
cd $BA_ALIGN_DIR/JUM_DRYvWET_multi/JUM_diff/AB_DRY_vs_AB_WET
$BA_TOOLS/JUM/JUM_B.sh --Folder $BA_TOOLS/JUM --Test pvalue --Cutoff 1 --TotalFileNum 69 --Condition1_fileNum_threshold 33 --Condition2_fileNum_threshold 34 --Condition1SampleName AB_DRY_21_691,AB_DRY_21_692,AB_DRY_21_720,AB_DRY_21_735,AB_DRY_29_702,AB_DRY_29_703,AB_DRY_29_719,AB_DRY_29_734,AB_DRY_30_688,AB_DRY_30_689,AB_DRY_30_699,AB_DRY_30_717,AB_DRY_30_730,AB_DRY_30_732,AB_DRY_37_685,AB_DRY_37_686,AB_DRY_37_697,AB_DRY_37_714,AB_DRY_37_728,AB_DRY_37_745,AB_DRY_53_683,AB_DRY_53_711,AB_DRY_53_726,AB_DRY_53_743,AB_DRY_60_681,AB_DRY_60_694,AB_DRY_60_709,AB_DRY_60_725,AB_DRY_60_740,AB_DRY_60_741,AB_DRY_63_706,AB_DRY_63_707,AB_DRY_63_722,AB_DRY_63_737 --Condition2SampleName AB_WET_21_678,AB_WET_21_679,AB_WET_21_690,AB_WET_21_704,AB_WET_29_674,AB_WET_29_675,AB_WET_29_700,AB_WET_29_701,AB_WET_29_718,AB_WET_29_733,AB_WET_30_687,AB_WET_30_698,AB_WET_30_715,AB_WET_30_716,AB_WET_30_729,AB_WET_30_731,AB_WET_37_684,AB_WET_37_696,AB_WET_37_712,AB_WET_37_713,AB_WET_37_727,AB_WET_37_744,AB_WET_53_695,AB_WET_53_710,AB_WET_53_742,AB_WET_60_680,AB_WET_60_708,AB_WET_60_723,AB_WET_60_724,AB_WET_60_738,AB_WET_60_739,AB_WET_63_693,AB_WET_63_705,AB_WET_63_721,AB_WET_63_736
cd FINAL_JUM_OUTPUT_pvalue_1
$BA_TOOLS/JUM/JUM_C.sh --Folder $BA_TOOLS/JUM --Test pvalue --Cutoff 1 --TotalCondition1FileNum 34 --TotalCondition2FileNum 35 --REF $BA_REFS_DIR/refFlat.txt
 
cd $BA_ALIGN_DIR/JUM_DRYvWET_multi/JUM_diff/TH_DRY_vs_TH_WET
$BA_TOOLS/JUM/JUM_B.sh --Folder $BA_TOOLS/JUM --Test pvalue --Cutoff 1 --TotalFileNum 70 --Condition1_fileNum_threshold 34 --Condition2_fileNum_threshold 34 --Condition1SampleName TH_DRY_21_619,TH_DRY_21_620,TH_DRY_21_648,TH_DRY_21_663,TH_DRY_29_604,TH_DRY_29_605,TH_DRY_29_630,TH_DRY_29_631,TH_DRY_29_647,TH_DRY_29_662,TH_DRY_30_616,TH_DRY_30_617,TH_DRY_30_627,TH_DRY_30_645,TH_DRY_30_658,TH_DRY_30_660,TH_DRY_37_613,TH_DRY_37_614,TH_DRY_37_625,TH_DRY_37_642,TH_DRY_37_656,TH_DRY_37_673,TH_DRY_53_611,TH_DRY_53_639,TH_DRY_53_654,TH_DRY_53_671,TH_DRY_60_609,TH_DRY_60_622,TH_DRY_60_653,TH_DRY_60_668,TH_DRY_60_669,TH_DRY_63_634,TH_DRY_63_635,TH_DRY_63_650,TH_DRY_63_665 --Condition2SampleName TH_WET_21_606,TH_WET_21_607,TH_WET_21_618,TH_WET_21_632,TH_WET_29_602,TH_WET_29_603,TH_WET_29_628,TH_WET_29_629,TH_WET_29_646,TH_WET_29_661,TH_WET_30_615,TH_WET_30_626,TH_WET_30_643,TH_WET_30_644,TH_WET_30_657,TH_WET_30_659,TH_WET_37_612,TH_WET_37_624,TH_WET_37_640,TH_WET_37_641,TH_WET_37_655,TH_WET_37_672,TH_WET_53_610,TH_WET_53_623,TH_WET_53_638,TH_WET_53_670,TH_WET_60_608,TH_WET_60_651,TH_WET_60_652,TH_WET_60_666,TH_WET_60_667,TH_WET_63_621,TH_WET_63_633,TH_WET_63_649,TH_WET_63_664
cd FINAL_JUM_OUTPUT_pvalue_1
$BA_TOOLS/JUM/JUM_C.sh --Folder $BA_TOOLS/JUM --Test pvalue --Cutoff 1 --TotalCondition1FileNum 35 --TotalCondition2FileNum 35 --REF $BA_REFS_DIR/refFlat.txt
 
cd $BA_ALIGN_DIR/JUM_DRYvWET_multi/JUM_diff/AB_DRY_vs_TH_DRY
$BA_TOOLS/JUM/JUM_B.sh --Folder $BA_TOOLS/JUM/ --Test pvalue --Cutoff 1 --TotalFileNum 69 --Condition1_fileNum_threshold 33 --Condition2_fileNum_threshold 34 --Condition1SampleName AB_DRY_21_691,AB_DRY_21_692,AB_DRY_21_720,AB_DRY_21_735,AB_DRY_29_702,AB_DRY_29_703,AB_DRY_29_719,AB_DRY_29_734,AB_DRY_30_688,AB_DRY_30_689,AB_DRY_30_699,AB_DRY_30_717,AB_DRY_30_730,AB_DRY_30_732,AB_DRY_37_685,AB_DRY_37_686,AB_DRY_37_697,AB_DRY_37_714,AB_DRY_37_728,AB_DRY_37_745,AB_DRY_53_683,AB_DRY_53_711,AB_DRY_53_726,AB_DRY_53_743,AB_DRY_60_681,AB_DRY_60_694,AB_DRY_60_709,AB_DRY_60_725,AB_DRY_60_740,AB_DRY_60_741,AB_DRY_63_706,AB_DRY_63_707,AB_DRY_63_722,AB_DRY_63_737 --Condition2SampleName TH_DRY_21_619,TH_DRY_21_620,TH_DRY_21_648,TH_DRY_21_663,TH_DRY_29_604,TH_DRY_29_605,TH_DRY_29_630,TH_DRY_29_631,TH_DRY_29_647,TH_DRY_29_662,TH_DRY_30_616,TH_DRY_30_617,TH_DRY_30_627,TH_DRY_30_645,TH_DRY_30_658,TH_DRY_30_660,TH_DRY_37_613,TH_DRY_37_614,TH_DRY_37_625,TH_DRY_37_642,TH_DRY_37_656,TH_DRY_37_673,TH_DRY_53_611,TH_DRY_53_639,TH_DRY_53_654,TH_DRY_53_671,TH_DRY_60_609,TH_DRY_60_622,TH_DRY_60_653,TH_DRY_60_668,TH_DRY_60_669,TH_DRY_63_634,TH_DRY_63_635,TH_DRY_63_650,TH_DRY_63_665
cd FINAL_JUM_OUTPUT_pvalue_1
$BA_TOOLS/JUM/JUM_C.sh --Folder $BA_TOOLS/JUM/ --Test pvalue --Cutoff 1 --TotalCondition1FileNum 34 --TotalCondition2FileNum 35 --REF $BA_REFS_DIR/refFlat.txt
/cerberus/programs/JUM/JUM_C.sh /cerberus/programs/JUM --Test pvalue --Cutoff 1 --TotalCondition1FileNum 34 --TotalCondition2FileNum 35 --REF $BA_REFS_DIR/refFlat.txt

cd $BA_ALIGN_DIR/JUM_DRYvWET_multi/JUM_diff/AB_WET_vs_TH_WET
$BA_TOOLS/JUM/JUM_B.sh --Folder $BA_TOOLS/JUM --Test pvalue --Cutoff 1 --TotalFileNum 70 --Condition1_fileNum_threshold 34 --Condition2_fileNum_threshold 34 --Condition1SampleName AB_WET_21_678,AB_WET_21_679,AB_WET_21_690,AB_WET_21_704,AB_WET_29_674,AB_WET_29_675,AB_WET_29_700,AB_WET_29_701,AB_WET_29_718,AB_WET_29_733,AB_WET_30_687,AB_WET_30_698,AB_WET_30_715,AB_WET_30_716,AB_WET_30_729,AB_WET_30_731,AB_WET_37_684,AB_WET_37_696,AB_WET_37_712,AB_WET_37_713,AB_WET_37_727,AB_WET_37_744,AB_WET_53_695,AB_WET_53_710,AB_WET_53_742,AB_WET_60_680,AB_WET_60_708,AB_WET_60_723,AB_WET_60_724,AB_WET_60_738,AB_WET_60_739,AB_WET_63_693,AB_WET_63_705,AB_WET_63_721,AB_WET_63_736 --Condition2SampleName TH_WET_21_606,TH_WET_21_607,TH_WET_21_618,TH_WET_21_632,TH_WET_29_602,TH_WET_29_603,TH_WET_29_628,TH_WET_29_629,TH_WET_29_646,TH_WET_29_661,TH_WET_30_615,TH_WET_30_626,TH_WET_30_643,TH_WET_30_644,TH_WET_30_657,TH_WET_30_659,TH_WET_37_612,TH_WET_37_624,TH_WET_37_640,TH_WET_37_641,TH_WET_37_655,TH_WET_37_672,TH_WET_53_610,TH_WET_53_623,TH_WET_53_638,TH_WET_53_670,TH_WET_60_608,TH_WET_60_651,TH_WET_60_652,TH_WET_60_666,TH_WET_60_667,TH_WET_63_621,TH_WET_63_633,TH_WET_63_649,TH_WET_63_664
cd FINAL_JUM_OUTPUT_pvalue_1
$BA_TOOLS/JUM/JUM_C.sh --Folder $BA_TOOLS/JUM --Test pvalue --Cutoff 1 --TotalCondition1FileNum 35 --TotalCondition2FileNum 35 --REF $BA_REFS_DIR/refFlat.txt
 
