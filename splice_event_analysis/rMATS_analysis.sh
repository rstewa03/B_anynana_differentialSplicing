#Step 1: make rMATS comparison text files with paths to sample bam Files
cd $BA_ALIGN_DIR/JUM_DRYvWET_multi/
ls AB_DRY* | grep -E '.bam' | sed -e 's?^?'$path'/''?'| tr '\r\n' ',' | sed 's/.$//' > $BA_HOME/diff_expr/08_rMATS_S/rMATS_AB_DRY.txt
ls AB_WET* | grep -E '.bam' | sed -e 's?^?'$path'/''?'| tr '\r\n' ',' | sed 's/.$//' > $BA_HOME/diff_expr/08_rMATS_S/rMATS_AB_WET.txt
ls TH_DRY* | grep -E '.bam' | sed -e 's?^?'$path'/''?'| tr '\r\n' ',' | sed 's/.$//' > $BA_HOME/diff_expr/08_rMATS_S/rMATS_TH_DRY.txt
ls TH_WET* | grep -E '.bam' | sed -e 's?^?'$path'/''?'| tr '\r\n' ',' | sed 's/.$//' > $BA_HOME/diff_expr/08_rMATS_S/rMATS_TH_WET.txt

ls AB_DRY_29* | grep -E '.bam' | sed -e 's?^?'$path'/''?'| tr '\r\n' ',' | sed 's/.$//'  > $BA_HOME/diff_expr/08_rMATS_S/rMATS_AB_DRY_29.txt
ls AB_WET_29* | grep -E '.bam' | sed -e 's?^?'$path'/''?'| tr '\r\n' ',' | sed 's/.$//'  > $BA_HOME/diff_expr/08_rMATS_S/rMATS_AB_WET_29.txt


#Step 2: Run rMATS
# test ith family 29
cd $BA_HOME/diff_expr/08_rMATS_S/

# test with family 29
python $BA_TOOLS/rmats-turbo/rmats.py --b1 rMATS_AB_DRY_29.txt --b2 rMATS_AB_WET_29.txt --gtf /cerberus/projects/racste/B_anynana/refs/GCF_900239965.1_Bicyclus_anynana_v1.2_genomic.2.gtf --od rMATS_AB_DRY_vs_AB_WET --tmp rMATS_AB_DRY_vs_AB_WET_temp -t paired --readLength 100 --cstat 0.0001 --libType fr-unstranded --nthread 8

# gtf details:
# There are 15845 distinct gene ID in the gtf file
# There are 23362 distinct transcript ID in the gtf file
# There are 12345 one-transcript genes in the gtf file
# There are 196681 exons in the gtf file
# There are 1916 one-exon transcripts in the gtf file
# There are 1908 one-transcript genes with only one exon in the transcript
# Average number of transcripts per gene is 1.474408
# Average number of exons per transcript is 8.418843
# Average number of exons per transcript excluding one-exon tx is 9.081647
# Average number of gene per geneGroup is 60.550649
# statistic: 0.00302886962891
# novel: 3496.40151715
# The splicing graph and candidate read have been saved into rMATS_AB_DRY_vs_AB_WET_temp/2020-06-23-15:53:45_449338.rmats
# save: 11.6561050415
# loadsg: 0.139969110489
#
# ==========
# Done processing each gene from dictionary to compile AS events
# Found 7019 exon skipping events
# Found 1014 exon MX events
# Found 2053 alt SS events
# There are 1152 alt 3 SS events and 901 alt 5 SS events.
# Found 191 RI events
# ==========
#
# ase: 0.575726032257
# count: 17.0569419861
# Processing count files.
# Done processing count files.


#output contains:
#     [AS_Event].MATS.JC.txt: Final output including only reads that span junctions defined by rmats (Junction Counts)
#     [AS_Event].MATS.JCEC.txt: Final output including both reads that span junctions defined by rmats (Junction Counts) and reads that do not cross an exon boundary (Exon Counts)
#     fromGTF.[AS_Event].txt: All identified alternative splicing (AS) events derived from GTF and RNA
#     fromGTF.novelJunction.[AS_Event].txt: Alternative splicing (AS) events which were identified only after considering the RNA (as opposed to analyzing the GTF in isolation). This does not include events with an unannotated splice site.
#     fromGTF.novelSpliceSite.[AS_Event].txt: This file contains only those events which include an unannotated splice site. Only relevant if --novelSS is enabled.
#     JC.raw.input.[AS_Event].txt: Event counts including only reads that span junctions defined by rmats (Junction Counts)
#     JCEC.raw.input.[AS_Event].txt: Event counts including both reads that span junctions defined by rmats (Junction Counts) and reads that do not cross an exon boundary (Exon Counts)
#     Shared columns:
#         ID: rMATS event id
#         GeneID: Gene id
#         geneSymbol: Gene name
#         chr: Chromosome
#         strand: Strand of the gene
#         IJC_SAMPLE_1: Inclusion counts for sample 1. Replicates are comma separated
#         SJC_SAMPLE_1: Skipping counts for sample 1. Replicates are comma separated
#         IJC_SAMPLE_2: Inclusion counts for sample 2. Replicates are comma separated
#         SJC_SAMPLE_2: Skipping counts for sample 2. Replicates are comma separated
#         IncFormLen: Length of inclusion form, used for normalization
#         SkipFormLen: Length of skipping form, used for normalization
#         PValue: Significance of splicing difference between the two sample groups. (Only available if the statistical model is on)
#         FDR: False Discovery Rate calculated from p-value. (Only available if statistical model is on)
#         IncLevel1: Inclusion level for sample 1. Replicates are comma separated. Calculated from normalized counts
#         IncLevel2: Inclusion level for sample 2. Replicates are comma separated. Calculated from normalized counts
#         IncLevelDifference: average(IncLevel1) - average(IncLevel2)
#     Event specific columns (event coordinates):
#         SE: exonStart_0base exonEnd upstreamES upstreamEE downstreamES downstreamEE
#             The inclusion form includes the target exon (exonStart_0base, exonEnd)
#         MXE: 1stExonStart_0base 1stExonEnd 2ndExonStart_0base 2ndExonEnd upstreamES upstreamEE downstreamES downstreamEE
#             If the strand is + then the inclusion form includes the 1st exon (1stExonStart_0base, 1stExonEnd) and skips the 2nd exon
#             If the strand is - then the inclusion form includes the 2nd exon (2ndExonStart_0base, 2ndExonEnd) and skips the 1st exon
#         A3SS, A5SS: longExonStart_0base longExonEnd shortES shortEE flankingES flankingEE
#             The inclusion form includes the long exon (longExonStart_0base, longExonEnd) instead of the short exon (shortES shortEE)
#         RI: riExonStart_0base riExonEnd upstreamES upstreamEE downstreamES downstreamEE
#             The inclusion form includes (retains) the intron (upstreamEE, downstreamES)
#
# --tmp contains the intermediate files generated by the prep step:
#
#     {datetime}.rmats: Summary generated from processing the BAM(s)

python $BA_TOOLS/rmats-turbo/rmats.py --b1 rMATS_AB_DRY.txt --b2 rMATS_AB_WET.txt --gtf /cerberus/projects/racste/B_anynana/refs/GCF_900239965.1_Bicyclus_anynana_v1.2_genomic.2.gtf --od rMATS_AB_DRY_vs_AB_WET --tmp rMATS_AB_DRY_vs_AB_WET_temp -t paired --readLength 100 --cstat 0.0001 --libType fr-unstranded --novelSS --variable-read-length --nthread 20


python $BA_TOOLS/rmats-turbo/rmats.py --b1 rMATS_TH_DRY.txt --b2 rMATS_TH_WET.txt --gtf /cerberus/projects/racste/B_anynana/refs/GCF_900239965.1_Bicyclus_anynana_v1.2_genomic.2.gtf --od rMATS_TH_DRY_vs_TH_WET --tmp rMATS_TH_DRY_vs_TH_WET_temp -t paired --readLength 100 --cstat 0.0001 --libType fr-unstranded --novelSS --variable-read-length --nthread 20
# statistic: 0.0028440952301
# novel: 7916.02216911
# The splicing graph and candidate read have been saved into rMATS_TH_DRY_vs_TH_WET_temp/2020-06-24-11:31:36_351954.rmats
# save: 1885.61783385
# loadsg: 1.83778619766
#
# ==========
# Done processing each gene from dictionary to compile AS events
# Found 38095 exon skipping events
# Found 45609 exon MX events
# Found 56494 alt SS events
# There are 28551 alt 3 SS events and 27943 alt 5 SS events.
# Found 3360 RI events
# ==========
#
# ase: 14.074488163
# count: 3984.72542

python $BA_TOOLS/rmats-turbo/rmats.py --b1 rMATS_AB_DRY.txt --b2 rMATS_TH_DRY.txt --gtf /cerberus/projects/racste/B_anynana/refs/GCF_900239965.1_Bicyclus_anynana_v1.2_genomic.2.gtf --od rMATS_AB_DRY_vs_TH_DRY --tmp rMATS_AB_DRY_vs_TH_DRY_temp -t paired --readLength 100 --cstat 0.0001 --libType fr-unstranded --novelSS --variable-read-length --nthread 20
# statistic: 0.00305891036987
# novel: 9914.7808919
# The splicing graph and candidate read have been saved into rMATS_AB_DRY_vs_TH_DRY_temp/2020-06-23-23:26:17_681523.rmats
# save: 2102.81114292
# loadsg: 3.43366503716
#
# ==========
# Done processing each gene from dictionary to compile AS events
# Found 49488 exon skipping events
# Found 30883 exon MX events
# Found 74783 alt SS events
# There are 37572 alt 3 SS events and 37211 alt 5 SS events.
# Found 4638 RI events
# ==========
#
# ase: 15.8693890572
# count: 2316.30162716
# Processing count files.
# Done processing count files.

python $BA_TOOLS/rmats-turbo/rmats.py --b1 rMATS_AB_WET.txt --b2 rMATS_TH_WET.txt --gtf /cerberus/projects/racste/B_anynana/refs/GCF_900239965.1_Bicyclus_anynana_v1.2_genomic.2.gtf --od rMATS_AB_WET_vs_TH_WET --tmp rMATS_AB_WET_vs_TH_WET_temp -t paired --readLength 100 --cstat 0.0001 --libType fr-unstranded --novelSS --variable-read-length --nthread 20
# statistic: 0.00272607803345
# novel: 8740.44043183
# The splicing graph and candidate read have been saved into rMATS_AB_WET_vs_TH_WET_temp/2020-06-24-03:38:40_768804.rmats
# save: 1621.72239399
# loadsg: 2.52596712112
#
# ==========
# Done processing each gene from dictionary to compile AS events
# Found 47358 exon skipping events
# Found 37703 exon MX events
# Found 76655 alt SS events
# There are 36632 alt 3 SS events and 40023 alt 5 SS events.
# Found 4639 RI events
# ==========
#
# ase: 12.3522670269
# count: 2457.12784982
# Processing count files.
# Done processing count files.
