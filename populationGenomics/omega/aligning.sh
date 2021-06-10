# align them using MASCE
# then point to that folder
# $Bicyclus_anynana=/path/to/project
# $TOOLS=/path/to/software

$Bicyclus_anynana/SCO_Parargae_aegeria

mkdir macse_alignments
# mv *fasta macse_alignments
# cd macse_alignments
ls SCO_combined_fastas/*.fasta > files2align
# align
cat files2align | parallel -P40 'java -jar  -Xmx1000m $tools/macse_v1.01b.jar -prog alignSequences -seq {} -out_NT macse_alignments/{/.}.masce.fa -out_AA macse_alignments/{/.}.masceAA.fa'

# assess CDS
$TOOLS/AMAS/amas/AMAS.py summary -f fasta -d dna -i macse_alignments/*.masce.fa -o summary_masce_CDS.tsv
# assess AA
$TOOLS/AMAS/amas/AMAS.py summary -f fasta -d dna -i macse_alignments/*.masceAA.fa -o summary_masce_AA.tsv
