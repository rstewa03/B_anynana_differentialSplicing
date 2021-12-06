### R script to average per-exon pi (as calculated by angsd)
### across all exons per gene, taking into account exon length

# read data
thetas.dir <- "angsd_thetas"
setwd(thetas.dir)
per_exon_thetas <- read.table("angsd_thetas_genes_CDS_exons.txt", header=T, comment.char="", sep="\t")

# how many exons per gene
dim(per_exon_thetas) # 103293 exons
levels(per_exon_thetas$locus) # 14402 loci 

## take length-weighted average for thetaW_rel and pi_rel across exons of same gene
# make new table per locus
gene.table <- data.frame(scaffold = NA, start = NA, end = NA, 
	locus = aggregate(per_exon_thetas$length~per_exon_thetas$locus, FUN=sum)[,1],
	length = aggregate(per_exon_thetas$length~per_exon_thetas$locus, FUN=sum)[,2],
	thetaW_rel = NA, pi_rel = NA,
	thetaW_rel.uncorr_ave = NA, pi_rel.uncorr_ave = NA,
	thetaW_rel.var = NA, pi_rel.var = NA)
head(gene.table)

# use length-weighted averages for thetas; also get variances
# and loop through table to fill with relevant field
for(ll in gene.table$locus){
	curr.loc <- ll	
	curr.exons <- per_exon_thetas[which(per_exon_thetas$locus == curr.loc),]
	nr.exons <- dim(curr.exons)[1]
	curr.length <- gene.table[which(gene.table$locus == curr.loc),]$length
	# print(ll)
	# print(curr.loc)
	# dim(curr.exons)
	# print(curr.length)
	gene.table[which(gene.table$locus == curr.loc),]$thetaW_rel <- sum(curr.exons$length / curr.length * curr.exons$thetaW_rel)
	gene.table[which(gene.table$locus == curr.loc),]$pi_rel <- sum(curr.exons$length / curr.length * curr.exons$pi_rel)
	gene.table[which(gene.table$locus == curr.loc),]$thetaW_rel.uncorr_ave <- mean(curr.exons$thetaW_rel, na.rm=T)
	gene.table[which(gene.table$locus == curr.loc),]$pi_rel.uncorr_ave <- mean(curr.exons$pi_rel, na.rm=T)
	gene.table[which(gene.table$locus == curr.loc),]$thetaW_rel.var <- var(curr.exons$thetaW_rel, na.rm=T)
	gene.table[which(gene.table$locus == curr.loc),]$pi_rel.var <- var(curr.exons$pi_rel, na.rm=T)
	gene.table[which(gene.table$locus == curr.loc),]$scaffold <- as.character(curr.exons[1,]$scaffold)
	gene.table[which(gene.table$locus == curr.loc),]$start <- curr.exons[1,]$start
	gene.table[which(gene.table$locus == curr.loc),]$end <- curr.exons[nr.exons,]$end
}

write.table(gene.table, "angsd_thetas_genes_CDS.txt", sep="\t", quote=F, row.names=F)

