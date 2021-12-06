### Headless R script to calculate pi on a per-gene (per exon actually) basis, using per-site thetas from angsd as input
## this script will be called with one parameter which is the scaffold name
## it is called from command line, from within a loop that loops over all scaffolds one at a time
## it reads in the angsd output file of per-site thetas, which has been subsetted per scaffold so as to only read in what's needed


## set environment variables
#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)

## test if there is precisely one argument: if not, return an error
if (length(args)==0) {
  stop("One argument must be supplied (input file).n", call.=FALSE)
} else if (length(args) > 1) {
  stop("One argument must be supplied (input file).n", call.=FALSE)
}


## read in scaffold
curr.scaffold <- args

# output here
angsd_temp_folder <- "angsd_thetas/angsd_per_scaffold_temp_folder/"
thetas_per_gene <- paste("angsd_thetas_CDS_",curr.scaffold, ".txt", sep="")

# per-site thetas per scaffold also here: (subsetted from main per-site-thetas angsd output file)
thetas.dir <- angsd_temp_folder

# genes here
genes.dir <- "angsd"

# read in gene names and locations; these are only the CDSs
setwd(genes.dir)
all.genes <- read.table("CDS_bed_and_names_for_R.txt", sep="\t")

# and only use those in current scaffold
genes <- all.genes[which(all.genes$V1==curr.scaffold),]


# stop if there is less than 1 gene in this scaffold
if (length(genes$V1)==0) {
  stop("no genes on this scaffold", call.=FALSE)
}

colnames(genes) <- c("scaffold", "start", "end", "locus")

# read in thetas for current scaffold
setwd(thetas.dir)
curr.thetas <- paste(curr.scaffold, "_angsd_thetas.txt", sep="")

thetas <- read.table(curr.thetas, header=F, comment.char="")
colnames(thetas) <- c("Chromo", "Pos", "Watterson", "Pairwise", "thetaSingleton", "thetaH", "thetaL")

N <- 10 # number of haploid sequences (5 individuals = 10 chromosomes)

#loop over genes and calculate thetas
genes$thetaW <- NA
genes$pi <- NA

date()
for (i in 1:nrow(genes))
{
  sites <- (thetas$Pos > genes[i,2] & thetas$Pos <= genes[i,3])
  genes$thetaW[i]  <- sum(exp(thetas$Watterson[sites]))
  genes$pi[i]      <- sum(exp(thetas$Pairwise[sites]))
  gc()
}
date()

# get gene length to correct thetaW and theta pi for number of sites across whole exon; add 1 because first and last position count for length
genes$length <- abs(genes$end - genes$start) + 1

genes$thetaW_rel <- genes$thetaW / genes$length
genes$pi_rel <- genes$pi / genes$length


# write output file

write.table(genes, file=thetas_per_gene, sep="\t", quote=F, row.names=F)
