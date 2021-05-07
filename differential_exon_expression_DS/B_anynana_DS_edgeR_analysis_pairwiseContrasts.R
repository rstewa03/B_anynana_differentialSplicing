
#' ---
#' title: "B. anynana pairwise differential splicing (exon expression)"
#' output:
#'   pdf_document:
#'     keep_tex: true
#' ---
#' 
#' 
#' *Season by Family Interaction edgeR DS analysis using STAR alignment 
#'  using diffSplice to only compare two groups at a time...

#### Set WD and load packages ####
setwd("/mnt/griffin/racste/B_anynana/diff_expr/edgeR")

library(edgeR)
# library(UpSetR)
# library(ggfortify)
# library(data.table)
# library(eulerr)
# library(matrixStats)
# library(gridExtra)
library(tidyverse)
library(ggpubr)

'%ni%' <- Negate('%in%')

#### load data and metadata ####
Metadata <- read.csv("/mnt/griffin/racste/B_anynana/B_anynana_metadataForSupplement.csv")%>% 
  dplyr::select(SRR, Season, Food.treatment, Family, Body.part) %>% 
  mutate(Run = factor(SRR), 
         Season = factor(Season), 
         Family = factor(Family), 
         Body.part = factor(Body.part),
         season_stress = paste0(Season, Food.treatment))

# load dgelist object
load(file = "star_fc_B_anynana_dgelist.Rdata")
y.all <- y.all[, Metadata$Run]

y.all$samples

#### Subset data ####
DS.data <- list()
#change body part
body.part <- "thorax"

if (body.part=="abdomen") {
  cols = c("#D95F02", "#1B9E77")
} else {
  cols =c("#FC8D62","#66C2A5")
}
pchs <- c(1,2,3,4,5,6,7)

DS.data$sample.sizes <- Metadata %>% filter(Body.part == body.part) %>% group_by(Season,Family) %>% summarize(samples = n())

run<- droplevels(Metadata$Run[Metadata$Body.part == body.part] )
seasons <- factor(Metadata$Season[Metadata$Body.part== body.part])
family <- factor(Metadata$Family[Metadata$Body.part== body.part])
season_family <- factor(paste0(Metadata$Season[Metadata$Body.part== body.part],
                                  Metadata$Family[Metadata$Body.part== body.part] ))
y.exons <-y.all[,colnames(y.all$counts) %in% run]

# filter by expression
# filterbyExpr calculates theminimum effect sample size library size as lib.size*norm.factors, 

keep.exons <- filterByExpr(y.exons, group = seasons)
DS.data$keep.table <- table(keep.exons)

# alternatively, Oostra et al. 2018:
# keep.exons <- rowSums(cpm(y.exons) > 0.25) >=3
# table(keep.ab)

y.exons.filt <- y.exons[keep.exons, keep.lib.sizes=FALSE]
y.exons.filt$samples$group <- seasons
y.exons.filt$samples$family<- family
y.exons.filt$samples$season_stress <- as.factor(paste0(Metadata$Season[Metadata$Body.part== body.part],
                                                       Metadata$Food.treatment[Metadata$Body.part== body.part] ))
# TMM normalization to eliminate composition biases between libraries
y.exons.filt<- calcNormFactors(y.exons.filt)

# Design matrices (skip ahead to differential splicing)
## Testing interaction effects following the statistical design described by Gordon Smyth
### https://support.bioconductor.org/p/56568/
## and the EdgeR Users Guide
### https://www.bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf
### e.g. 3.4.2, 4.2.8, etc. 

#### Set contrasts ####
design.exons.filt.full <- model.matrix(~ 0 +season_family) # interaction effects
SeasonMainContrasts <- makeContrasts((season_familywet21 + season_familywet29 + season_familywet30 + season_familywet37 +season_familywet53 +season_familywet60 +season_familywet63)/7 -
                                       (season_familydry21 + season_familydry29 + season_familydry30 + season_familydry37 +season_familydry53 +season_familydry60 +season_familydry63)/7, 
                                     levels = design.exons.filt.full)
colnames(SeasonMainContrasts) <- "SeasonMainContrasts"
rownames(SeasonMainContrasts) <- levels(season_family)
colSums(SeasonMainContrasts)
sum(SeasonMainContrasts[which(SeasonMainContrasts > 0),])
sum(SeasonMainContrasts[which(SeasonMainContrasts < 0),])

FamilyMainContrasts <- makeContrasts(
  f21v29=(season_familydry21 + season_familywet21)/2 - (season_familydry29 + season_familywet29)/2,
  f21v30=(season_familydry21 + season_familywet21)/2 - (season_familydry30 + season_familywet30)/2,
  f21v37=(season_familydry21 + season_familywet21)/2 - (season_familydry37 + season_familywet37)/2,
  f21v53=(season_familydry21 + season_familywet21)/2 - (season_familydry53 + season_familywet53)/2,
  f21v60=(season_familydry21 + season_familywet21)/2 - (season_familydry60 + season_familywet60)/2,
  f21v63=(season_familydry21 + season_familywet21)/2 - (season_familydry63 + season_familywet63)/2,
  f29v30=(season_familydry29 + season_familywet29)/2 - (season_familydry30 + season_familywet30)/2,
  f29v37=(season_familydry29 + season_familywet29)/2 - (season_familydry37 + season_familywet37)/2,
  f29v53=(season_familydry29 + season_familywet29)/2 - (season_familydry53 + season_familywet53)/2,
  f29v60=(season_familydry29 + season_familywet29)/2 - (season_familydry60 + season_familywet60)/2,
  f29v63=(season_familydry29 + season_familywet29)/2 - (season_familydry63 + season_familywet63)/2,
  f30v37=(season_familydry30 + season_familywet30)/2 - (season_familydry37 + season_familywet37)/2,
  f30v53=(season_familydry30 + season_familywet30)/2 - (season_familydry53 + season_familywet53)/2,
  f30v60=(season_familydry30 + season_familywet30)/2 - (season_familydry60 + season_familywet60)/2,
  f30v63=(season_familydry30 + season_familywet30)/2 - (season_familydry63 + season_familywet63)/2,
  f37v53=(season_familydry37 + season_familywet37)/2 - (season_familydry53 + season_familywet53)/2,
  f37v60=(season_familydry37 + season_familywet37)/2 - (season_familydry60 + season_familywet60)/2,
  f37v63=(season_familydry37 + season_familywet37)/2 - (season_familydry63 + season_familywet63)/2,
  f53v60=(season_familydry53 + season_familywet53)/2 - (season_familydry60 + season_familywet60)/2,
  f53v63=(season_familydry53 + season_familywet53)/2 - (season_familydry63 + season_familywet63)/2,
  f60v63=(season_familydry60 + season_familywet60)/2 - (season_familydry63 + season_familywet63)/2, 
  levels = design.exons.filt.full)
rownames(FamilyMainContrasts) <- levels(season_family)
colSums(FamilyMainContrasts)
for (ct in 1:dim(FamilyMainContrasts)[2]){
  print(sum(FamilyMainContrasts[which(FamilyMainContrasts[,ct] > 0),ct]))
  print(sum(FamilyMainContrasts[which(FamilyMainContrasts[,ct] < 0),ct]))
}  

InteractionContrasts <- makeContrasts(
  sxf21v29=(season_familydry21 - season_familywet21)/2 - (season_familydry29  - season_familywet29)/2,
  sxf21v30=(season_familydry21 - season_familywet21)/2 - (season_familydry30 - season_familywet30)/2,
  sxf21v37=(season_familydry21 - season_familywet21)/2 - (season_familydry37 - season_familywet37)/2,
  sxf21v53=(season_familydry21 - season_familywet21)/2 - (season_familydry53 - season_familywet53)/2,
  sxf21v60=(season_familydry21 - season_familywet21)/2 - (season_familydry60 - season_familywet60)/2,
  sxf21v63=(season_familydry21 - season_familywet21)/2 - (season_familydry63 - season_familywet63)/2,
  sxf29v30=(season_familydry29 - season_familywet29)/2 - (season_familydry30 - season_familywet30)/2,
  sxf29v37=(season_familydry29 - season_familywet29)/2 - (season_familydry37 - season_familywet37)/2,
  sxf29v53=(season_familydry29 - season_familywet29)/2 - (season_familydry53 - season_familywet53)/2,
  sxf29v60=(season_familydry29 - season_familywet29)/2 - (season_familydry60 - season_familywet60)/2,
  sxf29v63=(season_familydry29 - season_familywet29)/2 - (season_familydry63 - season_familywet63)/2,
  sxf30v37=(season_familydry30 - season_familywet30)/2 - (season_familydry37 - season_familywet37)/2,
  sxf30v53=(season_familydry30 - season_familywet30)/2 - (season_familydry53 - season_familywet53)/2,
  sxf30v60=(season_familydry30 - season_familywet30)/2 - (season_familydry60 - season_familywet60)/2,
  sxf30v63=(season_familydry30 - season_familywet30)/2 - (season_familydry63 - season_familywet63)/2,
  sxf37v53=(season_familydry37 - season_familywet37)/2 - (season_familydry53 - season_familywet53)/2,
  sxf37v60=(season_familydry37 - season_familywet37)/2 - (season_familydry60 - season_familywet60)/2,
  sxf37v63=(season_familydry37 - season_familywet37)/2 - (season_familydry63 - season_familywet63)/2,
  sxf53v60=(season_familydry53 - season_familywet53)/2 - (season_familydry60 - season_familywet60)/2,
  sxf53v63=(season_familydry53 - season_familywet53)/2 - (season_familydry63 - season_familywet63)/2,
  sxf60v63=(season_familydry60 - season_familywet60)/2 - (season_familydry63 - season_familywet63)/2, 
  
  levels = design.exons.filt.full)
rownames(InteractionContrasts) <- levels(season_family)
colSums(InteractionContrasts)
for (ct in 1:dim(InteractionContrasts)[2]){
  print(sum(InteractionContrasts[which(InteractionContrasts[,ct] > 0),ct]))
  print(sum(InteractionContrasts[which(InteractionContrasts[,ct] < 0),ct]))
}  

# Estimate dispersion
y.exons.filt.full <- estimateDisp(y.exons.filt, design.exons.filt.full, robust=TRUE)
DS.data$full.common.dispersion <- y.exons.filt.full$common.dispersion[1]
DS.data$full.trended.dispersion <- y.exons.filt.full$trended.dispersion[1] 

# Quasilikelihood fit object
fit.exons.filt.full <- glmQLFit(y.exons.filt.full,design.exons.filt.full, robust=TRUE)

pdf(file = paste0("DS/edgeR_DS_modelFit_",body.part,".pdf"), 
    height = 8, width = 8)
par(mfrow = c(2, 1))
plotBCV(y.exons.filt.full, main = "y.exons.filt.full") #plot dispersions
plotQLDisp(fit.exons.filt.full)
dev.off()

# save(fit.exons.filt.full, file = paste0("DS/fit.exons.",body.part,".full.Rdata"))

#### differential splicing ####
fit.exons.filt.full$genes$ExonID <- rownames(fit.exons.filt.full$genes)

#### DS between seasons ####
sp.exons.filt.full_season <- diffSpliceDGE(fit.exons.filt.full,contrast = SeasonMainContrasts, geneid="GeneID", exonid="ExonID")
ngenes <- length(sp.exons.filt.full_season$gene.genes$GeneID)
nexons <- length(sp.exons.filt.full_season$genes$ExonID)
simes_season <- topSpliceDGE(sp.exons.filt.full_season, test="simes", ngenes)  %>% dplyr::mutate(effect = "season") %>% 
  mutate(adjPvalBH = p.adjust(P.Value, method = "BH", n = length(P.Value)))
write_tsv(simes_season, paste0("DS/edgeR_DS_simesDataframe_Season_", body.part,".tsv"))

exons_season <- topSpliceDGE(sp.exons.filt.full_season, test="exon", n = Inf)  %>% dplyr::mutate(effect = "season") %>% 
  mutate(adjPvalBH = p.adjust(P.Value, method = "BH", n = length(P.Value))) %>% left_join(select(simes_season, GeneID,FDR), by = c("GeneID")) %>% 
  mutate(fill_fdr = if_else(FDR.y < .05 &  FDR.x < 0.05 & logFC<0, "seasonDry",if_else(FDR.y < .05 &  FDR.x < 0.05 & logFC > 0, "seasonWet", "rest")))

season_simes_pval <- ggplot(simes_season, aes(x = P.Value)) +
  geom_histogram() +
  theme_bw() +
  ggtitle("Simes p-value distribution")
season_exons_pval <- ggplot(exons_season, aes(x = P.Value)) +
  geom_histogram() +
  theme_bw()+
  ggtitle("Exon p-value distribution")
season_exons_volc <- ggplot(exons_season, aes(x = logFC, y = -log10(FDR.x))) +
  geom_point(aes(color = fill_fdr)) +
  scale_color_manual(values = c("black", cols[1], cols[2]), guide = FALSE) +
  annotate(geom = "text", x = c(-2,2), y = c(50, 50), label = c("Dry", "Wet")) +
  theme_bw()+
  ggtitle("Exon volcano plot")

panel1 <- ggarrange(season_simes_pval, season_exons_pval, season_exons_volc, ncol = 3, labels = "AUTO")
panel1
ggsave(plot = panel1, filename =  paste0("DS/edgeR_DS_simesPvalDist_season_",body.part,".pdf"), height = 6, width = 18)
write_tsv(exons_season,  paste0("DS/edgeR_DS_season_logFC_",body.part,".tsv"))


#### DS among families ####
sp.exons.filt.full_family <- list()
simes_family <- list()
exons_family <- list()
for (i in 1:dim(FamilyMainContrasts)[2]){
  sp.exons.filt.full_family[[i]]<- diffSpliceDGE(fit.exons.filt.full, contrast = FamilyMainContrasts[,i], geneid="GeneID", exonid="ExonID")
  names(sp.exons.filt.full_family)[i] <- colnames(FamilyMainContrasts)[i]
  ngenes <- length(sp.exons.filt.full_family[[i]]$gene.genes$GeneID)
  nexons <- length(sp.exons.filt.full_family[[i]]$genes$ExonID)
  simes_family[[i]] <- topSpliceDGE(sp.exons.filt.full_family[[i]], test = "simes",ngenes )
  names(simes_family)[i] <- colnames(FamilyMainContrasts)[i]
  simes_family[[i]]$contrast <- colnames(FamilyMainContrasts)[i]
  simes_family[[i]]$adjPvalBH <- p.adjust(p = simes_family[[i]]$P.Value, n = (ngenes*dim(FamilyMainContrasts)[2]), method = "BH")
  exons_family[[i]] <- topSpliceDGE(sp.exons.filt.full_family[[i]], test = "exon",nexons )
  names(exons_family)[i] <- colnames(FamilyMainContrasts)[i]
  exons_family[[i]]$contrast <- colnames(FamilyMainContrasts)[i]
  exons_family[[i]]$adjPvalBH <- p.adjust(p = exons_family[[i]]$P.Value, n = (ngenes*dim(FamilyMainContrasts)[2]), method = "BH")
}

simes_family_all <- data.frame(simes_family[[1]])%>% 
  mutate(contrast2 = contrast) %>% separate(contrast2, into= c("group1", "group2"), sep = "v")
for (i in 1:length(simes_family)){
  temp <- data.frame(simes_family[[i]])%>% 
    mutate(contrast2 = contrast) %>% separate(contrast2, into= c("group1", "group2"), sep = "v")
  simes_family_all<- union(simes_family_all, temp)
}

simes_pval_dist_family <- ggplot(simes_family_all, aes(x = P.Value)) +
  geom_histogram() +
  theme_bw() +
  facet_grid(group1 ~group2, scales = "free") 
ggsave(simes_pval_dist_family, filename =  paste0("DS/edgeR_DS_simesPvalDist_family_",body.part,".png"), height = 20, width = 30)

simes_family_summ <- simes_family_all %>% 
  group_by(GeneID) %>% 
  mutate(minpval = min(FDR), 
         nsig = sum(FDR <0.05)) %>% 
  select(GeneID, Chr, Strand, NExons, FDR, contrast, minpval, nsig) %>% 
  spread(key = contrast, value = FDR)
write_tsv(simes_family_summ, paste0("DS/edgeR_DS_simesDataframe_family_", body.part,".tsv"))


exons_family_all <- data.frame(exons_family[[1]])%>% 
  mutate(contrast2 = contrast) %>% separate(contrast2, into= c("group1", "group2"), sep = "v")
for (i in 1:length(exons_family)){
  temp <- data.frame(exons_family[[i]])%>% 
    mutate(contrast2 = contrast) %>% separate(contrast2, into= c("group1", "group2"), sep = "v")
  exons_family_all<- union(exons_family_all, temp)
}
exons_family_all$shape <- as.factor(if_else(exons_family_all$logFC < -1 & exons_family_all$FDR < 0.05, exons_family_all$group1, 
                                  if_else(exons_family_all$logFC > 1 & exons_family_all$FDR < 0.05, exons_family_all$group2, "rest")))
levels(exons_family_all$shape)

exon_pval_dist_family <- ggplot(exons_family_all, aes(x = P.Value)) +
  geom_histogram() +
  theme_bw() +
  facet_grid(group1 ~group2, scales = "free") 
ggsave(exon_pval_dist_family, filename =  paste0("DS/edgeR_DS_exonPvalDist_family_",body.part,".png"), height = 20, width = 30)

exon_volcano_family <- ggplot(exons_family_all, aes(x = logFC, y = -log10(FDR))) +
  geom_point(aes(shape = shape, color = shape)) +
  scale_shape_manual(values = c(2,3,4,5,6,7,1,2,3,4,5,6,1), guide = FALSE)+
scale_color_manual(values = c("black","black","black","black","black","black",
                              "black","black","black","black","black","black", "gray"), guide = FALSE) +
  theme_bw() +
  facet_grid(group1 ~group2, scales = "free") 
ggsave(exon_volcano_family, filename =  paste0("DS/edgeR_DS_exonVoclanoFDR_family_",body.part,".png"), height = 20, width = 30)

exons_family_summ <- exons_family_all %>% 
  group_by(GeneID, ExonID) %>% 
  summarise(minFDR = min(FDR), 
            nsig = sum(FDR <0.05),
            rangeFC = max(logFC) - min(logFC))
write_tsv(exons_family_summ,  paste0("DS/edgeR_DS_family_logFC_",body.part,".tsv"))

#### DS between seasons among families (SxF) ####
sp.exons.filt.full_SxF <- list()
simes_SxF <- list()
exons_SxF <- list()
for (i in 1:dim(InteractionContrasts)[2]){
  sp.exons.filt.full_SxF[[i]]<- diffSpliceDGE(fit.exons.filt.full, contrast = InteractionContrasts[,i], geneid="GeneID", exonid="ExonID")
  names(sp.exons.filt.full_SxF)[i] <- colnames(InteractionContrasts)[i]
  ngenes <- length(sp.exons.filt.full_SxF[[i]]$gene.genes$GeneID)
  nexons <- length(sp.exons.filt.full_SxF[[i]]$genes$ExonID)
  simes_SxF[[i]] <- topSpliceDGE(sp.exons.filt.full_SxF[[i]], test = "simes",ngenes )
  names(simes_SxF)[i] <- colnames(InteractionContrasts)[i]
  simes_SxF[[i]]$contrast <- colnames(InteractionContrasts)[i]
  simes_SxF[[i]]$adjPvalBH <- p.adjust(p = simes_SxF[[i]]$P.Value, n = (ngenes*dim(InteractionContrasts)[2]), method = "BH")
  exons_SxF[[i]] <- topSpliceDGE(sp.exons.filt.full_SxF[[i]], test = "exon",nexons )
  names(exons_SxF)[i] <- colnames(InteractionContrasts)[i]
  exons_SxF[[i]]$contrast <- colnames(InteractionContrasts)[i]
  exons_SxF[[i]]$adjPvalBH <- p.adjust(p = exons_SxF[[i]]$P.Value, n = (ngenes*dim(InteractionContrasts)[2]), method = "BH")
}

simes_SxF_all <- data.frame(simes_SxF[[1]])%>% 
  mutate(contrast2 = contrast) %>% separate(contrast2, into= c("group1", "group2"), sep = "v")
for (i in 1:length(simes_SxF)){
  temp <- data.frame(simes_SxF[[i]])%>% 
    mutate(contrast2 = contrast) %>% separate(contrast2, into= c("group1", "group2"), sep = "v")
  simes_SxF_all<- union(simes_SxF_all, temp)
}

simes_pval_dist_SxF <- ggplot(simes_SxF_all, aes(x = P.Value)) +
  geom_histogram() +
  theme_bw() +
  facet_grid(group1 ~group2, scales = "free") 
ggsave(simes_pval_dist_SxF, filename =  paste0("DS/edgeR_DS_simesPvalDist_SxF_",body.part,".png"), height = 20, width = 30)

simes_SxF_summ <- simes_SxF_all %>% 
  group_by(GeneID) %>% 
  mutate(minpval = min(FDR), 
         nsig = sum(FDR <0.05)) %>% 
  select(GeneID, Chr, Strand, NExons, FDR, contrast, minpval, nsig) %>% 
  spread(key = contrast, value = FDR) 
write_tsv(simes_SxF_summ, paste0("DS/edgeR_DS_simesDataframe_SxF_", body.part,".tsv"))


exons_SxF_all <- data.frame(exons_SxF[[1]])%>% 
  mutate(contrast2 = contrast) %>% separate(contrast2, into= c("group1", "group2"), sep = "v")
for (i in 1:length(exons_SxF)){
  temp <- data.frame(exons_SxF[[i]])%>% 
    mutate(contrast2 = contrast) %>% separate(contrast2, into= c("group1", "group2"), sep = "v")
  exons_SxF_all<- union(exons_SxF_all, temp)
}
exons_SxF_all$shape <- as.factor(if_else(exons_SxF_all$logFC < -1 & exons_SxF_all$FDR < 0.05, exons_SxF_all$group1, 
                                            if_else(exons_SxF_all$logFC > 1 & exons_SxF_all$FDR < 0.05, exons_SxF_all$group2, "yyy")))
levels(exons_SxF_all$shape)

exon_pval_dist_SxF <- ggplot(exons_SxF_all, aes(x = P.Value)) +
  geom_histogram() +
  theme_bw() +
  facet_grid(group1 ~group2, scales = "free") 
ggsave(exon_pval_dist_SxF, filename =  paste0("DS/edgeR_DS_exonPvalDist_SxF_",body.part,".png"), height = 20, width = 30)

exons_volcano_SxF <- ggplot(exons_SxF_all, aes(x = logFC, y = -log10(FDR))) +
  geom_point(aes(shape = shape, color = shape)) +
  scale_shape_manual(values = c(2,3,4,5,6,7,1,2,3,4,5,6,1), guide = FALSE)+
  scale_color_manual(values = c("black","black","black","black","black","black",
                                "black","black","black","black","black","black", "gray"), guide = FALSE) +
  theme_bw() +
  facet_grid(group1 ~group2, scales = "free") 
ggsave(exons_volcano_SxF, filename =  paste0("DS/edgeR_DS_exonVoclanoFDR_SxF_",body.part,".png"), height = 20, width = 30)

exons_SxF_summ <- exons_SxF_all %>% 
  group_by(GeneID, ExonID) %>% 
  summarise(minFDR = min(FDR), 
         nsig = sum(FDR <0.05),
         rangeFC = max(logFC) - min(logFC))
write_tsv(exons_SxF_summ,  paste0("DS/edgeR_DS_SxF_logFC_",body.part,".tsv"))
#### Combined DS ####
simes.season.pval <- length(unique(simes_season$GeneID[simes_season$P.Value < 0.05]))
simes.season.FDR <- length(unique(simes_season$GeneID[simes_season$FDR < 0.05]))
exons.season.total <- length(unique(exons_season$ExonID[(exons_season$GeneID %in% unique(simes_season$GeneID[simes_season$FDR < 0.05]))]))
exons.season.pval <- length(unique(exons_season$ExonID[exons_season$P.Value < 0.05 & (exons_season$GeneID %in% unique(simes_season$GeneID[simes_season$FDR < 0.05]))]))
exons.season.FDR <- length(unique(exons_season$ExonID[exons_season$FDR < 0.05 & (exons_season$GeneID %in% unique(simes_season$GeneID[simes_season$FDR < 0.05]))]))

simes.family.pval <- length(unique(simes_family_all$GeneID[simes_family_all$P.Value < 0.05]))
simes.family.FDR <- length(unique(simes_family_all$GeneID[simes_family_all$FDR < 0.05]))
exons.family.total <- length(unique(exons_family_all$ExonID[(exons_family_all$GeneID %in% unique(simes_family_all$GeneID[simes_family_all$FDR < 0.05]))]))
exons.family.pval <- length(unique(exons_family_all$ExonID[exons_family_all$P.Value < 0.05 & (exons_family_all$GeneID %in% unique(simes_family_all$GeneID[simes_family_all$FDR < 0.05]))]))
exons.family.FDR <- length(unique(exons_family_all$ExonID[exons_family_all$FDR < 0.05 & (exons_family_all$GeneID %in% unique(simes_family_all$GeneID[simes_family_all$FDR < 0.05]))]))

simes.SxF.pval <- length(unique(simes_SxF_all$GeneID[simes_SxF_all$P.Value < 0.05]))
simes.SxF.FDR <- length(unique(simes_SxF_all$GeneID[simes_SxF_all$FDR < 0.05]))
exons.SxF.total <- length(unique(exons_SxF_all$ExonID[(exons_SxF_all$GeneID %in% unique(simes_SxF_all$GeneID[simes_SxF_all$FDR < 0.05]))]))
exons.SxF.pval <- length(unique(exons_SxF_all$ExonID[exons_SxF_all$P.Value < 0.05 & (exons_SxF_all$GeneID %in% unique(simes_SxF_all$GeneID[simes_SxF_all$FDR < 0.05]))]))
exons.SxF.FDR <- length(unique(exons_SxF_all$ExonID[exons_SxF_all$FDR < 0.05 & (exons_SxF_all$GeneID %in% unique(simes_SxF_all$GeneID[simes_SxF_all$FDR < 0.05]))]))

DS_summary <- as.data.frame(t(cbind(test = c("Season main", "Family main", "SxF interaction"),
      simes.pval = c(simes.season.pval, simes.family.pval, simes.SxF.pval),
      simes.FDR = c(simes.season.FDR, simes.family.FDR, simes.SxF.FDR),
      exons.in.DSgenes = c(exons.season.total, exons.family.total, exons.SxF.total),
      sig.exons.in.DSgenes = c(exons.season.pval, exons.family.pval, exons.SxF.pval))))
write_tsv(DS_summary,  paste0("DS/edgeR_DS_summary_",body.part,".tsv"))

simes.exons_season_genes <- data.frame(GeneID = unique(simes_season$GeneID[simes_season$FDR < 0.05]),
                                       sig = 1, 
                                       analysis =  paste0("adjP.value_S_",body.part))
simes.exons_SxF_genes <- data.frame(GeneID = unique(simes_SxF_all$GeneID[simes_SxF_all$FDR < 0.05]),
                                    sig = 1, 
                                    analysis =  paste0("adjP.value_I_",body.part))
simes.exons_family_genes <- data.frame(GeneID = unique(simes_family_all$GeneID[simes_family_all$FDR < 0.05]),
                                       sig = 1, 
                                       analysis =  paste0("adjP.value_F_",body.part))

DS_joined <- simes.exons_season_genes %>% 
  union(simes.exons_family_genes) %>% 
  union(simes.exons_SxF_genes) %>% 
  spread(key = analysis, value = sig) %>% 
  mutate_at(vars(-GeneID), ~ if_else(is.na(.x),0,1)) %>% 
  mutate_at(vars(-GeneID), as.logical)
for (i in 1:length(DS_joined$GeneID)){
  DS_joined$all[i] <- paste(DS_joined[i,2], DS_joined[i,3], DS_joined[i,4], sep = "_")
}
DS_joined_summ <- DS_joined %>% group_by(all) %>% 
  summarise(n=n()) 
write_tsv(DS_joined_summ,  paste0("DS/edgeR_DS_joined_summ_",body.part,".tsv"))
write_tsv(DS_joined,  paste0("DS/edgeR_DS_joined_",body.part,".tsv"))

color2 <-c( "#BC3C29", "#BC3C29", "#BC3C29")
color_transparent2 <- adjustcolor(color2, alpha.f = 0.75)
pdf(file = paste0("DS/edgeR_DS_euler_",body.part,".pdf"))
DS_euler <- plot(euler(DS_joined[,c(4,2,3)], shape = "circle"),
     quantities = list(fontsize =20),
     fill = color_transparent2,
     lty = 1,
     lty = 1, labels =list(labels= c("Season\nmain effect", "Family\nmain effect", "SxF\ninteraction"), font = 2, fontsize = 18),
     adjust_labels = T
     )
DS_euler
dev.off()

ab <- read_tsv("DS/edgeR_DS_joined_abdomen.tsv")
th <- read_tsv("DS/edgeR_DS_joined_thorax.tsv")
all <- ab[1:4] %>% left_join(th[1:4]) %>% replace(., is.na(.), FALSE)

color2 <-c( "#BC3C29", "#BC3C29", "#BC3C29")
color_transparent2 <- adjustcolor(color2, alpha.f = 0.75)
pdf(file = paste0("DS/edgeR_DS_euler_abvth.pdf"), )
plot(euler(all[,c(2,5)], shape = "circle"), #family
     quantities = list(fontsize =20),
     fill = color_transparent2,
     lty = 1,
     lty = 1, labels =list(labels= c("Abdomen", "Thorax"), font = 2, fontsize = 18)
)
plot(euler(all[,c(3,6)], shape = "circle"), #interaction
     quantities = list(fontsize =20),
     fill = color_transparent2,
     lty = 1,
     lty = 1, labels =list(labels= c("Abdomen", "Thorax"), font = 2, fontsize = 18)
)
plot(euler(all[,c(4,7)], shape = "circle"), #season
      quantities = list(fontsize =20),
      fill = color_transparent2,
      lty = 1,
      lty = 1, labels =list(labels= c("Abdomen", "Thorax"), font = 2, fontsize = 18)
)
dev.off()

#### Heatmaps #### 

####
# make sure y.exons.filt is set for abdomen
resid.exons <- calcNormFactors(DGEList(y.exons.filt$counts), method = "TMM")
pseudo_resid.exons <- log2(cpm(resid.exons) + 1) %>% 
  data.frame() %>% 
  rownames_to_column(var = "exon_id") %>% 
  mutate(gene_id = exon_id, 
         gene_id = gsub(pattern = "gene-", replacement = "", gene_id)) %>% 
  separate(gene_id, into = c("gene_id", NA), sep = "_") %>% 
  gather(-c(gene_id, exon_id), key = "SRR", value = "tmm") %>% 
  group_by(SRR, gene_id) %>% 
  mutate(mean_tmm = mean(tmm),
         resid = tmm - mean_tmm) %>% 
  ungroup %>% 
  dplyr::select(exon_id, SRR, resid) %>% 
  spread(key = SRR, value = resid)

# TMM Heatmap : between Seasons
season_list  <- as.character(Metadata$Season[Metadata$Body.part == body.part])
names(season_list) <- Metadata$SRR[Metadata$Body.part == body.part]
season_col <- c("#D95F02", "#1B9E77")
names(season_col) <- c("dry", "wet") 

season_df <-as.data.frame(season_list) 

family_list <- as.character(Metadata$Family[Metadata$Body.part == body.part])
names(family_list) <- Metadata$SRR[Metadata$Body.part == body.part]
family_col <- brewer.pal(7, "Purples")
names(family_col) <- levels(Metadata$Family)

stress_list <- as.character(Metadata$Food.treatment[Metadata$Body.part == body.part])
names(stress_list) <- Metadata$SRR[Metadata$Body.part == body.part]
stress_col <- brewer.pal(n=3, "Set1") [1:2]
names(stress_col) <- levels(as.factor(Metadata$Food.treatment))


family_pch = c(1,2,3,4,5,6,7)
names(family_pch) <- levels(Metadata$Family)

exons_season <- if_else(body.part == "abdomen", read_tsv("DS/edgeR_DS_season_logFC_abdomen.tsv"), read_tsv("DS/edgeR_DS_season_logFC_thorax.tsv"))

Season_HM_exons <- exons_season %>% 
  filter(FDR.x < 0.05 & FDR.x != 0) %>% 
  arrange(FDR.x) %>% head(5000) %>% select(ExonID) %>% as.vector()

pseudo_resid.exons_DS_season<- pseudo_resid.exons[pseudo_resid.exons$exon_id %in% Season_HM_exons$ExonID,][,2:71] %>% as.matrix()
dim(pseudo_resid.exons_DS_season)
library(ComplexHeatmap)
seasonHT <- Heatmap(pseudo_resid.exons_DS_season, 
                    top_annotation = HeatmapAnnotation(Season = season_list,
                                                       Family = family_list,
                                                       `Food Stress` = stress_list,
                                                       simple_anno_size = unit(0.125, "in"),
                                                       height = unit(0.5, "in"),
                                                       col = list(Season = season_col, Family = family_col, `Food Stress` = stress_col)
                    ),
                    name = "TMM residual\nexon expression", column_dend_height = unit(2, "cm"),
                    show_row_dend = FALSE,
                    show_row_names = T,
                    column_title = NULL,
                    row_names_gp = gpar(fontsize = 1),
                    column_names_gp = gpar(fontsize = 5),
                    # show_column_names = FALSE,
                    # row_split= 6, 
                    column_split = season_df
                    )

pdf(file = paste0("DS/BA_exonExpression_seasonHeatmap_", body.part, ".pdf"), height = 6, width = 9) 
draw(seasonHT)
