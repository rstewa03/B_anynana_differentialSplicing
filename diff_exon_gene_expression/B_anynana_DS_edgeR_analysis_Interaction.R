setwd("/mnt/griffin/racste/B_anynana/diff_expr/edgeR")

library(edgeR)
library(UpSetR)
library(ggfortify)
library(data.table)
library(eulerr)
library(matrixStats)
library(gridExtra)
library(tidyverse)

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
load(file = "star_fc_genes_B_anynana_dgelist.Rdata")
y.genes.all <- y.genes.all[, Metadata$Run]



#### Subset data ####
DE.data <- list()
body.part <- "abdomen"

DE.data$sample.sizes <- Metadata %>% filter(Body.part == body.part) %>% group_by(Season,Family) %>% summarize(samples = n())

run<- droplevels(Metadata$Run[Metadata$Body.part == body.part] )
seasons <- factor(Metadata$Season[Metadata$Body.part== body.part])
family <- factor(Metadata$Family[Metadata$Body.part== body.part])

y.genes <-y.genes.all[,colnames(y.genes.all$counts) %in% run]

# filter by expression
# filterbyExpr calculates theminimum effect sample size library size as lib.size*norm.factors, 

keep.genes <- filterByExpr(y.genes, group = seasons)
DE.data$keep.table <- table(keep.genes)

# alternatively, Oostra et al. 2018:
# keep.genes <- rowSums(cpm(y.genes) > 0.25) >=3
# table(keep.genes)

y.genes.filt <- y.genes[keep.genes, keep.lib.sizes=FALSE]
y.genes.filt$samples$group <- seasons
y.genes.filt$samples$family<- family
y.genes.filt$samples$season_stress <- as.factor(paste0(Metadata$Season[Metadata$Body.part== body.part],
                                                       Metadata$Food.treatment[Metadata$Body.part== body.part] ))
# TMM normalization to eliminate composition biases between libraries
y.genes.filt<- calcNormFactors(y.genes.filt)

# Design matrices (skip ahead to differential splicing)
## Testing interaction effects following the statistical design described by Gordon Smyth
### https://support.bioconductor.org/p/56568/
## and the EdgeR Users Guide
### https://www.bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf
### e.g. 3.4.2, 4.2.8, etc. 

design.genes.filt.SxF<- model.matrix(~ 1 + family * seasons) # interaction effects
design.genes.filt.SandF<- model.matrix(~ 1 + family + seasons) # main effects

# Estimate dispersion
y.genes.filt.SxF <- estimateDisp(y.genes.filt, design.genes.filt.SxF, robust=TRUE)
DE.data$SxF.common.dispersion <- y.genes.filt.SxF$common.dispersion[1]
DE.data$SxF.trended.dispersion <- y.genes.filt.SxF$trended.dispersion[1] 

y.genes.filt.SandF <- estimateDisp(y.genes.filt, design.genes.filt.SandF, robust=TRUE)
DE.data$SandF.common.dispersion <- y.genes.filt.SandF$common.dispersion[1]
DE.data$SandF.trended.dispersion <- y.genes.filt.SandF$trended.dispersion[1] 

# Quasilikelihood fit object
fit.genes.filt.SxF <- glmQLFit(y.genes.filt.SxF,design.genes.filt.SxF, robust=TRUE)
fit.genes.filt.SandF <- glmQLFit(y.genes.filt.SandF,design.genes.filt.SandF, robust=TRUE)

pdf(file = paste0("DE/edgeR_DE_modelFit_",body.part,".pdf"), 
    height = 8, width = 8)
par(mfrow = c(2, 2))
plotBCV(y.genes.filt.SxF, main = "y.genes.filt.SxF") #plot dispersions
plotBCV(y.genes.filt.SandF, main = "y.genes.filt.SandF")
plotQLDisp(fit.genes.filt.SxF)
plotQLDisp(fit.genes.filt.SandF)
dev.off()

save(fit.genes.filt.SxF, file = paste0("DE/fit.genes.",body.part,".SxF.Rdata"))
save(fit.genes.filt.SandF, file = paste0("DE/fit.genes.",body.part,".SandF.Rdata"))

#### differential expression: filterByExpr ####
# load(file = paste0("DE/fit.genes.",body.part,".SxF.Rdata"))
# load(file = paste0("DE/fit.genes.",body.part,".SandF.Rdata"))

# ANOVA-like analysis of differential gene expression 

qlf.genes.ab.filt.SxF_interaction <- glmQLFTest(fit.genes.filt.SxF, coef = 9:14) # detects genes that respond differently in the wet season, relative to the dry season, in any of the families
qlf.genes.ab.filt.SandF_season <- glmQLFTest(fit.genes.filt.SandF, coef = 8) # Season main effect
qlf.genes.ab.filt.SandF_family <- glmQLFTest(fit.genes.filt.SandF, coef = 2:7)

DE.data$genes <- length(qlf.genes.ab.filt.SandF_season$genes$GeneID)

qlf.genes <- list()
qlf.genes$interaction <- (topTags(qlf.genes.ab.filt.SxF_interaction, DE.data$genes))$table %>% mutate(effect = "interaction")
qlf.genes$season <- (topTags(qlf.genes.ab.filt.SandF_season, DE.data$genes))$table %>% mutate(effect = "season")
qlf.genes$family <- (topTags(qlf.genes.ab.filt.SandF_family, DE.data$genes))$table %>% mutate(effect = "family")

# adjusted p-value for differentially spliced genes
all.genes.pvals <- rbind(dplyr::select(qlf.genes$interaction, GeneID, PValue, FDR, effect), 
                   dplyr::select(qlf.genes$season, GeneID, PValue, FDR, effect), 
                   dplyr::select(qlf.genes$family, GeneID, PValue, FDR, effect)) %>% 
  mutate(all.adjp = p.adjust(PValue, method="BH"))

qlf.genes$interaction$adjP.Value <- all.genes.pvals$all.adjp[all.genes.pvals$effect == "interaction"]
qlf.genes$season$adjP.Value <- all.genes.pvals$all.adjp[all.genes.pvals$effect == "season"]
qlf.genes$family$adjP.Value <- all.genes.pvals$all.adjp[all.genes.pvals$effect == "family"]

for (i in 1:3){
  write_tsv(qlf.genes[[i]], paste0("DE/edgeR_DE_qlfGenes_",body.part,"_", names(qlf.genes)[i],".tsv")) 
}

DE.data$top.genes$interaction <-qlf.genes$interaction$GeneID[qlf.genes$interaction$adjP.Value <= 0.05]
DE.data$top.genes$season <-qlf.genes$season$GeneID[qlf.genes$season$adjP.Value <= 0.05]
DE.data$top.genes$family <-qlf.genes$family$GeneID[qlf.genes$family$adjP.Value <= 0.05]

# only multiexon genes
qlf.genes.multiexon <- list()
qlf.genes.multiexon$interaction <- cbind(qlf.genes$interaction, multiexon = str_extract(qlf.genes$interaction$Start, ";+")) %>% 
  filter(multiexon == ";") %>% 
  dplyr::select(-multiexon)
qlf.genes.multiexon$season <- cbind(qlf.genes$season, multiexon = str_extract(qlf.genes$season$Start, ";+")) %>% 
  filter(multiexon == ";") %>% 
  dplyr::select(-multiexon)
qlf.genes.multiexon$family <- cbind(qlf.genes$family, multiexon = str_extract(qlf.genes$family$Start, ";+")) %>% 
  filter(multiexon == ";") %>% 
  dplyr::select(-multiexon)

DE.data$genes_multiexon <- length(qlf.genes.multiexon$family$GeneID)

# adjusted p-value for differentially spliced genes
all.genes.pvals.multiexon <- rbind(dplyr::select(qlf.genes.multiexon$interaction, GeneID, PValue, FDR, effect), 
                             dplyr::select(qlf.genes.multiexon$season, GeneID, PValue, FDR, effect), 
                             dplyr::select(qlf.genes.multiexon$family, GeneID, PValue, FDR, effect)) %>% 
  mutate(all.adjp = p.adjust(PValue, method="BH"))

qlf.genes.multiexon$interaction$adjP.Value <- all.genes.pvals.multiexon$all.adjp[all.genes.pvals.multiexon$effect == "interaction"]
qlf.genes.multiexon$season$adjP.Value <- all.genes.pvals.multiexon$all.adjp[all.genes.pvals.multiexon$effect == "season"]
qlf.genes.multiexon$family$adjP.Value <- all.genes.pvals.multiexon$all.adjp[all.genes.pvals.multiexon$effect == "family"]

family_logfc_genes <- cbind(qlf.genes.multiexon$family[8:13], rep(0,DE.data$genes_multiexon))
qlf.genes.multiexon$family <- left_join(qlf.genes.multiexon$family,
                                        cbind(GeneID = qlf.genes.multiexon$family$GeneID, 
                                              data.frame(Means = rowMeans(as.matrix(family_logfc_genes[,1:6]))),
                                              data.frame(Ranges = rowMaxs(as.matrix(family_logfc_genes))-rowMins(as.matrix(family_logfc_genes)))))

for (i in 1:3){
  write_tsv(qlf.genes.multiexon[[i]], paste0("DE/edgeR_DE_qlfMultiexonGenes_",body.part,"_", names(qlf.genes.multiexon)[i],".tsv")) 
}

interaction_logfc_genes <- cbind(qlf.genes.multiexon$interaction[8:13], rep(0,DE.data$genes_multiexon))
qlf.genes.multiexon$interaction <- left_join(qlf.genes.multiexon$interaction,
                                        cbind(GeneID = qlf.genes.multiexon$interaction$GeneID, 
                                              data.frame(Means = rowMeans(as.matrix(interaction_logfc_genes[,1:6]))),
                                              data.frame(Ranges = rowMaxs(as.matrix(interaction_logfc_genes))-rowMins(as.matrix(interaction_logfc_genes)))))

for (i in 1:3){
  write_tsv(qlf.genes.multiexon[[i]], paste0("DE/edgeR_DE_qlfMultiexonGenes_",body.part,"_", names(qlf.genes.multiexon)[i],".tsv")) 
}

DE.data$top.genes.multiexon$interaction <-qlf.genes.multiexon$interaction$GeneID[qlf.genes.multiexon$interaction$adjP.Value <= 0.05]
DE.data$top.genes.multiexon$season <-qlf.genes.multiexon$season$GeneID[qlf.genes.multiexon$season$adjP.Value <= 0.05]
DE.data$top.genes.multiexon$family <-qlf.genes.multiexon$family$GeneID[qlf.genes.multiexon$family$adjP.Value <= 0.05]

genes.coef.SvF <- qlf.genes.multiexon$family %>% 
  dplyr::select(GeneID, adjP.Value, Means, Ranges) %>% 
  dplyr::rename(adjP.Value_F = adjP.Value) %>% 
  left_join(dplyr::select(qlf.genes.multiexon$season, GeneID, adjP.Value, logFC)) %>% 
  mutate(DE_season = GeneID %in% DE.data$top.genes.multiexon$season,
         DE_family = GeneID %in% DE.data$top.genes.multiexon$family)

DE_summary <- list()
for (i in 1:3){
  DE_summary[[i]] <- c(
    length(qlf.genes[[i]]$GeneID),
    sum(qlf.genes[[i]]$PValue < 0.05),
    sum(qlf.genes[[i]]$FDR < 0.05),
    sum(qlf.genes[[i]]$adjP.Value < 0.05), 
    length(qlf.genes.multiexon[[i]]$GeneID),
    sum(qlf.genes.multiexon[[i]]$PValue< 0.05),
    sum(qlf.genes.multiexon[[i]]$FDR < 0.05), 
    sum(qlf.genes.multiexon[[i]]$adjP.Value < 0.05)
  )
}

DE.data$DE_summary <- data.frame(cbind(DE_summary[[2]],DE_summary[[3]], DE_summary[[1]])) %>% 
  dplyr::rename("Season" = X1, "Family" = X2, "Interaction" = X3) 
rownames(DE.data$DE_summary) <- c("Genes", "DE genes (P-value)", "DE genes (FDR)", "DE genes (adjusted)", "Multiexon Genes", "DE multiexon genes (P-value)", "DE multiexon genes (FDR)", "DE multiexon genes (adjusted)")
DE.data$DE_summary 
write_tsv(DE.data$DE_summary, paste0("DE/DE_summary_",body.part,".tsv"))


DE_joined<- list()
colnames1 <- c("GeneID", paste0("DE_pval_S_",body.part), paste0("DE_fdr_S_",body.part), paste0("DE_adjpval_S_",body.part))
colnames2 <- c(paste0("DE_pval_F_",body.part), paste0("DE_fdr_F_",body.part), paste0("DE_adjpval_F_",body.part))
colnames3 <- c(paste0("DE_pval_I_",body.part), paste0("DE_fdr_I_",body.part), paste0("DE_adjpval_I_",body.part))

DE_joined$numeric <- qlf.genes.multiexon$season %>% dplyr::select(GeneID, PValue ,FDR, adjP.Value) 
colnames(DE_joined$numeric) <- colnames1
DE_joined$numeric <- DE_joined$numeric %>% full_join(dplyr::select(qlf.genes.multiexon$family, GeneID, PValue ,FDR, adjP.Value))
colnames(DE_joined$numeric) = c(colnames1,colnames2)
DE_joined$numeric <- DE_joined$numeric %>% full_join(dplyr::select(qlf.genes.multiexon$interaction, GeneID, PValue ,FDR, adjP.Value))
colnames(DE_joined$numeric) = c(colnames1,colnames2,colnames3)

DE_joined$binary <- DE_joined$numeric %>% 
  mutate_at(vars(-GeneID), ~ if_else(.x<0.05,1,0)) %>%
  data.frame()

DE_joined$euler <- DE_joined$binary %>% mutate_at(vars(-GeneID), as.logical)

save(DE_joined, file = paste0("DE/DE_joined_",body.part))

for (i in 1:3){
  write_tsv(DE_joined[[i]], paste0("DE/DE_joined_",body.part,"_", names(DE_joined)[i],".tsv")) 
}

save(DE.data, file = paste0("DE/DE.data.", body.part))

#### Figures #### 

color2 <-c( "#0072B5", "#0072B5", "#0072B5")
color_transparent2 <- adjustcolor(color2, alpha.f = 0.75)
# Upset plots
pdf(file = paste0("DE/Overlap/edgeR_DE_upset_",body.part,".pdf"), 
    height = 4,
    width = 6)
upset(DE_joined$binary[,c(2,5,8)], order.by = "freq")
upset(DE_joined$binary[,c(3,6,9)], order.by = "freq")
upset(DE_joined$binary[,c(4,7,10)], order.by = "freq")
dev.off()

# Euler venn diagrams
pdf(file = paste0("DE/Overlap/edgeR_DE_euler_",body.part,".pdf"))
plot(euler(DE_joined$euler[, c(2,5,8)], shape = "ellipse"),
     quantities = list(fontsize =16),
     fill = color_transparent2,
     lty = 1, labels =list(labels= c("S", "F", "SxF"), font = 2, fontsize = 18))
DE_euler <- plot(euler(DE_joined$euler[, c(3,6,9)], shape = "ellipse"),
     quantities = list(fontsize =20),
     fill = color_transparent2,
     lty = 1,
     lty = 1, labels =list(labels= c("Season\nmain effect", "Family\nmain effect", "SxF\ninteraction"), font = 2, fontsize = 18),
     adjust_labels = T)
DE_euler

plot(euler(DE_joined$euler[, c(4,7,10)], shape = "ellipse"),
     quantities = list(fontsize =16),
     fill = color_transparent2,
     lty = 1, labels =list(labels= c("S", "F", "SxF"), font = 2, fontsize = 18))
dev.off() 

#MDS Plots
if (body.part=="abdomen") {
  cols = c("#D95F02", "#1B9E77")
} else {
  cols =c("#FC8D62","#66C2A5")
}
cols_season<- cols[y.genes.filt$samples$group]
pchs <- c(1,2,3,4,5,6,7)
pchs_family <- pchs[y.genes.filt$samples$family]

{pdf(file = paste0("DE/Clustering/edgeR_DE_MDSbySeasonStress_",body.part,".pdf"),
     height = 5, 
     width = 5)  
  plotMDS(y.genes.filt,  dim.plot = c(1,2), pch = pchs_family, lwd = 3, col = cols_season,
          main="Abdomen", gene.selection = "pairwise")
  legend("topleft", legend=c("Dry", "Wet"),
         col=cols, pch = 19, cex=0.8)            
  dev.off()} 

#comparison of dimensions
{pdf(file = paste0("DE/Clustering/edgeR_DE_MDSbySeasonStress_dims1to5_",body.part,".pdf"),
     height = 5.5, 
     width = 5)  
  plotMDS(y.genes.filt,  dim.plot = c(1,2), pch = pchs_family, lwd = 2.5, col = cols_season)
  plotMDS(y.genes.filt,  dim.plot = c(1,3), pch = pchs_family, lwd = 2.5, col = cols_season)
  plotMDS(y.genes.filt,  dim.plot = c(1,4), pch = pchs_family, lwd = 2.5, col = cols_season)
  plotMDS(y.genes.filt,  dim.plot = c(1,5), pch = pchs_family, lwd = 2.5, col = cols_season)
  plotMDS(y.genes.filt,  dim.plot = c(2,3), pch = pchs_family, lwd = 2.5, col = cols_season)
  plotMDS(y.genes.filt,  dim.plot = c(2,4), pch = pchs_family, lwd = 2.5, col = cols_season)
  plotMDS(y.genes.filt,  dim.plot = c(2,5), pch = pchs_family, lwd = 2.5, col = cols_season)
  plotMDS(y.genes.filt,  dim.plot = c(3,4), pch = pchs_family, lwd = 2.5, col = cols_season)
  plotMDS(y.genes.filt,  dim.plot = c(3,5), pch = pchs_family, lwd = 2.5, col = cols_season)
  plotMDS(y.genes.filt,  dim.plot = c(4,5), pch = pchs_family, lwd = 2.5, col = cols_season)   
  dev.off()} 

# PCA
tmm.genes <- calcNormFactors(DGEList(y.genes.filt$counts), method = "TMM")
pseudo_tmm.genes <- log2(cpm(tmm.genes) + 1)
top_var_pseudo_tmm.genes <-names(sort(apply(pseudo_tmm.genes, 1, var), decreasing = T))[1:2500]
t_pseudo_tmm.genes <- t(pseudo_tmm.genes[top_var_pseudo_tmm.genes,])
tmmPCA.genes <- prcomp(t_pseudo_tmm.genes)
tmmPCA_summary.genes <- as.data.frame(t((summary(tmmPCA.genes))$importance)) %>% rownames_to_column(var = "Principal Component")
tmmPCA.data.frame.genes <- as.data.frame(tmmPCA.genes$x) %>%
  rownames_to_column(var = "Run") %>% 
  left_join(dplyr::select(Metadata, Run, Season, Family, season_stress)) %>% 
  mutate(Family = as.factor(as.character(Family)),
         Season = as.factor(as.character(Season)))
write_tsv(tmmPCA.data.frame.genes,paste0("DE/edgeR_DE_tmmPCA_", body.part,".tsv"))

ggplot(tmmPCA_summary.genes[1:6,], aes( x= `Principal Component`, y = `Proportion of Variance`)) +
  geom_bar(stat = "identity") +
  labs(x = "Principal component", "Prop. of variance") +
  geom_label(stat = "identity", aes(label = round(`Proportion of Variance`, 3))) +
  theme_classic(base_size = 15)

ggsave(paste0("DE/Clustering/edgeR_DE_prcomp_propVariance_",body.part,".pdf"), height = 4.5, width = 6) 

PC_axes <- c(1,2)
PC_columns <- c(paste0("PC",PC_axes[1]), paste0("PC",PC_axes[2]))
PC_var <- c(paste0(PC_columns[1], " (", round(tmmPCA_summary.genes[PC_axes[1],3]*100, 1), "%)"), paste0(PC_columns[2], " (", round(tmmPCA_summary.genes[PC_axes[2],3]*100, 1), "%)"))

(pca.plot <- ggplot(tmmPCA.data.frame.genes, aes_(x = as.name(PC_columns[1]), y = as.name(PC_columns[2]))) +
    geom_point(size = 2, stroke = 1.25, aes(color = Season, shape= Family))+
    scale_color_manual(labels = c("Dry", "Wet"),
                       values =cols)+
    labs(x = PC_var[1], y =PC_var[2]) +
    scale_shape_manual(values = c(1,2,3,4,5,6,7), guide = FALSE) +
    theme_classic(base_size = 15) + theme(legend.position = "bottom"))
ggsave(plot = pca.plot, paste0("DE/Clustering/edgeR_DE_prcomp_",PC_columns[1], PC_columns[2],"_",body.part,".pdf"), height = 4.5, width = 4) 

## season
ngenes1 <- ggplot(qlf.genes.multiexon$season, aes(x = Length, -log10(adjP.Value))) +
  geom_point() +
  labs(x="Gene Length", y = expression(`-log`[10]* " adj. P-value"), title = "Season") +
  theme_classic()

volc1 <- ggplot(qlf.genes.multiexon$season, aes(x = logFC, y = -log10(adjP.Value))) +
  geom_point(aes(color = adjP.Value <= 0.05)) +
  scale_color_manual(values = c("gray", "black"), guide = FALSE) +
  labs(x=expression(`Relative log`[2]* " fold-change"), y = expression(`-log`[10]* " adj. P-value"), title = "Season fold change") +
  theme_classic()

bp1 <- ggplot(genes.coef.SvF, aes(x = DE_season, y = abs(logFC))) +
  geom_boxplot() +
  scale_x_discrete(labels = c("non-DE", "DE")) +
  labs(x = "Season", y = expression(`abs(Relative log`[2]* " fold-change)"), title = "Seasonal Log fold-change of genes within DE and non-DE genes")  +
  theme_classic()

bp2 <-ggplot(genes.coef.SvF, aes(x = DE_family, y = abs(logFC))) +
  geom_boxplot() +
  scale_x_discrete(labels = c("non-DE", "DE")) +
  labs(x = "Family", y = expression(`abs(Relative log`[2]* " fold-change)"), title = "Seasonal Log fold-change of genes within DE and non-DE genes")  +
  theme_classic()

## family
ngenes2 <- ggplot(qlf.genes.multiexon$family, aes(x = Length, -log10(adjP.Value))) +
  geom_point() +
  labs(x="Gene Length", y = expression(`-log`[10]* " adj. P-value"), title = "Family") +
  theme_classic()

volc2 <- ggplot(qlf.genes.multiexon$family, aes(x = Means, y = -log10(adjP.Value))) +
  geom_point(aes(color = adjP.Value <= 0.05)) +
  scale_color_manual(values = c("gray", "black"), guide = FALSE) +
  labs(x=expression(`Mean log`[2]* " fold-change"), y = expression(`-log`[10]* " adj. P-value"), title = "Family fold change (coef. means)") +
  theme_classic()

volc3 <- ggplot(qlf.genes.multiexon$family, aes(x = Ranges, y = -log10(adjP.Value))) +
  geom_point(aes(color = adjP.Value <= 0.05)) +
  scale_color_manual(values = c("gray", "black"), guide = FALSE) +
  labs(x=expression(`Range of log`[2]* " fold-changes"), y = expression(`-log`[10]* " adj. P-value"), title = "Family fold change (coef. ranges)") +
  theme_classic()

bp3 <- ggplot(genes.coef.SvF, aes(x = DE_season, y = abs(Ranges))) +
  geom_boxplot() +
  scale_x_discrete(labels = c("non-DE", "DE")) +
  labs(x = "Season", y = expression(`abs(Relative log`[2]* " fold-change)"), title = "Family Log fold-change of genes within DE and non-DE genes")  +
  theme_classic()

bp4 <- ggplot(genes.coef.SvF, aes(x = DE_family, y = abs(Ranges))) +
  geom_boxplot() +
  scale_x_discrete(labels = c("non-DE", "DE")) +
  labs(x = "Family", y = expression(`abs(Relative log`[2]* " fold-change)"), title = "Family Log fold-change of genes within DE and non-DE genes")  +
  theme_classic()

# interaction
ngenes3 <- ggplot(qlf.genes.multiexon$interaction, aes(x = Length, -log10(adjP.Value))) +
  geom_point() +
  labs(x="Gene Length", y = expression(`-log`[10]* " adj. P-value"), title = "Interaction") +
  theme_classic()

ngenes <- grid.arrange(ngenes1, ngenes2, ngenes3, nrow = 1)
ggsave(plot = ngenes, file = paste0("DE/edgeR_DE_NExonsSimesPlots_", body.part,".pdf"), height = 4, width = 12) 
volc <- grid.arrange(volc1, volc2, volc3, nrow = 1)
ggsave(plot = volc, file = paste0("DE/edgeR_DE_volcanoPlots_", body.part,".pdf"), height = 4, width = 12) 
bp <- grid.arrange( bp1, bp2, bp3, bp4, nrow = 2, ncol =2)
ggsave(plot = bp, file = paste0("DE/edgeR_DE_SandF_foldchange_", body.part,".pdf"), height = 8, width = 8) 



ab <- (read_tsv("DE/DE_joined_abdomen_euler.tsv"))[c(1,3, 6,9)]

th <- (read_tsv("DE/DE_joined_thorax_euler.tsv"))[c(1,3, 6,9)]
all <- ab%>% left_join(th) %>% replace(., is.na(.), FALSE)

color2 <-c( "#0072B5", "#0072B5", "#0072B5")
color_transparent2 <- adjustcolor(color2, alpha.f = 0.75)
pdf(file = paste0("DE/edgeR_DE_euler_abvth.pdf"))
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




