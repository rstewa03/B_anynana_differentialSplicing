setwd("/cerberus/projects/racste/B_anynana/diff_expr/")
library(tidyverse)
library(UpSetR)
library(eulerr)
#library(cowplot) 
#library(lme4)
# library(GeneOverlap)#
#abdomen dataframes
load("01_edgeR_exons_SxF/02_edgeR_exons/DS_dataframes/DS_joined_abdomen")
load("02_edgeR_genes_SxF/02_edgeR_genes/DE_dataframes/DE_joined_abdomen")
body.part <- "abdomen"


## thorax dataframes
# load("01_edgeR_exons_SxF/02_edgeR_exons/DS_dataframes/DS_joined_thorax")
# load("02_edgeR_genes_SxF/02_edgeR_genes/DE_dataframes/DE_joined_thorax")
# body.part <- "thorax"

adj <- list()
adj$exons_S <- DS_joined$binary$GeneID[DS_joined$binary[4] == 1]
adj$exons_F <- DS_joined$binary$GeneID[DS_joined$binary[7] == 1]
adj$exons_I <- DS_joined$binary$GeneID[DS_joined$binary[10] == 1]
adj$genes_S <- DE_joined$binary$GeneID[DE_joined$binary[4] == 1]
adj$genes_F <- DE_joined$binary$GeneID[DE_joined$binary[7] == 1]
adj$genes_I <- DE_joined$binary$GeneID[DE_joined$binary[10] == 1]

#Combined dataframe
edgeR_joined <- list()
edgeR_joined$numeric <- DS_joined$numeric %>% 
  full_join(DE_joined$numeric) %>% 
  data.frame()
edgeR_joined$binary <- DS_joined$binary %>% 
  full_join(DE_joined$binary) %>%
  mutate_at(vars(-GeneID), ~ if_else(is.na(.x),0,.x)) %>% 
  data.frame()
edgeR_joined$euler <- DS_joined$binary %>% 
  full_join(DE_joined$binary) %>% 
  mutate_at(vars(-GeneID), as.logical) %>% data.frame()

color2.edger <-c( "cornflowerblue", "firebrick")
color_transparent2.edger <- adjustcolor(color2.edger, alpha.f = 0.75)

#Upset
# Season main effects, DS vs DE comparisons
pdf(file = paste("03_edgeR_exonsVgenes/DSvDE_figures/Overlap/edgeR_DSvDE_season_upset", body.part,".pdf", sep = "_"),
    height = 2.66,
    width = 4)
upset(edgeR_joined$binary[c(2,11)], order.by = "freq")
upset(edgeR_joined$binary[c(3,12)], order.by = "freq")
upset(edgeR_joined$binary[c(4,13)], order.by = "freq")
dev.off()

# Family main effects, DS vs DE comparisons
pdf(file = paste("03_edgeR_exonsVgenes/DSvDE_figures/Overlap/edgeR_DSvDE_family_upset", body.part,".pdf", sep = "_"),
    height = 2.66,
    width = 4)
upset(edgeR_joined$binary[c(5,14)], order.by = "freq")
upset(edgeR_joined$binary[c(6,15)], order.by = "freq")
upset(edgeR_joined$binary[c(7,16)], order.by = "freq")
dev.off()

# Interaction main effects, DS vs DE comparisons
pdf(file = paste("03_edgeR_exonsVgenes/DSvDE_figures/Overlap/edgeR_DSvDE_interaction_upset", body.part,".pdf", sep = "_"),
    height = 2.66,
    width = 4)
upset(edgeR_joined$binary[c(8,17)], order.by = "freq")
upset(edgeR_joined$binary[c(9,18)], order.by = "freq")
upset(edgeR_joined$binary[c(10,19)], order.by = "freq")
dev.off()

pdf(file = paste("03_edgeR_exonsVgenes/DSvDE_figures/Overlap/edgeR_DSvDE_season_euler", body.part,".pdf", sep = "_"),
    height = 2.66,
    width = 4)
plot(euler(edgeR_joined$binary[, c(11,2)]), quantities = list(fontsize =15), lty = 1, fills = color_transparent2.edger, labels =c("DE", "DS"))
plot(euler(edgeR_joined$binary[, c(12,3)]),quantities = list(fontsize =15), lty = 1, fills = color_transparent2.edger, labels =c("DE", "DS"))
plot(euler(edgeR_joined$binary[, c(13,4)]),quantities = list(fontsize =15), lty = 1, fills = color_transparent2.edger, labels =c("DE", "DS"))
dev.off()

pdf(file = paste("03_edgeR_exonsVgenes/DSvDE_figures/Overlap/edgeR_DSvDE_family_euler", body.part,".pdf", sep = "_"),
    height = 2.66,
    width = 4)
plot(euler(edgeR_joined$binary[, c(14,5)]), quantities = list(fontsize =15), lty = 1, fills = color_transparent2.edger, labels =c("DE", "DS"))
plot(euler(edgeR_joined$binary[, c(15,6)]),quantities = list(fontsize =15), lty = 1, fills = color_transparent2.edger, labels =c("DE", "DS"))
plot(euler(edgeR_joined$binary[, c(16,7)]),quantities = list(fontsize =15), lty = 1, fills = color_transparent2.edger, labels =c("DE", "DS"))
dev.off()

pdf(file = paste("03_edgeR_exonsVgenes/DSvDE_figures/Overlap/edgeR_DSvDE_interaction_euler", body.part,".pdf", sep = "_"),
    height = 2.66,
    width = 4)
plot(euler(edgeR_joined$binary[, c(17,8)]), quantities = list(fontsize =15), lty = 1, fills = color_transparent2.edger, labels =c("DE", "DS"))
plot(euler(edgeR_joined$binary[, c(18,9)]),quantities = list(fontsize =15), lty = 1, fills = color_transparent2.edger, labels =c("DE", "DS"))
plot(euler(edgeR_joined$binary[, c(19,10)]),quantities = list(fontsize =15), lty = 1, fills = color_transparent2.edger, labels =c("DE", "DS"))
dev.off()


#### LogFC ####
library(lme4)
library(cowplot)
library(scales)

joined_logfc <- list()
exons_S_subset <- read_tsv(paste0("01_edgeR_exons_SxF/02_edgeR_exons/DS_dataframes/edgeR_DS_exons_subset_",body.part,"_season.tsv")) %>% 
  dplyr::select(GeneID, ExonID, logFC, adjP.Value) %>% 
  dplyr::rename(logFC_exon = logFC, adjP.Value_exon = adjP.Value)

joined_logfc$S_subset <- read_tsv(paste0("02_edgeR_genes_SxF/02_edgeR_genes/DE_dataframes/edgeR_DE_qlfMultiexonGenes_",body.part,"_season.tsv")) %>%
  right_join(exons_S_subset) %>% 
  mutate(DSG = GeneID %in% adj$exons_S,
         DEG = GeneID %in% adj$genes_S, 
         DEDS = factor(paste(DSG, DEG, sep = "_"), levels = c("TRUE_TRUE", "TRUE_FALSE", "FALSE_TRUE", "FALSE_FALSE"))) ; rm(exons_S_subset)

lmm1 <- lmer(logFC_exon ~ logFC + (1|GeneID), data = joined_logfc$S_subset)
lmm2 <- lmer(logFC_exon ~ 1 + (1|GeneID), data = joined_logfc$S_subset)
(lmm_anova <- anova(lmm1, lmm2))

mm<-model.matrix(~logFC, data = joined_logfc$S_subset)
y<-mm%*%fixef(lmm1)
pvar1 <- diag(mm %*% tcrossprod(vcov(lmm1),mm))
joined_logfc$S_subset <- data.frame(joined_logfc$S_subset,
                             y = y,
                             plo =y - 1.96*sqrt(pvar1),
                             phi = y + 1.96*sqrt(pvar1))

main <- ggplot(joined_logfc$S_subset, aes(x = logFC, y = logFC_exon), size = .5) +
  geom_point(aes(color = DEDS), alpha = 0.25) +
  geom_line(aes(x = logFC, y = y), color = "black", linetype = "solid") +
  geom_line(aes(x = logFC, y = phi), color = "black", linetype = "dashed") +
  geom_line(aes(x = logFC, y = plo), color = "black", linetype = "dashed") +
  annotate(geom = "text", x = 1.5, y = 1.5, label = paste( "P = ", scientific(lmm_anova$`Pr(>Chisq)`[2]))) +
  scale_color_manual( values = c("palevioletred4","firebrick"), guide = FALSE) +
  labs(x = expression(`Whole gene log`[2]* " fold-change"), y = expression(`Exon relative log`[2]* " fold-change")) +
  theme_classic(base_size = 20)

logfc_dens <- axis_canvas(main, axis = "x") +
  # geom_density(data=joined_foldchange_ab_season, aes(x=logFC), alpha=0.7, size=1, fill = "black") +
  geom_boxplot(data=joined_logfc$S_subset, aes(x=logFC, y = 1), width = 0.2, outlier.shape = 1) 
 # stat_summary(data=joined_logfc$S_subset,  aes(x = logFC, y=1), fun.data = "mean_cl_boot", colour = "red", size = 0.25)

logfc_exon_dens <- axis_canvas(main, axis = "y", coord_flip = TRUE) +
  # geom_density(data=joined_foldchange_ab_season, aes(x=logFC_exon), alpha=0.7, size=1, fill = "black") +
  geom_boxplot(data=joined_logfc$S_subset, aes(x=logFC_exon, y = 1), width = 0.2,  outlier.shape = 1) +
 # stat_summary(data=joined_logfc$S_subset,  aes(x = logFC_exon, y=1), fun.data = "mean_cl_boot", colour = "red", size = 0.25)+
  coord_flip()

p1 <- insert_xaxis_grob(main, logfc_dens, grid::unit(.1, "null"), position = "top")
p2 <- insert_yaxis_grob(p1, logfc_exon_dens, grid::unit(.1, "null"), position = "right")
(DEDSplot <- ggdraw(p2))
ggsave(plot = DEDSplot, paste0("03_edgeR_exonsVgenes/DSvDE_figures/edgeR_DSvDE_logFCscatterplot_S_subset", body.part,".pdf"), height = 5, width = 5)

ggplot(joined_logfc$S_subset, aes(x = DEG, y = logFC_exon)) +
  geom_boxplot() +
  labs(x = "Season DE", y = expression(`Exon relative log`[2]* " fold-change")) +
  scale_x_discrete(labels = c("non-DE", "DE")) +
  theme_classic()
ggsave(paste0("03_edgeR_exonsVgenes/DSvDE_figures/edgeR_DSvDE_boxplot_S_subset", body.part,".pdf"), height = 5, width = 5)

## all genes ##
exons_S_all <- read_tsv(paste0("01_edgeR_exons_SxF/02_edgeR_exons/DS_dataframes/edgeR_DS_exons_",body.part,"_season.tsv")) %>% 
  dplyr::select(GeneID, ExonID, logFC, adjP.Value) %>% 
  dplyr::rename(logFC_exon = logFC, adjP.Value_exon = adjP.Value)

joined_logfc$S_all <- read_tsv(paste0("02_edgeR_genes_SxF/02_edgeR_genes/DE_dataframes/edgeR_DE_qlfMultiexonGenes_",body.part,"_season.tsv")) %>%
  right_join(exons_S_all) %>% 
  mutate(DSG = GeneID %in% adj$exons_S,
         DEG = GeneID %in% adj$genes_S, 
         DEDS = factor(paste(DSG, DEG, sep = "_"), levels = c("TRUE_TRUE", "TRUE_FALSE", "FALSE_TRUE", "FALSE_FALSE"))) %>% 
  arrange(desc(DEDS)); rm(exons_S_all)

lmm1 <- lmer(logFC_exon ~ logFC + (1|GeneID), data = na.omit(joined_logfc$S_all))
lmm2 <- lmer(logFC_exon ~ 1 + (1|GeneID), data = na.omit(joined_logfc$S_all))
(lmm_anova <- anova(lmm1, lmm2))

summary(joined_logfc$S_all$logFC)
main <- ggplot(joined_logfc$S_all[joined_logfc$S_all$DEDS != "FALSE_FALSE",], aes(x = logFC, y = logFC_exon), size = 0.25) +
  # geom_point(data = joined_logfc$S_all[joined_logfc$S_all$DEDS == "FALSE_FALSE",], aes(x = logFC, y = logFC_exon), color = "black", alpha = 0.1) +
  geom_point(aes(color = DEDS), alpha = 0.25) +
  scale_color_manual(values = c("palevioletred4", "firebrick", "cornflowerblue"), guide = FALSE) +
  labs(x = expression(`Whole gene log`[2]* " fold-change"), y = expression(`Exon relative log`[2]* " fold-change")) +
  theme_classic(base_size = 20)

logfc_dens <- axis_canvas(main, axis = "x") +
  # geom_density(data=joined_foldchange_ab_season, aes(x=logFC), alpha=0.7, size=1, fill = "black") +
  geom_boxplot(data=joined_logfc$S_all, aes(x=logFC, y = 1), width = 0.2, outlier.shape = 1) 
  # stat_summary(data=joined_logfc$S_all,  aes(x = logFC, y=1), fun.data = "mean_cl_boot", colour = "red", size = 0.25)

logfc_exon_dens <- axis_canvas(main, axis = "y", coord_flip = TRUE) +
  # geom_density(data=joined_foldchange_ab_season, aes(x=logFC_exon), alpha=0.7, size=1, fill = "black") +
  geom_boxplot(data=joined_logfc$S_all, aes(x=logFC_exon, y = 1),  width = 0.2, outlier.shape = 1) +
  # stat_summary(data=joined_logfc$S_all,  aes(x = logFC_exon, y=1), fun.data = "mean_cl_boot", colour = "red", size = 0.25)+
  coord_flip()

p1 <- insert_xaxis_grob(main, logfc_dens, grid::unit(.1, "null"), position = "top")
p2 <- insert_yaxis_grob(p1, logfc_exon_dens, grid::unit(.1, "null"), position = "right")
(DEDSplot <- ggdraw(p2))

ggsave(plot = DEDSplot, paste0("03_edgeR_exonsVgenes/DSvDE_figures/edgeR_DSvDE_logFCscatterplot_S_all", body.part,".pdf"), height = 5, width = 5)


exons_F_subset <- read_tsv(paste0("01_edgeR_exons_SxF/02_edgeR_exons/DS_dataframes/edgeR_DS_exons_subset_",body.part,"_family.tsv")) %>% 
  dplyr::select(GeneID, ExonID, Ranges, adjP.Value) %>% 
  dplyr::rename(Ranges_exon = Ranges,  adjP.Value_exon = adjP.Value)

length(unique(exons_F_subset$GeneID))
joined_logfc$F_subset <- read_tsv(paste0("02_edgeR_genes_SxF/02_edgeR_genes/DE_dataframes/edgeR_DE_qlfMultiexonGenes_",body.part,"_family.tsv")) %>%
  right_join(exons_F_subset) %>% 
  mutate(DSG = GeneID %in% adj$exons_F,
         DEG = GeneID %in% adj$genes_F, 
         DEDS = factor(paste(DSG, DEG, sep = "_"), levels = c("TRUE_TRUE", "TRUE_FALSE", "FALSE_TRUE", "FALSE_FALSE"))); rm(exons_F_subset)

exons_F_all <- read_tsv(paste0("01_edgeR_exons_SxF/02_edgeR_exons/DS_dataframes/edgeR_DS_exons_",body.part,"_family.tsv")) %>% 
  dplyr::select(GeneID, ExonID, Ranges, adjP.Value) %>% 
  dplyr::rename(Ranges_exon = Ranges,adjP.Value_exon = adjP.Value)

joined_logfc$F_all <- read_tsv(paste0("02_edgeR_genes_SxF/02_edgeR_genes/DE_dataframes/edgeR_DE_qlfMultiexonGenes_",body.part,"_family.tsv")) %>%
  right_join(exons_F_all) %>% 
  mutate(DSG = GeneID %in% adj$exons_F,
         DEG = GeneID %in% adj$genes_F, 
         DEDS = factor(paste(DSG, DEG, sep = "_"), levels = c("TRUE_TRUE", "TRUE_FALSE", "FALSE_TRUE", "FALSE_FALSE"))) %>% 
  arrange(desc(DEDS)); rm(exons_F_all)

main <- ggplot(joined_logfc$F_subset, aes(x = Ranges, y = Ranges_exon), size = 0.5) +
  geom_point(aes(color = DEDS), alpha = 0.25) +
  scale_color_manual(values = c("firebrick" ,"palevioletred4", "cornflowerblue", "black"), guide = FALSE) +
  labs(x = expression(`Whole gene abs. log`[2]* " fold-change"), y = expression(`Exon relative abs. log`[2]* " fold-change")) +
  theme_classic(base_size =20)

logfc_dens <- axis_canvas(main, axis = "x") +
  # geom_density(data=joined_foldchange_ab_season, aes(x=logFC), alpha=0.7, size=1, fill = "black") +
  geom_boxplot(data=joined_logfc$F_subset, aes(x=Ranges, y = 1), width = 0.2, outlier.shape = 1) 
# stat_summary(data=joined_logfc$F_subset,  aes(x = logFC, y=1), fun.data = "mean_cl_boot", colour = "red", size = 0.25)

logfc_exon_dens <- axis_canvas(main, axis = "y", coord_flip = TRUE) +
  # geom_density(data=joined_foldchange_ab_season, aes(x=logFC_exon), alpha=0.7, size=1, fill = "black") +
  geom_boxplot(data=joined_logfc$F_subset, aes(x=Ranges_exon, y = 1), width = 0.2,  outlier.shape = 1) +
  # stat_summary(data=joined_logfc$F_subset,  aes(x = logFC_exon, y=1), fun.data = "mean_cl_boot", colour = "red", size = 0.25)+
  coord_flip()

p1 <- insert_xaxis_grob(main, logfc_dens, grid::unit(.1, "null"), position = "top")
p2 <- insert_yaxis_grob(p1, logfc_exon_dens, grid::unit(.1, "null"), position = "right")
DEDSplot <- ggdraw(p2)
ggsave(plot = DEDSplot, paste0("03_edgeR_exonsVgenes/DSvDE_figures/edgeR_DSvDE_logFCscatterplot_F_subset", body.part,".pdf"), height = 5, width = 5)

main <- ggplot(joined_logfc$F_all[joined_logfc$F_all$DEDS != "FALSE_FALSE",], aes(x = Ranges, y = Ranges_exon), size = 0.5) +
  geom_point(aes(color = DEDS), alpha = 0.25) +
   scale_color_manual( values = c("palevioletred4","firebrick", "cornflowerblue", "black"), guide = FALSE) +
  labs(x = expression(`Whole gene log`[2]* " fold-change"), y = expression(`Exon relative log`[2]* " fold-change")) +
  theme_classic(base_size = 20)

logfc_dens <- axis_canvas(main, axis = "x") +
  # geom_density(data=joined_foldchange_ab_season, aes(x=logFC), alpha=0.7, size=1, fill = "black") +
  geom_boxplot(data=joined_logfc$F_all, aes(x=Ranges, y = 1), width = 0.2, outlier.shape = 1) 
# stat_summary(data=joined_logfc$F_all,  aes(x = logFC, y=1), fun.data = "mean_cl_boot", colour = "red", size = 0.25)

logfc_exon_dens <- axis_canvas(main, axis = "y", coord_flip = TRUE) +
  # geom_density(data=joined_foldchange_ab_season, aes(x=logFC_exon), alpha=0.7, size=1, fill = "black") +
  geom_boxplot(data=joined_logfc$F_all, aes(x=Ranges_exon, y = 1),  width = 0.2, outlier.shape = 1) +
  # stat_summary(data=joined_logfc$F_all,  aes(x = logFC_exon, y=1), fun.data = "mean_cl_boot", colour = "red", size = 0.25)+
  coord_flip()


p1 <- insert_xaxis_grob(main, logfc_dens, grid::unit(.1, "null"), position = "top")
p2 <- insert_yaxis_grob(p1, logfc_exon_dens, grid::unit(.1, "null"), position = "right")
DEDSplot <- ggdraw(p2)

ggsave(plot = DEDSplot, paste0("03_edgeR_exonsVgenes/DSvDE_figures/edgeR_DSvDE_logFCscatterplot_F_all", body.part,".pdf"), height = 5, width = 5)

ggplot(joined_logfc$F_subset, aes(x = DEG, y = Ranges_exon)) +
  geom_boxplot() +
  labs(x = "Family DE", y = expression(`Exon abs. relative log`[2]* " fold-change")) +
  scale_x_discrete(labels = c("non-DE", "DE")) +
  theme_classic()
ggsave(paste0("03_edgeR_exonsVgenes/DSvDE_figures/edgeR_DSvDE_boxplot_F_subset", body.part,".pdf"), height = 5, width = 5)


#### Gene overlap ####



