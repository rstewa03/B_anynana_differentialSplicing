#### Tissue lists ####

library(tidyverse)
library(eulerr)
library(UpSetR)
#### THORAX v ABDOMEN: ####
# run the above script on abdomen and thorax before proceeding
load("05_edgeRVJUM/dataframes/star_fc_joined_genes_abdomen")
load("05_edgeRVJUM/dataframes/star_fc_joined_genes_thorax")

joined_genes_both <- list()
joined_genes_both$binary <- joined_genes_abdomen$binary %>% 
  full_join(joined_genes_thorax$binary) %>% 
  mutate_at(vars(-GeneID), ~ if_else(is.na(.x),0,.x)) %>% 
  data.frame()
joined_genes_both$euler <- joined_genes_abdomen$euler %>% 
  full_join(joined_genes_thorax$euler) %>% 
  mutate_at(vars(-GeneID), ~ if_else(is.na(.x),FALSE,.x)) %>% 
  data.frame()

#1# $ GeneID              : chr  "gene14013" "gene9416" "gene6427" "gene12579" ...
# $ DE_pval_S_abdomen   : num  1 1 1 1 1 1 1 1 1 1 ...
# $ DE_fdr_S_abdomen    : num  1 1 1 1 1 1 1 1 1 1 ...
# $ DE_adjpval_S_abdomen: num  1 1 1 1 1 1 1 1 1 1 ...
#5 # $ DE_pval_F_abdomen   : num  1 1 1 1 1 0 0 1 1 1 ...
# $ DE_fdr_F_abdomen    : num  1 1 1 1 1 0 0 1 1 1 ...
# $ DE_adjpval_F_abdomen: num  1 1 1 1 1 0 0 1 1 1 ...
# $ DE_pval_I_abdomen   : num  1 0 0 0 0 0 0 0 0 0 ...
# $ DE_fdr_I_abdomen    : num  0 0 0 0 0 0 0 0 0 0 ...
#10 # $ DE_adjpval_I_abdomen: num  0 0 0 0 0 0 0 0 0 0 ...
# $ P.value_S_abdomen   : num  0 0 1 0 0 0 0 0 0 0 ...
# $ FDR_S_abdomen       : num  0 0 1 0 0 0 0 0 0 0 ...
# $ adjP.value_S_abdomen: num  0 0 1 0 0 0 0 0 0 0 ...
# $ P.value_F_abdomen   : num  0 0 0 0 0 0 0 0 0 0 ...
#15# $ FDR_F_abdomen       : num  0 0 0 0 0 0 0 0 0 0 ...
# $ adjP.value_F_abdomen: num  0 0 0 0 0 0 0 0 0 0 ...
# $ P.value_I_abdomen   : num  0 0 0 0 0 0 0 0 0 0 ...
# $ FDR_I_abdomen       : num  0 0 0 0 0 0 0 0 0 0 ...
# $ adjP.value_I_abdomen: num  0 0 0 0 0 0 0 0 0 0 ...
#20# $ sigJUM_abdomen       : num  0 0 0 0 0 0 0 0 0 0 ...
#  # $ sigJUM_strict_abdomen: num  0 0 0 0 0 0 0 0 0 0 ...
# $ DE_pval_S_thorax    : num  0 1 1 1 1 1 1 1 1 0 ...
# $ DE_fdr_S_thorax     : num  0 1 1 1 1 1 1 1 1 0 ...
# $ DE_adjpval_S_thorax : num  0 1 1 1 1 1 1 1 1 0 ...
#25# $ DE_pval_F_thorax    : num  0 1 1 1 1 1 1 0 1 0 ...
# $ DE_fdr_F_thorax     : num  0 1 1 1 1 1 0 0 1 0 ...
# $ DE_adjpval_F_thorax : num  0 1 1 1 1 1 0 0 1 0 ...
# $ DE_pval_I_thorax    : num  0 0 1 0 0 0 0 0 0 0 ...
# $ DE_fdr_I_thorax     : num  0 0 0 0 0 0 0 0 0 0 ...
#30# $ DE_adjpval_I_thorax : num  0 0 1 0 0 0 0 0 0 0 ...
# $ P.value_S_thorax    : num  0 0 0 0 0 0 0 0 0 0 ...
# $ FDR_S_thorax        : num  0 0 0 0 0 0 0 0 0 0 ...
# $ adjP.value_S_thorax : num  0 0 0 0 0 0 0 0 0 0 ...
# $ P.value_F_thorax    : num  0 0 0 0 0 0 0 0 0 0 ...
#35# $ FDR_F_thorax        : num  0 0 0 0 0 0 0 0 0 0 ...
# $ adjP.value_F_thorax : num  0 0 0 0 0 0 0 0 0 0 ...
# $ P.value_I_thorax    : num  0 0 0 0 0 0 0 0 0 0 ...
# $ FDR_I_thorax        : num  0 0 0 0 0 0 0 0 0 0 ...
# $ adjP.value_I_thorax : num  0 0 0 0 0 0 0 0 0 0 ...
#40# $ sigJUM_thorax        : num  0 0 0 0 0 0 0 0 0 0 ...
# $ sigJUM_strict_thorax : num  0 0 0 0 0 0 0 0 0 0 ...

color_transparent.edgerDE <- adjustcolor(c( "cornflowerblue", "cornflowerblue"), alpha.f = 0.75)
color_transparent.edgerDS <- adjustcolor(c( "firebrick", "firebrick"), alpha.f = 0.75)
color_transparent.JUM <- adjustcolor(c( "goldenrod", "goldenrod"), alpha.f = 0.75)


pdf(file = "05_edgeRVJUM/figures/edgeRDSvJUM_season_both.pdf",
    height = 3,
    width = 4.5)
upset(joined_genes_both$binary[,c(13, 21, 33,41)], order.by = "freq")
plot(euler(joined_genes_both$euler[, c(13, 33)],  shape= "circle"), 
     quantities = list(fontsize =16),
     fill = color_transparent.edgerDS,
     lty = 1, labels =list(labels= c("A", "T"), font = 2, fontsize = 18),
     main = "edgeR DS")

plot(euler(joined_genes_both$euler[, c(21,41)],  shape= "circle"), 
     quantities = list(fontsize =16),
     fill = color_transparent.JUM,
     lty = 1, labels =list(labels= c("A", "T"), font = 2, fontsize = 18),
     main = "JUM DS")

plot(euler(joined_genes_both$euler[, c(4,24)],  shape= "circle"), 
     quantities = list(fontsize =16),
     fill = color_transparent.edgerDE,
     lty = 1, labels =list(labels= c("A", "T"), font = 2, fontsize = 18),
     main = "edgeR DE")

dev.off()


#### JUM Overlap: abdomen v. thorax ####
load("/cerberus/projects/racste/B_anynana/diff_expr/05_edgeRVJUM/dataframes/star_fc_DS_JUM_abdomen")
load("/cerberus/projects/racste/B_anynana/diff_expr/05_edgeRVJUM/dataframes/star_fc_DS_JUM_thorax")

DS_JUM_ab_event <- DS_JUM_abdomen$all_events %>% 
  dplyr::select(GeneID, AS_event_ID, event, sigJUM_strict_abdomen) %>% 
  mutate(isAB = 1)

DS_JUM_th_event <- DS_JUM_thorax$all_events %>% 
  dplyr::select(GeneID, AS_event_ID, event, sigJUM_strict_thorax) %>% 
  mutate(isTH = 1)

DS_JUM_ab_genes <- DS_JUM_abdomen$unique_genes %>% 
  select(GeneID, sigJUM_strict_abdomen) %>% 
  mutate(isAB = 1) 

DS_JUM_th_genes <- DS_JUM_thorax$unique_genes %>% 
  select(GeneID, sigJUM_strict_thorax) %>% 
  mutate(isTH = 1) 

DS_JUM_both_events <- DS_JUM_ab_event %>%  full_join(DS_JUM_th_event) %>% 
#sum(is.na(DS_JUM_both_events$isAB)) == 82
#sum(is.na(DS_JUM_both_events$isTH)) == 137
  mutate_at(vars(-c(GeneID, AS_event_ID, event)), ~ if_else(is.na(.x),0,.x)) %>% 
  data.frame()%>% 
  mutate_at(vars(-c(GeneID, AS_event_ID, event)), as.logical)

DS_JUM_both_genes <- DS_JUM_ab_genes %>% 
  mutate(isAB = 1) %>% 
  full_join(DS_JUM_th_genes) %>% 
  mutate_at(vars(-c(GeneID)), ~ if_else(is.na(.x),0,.x)) %>% 
  data.frame()%>% 
  mutate_at(vars(-GeneID), as.logical)

pdf(file = "05_edgeRVJUM/figures/star_B_anynana_both_JUM_events_genes.pdf",
    height = 2,
    width = 4) 
plot(euler(DS_JUM_both_events[, c(4,6)],  shape= "ellipse"), 
     quantities = list(fontsize = 15),
     labels = c("Abdomen", "Thorax"),
     main = "Significant splice events")

plot(euler(DS_JUM_both_events[, c(5,7)],  shape= "ellipse"), 
     quantities = list(fontsize = 15),
     labels = c("Abdomen", "Thorax"),
     main = "All splice events")

plot(euler(DS_JUM_both_genes[, c(2,4)],  shape= "ellipse"), 
     quantities = list(fontsize = 15),
     labels = c("Abdomen", "Thorax"), 
     main = "Genes containing significant splice events")

plot(euler(DS_JUM_both_genes[, c(3,5)],  shape= "ellipse"), 
     quantities = list(fontsize = 15),
     labels = c("Abdomen", "Thorax"),
     main = "Genes containing splice events")
dev.off()



