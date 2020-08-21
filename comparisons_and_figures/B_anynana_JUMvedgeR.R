library(tidyverse)

# JUM v. EdgeR
setwd("/cerberus/projects/racste/B_anynana/diff_expr/")

body.part <- "thorax" #abdomen, thorax
comparison <- "TH_DRY_vs_TH_WET" #AB_DRY_vs_AB_WET, TH_DRY_vs_TH_WET

#### LOAD DATA ####
DS_JUM <- list()
events <- c("A3SS", "A5SS", "SE", "Comp", "IR", "MXE")

DS_JUM$all_events <- read_tsv(paste0("04_JUM_S/", comparison,"/AS_differential_JUM_all_splice_events_simplified_", comparison,".tsv")) %>% 
  dplyr::rename("GeneID" = "Gene") %>% 
  filter(GeneID != 'NONE') %>% 
  mutate(sigJUM = if_else(pvalue < 0.05, 1, 0), 
         sigJUM_strict = if_else(qvalue < 0.05 & max_dpsi > 0.05, 1,0),
         event = factor(event, levels = events)) 

DS_JUM$unique_genes <-  DS_JUM$all_events %>% 
  dplyr::select(GeneID, sigJUM, sigJUM_strict) %>% 
  group_by(GeneID) %>% 
  summarise(sigJUM= max(sigJUM),
            sigJUM_strict= max(sigJUM_strict))

load(paste0("01_edgeR_exons_SxF/02_edgeR_exons/DS_dataframes/DS_joined_",body.part))
load(paste0("02_edgeR_genes_SxF/02_edgeR_genes/DE_dataframes/DE_joined_",body.part))


joined_genes <- list()
joined_genes$binary <- full_join(DE_joined$binary, DS_joined$binary) %>% 
  full_join(DS_JUM$unique_genes) %>% 
  mutate_at(vars(-GeneID), ~ if_else(is.na(.x),0,.x)) %>% 
  data.frame()
joined_genes$euler <-  joined_genes$binary %>% 
  mutate_at(vars(-GeneID), as.logical)


#1  $ GeneID         : chr  "gene14013" "gene9416" "gene6427" "gene12579" ...
#2  $ P.value_S_DE: num  1 1 1 1 1 1 1 1 1 1 ...
#3  $ FDR_S_DE    : num  1 1 1 1 1 1 1 1 1 1 ...
#4  $ adjPVal_S_DE: num  1 1 1 1 1 1 1 1 1 1 ...
#5  $ P.value_F_DE: num  1 1 1 1 1 0 0 1 1 1 ...
#6  $ FDR_F_DE    : num  1 1 1 1 1 0 0 1 1 1 ...
#7  $ adjPVal_F_DE: num  1 1 1 1 1 0 0 1 1 1 ...
#8  $ P.value_I_DE: num  1 0 0 0 0 0 0 0 0 0 ...
#9  $ FDR_I_DE    : num  0 0 0 0 0 0 0 0 0 0 ...
#10 $ adjPVal_I_DE: num  0 0 0 0 0 0 0 0 0 0 ...
#11 $ P.value_S   : num  0 0 1 0 0 0 0 0 0 0 ...
#12 $ FDR_S       : num  0 0 1 0 0 0 0 0 0 0 ...
#13 $ adjP.value_S: num  0 0 1 0 0 0 0 0 0 0 ...
#14 $ P.value_F   : num  0 0 0 0 0 0 0 0 0 0 ...
#15 $ FDR_F       : num  0 0 0 0 0 0 0 0 0 0 ...
#16 $ adjP.value_F: num  0 0 0 0 0 0 0 0 0 0 ...
#17 $ P.value_I   : num  0 0 0 0 0 0 0 0 0 0 ...
#18 $ FDR_I       : num  0 0 0 0 0 0 0 0 0 0 ...
#19 $ adjP.value_I: num  0 0 0 0 0 0 0 0 0 0 ...
#20 $ sigJUM         : num  0 0 0 0 0 0 0 0 0 0 ...
#21 $ sigJUM_strict  : num  0 0 0 0 0 0 0 0 0 0 ...

JUM_vars <- c(paste0("sigJUM_", body.part), paste0("sigJUM_strict_", body.part))
names(DS_JUM$all_events)[c(10,11)] <- JUM_vars
names(DS_JUM$unique_genes)[c(2,3)] <- JUM_vars
names(joined_genes$binary)[c(20,21)] <- JUM_vars
names(joined_genes$euler)[c(20,21)] <- JUM_vars

if (body.part == "abdomen") {
  DS_JUM_abdomen <- list()
  DS_JUM_abdomen <- DS_JUM
  save(DS_JUM_abdomen, file = paste0("05_edgeRVJUM/dataframes/star_fc_DS_JUM_abdomen"))
  
  joined_genes_thorax <- list()
  joined_genes_abdomen <- joined_genes
  save(joined_genes_abdomen, file = paste0("05_edgeRVJUM/dataframes/star_fc_joined_genes_abdomen"))
  
} else if (body.part == "thorax"){
  DS_JUM_thorax <- list()
  DS_JUM_thorax <- DS_JUM
  save(DS_JUM_thorax, file = paste0("05_edgeRVJUM/dataframes/star_fc_DS_JUM_thorax"))
  joined_genes_thorax <- list()
  joined_genes_thorax <- joined_genes
  save(joined_genes_thorax, file = paste0("05_edgeRVJUM/dataframes/star_fc_joined_genes_thorax"))
}

#### Overlap Figures  ####

color3 <-c( "cornflowerblue", "goldenrod", "firebrick")
color_transparent3 <- adjustcolor(color3, alpha.f = 0.75)

pdf(file = paste0("05_edgeRVJUM/figures/edgeRDEvedgeRDSvJUM_season_", comparison,".pdf"),
    height = 4,
    width = 6)
upset(joined_genes$binary[,c(4,13,21)], order.by = "freq")
plot(euler(joined_genes$euler[, c(4,21,13)],  shape= "ellipse"), 
     quantities = list(fontsize =16),
     fill = color_transparent3,
     lty = 1, labels =list(labels= c("DE", "DS JUM", "DS edgeR"), font = 2, fontsize = 18))
dev.off()


#### Splicing event figures ####
event_overlap <- list()
event_overlap$all_events <- DS_JUM$all_events[,c(1,2,11)] %>% 
  full_join(DS_joined$binary[,c(1, 3)]) %>% 
  mutate(event = as.factor(event)) %>% 
  filter(GeneID != "NONE" & event != "NA") 
colnames(event_overlap$all_events) <- c("event", "GeneID", "JUM", "edgeR")

event_overlap$summary<- event_overlap$all_events%>% 
  group_by(event) %>% 
  mutate(n_events = n()) %>% 
  ungroup() %>% 
  group_by(event, JUM, edgeR) %>% 
  summarise(n_JUM_EdgeR = n()) %>% ungroup() %>% na.omit() %>% 
  mutate(type = factor(paste0(JUM, "_", edgeR), levels = c("0_0", "0_1", "1_1", "1_0")),
         type2 = as.factor(if_else(type == "0_1", "0_0", as.character(type))))
event_overlap$summary$type2 <- factor(event_overlap$summary$type2, levels = c("0_0", "1_1", "1_0"))

ggplot(event_overlap$summary, aes(x = event, y = n_JUM_EdgeR, fill = type)) +
  geom_bar(position = "stack", stat= "identity", alpha = 0.75) +
  labs(x = "Splice event", y = "Count") +
  scale_x_discrete(labels= c("A3SS", "A5SS", "SE", "Comp.", "IR", "MXE")) +
  scale_fill_manual( name = "DS Results", 
                     values = c("lightgray", "firebrick", "darkorange3", "goldenrod"),
                     labels =c("NS",  "NS JUM, Sig. edgeR", "Sig. JUM, Sig. edgeR", "Sig. JUM, NS edgeR")) +
  theme_classic(base_size = 15)+
  theme(panel.grid = element_blank())
ggsave(paste0("05_edgeRVJUM/figures/JUMspliceEvents_bar_", comparison,".pdf"), width = 6, height = 4)

