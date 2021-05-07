library(tidyverse) 
library(matrixStats)
library(scales)
library(ggrepel)

'%ni%' <- Negate('%in%')
dir.list <- c('AB_DRY_vs_AB_WET',
              'TH_DRY_vs_TH_WET',
              'AB_DRY_vs_TH_DRY',
              'AB_WET_vs_TH_WET',
              'AB_vs_TH')

for (x in 1:length(dir.list)) {
  dir <- dir.list[x]
  if (dir == "AB_DRY_vs_AB_WET" | dir == "AB_DRY_vs_TH_DRY"){
    body.n = 69; t1 = 34; t2 = 68; t3 = 103; t4 = 138} else if (dir == "TH_DRY_vs_TH_WET" | dir == "AB_WET_vs_TH_WET") {
      body.n = 70; t1 = 35;  t2 = 70; t3 = 105; t4 = 140} else {
        body.n = 139; t1 = 69;  t2 = 138; t3 = 208; t4 = 278}
  print(dir)
  wd <- paste0("/cerberus/projects/racste/B_anynana/diff_expr/08_rMATS_S/rMATS_",dir, "/")
  setwd(wd[[1]])
  temp1 = list.files(pattern="*.MATS.JC.txt")
  rMATS_all = lapply(temp1, read_tsv)
  events <- c("A3SS", "A5SS", "MXE", "IR", "SE")
  names(rMATS_all) <- events
  
  for (i in 1:5) {
    rMATS_all[[i]] <- rMATS_all[[i]] %>% mutate(event = events[i], ID = paste0(event,ID), ID_1 = paste0(event,ID_1), comp = dir) %>% data.frame()
  }
  threshold <- seq(from = 0, to = 200, by = 5)
  threshold[1] <- 1
  individual <- seq(from = 0, to = 33, by = 3)
  individual[1] <- 1
  keep.genes <- list()
  keep.genes.threshold <- data.frame(threshold)
  keep.genes.individual <- data.frame(individual)
  rMATS_keep.threshold <- list()
  rMATS_keep.individual <- list()
  for (i in 1:length(rMATS_all)){
    all_sep <- rMATS_all[[i]] %>% select(IJC_SAMPLE_1, IJC_SAMPLE_2, SJC_SAMPLE_1, SJC_SAMPLE_2) %>% 
      separate(IJC_SAMPLE_1, sep = ",", into = c(paste0("S1_I_",seq(1:t1)))) %>% 
      separate(SJC_SAMPLE_1, sep = ",", into = c(paste0("S1_S_",seq(1:t1)))) %>%  
      separate(IJC_SAMPLE_2, sep = ",", into = c(paste0("S2_I_",seq(1:(t3-t2))))) %>% 
      separate(SJC_SAMPLE_2, sep = ",", into = c(paste0("S2_S_",seq(1:(t3-t2))))) %>% 
      mutate(across(everything(), as.numeric))
    
    temp1 <- list()
    temp2<- list()
    temp3 <- list()
    temp4 <-list()
    for (j in 1:length(threshold)) { 
      temp1[j] <- sum((rowSums(all_sep[1:t1]) >= threshold[j] | rowSums(all_sep[(1+t1):t2]) >= threshold[j]) &
                        (rowSums(all_sep[(1+t2):t3]) >= threshold[j] | rowSums(all_sep[(1+t3):t4]) >= threshold[j]))
      temp2[[j]] <- sum((rowSums(all_sep[1:t1]) >= threshold[j] | rowSums(all_sep[(1+t1):t2]) >= threshold[j]) &
        (rowSums(all_sep[(1+t2):t3]) >= threshold[j] | rowSums(all_sep[(1+t3):t4]) >= threshold[j]))
    }
    temp3 <- list()
    temp4 <-list()
    
    for (k in 1:length(individual)) { 
      temp3[k] <- sum((rowSums(all_sep[1:t1]>=5) >= individual[k] | rowSums(all_sep[(1+t1):t2]>=5) >= individual[k]) &
                        (rowSums(all_sep[(1+t2):t3]>=5) >= individual[k] | rowSums(all_sep[(1+t3):t4]>=5) >= individual[k] ))
      temp4[[k]] <- (rowSums(all_sep[1:t1]>=5) >= individual[k] | rowSums(all_sep[(1+t1):t2]>=5) >= individual[k]) &
        (rowSums(all_sep[(1+t2):t3]>=5) >= individual[k] | rowSums(all_sep[(1+t3):t4]>=5) >= individual[k] )
    }
    names(temp2) <- as.character(threshold)
    names(temp4) <- as.character(individual)
    
    keep.genes.threshold <- cbind(keep.genes.threshold,unlist(temp1))
    keep.genes.individual<- cbind(keep.genes.individual,unlist(temp3))
    
    rMATS_keep.threshold[[i]] <- temp2[c(1,5,11,21,41)]
    rMATS_keep.individual[[i]] <- temp4[c(1,2,6,12)]
  }
  
  names(keep.genes.threshold)[2:6] <- events
  gather(keep.genes.threshold,-threshold, key = event, value = count) %>% 
    mutate(event = factor(event, levels = c("A3SS", "A5SS", "SE", "IR", "MXE"))) %>% 
    ggplot(aes(x = threshold, y = count, color = event)) + 
    geom_line(size = 1.5) +
    geom_vline(xintercept = 20, linetype = "dashed") +
    labs(x = "Supporting read threshold", y = "Number of events", title = dir) +
    scale_color_manual(values = c("#1B9E77", "#D95F02", "#7570B3", "#66A61E", "#E6AB02"), 
                       name = "Event") +
    scale_x_log10() + #scale_y_log10() +
    theme_classic(base_size = 15) +
    theme(legend.position = "bottom")
  # ggsave(file = paste0("/cerberus/projects/racste/B_anynana/diff_expr/08_rMATS_S/figures/effect_of_read_threshold_",dir,".pdf"), height = 5, width = 5)
  
  names(keep.genes.individual)[2:6] <- events
  gather(keep.genes.individual,-individual, key = event, value = count) %>% 
    mutate(event = factor(event, levels = c("A3SS", "A5SS", "SE", "IR", "MXE"))) %>% 
    ggplot(aes(x = individual, y = count, color = event)) + 
    geom_line(size = 1.5) +
    geom_vline(xintercept = 33, linetype = "dashed") +
    labs(x = "Sample threshold", y = "Number of events", title = dir) +
    scale_color_manual(values = c("#1B9E77", "#D95F02", "#7570B3", "#66A61E", "#E6AB02"), 
                       name = "") +
    scale_x_log10() + #scale_y_log10() +
    theme_classic(base_size = 15) +
    theme(legend.position = "bottom")
  # ggsave(file = paste0("/cerberus/projects/racste/B_anynana/diff_expr/08_rMATS_S/figures/effect_of_sample_threshold_",dir,"2.pdf"), height = 5, width = 5)
  
  rMATS_filt20 <- list()
  for (i in 1:5) {
    temp <- rMATS_keep.threshold[[i]]$`20`
    rMATS_temp <- rMATS_all[[i]]
    rMATS_filt20[[i]] <- rMATS_temp[temp,]
  }
  names(rMATS_filt20) <- events
  
  rMATS_filt3 <- list()
  for (i in 1:5) {
    temp <- rMATS_keep.individual[[i]]$`3`
    rMATS_temp <- rMATS_all[[i]]
    rMATS_filt3[[i]] <- rMATS_temp[temp,]
  }
  names(rMATS_filt3) <- events
  save(rMATS_filt3, file = "rMATS_filt3ind_TH_DRY_vs_TH_WET.Rdata")
  
  # if (dir == "AB_DRY_vs_AB_WET") {
  #   rMATS_all_AB_DRY_vs_AB_WET <- rMATS_all
  #   save(rMATS_all_AB_DRY_vs_AB_WET, file = paste0("/cerberus/projects/racste/B_anynana/diff_expr/08_rMATS_S/rMATS_",dir, "/rMATS_all_",dir,".Rdata"))
  #   rMATS_filt20_AB_DRY_vs_AB_WET <- rMATS_filt20
  #   save(rMATS_filt20_AB_DRY_vs_AB_WET, file = paste0("/cerberus/projects/racste/B_anynana/diff_expr/08_rMATS_S/rMATS_",dir, "/rMATS_filt20_",dir,".Rdata"))
  #   
  # } else if (dir == "TH_DRY_vs_TH_WET") {
  #   rMATS_all_TH_DRY_vs_TH_WET <- rMATS_all
  #   save(rMATS_all_TH_DRY_vs_TH_WET, file = paste0("/cerberus/projects/racste/B_anynana/diff_expr/08_rMATS_S/rMATS_",dir, "/rMATS_all_",dir,".Rdata"))
  #   rMATS_filt20_TH_DRY_vs_TH_WET <- rMATS_filt20
  #   save(rMATS_filt20_TH_DRY_vs_TH_WET, file = paste0("/cerberus/projects/racste/B_anynana/diff_expr/08_rMATS_S/rMATS_",dir, "/rMATS_filt20_",dir,".Rdata"))
  #   
  # } else if (dir == "AB_DRY_vs_TH_DRY") {
  #   rMATS_all_AB_DRY_vs_TH_DRY <- rMATS_all
  #   save(rMATS_all_AB_DRY_vs_TH_DRY, file = paste0("/cerberus/projects/racste/B_anynana/diff_expr/08_rMATS_S/rMATS_",dir, "/rMATS_all_",dir,".Rdata"))
  #   rMATS_filt20_AB_DRY_vs_TH_DRY <- rMATS_filt20
  #   save(rMATS_filt20_AB_DRY_vs_TH_DRY, file = paste0("/cerberus/projects/racste/B_anynana/diff_expr/08_rMATS_S/rMATS_",dir, "/rMATS_filt20_",dir,".Rdata"))
  #   
  # } else if (dir == "AB_WET_vs_TH_WET") {
  #   rMATS_all_AB_WET_vs_TH_WET <- rMATS_all
  #   save(rMATS_all_AB_WET_vs_TH_WET, file = paste0("/cerberus/projects/racste/B_anynana/diff_expr/08_rMATS_S/rMATS_",dir, "/rMATS_all_",dir,".Rdata"))
  #   rMATS_filt20_AB_WET_vs_TH_WET <- rMATS_filt20
  #   save(rMATS_filt20_AB_WET_vs_TH_WET, file = paste0("/cerberus/projects/racste/B_anynana/diff_expr/08_rMATS_S/rMATS_",dir, "/rMATS_filt20_",dir,".Rdata"))}
  # 
  all_splice_events_filt20 <- list()
  for (i in 1:5) {
    all_splice_events_filt20 <- rbind(all_splice_events_filt20, select(rMATS_filt20[[i]], comp, event, ID, GeneID, FDR, IncLevelDifference))
  }
  splice_summary_filt20 <- all_splice_events_filt20 %>% group_by(event) %>%
    summarize(n_spl_events = n(),
              genes = length(unique(GeneID)),
              spl_event_per_gene = n_spl_events/genes,
              n_sig_spl_events = sum(as.numeric(FDR) < 0.05 & IncLevelDifference > 0.1)) %>% ungroup %>%
    mutate(total_spl_events = sum(n_spl_events) ,
           total_sig_spl_events = sum(n_sig_spl_events),
           expected_spl_events = (n_spl_events/total_spl_events)*total_sig_spl_events,
           event = factor(event, levels = c("A3SS", "A5SS", "SE", "IR", "MXE")))
  
  for (i in 1:nrow(splice_summary_filt20)){
    a <- (splice_summary_filt20$n_sig_spl_events)
    b <- (splice_summary_filt20$n_spl_events - a)
    c <- (round(splice_summary_filt20$expected_spl_events, 0))
    d <- (splice_summary_filt20$n_spl_events[i] -c)
    sig_exp <- matrix(c((a[i]), (b[i]),
                        (c[i]), (d[i])),
                      nrow = 2, byrow = TRUE,
                      dimnames = list(c("Sig", "NS"),
                                      c("Actual", "Expected")))
    sig_exp_ft <- fisher.test(sig_exp)
    splice_summary_filt20$ftest[i] <- as.numeric((fisher.test(sig_exp))$p.value)
    splice_summary_filt20$ftest_adj[i] <-  as.numeric(p.adjust((fisher.test(sig_exp))$p.value, method = "fdr", n = 5))
  }
  # write_tsv(splice_summary_filt20, paste0("/cerberus/projects/racste/B_anynana/diff_expr/08_rMATS_S/rMATS_",dir,"/read_threshold_20_summary_",dir,".tsv"))
  
  ggplot() +
    geom_bar(data = splice_summary_filt20, aes(x = event, y = n_spl_events), stat= "identity", fill = "grey") +
    geom_bar(data = splice_summary_filt20, aes(x = event, y = n_sig_spl_events, fill = event), stat= "identity", color = NA) +
    geom_boxplot(data = splice_summary_filt20, aes (x = event, y = expected_spl_events))+
    theme_classic(base_size = 15) +
    labs(x = "Splice Event", y = "Number of events") +
    theme(legend.position = "bottom") +
    scale_fill_manual(values = c("#1B9E77", "#D95F02", "#7570B3", "#66A61E", "#E6AB02"), 
                      guide = FALSE) +
    scale_x_discrete(label = c("A3SS", "A5SS", "SE", "IR.", "MXE")) 
  # ggsave(file = paste0("/cerberus/projects/racste/B_anynana/diff_expr/08_rMATS_S/figures/read_threshold_20_bar_",dir,"2.pdf"), height = 10, width = 4)
  
  
  
  all_splice_events_filt33 <- list()
  for (i in 1:5) {
    all_splice_events_filt33 <- rbind(all_splice_events_filt33, select(rMATS_filt33[[i]], comp, event, ID, GeneID, FDR, IncLevelDifference))
  }
  splice_summary_filt33 <- all_splice_events_filt33 %>% group_by(event) %>%
    summarize(n_spl_events = n(),
              genes = length(unique(GeneID)),
              spl_event_per_gene = n_spl_events/genes,
              n_sig_spl_events = sum(as.numeric(FDR) < 0.05 & IncLevelDifference > 0.1)) %>% ungroup %>%
    mutate(total_spl_events = sum(n_spl_events) ,
           total_sig_spl_events = sum(n_sig_spl_events),
           expected_spl_events = (n_spl_events/total_spl_events)*total_sig_spl_events,
           event = factor(event, levels = c("A3SS", "A5SS", "SE", "IR", "MXE"))) %>% arrange(event)
  
  for (i in 1:nrow(splice_summary_filt33)){
    a <- (splice_summary_filt33$n_sig_spl_events)
    b <- (splice_summary_filt33$n_spl_events - a)
    c <- (round(splice_summary_filt33$expected_spl_events, 0))
    d <- (splice_summary_filt33$n_spl_events[i] -c)
    sig_exp <- matrix(c((a[i]), (b[i]),
                        (c[i]), (d[i])),
                      nrow = 2, byrow = TRUE,
                      dimnames = list(c("Sig", "NS"),
                                      c("Actual", "Expected")))
    sig_exp_ft <- fisher.test(sig_exp)
    splice_summary_filt33$ftest[i] <- as.numeric((fisher.test(sig_exp))$p.value)
    splice_summary_filt33$ftest_adj[i] <-  as.numeric(p.adjust((fisher.test(sig_exp))$p.value, method = "fdr", n = 5))
  }
  # write_tsv(splice_summary_filt33, paste0("/cerberus/projects/racste/B_anynana/diff_expr/08_rMATS_S/rMATS_",dir,"/sample_threshold_33_summary_",dir,".tsv"))
  
  
  ggplot() +
    geom_bar(data = splice_summary_filt33, aes(x = event, y = n_spl_events), stat= "identity", fill = "grey") +
    geom_bar(data = splice_summary_filt33, aes(x = event, y = n_sig_spl_events, fill = event), stat= "identity", color = NA) +
  #  geom_boxplot(data = splice_summary_filt33, aes (x = event, y = expected_spl_events))+
    theme_classic(base_size = 15) +
    labs(x = "Splice Event", y = "Number of events") +
    theme(legend.position = "bottom") +
    scale_fill_manual(values = c("#1B9E77", "#D95F02", "#7570B3", "#66A61E", "#E6AB02"),
                      guide = FALSE) +
    scale_x_discrete(label = c("A3SS", "A5SS", "SE", "IR.", "MXE")) +
    scale_y_continuous(trans = "log", breaks = c(1,10, 100, 1000, 10000))
  # ggsave(file = paste0("/cerberus/projects/racste/B_anynana/diff_expr/08_rMATS_S/figures/sample_threshold_33_bar_",dir,".pdf"), height = 10, width = 4)
  
}

pie_data_20  <- splice_summary_filt20 %>% select(event, n_spl_events, n_sig_spl_events) %>% 
  gather(-c(event), key = type, value = n_events) %>% 
  mutate(event = factor(event, levels = c("IR", "A3SS", "A5SS", "SE", "MXE"))) %>% 
  arrange(type, event) %>% 
  mutate(total = if_else(pie_data$type == "n_spl_events",  1990, 90728),
         percent = n_events/total,
         label = if_else(pie_data$type == "n_spl_events", 1 - cumsum(percent) + percent/2, 1 +(1 - + cumsum(percent) + percent/2 )))

ggplot(pie_data_20, aes( y=percent, x = (type)))+
  geom_bar(aes(fill = event), stat = "identity", width = 1, alpha = c(1,1,1,1,1,0.5,0.5,0.5,0.5,0.5), color = "lightgray", ) +
  geom_text_repel(aes(y = label, label = percent(percent, accuracy = 0.1)), fontface = "bold",
                  size = 5,
                  nudge_x = .25, 
                  show.legend = FALSE,
                  segment.alpha = 0) +#
  coord_polar("y") + 
  scale_fill_manual(values = c("#66A61E","#1B9E77", "#D95F02",  "#7570B3", "#E6AB02"), guide = FALSE) +
  theme_void()+
  theme(axis.text.x=element_blank())

ggsave(filename = "/cerberus/projects/racste/B_anynana/diff_expr/08_rMATS_S/figures/read_threshold_20_pie_AB_vs_TH.pdf")



pie_data_33  <- splice_summary_filt33 %>% select(event, n_spl_events, n_sig_spl_events) %>% 
  gather(-c(event), key = type, value = n_events) %>% 
  mutate(event = factor(event, levels = c("IR", "A3SS", "A5SS", "SE", "MXE"))) %>% 
  arrange(type, event) %>% 
  mutate(total = if_else(pie_data$type == "n_spl_events", 1841, 15553),
         percent = n_events/total,
        label = if_else(pie_data$type == "n_spl_events", 1 - cumsum(percent) + percent/2, 1 +(1 - + cumsum(percent) + percent/2 )))


ggplot(pie_data_33, aes( y=percent, x = type))+
  geom_bar(aes(fill = event), stat = "identity", width = 1, alpha = c(1,1,1,1,1,0.5,0.5,0.5,0.5,0.5), color = "lightgray", ) +
  geom_text_repel(aes(y = label, label = percent(percent, accuracy = 0.1)), fontface = "bold",
                  size = 5,
                  nudge_x = .25, 
                  show.legend = FALSE,
                  segment.alpha = 0) +#
  coord_polar("y") + 
  scale_fill_manual(values = c("#66A61E","#1B9E77", "#D95F02",  "#7570B3", "#E6AB02"), guide = FALSE) +
  theme_void()+
  theme(axis.text.x=element_blank()) 
ggsave(filename = "/cerberus/projects/racste/B_anynana/diff_expr/08_rMATS_S/figures/sample_threshold_33_pie_AB_vs_TH.pdf")

