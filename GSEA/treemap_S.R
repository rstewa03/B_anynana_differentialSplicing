#' ---
#' title: "REVIGO analyses for between-season DS and DE genesets"
#' output: html_document
#' author: "Rachel Steward"
#' date: '2021-05-07'
#' ---
#' 
# A treemap R script produced by the REVIGO server at http://revigo.irb.hr/
# If you found REVIGO useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800

# author: Anton Kratz <anton.kratz@gmail.com>, RIKEN Omics Science Center, Functional Genomics Technology Team, Japan
# created: Fri, Nov 02, 2012  7:25:52 PM
# last change: Fri, Nov 09, 2012  3:20:01 PM

# -----------------------------------------------------------------------------
# If you don't have the treemap package installed, uncomment the following line:
# install.packages( "treemap" );
library(treemap) 								# treemap package by Martijn Tennekes
library(tidyverse)
library(ggrepel)
library(ggnewscale)
library(cowplot)
library(eulerr)
# --------------------------------------------------------------------------
# Data from REVIGO. Scroll down for plot configuration options.

# Prep data from Revigo for analysis by removing 4 lines above header
# tail -n +5 DEDS_ab_BP_parchild_REVIGO_treemap.csv > file.txt.new && mv file.txt.new DEDS_ab_BP_parchild_REVIGO_treemap.csv


#set working directory
setwd("~/OneDrive/B_anynana_altSplicing/PanicTime/genelists/REVIGO_S/treemap/")
#setwd("~/OneDrive/B_anynana_altSplicing/B_anynana_GSEA/edgeRDEvsDS/Treemap/")
revigo.names <- c("term_ID","description","freqInDbPercent","abslog10pvalue","uniqueness","dispensability","representative");

temp_revigo <- list.files(pattern="*.csv")
my_names <- str_sub(temp_revigo, end = -5)
my_revigo <- lapply(temp_revigo, read.csv)
names(my_revigo) <- my_names
names(my_revigo[1])
for (i in 1:length(my_revigo)) {
  revigo.data <- my_revigo[[i]]
  revigo.data$abslog10pvalue <- abs(revigo.data$log10pvalue)
  
  # by default, outputs to a PDF file
  pdf(file= paste(names(my_revigo[i]),"revigo_treemap_simple.pdf", sep = "_"), width=8, height=4.5 )# width=16, height=9 ) # width and height are in inches
  
  # check the tmPlot command documentation for all possible parameters - there are a lot more
  treemap(
    revigo.data,
    index = c("representative"), #,"description"),
    vSize = "abslog10pvalue",
    type = "categorical",
    vColor = "representative",
    title =  names(my_revigo[i]),
    inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
    lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
    bg.labels = "#CCCCCCAA",     # define background color of group labels
    # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
    position.legend = "none"
  )
  dev.off()
  
  pdf(file= paste(names(my_revigo[i]),"revigo_treemap.pdf", sep = "_"), width=8, height=4.5 )# width=16, height=9 ) # width and height are in inches
  
  # check the tmPlot command documentation for all possible parameters - there are a lot more
  treemap(
    revigo.data,
    index = c("representative","description"),
    vSize = "abslog10pvalue",
    type = "categorical",
    vColor = "representative",
    title =  names(my_revigo[i]),
    inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
    lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
    bg.labels = "#CCCCCCAA",     # define background color of group labels
    # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
    position.legend = "none"
  )
  dev.off()
}

#### Semantic space ####
# setwd("~/OneDrive/B_anynana_altSplicing/B_anynana_GSEA/DEvsAllDS/parentChild/REVIGO/SemSpace/")
setwd("/Users/rachelsteward/OneDrive/B_anynana_altSplicing/PanicTime/genelists/REVIGO_S/SemSpace/")
semspace.names <- c("term_ID","description","frequency","plot_X","plot_Y","plot_size","log10 p-value","uniqueness","dispensability","representative", "eliminated")

temp_semspace <- list.files(pattern="*.csv")
my_names <- str_sub(temp_semspace, end = -5)
my_semspace <- lapply(temp_semspace, read.csv)
names(my_semspace) <- my_names
names(my_semspace[1])
plots <- vector('list', length(my_semspace))

for (i in 1:length(my_semspace)) {
  my_semspace[[i]] <- my_semspace[[i]] [(my_semspace[[i]]$plot_X != "null" & my_semspace[[i]]$plot_Y != "null"), ];
  my_semspace[[i]]$plot_X <- as.numeric( as.character(my_semspace[[i]]$plot_X) );
  my_semspace[[i]]$plot_Y <- as.numeric( as.character(my_semspace[[i]]$plot_Y) );
  my_semspace[[i]]$plot_size <- as.numeric( as.character(my_semspace[[i]]$plot_size) );
  my_semspace[[i]]$log10.p.value <- as.numeric( as.character(my_semspace[[i]]$log10.p.value) );
  my_semspace[[i]]$frequency <- as.numeric(as.numeric(sub("%", "",my_semspace[[i]]$frequency,fixed=TRUE))/100); 
  my_semspace[[i]]$uniqueness <- as.numeric( as.character(my_semspace[[i]]$uniqueness) );
  my_semspace[[i]]$dispensability <- as.numeric( as.character(my_semspace[[i]]$dispensability) );
  my_semspace[[i]]$name <- names(my_semspace[i])
  
  ex <- my_semspace[[i]] [ my_semspace[[i]]$dispensability < 0.15 & my_semspace[[i]]$log10.p.value < log10(0.05), ]
  one.x_range = max(my_semspace[[i]]$plot_X) - min(my_semspace[[i]]$plot_X);
  one.y_range = max(my_semspace[[i]]$plot_Y) - min(my_semspace[[i]]$plot_Y);
  
  p1 <- ggplot( data = my_semspace[[i]][ my_semspace[[i]]$log10.p.value < log10(0.05), ] ) +
    geom_point( aes( plot_X, plot_Y, colour = log10.p.value, size = plot_size), alpha = I(0.6) ) + 
    scale_colour_gradientn( colours = c("blue", "green", "yellow", "red") , limits = c( -6, 0))+
    geom_point( aes(plot_X, plot_Y, size = plot_size), shape = 21, fill = "transparent", colour = I (alpha ("black", 0.6) )) +
    scale_size( range=c(5, 30), guide = FALSE) + 
    theme_bw() +
    geom_text( data = ex, aes(plot_X, plot_Y, label = description), colour = I(alpha("black", 0.85)), size = 3 ) + 
    labs (y = "Semantic space x", x = "Semantic space y") + theme(legend.key = element_blank()) +
    xlim(min(my_semspace[[i]]$plot_X)-one.x_range/10,max(my_semspace[[i]]$plot_X)+one.x_range/10) +
    ylim(min(my_semspace[[i]]$plot_Y)-one.y_range/10,max(my_semspace[[i]]$plot_Y)+one.y_range/10);
  ggsave(plot = p1, file= paste0(names(my_semspace[i]),".pdf"), width=9, height=9 ) # width and height are in inches
  
  semspace.colors <- c("#0072B5","#0072B5",
                       "#7876B1","#7876B1",
                       "#BC3C29","#BC3C29",
                       "#0072B5","#0072B5",
                       "#7876B1","#7876B1",
                       "#BC3C29","#BC3C29")
  message(i)
  plots[[i]] <- local({
    p1 <- ggplot( ) +
      geom_point(data = my_semspace[[i]][ my_semspace[[i]]$log10.p.value < log10(0.05), ], aes( plot_X, plot_Y, size = plot_size,  alpha = -(log10.p.value)), colour = semspace.colors[i] ) + 
      geom_point( data = my_semspace[[i]][ my_semspace[[i]]$log10.p.value < log10(0.05), ], aes(plot_X, plot_Y, size = plot_size), shape = 21, fill = "transparent", colour = I (alpha ("black", 0.6) )) +
      scale_colour_gradient(limits = c( -6, 0))+
      scale_size( range=c(5, 20), guide = FALSE) + 
      theme_bw() +
      geom_text_repel( data = my_semspace[[i]][ my_semspace[[i]]$dispensability < 0.1 &  my_semspace[[i]]$log10.p.value < log10(0.05), ], 
                       aes(plot_X, plot_Y, label = description), colour = I(alpha("black", 0.85)), size = 4, ) + 
      labs (y = "Semantic space x", x = "Semantic space y") + theme(legend.key = element_blank()) +
      xlim(-9, 9) +
      ylim(-9,9) +
      theme(legend.position = "bottom")
    print(p1)
  })
}
names(my_semspace)
ab_BP <- plot_grid(plots[[1]]+ theme(legend.position = "none"),
                   plots[[3]]+ theme(legend.position = "none"), 
                   plots[[5]]+ theme(legend.position = "none"), 
                   NULL,
                   cowplot::get_legend(plots[[9]]),
                   NULL,
                   ncol = 3, nrow =2, labels = c("A", "B", "C", ""), 
                   rel_widths = c(1,1,1), rel_heights = c(1,.1))
ggsave(ab_BP, filename = "ab_BP_parchild_REVIGO_semspace.pdf", height = 7, width = 20)

th_BP <- plot_grid(plots[[7]]+ theme(legend.position = "none"), #DE 
                   plots[[9]]+ theme(legend.position = "none"), #DEDS
                   plots[[11]]+ theme(legend.position = "none"), #DS
                   NULL,
                   cowplot::get_legend(plots[[11]]),
                   NULL,
                   ncol = 3, nrow =2, labels = c("D", "E", "F", ""), 
                   rel_widths = c(1,1,1), rel_heights = c(1,.1))
ggsave(th_BP, filename = "th_BP_parchild_REVIGO_semspace.pdf", height = 7, width = 20)


ab_MF <- plot_grid(plots[[2]]+ theme(legend.position = "none"),
                   plots[[4]]+ theme(legend.position = "none"), 
                   plots[[6]]+ theme(legend.position = "none"), 
                   NULL,
                   cowplot::get_legend(plots[[10]]),
                   NULL,
                   ncol = 3, nrow =2, labels = c("A", "B", "C", ""), 
                   rel_widths = c(1,1,1), rel_heights = c(1,.1))
ggsave(ab_MF, filename = "ab_MF_parchild_REVIGO_semspace.pdf", height = 7, width = 20)

th_MF <- plot_grid(plots[[8]]+ theme(legend.position = "none"),
                   plots[[10]]+ theme(legend.position = "none"), 
                   plots[[12]]+ theme(legend.position = "none"), 
                   NULL,
                   cowplot::get_legend(plots[[12]]),
                   NULL,
                   ncol = 3, nrow =2, labels = c("D", "E", "F", ""), 
                   rel_widths = c(1,1,1), rel_heights = c(1,.1))
ggsave(th_MF, filename = "th_MF_parchild_REVIGO_semspace.pdf", height = 7, width = 20)

my_names <- names(my_semspace)
for (i in 1:length(my_semspace)) {
  my_revigo[[i]] <- my_semspace[[i]] %>% select(term_ID, log10.p.value)
  colnames(my_revigo[[i]]) <- c("term_ID", my_names[i])
}
all_terms <- data.frame(my_revigo[[1]]) 
for (i in 2:length(temp_revigo)){
  all_terms<- full_join(all_terms, my_revigo[[i]])
}

all_terms <- all_terms %>% 
  mutate_at(vars(-term_ID), ~ if_else(.x<0.05,1,0)) %>%
  mutate_at(vars(-term_ID), ~ if_else(is.na(.x),0,.x))
semspace.colors2 <- c("#7876B1","#0072B5","#BC3C29")
# 
# pdf(file ="eulerr_overlapping_GO_terms.pdf", height = 4, width = 4) 
# plot(euler(all_terms[, c(3,1,5)],  shape= "ellipse"),
#      quantities = list(fontsize =16),
#      fill = semspace.colors2, alpha = 0.85,
#      labels = c( "DE-DS", "DE only","DS only") ,
#      lty = 1, title = "AB_BP")
# # plot(venn(all_terms[, c(2,6,10)], labels = c( "DE-DS", "DE only","DS only"),fill = semspace.colors2, alpha = 0.85))
# 
# 
# plot(euler(all_terms[, c(4,2,6)],  shape= "ellipse"), 
#      quantities = list(fontsize =16),
#      fill = semspace.colors2, alpha = 0.85,
#      labels = c( "DE-DS", "DE only","DS only") ,
#      lty = 1, title = "AB_MF")  
# # plot(venn(all_terms[, c(3,7,11)]),labels = c( "DE-DS", "DE only","DS only") , fill = semspace.colors2, alpha = 0.85)
# 
# plot(euler(all_terms[, c(9,7,11)],  shape= "ellipse"), 
#      quantities = list(fontsize =16),
#      fill = semspace.colors2, alpha = 0.85,
#      labels = c( "DE-DS", "DE only","DS only") ,
#      lty = 1, title = "TH_BP")  
# 
# # plot(venn(all_terms[, c(4,8,12)]),labels = c( "DE-DS", "DE only","DS only") , fill = semspace.colors2, alpha = 0.85)
# 
# plot(euler(all_terms[, c(10,8,12)],  shape= "ellipse"), 
#      quantities = list(fontsize =16),
#      fill = semspace.colors2, alpha = 0.85,
#      labels = c( "DE-DS", "DE only","DS only") ,
#      lty = 1)  
# # plot(venn(all_terms[, c(5,9,13)]),labels = c( "DE-DS", "DE only","DS only") , fill = semspace.colors2, alpha = 0.85)
# 
# dev.off()

ab_BP_overlap <- all_terms[, c(1,2,4,6)]%>% filter_at(vars(-term_ID), any_vars(. > 0)) %>% mutate(sum = rowSums(.[2:4]), group = "ab_BP") %>% filter(sum >1)
colnames(ab_BP_overlap) <- c("termID", "DE", "DEDS", "DS", "sum", "group")
ab_MF_overlap <- all_terms[, c(1,3,5,7)] %>% filter_at(vars(-term_ID), any_vars(. > 0))%>% mutate(sum = rowSums(.[2:4]), group = "ab_MF") %>% filter(sum >1)
colnames(ab_MF_overlap) <- c("termID", "DE", "DEDS", "DS", "sum", "group")
th_BP_overlap <- all_terms[, c(1,8,10,12)]%>% filter_at(vars(-term_ID), any_vars(. > 0))%>% mutate(sum = rowSums(.[2:4]), group = "th_BP") %>% filter(sum >1) 
colnames(th_BP_overlap) <- c("termID", "DE", "DEDS", "DS", "sum", "group")
th_MF_overlap <- all_terms[, c(1,9,11,13)]%>% filter_at(vars(-term_ID), any_vars(. > 0))%>% mutate(sum = rowSums(.[2:4]), group = "th_MF") %>% filter(sum >1) 
colnames(th_MF_overlap) <- c("termID", "DE", "DEDS", "DS", "sum", "group")

overlap <- union(ab_BP_overlap, ab_MF_overlap) %>% union(th_BP_overlap)%>% union(th_MF_overlap)
write.csv(overlap, "/Users/rachelsteward/OneDrive/B_anynana_altSplicing/PanicTime/genelists/REVIGO_S/overlappping_go_terms.csv")

setwd("/Users/rachelsteward/OneDrive/B_anynana_altSplicing/PanicTime/genelists/REVIGO_S/SemSpace/all/")
temp_combined <- list.files(pattern="*csv")
my_names <- str_sub(temp_combined, end = -12)
my_combined <- lapply(temp_combined, read.csv)
names(my_combined) <- my_names
names(my_combined)
for (i in 1:length(my_combined)) {
  my_combined[[i]] <- my_combined[[i]][,c(1,4,5)] %>% 
  mutate(plot_X = as.numeric(as.character(PlotX)),
         plot_Y = as.numeric(PlotY)) %>%
    select(TermID, plot_X, plot_Y) %>% rename(term_ID = TermID)
}
names(my_combined[[i]])
my_semspace_ab_BP <- my_semspace[grepl("ab_DE_S_BP|ab_DEDS_S_BP|ab_DS_S_BP", names(my_semspace))]
my_semspace_th_BP <- my_semspace[grepl("th_DE_S_BP|th_DEDS_S_BP|th_DS_S_BP", names(my_semspace))] 
my_semspace_ab_MF <- my_semspace[grepl("ab_DE_S_MF|ab_DEDS_S_MF|ab_DS_S_MF", names(my_semspace))]
my_semspace_th_MF <- my_semspace[grepl("th_DE_S_MF|th_DEDS_S_MF|th_DS_S_MF", names(my_semspace))] 

combined_list <- list()
factorNum <- c("DE_only", "DE_DS", "DS_only")
for (i in 1:length(my_semspace_ab_BP)) {
  my_semspace_ab_BP[[i]] <- my_semspace_ab_BP[[i]] %>% dplyr::select(-c(plot_X, plot_Y)) %>% filter(dispensability < 0.3) %>% left_join(my_combined$ab_S_BP) %>% mutate(description = as.character(description))
  my_semspace_ab_BP[[i]]$name <- factorNum[[i]]
  combined_list$my_combined_ab_BP <-rbind(combined_list$my_combined_ab_BP, my_semspace_ab_BP[[i]])
}
combined_list$my_combined_ab_BP$textBreaks <- sapply(strsplit(combined_list$my_combined_ab_BP$description, " "), function(x) {
  spacePosition <- cumsum(nchar(x))
  placeBreak <- spacePosition[which(diff(spacePosition %/% 50) == 1)] + 1
  result <- paste(x, collapse = " ")
  for(i in placeBreak) {
    substring(result, i, i) <- "\n"
  }
  result
})

for (i in 1:length(my_semspace_th_BP)) {
  my_semspace_th_BP[[i]] <- my_semspace_th_BP[[i]] %>% dplyr::select(-c(plot_X, plot_Y)) %>% filter(dispensability < 0.3) %>% left_join(my_combined$th_S_BP)%>% mutate(description = as.character(description))
  my_semspace_th_BP[[i]]$name <- factorNum[[i]]
  combined_list$my_combined_th_BP <-rbind(combined_list$my_combined_th_BP, my_semspace_th_BP[[i]])
}
combined_list$my_combined_th_BP$textBreaks <- sapply(strsplit(combined_list$my_combined_th_BP$description, " "), function(x) {
  spacePosition <- cumsum(nchar(x))
  placeBreak <- spacePosition[which(diff(spacePosition %/% 50) == 1)] + 1
  result <- paste(x, collapse = " ")
  for(i in placeBreak) {
    substring(result, i, i) <- "\n"
  }
  result
})

for (i in 1:length(my_semspace_ab_MF)) {
  my_semspace_ab_MF[[i]] <- my_semspace_ab_MF[[i]] %>% dplyr::select(-c(plot_X, plot_Y)) %>%filter(dispensability < 0.3)  %>% left_join(my_combined$ab_S_MF) %>% mutate(description = as.character(description))
  my_semspace_ab_MF[[i]]$name <- factorNum[[i]]
  combined_list$my_combined_ab_MF <-rbind(combined_list$my_combined_ab_MF, my_semspace_ab_MF[[i]])
}
combined_list$my_combined_ab_MF$textBreaks <- sapply(strsplit(combined_list$my_combined_ab_MF$description, " "), function(x) {
  spacePosition <- cumsum(nchar(x))
  placeBreak <- spacePosition[which(diff(spacePosition %/% 50) == 1)] + 1
  result <- paste(x, collapse = " ")
  for(i in placeBreak) {
    substring(result, i, i) <- "\n"
  }
  result
})

for (i in 1:length(my_semspace_th_MF)) {
  my_semspace_th_MF[[i]] <- my_semspace_th_MF[[i]] %>% dplyr::select(-c(plot_X, plot_Y)) %>%filter(dispensability < 0.3)  %>% left_join(my_combined$th_S_MF)%>% mutate(description = as.character(description))
  my_semspace_th_MF[[i]]$name <- factorNum[[i]]
  combined_list$my_combined_th_MF <-rbind(combined_list$my_combined_th_MF, my_semspace_th_MF[[i]])
}
combined_list$my_combined_th_MF$textBreaks <- sapply(strsplit(combined_list$my_combined_th_MF$description, " "), function(x) {
  spacePosition <- cumsum(nchar(x))
  placeBreak <- spacePosition[which(diff(spacePosition %/% 50) == 1)] + 1
  result <- paste(x, collapse = " ")
  for(i in placeBreak) {
    substring(result, i, i) <- "\n"
  }
  result
})

semspace.colors2 <- c("#0072B5","#7876B1","#BC3C29")
semspace.colors3 <-  c("dodgerblue4","orchid4","darkred")


plots2 <- vector('list', length(combined_list))
plots3 <- vector('list', length(combined_list))

for (i in 1:length(combined_list)){
  message(i)
  plots2[[i]] <- local({
    p1 <- ggplot( ) +
      geom_point(data = combined_list[[i]], aes( plot_X, plot_Y, size = plot_size, colour = name ), alpha =0.75) + 
      geom_point( data = combined_list[[i]], aes(plot_X, plot_Y, size = plot_size), shape = 21, fill = "transparent", colour = I (alpha ("black", 0.6) )) +
      scale_color_manual(values = semspace.colors2,
                         labels = c("DE only","DS only"))+
      scale_size( range=c(5, 20), guide = FALSE) + 
      theme_minimal(base_size = 15) +
      new_scale_color() +
      #geom_text_repel( data = combined_list[[i]][combined_list[[i]]$dispensability < 0.001, ], 
      #aes(plot_X, plot_Y, label = description, colour = name), size = 4 ) +
      geom_text_repel(
        data = na.omit(combined_list[[i]][combined_list[[i]]$plot_X > 0 & combined_list[[i]]$dispensability < 0.1,]),
        aes(plot_X, plot_Y, label = description, colour = name), size = 3.5,
       nudge_x       = 10 - na.omit(combined_list[[i]][combined_list[[i]]$plot_X > 0 & combined_list[[i]]$dispensability < 0.1,]$plot_X),
        segment.size  = 0.2,
        direction     = "y",
        hjust         = 0)  +
      geom_text_repel(
        data = na.omit(combined_list[[i]][combined_list[[i]]$plot_X < 0 & combined_list[[i]]$dispensability < 0.1,]),
        aes(plot_X, plot_Y, label = description, colour = name, size = log10.p.value), size =3.5,
         nudge_x       = -10 - na.omit(combined_list[[i]][combined_list[[i]]$plot_X < 0 & combined_list[[i]]$dispensability < 0.1,]$plot_X),
        segment.size  = 0.2,
        direction     = "y",
        hjust         = 1) +
      scale_color_manual(values = semspace.colors3,
                         labels = c("DE only","DS only"))+
      labs (x = "Semantic space x", y = "Semantic space y") + theme(legend.key = element_blank()) +
      xlim(-20,20) +
      ylim(-9,9) +
      theme(legend.position = "bottom") # This widens the right margin
    print(p1)
  })
  plots3[[i]] <- local({
    p1 <- ggplot( ) +
      geom_point(data = combined_list[[i]][combined_list[[i]]$log10.p.value < log10(0.01), ], aes( plot_X, plot_Y, size = plot_size, colour = name ), alpha =0.75) +
      geom_point( data = combined_list[[i]][combined_list[[i]]$log10.p.value < log10(0.01), ], aes(plot_X, plot_Y, size = plot_size), shape = 21, fill = "transparent", colour = I (alpha ("black", 0.6) )) +
      scale_color_manual(values = semspace.colors2,
                         labels = c("DE only","DS only"))+
      scale_size( range=c(5, 20), guide = FALSE) + 
      theme_minimal(base_size = 15) +
      new_scale_color() +
      #geom_text_repel( data = combined_list[[i]][combined_list[[i]]$dispensability < 0.001, ], 
      #aes(plot_X, plot_Y, label = description, colour = name), size = 4 ) +
      geom_text_repel(
        data = na.omit(combined_list[[i]][combined_list[[i]]$plot_X > 0 & combined_list[[i]]$dispensability < 0.1 & combined_list[[i]]$log10.p.value < log10(0.01),]),
        aes(plot_X, plot_Y, label = description, colour = name), size = 3.5,
        nudge_x       = 10 - na.omit(combined_list[[i]][combined_list[[i]]$plot_X > 0 & combined_list[[i]]$dispensability < 0.1& combined_list[[i]]$log10.p.value < log10(0.01),]$plot_X) ,
        segment.size  = 0.2,
        direction     = "y",
        hjust         = 0)  +
      geom_text_repel(
        data = na.omit(combined_list[[i]][combined_list[[i]]$plot_X < 0 & combined_list[[i]]$dispensability < 0.1 & combined_list[[i]]$log10.p.value < log10(0.01),]),
        aes(plot_X, plot_Y, label = description, colour = name, size = log10.p.value), size =3.5,
        nudge_x       = -10 - na.omit(combined_list[[i]][combined_list[[i]]$plot_X < 0 & combined_list[[i]]$dispensability < 0.1 & combined_list[[i]]$log10.p.value < log10(0.01),]$plot_X),
        segment.size  = 0.2,
        direction     = "y",
        hjust         = 1) +
      scale_color_manual(values = semspace.colors3,
                         labels = c("DE only","DS only"))+
      labs (x = "Semantic space x", y = "Semantic space y") + theme(legend.key = element_blank()) +
      xlim(-15,15) +
      ylim(-9,9) +
      theme(legend.position = "bottom")
    print(p1)
  })
}


all_overlap <- plot_grid(plots2[[1]]+ theme(legend.position = "none"),
                         plots2[[2]]+ theme(legend.position = "none"), 
                         plots2[[3]]+ theme(legend.position = "none"),  
                         plots2[[4]]+ theme(legend.position = "none"), 
                         ncol = 2, nrow =2, labels = c("A", "B", "C", "D"), 
                         rel_widths = c(1,1), rel_heights = c(1,1))
ggsave(all_overlap, filename = "all_overlap_parchild_REVIGO_semspace.pdf", height = 12, width = 20)


all_overlap_0.01 <- plot_grid(plots3[[1]]+ theme(legend.position = "none"),
                              plots3[[2]]+ theme(legend.position = "none"), 
                              plots3[[3]]+ theme(legend.position = "none"),  
                              plots3[[4]]+ theme(legend.position = "none"), 
                              ncol = 2, nrow =2, labels = c("A", "B", "C", "D"), 
                              rel_widths = c(1,1), rel_heights = c(1,1))
ggsave(all_overlap_0.01, filename = "all_overlap_0.01_parchild_REVIGO_semspace.pdf", height = 12, width = 20)



plots4 <- vector('list', length(combined_list))
plots5 <- vector('list', length(combined_list))
for (i in 1:length(combined_list)){
  message(i)
  plots4[[i]] <- local({
    p1 <- ggplot( ) +
      geom_point(data = combined_list[[i]], aes( plot_X, plot_Y, size = plot_size, colour = name ), alpha =0.75) + 
      geom_point( data = combined_list[[i]], aes(plot_X, plot_Y, size = plot_size), shape = 21, fill = "transparent", colour = I (alpha ("black", 0.6) )) +
      scale_color_manual(values = semspace.colors2[c(2,1,3)],
                         labels = c("DE & DS","DE only", "DS only"))+
      scale_size( range=c(2, 10), guide = FALSE) + 
      labs(x = "Semantic Space X", y = "Semantic Space Y") +
      theme_bw(base_size = 15)
    print(p1)
  })

  
  plots5[[i]] <- local({
    p1 <- ggplot( ) +
      geom_point(data = combined_list[[i]][combined_list[[i]]$log10.p.value < log10(0.01), ], aes( plot_X, plot_Y, size = plot_size, colour = name ), alpha =0.75) +
      geom_point( data = combined_list[[i]][combined_list[[i]]$log10.p.value < log10(0.01), ], aes(plot_X, plot_Y, size = plot_size), shape = 21, fill = "transparent", colour = I (alpha ("black", 0.6) )) +
      scale_color_manual(values = semspace.colors2[c(2,1,3)],
                         labels = c("DE & DS","DE only", "DS only"))+
      scale_size( range=c(2, 10), guide = FALSE) + 
      labs(x = "Semantic Space X", y = "Semantic Space Y") +
      theme_bw(base_size = 15)
     
    print(p1)
  })
 
}

save(plots4, file = paste0("Fig2_parchild_REVIGO_semspacePlots_season.Rdata"))
save(plots5, file = paste0("Fig2_0.01_parchild_REVIGO_semspacePlots_season.Rdata"))

