#compare the correlation of state 1 vs pre surgery and post surgery vs pre surgery

#load libraries 

library(tidyverse)
library(ggpubr)

#load data 

#state 1 vs pre surgery
state1vsCCT0 <- read_csv("Results_csv/IMMERSE paper/SCRIPT 4/state1vsCCT0.csv")%>% 
  dplyr::rename("log2FC_State1_vs_Presurgery" = "log2FoldChange", "adj.P.Valstate1_vs_pre" = "padj")%>% 
  select(gene_ID,hgnc_symbol,log2FC_State1_vs_Presurgery,adj.P.Valstate1_vs_pre) 
# add a column of NAs
state1vsCCT0$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
state1vsCCT0$diffexpressed[state1vsCCT0$log2FC_State1_vs_Presurgery > 1 & state1vsCCT0$adj.P.Valstate1_vs_pre < 0.01] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
state1vsCCT0$diffexpressed[state1vsCCT0$log2FC_State1_vs_Presurgery < -1 & state1vsCCT0$adj.P.Valstate1_vs_pre < 0.01] <- "DOWN"
state1vsCCT0 <- state1vsCCT0 %>% dplyr::rename("diffexpressed_State1vsPresurgery" = "diffexpressed")
table(state1vsCCT0$diffexpressed_sepsisvsPresurgery)
colnames(state1vsCCT0)


#read in post surgery vs pre surgery
PostsurgeryvsPresurgery <- read_csv("Results_csv/IMMERSE paper/SCRIPT 1/Cohort_time_Cardiac_T1_vs_Cardiac_T0.csv")%>% 
  dplyr::rename("log2FC_postsurgeryvsPresurgery" = "log2FoldChange", "adj.P.Valpostsurgeryvspresurgery" = "padj") %>% 
  select(gene_ID,hgnc_symbol,log2FC_postsurgeryvsPresurgery,adj.P.Valpostsurgeryvspresurgery) 
# add a column of NAs
PostsurgeryvsPresurgery$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
PostsurgeryvsPresurgery$diffexpressed[PostsurgeryvsPresurgery$log2FC_postsurgeryvsPresurgery > 1 & PostsurgeryvsPresurgery$adj.P.Valpostsurgeryvspresurgery < 0.01] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
PostsurgeryvsPresurgery$diffexpressed[PostsurgeryvsPresurgery$log2FC_postsurgeryvsPresurgery < -1 & PostsurgeryvsPresurgery$adj.P.Valpostsurgeryvspresurgery < 0.01] <- "DOWN"
PostsurgeryvsPresurgery <- PostsurgeryvsPresurgery %>% dplyr::rename("diffexpressed_PostsurgeryvsPresurgery" = "diffexpressed")




compare1 <- merge(state1vsCCT0,PostsurgeryvsPresurgery)
compare1 <- compare1 %>% mutate(hgnc_symbol = coalesce(hgnc_symbol,gene_ID)) 
colnames(compare1)

compare1 <- compare1 %>% filter(between(log2FC_postsurgeryvsPresurgery, -15, -1)|
                                  between(log2FC_postsurgeryvsPresurgery, 1, 15)|
                                  between(log2FC_State1_vs_Presurgery, -15, -1)|
                                  between(log2FC_State1_vs_Presurgery, 1, 15))
compare1$Significance <- "NO"



# if log2Foldchange >  and pvalue < 0.05, set as "UP" 


compare1$Significance[compare1$diffexpressed_State1vsPresurgery == "UP" & compare1$diffexpressed_PostsurgeryvsPresurgery == "UP"] <- "State 1 and sterile"
compare1$Significance[compare1$diffexpressed_State1vsPresurgery == "DOWN" & compare1$diffexpressed_PostsurgeryvsPresurgery == "DOWN"] <- "State 1 and sterile"
compare1$Significance[compare1$diffexpressed_State1vsPresurgery == "UP" & compare1$diffexpressed_PostsurgeryvsPresurgery == "DOWN"] <- "State 1 and sterile"
compare1$Significance[compare1$diffexpressed_State1vsPresurgery == "DOWN" & compare1$diffexpressed_PostsurgeryvsPresurgery == "UP"] <- "State 1 and sterile"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
compare1$Significance[compare1$diffexpressed_State1vsPresurgery == "UP" & compare1$diffexpressed_PostsurgeryvsPresurgery == "NO"] <- "State 1 only"
compare1$Significance[compare1$diffexpressed_State1vsPresurgery == "DOWN" & compare1$diffexpressed_PostsurgeryvsPresurgery == "NO"] <- "State 1 only"
compare1$Significance[compare1$diffexpressed_State1vsPresurgery == "NO" & compare1$diffexpressed_PostsurgeryvsPresurgery == "UP"] <- "Sterile only"
compare1$Significance[compare1$diffexpressed_State1vsPresurgery == "NO" & compare1$diffexpressed_PostsurgeryvsPresurgery == "DOWN"] <- "Sterile only"
compare1$Significance[compare1$diffexpressed_State1vsPresurgery == "NO" & compare1$diffexpressed_PostsurgeryvsPresurgery == "NO"] <- "No change"



# Create a new column "delabel" to de, that will contain the name of genes differentially expressed (NA in case they are not)
compare1$delabel <- NA
compare1$delabel[compare1$Significance != "No change"] <- compare1$hgnc_symbol[compare1$Significance != "No change"]

colnames(compare1)
table(compare1$Significance)
compare1 <- compare1 %>% filter(!Significance == "No change",
                              !Significance == "NO")

 state1vssterile <- ggplot(compare1, aes(x = log2FC_State1_vs_Presurgery, y = log2FC_postsurgeryvsPresurgery))+
 # annotate(geom = "rect", xmin = -Inf, xmax = 0,   ymin = -Inf, ymax = 0,
  #         fill = "#e23c52", alpha = 0.1)+
#  annotate(geom = "rect", xmin = +Inf, xmax = 0,   ymin = +Inf, ymax = 0,
 #          fill = "#e23c52", alpha = 0.1)+
#  annotate(geom = "rect", xmin = 0, xmax = -Inf,   ymin = 0, ymax = +Inf,
 #          fill = "#fffdd0", alpha = 0.3)+
#  annotate(geom = "rect", xmin = 0, xmax = +Inf,   ymin = 0, ymax = -Inf,
 #          fill = "#fffdd0", alpha = 0.3)+
  geom_point(pch = 21, aes(fill = Significance), size = 2)+
  scale_fill_manual(values = c("State 1 and sterile" = "#CF3721", 
                               "State 1 only" = "#A1D6E2",
                               "Sterile only" = "#F5BE41",
                               "No change" = "#BCBABE"))+
  geom_hline(yintercept=0, linetype="dashed", color = "black")+
  geom_vline(xintercept=0, linetype="dashed", color = "black")+
  ylim(-10,10)+
  xlim(-10,10)+
  theme_classic()+
  geom_smooth(method='lm', formula= y~x, colour = "black")+ 
  stat_cor(method = "spearman", size =8)+
  guides(fill=guide_legend(override.aes=list(size = 4)))+
  theme(legend.position = "bottom")+
  theme(legend.text=element_text(size=18),
        legend.title = element_text(size=18))+
  theme(axis.text.x = element_text(color = "black", size = 18,  face = "plain"),
        axis.text.y = element_text(color = "black", size = 18, face = "plain"),  
        axis.title.x = element_text(color = "black", size = 18,  face = "plain"),
        axis.title.y = element_text(color = "black", size = 18,  face = "plain"))+
   ylab("Post- vs Pre-surgery Log2fold change")+
   xlab("State 1 vs Pre-surgery Log2fold change")


ggsave(plot = state1vssterile, "FIGURES/IMMERSE paper/SCRIPT 5/state1vssterile.png", dpi = 300, height =8, width = 8)

# get the  number of genes going up and down in each comparison
compare1$diffexpressed_State1vsPresurgery
counts <- compare1 %>% group_by(diffexpressed_PostsurgeryvsPresurgery,diffexpressed_State1vsPresurgery, Significance) %>% 
  summarise(count = n())
#save this for making figure later
write_csv(counts, "Results_csv/IMMERSE paper/SCRIPT 5/State1 vs sterile.csv")


#######################################################################################################################################
#######################################################################################################################################
# make the plot for the State 2 vs sterile 


#state 1 vs pre surgery
state2vsCCT0 <- read_csv("Results_csv/IMMERSE paper/SCRIPT 4/state2vsCCT0.csv")%>% 
  dplyr::rename("log2FC_State2_vs_Presurgery" = "log2FoldChange", "adj.P.Valstate2_vs_pre" = "padj")%>% 
  select(gene_ID,hgnc_symbol,log2FC_State2_vs_Presurgery,adj.P.Valstate2_vs_pre) 
# add a column of NAs
state2vsCCT0$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
state2vsCCT0$diffexpressed[state2vsCCT0$log2FC_State2_vs_Presurgery > 1 & state2vsCCT0$adj.P.Valstate2_vs_pre < 0.01] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
state2vsCCT0$diffexpressed[state2vsCCT0$log2FC_State2_vs_Presurgery < -1 & state2vsCCT0$adj.P.Valstate2_vs_pre < 0.01] <- "DOWN"
state2vsCCT0 <- state2vsCCT0 %>% dplyr::rename("diffexpressed_State2vsPresurgery" = "diffexpressed")
table(state2vsCCT0$diffexpressed_sepsisvsPresurgery)
colnames(state2vsCCT0)




compare2 <- merge(state2vsCCT0,PostsurgeryvsPresurgery)
compare2 <- compare2 %>% mutate(hgnc_symbol = coalesce(hgnc_symbol,gene_ID)) 
colnames(compare2)

compare2 <- compare2 %>% filter(between(log2FC_postsurgeryvsPresurgery, -15, -1)|
                                  between(log2FC_postsurgeryvsPresurgery, 1, 15)|
                                  between(log2FC_State2_vs_Presurgery, -15, -1)|
                                  between(log2FC_State2_vs_Presurgery, 1, 15))
compare2$Significance <- "NO"



# if log2Foldchange >  and pvalue < 0.05, set as "UP" 


compare2$Significance[compare2$diffexpressed_State2vsPresurgery == "UP" & compare2$diffexpressed_PostsurgeryvsPresurgery == "UP"] <- "State 2 and sterile"
compare2$Significance[compare2$diffexpressed_State2vsPresurgery == "DOWN" & compare2$diffexpressed_PostsurgeryvsPresurgery == "DOWN"] <- "State 2 and sterile"
compare2$Significance[compare2$diffexpressed_State2vsPresurgery == "UP" & compare2$diffexpressed_PostsurgeryvsPresurgery == "DOWN"] <- "State 2 and sterile"
compare2$Significance[compare2$diffexpressed_State2vsPresurgery == "DOWN" & compare2$diffexpressed_PostsurgeryvsPresurgery == "UP"] <- "State 2 and sterile"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
compare2$Significance[compare2$diffexpressed_State2vsPresurgery == "UP" & compare2$diffexpressed_PostsurgeryvsPresurgery == "NO"] <- "State 2 only"
compare2$Significance[compare2$diffexpressed_State2vsPresurgery == "DOWN" & compare2$diffexpressed_PostsurgeryvsPresurgery == "NO"] <- "State 2 only"
compare2$Significance[compare2$diffexpressed_State2vsPresurgery == "NO" & compare2$diffexpressed_PostsurgeryvsPresurgery == "UP"] <- "Sterile only"
compare2$Significance[compare2$diffexpressed_State2vsPresurgery == "NO" & compare2$diffexpressed_PostsurgeryvsPresurgery == "DOWN"] <- "Sterile only"
compare2$Significance[compare2$diffexpressed_State2vsPresurgery == "NO" & compare2$diffexpressed_PostsurgeryvsPresurgery == "NO"] <- "No change"


colnames(compare2)
table(compare2$Significance)
compare2 <- compare2 %>% filter(!Significance == "No change",
                                !Significance == "NO")

state2vssterile <- ggplot(compare2, aes(x = log2FC_State2_vs_Presurgery, y = log2FC_postsurgeryvsPresurgery))+
  # annotate(geom = "rect", xmin = -Inf, xmax = 0,   ymin = -Inf, ymax = 0,
  #         fill = "#e23c52", alpha = 0.1)+
  #  annotate(geom = "rect", xmin = +Inf, xmax = 0,   ymin = +Inf, ymax = 0,
  #          fill = "#e23c52", alpha = 0.1)+
  #  annotate(geom = "rect", xmin = 0, xmax = -Inf,   ymin = 0, ymax = +Inf,
  #          fill = "#fffdd0", alpha = 0.3)+
  #  annotate(geom = "rect", xmin = 0, xmax = +Inf,   ymin = 0, ymax = -Inf,
  #          fill = "#fffdd0", alpha = 0.3)+
  geom_point(pch = 21, aes(fill = Significance), size = 2)+
  scale_fill_manual(values = c("State 2 and sterile" = "#CF3721", 
                               "State 2 only" = "#A1D6E2",
                               "Sterile only" = "#F5BE41",
                               "No change" = "#BCBABE"))+
  geom_hline(yintercept=0, linetype="dashed", color = "black")+
  geom_vline(xintercept=0, linetype="dashed", color = "black")+
  ylim(-10,10)+
  xlim(-10,10)+
  theme_classic()+
  geom_smooth(method='lm', formula= y~x, colour = "black")+ 
  stat_cor(method = "spearman", size =8)+
  guides(fill=guide_legend(override.aes=list(size = 4)))+
  theme(legend.position = "bottom")+
  theme(legend.text=element_text(size=18),
        legend.title = element_text(size=18))+
  theme(axis.text.x = element_text(color = "black", size = 18,  face = "plain"),
        axis.text.y = element_text(color = "black", size = 18, face = "plain"),  
        axis.title.x = element_text(color = "black", size = 18,  face = "plain"),
        axis.title.y = element_text(color = "black", size = 18,  face = "plain"))+
  ylab("Post- vs Pre-surgery Log2fold change")+
  xlab("State 2 vs Pre-surgery Log2fold change")


ggsave(plot = state2vssterile, "FIGURES/IMMERSE paper/SCRIPT 5/state2vssterile.png", dpi = 300, height =8, width = 8)

# get the  number of genes going up and down in each comparison
compare2$diffexpressed_State1vsPresurgery
counts <- compare2 %>% group_by(diffexpressed_PostsurgeryvsPresurgery,diffexpressed_State2vsPresurgery, Significance) %>% 
  summarise(count = n())
#save this for making figure later
write_csv(counts, "Results_csv/IMMERSE paper/SCRIPT 5/State2 vs sterile.csv")
#######################################################################################################################################
#######################################################################################################################################
# make the plot for the State 3 vs sterile 


#state 3 vs pre surgery
state3vsCCT0 <- read_csv("Results_csv/IMMERSE paper/SCRIPT 4/state3vsCCT0.csv")%>% 
  dplyr::rename("log2FC_state3_vs_Presurgery" = "log2FoldChange", "adj.P.Valstate3_vs_pre" = "padj")%>% 
  select(gene_ID,hgnc_symbol,log2FC_state3_vs_Presurgery,adj.P.Valstate3_vs_pre) 
# add a column of NAs
state3vsCCT0$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
state3vsCCT0$diffexpressed[state3vsCCT0$log2FC_state3_vs_Presurgery > 1 & state3vsCCT0$adj.P.Valstate3_vs_pre < 0.01] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
state3vsCCT0$diffexpressed[state3vsCCT0$log2FC_state3_vs_Presurgery < -1 & state3vsCCT0$adj.P.Valstate3_vs_pre < 0.01] <- "DOWN"
state3vsCCT0 <- state3vsCCT0 %>% dplyr::rename("diffexpressed_state3vsPresurgery" = "diffexpressed")
table(state3vsCCT0$diffexpressed_sepsisvsPresurgery)
colnames(state3vsCCT0)




compare3 <- merge(state3vsCCT0,PostsurgeryvsPresurgery)
compare3 <- compare3 %>% mutate(hgnc_symbol = coalesce(hgnc_symbol,gene_ID)) 
colnames(compare3)

compare3 <- compare3 %>% filter(between(log2FC_postsurgeryvsPresurgery, -15, -1)|
                                  between(log2FC_postsurgeryvsPresurgery, 1, 15)|
                                  between(log2FC_state3_vs_Presurgery, -15, -1)|
                                  between(log2FC_state3_vs_Presurgery, 1, 15))
compare3$Significance <- "NO"



# if log2Foldchange >  and pvalue < 0.05, set as "UP" 


compare3$Significance[compare3$diffexpressed_state3vsPresurgery == "UP" & compare3$diffexpressed_PostsurgeryvsPresurgery == "UP"] <- "State 3 and sterile"
compare3$Significance[compare3$diffexpressed_state3vsPresurgery == "DOWN" & compare3$diffexpressed_PostsurgeryvsPresurgery == "DOWN"] <- "State 3 and sterile"
compare3$Significance[compare3$diffexpressed_state3vsPresurgery == "UP" & compare3$diffexpressed_PostsurgeryvsPresurgery == "DOWN"] <- "State 3 and sterile"
compare3$Significance[compare3$diffexpressed_state3vsPresurgery == "DOWN" & compare3$diffexpressed_PostsurgeryvsPresurgery == "UP"] <- "State 3 and sterile"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
compare3$Significance[compare3$diffexpressed_state3vsPresurgery == "UP" & compare3$diffexpressed_PostsurgeryvsPresurgery == "NO"] <- "State 3 only"
compare3$Significance[compare3$diffexpressed_state3vsPresurgery == "DOWN" & compare3$diffexpressed_PostsurgeryvsPresurgery == "NO"] <- "State 3 only"
compare3$Significance[compare3$diffexpressed_state3vsPresurgery == "NO" & compare3$diffexpressed_PostsurgeryvsPresurgery == "UP"] <- "Sterile only"
compare3$Significance[compare3$diffexpressed_state3vsPresurgery == "NO" & compare3$diffexpressed_PostsurgeryvsPresurgery == "DOWN"] <- "Sterile only"
compare3$Significance[compare3$diffexpressed_state3vsPresurgery == "NO" & compare3$diffexpressed_PostsurgeryvsPresurgery == "NO"] <- "No change"


colnames(compare3)
table(compare3$Significance)
compare3 <- compare3 %>% filter(!Significance == "No change",
                                !Significance == "NO")

state3vssterile <- ggplot(compare3, aes(x = log2FC_state3_vs_Presurgery, y = log2FC_postsurgeryvsPresurgery))+
  # annotate(geom = "rect", xmin = -Inf, xmax = 0,   ymin = -Inf, ymax = 0,
  #         fill = "#e23c52", alpha = 0.1)+
  #  annotate(geom = "rect", xmin = +Inf, xmax = 0,   ymin = +Inf, ymax = 0,
  #          fill = "#e23c52", alpha = 0.1)+
  #  annotate(geom = "rect", xmin = 0, xmax = -Inf,   ymin = 0, ymax = +Inf,
  #          fill = "#fffdd0", alpha = 0.3)+
  #  annotate(geom = "rect", xmin = 0, xmax = +Inf,   ymin = 0, ymax = -Inf,
  #          fill = "#fffdd0", alpha = 0.3)+
  geom_point(pch = 21, aes(fill = Significance), size = 2)+
  scale_fill_manual(values = c("State 3 and sterile" = "#CF3721", 
                               "State 3 only" = "#A1D6E2",
                               "Sterile only" = "#F5BE41",
                               "No change" = "#BCBABE"))+
  geom_hline(yintercept=0, linetype="dashed", color = "black")+
  geom_vline(xintercept=0, linetype="dashed", color = "black")+
  ylim(-10,10)+
  xlim(-10,10)+
  theme_classic()+
  geom_smooth(method='lm', formula= y~x, colour = "black")+ 
  stat_cor(method = "spearman", size =8)+
  guides(fill=guide_legend(override.aes=list(size = 4)))+
  theme(legend.position = "bottom")+
  theme(legend.text=element_text(size=18),
        legend.title = element_text(size=18))+
  theme(axis.text.x = element_text(color = "black", size = 18,  face = "plain"),
        axis.text.y = element_text(color = "black", size = 18, face = "plain"),  
        axis.title.x = element_text(color = "black", size = 18,  face = "plain"),
        axis.title.y = element_text(color = "black", size = 18,  face = "plain"))+
  ylab("Post- vs Pre-surgery Log2fold change")+
  xlab("State 3 vs Pre-surgery Log2fold change")


ggsave(plot = state3vssterile, "FIGURES/IMMERSE paper/SCRIPT 5/state3vssterile.png", dpi = 300, height =8, width = 8)


# get the  number of genes going up and down in each comparison
compare1$diffexpressed_State1vsPresurgery
counts <- compare3 %>% group_by(diffexpressed_PostsurgeryvsPresurgery,diffexpressed_state3vsPresurgery, Significance) %>% 
  summarise(count = n())
#save this for making figure later
write_csv(counts, "Results_csv/IMMERSE paper/SCRIPT 5/State3 vs sterile.csv")


