# analysis proportions of other endotypes in the in the states 

library(tidyverse)


#import the metadata with states. 

sweeney_groups <- read_csv("Results_csv/IMMERSE paper/SCRIPT 1/sweeney_groups.csv")
metadata_states_Phate <- read_csv("Metadata/metadata_states_Phate.csv")
mars_groupings <- read_csv("Results_csv/IMMERSE paper/SCRIPT 1/mars_groupings.csv")

df <- merge(metadata_states_Phate, sweeney_groups)
df <- merge(df,mars_groupings )

colnames(df)
df$Cohort
#SRS davenport
df <- df %>% 
  mutate(SRS_davenport = case_when(SRS_davenport == "SRS3" ~ "SRS2",
                                   SRS_davenport == "SRS1" ~ "SRS1",
                                   SRS_davenport == "SRS2" ~ "SRS2")) 


colnames(df)
df$SRS_davenport <- as.factor(df$SRS_davenport)

df_summary <- df %>%  group_by(State,SRS_davenport) %>% summarise(Count=n())%>%
  mutate(freq = Count / sum(Count)*100) %>% 
  mutate_at(vars(freq), funs(round(., 1)))


SRS_grouping <- ggplot(df_summary, aes(x = State, y = Count, fill = SRS_davenport, label = freq))+
  geom_bar(stat="identity", colour = "black")+
  theme_classic()+
  scale_fill_manual(values = c("SRS1" = "#810103",
                               "SRS2" ="#4d7faf",
                               "SRS3" = "#010089"))+ 
  theme(text = element_text(size = 16))+
  geom_text(size = 6, position = position_stack(vjust = 0.5))

ggsave(plot = SRS_grouping, "FIGURES/IMMERSE paper/SCRIPT 6/SRS_dp_state_barplot.png", height = 4, width = 5, dpi = 300)



#SRS extended
df_summary <- df %>% group_by(State,SRS_extended) %>% summarise(Count=n())%>%
  mutate(freq = Count / sum(Count)*100) %>% 
  mutate_at(vars(freq), funs(round(., 1)))


SRS_extended_grouping <- ggplot(df_summary, aes(x = State, y = Count, fill = SRS_extended, label = freq))+
  geom_bar(stat="identity", colour = "black")+
  theme_classic()+
  scale_fill_manual(values = c("SRS1" = "#810103",
                               "SRS2" ="#4d7faf",
                               "SRS3" = "#010089"))+ 
  theme(text = element_text(size = 16))+
  geom_text(size = 6, position = position_stack(vjust = 0.5))

ggsave(plot = SRS_extended_grouping, "FIGURES/IMMERSE paper/SCRIPT 6/SRS_extended_barplot.png", height = 4, width = 5, dpi = 300)




#Sweeney
df_summary <- df %>% group_by(State,`Sweeney 2019`) %>% summarise(Count=n())%>%
  mutate(freq = Count / sum(Count)*100) %>% 
  mutate_at(vars(freq), funs(round(., 1)))
df_summary$fre

Sweeney_grouping <- ggplot(df_summary, aes(x = State, y = Count, fill = `Sweeney 2019`, label = freq))+
  geom_bar(stat="identity", colour = "black")+
  theme_classic()+
  scale_fill_manual(values = c("adaptive" = "#6cc0e5",
                               "coagulopathic" ="#fbc93d",
                               "inflammopathic" = "#fb4f4f"))+ 
  theme(text = element_text(size = 16))+
  geom_text(size = 6, position = position_stack(vjust = 0.5))

ggsave(plot = Sweeney_grouping, "FIGURES/IMMERSE paper/SCRIPT 6/Sweeney_grouping.png", height = 4, width = 5, dpi = 300)


#mars
colnames(df)
df_summary <- df %>% group_by(State,`MARS`) %>% summarise(Count=n())%>%
  mutate(freq = Count / sum(Count)*100) %>% 
  mutate_at(vars(freq), funs(round(., 1)))

mars_grouping <- ggplot(df_summary, aes(x = State, y = Count, fill = `MARS`, label = freq))+
  geom_bar(stat="identity", colour = "black")+
  theme_classic()+
  scale_fill_manual(values = c("Mars1" = "#c7deb2",
                               "Mars2" ="#2e64aa",
                               "Mars3" = "#9dacd5",
                               "Mars4" = "#82ba55"))+ 
  theme(text = element_text(size = 16))+
  geom_text(size = 6, position = position_stack(vjust = 0.5))

ggsave(plot = mars_grouping, "FIGURES/IMMERSE paper/SCRIPT 6/mars_grouping.png", height = 4, width = 5, dpi = 300)


