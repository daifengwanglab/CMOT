# library
library(ggplot2)
library(plotly)
library(dplyr)
library(tidyverse)
library(hrbrthemes)
library(viridis)

df <-read.csv("./CMOT_Poster_Plots/LiuCancer/ATAC_TO_RNA_SilhouetteScore_Boxplot.csv")
colnames(df)[1] <- "MOFA+"
colnames(df)[2] <- "Measured \nGene Expression"
g3 <- melt(df)
# Boxplot basic

p<-  ggplot( g3,aes(x=variable, y=value)) +
  geom_boxplot(fill = "#b2b2b2", colour = "black") +
  #scale_fill_viridis(discrete = TRUE, alpha=0.6, option="E") +
  #theme_ipsum() +
  theme(axis.text=element_text(size=9,color="black", font),axis.title=element_text(size=12),
          #panel.border = element_rect(colour = "black", fill=NA, size=0.5),
          panel.background = element_blank(),
          axis.line = element_line(),
        legend.title = element_blank(), legend.position = "top") +
  #ggtitle("Basic boxplot") +
  xlab("") + ylab("Silhouette Score")

ggsave("./CMOT_Poster_Plots/LiuCancer/LiuCancer_SilhouetteScore_60_40_split_2.png", width=4, height=3)  

