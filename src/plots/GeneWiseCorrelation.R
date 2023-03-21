# code was adapted from: https://github.com/Noble-Lab/Polarbear
library(data.table)
library(readxl)
library(ROCR)
library(Matrix)
library(scran)
library(ggplot2)
library(ggallin)

df <- read.csv("./GeneWiseCorrelations.csv") #columns: methods, rows: inferred vs measured correlation of a gene 

#colnames(df)[2]<- "Polarbear-coassay" #for peak-wise
#colnames(df)[3]<- "Polarbear" #for peak-wise

#colnames(df)[2]<- "MOFA+" #for gene-wise
colnames(df)[4]<- "MOFA+" #for peak-wise

count_diagonal = table(df$`CMOT`>df$`MOFA+`)
pval = wilcox.test(df$`MOFA+`, df$CMOT, alternative = "less")$p.value

p2_cor_compare <- ggplot(df, aes(x = `MOFA+`, y = CMOT)) + 
  #geom_point(size = 1, alpha=.5, color="#2c7bb6")+ ggtitle("peak-wise AUROC")+
  geom_point(size = 1, alpha=.5, color="black")+ ggtitle("gene-wise correlation")+
  theme_classic()+
  theme(panel.background = element_rect(fill = "white", colour = "white"), panel.border = element_rect(colour = "black", fill=NA, size=0.8)) +
  theme(axis.text=element_text(size=12,colour="black"), axis.title=element_text(size=15, colour="black")) + ylab("CMOT")

p2_cor_compare <- p2_cor_compare +  geom_abline(intercept =0 , slope = 1, linetype = "dashed")+
  annotate("text", y = min(df$`CMOT`, na.rm = TRUE)+0.01, x = max(df$`MOFA+`, na.rm = TRUE) -0.2 , label = count_diagonal[1])+
  annotate("text", y = max(df$`CMOT`, na.rm = TRUE)+0.01, x = min(df$`MOFA+`, na.rm = TRUE) +0.2 , label = count_diagonal[2])+
  annotate("text", x = max(df$`MOFA+`, na.rm = TRUE)-0.2, y = min(df$`CMOT`, na.rm = TRUE)+0.1 , label = paste0("P = ",formatC(pval, format = "e", digits = 2)))

png(paste0("./GeneWiseCorrelations.png"), width = 330, height = 350, res=110)
p2_cor_compare
dev.off()
