# library
library(ggplot2)
library(plotly)
library(dplyr)
library(hrbrthemes)

p<- read.csv("./LiuCancer_RNA_TO_ATAC_Correlations_10000genes_BinNorm_60_40_split_2.csv")

#p$Method <- factor(p$Method, levels = c('CMOT(p=25%)', 'CMOT(p=50%)', 'CMOT(p=75%)','CMOT(p=100%)', 'MOFA+', 'Seurat', 'TotalVI'))

p$Method <- factor(p$Method, levels = c('CMOT(p=25%)', 'CMOT(p=50%)', 'CMOT(p=75%)','CMOT(p=100%)', 'MOFA+', 'Seurat'))

#p$Method <- factor(p$Method, levels = c('CMOT(p=25%)', 'CMOT(p=50%)', 'CMOT(p=75%)','CMOT(p=100%)', 'Polarbear-coassay','Polarbear-semi'))

my_labels=c('CMOT\n(p=25%)', 'CMOT\n(p=50%)', 'CMOT\n(p=75%)','CMOT\n(p=100%)', 'MOFA+', 'Seurat')

#my_labels=c('CMOT\n(p=25%)', 'CMOT\n(p=50%)', 'CMOT\n(p=75%)','CMOT\n(p=100%)', 'Polarbear-\ncoassay','Polarbear')


p10 <-ggplotly( ggplot(p, aes(x = Method, y = Pearson.Correlation)) +
  geom_boxplot(fill="white")
  +labs(x = "Methods", y = "Cell-wise correlation between inferred and \nmeasured chromatin accessiblity") +
 theme(axis.text.x = element_text(face="bold", color="#000000", 
                                     size=18, angle=0),
          axis.text.y = element_text(face="bold", color="#000000", 
                                     size=18, angle=90),
          #panel.border = element_rect(linetype = "dashed", fill="NA"),
       #panel.background = element_rect(fill = "white", colour = "grey50"),
         axis.title.y = element_text(size = 16),
         axis.title.x = element_text(size = 16)
         )+   scale_x_discrete(labels= my_labels)
  )
p10
