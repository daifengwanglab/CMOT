library(plotly)

df <- read.csv("./CMOT_Final/PBMC/Cell_TruevsPred_22642.csv")

scaleFUN <- function(x) format(x,scientific = TRUE, digits=3)
#df$Name <- c('CD3','CD4','CD8a','CD14','CD15','CD16','CD56','CD19','CD25','CD45RA','CD45RO','PD','TIGIT','CD127')

p<-ggplot(df, aes(x = TrueGE, y = PredictedGE )) +geom_point(color = "black", alpha=.5, size = 4)+ ggtitle("Cell2") 
   #+geom_point(size = 5,  color="#7e6148")+geom_text(x = 0.08, y = 0.04, check_overlap = T)
   #+geom_smooth(method = 'lm',color = "black")
p<-p+geom_abline(intercept =0 , slope = 1, linetype = "dashed") 
    #labs(x = "True Gene Expression", y = "Predicted Gene Expression") +
p<-p+annotate("text", 
             label=paste0("r = ", round(with(df, cor.test(TrueGE, PredictedGE))$estimate, 3),
                          "\np = ", format(with(df, cor.test(TrueGE, PredictedGE))$p.value,scientific = TRUE, digits = 2)),
             x = max(df$TrueGE)-0.6, y = max(df$PredictedGE)-0.1, color="black", size=10)+theme_classic()+theme(plot.title = element_text(size = 26,  vjust = -0.5))
p<-p+theme(panel.background = element_rect(fill = "white", colour = "white"), panel.border = element_rect(colour = "black", fill=NA, size=0.8)) 
p<-p+theme(axis.text=element_text(size=20,colour="black"), axis.title=element_text(size=15, colour="black"))

ggsave("./CMOT_Final/PBMC/Cell_TruevsPred_22642.png", width = 5, height = 5)
p
#dev.off()

