library(ggpubr)


df <- read.csv("./CMOT_Final/A549/A549_DisPlusenrichment_CMOTvsMOFA2.csv")
df$MOFA. <- as.numeric(df$MOFA.)
df$X = str_wrap(df$X, width = 30)


ggbarplot(df, x = "X", y = "MOFA.",
          fill = "#7dbb41",               # change fill color by cyl
          color = "black",            # Set bar border colors to white
          palette = "jco",            # jco journal color palett. see ?ggpar
          width = 0.5,
          # sort.val = "asc",           # Sort the value in dscending order
         # sort.by.groups = TRUE,      # Sort inside each group
          x.text.angle = 0,           # Rotate vertically x axis texts
         xlab = "",
         ylab = "-log10(adj. pval)"
         
) + coord_flip()

ggsave("./CMOT_Final/A549/A549_DisPlusenrichment_MOFA2.png", width=5, height=4)