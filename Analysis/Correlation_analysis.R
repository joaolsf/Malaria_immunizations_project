library(RColorBrewer)
library(ggplot2)
library(ggthemes)
library(ggrepel)
library(viridis)
library(longitudinal)
library(ggcorrplot)
library(Hmisc)
library(corrplot)
library(missMDA)

# Preprocssind df
multiplex <- read.csv("Correlation_analysis_inoclevel4.csv")
df_to_cluster <- multiplex[,-1]
rownames(df_to_cluster) <- multiplex$SAMPLE_ID

# Scale data
mydata <- scale(df_to_cluster)
mydata2 <- mydata[ , colSums(is.na(mydata))==0]

library(Hmisc)
cor <- rcorr(mydata2, type = "spearman") # rcorr Calcula o p-value das correlacoes

library(tabletools)
#calculate adjusted p-values of the correlations: https://rdrr.io/github/JMLuther/tabletools/man/rcorr_padjust.html
cor2 <- rcorr_padjust(cor) # BH by default

#Extract r and p values to plot a corr plot for the correlation of specific variables
cor$r 
cor2$r
write.csv(cor$r , "Correlation_person_rcorr_inoclevel4_rvalues.csv")

cor$P #unadjusted p-values
cor2$P #adjusted p-values using the BH correction method
write.csv(cor$P , "Correlation_person_rcorr_inoclevel4_pvalues.csv")
write.csv(cor2$P , "Correlation_person_rcorr_inoclevel4_adjustpvalues.csv")

# reading file (importar o arquivo)
multiplex2 <- read.csv("Correlation_spearman_rcorr_inoclevel4_rvalues_edited.csv")
multiplex3 <- read.csv("Correlation_spearman_rcorr_inoclevel4_adjustpvalues_edited.csv")

rownames(multiplex2) <- multiplex2$Features
rownames(multiplex3) <- multiplex3$Features

multiplex2 <- multiplex2[,-1]
multiplex2 <- as.matrix.data.frame(multiplex2)

multiplex3 <- multiplex3[,-1]
multiplex3 <- as.matrix.data.frame(multiplex3)

multiplex4 <- t(multiplex2)
multiplex5 <- t(multiplex3)

write.csv(multiplex4 , "Correlation_spearman_rcorr_inoclevel4_rvalues_transpose.csv")
write.csv(multiplex5 , "Correlation_spearman_rcorr_inoclevel4_adjustpvalues_transpose.csv")

# reading file (importar o arquivo)
multiplex2 <- read.csv("Correlation_spearman_rcorr_inoclevel4_rvalues_transpose.csv", header=T)
multiplex3 <- read.csv("Correlation_spearman_rcorr_inoclevel4_adjustpvalues_transpose.csv", header=T)

rownames(multiplex2) <- multiplex2$Features
rownames(multiplex3) <- multiplex3$Features

multiplex2 <- multiplex2[,-1]
multiplex2 <- as.matrix.data.frame(multiplex2)

multiplex3 <- multiplex3[,-1]
multiplex3 <- as.matrix.data.frame(multiplex3)

#plot the correlation of specific variables using the files above
col4 <- colorRampPalette(c("#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "#FFFFFF", "#FDDBC7", "#F4A582","#D6604D",  "#B2182B", "#67001F"))
library(corrplot)
cor_matrix <- corrplot(multiplex2, method = "square", p.mat = multiplex3, 
                        type = "full", is.corr=T, insig = 'label_sig', sig.level = c(0.001, 0.01, 0.05), pch.cex = 0.25,pch.col = 'white',
                        tl.cex=0.25, tl.col = "black", number.cex=0.2, number.font=0.7, diag=F,cl.cex = 0.3, cl.ratio = 0.9,
                        mar=c(0.05,0.05,0.05,0.05), col = col4(200), addgrid.col="light grey", add=F)



