
#We assume that the images have been analyzed using GenePix to produce a .gpr ﬁle for each array
#and that a targets ﬁle targets.txt has been prepared with a column containing the names of the .gpr ﬁles.
#The following protocol is for Two-channel arrays with common reference designs
library(tidyverse)
library(RColorBrewer)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(EnhancedVolcano)
library(viridis)
library(ggVennDiagram)
library(limma)

targets <- readTargets("Targets3.txt")

#Set up a ﬁlter so that any spot with a ﬂag of −99 or less gets zero weight.
#f <- function(x) as.numeric(x$Flags > -99)

#1- Read in the data.
#RG <- read.maimages(targets, source="genepix",wt.fun=f, columns=list(R="F635 Median" ,G="F532 Median",Rb="B635",Gb="B532" ))

#USE THIS
Pv_monkey <- read.maimages(targets, source="genepix", wt.fun=wtflags(weight=0,cutoff=-1))
                    
                   # columns=list(R="F635 Median" ,G="F532 Median",Rb="B635",Gb="B532" ))


#It is a good idea to look at your data to check that it has been read in correctly.
show(Pv_monkey)
#to see a print out of the ﬁrst few lines of data.
summary(Pv_monkey$R)
summary(Pv_monkey$F)
names(Pv_monkey$genes)

#2- Quality Assessment
plotMD(Pv_monkey)
plotMDS(Pv_monkey)

#Boxplots of the background intensities from each array
#will highlight any arrays unusually with high background intensities.
#If the plots suggest that some arrays are of lesser quality than others, 
#it may be useful to estimate array quality weights to be used in the linear model analysis, see Section 14.

boxplot(data.frame(log2(Pv_monkey$Gb)),main="Green background")
boxplot(data.frame(log2(Pv_monkey$Rb)),main="Red background")

#3- Background correction
#The following command implements a type of adaptive background correction. 
#This is optional but recommended for GenePix data and assessging differential gene expression
RG <- backgroundCorrect(Pv_monkey, method="normexp", offset=50)

#4- Within-Array Normalization
#If the RGList object has not already been background corrected, then normalizeWithinArrays will do this by default.
#Print-tip loess normalization:
#MA <- normalizeWithinArrays(RG)
#Global loess normalization:
MA <- normalizeWithinArrays(RG, method="loess")
plotDensities(MA)

cols <- MA$targets$Cy5
cols[cols=="PvSal_d-1_1st_inoc"] <- "blue"
cols[cols=="PvSal_d14_1st_inoc"] <- "blue"
cols[cols=="PvSal_d-1_2nd_inoc"] <- "yellow"
cols[cols=="PvSal_d14_2nd_inoc"] <- "yellow"
cols[cols=="PvSal_d28_2nd_inoc"] <- "yellow"
cols[cols=="PvSal_d-1_3rd_inoc"] <- "orange"
cols[cols=="PvSal_d14_3rd_inoc"] <- "orange"
cols[cols=="PvSal_d28_3rd_inoc"] <- "orange"
cols[cols=="PvSal_d-1_4th_inoc"] <- "red"
cols[cols=="PvSal_d14_4th_inoc"] <- "red"
cols[cols=="PvSal_d28_4th_inoc"] <- "red"
cols[cols=="PvSal_d-1_1st_control"] <- "brown"
cols[cols=="PvSal_d14_1st_control"] <- "brown"
boxplot(MA$M~col(MA$M),names=rownames(MA$targets),col=cols,xlab="Sample",ylab="M-values")

#5- Between-Array Normalization
MA.q <- normalizeBetweenArrays(MA, method="quantile")
plotDensities(MA.q)
plotMDS(MA.q)

show(MA.q$M)
show(MA.q$A)
show(MA.q$other)
show(MA.q$genes)
p# <- as.matrix.MAList(MA.q)
#write.csv((broom::tidy(p)), "PV_MN_Microarray_Results_Quantile_Loess_Norm_JL_run2_2.csv")

save(MA.q, file = "MA.q_run2")
load("Ma.q_run2")

#filtering
F>2B
quanlity_filter=function(gpr)
{
  s_filter = matrix(data=NA,ncol=ncol(gpr$R),nrow=nrow(gpr$R))
  n=ncol(gpr$R)
  for(i in 1:n)
  {
    tmp1 = cbind(gpr$R[,i],gpr$Rb[,i],gpr$G[,i],gpr$Gb[,i])
    
    w1 = apply(tmp1,1,function(x) x[1]<cutoff*x[2]&x[3]<cutoff*x[4])
    
    s_filter[,i]=as.numeric(w1)
  }
  return(s_filter)
}

set_NA= function(x)
{
  l = length(x)/2
  y1 = x[1:l]
  l = l+1
  y2 = x[l:length(x)]
  
  y1[y2>1]=NA
  return(y1)
}

cutoff = 2 #raw data F>2B
gf.filter = quanlity_filter(Pv_monkey)
gf.w = abs(Pv_monkey$weight-1)
gf.mask = gf.filter+gf.w
gf.final = t(apply(cbind(MA$M,gf.mask),1,function(x) set_NA(x)))

data.gf = cbind(Pv_monkey$genes,gf.final)
write.table(data.gf,file="PV_MN_Microarray_Results_Quantile_Loess_Norm_filtered_JL_run3.txt",quote=F,sep="\t",row.names=F)
MA.q$M <- gf.final
save(MA.q, file = "MA.q_run3")
load("Ma.q_run3")

#Weighted arrays
arrayw <- arrayWeights(MA.q)
barplot(arrayw, xlab="Array", ylab="Weight", col="white", las=2)
abline(h=1, lwd=1, lty=2)


#9.2- Interaction Models: 2 × 2 Factorial Designs or 3 x 3 factorial design
#Factorial designs are those where more than one experimental dimension is being varied and each combination of treatment conditions is observed.
#The two experimental dimensions or factors here are Inoculation level(1st, 2nd or 3rd) and Time point after inoculation (day -1, day 14, day 28).
#It will be convenient for us to collect the Strain/Treatment combinations into one vector as follows:
TS <- paste(targets$InocLevel_Strain, targets$PostInocDay, sep="")
TS

#It is especially important with a factorial design to decide what are the comparisons of interest. Analysing as for a Single Factor:
#We will assume here that the experimenter is interested in:
#1. which genes respond to stimulation in wild-type cells - in our case, which genes change expression across time during the first inoculation; second, third,etc
#Level1day28-level1day14; level1day14-level1day-1, etc
#2. which genes respond diﬀerently between each inoculation (2nd vs 1st; 3rd vs 2nd). It relates to the diﬀerence of diﬀerences, which is called the interaction term.
#(Level2day28-level2day14)-(Level1day28-level1day14)
#(level2day14-level2day-1)-(level1day14-level1day-1)
#(Level3day28-level3day14)-(Level2day28-level2day14)
#(Level3day14-level3day-1)-(Level2day14-level2day-1)

#All the approaches considered are equivalent and yield identical bottom-line results.
#The most basic approach is to ﬁt a model with a coeﬃcient for each of the four factor combinations a(in our case will be 9) nd then to extract the comparisons of interest as contrasts:
TS <- factor(TS, levels=c("NoNo", "Sal1d14","Sal2d14", "AMRU4d14", "AMRU4d28", "AMRU2d1", "AMRU2d14", "AMRU1d1", "AMRU1d14"))
design <- model.matrix(~0+TS)
colnames(design) <- levels(TS)
fit <- lmFit(MA.q, design)

cont.matrix <- makeContrasts(
  #Differences in the average gene expression in AMRU versus averaged Sal gene expression across all 3 inoculations:
  Diff1 = (AMRU4d28+AMRU4d14)/2-(Sal2d14+Sal1d14)/2,
  #Differences in average AMRU gene expression due to heterologous challenge (heterologous-challenged-monkeys vs non-immunized monkey):
  Diff2 = (AMRU4d28+AMRU4d14)/2-(AMRU1d14+AMRU1d1)/2,
  #compare the effect of heterologous challenge with previous 2 homologous challenges (3 Sal infections) vs heterologous challenge without previous homologous challenge (only 1 Sal infection).
  Diff3 = (AMRU4d28+AMRU4d14)/2-(AMRU2d14+AMRU2d1)/2,
  #compare the effect of AMRU with 1Sal infection versus Non Immmune
  Diff4 = (AMRU2d14+AMRU2d1)/2-(AMRU1d14+AMRU1d1)/2,
  #compare Sal 2nd inoculation vs Sal 1st inoculation
  Diff5 = Sal2d14-Sal1d14,
  levels=design)

fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
save(fit2, file = "fit2_run3")
load("fit2_run3")

#We can use topTable() to look at lists of diﬀerentially expressed genes for each of three contrasts, or else to look at all three contrasts simultaneously.
#This approach is recommended for most users, because the contrasts that are being tested are formed explicitly.

DGE1 = topTable(fit2, coef = "Diff4", adjust.method = "BH", number = "8208", sort.by = "p")
write.csv(topTable(fit2, coef = "Diff4", adjust.method = "BH", number = "8208", sort.by = "p"), "DGE_AMRU_1SalvsNonImmune.csv")

DGE1 <- read.csv("DGE_AMRU_averageSal_v2.csv")
options(ggrepel.max.overlaps = 50000)
EnhancedVolcano(toptable = DGE1,
                x = "logFC",
                y = "adj.P.Val",
                lab = DGE1$Name,
                xlim = c(-5, +5),
                ylim = c(0,5),
                pCutoff = 0.001,
                FCcutoff = 1.0,
                pointSize = 3,
                labSize = 3.0,
                labCol = 'black',
                boxedLabels = FALSE,
                col = c(brewer.pal(n = 4, name = "Spectral"), space = "RGB"),
                colAlpha = 4/5,
                title = "DGE_AMRU_averageSal_v2",
                legendLabels = c(
                  'Not significant',
                  'Fold change (but do not pass padj cutoff)',
                  'Pass padj cutoff',
                  'Pass both padj & fold change'),
                legendPosition = 'right',
                legendLabSize = 6.0,
                legendIconSize = 2.0,
                drawConnectors = TRUE,
                widthConnectors = 0.2,
                typeConnectors = "open",
                colConnectors = 'black'
) + theme_light()

DGE2 = topTable(fit2, coef = "Diff2", adjust.method = "BH", number = "8208", sort.by = "p")
write.csv(topTable(fit2, coef = "Diff2", adjust.method = "BH", number = "8208", sort.by = "p"), "DGE_AMRU_HeterovsNonImmune.csv")


#Over-ride colouring scheme with custom key-value pairs
#In certain situations, one may wish to over-ride the default colour scheme with their own colour-scheme, 
#such as colouring variables by pathway, cell-type or group. 
#This can be achieved by supplying a named vector as ‘colCustom’

DGE2 <- read.csv("Var_genes_DGE_AMRU_merged.csv")

# create custom key-value pairs for 'high', 'low', 'mid' expression by fold-change
# this can be achieved with nested ifelse statements
keyvals <- ifelse(
  DGE2$cond < 1, 'royalblue',
  ifelse(DGE2$cond < 2, 'gold',
         'black'))
keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'royalblue'] <- 'AMRU1Sal_vs_AMRUNaive'
names(keyvals)[keyvals == 'gold'] <- 'AMRU3Sal_vs_AMRUNaive'
names(keyvals)[keyvals == 'black'] <- 'AMRU3Sal_vs_AMRU1Sal'

EnhancedVolcano(toptable = DGE2,
                x = "logFC",
                y = "adj.P.Val",
                lab = DGE2$ID,
                #selectLab = rownames(cond)[which(names(keyvals) %in% c('AMRU1Sal_vs_AMRUNaive', 'AMRU3Sal_vs_AMRUNaive', 'AMRU3Sal_vs_AMRU1Sal'))],
                xlim = c(-5, +5),
                ylim = c(0,4),
                pCutoff = 0.001,
                FCcutoff = 1.0,
                pointSize = 3,
                labSize = 3.0,
                labCol = 'black',
                boxedLabels = FALSE,
                colCustom = keyvals,
                #col = c(brewer.pal(n = 4, name = "Spectral"), space = "RGB"),
                colAlpha = 4/5,
                title = "Var_genes_DGE_AMRU_merged",
                legendLabels = c(
                  'Not significant',
                  'Fold change (but do not pass padj cutoff)',
                  'Pass padj cutoff',
                  'Pass both padj & fold change'),
                legendPosition = 'right',
                legendLabSize = 6.0,
                legendIconSize = 2.0,
                drawConnectors = TRUE,
                widthConnectors = 0.2,
                typeConnectors = "open",
                colConnectors = 'black'
) + theme_light()

DGE3 <- read.csv("Var_genes_DGE_AMRU_Sal_merged.csv")

# create custom key-value pairs for 'high', 'low', 'mid' expression by fold-change
# this can be achieved with nested ifelse statements
keyvals <- ifelse(
  DGE3$cond < 1, 'black', 'royalblue')

#keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'royalblue'] <- 'AMRU_vs_Sal' 
names(keyvals)[keyvals == 'black'] <- 'Sal_level2vs1'
#names(keyvals)[keyvals == 'black'] <- 'AMRU3Sal_vs_AMRU1Sal'

EnhancedVolcano(toptable = DGE3,
                x = "logFC",
                y = "adj.P.Val",
                lab = DGE3$ID,
                #selectLab = rownames(cond)[which(names(keyvals) %in% c('AMRU1Sal_vs_AMRUNaive', 'AMRU3Sal_vs_AMRUNaive', 'AMRU3Sal_vs_AMRU1Sal'))],
                xlim = c(-5, +5),
                ylim = c(0,4),
                pCutoff = 0.001,
                FCcutoff = 1.0,
                pointSize = 3,
                labSize = 3.0,
                labCol = 'black',
                boxedLabels = FALSE,
                colCustom = keyvals,
                #col = c(brewer.pal(n = 4, name = "Spectral"), space = "RGB"),
                colAlpha = 4/5,
                title = "Var_genes_DGE_AMRU_Sal_merged",
                legendLabels = c(
                  'Not significant',
                  'Fold change (but do not pass padj cutoff)',
                  'Pass padj cutoff',
                  'Pass both padj & fold change'),
                legendPosition = 'right',
                legendLabSize = 6.0,
                legendIconSize = 2.0,
                drawConnectors = TRUE,
                widthConnectors = 0.2,
                typeConnectors = "open",
                colConnectors = 'black'
) + theme_light()

#Load lists and re-calculate adj p-values for the lsit of filtered genes (after removing hypothetical proteins)
E <- read.csv("DGE_Sal_level2vs1_v2.csv")
E2 <- p.adjust(E$P.Value, method = "BH")
E3 <- cbind(E, E2)
write.csv(E3, "DGE_Sal_level2vs1_v2.csv")

#Filter genes with logFC cutoff > 1.0
DGE <- read.csv("DGE_Sal_level2vs1_v2.csv")
diff = DGE %>% 
  as.data.frame() %>% 
  rownames_to_column("genes") %>% 
  filter(abs(logFC) > 1.0 ) %>% 
  arrange(desc(logFC), 
          desc(adj.P.Val))
head(diff)
write.csv(diff, "DGE_Sal_level2vs1_filtered.csv")


#Venn Diagrams for similarly DEGs
#1- Subset Diff expressed genes based on logFC cutoff = 1.0 and p.adj. values = 0.05
DGE <- read.csv("DGE_Sal_level2vs1_filtered.csv")
diff = DGE %>% 
  as.data.frame() %>% 
  rownames_to_column("Genes") %>% 
  filter(adj.P.Val < 0.05 & abs(logFC) > 1.0) %>% 
  arrange(desc(logFC), 
          desc(adj.P.Val))
head(diff)
write.csv(diff, "DGE_Sal_level2vs1_FC_adjp_filtered.csv")


#Load lists
E <- read.csv("List_AMRU_conditions_v2.csv")
E <- as.list(E, sep="")

library(ggVennDiagram)
# Default plot
ggVennDiagram(E)
# Remove labels background color
ggVennDiagram(E, label_alpha = 0)

# Change category names
# Change the gradient fill color
ggVennDiagram(E, color="black", size=2, lwd = 0.8, lty = 1, label_alpha = 0,
  category.names = c("AMRU_3SalvsNonImmune", "AMRU_1SalvsNonImmune", "AMRU_3Salvs1Sal")) +
  ggplot2::scale_fill_gradient(low = "#F4FAFE", high = "#4981BF")

library(ggvenn)
ggvenn(E)
# Change category names
# Change the fill color
names(E) <- c("AMRU_3SalvsNonImmune", "AMRU_1SalvsNonImmune", "AMRU_3Salvs1Sal")
ggvenn(
  E, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF"),
  stroke_size = 0.5, set_name_size = 4
)

"#868686FF", "#CD534CFF"

library(VennDiagram)
venn.diagram(E, filename = 'NULL')
overlap <- calculate.overlap(E)
write.csv(overlap[["a6"]], "Overlap_AMRU_1SalvsNonImmune_3Salvs1Sal.csv")











