BiocManager::install("ComplexHeatmap", force = T)

library(rnaseqDTU)
library(tximport)
library(GenomicFeatures)
library(DRIMSeq)
library(stageR)
library(ggbeeswarm)
library(ggrepel)
library(scales)
library("ComplexHeatmap")


metabolites_all <- metabolites_all[!duplicated(metabolites_all$No),]
metabolites_all <- metabolites_all[,c("HMDB.ID.x","X1_1.x", "X1_2.x","X1_3.x","X1_4.x","X1_5.x",
                                      "X2_1","X2_2","X2_3","X2_4","X2_5",
                                      "X3_1","X3_2","X3_3","X3_4","X3_5",
                                      "X4_1","X4_2","X4_3","X4_4","X4_5",
                                      "X5_1","X5_2","X5_3","X5_4","X5_5")]
head(metabolites_all)
samps <- data.frame(0)
samps$sample_id <- NA
sample_id <- colnames(metabolites_all)[2:26]
samps <- as.data.frame(sample_id)
samps$condition <- NA
samps$condition[1:5] <- "CONTROL"
samps$condition[6:10] <- "COVID1"
samps$condition[11:15] <- "COVID2"
samps$condition[16:20] <- "COVID3"
samps$condition[21:25] <- "COVID4"

Boxplot_metabolites <- metabolites_all
Boxplot_metabolites$tid <- NA
Boxplot_metabolites$tid <- Boxplot_metabolites$HMDB.ID.x
Boxplot_metabolites$gid <- Boxplot_metabolites$HMDB.ID.x

Boxplot_metabolites$HMDB.ID.x <- NULL
Boxplot_metabolites <- Boxplot_metabolites[,c("tid","gid","X1_1.x", "X1_2.x","X1_3.x","X1_4.x","X1_5.x",
                                              "X2_1","X2_2","X2_3","X2_4","X2_5",
                                              "X3_1","X3_2","X3_3","X3_4","X3_5",
                                              "X4_1","X4_2","X4_3","X4_4","X4_5",
                                              "X5_1","X5_2","X5_3","X5_4","X5_5")]

plotExpression <- function(expData = NULL, geneID = NULL, samps = NULL, isProportion = FALSE) {
  colnames(expData)[1:2] = c("gid","tid")
  sub = subset(expData, gid == geneID)
  sub = reshape2::melt(sub, id = c("gid", "tid"))
  sub = merge(samps, sub, by.x = "sample_id", by.y = "variable")
  group = sub[,2]
  if(!isProportion) {
    sub$value = log(sub$value)
  }
  
  clrs = c("dodgerblue3", "maroon2",  "forestgreen", "darkorange1", "blueviolet", "firebrick2",
           "deepskyblue", "orchid2", "chartreuse3", "gold", "slateblue1", "tomato" , "blue", "magenta", "green3",
           "yellow", "purple3", "red" ,"darkslategray1", "lightpink1", "lightgreen", "khaki1", "plum3", "salmon")
  
  p = ggplot(sub, aes(dplyr::desc(tid), value, color = group, fill = group)) +
    geom_boxplot(alpha = 0.4, outlier.shape = NA, width = 0.8, lwd = 0.5) +
    stat_summary(fun = mean, geom = "point", shape = 5, size = 3, position=position_dodge(width = 0.8)) +
    scale_color_manual(values = clrs) + scale_fill_manual(values = clrs) +
    geom_quasirandom( size = 1, dodge.width = 0.8) + theme_bw() +
    ggtitle(geneID) + xlab("Comparison") +theme(axis.text.x = element_blank())
  
  if(!isProportion) {
    p = p + ylab("log(Expression)")
  } else {
    p = p + ylab("Counts")
  }
  p
}

head(Boxplot_metabolites)

plotExpression(Boxplot_metabolites,geneID = "HMDB0000714",samps,isProportion = TRUE)



unique(lista_sigów)

Boxplot_metabolites$tid[25] <- "HMDB0011507_1"
Boxplot_metabolites$tid[24] <- "HMDB0011509_1"

Boxplot_metabolites$gid[25] <- "HMDB0011507_1"
Boxplot_metabolites$gid[24] <- "HMDB0011509_1"

lista_sigów <- Boxplot_metabolites

lista_sigów <- lista_sigów$tid
i <- lista_sigów[1]


for (i in 1:length(lista_sigów)) {
  png(paste0("/dane/Marta_Majewska/Metabolom/Boxplots/",lista_sigów[i],".png"), width=6, height=6, units = "in", res = 300)
  p <- plotExpression(Boxplot_metabolites, lista_sigów[i] , samps, isProportion = TRUE)
  print(p)
  dev.off()
}

