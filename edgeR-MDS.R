Data <- read.csv("cm.csv",header=T)
CountData <- read.csv("cm.csv",header=T, row.names = 1)
GeneIDs <- Data[,1]
head(GeneIDs)
library(edgeR)
Condition <- c("Ctrl","Ctrl","Ctrl",
               "LIN28AKO","LIN28AKO","LIN28AKO",
               "LIN28BKO","LIN28BKO","LIN28BKO")

SampleID <- colnames(CountData)
MetaData <- data.frame(SampleID ,Condition)
str(MetaData)

targets <- MetaData ##Metadata
targets$Condition <- factor(targets$Condition)
str(targets)
levels(targets$Condition)
Condition <-factor(Condition,levels=c("Ctrl","LIN28AKO","LIN28BKO"))
targets$Condition <- factor(targets$Condition,levels=c("Ctrl","LIN28AKO","LIN28BKO"))

str(Condition)
head(CountData)
targets

design <- model.matrix(~Condition, data=targets)
y <- DGEList(counts=CountData, gene=GeneIDs)
y <- calcNormFactors(y)
y <- estimateGLMRobustDisp(y,design)

MDS_data <- plotMDS(y, top=20000)