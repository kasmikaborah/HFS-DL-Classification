setwd("F:/PHD/chapter2/dataset/stomach_cancer/TCGA_GDC/Gene_expression/Normalization/")

count1 <- read.csv("preprocess_threshold_15.csv", row.names = 1)
count1
dim(count1)

dd1_t<-t(count1)
dd1_t
dim(dd1_t)

min_max_normalization <- function(x) {
  return((x - min(x)) / (max(x) - min(x)))
}


normalized_data <- apply(dd1_t, 2, min_max_normalization)
head(normalized_data)
dim(normalized_data)
View (normalized_data)
#abline(normalized_data=0.1)


counts<-t(normalized_data)
dim(counts)
head(counts)
class(counts)


write.table(counts,file="normalized_data_R.txt",quote=F,sep="\t",row.names=TRUE,col.name=T)


library(edgeR)
library(RColorBrewer)
library(gplots)
library(limma)

d0 <- DGEList(counts)
d0

metadata <- read.csv("metadata.csv")
head(metadata)
dim(metadata)

identical(metadata$Run, colnames(counts))
#counts <- counts[,metadata$Run](if rows are not same than run this line)

#preprocessing
d0 <- calcNormFactors(d0)
d0$samples

plotMDS(d0, col = as.numeric(metadata))

#write.table(d0,"normalizationd0.txt",sep = "\t")


mm <- model.matrix(~0 + Stage, data = metadata)
mm

#Voom transformation

y <- voom(d0, mm, plot = T)

#Fitting linear models in limma

fit <- lmFit(y, mm)
head(coef(fit))

#Specify which groups to compare using contrasts:
contr <- makeContrasts(StageNormal - StageTumor, levels = colnames(coef(fit)))

#Estimate contrast for each gene
tmp <- contrasts.fit(fit, contr)
tmp
tmp <- eBayes(tmp)
tmp


#toptable with adjustment like BH...
top.table <- topTable(tmp, adjust.method = "BH", sort.by = "P", n = Inf)
head(top.table,20)

length(which(top.table$adj.P.Val < 0.001))#(6298)

DEGlist<-top.table[which(top.table$adj.P.Val<0.001),]
head(DEGlist)

write.table(DEGlist,file="adjpval_0.001__minmax.txt",quote=F,sep="\t",row.names=TRUE,col.name=T)
