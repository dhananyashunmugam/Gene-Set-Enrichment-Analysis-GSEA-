setwd("H:/GSEA")
##ref:http://yulab-smu.top/biomedical-knowledge-mining-book/reactomepa.html
##GSEA using clusterprofiler
library(readxl)
library(org.Hs.eg.db)
library(clusterProfiler)
library(ReactomePA)
library(ggplot2)

##DR_TB_vs_DS_TB
DR_TB_vs_DS_TB <- as.data.frame(read_excel("fold_change_DR-TB_vs_DS-TB.xlsx",sheet=2))
i=which(is.na(DR_TB_vs_DS_TB$Entrez_ID))
DR_TB_vs_DS_TB=DR_TB_vs_DS_TB[-i,]


genelist=-log10(DR_TB_vs_DS_TB$pvalue) * sign(DR_TB_vs_DS_TB$log2FoldChange)
names(genelist)=DR_TB_vs_DS_TB$Entrez_ID

j=which(is.na(genelist))
genelist=genelist[-j]


##sorting in decreasing order
genelist<-sort(genelist, decreasing =TRUE)
##GSEA
set.seed(123)
DR_TB_vs_DS_TB_GSEA <- gsePathway(genelist, 
                pvalueCutoff = 0.05,
                pAdjustMethod = "BH",
                organism = "human",
                verbose = FALSE,
                minGSSize =15,
                maxGSSize =500,
                seed = TRUE)
DR_TB_vs_DS_TB_GSEA_symbol <- setReadable(DR_TB_vs_DS_TB_GSEA, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
write.csv(as.data.frame(DR_TB_vs_DS_TB_GSEA_symbol),"DR_TB_vs_DS_TB_GSEA.csv")

##pvalue less than 10
set.seed(123)
DR_TB_vs_DS_TB_GSEA_pvalue10 <- gsePathway(genelist, 
                                  pvalueCutoff = 10,
                                  pAdjustMethod = "BH",
                                  organism = "human",
                                  verbose = FALSE,
                                  minGSSize =15,
                                  maxGSSize =500,
                                  seed=FALSE)
DR_TB_vs_DS_TB_GSEA_pvalue10_symbol <- setReadable(DR_TB_vs_DS_TB_GSEA_pvalue10, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
write.csv(as.data.frame(DR_TB_vs_DS_TB_GSEA_pvalue10_symbol),"DR_TB_vs_DS_TB_GSEA_pval_lessthan10.csv")

###ridge plot
ridge_DR_TB_vs_DS_TB_GSEA<-ridgeplot(DR_TB_vs_DS_TB_GSEA,showCategory = 50,orderBy = "pvalue")+ labs(title = "GSEA - DR-TB vs DS-TB")




###########################fgsea###########################
##DR_TB_vs_DS_TB
library(data.table)
setwd("H:/fgsea")
library(fgsea)
pathways<-gmtPathways("c2.cp.reactome.v2023.1.Hs.entrez.gmt")
DR_TB_vs_DS_TB <- as.data.frame(read_excel("fold_change_DR-TB_vs_DS-TB.xlsx",sheet=2))
i=which(is.na(DR_TB_vs_DS_TB$Entrez_ID))
DR_TB_vs_DS_TB=DR_TB_vs_DS_TB[-i,]


genelist=-log10(DR_TB_vs_DS_TB$pvalue) * sign(DR_TB_vs_DS_TB$log2FoldChange)
names(genelist)=DR_TB_vs_DS_TB$Entrez_ID

j=which(is.na(genelist))
genelist=genelist[-j]


##sorting in decreasing order
genelist<-sort(genelist, decreasing =TRUE)
stats<-genelist
fgsea_DR_TB_vs_DS_TB<-fgsea(pathways,stats,minSize=4,maxSize = 800)
fwrite(fgsea_DR_TB_vs_DS_TB, file="fgsea_DR_TB_vs_DS_TB.tsv", sep="\t", sep2=c("", " ", ""))

##DR_TB_vs_HC
DR_TB_vs_HC <- as.data.frame(read_excel("fold_change_DR-TB_vs_HC.xlsx",sheet=2))
i=which(is.na(DR_TB_vs_HC$Entrez_ID))
DR_TB_vs_HC=DR_TB_vs_HC[-i,]


genelist=-log10(DR_TB_vs_HC$pvalue) * sign(DR_TB_vs_HC$log2FoldChange)
names(genelist)=DR_TB_vs_HC$Entrez_ID

j=which(is.na(genelist))
genelist=genelist[-j]


##sorting in decreasing order
genelist<-sort(genelist, decreasing =TRUE)
stats<-genelist
fgsea_DR_TB_vs_HC<-fgsea(pathways,stats,minSize=4,maxSize = 800)
fwrite(fgsea_DR_TB_vs_HC, file="fgsea_DR_TB_vs_HC.tsv", sep="\t", sep2=c("", " ", ""))


########Disease - DO - GSEA###########
library(DOSE)
##DR_TB_vs_DS_TB
DR_TB_vs_DS_TB <- as.data.frame(read_excel("fold_change_DR-TB_vs_DS-TB.xlsx",sheet=2))
i=which(is.na(DR_TB_vs_DS_TB$Entrez_ID))
#i=which(DR_TB_vs_DS_TB$Entrez_ID=="Novel Transcript")
DR_TB_vs_DS_TB=DR_TB_vs_DS_TB[-i,]


genelist=-log10(DR_TB_vs_DS_TB$pvalue) * sign(DR_TB_vs_DS_TB$log2FoldChange)
names(genelist)=DR_TB_vs_DS_TB$Entrez_ID

j=which(is.na(genelist))
genelist=genelist[-j]


##sorting in decreasing order
genelist<-sort(genelist, decreasing =TRUE)
##GSEA
set.seed(123)
DR_TB_vs_DS_TB_GSEA <- gseDO(genelist, 
                                  pvalueCutoff = 0.05,
                                  pAdjustMethod = "BH",
                                  verbose = FALSE,
                                  minGSSize =15,
                                  maxGSSize =500,
                                  seed = TRUE)
DR_TB_vs_DS_TB_GSEA_symbol <- setReadable(DR_TB_vs_DS_TB_GSEA, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
write.csv(as.data.frame(DR_TB_vs_DS_TB_GSEA_symbol),"DR_TB_vs_DS_TB_GSEA.csv")


###ridge plot
ridge_DR_TB_vs_DS_TB_GSEA<-ridgeplot(DR_TB_vs_DS_TB_GSEA,showCategory = 50,orderBy = "pvalue")+ labs(title = "GSEA - DR-TB vs DS-TB")


####################Disease enrichment - Disgenet##########################
library(DOSE)
##DR_TB_vs_DS_TB
DR_TB_vs_DS_TB <- as.data.frame(read_excel("fold_change_DR-TB_vs_DS-TB.xlsx",sheet=2))
i=which(is.na(DR_TB_vs_DS_TB$Entrez_ID))
#i=which(DR_TB_vs_DS_TB$Entrez_ID=="Novel Transcript")
DR_TB_vs_DS_TB=DR_TB_vs_DS_TB[-i,]


genelist=-log10(DR_TB_vs_DS_TB$pvalue) * sign(DR_TB_vs_DS_TB$log2FoldChange)
names(genelist)=DR_TB_vs_DS_TB$Entrez_ID

j=which(is.na(genelist))
genelist=genelist[-j]


##sorting in decreasing order
genelist<-sort(genelist, decreasing =TRUE)
##GSEA
set.seed(123)
DR_TB_vs_DS_TB_GSEA <- gseDGN(genelist, 
                             pvalueCutoff = 0.05,
                             pAdjustMethod = "BH",
                             verbose = FALSE,
                             minGSSize =15,
                             maxGSSize =500,
                             seed = TRUE)
DR_TB_vs_DS_TB_GSEA_symbol <- setReadable(DR_TB_vs_DS_TB_GSEA, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
write.csv(as.data.frame(DR_TB_vs_DS_TB_GSEA_symbol),"DR_TB_vs_DS_TB_GSEA.csv")


###ridge plot
ridge_DR_TB_vs_DS_TB_GSEA<-ridgeplot(DR_TB_vs_DS_TB_GSEA,showCategory = 50,orderBy = "pvalue")+ labs(title = "GSEA - DR-TB vs DS-TB")



