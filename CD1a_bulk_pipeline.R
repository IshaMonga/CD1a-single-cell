library( "DESeq2" )
library(ggplot2)

#https://lashlock.github.io/compbio/R_presentation.html
setwd("~/Documents/CD1a-work/New_Folder_With_Items/Bulk_RNAseq/hit-counts")
# countsName<-read.csv("./CD1a_counts_files.csv")

########Step1. Improt count file of 3stima nd 3 unstim samples#######
countData <- read.csv('./CD1a_counts_files.csv', header = TRUE, sep = ",")
head(countData)
names(countData)
names(countData)[2]<-gsub("[^921.AdJ.stim]","",names(countData)[2])
names(countData)[3]<-gsub("[^921.stim]","",names(countData)[3])
names(countData)[4]<-gsub("[^966.stim]","",names(countData)[4])
names(countData)[5]<-gsub("[^921.AdJ.unstim]","",names(countData)[5])
names(countData)[6]<-gsub("[^921.unstim]","",names(countData)[6])
names(countData)[7]<-gsub("[^966.unstim]","",names(countData)[7])
countData[1]

######Step2. Improt meta data file having condition information of the samples, this will be used for Differential Expression(DE) testing#####
metaData <- read.csv("./metadata_Bulk.csv")
head(metaData)
names(metaData)
names(metaData)[1]<-gsub("[^id]","",names(metaData)[1])



#####Step3: Create DEseq object using DESeqDataSetFromMatrix function####
dds <- DESeqDataSetFromMatrix(countData=countData, 
                              colData=metaData, 
                              design=~condition, tidy = TRUE)
##check it##
dds


#####Step4: Converting ENSG to gene symbols using biomaRt on the dds DEseq object ####
library(biomaRt)
ensg<-rownames(dds)
head(ensg)
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
ensg.biomart.out<-getBM(attributes=c('ensembl_gene_id','ensembl_transcript_id','hgnc_symbol'),filters = 'ensembl_gene_id', values = ensg, mart = ensembl)
write.csv(ensg.biomart.out, "./bulk_isha_analysis/ensg_gene_symbol.csv")

idx<-match(ensg,ensg.biomart.out$ensembl_gene_id)
ensg<-ensg.biomart.out$ensembl_gene_id[idx]
rownames(dds)<-ensg.biomart.out$hgnc_symbol[idx]
dds


              ##Pre-filtering#########
# While it is not necessary to pre-filter low count genes before running the DESeq2 functions, 
# there are two reasons which make pre-filtering useful: 
# by removing rows in which there are very few reads, we reduce the memory size of the dds data object, and we increase the speed of the transformation and testing functions within DESeq2. 
# Here we perform a minimal pre-filtering to keep only rows 
# that have at least 10 reads total. Note that more strict filtering to increase power is automatically applied via independent filtering on the mean of normalized counts within the results function.
keep <- rowSums(counts(dds)) >= 10

dds <- dds[keep,]
rownames(dds)
x<-na.omit(rownames(dds))
dds<-dds[x,]
rownames(dds)
dds$condition <- factor(dds$condition, levels = c("unstimulated","Stimulated"))
dds$condition <- relevel(dds$condition, ref = "unstimulated")
dds$condition <- droplevels(dds$condition)










##Now we’re ready to run DESEQ function to do DE testing##
dds <- DESeq(dds)
res <- results(dds)
res
res <- results(dds, name="condition_Stimulated_vs_unstimulated")
res <- results(dds, contrast=c("condition","Stimulated","unstimulated"))
head(res)

####Take a look at the results table##
resultsNames(dds)

####Log fold change shrinkage for visualization and ranking####
resLFC <- lfcShrink(dds, coef="condition_Stimulated_vs_unstimulated", type="apeglm")
resLFC
head(resLFC)
##Speed-up and parallelization thoughts##3

library("BiocParallel")
register(MulticoreParam(4))
resOrdered <- res[order(res$pvalue),]
summary(res)
sum(res$padj < 0.1, na.rm=TRUE)
res05 <- results(dds, alpha=0.05)
summary(res05)
sum(res05$padj < 0.05, na.rm=TRUE)

#####Independent hypothesis weighting######
library("IHW")
resIHW <- results(dds, filterFun=ihw)
summary(resIHW)
sum(resIHW$padj < 0.1, na.rm=TRUE)
metadata(resIHW)$ihwResult
plotMA(res, ylim=c(-2,2))
plotMA(resLFC, ylim=c(-2,2))
res[1]

# idx <- identify(res$baseMean, res$log2FoldChange)
# rownames(res)[idx]



####Data transformations and visualization######

####Data transformations######
#Extracting transformed values
vsd <- vst(dds, blind=FALSE)
rownames(vsd)
colnames(vsd)
head(assay(vsd), 3)
vst_counts<-assay(vsd)

write.csv(vst_counts,"./bulk_isha_analysis/vst_counts.csv")
# ####making this ENSG count to gene names at VST level (not necessary)###
# library(biomaRt)
# ensg<-rownames(vst_counts)
# head(ensg)
# ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
# ensg.biomart.out<-getBM(attributes=c('ensembl_gene_id','ensembl_transcript_id','hgnc_symbol'),filters = 'ensembl_gene_id', values = ensg, mart = ensembl)
# write.csv(ensg.biomart.out, "./bulk_isha_analysis/ensg_gene_symbol.csv")
# 
# idx<-match(ensg,ensg.biomart.out$ensembl_gene_id)
# ensg<-ensg.biomart.out$ensembl_gene_id[idx]
# rownames(vst_counts)<-ensg.biomart.out$hgnc_symbol[idx]
# vst_counts
# write.csv(vst_counts,"./bulk_isha_analysis/vst_counts_gene_Symbol.csv")



ntd <- normTransform(dds)
# BiocManager::install("vsn")
library("vsn")
meanSdPlot(assay(ntd))
meanSdPlot(assay(vsd)) #variance stabilized counts#
##Pheatmap:Data quality assessment by sample clustering and visualization###

####Data visualization######
library("pheatmap")
# dds$id
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:100]
df <- as.data.frame(colData(dds)[,c("condition","id")])
pheatmap(assay(vsd)[select,], cluster_rows=TRUE, show_rownames=TRUE,fontsize_row = 12,
         cluster_cols=FALSE, annotation_col=df)


##Heatmap of the sample-to-sample distances##
sampleDists <- dist(t(assay(vsd)))
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$id, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

##Principal component plot of the samples##
plotPCA(vsd, intgroup=c("condition", "id"))


par(mar=c(8,5,2,2))
boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2)
