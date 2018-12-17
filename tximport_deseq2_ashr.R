## Analysis of Whitmore human neutrophils timecourse w/ h. pylori infection
## Date: 10.30.2018
## Author: Michael Chimenti
## Organism: hg38
## Aligners: hisat2 / salmon
## Design: neutrophils at 0, 6, 24 hrs.  Plus/minus h.pylori infection
## Reps: 3

##########
## Imports
##########

#source("https://bioconductor.org/biocLite.R")
#biocLite("DEGreport")

#negative binomial GLM and related
library('DESeq2')
library('calibrate')
library('tximport')
library('readr')
#annotation
library('biomaRt')
library("AnnotationDbi")
library("org.Hs.eg.db")
#Exploratory analysis
library('tidyverse')
library('pcaExplorer')
#pathway and gene clusters
library('DEGreport')
#library(pathview)
#library(gage)
#library(gageData)
#library(ggplot2)

setwd("~/iihg/RNA_seq/whitmore/neutrophil_reanalysis_oct2018/")

###########
##Function Defs
###########

get_annotation <- function(dds, biomart_dataset, idtype){
  if(is.null(biomart_dataset))
    stop("Select a species to generate the corresponding annotation.
         To obtain a list, type mart = useMart('ensembl'), followed by listDatasets(mart).")
  
  mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL",
                  host="www.ensembl.org",
                  dataset=biomart_dataset,
                  ensemblRedirect = FALSE)  ## Ensembl biomart redirects to [web mirror]; [[force NO]]
  
  anns <- getBM(attributes = c(idtype, "external_gene_name", "description"),
                filters = idtype,
                values = rownames(dds),
                mart = mart)
  
  # keep and match with the ones that are actually there
  anns2 <- anns[match(rownames(dds), anns[, 1]), ]
  rownames(anns2) <- rownames(dds)
  # rename the columns rsp. add row names to be consistent with other function
  colnames(anns2) <- c("gene_id","gene_name","description")
  
  return(anns2)
}

## Volcano Plot function 
volcanoplot <- function (res, lfcthresh=2, sigthresh=0.05, main="Volcano Plot", legendpos="topright", labelsig=TRUE, textcx=1, ...) {
  with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main=main, ...))
  with(subset(res, padj<sigthresh ), points(log2FoldChange, -log10(pvalue), pch=20, col="red", ...))
  with(subset(res, abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="orange", ...))
  with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="green", ...))
  if (labelsig) {
    require(calibrate)
    with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), textxy(log2FoldChange, -log10(pvalue), labs=ext_gene, cex=textcx, offset=0.3, ...))
  }
  #legend(legendpos, xjust=1, yjust=1, legend=c(paste("FDR<",sigthresh,sep=""), paste("|LogFC|>",lfcthresh,sep=""), "both"), pch=20, col=c("red","orange","green"))
}

#######################################
## tximport > DESeq2 
#######################################
samples <- read.table("samples.csv", sep=',', header=TRUE)
rownames(samples) <- samples$sname
samples$group <- paste0(samples$time, samples$infect)

files <- file.path(getwd(), samples$sname, 'salmon', 'quant.sf')
names(files) <- samples$sname

tx2gene <- read_csv(file.path(getwd(), "tx2gene.csv"), col_names = FALSE)
tx2gene$X1 <- tx2gene$X1 %>%
  strsplit(split = '.', fixed = TRUE) %>%
  sapply( "[", 1)  ## obtuse sapply statement needed b/c of annoying way strsplit returns list of lists


txi <- tximport(files, type="salmon", tx2gene=tx2gene)
ddsTxi <- DESeqDataSetFromTximport(txi,
                                   colData = samples,
                                   design = ~ donor + prep + group)

ddsTxi <- ddsTxi[ rowSums(counts(ddsTxi)) > 5, ]
ddsTxi <- DESeq(ddsTxi)

anno <- get_annotation(ddsTxi, 'hsapiens_gene_ensembl','ensembl_gene_id')
anno <- na.omit(anno)

rldTxi <- rlog(ddsTxi, blind=FALSE)
pcaExplorer(dds=ddsTxi,annotation=anno,rlt=rldTxi)

## look at dispersion estimates 
plotDispEsts(ddsTxi)
plotMA(object = ddsTxi, alpha = 0.05)
plotPCA(object = rldTxi, intgroup = 'group')

#looking at PCs 3 and 4
rld_mat <- assay(rldTxi)
pca <- prcomp(t(rld_mat))
df <- cbind(samples, pca$x)
ggplot(df) + geom_point(aes(x=PC3,y=PC4, color = prep))

#############################################################################
### sample SWAP!  Confirmed with Laura Whitmore
### sample 6 and 7 are swapped (0 vs 6 hr)

samples <- read.table("samples_SWAP.csv", sep=',', header=TRUE)
rownames(samples) <- samples$sname
samples$group <- paste0(samples$time, samples$infect)

files <- file.path(getwd(), samples$sname, 'salmon', 'quant.sf')
names(files) <- samples$sname

tx2gene <- read_csv(file.path(getwd(), "tx2gene.csv"), col_names = FALSE)
tx2gene$X1 <- tx2gene$X1 %>%
  strsplit(split = '.', fixed = TRUE) %>%
  sapply( "[", 1)  ## obtuse sapply statement needed b/c of annoying way strsplit returns list of lists
txi <- tximport(files, type="salmon", tx2gene=tx2gene)

ddsTxi <- DESeqDataSetFromTximport(txi,
                                   colData = samples,
                                   design = ~ donor + prep + group)

ddsTxi <- ddsTxi[ rowSums(counts(ddsTxi)) > 5, ]
ddsTxi <- DESeq(ddsTxi)

anno <- get_annotation(ddsTxi, 'hsapiens_gene_ensembl','ensembl_gene_id')
anno <- na.omit(anno)

rldTxi <- rlog(ddsTxi, blind=FALSE)
pcaExplorer(dds=ddsTxi,annotation=anno,rlt=rldTxi)

###MNDA boxplot  ENSG00000163563

#normalized counts
mnda_counts <- counts(ddsTxi['ENSG00000163563',], normalized = TRUE)
m <- list(counts = as.numeric(mnda_counts), group = as.factor(samples$group))
m <- as.tibble(m)
p <- ggplot(m, aes(group, counts)) + geom_boxplot() + geom_jitter(width = 0.1)
p <- p + labs(x = "Experimental Group", y = "Normalized Counts ", title = "Normalized Expression of MNDA")
p <- p + scale_x_discrete(labels=c("00hn" = "PMN, 0hrs", "06hn" = "PMN, 6hrs",
                                   "06hy" = "PMN+Hp, 6hrs", "24hn" = "PMN, 24hrs", "24hy" = "PMN+Hp, 24hrs"))
p


traf1_counts <- counts(ddsTxi['ENSG00000056558',], normalized = TRUE)
m <- list(counts = as.numeric(traf1_counts), group = as.factor(samples$group))
m <- as.tibble(m)
q <- ggplot(m, aes(group, counts)) + geom_boxplot() + geom_jitter(width = 0.1)
q <- q + labs(x = "Experimental Group", y = "Normalized Counts ", title = "Normalized Expression of TRAF1")
q <- q + scale_x_discrete(labels=c("00hn" = "PMN, 0hrs", "06hn" = "PMN, 6hrs",
                                   "06hy" = "PMN+Hp, 6hrs", "24hn" = "PMN, 24hrs", "24hy" = "PMN+Hp, 24hrs"))
q

#rlog counts
mnda_rld <- assay(rldTxi["ENSG00000163563",])
m <- list(counts = as.numeric(mnda_rld), group = as.factor(samples$group))
m <- as.tibble(m)
p <- ggplot(m, aes(group, counts)) + geom_boxplot() + geom_jitter(width = 0.1)
p <- p + labs(x = "Experimental Group", y = "Rlog Variance Stabilized Expression ")
p

########




plotPCA(object = rldTxi, intgroup = c("infect","time"))

##############################################################################
##DE testing 

## Volcano Plot function 
volcanoplot2 <- function (res, lfcthresh=2, sigthresh=0.05, main="Volcano Plot", legendpos="topright", labelsig=TRUE, textcx=1, ...) {
  with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main=main, col="gray", ...))
  with(subset(res, padj<sigthresh ), points(log2FoldChange, -log10(pvalue), pch=20, col="gray", ...))
  with(subset(res, abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="gray", ...))
  with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="black", ...))
  if (labelsig) {
    require(calibrate)
    with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), textxy(log2FoldChange, -log10(pvalue), labs=ext_gene, cex=textcx, offset=0.7, ...))
  }
}

###  24 hr uninfected vs. 0 hr uninfected
res_24n_0n <- lfcShrink(ddsTxi, contrast = c("group","24hn","00hn"), type = 'ashr')
res_24n_0n <- na.omit(res_24n_0n)  #drop NA rows
res_24n_0n_sig <- res_24n_0n[res_24n_0n$padj < 0.001 & res_24n_0n$baseMean > 5.0 & abs(res_24n_0n$log2FoldChange) > 2,]
res_24n_0n_ord <- res_24n_0n_sig[order(res_24n_0n_sig$padj),]
res_24n_0n_ord$ext_gene <- anno[row.names(res_24n_0n_ord), "gene_name"]

png("volcano_24n_0n_DEgenes.png", 1200, 1500, pointsize=20, res=100)
volcanoplot2(res_24n_0n_ord, lfcthresh=2, 
            sigthresh=1e-15, main="24h PMN vs. 0h PMN, p-adjusted < 1e-15 (black circles)", textcx=.45, xlim=c(-10,10), ylim = c(4.5,50))
dev.off()

### 6 hr infected vs 6 hr uninfected

res_6y_6n <- lfcShrink(ddsTxi, contrast = c("group","06hy","06hn"), type = 'ashr')
res_6y_6n <- na.omit(res_6y_6n)  #drop NA rows
res_6y_6n_sig <- res_6y_6n[res_6y_6n$padj < 0.001 & res_6y_6n$baseMean > 5.0 & abs(res_6y_6n$log2FoldChange) > 2,]
res_6y_6n_ord <- res_6y_6n_sig[order(res_6y_6n_sig$padj),]
res_6y_6n_ord$ext_gene <- anno[row.names(res_6y_6n_ord), "gene_name"]

png("volcano_6hrinfect_6hourcntrl_DEgenes.png", 1200, 1500, pointsize=20, res=100)
volcanoplot2(res_6y_6n_ord, lfcthresh=2, 
            sigthresh=1e-20, main="6h PMN + HP vs. 6h PMN, p-adjusted < 1e-20 (black circles)", 
            textcx=.45, xlim=c(-8,8), ylim = c(4.5,160))
dev.off()

### 6 hour infected vs 0 hr un-infected

res_6y_0n <- lfcShrink(ddsTxi, contrast = c("group","06hy","00hn"), type = 'ashr')
res_6y_0n <- na.omit(res_6y_0n)  #drop NA rows
res_6y_0n_sig <- res_6y_0n[res_6y_0n$padj < 0.001 & res_6y_0n$baseMean > 5.0 & abs(res_6y_0n$log2FoldChange) > 2,]
res_6y_0n_ord <- res_6y_0n_sig[order(res_6y_0n_sig$padj),]
res_6y_0n_ord$ext_gene <- anno[row.names(res_6y_0n_ord), "gene_name"]

png("volcano_6hrinfect_0hourcntrl_DEgenes.png", 1200, 1500, pointsize=20, res=100)
volcanoplot2(res_6y_0n_ord, lfcthresh=2, 
            sigthresh=1e-20, main="6h PMN + HP vs. 0h PMN, p-adjusted < 1e-20 (black circles)", textcx=.45, xlim=c(-12,12), ylim = c(4.5,170))
dev.off()

###  24hr infected vs. 0 hr uninfected
res_24y_0n <- lfcShrink(ddsTxi, contrast = c("group","24hy","00hn"), type = 'ashr')
res_24y_0n <- na.omit(res_24y_0n)  #drop NA rows
res_24y_0n_sig <- res_24y_0n[res_24y_0n$padj < 0.001 & res_24y_0n$baseMean > 5.0 & abs(res_24y_0n$log2FoldChange) > 2,]
res_24y_0n_ord <- res_24y_0n_sig[order(res_24y_0n_sig$padj),]
res_24y_0n_ord$ext_gene <- anno[row.names(res_24y_0n_ord), "gene_name"]

png("volcano_24hrinfect_0hourcntrl_DEgenes.png", 1200, 1500, pointsize=20, res=100)
volcanoplot2(res_24y_0n_ord, lfcthresh=2, 
            sigthresh=1e-50, main="24h PMN + HP vs. 0h PMN, p-adjusted < 1e-50 (black circles)", textcx=.45, xlim=c(-12,12), ylim = c(4.5,200))
dev.off()



###  24hr infected vs. 24 hr uninfected

res_24y_24n <- lfcShrink(ddsTxi, contrast = c("group","24hy","24hn"), type = 'ashr')
res_24y_24n <- na.omit(res_24y_24n)  #drop NA rows
res_24y_24n_sig <- res_24y_24n[res_24y_24n$padj < 0.001 & res_24y_24n$baseMean > 5.0 & abs(res_24y_24n$log2FoldChange) > 2,]
res_24y_24n_ord <- res_24y_24n_sig[order(res_24y_24n_sig$padj),]
res_24y_24n_ord$ext_gene <- anno[row.names(res_24y_24n_ord), "gene_name"]

png("volcano_24hrinfect_24hourcntrl_DEgenes.png", 1200, 1500, pointsize=20, res=100)
volcanoplot2(res_24y_24n_ord, lfcthresh=2, 
            sigthresh=1e-40, main="24h PMN + HP vs. 24h PMN, p-adjusted < 1e-20 (black circles)", textcx=.4, xlim=c(-8,10), ylim = c(4.5,210))
dev.off()

###  24hr infected vs. 6 hr infected

res_24y_6y <- lfcShrink(ddsTxi, contrast = c("group","24hy","06hy"), type = 'ashr')
res_24y_6y <- na.omit(res_24y_6y)  #drop NA rows
res_24y_6y_sig <- res_24y_6y[res_24y_6y$padj < 0.001 & res_24y_6y$baseMean > 5.0 & abs(res_24y_6y$log2FoldChange) > 2,]
res_24y_6y_ord <- res_24y_6y_sig[order(res_24y_6y_sig$padj),]
res_24y_6y_ord$ext_gene <- anno[row.names(res_24y_6y_ord), "gene_name"]

png("volcano_24hrinfect_6hourinfect_DEgenes.png", 1200, 1500, pointsize=20, res=100)
volcanoplot2(res_24y_6y_ord, lfcthresh=2, 
            sigthresh=1e-15, main="24h PMN + HP vs. 6h PMN + HP, p-adjusted < 1e-15 (black circles)", textcx=.4, xlim=c(-7,7), ylim = c(4.5,50))
dev.off()

### gene lists 
mycols <- c("baseMean", "log2FoldChange", "padj", "ext_gene")

write.csv(x = res_24n_0n_ord[,mycols], file = "DE_genes_24hrcntrl_0hrcntrl_padj_0p001_log2FC_abs_2.csv")
write.csv(x = res_6y_6n_ord[,mycols], file = "DE_genes_6hrHP_6hrcntrl_padj_0p001_log2FC_abs_2.csv")
write.csv(x = res_6y_0n_ord[,mycols], file = "DE_genes_6hrHP_0hrcntrl_padj_0p001_log2FC_abs_2.csv")
write.csv(x = res_24y_0n_ord[,mycols], file = "DE_genes_24hrHP_0hrcntrl_padj_0p001_log2FC_abs_2.csv")
write.csv(x = res_24y_24n_ord[,mycols], file = "DE_genes_24hrHP_24hrcntrl_padj_0p001_log2FC_abs_2.csv")
write.csv(x = res_24y_6y_ord[,mycols], file = "DE_genes_24hrHP_6hrHP_padj_0p001_log2FC_abs_2.csv")


##################
## Heatmaps 


apop <- read.csv('apop_genes_PMN.csv', header = FALSE)

## this section queries Ensemble online database for gene names associated with transcript IDs
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host="www.ensembl.org", ensemblRedirect = FALSE)
t2g <- getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id", "external_gene_name"), mart = mart)
t2g <- rename(t2g, 'target_id' = 'ensembl_transcript_id', 'ens_gene' = 'ensembl_gene_id', 'ext_gene' = 'external_gene_name')

names(apop) <- "ext_gene"

apop <- left_join(apop, t2g, on = 'ext_gene')
apop <- unique(apop[,c("ens_gene","ext_gene")])

library("pheatmap")
library("RColorBrewer")
library("viridis")
df <- as.data.frame(colData(ddsTxi)[,c("infect","time")])

apop_tib <- as.data.frame(assay(rldTxi)) %>%
  rownames_to_column(var = "ens_gene") %>%
  as.tibble() %>%
  filter(ens_gene %in% apop$ens_gene)

apop_tib <- left_join(apop_tib, apop, by = "ens_gene")
apop_mat <- as.matrix(apop_tib[,2:16])
row.names(apop_mat) <- apop_tib$ext_gene


quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}

mat_breaks <- quantile_breaks(apop_mat, n = 11)

pheatmap(apop_mat, 
         breaks = mat_breaks,
         cluster_rows=TRUE, 
         show_rownames=TRUE,
         cluster_cols=TRUE, 
         annotation_col=df,
         color = viridis(10),
         fontsize = 9,
         drop_levels = TRUE,
         show_colnames = FALSE,
         cutree_rows = 5,
         treeheight_row = 20,
         treeheight_col = 20,
         fontsize_row = 14
         )


##################
## Timecourse analysis: DEGpatterns


####  Non infected gene clusters

genes <- row.names(res_24n_0n_ord)
rldTxi_noHP <- rldTxi[ , rldTxi$infect %in% 'n']
mat <- assay(rldTxi_noHP[genes,])
clusters <- degPatterns(mat, metadata = as.data.frame(colData(rldTxi_noHP)), time = 'time',
                        reduce = TRUE, cutoff = 0.85, minc = 20)


cluster1_genes_noHP <- clusters$df[clusters$df$cluster %in% "1",]$genes
cluster2_genes_noHP <- clusters$df[clusters$df$cluster %in% "2",]$genes

library(clusterProfiler)
clust24 <- unique(clusters$df$cluster)
for(i in clust24){
  group <- filter(clusters$df, cluster %in% as.character(i))
  group.df <- bitr(group$genes, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  ego <- enrichGO(gene = group.df$ENTREZID,
                  OrgDb = org.Hs.eg.db,
                  ont = "CC",
                  readable = TRUE)
  print("Enriched GO terms for cluster")
  print(i)
  print(head(ego))
}

for(i in clust24){
  group <- filter(clusters$df, cluster %in% as.character(i))
  group.df <- suppressMessages(bitr(group$genes, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db))
  ego <- enrichKEGG(gene = group.df$ENTREZID,
                    organism = 'hsa',
                    pvalueCutoff = 0.05)
  print(paste0("Enriched KEGG terms for cluster #:", i))
  print(head(ego))
}


####  Infected gene clusters
## subsetting to just the infected samples and the 00h control 
## want to track clusters of expression just in infected


genes <- row.names(res_24y_0n_ord[1:2000,])
rldTxi_HP <- rldTxi[ , rldTxi$infect %in% 'y' | rldTxi$time %in% '00h']
mat <- assay(rldTxi_HP[genes,])
clusters <- degPatterns(mat, metadata = as.data.frame(colData(rldTxi_HP)), time = 'time',
                        reduce = TRUE, cutoff = 0.85, minc = 20)


library(clusterProfiler)
clust24 <- unique(clusters$df$cluster)
for(i in clust24){
  group <- filter(clusters$df, cluster %in% as.character(i))
  group.df <- bitr(group$genes, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  ego <- enrichGO(gene = group.df$ENTREZID,
                  OrgDb = org.Hs.eg.db,
                  ont = "CC",
                  readable = TRUE)
  print("Enriched GO terms for cluster")
  print(i)
  print(head(ego))
}

for(i in clust24){
  group <- filter(clusters$df, cluster %in% as.character(i))
  group.df <- suppressMessages(bitr(group$genes, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db))
  ego <- enrichKEGG(gene = group.df$ENTREZID,
                    organism = 'hsa',
                    pvalueCutoff = 0.05)
  print(paste0("Enriched KEGG terms for cluster #:", i))
  print(head(ego))
}


