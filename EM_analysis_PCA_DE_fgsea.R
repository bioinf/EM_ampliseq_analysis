stringsAsFactors=F
setwd("/set/to/your/work/dir")

library(DESeq2)
library(ggplot2)
library(ggrepel)
library(reshape)
library(fgsea)

## whole dataset DESeq-normalized 

tt <- read.table("Heterotopies.all_gene.rsem.counts",header=T,row.names=1)
tt2 <- tt[,!names(tt) %in% c("Symbol", "Gene_type")]
ann <- tt[,names(tt) %in% c("Symbol", "Gene_type")]
cond <- read.table("Conditions2.txt",row.names=1,header=T)

dds <- DESeqDataSetFromMatrix(countData=round(tt2,0),colData=cond,design = ~ Donor + Condition)
dds <- estimateSizeFactors(dds)
exp_norm <- counts(dds, normalized=TRUE)
exp_rlog <- assay(rlog(dds,blind=F,fitType="mean"))

ann_norm <- merge(ann, exp_norm, by = "row.names", all.x = T)
ann_rlog <- merge(ann, exp_rlog, by = "row.names", all.x = T)

write.table(ann_norm,file="EM_ampliseq.all_gene_DEseq2_norm.tsv",quote=F,sep="\t",row.names=F)
write.table(ann_rlog,file="EM_ampliseq.all_gene_rlog_norm.tsv",quote=F,sep="\t",row.names=F)



## Processing typical RNA-seq data, with DESeq2 based normalization and/or rlog transform for PCA. 

tt <- read.table("Ampliseq_heterotopy.filtered.rsem.counts",header=T,row.names=1)
cond <- read.table("Conditions2.txt",row.names=1,header=T)
head(tt)
cond

tt2 <- tt[,!names(tt) %in% c("Symbol", "Gene_type")]
ann <- tt[,names(tt) %in% c("Symbol", "Gene_type")]
ttn <- sweep(tt2, 2, colSums(tt2)/mean(colSums(tt2)), FUN="/")
colSums(ttn)
exp <- ttn[apply(ttn, 1, mean) > 10,]
dim(exp)

expl <- log2(exp+1)
boxplot(expl)
pca <- prcomp(t(expl))
pcat <- as.data.frame(pca$x[, 1:2])

pcat$Condition <- as.character(sapply(rownames(pcat), function(elt) cond[elt, 1])) 
pcat$Donor     <- as.character(sapply(rownames(pcat), function(elt) cond[elt, 2])) 

ggplot(pcat, aes(PC1, PC2)) + geom_point(aes(color = Condition), size =5) + geom_text_repel(data = pcat, aes(label = rownames(pcat)))
ggplot(pcat, aes(PC1, PC2)) + geom_point(aes(color = Donor), size =5) + geom_text_repel(data = pcat, aes(label = rownames(pcat)))

##################################################################################
###                    DESeq2 with and w/o doner effect                        ###
##################################################################################

tt2 <- tt2[, rownames(cond)]
dds_cond  <- DESeqDataSetFromMatrix(countData=round(tt2,0),colData=cond,design = ~ Condition)
dds_donor <- DESeqDataSetFromMatrix(countData=round(tt2,0),colData=cond,design = ~ Donor + Condition)
de_cond   <- DESeq(dds_cond)
de_donor  <- DESeq(dds_donor)

res_cond  <- results(de_cond,contrast=c("Condition","lesion","eutopic"))
res_donor <- results(de_donor,contrast=c("Condition","lesion","eutopic"))
resO_cond  <- as.data.frame(res_cond[order(res_cond$padj),])
resO_donor <- as.data.frame(res_donor[order(res_donor$padj),])

xx <- merge(resO_cond, ann, by = "row.names", all.x = T)
xx2 <- xx[, c("Row.names","Symbol","Gene_type","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj" )]
row.names(xx2) <- xx2$Row.names
xx2$Row.names <- NULL
resO_cond <- xx2[order(xx2$padj), ]

xx <- merge(resO_donor, ann, by = "row.names", all.x = T)
xx2 <- xx[, c("Row.names","Symbol","Gene_type","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj" )]
row.names(xx2) <- xx2$Row.names
xx2$Row.names <- NULL
resO_donor <- xx2[order(xx2$padj), ]

summary(res_cond)
summary(res_donor)

deg_cond <- resO_cond[complete.cases(resO_cond),]
deg_cond <- deg_cond[deg_cond$padj < 0.05,]
deg_cond <- deg_cond[order(deg_cond$log2FoldChange,decreasing = T),]

deg_donor <- resO_donor[complete.cases(resO_donor),]
deg_donor <- deg_donor[deg_donor$padj < 0.05,]
deg_donor <- deg_donor[order(deg_donor$log2FoldChange,decreasing = T),]


write.table(deg_cond,file="New_eutopic_vs_lesion.cond_only.tsv",quote=F,sep="\t")
write.table(deg_donor,file="New_eutopic_vs_lesion.with_donor.tsv",quote=F,sep="\t")

## let's get some plots - MA, volcano

plotMA(res_donor, main="DESeq2", alpha=0.05, ylim=c(-10,10))

set.seed(42)
genes <- resO_donor
genes$Significant <- ifelse(genes$padj < 0.05, "FDR < 0.05", "Not Sig")

ggplot(genes, aes(x = log2FoldChange, y = -log10(padj))) + geom_point(aes(color = Significant)) + scale_color_manual(values = c("red", "grey")) +
  theme_bw(base_size = 16) + geom_text_repel(data = subset(genes, abs(log2FoldChange) > 5),
  aes(label = Symbol),size = 5,box.padding = unit(0.35, "lines"),point.padding = unit(0.3, "lines"))

##################################################################################
###                                fgsea on both                               ###
##################################################################################

##gmt1 <- "/home/apredeus/Dropbox/GSEA/Msigdb_gmts_v5.1/c2.cp.kegg.v5.1.symbols.gmt"
##gmt2 <- "/home/apredeus/Dropbox/GSEA/Msigdb_gmts_v5.1/c2.cp.v5.1.symbols.gmt"
gmt1 <- "C:/users/apredeus/Dropbox/GSEA/Msigdb_gmts_v5.1/c2.cp.kegg.v5.1.symbols.gmt"
gmt2 <- "C:/users/apredeus/Dropbox/GSEA/Msigdb_gmts_v5.1/c2.cp.v5.1.symbols.gmt"

kegg   <- gmtPathways(gmt1)
cp     <- gmtPathways(gmt2)

rnk_cond  <- aggregate(stat ~ Symbol, data=resO_cond, function(x) ifelse(mean(x)>0,max(x),min(x)))
rnk_donor <- aggregate(stat ~ Symbol, data=resO_donor, function(x) ifelse(mean(x)>0,max(x),min(x)))
rnk_cond  <- setNames(rnk_cond$stat,rnk_cond$Symbol)
rnk_donor <- setNames(rnk_donor$stat,rnk_donor$Symbol)

fg_cp_cond   <- fgsea(cp, rnk_cond, minSize=15, maxSize=500, nperm=1000000)
fg_cp_donor  <- fgsea(cp, rnk_donor, minSize=15, maxSize=500, nperm=1000000)

source('C:/users/apredeus/Dropbox/R_scripts/My_functions.R')
gsea_table(fg_cp_cond,rnk_cond,cp)
gsea_table(fg_cp_donor,rnk_donor,cp)

cp_cond <- fg_cp_cond[order(NES,decreasing = T),]
cp_cond$leadingEdge <- vapply(cp_cond$leadingEdge, paste, collapse = ", ", character(1L))
write.table(cp_cond,"cp_eutopic_vs_lesion.cond_only.significant_pathways.tsv",sep = "\t",quote=F,row.names=F)

cp_donor <- fg_cp_donor[order(NES,decreasing = T),]
cp_donor$leadingEdge <- vapply(cp_donor$leadingEdge, paste, collapse = ", ", character(1L))
write.table(cp_donor,"cp_eutopic_vs_lesion.with_donor.significant_pathways.tsv",sep = "\t",quote=F,row.names=F)

## individual figures GSEA


plotEnrichment(cp[["KEGG_REGULATION_OF_ACTIN_CYTOSKELETON"]], rnk_donor) + labs(title="KEGG: Actin cytoskeleton")
plotEnrichment(cp[["REACTOME_METABOLISM_OF_LIPIDS_AND_LIPOPROTEINS"]], rnk_donor) + labs(title="Reactome: Metabolism of lipids and lipoproteins")
plotEnrichment(cp[["NABA_CORE_MATRISOME"]], rnk_donor) + labs(title="Naba: Core matrisome")
plotEnrichment(cp[["REACTOME_CELL_CYCLE"]], rnk_donor) + labs(title="Reactome: Cell cycle")
plotEnrichment(cp[["PID_P53_DOWNSTREAM_PATHWAY"]], rnk_donor) + labs(title="PID: p53 downstream pathway")
plotEnrichment(cp[["REACTOME_CELL_CYCLE_MITOTIC"]], rnk_donor) + labs(title="Reactome: Cell cycle mitotic")

## now some more bar-charts for the figures, and we will be all set!

tt1 <- read.table("Ampliseq.hallmark.tsv",row.names=1)
tt2 <- read.table("Ampliseq.non_hallmark.tsv",row.names=1)
tt3 <- read.table("DE_gene_counts.tsv",row.names=1)
tt4 <- read.table("Endometriosis - Differential genes  - donor_up_MsigDB.tsv",row.names=1)
tt5 <- read.table("Endometriosis - Differential genes  - donor_dn_MsigDB.tsv",row.names=1)
tt6 <- read.table("Endometriosis - Differential genes  - Go_bio_up.tsv",row.names=1,sep="\t")
tt7 <- read.table("Endometriosis - Differential genes  - Go_bio_dn.tsv",row.names=1,sep="\t")

tt6$log <- -log10(tt6$V6)
tt7$log <- -log10(tt7$V6)


ggplot(data=tt1, aes(x=reorder(row.names(tt1),tt1$V2),y=V2)) + geom_bar(stat="identity")+coord_flip()
ggplot(data=tt2, aes(x=reorder(row.names(tt2),tt2$V2),y=V2)) + geom_bar(stat="identity")+coord_flip()
ggplot(data=tt3, aes(x=reorder(row.names(tt3),tt3$V3),y=V3)) + geom_bar(stat="identity")+coord_flip()
ggplot(data=tt4, aes(x=reorder(row.names(tt4),tt4$V6),y=V6)) + geom_bar(stat="identity")+coord_flip()
ggplot(data=tt5, aes(x=reorder(row.names(tt5),tt5$V6),y=V6)) + geom_bar(stat="identity")+coord_flip()
ggplot(data=tt6, aes(x=reorder(row.names(tt6),tt6$log),y=log)) + geom_bar(stat="identity")+coord_flip()
ggplot(data=tt7, aes(x=reorder(row.names(tt7),tt7$log),y=log)) + geom_bar(stat="identity")+coord_flip()




## fgsea has no stat power - let's try and redo with all the genes. 
### ATTENTION: THE APPROACH BELOW DOES NOT WORK
### It only results in smallest pathways being selected - clearly false positive. 

exp_full <- read.table("Heterotopies.all_gene.rsem.counts",header=T,row.names=1)
head(exp_full)

yy <- exp_full[exp_full$Gene_type=="protein_coding",]
dim(yy)
yy2 <- yy[,!names(yy) %in% c("Symbol", "Gene_type")]
yy3 <- sweep(yy2, 2, colSums(yy2)/mean(colSums(yy2)), FUN="/")
colSums(yy3)
yy4 <- yy3
yy4[yy4 == 0] <- 1 ## neat
colSums(yy4)
ann2 <- yy[,names(yy) %in% c("Symbol", "Gene_type")]

yy4 <- yy4[, rownames(cond)]

dds_full  <- DESeqDataSetFromMatrix(countData=round(yy4,0),colData=cond,design = ~ Condition)
de_full   <- DESeq(dds_full)
res_full  <- results(de_full,contrast=c("Condition","lesion","eutopic"))
resO_full  <- as.data.frame(res_full[order(res_full$padj),])

xx <- merge(resO_full, ann2, by = "row.names", all.x = T)
xx2 <- xx[, c("Row.names","Symbol","Gene_type","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj" )]
row.names(xx2) <- xx2$Row.names
xx2$Row.names <- NULL
resO_full <- xx2[order(xx2$padj), ]

summary(res_full)

rnk_full   <- aggregate(stat ~ Symbol, data=resO_full, function(x) ifelse(mean(x)>0,max(x),min(x)))
rnk_full   <- setNames(rnk_full$stat,rnk_full$Symbol)

fg_cp_full <- fgsea(cp, rnk_full, minSize=15, maxSize=500, nperm=1000000)

gsea_table(fg_cp_full,rnk_full,cp)

cp_full <- fg_cp_full[order(NES,decreasing = T),]
cp_full$leadingEdge <- vapply(cp_full$leadingEdge, paste, collapse = ", ", character(1L))
write.table(cp_full,"full_eutopic_vs_lesion.cond_only.significant_pathways.tsv",sep = "\t",quote=F,row.names=F)

## just for the record: this approach fails miserably :)
