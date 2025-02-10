#comandi da fare nel terminale di rstudio

pwd

ls -l

cd datiesame

mkdir -p rawdata

tar -xzvf data_rnaseq.tar.gz -C rawdata

cd rawdata

export PATH=${PATH}:/usr/local/bin

salmon

for sample in `ls *_1.fasta.gz`
do
index="/home/gitpod/datiesame/datasets_reference_only/trascriptome/chr21_transcripts_index"
name=${sample%_1.fasta.gz}
echo "quantifying $name"
salmon quant \
 -p 2 \
 -i $index \
 -l IU \
 -1 "${name}_1.fasta.gz" -2 "${name}_2.fasta.gz" \
 --validateMappings \
 -o "${name}.quant"
echo -e "$name done now\n"
done

#spostarsi adesso sulla console e settare la working directory su datiesame

library(DESeq2)
library(tximport)
library(tidyverse)
library(pheatmap)
library(clusterProfiler)
library(DOSE)
library(org.Hs.eg.db)

###################################
## PREPARE DATASET CONDITIONS #####
###################################

dataset <- tibble(
  sample = c("sample_01",
             "sample_02",
             "sample_03",
             "sample_04",
             "sample_05",
             "sample_06"),
  condition = c(rep("control", 3),
                rep("case", 3))
)
tx2gene <- read_tsv("/home/gitpod/datiesame/datasets_reference_only/trascriptome/gencode.v29.transcripts_no-vers_chr21_tx2gene.txt")


###################################
#### READ LOCAL FILES IN ##########
###################################

files <- file.path("/home/gitpod/datiesame/rawdata/", paste0(dataset$sample,".quant"), "quant.sf")
names(files) <- dataset$sample

txi <- tximport(files, type = "salmon", tx2gene = tx2gene)

colnames(txi$counts)
rownames(dataset) <- colnames(txi$counts)

dds <- DESeqDataSetFromTximport(txi, dataset, ~condition)

###################################
## PREFILTER MIN COUNTS >10 #####
###################################

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

### make sure base level is control
dds$condition <- relevel(dds$condition, ref = "control")


###################################
##### DIFFERENTIAL EXPRESSION #####
###################################

dds <- DESeq(dds)


###################################
## EXTRACT ANALYSIS RESULTS #####
###################################

res <- results(dds)
resOrdered <- res[order(res$pvalue),]

## writeLines(summary(res), "differential_expression_summary.txt")

plotMA(res, ylim=c(-3,3))

plotDispEsts(dds)

plotCounts(dds, gene=which.min(res$padj), intgroup="condition")

###################################
## WRITE RESULTS OF ANALYSIS #####
###################################

resdata <- as_tibble(resOrdered)
resdata$gene <- rownames(resOrdered)
write_tsv(resdata, "analysis_results.tsv")


############################################
## CLUSTERING ##############################
############################################

ntd <- normTransform(dds)
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("condition")])

pheatmap(assay(ntd)[select,],
         cluster_cols=FALSE, annotation_col=df$condition)

plotPCA(ntd, intgroup=c("condition"))



###################################
## EXTRACT SIGNIFICANT GENES #####
###################################

universe <- AnnotationDbi::select(org.Hs.eg.db,
                                  keys = keys(org.Hs.eg.db),
                                  columns = c('ENTREZID','SYMBOL','ENSEMBL','ENSEMBLTRANS'),
                                  keytype = 'ENTREZID')

sig_genes <- resdata$gene[which(resdata$padj<0.05)]
entrez_genes_sig <- unique(universe[which(universe$ENSEMBL %in% sig_genes),]$ENTREZID)

pvalue_ens_genes <- resdata$padj[which(resdata$padj<0.05)]
names(pvalue_ens_genes)<-sig_genes

pvalue_entrez_genes <- resdata$padj[which(resdata$padj<0.05)]
names(pvalue_entrez_genes) <- entrez_genes_sig


###################################
## ENRICH GO ANALYSIS #####
###################################

ego <- enrichGO( gene = sig_genes,
                 universe = unique(tx2gene$GENEID),
                 OrgDb = org.Hs.eg.db,
                 keyType = 'ENSEMBL',
                 ont = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.05)


ego

dotplot(ego, showCategory=10)

cnetplot(ego, foldChange=resdata$log2FoldChange[which(resdata$padj<0.5)])

ego_MF <- enrichGO(gene = sig_genes,
                   universe = unique(tx2gene$GENEID),
                   OrgDb = org.Hs.eg.db,
                   keyType = 'ENSEMBL',
                   ont = "MF",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05)

ego_MF

dotplot(ego_MF, showCategory = 10)

cnetplot(ego_MF, foldChange = resdata$log2FoldChange[which(resdata$padj<0.5)])

ego_CC <- enrichGO(gene = sig_genes,
                   universe = unique(tx2gene$GENEID),
                   OrgDb = org.Hs.eg.db,
                   keyType = 'ENSEMBL',
                   ont = "CC",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05)

ego_CC

dotplot(ego_CC, showCategory = 10)

cnetplot(ego_CC, foldChange = resdata$log2FoldChange[which(resdata$padj<0.5)])


###################################
## DISGNET ANALYSIS #####
###################################


gda <- read_tsv(gzfile("/home/gitpod/datiesame/datasets_reference_only/trascriptome/all_gene_disease_associations.tsv.gz"))

disease2gene=gda[, c("diseaseId", "geneId")]
disease2name=gda[, c("diseaseId", "diseaseName")]

disgnet = enricher(entrez_genes_sig, TERM2GENE=disease2gene, TERM2NAME=disease2name)


cnetplot(disgnet, foldChange=resdata$log2FoldChange[which(resdata$padj<0.5)])


save.image("deseq2_analysis.RData")