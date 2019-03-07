
## get setup

.libPaths("")
.libPaths('/cluster/kappa/90-days-archive/bio/tools/R_libs/3.5/')
library(DESeq2)
library(vsn) 
library(ggplot2)
library(dplyr)
library(tidyverse)
library(ggrepel)
# library(DEGreport) not installed
library(RColorBrewer)
library(pheatmap)
library(org.Sc.sgd.db)
library(clusterProfiler)

# ## loading and exploring data
# readcounts<-read.table("/cluster/tufts/rt/rbator01/work/rnaseq_course_spring_2019/quantification/sacCerfeatureCounts_gene_results.txt", header = TRUE)
# row.names(readcounts)<-readcounts$Geneid
# readcounts <- readcounts[ , -c(1:6)]
# orig_names <- names(readcounts)
# names(readcounts) <- gsub(".*(WT|SNF2)(_rep[0-9]+).*", "\\1\\2", orig_names)
# head(readcounts)
# sample_info <- data.frame(condition = gsub( "_rep[0-9]+", "", names(readcounts)), row.names = names(readcounts) )
# print(sample_info)                          
# 
# ## instead, just load up the final matrix
# write.table(readcounts,"/cluster/tufts/rt/rbator01/work/rnaseq_course_spring_2019/sacCerfeatureCounts_gene_results.formatted.txt",sep="\t")
# write.table(sample_info,"/cluster/tufts/rt/rbator01/work/rnaseq_course_spring_2019/sample_info.txt",sep="\t")

## read in preprocessed count and meta data
readcounts<-read.table("/cluster/tufts/rt/rbator01/work/rnaseq_course_spring_2019/sacCerfeatureCounts_gene_results.formatted.txt",header=TRUE)
meta <- read.table("/cluster/tufts/rt/rbator01/work/rnaseq_course_spring_2019/sample_info.txt", header=TRUE)

## run analysis
dds  <- DESeqDataSetFromMatrix(countData = readcounts, colData = meta, design = ~ condition)
dds <- DESeq(dds)

## go through step by step
sizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)

## how to number correlate with size factor
colSums(counts(dds))

## Total number of normalized counts per sample
colSums(counts(dds, normalized=T))

## Plot dispersion estimates
plotDispEsts(dds)

## creating contrasts
#contrast <- c("condition", "level_to_compare", "base_level")
#results(dds, contrast = contrast, alpha = alpha_threshold)

contrast <- c("condition", "SNF2", "WT")
res_unshrunken <- results(dds, contrast=contrast, alpha = 0.05)
res <- lfcShrink(dds, contrast=contrast, res=res_unshrunken)

plotMA(res_unshrunken, ylim=c(-2,2))
plotMA(res, ylim=c(-2,2))

## explore the results table
mcols(res, use.names=T)
res %>% data.frame() %>% View()

## Summarize results
summary(res)

### Set thresholds
padj.cutoff <- 0.05
lfc.cutoff <- 0.58

### subsetting results by filtering

res_tb <- res %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>%
  filter(is.na(log2FoldChange) == FALSE ) %>%
  as_tibble()

res_tb %>% arrange(padj) %>% View()

sig_tb <- res_tb %>%
  filter(padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff)


## visualizing results

meta_tb <- meta %>% 
  rownames_to_column(var="samplename") %>% 
  as_tibble()

normalized_counts <- normalized_counts %>% 
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

## plotting a single gene


## add something here about finding the gene name

keytypes(org.Sc.sgd.db)

anno <- select(org.Sc.sgd.db,
               keys = res_tb$gene, keytype = "ORF", # to retrieve all genes: keys = keys(org.Sc.sgd.db)
               columns = c("GENENAME","ENSEMBL","ENTREZID")) %>% as_tibble()
subset(anno, GENENAME == "SNF2")

# simple
plotCounts(dds, gene="YOR290C", intgroup="condition") 

# Plotting the YOR290C normalized counts, using the samplenames (rownames of d as labels)

d <- plotCounts(dds, gene="YOR290C", intgroup="condition", returnData=TRUE)

ggplot(d, aes(x = condition, y = count, color = condition)) + 
  geom_point(position=position_jitter(w = 0.1,h = 0)) +
  geom_text_repel(aes(label = rownames(d))) + 
  theme_bw() +
  ggtitle("SNF2") +
  theme(plot.title = element_text(hjust = 0.5))

## plotting multiple genes

# Order results by padj values
top20_genes <- res_tb %>% 
  arrange(padj) %>% 	#Arrange rows by padj values
  pull(gene) %>% 		#Extract character vector of ordered genes
  head(n=20) 		#Extract the first 20 genes

# Then, we can extract the normalized count values for these top 20 genes:
# normalized counts for top 20 significant genes
top20_genes_counts <- normalized_counts %>%
  filter(gene %in% top20_genes)


### Extract normalized expression for significant genes from the OE and control samples (4:9), and set the gene column (1) to row names

norm_sig <- normalized_counts %>% 
  filter(gene %in% sig_tb$gene) %>% 
  data.frame() %>%
  column_to_rownames(var = "gene") 


#Now let’s draw the heatmap using pheatmap:
### Annotate our heatmap (optional)

### Set a color palette
heat_colors <- brewer.pal(6, "YlOrRd")

### Run pheatmap
pheatmap(norm_sig, 
         color = heat_colors, 
         cluster_rows = T, 
         show_rownames = F,
         annotation = meta, 
         border_color = NA, 
         fontsize = 10, 
         scale = "row", 
         fontsize_row = 10, 
         height = 20)

## testing for functional enrichment
res_with_ens <- inner_join(res_tb, anno, by=c("gene"="ORF"))     
sig_with_ens <- inner_join(sig_tb, anno, by=c("gene"="ORF"))    

sig_with_ens_only <- subset(sig_with_ens, is.na("ENTREZID") == FALSE)$ENTREZID
res_with_ens_only <- subset(res_with_ens, is.na("ENTREZID") == FALSE)$ENTREZID

ggo <- groupGO(gene         = sig_with_ens_only,
                OrgDb        = org.Sc.sgd.db,
                keyType      = "ENTREZID",
                ont          = "BP",
                level        = 3, readable = TRUE) ## what is level?


ego <- enrichGO(gene = sig_with_ens_only, 
                universe = res_with_ens_only,
                keyType = "ENTREZID",
                OrgDb = org.Sc.sgd.db, 
                ont = "BP", readable = TRUE)

keytypes(org.Sc.sgd.db)

cluster_summary <- data.frame(ego)
write.csv(cluster_summary, "ego.csv")

dotplot(ego, showCategory=8)
emapplot(ego, showCategory = 50)

head(cluster_summary)

foldchanges <- sig_tb$log2FoldChange

cnet<-cnetplot(ego, 
               categorySize="pvalue", 
               showCategory = 10, 
               foldChange=foldchanges, 
               vertex.label.font=6)
show(cnet)

