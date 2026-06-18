
list.files("data_raw")
list.files()
setwd("Autism_GWAS")
list.files("data_raw")
gwas <- read.delim("data_raw/ASD_SPARK_EUR_iPSYCH_PGC.tsv")
head(gwas)
dim(gwas)
colnames(gwas)
min(gwas$P.value, na.rm = TRUE)        #minimum P value
sum(gwas$P.value < 5e-8, na.rm = TRUE) #number of GW Significant variants
sig_snps <- subset(gwas, P.value < 5e-8)
write.csv(sig_snps,
          "data_processed/significant_variants.csv",
          row.names = FALSE)
head(sig_snps$MarkerName)
sig_snps$rsid <- sub(".*_", "", sig_snps$MarkerName)
head(sig_snps$rsid, 10)
write.csv(
  sig_snps,
  "data_processed/significant_variants_with_rsid.csv",
  row.names = FALSE
)
length(unique(sig_snps$rsid))
table(sig_snps$Chromosome) #table of no. of significant SNPs corresponding to the chromosome
chr20 <- sig_snps[sig_snps$Chromosome == 20, ]
range(chr20$Position)
colnames(chr20)
class(chr20$Position)
nrow(chr20)
str(sig_snps$Chromosome)
head(sig_snps$Chromosome)
sig_snps[which.min(sig_snps$P.value), ]
library(BiocManager)

BiocManager::install("biomaRt")
install.packages("qqman")
colnames(gwas)
library(qqman)
manhattan_data <- data.frame(
  SNP = gwas$MarkerName,
  CHR = gwas$Chromosome,
  BP  = gwas$Position,
  P   = gwas$P.value
)
manhattan(
  manhattan_data,
  main = "Autism GWAS",
  genomewideline = -log10(5e-8),
  suggestiveline = -log10(1e-5)
)
head(manhattan_data$CHR) #To check numerical value or string
unique(manhattan_data$CHR)
class(gwas$Chromosome)
class(gwas$P.value)
class(gwas$Position)
exists("manhattan_data")
str(manhattan_data)
dim(manhattan_data)
test_data <- manhattan_data[1:100000, ]
library(qqman)

manhattan(test_data)
sum(is.na(manhattan_data$P))
sig_data <- subset(manhattan_data, P < 5e-8)

nrow(sig_data)
library(qqman)
manhattan(sig_data)  #plots only genome-wide significant variants
set.seed(123)

plot_data <- manhattan_data[sample(nrow(manhattan_data), 500000), ] #plotting only the 1st 500000 points

manhattan(plot_data,
          main = "Autism GWAS Overview")
manhattan(sig_data,
          main = "Genome-wide Significant Autism Loci")             #plotting only the significant points
top10 <- sig_data[order(sig_data$P), ][1:10, ]

top10[, c("SNP", "CHR", "BP", "P")]
range(sig_data$BP[sig_data$CHR == 20])
sig_data[which.min(sig_data$P), ]
library(biomaRt)

ensembl <- useEnsembl(
  biomart = "genes",
  dataset = "hsapiens_gene_ensembl"
)

genes_20 <- getBM(
  attributes = c(
    "hgnc_symbol",
    "chromosome_name",
    "start_position",
    "end_position",
    "gene_biotype"
  ),
  filters = c(
    "chromosome_name",
    "start",
    "end"
  ),
  values = list(
    20,
    21137199,
    21433399
  ),
  mart = ensembl
)
head(genes_20)
gene_list <- unique(genes_20$hgnc_symbol)
gene_list <- gene_list[gene_list != ""]
gene_list <- na.omit(gene_list)
length(gene_list)
head(gene_list)
range(sig_data$BP[sig_data$CHR == 1])
genes_1 <- getBM(
  attributes = c(
    "hgnc_symbol",
    "chromosome_name",
    "start_position",
    "end_position",
    "gene_biotype"
  ),
  filters = c(
    "chromosome_name",
    "start",
    "end"
  ),
  values = list(
    1,
    96042484,
    96136829
  ),
  mart = ensembl
)
head(genes_1)
gene_list <- unique(genes_1$hgnc_symbol)
gene_list <- gene_list[gene_list != ""]
gene_list <- na.omit(gene_list)
length(gene_list)
head(gene_list)
range(sig_data$BP[sig_data$CHR == 8])
genes_8 <- getBM(
  attributes = c(
    "hgnc_symbol",
    "chromosome_name",
    "start_position",
    "end_position",
    "gene_biotype"
  ),
  filters = c(
    "chromosome_name",
    "start",
    "end"
  ),
  values = list(
    8,
    10660810,
    10725996
  ),
  mart = ensembl
)
head(genes_8)
gene_list <- unique(genes_8$hgnc_symbol)
gene_list <- gene_list[gene_list != ""]
gene_list <- na.omit(gene_list)
length(gene_list)
head(gene_list)

range(sig_data$BP[sig_data$CHR == 17])
genes_17 <- getBM(
  attributes = c(
    "hgnc_symbol",
    "chromosome_name",
    "start_position",
    "end_position",
    "gene_biotype"
  ),
  filters = c(
    "chromosome_name",
    "start",
    "end"
  ),
  values = list(
    17,
    45815418,
    46273553
  ),
  mart = ensembl
)
head(genes_17)
gene_list <- unique(genes_17$hgnc_symbol)
gene_list <- gene_list[gene_list != ""]
gene_list <- na.omit(gene_list)
length(gene_list)
View(gene_list)
head(gene_list)
all_genes <- rbind(
  genes_1,
  genes_8,
  genes_17,
  genes_20
)
length(all_genes)
nrow(all_genes)
gene_list <- unique(all_genes$hgnc_symbol)

gene_list <- gene_list[!is.na(gene_list)]
gene_list <- gene_list[gene_list != ""]
head(gene_list, 20)

library(BiocManager)
library(org.Hs.eg.db)
library(clusterProfiler)
gene_df <- bitr(
  gene_list,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

head(gene_df)

go_bp <- enrichGO(
  gene = gene_df$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)
head(go_bp)
go_bp    #zero enrichment terms found due to very small gene set

save.image("autism_gwas_workspace.RData")
