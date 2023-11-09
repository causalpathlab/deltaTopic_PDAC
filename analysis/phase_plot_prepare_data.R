library(reticulate)
library(Matrix)
library(data.table)
library(ggplot2)
library(argparse)
library(magrittr)
library(ggpubr)
library(data.table)
source("util.R")
source("util_analysis.R")
use_virtualenv("/home/BCCRC.CA/yzhang/DisNet", require = TRUE)

sc <- import("scanpy")
os <- import("os")
scipy <- import("scipy")
parser <- ArgumentParser()

parser$add_argument("--SavePath",
    help = "relative path to save folder")


args <- parser$parse_args()

SaveFolderPath = "models/BDeltaTopic_allgenes_ep2000_nlv32_bs1024_combinebyadd_lr0.01_train_size1.0_pip0rho_0.1_pip0delta_0.1_klbeta_10.0v1"

DataDIR <- os$path$join(os$path$expanduser('~'), "projects/data")

adata_spliced <- sc$read_h5ad(os$path$join(DataDIR,'CRA001160/final_CRA001160_spliced_allgenes.h5ad'))
adata_unspliced <- sc$read_h5ad(os$path$join(DataDIR,'CRA001160/final_CRA001160_unspliced_allgenes.h5ad'))

sc$pp$normalize_per_cell(adata_spliced)
sc$pp$normalize_per_cell(adata_unspliced)
sc$pp$log1p(adata_spliced)
sc$pp$log1p(adata_unspliced)

S = Matrix(adata_spliced$X$toarray(), sparse=TRUE) 
colnames(S) <- adata_spliced$var$gene
rownames(S) <- adata_spliced$obs$sample_id 
U = Matrix(adata_unspliced$X$toarray(), sparse=TRUE) 
colnames(U) <- adata_unspliced$var$gene
rownames(U) <- adata_unspliced$obs$sample_id

rm(adata_spliced, adata_unspliced)

#my_gene_list_df = fread(paste0(SaveFolderPath, paste0("/","delta","_topK_genes.csv.gz"))) %>% as.data.frame()

# construct the gene list from
## 1. top genes from gsea
result_IMMUNESIG_delta <- readRDS(paste0(SaveFolderPath, "/gsea/IMMUNESIG_delta.rds"))
result_HALLMARK <- readRDS(paste0(SaveFolderPath, "/gsea/HALLMARK_delta.rds"))
result_KEGG <- readRDS(paste0(SaveFolderPath, "/gsea/KEGG_delta.rds"))

fgsea.results.aggregata <- tibble()

fgsea.results.aggregata <- bind_rows(fgsea.results.aggregata,present.fgsea.result.alltopics(result_KEGG, N_pathways = 1, N_genes = 10))
fgsea.results.aggregata <- bind_rows(fgsea.results.aggregata,present.fgsea.result.alltopics(result_HALLMARK, N_pathways = 1, N_genes = 10))
fgsea.results.aggregata <- bind_rows(fgsea.results.aggregata,present.fgsea.result.alltopics(result_IMMUNESIG_delta, N_pathways = 1, N_genes = 10))

fgsea.results.aggregata <- fgsea.results.aggregata %>% as.data.table
fgsea.results.aggregata[,.SD[which.min(padj)], by = topic] %>%  filter(padj < 0.1) -> gsea_top
gsea_top$topGenes
my_gene_list = c()
for(j in 1:length(gsea_top$topGenes)){
    gsea_top$topGenes[j] %>% strsplit(", ") %>% unlist() -> tmp
    my_gene_list = c(my_gene_list, tmp)
}

## 2. top genes from loadings
dt.delta <- readRDS(paste0(SaveFolderPath, "/dt_delta.rds"))

dt.delta.top <- dt.delta %>% arrange(desc(weight)) %>% 
  group_by(Var1) %>%
  slice(1:10) %>%
  select(c("Var1", "weight", "variable"))

my_gene_list <- c(my_gene_list, dt.delta.top$variable)
my_gene_list <- my_gene_list %>% unique()
## 3. top genes from DE ?

U_top_genes = U[,my_gene_list]
saveRDS(U_top_genes, file = 'figures/fig1b/data_top_u.rds')
S_top_genes = S[,my_gene_list]
saveRDS(S_top_genes, file = 'figures/fig1b/data_top_s.rds' )
