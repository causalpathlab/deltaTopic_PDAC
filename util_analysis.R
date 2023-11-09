require(msigdbr)
require(goseq)
require(fgsea)
require(data.table)
require(magrittr)
require(reshape2)
require(dplyr)
require(patchwork)
require(ggplot2)
require(grkmisc) # Truncates, trims, and wraps strings. Built for ggplot2 plots with long string labels.

# hypergeom test with set cutoffs
run.hyper.test <- function(gene.weights, geneset.dt, k=500, min.weight=1e-9) {
	.genes <- gene.weights$variable
    .topk <- gene.weights %>% arrange(desc(weight)) %>% dplyr::filter(weight >= min.weight)
	print(.topk %>% tail(n = 1) %>% select(weight))
	.topk.genes <- .topk %>% pull(variable)
	.dt <- as.data.table(geneset.dt)
	.dt <- .dt[gene_symbol %in% .genes, .(gene_symbol, gs_name)]

	ntot = length(unique(.genes))
	nthresh = length(.topk.genes)
	gs.size <- .dt[, .(m = length(unique(gene_symbol))), by = .(gs_name)]

	overlap.size <- .dt[gene_symbol %in% .topk.genes,
											.(q = length(unique(gene_symbol))),
											by = .(gs_name)]

	left_join(gs.size, overlap.size, by = "gs_name") %>% 
		mutate(`q` = if_else(is.na(`q`), 0, as.numeric(`q`))) %>% 
		mutate(n = `ntot` - `m`) %>% 
		mutate(k = nthresh) %>%  
		mutate(p.val = phyper(`q`, `m`, `n`, `k`, lower.tail=FALSE)) %>% 
		arrange(p.val) %>% 
		as.data.table
}

# rank-based test
make.gs.lol <- function(.dt) {
    .dt <- as.data.table(.dt) %>% unique()
    .list <-
        .dt[, .(gene = .(gene_symbol)), by = .(gs_name)] %>%
        as.list()
    .names <- .list$gs_name
    .ret <- .list$gene
    names(.ret) <- .names
    return(.ret)
}

present.fgsea.result <- function(result, N_pathways=10, N_genes=3){
    result[,
           topGenes := paste0(head(unlist(`leadingEdge`), N_genes), collapse=", "),
           by = .(pathway)]

    result %>%
        arrange(pval) %>%
        head(N_pathways) %>% 
        select(-c("leadingEdge","ES","NES"))
} 


present.fgsea.result.alltopics <- function(result, N_pathways=1, N_genes=3){
    result[,
           topGenes := paste0(head(unlist(`leadingEdge`), N_genes), collapse=", "),
           by = .(topic, pathway)]

    result %>% 
        group_by(topic) %>%
        top_n(-N_pathways, wt = pval) %>%
        select(-c("leadingEdge","ES","NES")) %>% 
        arrange(topic, pval)
}

geneset_heatmap <- function(param.dt, .db, weight_cutoff = 5){
    .show.param <-
        param.dt %>%
        dplyr::filter(variable %in% .db$gene) %>%
        mutate(weight = pmin(pmax(weight, -weight_cutoff), weight_cutoff)) %>% 
        order.pair(ret.tab = TRUE) %>%
        as.data.table

    .show.db <-
        .db[gene %in% param.dt$variable] %>%
        mutate(row = gene, col = `pathway`, weight = 1) %>%
        col.order(.ro = sort(unique(.show.param$col)), ret.tab=TRUE) %>%
        mutate(gene = `row`, `gene set` = `col`)

    theme_set(theme_classic() +
            theme(legend.key.height = unit(.2, "lines")) +
            theme(axis.text.x = element_blank()) +
            theme(axis.ticks.x = element_blank()))

    p1 <-
        ggplot(.show.param, aes(col, row, fill = weight)) +
        theme(legend.position = "top") +
        geom_tile() + xlab("genes") + ylab("topics") +
        scale_fill_distiller("topic-specific\ngene activities", palette = "RdBu", direction=-1)

    p2 <-
        ggplot(.show.db, aes(`gene`, `gene set`)) +
        ggtitle("marker genes") +
        geom_tile() + xlab("genes") +
        scale_y_discrete(
  label = format_pretty_string(truncate_at = 30)
)
    plt <- p1/p2
    plt
}

celltype_heatmap <- function(param.dt, .db, weight_cutoff = 5){
    .show.param <-
        param.dt %>%
        dplyr::filter(variable %in% .db$gene) %>%
        mutate(weight = pmin(pmax(weight, -weight_cutoff), weight_cutoff)) %>% 
        order.pair(ret.tab = TRUE) %>%
        as.data.table

    .show.db <-
        .db[gene %in% param.dt$variable] %>%
        mutate(row = gene, col = `cell type`, weight = 1) %>%
        col.order(.ro = sort(unique(.show.param$col)), ret.tab=TRUE) %>%
        mutate(gene = `row`, `cell type` = `col`)

    theme_set(theme_classic() +
            theme(legend.key.height = unit(.2, "lines")) +
            theme(axis.text.x = element_blank()) +
            theme(axis.ticks.x = element_blank()))

    p1 <-
        ggplot(.show.param, aes(col, row, fill = weight)) +
        theme(legend.position = "top") +
        geom_tile() + xlab("genes") + ylab("topics") +
        scale_fill_distiller("topic-specific\ngene activities", palette = "RdBu", direction=-1)

    p2 <-
        ggplot(.show.db, aes(`gene`, `cell type`)) +
        ggtitle("marker genes") +
        geom_tile() + xlab("genes")
    plt <- p1/p2 + plot_layout(heights = c(2.5, 1))
    plt
}

get_geneset.dt <- function(gsea_result, N_genes = 10, p_cutoff = 0.05) {
   a <- gsea_result %>% dplyr::filter(padj < p_cutoff) %>% present.fgsea.result.alltopics(N_pathways = 1, N_genes = N_genes) 
    sig_pathways <- a$pathway %>% unique()
    df <- data.frame(pathway = character(), gene = character())
    colnames(df) <- c("pathway", "gene")
    for(sig_path in sig_pathways){
        a %>% dplyr::filter(pathway == sig_path) %>% select(topGenes) -> b
        for(gene in strsplit(b$topGenes, ", ") %>% unlist()){
            new_df <- data.frame(pathway = sig_path, gene = gene)
            df <- rbind(df, new_df)
        }
    }
    .db <- unique(df) %>% data.table()
    return(.db)
}


fgsea_all_topics <- function(.show.param.full, geneset){
    n_topics <- .show.param.full$Var1 %>% unique() %>% length()
    result <- data.table()
    for (i in 1:n_topics) {
        subset <- .show.param.full %>% dplyr::filter(Var1 == i)
        scores <- subset$weight
        names(scores) <- subset$variable
        fgsea <- fgsea::fgsea(pathways = geneset, stats = scores[scores > 0], scoreType = "pos")
        result <- rbind(result, fgsea[, topic := i])
    }
    return(result)
}

