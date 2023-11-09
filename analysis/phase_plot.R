# the scatter plot of genes (spliced versus unspliced) topic specific
library(data.table)
library(ggplot2)
library(dplyr)
library(patchwork)
library(ComplexHeatmap)
library(stringr)

# fix the range of axis and y axis for each gene
# only remove （0,0） points, showing topic-specific genes with at least 100 cell
# geom_hex f colored

# grab the topic from the saved model
SaveFolderPath = "models/BDeltaTopic_allgenes_ep2000_nlv32_bs1024_combinebyadd_lr0.01_train_size1.0_pip0rho_0.1_pip0delta_0.1_klbeta_10.0v1"
topics_df <- fread(paste0(SaveFolderPath, "/topics.csv"))
rownames(topics_df) <- topics_df$V1
topics_df = topics_df[, -1]
topics_df$topics <- colnames(topics_df)[apply(topics_df, 1, which.max)]

dt.delta <- readRDS(paste0(SaveFolderPath, "/dt_delta.rds"))

# top 10 genes
dt.delta.top <- dt.delta %>% arrange(desc(weight)) %>%
  group_by(Var1) %>%
  slice(1:10) %>%
  select(c("Var1", "weight", "variable"))

# survival-significant topics
topics_to_show <- c(4, 7, 11, 26, 29, 6, 10)

topics_to_show_df <- dt.delta %>%
                  arrange(desc(weight)) %>%
                  group_by(Var1) %>%
                  slice(1:10) %>%
                  filter(Var1 %in% topics_to_show)
top_genes <- topics_to_show_df$variable %>% unique()

# weight matrix
heatmap_df_delta  <- dt.delta %>% filter(variable %in% top_genes) %>% filter(Var1 %in% topics_to_show) %>% select(c("Var1", "weight", "variable"))
heatmap_df_delta <- heatmap_df_delta %>% reshape(idvar = "variable", timevar = "Var1", direction = "wide")
heatmap_mat_delta <- heatmap_df_delta %>% select(-c("variable")) %>% data.matrix()
rownames(heatmap_mat_delta) <- heatmap_df_delta$variable
colnames(heatmap_mat_delta) <- str_replace(colnames(heatmap_mat_delta), "weight.", "topic")

str(heatmap_mat_delta)

pdf(paste0("survival/figures/heatmap/topgenes_delta.pdf"))
ht_delta <- Heatmap(heatmap_mat_delta, show_column_dend = FALSE, show_row_dend = FALSE,name = "weight")
dev.off()

# top genes in ETM
dt.ETM<- readRDS(paste0("models/BETM_spliced_ep2000_nlv32_bs1024_lr0.01_train_size0.9_pip0rho_0.1_klbeta_10.0v1", "/dt_rho.rds"))

#topics_to_show <- c(9, 10, 17, 18, 23) # 4, 7, 11, 26, 29

# BETM topics matched with significant deltaTopics
topics_to_show <- c(4, 21, 8)
topics_to_show_df <- dt.ETM %>%
                  arrange(desc(weight)) %>%
                  group_by(Var1) %>%
                  slice(1:10) %>%
                  filter(Var1 %in% topics_to_show)
#top_genes <- topics_to_show_df$variable %>% unique()
# use the same top genes as DeltaTopic without recomputing
# weight matrix
heatmap_df_ETM  <- dt.ETM %>% filter(variable %in% top_genes) %>% filter(Var1 %in% topics_to_show) %>% select(c("Var1", "weight", "variable"))
heatmap_df_ETM <- heatmap_df_ETM %>% reshape(idvar = "variable", timevar = "Var1", direction = "wide")
heatmap_mat_ETM <- heatmap_df_ETM %>% select(-c("variable")) %>% data.matrix()
rownames(heatmap_mat_ETM) <- heatmap_df_ETM$variable
colnames(heatmap_mat_ETM) <- str_replace(colnames(heatmap_mat_ETM), "weight.", "topic")

pdf(paste0("survival/figures/heatmap/topgenes_ETM.pdf"))
ht_ETM <- Heatmap(heatmap_mat_ETM, cluster_rows = F, cluster_columns = FALSE, show_column_dend = FALSE, show_row_dend = FALSE,name = "weight")
dev.off()

pdf(paste0("survival/figures/heatmap/matched_topgenes.pdf"))
ht_delta + ht_ETM
dev.off()

##### scatter plot
U <- readRDS('figures/fig2a_hex_scatter/data_top_u.rds')
S <- readRDS('figures/fig2a_hex_scatter/data_top_s.rds')


# i topic index
for(i in 1:32){
  subset_exp <- dt.delta.top %>% filter(Var1 == i)
  Var1 <- i
  print(paste0("topic-", Var1))
  plots <- list()
  # loop over all top genes in each topic
  df_topic <- data.frame()
  for(j in 1:nrow(subset_exp)){
    my_gene <- subset_exp$variable[j]
    print(my_gene)
    df <- data.frame(S = S[, my_gene], U = U[, my_gene], gene = my_gene)
    df <- cbind(df, topics_df)
    df <- df %>% filter(topics == paste0('topic_', Var1 - 1)) %>% select(c("S","U","gene"))
    # filter out (S,U) = (0,0)
    df <- df %>% filter(S > 0 | U > 0)
    # 95 quantile of S and U
    df <- df %>% filter(U >= quantile(df$U, 0.025) & U <= quantile(df$U, 0.975)) %>% filter(S >= quantile(df$S, 0.025) & S <= quantile(df$S, 0.975))
    if(nrow(df) > 0){
      df$topic <- paste0('topic_', Var1)
      df_topic <- rbind(df_topic, df)
    }
  }
  DT_topic <- df_topic %>% data.table()
  if(nrow(DT_topic)>0){
    DT_topic[, .SD[.N > 100], by = .(topic, gene)]-> DT_topic_plot
  }
  if(nrow(DT_topic_plot)>0){
      p <- ggplot(DT_topic_plot, aes(x = S, y = U)) +
            geom_hex(bins = 50) +#35393a
            #geom_density_2d(mapping = aes(x = S, y = U), data = DT_topic_plot[(S>0)&(U>0),,], color = "black", alpha = 0.25) +
            facet_wrap(~gene, ncol = 3) +
            #xlim(0,1) + ylim(0,1) +
            scale_fill_gradient(low="lightblue1",high="darkblue",trans="log10") + 
            geom_abline(intercept = 0, slope = 1, colour = "red") +
            theme_classic() +
            ggtitle(paste0("Top Genes In Topic ", Var1))
      ggsave(paste0('figures/fig2a_hex_scatter/Topic',Var1,"_",'TopGenes.pdf'), p)
    }
    topgenes <- DT_topic_plot$gene %>% unique()
    for(mygene in topgenes){
      DT_topic_gene_plot <- DT_topic_plot %>% filter(gene == mygene)
      if(nrow(DT_topic_gene_plot)>0){
        max_val <- max(max(DT_topic_gene_plot$U), max(DT_topic_gene_plot$S))
        p0 <- ggplot(DT_topic_gene_plot, aes(x = S, y = U)) +
              geom_hex(bins = 50) +#35393a
              scale_fill_gradient(low="lightblue1",high="darkblue",trans="log10") + 
              geom_abline(intercept = 0, slope = 1, colour = "red") +
              theme_classic() +
              xlim(0, max_val) + ylim(0, max_val) + # U and S are in the same range
              ggtitle(paste0("Topic ", Var1, "-", mygene))
        if(nrow(DT_topic_gene_plot[(S>0)&(U>0),,]) > 20){
          p0 <- p0 + geom_density_2d(mapping = aes(x = S, y = U), data = DT_topic_gene_plot[(S>0)&(U>0),,], color = "black", alpha = 0.25) 
        }
        ggsave(paste0('figures/fig2a_hex_scatter_same_range/Topic',Var1,"_",mygene,'.pdf'), p0)
      }
    }
  }
mygene <- topgenes[1]
mygene
#my_gene = "COL1A1"
dt.delta.top5 <- dt.delta %>% arrange(desc(weight)) %>%
  group_by(Var1) %>%
  slice(1:5) %>%
  select(c("Var1", "weight", "variable"))


my_gene_list = dt.delta.top5$variable  %>% unique()
for(my_gene in my_gene_list){
    df = data.frame(S = S[,my_gene], U = U[,my_gene])
    DT = cbind(df, topics_df) %>% as.data.table()
    
    DT_p2 <- DT[,.SD[U>0 & S >0 & (U>=quantile(U, 0.025) & U <= quantile(U, 0.975)) & (S>=quantile(U, 0.025) & S <= quantile(U, 0.975))],by = topics]
    
    N_topic_cutoff <- 50 
    DT_p2[, .SD[.N < N_topic_cutoff], by = .(topics) ] %>% select(topics)  %>% unique() -> small_topics
    DT_p2[topics %in% small_topics$topics,topics := paste0("topics N < ", N_topic_cutoff ),]
    p2 <- ggplot(DT_p2, aes(x = S, y = U, color = topics)) +
    geom_point(alpha = 0.5) +
    #geom_jitter() + 
    #geom_density2d(aes(color = topics)) +
    #geom_hex(bins = 100) +
    theme_bw() +
    geom_abline(intercept = 0, slope = 1, colour = "black", linetype = "dashed") + 
    ggtitle(my_gene)

    ggsave(paste0('figures/fig1c/top5/',my_gene,'_scatter.pdf'), p2)
}

#top2
dt.delta.top2 <- dt.delta %>% arrange(desc(weight)) %>%
  group_by(Var1) %>%
  slice(1:2) %>%
  select(c("Var1", "weight", "variable"))


my_gene_list = dt.delta.top2$variable  %>% unique()
for(my_gene in my_gene_list){
    df = data.frame(S = S[,my_gene], U = U[,my_gene])
    DT = cbind(df, topics_df) %>% as.data.table()
    
    DT_p2 <- DT[,.SD[U>0 & S >0 & (U>=quantile(U, 0.025) & U <= quantile(U, 0.975)) & (S>=quantile(U, 0.025) & S <= quantile(U, 0.975))],by = topics]
    
    N_topic_cutoff <- 50 
    DT_p2[, .SD[.N < N_topic_cutoff], by = .(topics) ] %>% select(topics)  %>% unique() -> small_topics
    DT_p2[topics %in% small_topics$topics,topics := paste0("topics N < ", N_topic_cutoff ),]
    p2 <- ggplot(DT_p2, aes(x = S, y = U, color = topics)) +
    geom_point(alpha = 0.5) +
    geom_jitter() + 
    #geom_density2d(aes(color = topics)) +
    #geom_hex(bins = 100) +
    theme_bw() +
    geom_abline(intercept = 0, slope = 1, colour = "black", linetype = "dashed") + 
    ggtitle(my_gene)

    ggsave(paste0('figures/fig1c/top2/',my_gene,'_scatter.pdf'), p2)
}


#top1
dt.delta.top1 <- dt.delta %>% arrange(desc(weight)) %>%
  group_by(Var1) %>%
  slice(1) %>%
  select(c("Var1", "weight", "variable"))


my_gene_list = dt.delta.top1$variable  %>% unique()
for(my_gene in my_gene_list){
    df = data.frame(S = S[,my_gene], U = U[,my_gene])
    DT = cbind(df, topics_df) %>% as.data.table()
    
    DT_p2 <- DT[,.SD[U>0 & S >0 & (U>=quantile(U, 0.025) & U <= quantile(U, 0.975)) & (S>=quantile(U, 0.025) & S <= quantile(U, 0.975))],by = topics]
    
    N_topic_cutoff <- 50 
    DT_p2[, .SD[.N < N_topic_cutoff], by = .(topics) ] %>% select(topics)  %>% unique() -> small_topics
    DT_p2[topics %in% small_topics$topics,topics := paste0("topics N < ", N_topic_cutoff ),]
    
    
    p2 <- ggplot(DT_p2, aes(x = S, y = U, color = topics)) +
    geom_point(alpha = 0.5) +
    geom_jitter() + 
    #geom_density2d(aes(color = topics)) +
    #geom_hex(bins = 100) +
    theme_bw() +
    geom_abline(intercept = 0, slope = 1, colour = "black", linetype = "dashed") + 
    ggtitle(my_gene)

    ggsave(paste0('figures/fig1c/top1/',my_gene,'_scatter.pdf'), p2)
}


featured_DT <- data.table()
my_gene_list = c("AFF3","B2M", "BTBD2", "COL1A1", "GSN", "IGFBP7","MSR1","SPARCL1","SYPL1")
for(my_gene in my_gene_list){
    df = data.frame(S = S[,my_gene], U = U[,my_gene])
    DT = cbind(df, topics_df) %>% as.data.table()
    
    DT_p2 <- DT[,.SD[U>0 & S >0 & (U>=quantile(U, 0.025) & U <= quantile(U, 0.975)) & (S>=quantile(U, 0.025) & S <= quantile(U, 0.975))],by = topics]
    
    N_topic_cutoff <- 50 
    DT_p2[, .SD[.N < N_topic_cutoff], by = .(topics) ] %>% select(topics)  %>% unique() -> small_topics
    DT_p2[topics %in% small_topics$topics,topics := paste0("topics N < ", N_topic_cutoff ),]
    DT_p2$gene <- my_gene

    featured_DT <- rbind(featured_DT, DT_p2)

    p2 <- ggplot(featured_DT, aes(x = S, y = U, color = topics)) +
    geom_point(alpha = 0.05, size = 1) +
    #geom_jitter() + 
    #geom_density2d(aes(color = topics)) +
    #geom_hex(bins = 100) +
    theme_bw() +
    facet_wrap(~gene, ncol =3) +
    geom_abline(intercept = 0, slope = 1, colour = "black", linetype = "dashed") 

}

p2

ggsave(paste0('figures/fig1c/featured_genes.pdf'), p2)
