source("util.R")
source("util_analysis.R")
library(ggplot2)
library(data.table)
library(magrittr)
library(reshape2)
library(dplyr)
library(patchwork)
library(RColorBrewer)
library(grDevices)
library(ggrastr)

Save_Path <- "models/BDeltaTopic_allgenes_ep2000_nlv32_bs1024_combinebyadd_lr0.01_train_size1.0_pip0rho_0.1_pip0delta_0.1_klbeta_10.0v1"
topics <- fread(paste0(Save_Path, "/topics.csv"))
K <- ncol(topics)-1
colnames(topics) <- c("V1", 1:K)

.dt.deltaTopic <- as.data.table(reshape2::melt(topics, id.vars=1)) %>%
    mutate(col = V1, row = variable, weight = value) %>%
    col.order(rev(1:K), ret.tab=TRUE) %>%
    as.data.table

Save_Path1 <- "models/BETM_spliced_ep2000_nlv32_bs1024_lr0.01_train_size0.9_pip0rho_0.1_klbeta_1.0v1"
topics <- fread(paste0(Save_Path1, "/topics.csv"))
K <- ncol(topics)-1
colnames(topics) <- c("V1", 1:K)

.dt.BETM_spliced <- as.data.table(reshape2::melt(topics, id.vars=1)) %>%
    mutate(col = V1, row = variable, weight = value) %>%
    col.order(rev(1:K), ret.tab=TRUE) %>%
    as.data.table

Save_Path2 <- "models/BETM_unspliced_ep2000_nlv32_bs1024_lr0.01_train_size0.9_pip0rho_0.1_klbeta_1.0v1"
topics <- fread(paste0(Save_Path2, "/topics.csv"))
K <- ncol(topics)-1
colnames(topics) <- c("V1", 1:K)

.dt.BETM_unspliced <- as.data.table(reshape2::melt(topics, id.vars=1)) %>%
    mutate(col = V1, row = variable, weight = value) %>%
    col.order(rev(1:K), ret.tab=TRUE) %>%
    as.data.table

#check if the cell names match
all(.dt.BETM_spliced$col == .dt.deltaTopic$col) 
all(.dt.BETM_unspliced$col == .dt.deltaTopic$col)

all(.dt.BETM_unspliced$weight == .dt.deltaTopic$weight)
all(.dt.BETM_spliced$weight == .dt.deltaTopic$weight)
all(.dt.BETM_unspliced$weight == .dt.BETM_spliced$weight)

mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(K)

# structure plot
plt.deltaTopic <- 
    ggplot(.dt.deltaTopic, aes(x=`col`, y=`weight`, fill=as.factor(`variable`), color=as.factor(`variable`))) +
    xlab("cells") + ylab("topic proportion") +
    theme(legend.position = "top", 
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
    rasterise(geom_bar(stat="identity", position="stack", size=0), dpi=300) +
    scale_fill_manual(values=mycolors) +
    scale_color_manual(values=mycolors)

.ggsave(paste0("models/structure_overlay_deltaTopic.pdf"), plot=plt.deltaTopic,
        width=7, height=3)

plt.BETM_spliced <- 
    ggplot(.dt.BETM_spliced, aes(x=`col`, y=`weight`, fill=as.factor(`variable`), color=as.factor(`variable`))) +
    xlab("cells") + ylab("topic proportion") +
    theme(legend.position = "none", 
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
    rasterise(geom_bar(stat="identity", position="stack", size=0), dpi=300) +
    scale_fill_manual(values=mycolors) +
    scale_color_manual(values=mycolors)

.ggsave(paste0("models/structure_overlay_BETM_spliced.pdf"), plot=plt.BETM_spliced,
        width=7, height=3)

plt.BETM_unspliced <- 
    ggplot(.dt.BETM_unspliced, aes(x=`col`, y=`weight`, fill=as.factor(`variable`), color=as.factor(`variable`))) +
    xlab("cells") + ylab("topic proportion") +
    theme(legend.position = "none", 
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
    rasterise(geom_bar(stat="identity", position="stack", size=0), dpi=300) +
    scale_fill_manual(values=mycolors) +
    scale_color_manual(values=mycolors)

.ggsave(paste0("models/structure_overlay_BETM_unspliced.pdf"), plot=plt.BETM_unspliced,
        width=7, height=3)

plt <- plt.deltaTopic/plt.BETM_spliced/plt.BETM_unspliced

#.ggsave(paste0("models/structure_overlay.png"), plot=plt,
#        width=7, height=9)

.ggsave(paste0("models/structure_overlay.pdf"), plot=plt,
        width=7, height=9)



.dt.deltaTopic[, .SD[which.max(weight)], by = V1] -> x1
.dt.BETM_spliced[, .SD[which.max(weight)], by = V1] -> x2
.dt.BETM_unspliced[, .SD[which.max(weight)], by = V1] -> x3
out <- table(x1$row, x2$row)

rownames(out) <- paste0("D", rownames(out))
colnames(out) <- paste0("E", colnames(out))

out2 <- table(x1$row, x3$row)

rownames(out2) <- paste0("D", rownames(out2))
colnames(out2) <- paste0("E", colnames(out2))

# Chord diagram
library(circlize)
pdf(paste0(Save_Path,"/figures/chordDiag_deltaETMspliced.pdf"))
set.seed(123)
circos.par(start.degree = 0)
chordDiagram(out, reduce = 0.005, big.gap = 20)
title("Topic Correspondence\ndeltaTopic and ETM(spliced)")
abline(h = 0, lty = 2, col = "#00000080")
dev.off()
circos.clear()

pdf(paste0(Save_Path,"/figures/chordDiag_deltaETMunspliced.pdf"))
set.seed(123)
circos.par(start.degree = 0)
chordDiagram(out2, reduce = 0.005, big.gap = 20)
title("Topic Correspondence\ndeltaTopic and ETM(unspliced)")
abline(h = 0, lty = 2, col = "#00000080")
dev.off()
circos.clear()

# heatmap
heatmap(out)
