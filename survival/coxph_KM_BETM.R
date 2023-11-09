library(survival)
library(survminer)
library(data.table)
library(dplyr)
library(ComplexHeatmap)
library(EnsDb.Hsapiens.v79)
library(msigdbi)
source("util.R")

location <- "CA"
#CA
#US
target <- "rho"

model.path <- "models/BETM_spliced_ep2000_nlv32_bs1024_lr0.01_train_size0.9_pip0rho_0.1_klbeta_10.0v1"
# ETM on spliced
TCGA_data_DIR <- "/data/PDAC/bulk/"

if(location == "US"){
    exp_seq <- fread(paste0(TCGA_data_DIR, "processed/exp_seq.PAAD-US.tsv.gz"))
    donor <- fread(paste0(TCGA_data_DIR, "PAAD-US/donor.PAAD-US.tsv.gz"))
}

if(location == "AU"){
    exp_seq <- fread(paste0(TCGA_data_DIR, "processed/exp_seq.PACA-AU.tsv.gz"))
    donor <- fread(paste0(TCGA_data_DIR, "PACA-AU/donor.PACA-AU.tsv.gz"))

    # create mapping betwwen ensembl id and gene symbol
    ensembl.genes = exp_seq$gene_id
    geneID <- ensembldb::select(EnsDb.Hsapiens.v79, keys= ensembl.genes, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
    colnames(geneID) <- c("gene", "gene_id")
    exp_seq <- dplyr::left_join(exp_seq, geneID,  by = 'gene_id')
    exp_seq <- exp_seq[, gene_id := gene]
}

if(location == "CA"){
    exp_seq <- fread(paste0(TCGA_data_DIR, "processed/exp_seq.PACA-CA.tsv.gz"))
    donor <- fread(paste0(TCGA_data_DIR, "PACA-CA/donor.PACA-CA.tsv.gz"))

    # create mapping betwwen ensembl id and gene symbol
    ensembl.genes = exp_seq$gene_id
    geneID <- ensembldb::select(EnsDb.Hsapiens.v79, keys= ensembl.genes, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
    colnames(geneID) <- c("gene", "gene_id")
    exp_seq <- dplyr::left_join(exp_seq, geneID,  by = 'gene_id')
    exp_seq <- exp_seq[, gene_id := gene]
}

# creat time-to-event and status
donor[is.na(donor),] <- 0
donor[, time := donor_survival_time + donor_interval_of_last_followup]
donor$status  <- ifelse(donor$donor_vital_status == "alive", 0,
        ifelse(donor$donor_vital_status == "deceased", 1, NA))

donor <- donor[!is.na(status),]

# estimate dynamic (delta) topic proporitons

# choose which weight to run analysis for
if(target == "delta"){
    weights <- readRDS(paste0(model.path, "/dt_delta.rds"))
    weights <- weights %>% mutate(topic = row, gene = col) %>% dcast(gene ~ topic, value.var = "weight")
}else if (target == "rho"){
    weights <- readRDS(paste0(model.path, "/dt_rho.rds"))
    weights <- weights %>% mutate(topic = row, gene = col) %>% dcast(gene ~ topic, value.var = "weight")
}

# inner join on the common genes
exp_seq_wide <- exp_seq[, gene:= gene_id] %>% dcast(gene ~ icgc_donor_id, value.var = "raw_read_count", fun.aggregate = sum) %>% data.table()
merged.dt <- exp_seq_wide[weights, on = .(gene), nomatch = NULL]

# create weight matrix
weight.mat <- merged.dt %>% dplyr::select(paste0(1:32)) %>% as.matrix()
colnames(weight.mat) <- paste0("topic", 1:32)
# create expression matrix
exp.mat <- merged.dt %>% dplyr::select(-c("gene", paste0(1:32))) %>% as.matrix()
exp.mat <- scale(exp.mat) # zscore expression matrix
colnames(exp.mat) <- merged.dt %>% dplyr::select(-c("gene", paste0(1:32))) %>% colnames()
#heatmap(exp.mat)
##  Sample X gene * gene X Topic
weight.prop.hat <- t(exp.mat) %*% weight.mat
#Heatmap(weight.prop.hat)
# Create the heatmap annotation

weight.prop.hat.2 <- weight.prop.hat %>% as.data.table() %>% mutate(icgc_donor_id := rownames(weight.prop.hat))
merged.meta.hat <- weight.prop.hat.2[donor, on = .(icgc_donor_id), nomatch = NULL] 


mat_to_plt_heatmap <- merged.meta.hat %>% dplyr::select(paste0("topic", 1:32)) %>% as.matrix()
rownames(mat_to_plt_heatmap) <- merged.meta.hat$icgc_donor_id

#Heatmap(scale(mat_to_plt_heatmap))

ha <- HeatmapAnnotation(
  vital_status = merged.meta.hat$donor_vital_status
)


pdf(paste0("survival/figures/heatmap/", "heatmap_BETM_", location,"_", target,".pdf"))
Heatmap(t(mat_to_plt_heatmap),
    name = "deltaTopics",
    top_annotation = ha,
    column_names_gp = gpar(fontsize = 3), km = 4)
dev.off()

pdf(paste0("survival/figures/heatmap/", "heatmap_zscore_BETM_", location, "_", target,".pdf"))
Heatmap(scale(mat_to_plt_heatmap),
    name = "deltaTopics",
    split = merged.meta.hat$status,
    row_names_gp = gpar(fontsize = 5),
    k = 4)
dev.off()

# fit coxph model
res_cox_summary = data.frame()
for(i in 1:32) {
    f <- as.formula(paste0("Surv(time, status) ~ topic", i))
    res.cox <- coxph(f, data =  merged.meta.hat)
    print(paste("Topic", i, summary(res.cox)$coefficients[1,5]))
    res_cox_summary <- rbind(res_cox_summary, data.frame(summary(res.cox)$coefficients))
}
colnames(res_cox_summary)[5] <- "wald_p"
res_cox_summary$topic <- rownames(res_cox_summary)
res_cox_summary <- res_cox_summary %>% data.table()
res_cox_summary$location <- location
res_cox_summary[wald_p < 0.1]$topic

file = paste0("survival/", "coxph_res_BETM_", location, "_", target, ".txt")
if(!exists(file)){
    fwrite(res_cox_summary, file = file)
}

# fit survival model
for(topic in paste0("topic", 1:32)) {
    # set Z-scale cut-offs for high and low and (mid) expression
    highExpr <- 1.0
    lowExpr <- -1.0
    
    merged.meta.hat[[paste0(topic,"_3group")]] <- ifelse(scale(merged.meta.hat[[topic]]) >= highExpr, 'High',
        ifelse(scale(merged.meta.hat[[topic]]) <= lowExpr, 'Low', 'Mid'))

    merged.meta.hat[[paste0(topic,"_2group")]] <- ifelse(scale(merged.meta.hat[[topic]]) >= 0, 'UP','DOWN')

    # fit survival model
    f_2group <- as.formula(paste0("Surv(time, status) ~ ", topic, "_2group"))
    fit_2group <- surv_fit(f_2group, data = merged.meta.hat)
   
    # Drawing survival curves
    #ggsurvplot(fit)
    plt_2group <- ggsurvplot(fit_2group, 
            conf.int=TRUE, 
            pval=TRUE, risk.table=TRUE, 
            legend.labs=c("DOWN", "UP"), 
            legend.title=topic,  
            title="Kaplan-Meier Curve for PDAC Survival", 
            risk.table.height=.25)
    pdf(paste0("survival/figures/KM/", "KM_BETM_", topic,"_2group_",location,"_", target,".pdf"))
    print(plt_2group, newpage = FALSE)
    dev.off()       

    # fit survival model
    f_3group <- as.formula(paste0("Surv(time, status) ~ ", topic, "_3group"))
    fit_3group <- surv_fit(f_3group, data = merged.meta.hat)
   
    # Drawing survival curves
    #ggsurvplot(fit)
    plt_3group <- ggsurvplot(fit_3group, 
            conf.int=TRUE, 
            pval=TRUE, risk.table=TRUE, 
            #legend.labs=c("High", "Low", "Mid"), 
            legend.title=topic,  
            title="Kaplan-Meier Curve for PDAC Survival", 
            risk.table.height=.25)
    pdf(paste0("survival/figures/KM/", "KM_BETM_", topic,"_3group_",location,"_", target,".pdf"))
    #.ggsave(paste0("survival/figures/KM/", "KM_", topic,"_US.pdf"), plot=plt,
    #    width=8, height=6)
    print(plt_3group, newpage = FALSE)
    dev.off()        
}
