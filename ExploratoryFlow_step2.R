#This script transforms the flowset in a singlecellexperiment object adding metadata
#Then the script performs unsupervised clustering by using flowSOM, dimensionality reduction (UMAP)
#At the end of the script by using diffcyt it performs differential abundance analysis

#Here are the workflows:
#http://127.0.0.1:26337/library/CATALYST/doc/differential.html
#http://127.0.0.1:26337/library/diffcyt/doc/diffcyt_workflow.html

rm(list = ls())

library(rstudioapi)
library(flowCore)
#library(cytofCore)
#library(cytofkit2)
#library(FlowSOM)
#library(cluster)
#library(Rtsne)
library(ggplot2)
library(dplyr)
#library(ggthemes)
#library(RColorBrewer)
#library(uwot)
library(CATALYST)
library(diffcyt)
library(stringr)
#library(scran)
#library(scater)
#library(ggcyto)
library(SingleCellExperiment)
#library(flowWorkspace)
#library(reshape2)
#library(ggrepel)
#library(slingshot)
#library(knn.covertree)
#library(destiny)
#library(readxl)
#library(scDataviz)
#library(ggpubr)
browseVignettes("diffcyt")
# Set PrimaryDirectory where this script is located
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()
PrimaryDirectory <- getwd()
PrimaryDirectory

# Define directory where transformed FCS files are located
fcs_t <- "fcs_t"
FCSDirectory <- paste(PrimaryDirectory, fcs_t, sep = "/")

# List FCS files and create flowset
FCSfiles <- list.files(path = FCSDirectory, pattern = ".fcs", full.names= FALSE)
fs <- read.flowSet(files = FCSfiles, path = FCSDirectory, transformation = FALSE, truncate_max_range = FALSE)

# Keyword ($CYT = "FACS"): CATALYST is a mass spectrometry package. 
# We need to change the keyword in order to adapt it to flow cytometry data
names(keyword(fs[[1]]))[15] <- "$CYT"
ds <- keyword(fs[[1]])
l <- list(cyt = "\\$CYT$")
keep <- lapply(l, grep, names(ds))
ds[[keep$cyt]] <- "FACS"
keyword(fs[[1]])[[keep$cyt]] <- "FACS"

##Building panel dataframe
# Define channels of interest and marker_class (look at the workflow)
fcs_colname <- colnames(fs)[c(7:9, 11:18)]
marker_class <- c(rep("type", 8), rep("state",2), "type")
antigen <- fcs_colname
length(marker_class) == length(fcs_colname)

#Panel
panel <- cbind(fcs_colname, antigen, marker_class)
panel <- as.data.frame(panel)
all(panel$fcs_colname %in% colnames(fs))

##Building metadata dataframe
condition <- FCSfiles
condition <- word(condition, 2,4, sep = "_")
condition[grepl("HC", condition)] <- "HC"
condition <- gsub("DG", "base", condition)

patient_id <- FCSfiles
patient_id <- word(patient_id, 1, sep = "_")

sample_id <- paste(patient_id, condition, sep = "_")

file_name <- FCSfiles

md <- cbind(file_name, sample_id, condition, patient_id)
md <- data.frame(md)

##SCE object
sce <- CATALYST::prepData(fs, panel = panel, md = md, transform = FALSE)
sce@assays@data$exprs <- sce@assays@data$counts

## QC
# Density
p <- plotExprs(sce, color_by = "condition")
p$facet$params$ncol <- 4
p

# Counts
n_cells(sce)
plotCounts(sce, group_by = "sample_id", color_by = "condition")

# Multidimensional scaling (MDS)
CATALYST::pbMDS(sce, color_by = "condition", label = NULL)

# Non Redundancy Score (NRS)
plotNRS(sce, features = type_markers(sce), color_by = "condition")

## Clustering
#FlowSOM
set.seed(1234)
sce <- cluster(sce,
               xdim = 10, ydim = 10, maxK = 20, seed = 1234)

delta_area(sce)
plotExprHeatmap(sce, features = "type",
                by = "cluster_id", k = "meta7", 
                bars = TRUE, perc = TRUE, col_anno = TRUE, scale = "last")

#Filter out subset <1%
sce <- filterSCE(sce, cluster_id %in% c(1:4,6), k = "meta7")

#Set parameters for UMAP
set.seed(1234)
n_cells <- 1000
exaggeration_factor <- 12.0
eta <- n_cells/exaggeration_factor
sce <- runDR(sce, dr =  "UMAP", cells = n_cells, features = "type")

##Add annotations
#Read annotation file:3 main regions Resting, Activated, Terminal differentiated
annotation_table <- readxl::read_excel("annotation_1.xlsx")

# convert to factor with merged clusters in desired order
annotation_table$new_cluster <- factor(annotation_table$new_cluster,
                                       levels = c("senescent-like CD8+", "PD1+CD28+ CD8+", "activated CD8+",
                                                  "naive CD8+"))
# apply manual annotation
sce <- mergeClusters(sce, k = "meta7", 
                     table = annotation_table, id = "cluster_annotation", overwrite = TRUE)
sce$condition <- ifelse(sce$condition == "CR_BM_base", "Res_bas", 
                        ifelse(sce$condition == "NR_BM_base", "NonRes_bas", 
                               ifelse(sce$condition == "CR_BM_post", "Res_post",
                                      ifelse(sce$condition == "NR_BM_post", "NonRes_post", "HD"))))

#Heatmap with annotations
plotExprHeatmap(sce, features = "type", 
                by = "cluster_id", k = "meta7", m = "cluster_annotation", scale = "last")
tiff("./plots/heat.tiff", width = 5*900, height = 5*300, res = 300, pointsize = 5)     
p <- plotExprHeatmap(sce, features = "type", 
                     by = "cluster_id", k = "cluster_annotation", scale = "last")
p
dev.off()

#UMAP with annotations
tiff("./plots/umap.tiff", width = 5*900, height = 5*900, res = 300, pointsize = 5)     
plot <- plotDR(sce, "UMAP", color_by = "cluster_annotation") + facet_wrap("condition")
plot + guides(color = guide_legend(ncol = 1, override.aes = list(size = 10))) + geom_density2d(binwidth = 0.006, colour = "black") + 
  theme(strip.text = element_text(size=15)) +theme(legend.text=element_text(size=15), legend.title=element_text(size=15))
dev.off()

# Statistical Analysis (diffcyt) To be reviewed!!!!
ei <- sce@metadata$experiment_info

design <- createDesignMatrix(
  ei, cols_design = c("condition", "patient_id")
)

design <- design[,-c(26, 32)] # not a full rank matrix, drop column 19

#NR base vs CR base
colnames(design)[1:5]
contrast <- createContrast(c(0, 0, 0, 1, 0, rep(0, 25)))
out_DA1 <- diffcyt(sce, design = design, contrast = contrast, method_DA = "diffcyt-DA-voom",
                   clustering_to_use = "cluster_annotation", verbose = TRUE, subsampling = TRUE, transform = FALSE, min_samples = 6)
topTable(out_DA1, format_vals = TRUE)
da_1 <- rowData(out_DA1$res)

# CR post NR post
colnames(design)[1:5]
contrast <- createContrast(c(0, -1, 0, 0, 1, rep(0, 25)))
out_DA2 <- diffcyt(sce, design = design, contrast = contrast, method_DA = "diffcyt-DA-voom",
                   clustering_to_use = "cluster_annotation", verbose = TRUE, subsampling = TRUE, transform = FALSE, min_samples = 6)
topTable(out_DA2, format_vals = TRUE)
da_2 <- rowData(out_DA2$res)

#NR base vs HC
colnames(design)[1:5]
contrast <- createContrast(c(0, 0, -1, 1, 0, rep(0, 25)))
out_DA3 <- diffcyt(sce, design = design, contrast = contrast, method_DA = "diffcyt-DA-voom",
                   clustering_to_use = "cluster_annotation", verbose = TRUE, subsampling = TRUE, transform = FALSE, min_samples = 6)
topTable(out_DA3, format_vals = TRUE)
da_3 <- rowData(out_DA3$res)

#NR post vs HC
colnames(design)[1:5]
contrast <- createContrast(c(0, 0, -1, 0, 1, rep(0, 25)))
out_DA4 <- diffcyt(sce, design = design, contrast = contrast, method_DA = "diffcyt-DA-voom",
                   clustering_to_use = "cluster_annotation", verbose = TRUE, subsampling = TRUE, transform = FALSE, min_samples = 6)
topTable(out_DA4, format_vals = TRUE)
da_4 <- rowData(out_DA4$res)

#Add bars to boxplot
#base CR vs base NR
tiff("./plots/box.tiff", width = 5*1000, height = 5*300, res = 300, pointsize = 5)     
bxp <- plotAbundances(sce, k = "cluster_annotation", by = "cluster_id")
bxp  + theme(strip.text = element_text(size=25)) +theme(legend.text=element_text(size=15), legend.title=element_text(size=15))

dev.off()


stat.test_1 <- as_tibble(da_1)
stat.test_1$p_adj
p.adj.signif <- c(rep("**",3))
y.position <- c(65, 58, 82)
group1 <- (rep("CR_BM_base",3))
group2 <- (rep("NR_BM_base", 3))
stat.test_1 <- cbind(stat.test_1, group1, group2, p.adj.signif, y.position)

bxp_1 <- bxp + stat_pvalue_manual(stat.test_1,  hide.ns = TRUE, label = "p.adj.signif", tip.length = 0, size = 2.5) 
bxp_1

stat.test_2 <- as_tibble(da_2)
stat.test_2$p_adj
p.adj.signif <- c("***", "ns", "***")
y.position <- c(62.5, 40, 79.8)
group1 <- (rep("CR_BM_post",3))
group2 <- (rep("NR_BM_post", 3))
stat.test_2 <- cbind(stat.test_2, group1, group2, p.adj.signif, y.position)

bxp_2 <- bxp_1 + stat_pvalue_manual(stat.test_2, hide.ns = TRUE, label = "p.adj.signif", tip.length = 0, size = 2.5)

bxp_2

stat.test_3 <- as_tibble(da_3)
stat.test_3$p_adj
p.adj.signif <- c("**", "ns", "***")
y.position <- c(60, 39, 77)
group1 <- (rep("NR_BM_base",3))
group2 <- (rep("HC", 3))
stat.test_3 <- cbind(stat.test_3, group1, group2, p.adj.signif, y.position)

bxp_3 <- bxp_2 + stat_pvalue_manual(stat.test_3, hide.ns = TRUE, label = "p.adj.signif", tip.length = 0, size = 2.5)

bxp_3

stat.test_4 <- as_tibble(da_4)
stat.test_4$p_adj
p.adj.signif <- c("*", "ns", "**")
y.position <- c(67, 45, 85)
group1 <- (rep("NR_BM_post",3))
group2 <- (rep("HC", 3))
stat.test_4 <- cbind(stat.test_4, group1, group2, p.adj.signif, y.position)

bxp_4 <- bxp_3 + stat_pvalue_manual(stat.test_4, hide.ns = TRUE, label = "p.adj.signif", tip.length = 0, size = 2.5)

bxp_4

# save workspace 
save(list = ls(), file = "explflow_step2.rds")
