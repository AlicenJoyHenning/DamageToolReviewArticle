
## Install scDesign2
#library(devtools)
#devtools::install_github("JSB-UCLA/scDesign2")
#install.packages("copula")


# Load all packages 
packages <- c("scDesign2", "copula", "plyr", "reshape2", "gridExtra", 
              "Seurat", "SeuratData", "ggpubr", "cowplot", "ggplot2", "scRNAseq")

for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    library(pkg)
  }
}

# SCE available data ----
# Read single cell data for testing 
help(package = "scRNAseq")
view_data <- as.data.frame(surveyDatasets())
View(view_data)


# Sample Prep 
# PBMC multicellular
pbmc <- fetchDataset("kotliarov-pbmc-2020", "2024-04-18")
metadata <- colData(pbmc)
metadata$phenotype_sample <- paste0(metadata$adjmfc.time, "_", metadata$sample)
counts <- pbmc@assays@data$counts
rownames(counts) <- rownames(pbmc)
colnames(counts) <- colnames(pbmc)
counts  <- as(counts, "sparseMatrix")
high_sample <- subset(metadata, phenotype_sample == "d0 high_207_d0")
high_sample <- rownames(high_sample)
high_counts  <- counts[, high_sample]
low_sample <- subset(metadata, phenotype_sample == "d0 low_277_d0")
low_sample <- rownames(low_sample)
low_counts  <- counts[, low_sample]

# One cell type ("unicellular") 
tcell <- fetchDataset("bacher-tcell-2020", "2023-12-21")
counts <- tcell@assays@data$counts
rownames(counts) <- rownames(tcell)
colnames(counts) <- colnames(tcell)
metadata <- data.frame(colData(tcell))
covid <- paste0(metadata$barcode[metadata$diagnosis == "COVID19"], "-", metadata$sample[metadata$diagnosis == "COVID19"])
healthy <- paste0(metadata$barcode[metadata$diagnosis == "Healthy"], "-", metadata$sample[metadata$diagnosis == "Healthy"])
covid19 <- counts[, covid]
healthy <- counts[, healthy]
covid19  <- as(covid19 , "sparseMatrix")
selected_covid19  <- sample(colnames(covid19 ), 10000)
covid19  <- covid19[, selected_covid19]
healthy  <- as(healthy, "sparseMatrix")
selected_healthy  <- sample(colnames(healthy), 10000)
healthy  <- healthy[, selected_healthy]


# TUMOUR 
cancer <- fetchDataset("zilionis-lung-2019", "2023-12-20", "human")

# Filter 
mito_genes <- grepl("^MT-", rownames(cancer))
percent_mito <- ((colSums(counts(cancer)[mito_genes, , drop = FALSE])) / (colSums(counts(cancer)))) * 100
keep_cells <- which(percent_mito <= 15)
cancer <- cancer[, keep_cells]
counts <- cancer@assays@data$counts
rownames(counts) <- rownames(cancer)
colnames(counts) <- colData(cancer)[, "Barcode"]
metadata <- data.frame(colData(cancer))
table(metadata$Patient)
# p1    p2    p3    p4    p5    p6    p7 
# 24497 17071 20695 10694  8554  9408 14498 

p1 <- metadata$Barcode[metadata$Patient == "p5"]
p2 <- metadata$Barcode[metadata$Patient == "p7"]

# divide into two categories 
p1_counts <- counts[, p1]
p2_counts <- counts[, p2]

p1_counts  <- as(p1_counts, "sparseMatrix")
selected_p1_counts  <- sample(colnames(p1_counts), 5000)
p1_counts  <- p1_counts[, selected_p1_counts]

p2_counts  <- as(p2_counts, "sparseMatrix")
selected_p2_counts  <- sample(colnames(p2_counts), 5000)
p2_counts  <- p2_counts[, selected_p2_counts]


p1_results <- simulate_damage(p1_counts) 
p2_results <- simulate_damage(p2_counts)

results <- p1_results
results <- p2_results


# Plot to see
mito_ribo_old <- ggplot(results$qc_summary, aes(x = Original_RiboProp, y = Original_MitoProp, colour = Damaged_Status)) + 
  scale_color_manual(values = c("damaged" = "#FE3A55", "cell" = "#C1C1C2")) + 
  scale_y_continuous(labels = number_format(accuracy = 0.1)) + 
  geom_point() + 
  theme_classic() + 
  scale_x_continuous(limits = c(0, 1), labels = number_format(accuracy = 0.1)) +
  scale_y_continuous(limits = c(0, 1), labels = number_format(accuracy = 0.1)) +
  theme(plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(fill = NA, color = "black"),
        legend.position = "none")

mito_ribo_new <- ggplot(results$qc_summary, aes(x = New_RiboProp, y = New_MitoProp, colour = Damaged_Status)) + 
  scale_color_manual(values = c("damaged" = "#FE3A55", "cell" = "#C1C1C2")) + 
  scale_y_continuous(labels = number_format(accuracy = 0.1)) + 
  geom_point() + 
  scale_x_continuous(limits = c(0, 1), labels = number_format(accuracy = 0.1)) +
  scale_y_continuous(limits = c(0, 1), labels = number_format(accuracy = 0.1)) +
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(fill = NA, color = "black"),
        legend.position = "none")

mito_features_old <- ggplot(results$qc_summary, aes(x = Original_Features, y = Original_MitoProp, colour = Damaged_Status)) + 
  scale_color_manual(values = c("damaged" = "#FE3A55", "cell" = "#C1C1C2")) + 
  scale_y_continuous(labels = number_format(accuracy = 0.1)) + 
  geom_point() + 
  scale_y_continuous(limits = c(0, 1), labels = number_format(accuracy = 0.1)) +
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(fill = NA, color = "black"),
        legend.position = "none")

mito_features_new <- ggplot(results$qc_summary, aes(x = New_Features, y = New_MitoProp, colour = Damaged_Status)) + 
  scale_color_manual(values = c("damaged" = "#FE3A55", "cell" = "#C1C1C2")) + 
  scale_y_continuous(labels = number_format(accuracy = 0.1)) + 
  geom_point() + 
  scale_y_continuous(limits = c(0, 1), labels = number_format(accuracy = 0.1)) +
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(fill = NA, color = "black"),
        legend.position = "none")

((mito_features_old / mito_features_new) | (mito_ribo_old / mito_ribo_new))

# Seurat from output 
seurat <- CreateSeuratObject(counts = results$new_matrix, assay = "RNA") # will convert to dgCMatrix
seurat <- NormalizeData(seurat) %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(dims = 1:30) %>%
  FindClusters() %>%
  RunUMAP(dims = 1:30)
seurat$orig.ident <- results$qc_summary$Damaged_Status[match(rownames(seurat@meta.data), results$qc_summary$Cell)]


umap <- DimPlot(seurat, group.by = "orig.ident", cols = c("damaged" = "#FE3A55", "cell" = "#C1C1C2")) + NoAxes() + theme(panel.background = element_rect(fill = NULL, color = "black"))
umap | ((mito_features_old / mito_features_new) | (mito_ribo_old / mito_ribo_new))





# My processed -----

test <- readRDS("Projects/ReviewArticle/Zenodo/C_Test_Strategies/data/preprocess_output/HCT116_processed.rds")
test <- subset(test, mt.percent <= 20) # Ensure clean 
test <- as.sparse(test@assays$RNA$counts) # 16876 features 4338 cells
dim(test)
selected <- sample(colnames(test), 2500)
test <- test[, selected]
dim(test) # 16878 2500

mito_idx <- grep("^MT-", rownames(test), ignore.case = FALSE) # YES expression
test[mito_idx, 1]
results <- simulate_damage(count_matrix = test, organism = "Hsap") 


# Plot to see
mito_ribo_old <- ggplot(results$qc_summary, aes(x = Original_RiboProp, y = Original_MitoProp, colour = Damaged_Status)) + 
  scale_color_manual(values = c("damaged" = "#FE3A55", "cell" = "#C1C1C2")) + 
  scale_y_continuous(labels = number_format(accuracy = 0.1)) + 
  geom_point() + 
  theme_classic() + 
  scale_x_continuous(limits = c(0, 1), labels = number_format(accuracy = 0.1)) +
  scale_y_continuous(limits = c(0, 1), labels = number_format(accuracy = 0.1)) +
  theme(plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(fill = NA, color = "black"),
        legend.position = "none")

mito_ribo_new <- ggplot(results$qc_summary, aes(x = New_RiboProp, y = New_MitoProp, colour = Damaged_Status)) + 
  scale_color_manual(values = c("damaged" = "#FE3A55", "cell" = "#C1C1C2")) + 
  scale_y_continuous(labels = number_format(accuracy = 0.1)) + 
  geom_point() + 
  scale_x_continuous(limits = c(0, 1), labels = number_format(accuracy = 0.1)) +
  scale_y_continuous(limits = c(0, 1), labels = number_format(accuracy = 0.1)) +
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(fill = NA, color = "black"),
        legend.position = "none")

mito_features_old <- ggplot(results$qc_summary, aes(x = Original_Features, y = Original_MitoProp, colour = Damaged_Status)) + 
  scale_color_manual(values = c("damaged" = "#FE3A55", "cell" = "#C1C1C2")) + 
  scale_y_continuous(labels = number_format(accuracy = 0.1)) + 
  geom_point() + 
  scale_y_continuous(limits = c(0, 1), labels = number_format(accuracy = 0.1)) +
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(fill = NA, color = "black"),
        legend.position = "none")

mito_features_new <- ggplot(results$qc_summary, aes(x = New_Features, y = New_MitoProp, colour = Damaged_Status)) + 
  scale_color_manual(values = c("damaged" = "#FE3A55", "cell" = "#C1C1C2")) + 
  scale_y_continuous(labels = number_format(accuracy = 0.1)) + 
  geom_point() + 
  scale_y_continuous(limits = c(0, 1), labels = number_format(accuracy = 0.1)) +
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(fill = NA, color = "black"),
        legend.position = "none")

((mito_features_old / mito_features_new) | (mito_ribo_old / mito_ribo_new))

# Seurat from output 
seurat <- CreateSeuratObject(counts = results$new_matrix, assay = "RNA") # will convert to dgCMatrix
seurat <- NormalizeData(seurat) %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(dims = 1:30) %>%
  FindClusters() %>%
  RunUMAP(dims = 1:30)
seurat$orig.ident <- results$qc_summary$Damaged_Status[match(rownames(seurat@meta.data), results$qc_summary$Cell)]


umap <- DimPlot(seurat, group.by = "orig.ident", cols = c("damaged" = "#FE3A55", "cell" = "#C1C1C2")) + NoAxes() + theme(panel.background = element_rect(fill = NULL, color = "black"))
umap | ((mito_features_old / mito_features_new) | (mito_ribo_old / mito_ribo_new))

# SeuratData publicly available data ----
available_data <- AvailableData()


# IFNB-Stimulated & Control PBMCs
options(timeout = 600) # allow longer than 60 second to load
InstallData("pbmcsca")
data("pbmcsca")
pbmcsca <- UpdateSeuratObject(pbmcsca)

# Recollect celltypes  (don't need fine annotations, only coarse)
pbmcsca$mt.percent <- PercentageFeatureSet(pbmcsca, "^MT-")
pbmcsca <- subset(pbmcsca, mt.percent <= 25) # Ensure all definitely high quality (may lose some high quality but not concern here)
pbmcsca <- subset(pbmcsca, Method %in% c("10x Chromium (v3)", "10x Chromium (v2)"))
pbmcsca <- subset(pbmcsca, CellType != "Unassigned")

# Isolate stimulated and control cells (used to simulate datasets in isolation- mimic reality)
sample_10x_v3 <- subset(pbmcsca, Method == "10x Chromium (v3)")  
sample_10x_v3 <- subset(sample_10x, nFeature_RNA >= 700)
sample_10x_v2 <- subset(pbmcsca, Method == "10x Chromium (v2)") 

# Select random subset of each 
selected_sample_10x_v3 <- sample(colnames(sample_10x), 2500)
selected_sample_10x_v2 <- sample(colnames(sample_10x_v2), 2500)

sample_10x  <- as.sparse(sample_10x@assays$RNA$counts) 
sample_10x  <- sample_10x[, selected_sample_10x]
dim(sample_10x ) # 33694  2500

sample_10x_v2  <- as.sparse(sample_10x_v2@assays$RNA$counts) 
sample_10x_v2  <- sample_10x_v2[, selected_sample_10x_v2]
dim(sample_seqwell) # 33694  2500

# Run damage perturbation
sample_10x_perturbed <- simulate_damage(count_matrix = sample_10x)
sample_10x_v2_peturbed <- simulate_damage(count_matrix = sample_10x_v2)

test <- sample_10x_perturbed
test <- sample_10x_v2_peturbed
  
# Plot to see
mito_ribo_old <- ggplot(test$qc_summary, aes(x = Original_RiboProp, y = Original_MitoProp, colour = Damaged_Status)) + 
  scale_color_manual(values = c("damaged" = "#FE3A55", "cell" = "#C1C1C2")) + 
  scale_y_continuous(labels = number_format(accuracy = 0.1)) + 
  geom_point() + 
  theme_classic() + 
  scale_x_continuous(limits = c(0, 1), labels = number_format(accuracy = 0.1)) +
  scale_y_continuous(limits = c(0, 1), labels = number_format(accuracy = 0.1)) +
  theme(plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(fill = NA, color = "black"),
        legend.position = "none")

mito_ribo_new <- ggplot(test$qc_summary, aes(x = New_RiboProp, y = New_MitoProp, colour = Damaged_Status)) + 
  scale_color_manual(values = c("damaged" = "#FE3A55", "cell" = "#C1C1C2")) + 
  scale_y_continuous(labels = number_format(accuracy = 0.1)) + 
  geom_point() + 
  scale_x_continuous(limits = c(0, 1), labels = number_format(accuracy = 0.1)) +
  scale_y_continuous(limits = c(0, 1), labels = number_format(accuracy = 0.1)) +
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(fill = NA, color = "black"),
        legend.position = "none")

mito_features_old <- ggplot(test$qc_summary, aes(x = Original_Features, y = Original_MitoProp, colour = Damaged_Status)) + 
  scale_color_manual(values = c("damaged" = "#FE3A55", "cell" = "#C1C1C2")) + 
  scale_y_continuous(labels = number_format(accuracy = 0.1)) + 
  geom_point() + 
  scale_y_continuous(limits = c(0, 1), labels = number_format(accuracy = 0.1)) +
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(fill = NA, color = "black"),
        legend.position = "none")

mito_features_new <- ggplot(test$qc_summary, aes(x = New_Features, y = New_MitoProp, colour = Damaged_Status)) + 
  scale_color_manual(values = c("damaged" = "#FE3A55", "cell" = "#C1C1C2")) + 
  scale_y_continuous(labels = number_format(accuracy = 0.1)) + 
  geom_point() + 
  scale_y_continuous(limits = c(0, 1), labels = number_format(accuracy = 0.1)) +
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(fill = NA, color = "black"),
        legend.position = "none")

((mito_features_old / mito_features_new) | (mito_ribo_old / mito_ribo_new))

# Seurat from output 
seurat <- CreateSeuratObject(counts = test$new_matrix, assay = "RNA") # will convert to dgCMatrix
seurat <- NormalizeData(seurat) %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(dims = 1:30) %>%
  FindClusters() %>%
  RunUMAP(dims = 1:30)
seurat$orig.ident <- test$qc_summary$Damaged_Status[match(rownames(seurat@meta.data), test$qc_summary$Cell)]


umap <- DimPlot(seurat, group.by = "orig.ident", cols = c("damaged" = "#FE3A55", "cell" = "#C1C1C2")) + NoAxes() + theme(panel.background = element_rect(fill = NULL, color = "black"))
umap | ((mito_features_old / mito_features_new) | (mito_ribo_old / mito_ribo_new))


# IFNB-Stimulated & Control PBMCs
options(timeout = 600) # allow longer than 60 second to load
InstallData("hcabm40k")
data("hcabm40k")
hcabm40k <- UpdateSeuratObject(hcabm40k)

# Recollect celltypes  (don't need fine annotations, only coarse)
hcabm40k$mt.percent <- PercentageFeatureSet(hcabm40k, "^MT-")
hcabm40k <- subset(hcabm40k, mt.percent <= 25) # Ensure all definitely high quality (may lose some high quality but not concern here)
selected_hcabm40k <- sample(colnames(hcabm40k), 2500)

hcabm40k  <- as.sparse(hcabm40k@assays$RNA$counts) 
hcabm40k  <- hcabm40k[, selected_hcabm40k]
dim(hcabm40k) 

# Run damage perturbation
hcabm40k_perturbed <- simulate_damage(count_matrix = hcabm40k)

test <- hcabm40k_perturbed

# Plot to see
mito_ribo_old <- ggplot(test$qc_summary, aes(x = Original_RiboProp, y = Original_MitoProp, colour = Damaged_Status)) + 
  scale_color_manual(values = c("damaged" = "#FE3A55", "cell" = "#C1C1C2")) + 
  scale_y_continuous(labels = number_format(accuracy = 0.1)) + 
  geom_point() + 
  theme_classic() + 
  scale_x_continuous(limits = c(0, 1), labels = number_format(accuracy = 0.1)) +
  scale_y_continuous(limits = c(0, 1), labels = number_format(accuracy = 0.1)) +
  theme(plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(fill = NA, color = "black"),
        legend.position = "none")

mito_ribo_new <- ggplot(test$qc_summary, aes(x = New_RiboProp, y = New_MitoProp, colour = Damaged_Status)) + 
  scale_color_manual(values = c("damaged" = "#FE3A55", "cell" = "#C1C1C2")) + 
  scale_y_continuous(labels = number_format(accuracy = 0.1)) + 
  geom_point() + 
  scale_x_continuous(limits = c(0, 1), labels = number_format(accuracy = 0.1)) +
  scale_y_continuous(limits = c(0, 1), labels = number_format(accuracy = 0.1)) +
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(fill = NA, color = "black"),
        legend.position = "none")

mito_features_old <- ggplot(test$qc_summary, aes(x = Original_Features, y = Original_MitoProp, colour = Damaged_Status)) + 
  scale_color_manual(values = c("damaged" = "#FE3A55", "cell" = "#C1C1C2")) + 
  scale_y_continuous(labels = number_format(accuracy = 0.1)) + 
  geom_point() + 
  scale_y_continuous(limits = c(0, 1), labels = number_format(accuracy = 0.1)) +
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(fill = NA, color = "black"),
        legend.position = "none")

mito_features_new <- ggplot(test$qc_summary, aes(x = New_Features, y = New_MitoProp, colour = Damaged_Status)) + 
  scale_color_manual(values = c("damaged" = "#FE3A55", "cell" = "#C1C1C2")) + 
  scale_y_continuous(labels = number_format(accuracy = 0.1)) + 
  geom_point() + 
  scale_y_continuous(limits = c(0, 1), labels = number_format(accuracy = 0.1)) +
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(fill = NA, color = "black"),
        legend.position = "none")

((mito_features_old / mito_features_new) | (mito_ribo_old / mito_ribo_new))

# Seurat from output 
seurat <- CreateSeuratObject(counts = test$new_matrix, assay = "RNA") # will convert to dgCMatrix
seurat <- NormalizeData(seurat) %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(dims = 1:30) %>%
  FindClusters() %>%
  RunUMAP(dims = 1:30)
seurat$orig.ident <- test$qc_summary$Damaged_Status[match(rownames(seurat@meta.data), test$qc_summary$Cell)]


umap <- DimPlot(seurat, group.by = "orig.ident", cols = c("damaged" = "#FE3A55", "cell" = "#C1C1C2")) + NoAxes() + theme(panel.background = element_rect(fill = NULL, color = "black"))
umap | ((mito_features_old / mito_features_new) | (mito_ribo_old / mito_ribo_new))





