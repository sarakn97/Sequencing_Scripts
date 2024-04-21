# scRNA Analysis Bone Marrow (BM)

# initialization
library(dplyr)
library(Seurat)
library(SeuratData)
library(patchwork)
library(ggplot2)
library(DoubletFinder)

# Load Seurat Object
setwd("/storage1/fs1/leyao.wang/Active/saran/scRNA_BM_sara/")
combined.all <- readRDS("./combined_all_ShanLab_BM.rds")

############################ Read Input Data #######################################
# Read in Cell Calls
BM_Infected <- read.csv("/storage1/fs1/leyao.wang/Active/saran/scRNA_BM_sara/BM_I/gem_classification.csv")
BM_Infected_mouse <- subset(BM_Infected, BM_Infected$call=="mm10")

BM_NOTinfected <- read.csv("/storage1/fs1/leyao.wang/Active/saran/scRNA_BM_sara/BM_NI/gem_classification.csv")
BM_NOTinfected_mouse <- subset(BM_NOTinfected, BM_NOTinfected$call=="mm10")

################## read in sample 1 (Infected)
BM_HIV = Read10X("/storage1/fs1/leyao.wang/Active/saran/scRNA_BM_sara/BM_I/filtered_feature_bc_matrix/")
dim(BM_HIV)
#c6888c6  #5356

# remove mouse cells (BM_HIV)  
BM_HIV <- BM_HIV[,!colnames(BM_HIV)%in%BM_Infected_mouse$barcode]
dim(BM_HIV)
#68886   #5318

# Create Sample1 Seurat Object (Infected)
infected_BM <- CreateSeuratObject(counts = BM_HIV, project = "BM_infected", min.cells = 3, min.features = 200)
dim(infected_BM)
#24286 5286
head(infected_BM)

################# read sample 2 (Control)
BM_NI = Read10X("/storage1/fs1/leyao.wang/Active/saran/scRNA_BM_sara/BM_NI/filtered_feature_bc_matrix/")
dim(BM_NI)
#68886  7553

# remove mouse cells (BM_NI) 
BM_NI <- BM_NI[,!colnames(BM_NI)%in%BM_NOTinfected_mouse$barcode]
dim(BM_NI)
#c6888c6 7529

# create sample2 suerat object
NotInf_BM <- CreateSeuratObject(counts = BM_NI, project = "NotInf_BM", min.cells = 3, min.features = 200)
dim(NotInf_BM)
# 29325  7409
head(NotInf_BM)

############################################# Filtering ###########################
# Identify Condition
infected_BM[["sample"]] <- "Infected"
NotInf_BM[["sample"]] <- "Not"

# Merge Samples
merged_seurat <- merge(x=infected_BM, y=NotInf_BM, add.cell.id = c("Infected", "Not"))
# View(merged_seurat@meta.data)

sum(grepl("MT",rownames(merged_seurat)))
# 254
rownames(merged_seurat)[grepl("MT-",rownames(merged_seurat))]

# calculate the proportion of transcripts mapping to mitochondrial genes
merged_seurat[["percent.mt"]] <- PercentageFeatureSet(merged_seurat, pattern = "MT-")

max(merged_seurat@meta.data$percent.mt)
#84.3038

# PLot Quality
VlnPlot(merged_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(object=merged_seurat, feature1="nCount_RNA", feature2= "percent.mt", cols=c("red", "blue"))
FeatureScatter(object=merged_seurat, feature1="nCount_RNA", feature2= "nFeature_RNA", cols=c("red", "blue"))

# Filter based on quality plots 
Filter_BM <- subset(x=merged_seurat, subset=nFeature_RNA > 200 & nFeature_RNA < 9000 & percent.mt > -Inf & percent.mt < 6)
# Reassess Plots
VlnPlot(Filter_BM, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)
FeatureScatter(object=Filter_BM, feature1="nCount_RNA", feature2= "percent.mt", cols=c("red", "blue"))
FeatureScatter(object=merged_seurat, feature1="nCount_RNA", feature2= "nFeature_RNA", cols=c("red", "blue"))

# Sums post Filtering
sum(Filter_BM$sample=="Infected")
# 5062
sum(Filter_BM$sample=="Not")
# 6842

# Fix Gene Names
Filter_BM@assays$RNA@counts@Dimnames[[1]] <- gsub("GRCh38-","",Filter_BM@assays$RNA@counts@Dimnames[[1]])
Filter_BM@assays$RNA@data@Dimnames[[1]] <- gsub("GRCh38-","",Filter_BM@assays$RNA@counts@Dimnames[[1]])

split_seurat <- SplitObject(Filter_BM, split.by = "sample")
split_seurat <- split_seurat[c("Infected", "Not")]

############################################################## DOUBLET FINDER #######
# Pre-Process Seurat Object
for (i in 1:length(split_seurat)) {
  split_seurat[[i]] <- SCTransform(split_seurat[[i]], vars.to.regress = c("percent.mt"))
  split_seurat[[i]] <- RunPCA(split_seurat[[i]])
  split_seurat[[i]] <- RunUMAP(split_seurat[[i]], dims=1:20)
}

# Doublet Finder Infected
for (i in 1:length(split_seurat)) {
  split_seurat[[i]] <- split_seurat[[i]] %>%
    FindNeighbors(dims = 1:20)%>%
    FindClusters(resolution=0.6)
}

###################### Function to Find Doublets
findDoublets <- function(seu, dfr, cluster) {
  # seu = seurat object
  # dfr = doublet formation rate, adjust for cell count
  # cluster = seurat metadata column containing clusters
  
  sweep.res.list_1 <- paramSweep_v3(seu, PCs= 1:10, sct=T)
  sweep.stats_1 <- summarizeSweep(sweep.res.list_1, GT=F)
  bcmvn_1 <- find.pK(sweep.stats_1)
  
  pK=as.numeric(as.character(bcmvn_1$pK))
  BCmetric=bcmvn_1$BCmetric
  pK_chosen = pK[which(BCmetric %in% max(BCmetric))]
  
  annotations <- seu@meta.data$cluster
  homotypic.prop <- modelHomotypic(annotations)
  nExp_poi <- round(dfr*nrow(seu@meta.data)) # Adjust according to Doublet Rate (4000 cells) 
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  
  # find doublets, will create metadata column specifying cells as doublet or singlet
  seu <- doubletFinder_v3(seu, PCs=1:10, pN = 0.25, pK = pK_chosen, nExp = nExp_poi, reuse.pANN = FALSE, sct=T)
}

# Run Doublet Finder
findDoublets(split_seurat[[1]], 0.039, "SCT_snn_res.0.6")
findDoublets(split_seurat[[2]], 0.0524, "SCT_snn_res.0.6")


# REMOVE DOUBLETS #
split_seurat[[1]] <- subset(x = split_seurat[[1]], 
                            subset = DF.classifications_0.25_0.15_197!="Doublet") 
dim(split_seurat[[1]])
#4865
split_seurat[[2]] <- subset(x = split_seurat[[2]], 
                            subset =DF.classifications_0.25_0.29_359 !="Doublet") 
dim(split_seurat[[2]])
#6483

########################################################## Downsample ###########
set.seed(111)
split_seurat[[1]] <- split_seurat[[1]][, sample(colnames(split_seurat[[1]]), size = 4850, replace=F)] # set to size of smallest sample
split_seurat[[2]] <- split_seurat[[2]][, sample(colnames(split_seurat[[2]]), size = 4850, replace=F)]

for (i in 1:length(split_seurat)){
  split_seurat[[i]] <- RunTSNE(split_seurat[[i]], dims = 1:20)
}

########################################################### Integration ##########
# select most variable features
integ_features <- SelectIntegrationFeatures(object.list=split_seurat, nfeatures=6000)

# prepare SCT list object for integration
split_seurat <- PrepSCTIntegration(object.list=split_seurat, anchor.features=integ_features)

# Find Anchors
integ_anchors <- FindIntegrationAnchors(object.list=split_seurat, normalization.method = "SCT",
                                        anchor.features=integ_features)

# Integrate
combined.all <- IntegrateData(anchorset=integ_anchors, normalization.method="SCT")

# Save integrated seurat object & workspace
#saveRDS(combined.all, "/storage1/fs1/leyao.wang/Active/scRNA_BM_sara/SCT_integrated_BM_seurat.rds")

########################## Pre-Processing ########################################

DefaultAssay(combined.all) <- "integrated"
combined.all <- combined.all %>% RunPCA(npcs = 30, verbose = FALSE) %>%
  RunUMAP(dims = 1:20) %>% FindNeighbors(reduction = "pca", dims = 1:20) %>% 
  FindClusters(resolution = c(0.6, 0.8)) %>% RunTSNE() %>% FindVariableFeatures(selection.method = "vst", nfeatures = 2000)

######## Save files
# write expression counts matrix
library(Matrix)
counts_matrix <- GetAssayData(combined.all, assay='RNA', slot='counts')
writeMM(counts_matrix, file='/storage1/fs1/leyao.wang/Active/scRNA_sara/scRNA_BM_sara/counts_BM_organoid.mtx')
rownames(counts_matrix)

########## Make tSNE and UMAP plots ###############################################
DefaultAssay(combined.all) <- "integrated"
DimPlot(combined.all, reduction = "pca", split.by = "sample",label = TRUE,repel = TRUE)

DimPlot(combined.all, reduction = "umap", label = TRUE,repel = TRUE)
DimPlot(combined.all, reduction = "umap", split.by = "sample",label = TRUE,repel = TRUE)
DimPlot(combined.all, reduction = "tsne", label = TRUE,repel = TRUE)

# Use RNA Assay for Visualization, scale and normalize raw counts
DefaultAssay(combined.all) <- "RNA"
combined.all <- NormalizeData(assay = "RNA") %>%
  ScaleData(assay = "RNA")

# Manually Annotate Clusters
combined_markers <- FindAllMarkers(object = combined.all, 
                                   only.pos = TRUE,
                                   logfc.threshold = 0.25) 

# create columns for annotations
combined.all@meta.data$basic <- "NA"
combined.all@meta.data$zoom <- "NA"

########################### Manual Annotation ####################################
#### Feature Plots
DefaultAssay(combined.all) <- "integrated"
# Feature Plots
DefaultAssay(combined.all) <- "RNA"
FeaturePlot(combined.all, features = c("CD34", "CD38", "AVP", "MPO"), min.cutoff = 0, max.cutoff = 3, slot = "counts", label.size = 2.5, label = F, repel=T, pt.size = 0.7)
FeaturePlot(combined.all, features = c("DNTT", "VPREB3","RAG1", "MKI67"), min.cutoff = 0, max.cutoff = 3, slot = "counts", label.size = 2.5, label = F, repel=T, pt.size = 0.7)
FeaturePlot(combined.all, features = c("IRF8", "BATF3", "LILRA4", "S100A8"), min.cutoff = 0, max.cutoff = 3, slot = "counts", label.size = 2.5, label = F, repel=T, pt.size = 0.7)
FeaturePlot(combined.all, features = c("ELANE", "CEBPD", "CYGB", "HVCN1"), min.cutoff = 0, max.cutoff = 3, slot = "counts", label.size = 2.5, label = F, repel=T, pt.size = 0.7)
FeaturePlot(combined.all, features = c("APOC1", "GATA1","GATA2", "HVCN1"), min.cutoff = 0, max.cutoff = 3, slot = "counts", label.size = 2.5, label = F, repel=T, pt.size = 0.7)
FeaturePlot(combined.all, features = c("VCAN", "ACY3", "FCER1A", "HOPX"), min.cutoff = 0, max.cutoff = 3, slot = "counts", label.size = 2.5, label = F, repel=T, pt.size = 0.7)

# Apply Annotations
gmp <- WhichCells(combined.all, idents = 8)
memp <- WhichCells(combined.all, idents = 10)
pre_mono <- WhichCells(combined.all, idents = 9)
dc <- WhichCells(combined.all, idents = c(6,11))
clp1 <- WhichCells(combined.all, idents = 13)
clp2 <- WhichCells(combined.all, idents = 0)
clp_cyc <- WhichCells(combined.all, idents = 5)
preproB <- WhichCells(combined.all, idents = c(7))
proB <- WhichCells(combined.all, idents = 4)
preB <- WhichCells(combined.all, idents = c(1,2))

combined.all@meta.data[gmp, "annotate"] <- "GMP"
combined.all@meta.data[memp, "annotate"] <- "MEMP"
combined.all@meta.data[pre_mono, "annotate"] <- "PRE-MONO"
combined.all@meta.data[dc, "annotate"] <- "DC"
combined.all@meta.data[clp1, "annotate"] <- "CLP.1"
combined.all@meta.data[clp2, "annotate"] <- "CLP.2"
combined.all@meta.data[clp_cyc, "annotate"] <- "CLP-cycling"
combined.all@meta.data[preproB, "annotate"] <- "Pre-Pro-B"
combined.all@meta.data[proB, "annotate"] <- "Pro-B"
combined.all@meta.data[preB, "annotate"] <- "Pre-B"

# Plot w/ Annotations
combined.all <- SetIdent(combined.all, value = "annotate")
DimPlot(combined.all, label = T)

# Function to process any clusters you decide to subcluster 
process_subcluster <- function(x) {
  DefaultAssay(x) <- "integrated"
  x <- x %>% RunPCA(npcs = 30, verbose = FALSE) %>%
    RunUMAP(dims = 1:20) %>% FindNeighbors(reduction = "pca", dims = 1:20) %>% 
    FindClusters(resolution = c(0.6, 0.8))
  
}

# Sub-Cluster C.3
hsc_mpp <- subset(combined.all, ident = 3)
hsc_mpp <- process_subcluster(hsc_mpp)
DimPlot(hsc_mpp)
FeaturePlot(hsc_mpp, features = c("CD38", "CD34", "PROM1", "MLLT3", "MEIS1", "MSI2", "SPINK2", "FLT3", "HOPX"))
hsc <- WhichCells(hsc_mpp, idents = 1)
mpp1 <- WhichCells(hsc_mpp, idents = c(2,4))
mpp2 <- WhichCells(hsc_mpp, idents = 3)
mpp3 <- WhichCells(hsc_mpp, idents = c(5,6))
combined.all@meta.data[hsc, "annotate"] <- "HSC"
combined.all@meta.data[mpp1, "annotate"] <- "MPP1"
combined.all@meta.data[mpp2, "annotate"] <- "MPP2"
combined.all@meta.data[mpp3, "annotate"] <- "MPP3"

###################################### Analysis #################################
####################### DESeq Differential Expression ###########################
table(Idents(combined.all),combined.all$sample)
library(edgeR)

deseq <- function(seu, x) {
  # seu = seurat object
  # x = cluster name
  sc = subset(seu,idents= x)
  Idents(sc) = sc$sample
  DefaultAssay(sc) = 'RNA'
  sc = NormalizeData(sc) %>% ScaleData()
  mtx = as.matrix(GetAssayData(sc,slot='counts'))
  filter = rowSums(mtx>0) > 10
  dge = DGEList(mtx[filter,],norm.factors=rep(1,length(mtx[1,])),group=sc$sample)
  sc$sample = factor(sc$sample,levels=c('Not','Infected'))
  design = model.matrix(~0 + sc$sample)
  colnames(design) = c('Not','Infected')
  contrast.matrix = makeContrasts(Infected-Not,levels=design)
  dge = calcNormFactors(dge,method='TMM')
  dge = estimateDisp(dge,design=design)
  
  fit = glmQLFit(dge,design=design)
  res = glmQLFTest(fit,contrast=contrast.matrix)
  pval = res$table
  pval$fdr = p.adjust(pval$PValue,method='fdr')
  pval = arrange(pval,fdr)
  
  write.csv(pval, paste('/storage1/fs1/leyao.wang/Active/saran/scRNA_BM_sara/DEGs', x, 'csv', sep = '.'))
}

# Run DESeq Differential Expression for All Clusters
for(i in list(unique(combined.all$main))) {
  deseq(combined.all, i)
}

########################## Seurat Differential Expression #########################
# Default : Wilcoxon Rank Sum Test
combined.all$celltype.stim <- paste(Idents(combined.all), combined.all$sample, sep = "_")
Idents(combined.all) <- "celltype.stim"
table(combined.all@meta.data$celltype.stim)

proGran_DEgenes <- FindMarkers(combined.all, ident.1 = "pro-Granulocytes_Infected", ident.2 = "pro-Granulocytes_Not", verbose = FALSE)

myeDC_DEgenes <- FindMarkers(combined.all, ident.1 = "Mye-DC_Infected", ident.2 = "Mye-DC_Not", verbose = FALSE)

premono_DEgenes <- FindMarkers(combined.all, ident.1 = "pre-Monocytes_Infected", ident.2 = "pre-Monocytes_Not", verbose = FALSE)

mono_DEgenes <- FindMarkers(combined.all, ident.1 = "Monocytes_Infected", ident.2 = "Monocytes_Not", verbose = FALSE)

predcs_DEgenes <- FindMarkers(combined.all, ident.1 = "pre-DC Cycle_Infected", ident.2 = "pre-DC Cycle_Not", verbose = FALSE)

predcsCYC_DEgenes <- FindMarkers(combined.all, ident.1 = "pre-DCs_Infected", ident.2 = "pre-DCs_Not", verbose = FALSE)

MEMP_DEgenes <- FindMarkers(combined.all, ident.1 = "MEMP_Infected", ident.2 = "MEMP_Not", verbose = FALSE)

seven_DEgenes <- FindMarkers(combined.all, ident.1 = "7_Infected", ident.2 = "7_Not", verbose = FALSE)

mpp1_DEgenes <- FindMarkers(combined.all, ident.1 = "MPP1_Infected", ident.2 = "MPP1_Not", verbose = FALSE)

mpp2_DEgenes <- FindMarkers(combined.all, ident.1 = "MPP2_Infected", ident.2 = "MPP2_Not", verbose = FALSE)

mpp3_DEgenes <- FindMarkers(combined.all, ident.1 = "MPP3_Infected", ident.2 = "MPP3_Mock", verbose = FALSE)

mpp4_DEgenes <- FindMarkers(combined.all, ident.1 = "MPP4_Infected", ident.2 = "MPP4_Mock", verbose = FALSE)

hsc_DEgenes <- FindMarkers(combined.all, ident.1 = "HSC_Infected", ident.2 = "HSC_Mock", verbose = FALSE)

cdp_DEgenes <- FindMarkers(combined.all, ident.1 = "CDP_Infected", ident.2 = "CDP_Mock", verbose = FALSE)

lmpp_DEgenes <- FindMarkers(combined.all, ident.1 = "LMPP_Infected", ident.2 = "LMPP_Mock", verbose = FALSE)

clp1_DEgenes <- FindMarkers(combined.all, ident.1 = "CLP.1_Infected", ident.2 = "CLP.1_Not", verbose = FALSE)

clp2_DEgenes <- FindMarkers(combined.all, ident.1 = "CLP.2_Infected", ident.2 = "CLP.2_Not", verbose = FALSE)

proB_DEgenes <- FindMarkers(combined.all, ident.1 = "Pro-B_Infected", ident.2 = "Pro-B_Not", verbose = FALSE)

preProB_DEgenes <- FindMarkers(combined.all, ident.1 = "pre-Pro-B_Infected", ident.2 = "pre-Pro-B_Not", verbose = FALSE)

unident_DEgenes <- FindMarkers(combined.all, ident.1 = "unidentified_Infected", ident.2 = "unidentified_Mock", verbose = FALSE)
unident_marks <- FindMarkers(combined.all, ident.1 = "unidentified")
saveRDS(unident_marks, "./unident_cluster_marks.rds")


######################################### Volcano Plots ###########################
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())

preM <- subset(combined.all, idents = "pre-Monocytes")
Idents(preM) <- "sample"
avg.preM <- as.data.frame(AverageExpression(preM, verbose = FALSE, slot = "data")$RNA)
avg.preM$gene <- rownames(avg.preM)

Mono <- subset(combined.all, idents = "Monocytes")
Idents(Mono) <- "sample"
avg.Mono <- as.data.frame(AverageExpression(Mono, verbose = FALSE, slot="data")$RNA)
avg.Mono$gene <- rownames(avg.Mono)

### Volcano Plot function 
volcano <- function(df, title, genes.to.label) {
  p1 <- ggplot(df, aes(Infected, Not)) + geom_point() + ggtitle(title)
  p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, xnudge = 0, ynudge = 0)
  p1
}
genes.to.label = c("CXCL10", "IFI27", "B2M", "ISG15", "ISG20", 'IRF8', "GBP1", "CXCL9", "SCT", "PSME2")
volcano(avg.Mono, "Monocytes", genes.to.label)
volcano(avg.preM, "pre-Monocytes", genes.to.label)

###################################################################################
# 10/3
saveRDS(combined.all, "./combinedALL.rds")
combined.all <- readRDS("./combinedALL.rds")
#################################################################################
# SCpubR Plots
library(SCpubr)
SCpubr::do_DimPlot(sample = combined.all, pt.size =1.5, label = T, label.size = 3, label.box = T, split.by = "sample", group.by = "main")
SCpubr::do_FeaturePlot(sample = combined.all,
                       features = c("AVP", "MPO"), label.color = "black", label = T)
combined.all <- SetIdent(combined.all, value = "sample")
de_genes <- FindMarkers(combined.all, ident.1 = "Infected", ident.2 = "Mock")
SCpubr::do_VolcanoPlot(sample = combined.all,
                            de_genes = de_genes,
                            n_genes = 15,
                            pval_cutoff = 1.3,
                            FC_cutoff = 1.5,
                            colors.use = c("red"),
                            font.size=12)
combined.all <- SetIdent(combined.all, value = combined.all$celltype.stim)
SCpubr::do_ExpressionHeatmap(sample = combined.all,
                                  features = c("MKI67"), 
                                  viridis_direction = -1,
                                  # group.by = "sample",
                                  flip= T, cell_size = 7)

######################## Trajectory Analysis ####################################
########################## TRAJ ANALYSIS MARCH 30, 2023 - Apr 10
## DE TRAJECTORY 
library(slingshot); library(SingleCellExperiment)
library(RColorBrewer); library(scales)
library(viridis); library(UpSetR)
library(pheatmap); library(msigdbr)
library(fgsea); library(knitr)
library(ggplot2); library(gridExtra)
library(tradeSeq); library(condiments)

DefaultAssay(combined.all) <- "RNA"
# combined.all <- FindVariableFeatures(combined.all, nfeatures = 500)
# var.features <- combined.all@assays$RNA@var.features
# combined.all <- subset(combined.all, features = combined.all@assays$RNA@var.features)

sce <- as.SingleCellExperiment(combined.all, assay = "RNA")
sce <- sce[var.features,]

# Plot Differential Topology
shuffle <- sample(ncol(sce))
layout(matrix(1, nrow = 1))
par(mar = c(4.5,4,1,1))

plot(reducedDims(sce)$UMAP[shuffle, ],
     asp = 1, pch = 16, xlab = "UMAP-1", ylab = "UMAP-2",
     col = alpha(c(1:2)[factor(colData(sce)$sample)][shuffle], alpha = .5))
legend("topright", pch = 16, col = 1:2, bty = "n", 
       legend = levels(factor(colData(sce)$sample)))

# Regions with a high score indicate that the local cell distribution according to treatment label 
# is unbalanced compared the overall distribution. Here, we see that, while there are some small regions of imbalance, 
# the global path along the development axis is well-balanced. 
# This means that we can fit a global trajectory to the full dataset. 

scores <- condiments::imbalance_score(
  reducedDims(sce)$UMAP, 
  condition = colData(sce)$sample,
  k = 20, smooth = 40)

grad <- viridis::plasma(10, begin = 0, end = 1)
names(grad) <- levels(cut(scores$scaled_scores, breaks = 10))
plot(reducedDims(sce)$UMAP, col = grad[cut(scores$scaled_scores, breaks = 10)],
     asp = 1, pch = 16, xlab = "UMAP-1", ylab = "UMAP-2", cex = .8)
legend(10,8, legend = names(grad), col = grad, pch = 20, bty = "n", cex = 2.5 / 3)

# Run Slingshot Trajectory Analysis
sce <- slingshot(sce, reducedDim = 'UMAP', clusterLabels = colData(sce)$integrated_snn_res.2, 
                 start.clus = '12', approx_points = 300)
sce <- slingshot(sce, reducedDim = 'UMAP', clusterLabels = colData(sce)$traj.clusters, 
                 start.clus = '15', approx_points = 300)
sce <- slingshot(sce, reducedDim = 'UMAP', clusterLabels = colData(sce)$main, dist.method = "mnn",
                 start.clus = 'HSC', approx_points = 300, extend = 'n')
sce <- slingshot(sce, reducedDim = 'UMAP', clusterLabels = colData(sce)$trajectory_clusters,
                 start.clus = '3', approx_points = 300,  dist.method = "mnn")

# create list of colors
color_list <- as.character(combined.all$colors)
names(color_list) <- combined.all$main

colors = c('#87CEFA', '#20B2AA', '#00FFFF', '#DB7093',  '#FFB6C1', '#1E90FF',  '#9ACD32',  '#A0522D',
           '#CC9966',  '#0000CD', '#98FB98', '#DC143C', '#E0E0E0',  '#FF69B4', '#FF7F50', '#7B68EE',
           '#6B8E23',  '#BC8F8F', '#8A2BE2',  '#FFFF00', '#FFA500',  '#BA55D3', '#B0C4DE')

vars <- as.data.frame(table(combined.all$main))
levels(vars$Var1)

# Plot Trajectories
plot(reducedDim(SlingshotDataSet(sce)), col = color_list, pch = 10, cex = 0.4)
legend("topright", legend = vars$Var1 , col = colors, pch = 16, bty = "n", cex = 2 / 3)
lines(SlingshotDataSet(sce), lwd = 2, type = 'lineages', col = 'black', show.constraints = T)

layout(matrix(1:2, nrow = 1))
par(mar = c(3,4,3,6))
windows(width = 5, height = 6)

# CUrve Trajectories
crv1 <- getCurves(sce)
crv1
crv <- SlingshotDataSet(crv1)
crv@lineages <- crv@lineages[c(1,2,3,4,5,6,7,9,10,11,12,14)]
crv@curves <- crv@curves[c(1,2,3,4,5,6,7,9,10,11,12,14)]
par(mar = c(3,4,3,10))

# Plot Trajectory Curves
plot(reducedDim(SlingshotDataSet(sce)), col = cell_colors_clust, pch = 10, cex = 0.4)
legend("topright", legend = vars$Var1 , col = cols, pch = 20, bty = "n", cex = 1,  inset = c(-0.23, 0), xpd = T, text.font = 2)
lines(crv, lwd = 2,  col = 'black')

##### Kolmogorov-Smirnov Test to assess whether the two groups of pseudotime values are derived from the same distribution
ks.test(slingPseudotime(sce)[colData(sce)$sample == "Infected", 1:5],
        slingPseudotime(sce)[colData(sce)$sample == "Mock", 1:5])

ks.test(slingPseudotime(sce)[colData(sce)$sample == "Infected", 6:7],
        slingPseudotime(sce)[colData(sce)$sample == "Mock", 6:7])

ks.test(slingPseudotime(sce)[colData(sce)$sample == "Infected", c(12,14)],
        slingPseudotime(sce)[colData(sce)$sample == "Mock", c(12,14)])

ks.test(slingPseudotime(sce)[colData(sce)$sample == "Infected", 9],
        slingPseudotime(sce)[colData(sce)$sample == "Mock", 9])

ks.test(slingPseudotime(sce)[colData(sce)$sample == "Infected", 10:11],
        slingPseudotime(sce)[colData(sce)$sample == "Mock", 10:11])

ks.test(slingPseudotime(sce)[colData(sce)$sample == "Infected", 6],
        slingPseudotime(sce)[colData(sce)$sample == "Mock", 6])

ks.test(slingPseudotime(sce)[colData(sce)$sample == "Infected", 7],
        slingPseudotime(sce)[colData(sce)$sample == "Mock", 7])

ks.test(slingPseudotime(sce)[colData(sce)$sample == "Infected", 8],
        slingPseudotime(sce)[colData(sce)$sample == "Mock", 8])

ks.test(slingPseudotime(sce)[colData(sce)$sample == "Infected", 9],
        slingPseudotime(sce)[colData(sce)$sample == "Mock", 9])

ks.test(slingPseudotime(sce)[colData(sce)$sample == "Infected", 10],
        slingPseudotime(sce)[colData(sce)$sample == "Mock", 10])

ks.test(slingPseudotime(sce)[colData(sce)$sample == "Infected", 11],
        slingPseudotime(sce)[colData(sce)$sample == "Mock", 11])

ks.test(slingPseudotime(sce)[colData(sce)$sample == "Infected", 12],
        slingPseudotime(sce)[colData(sce)$sample == "Mock", 12])

ks.test(slingPseudotime(sce)[colData(sce)$sample == "Infected", 13],
        slingPseudotime(sce)[colData(sce)$sample == "Mock", 13])

ks.test(slingPseudotime(sce)[colData(sce)$sample == "Infected", 14],
        slingPseudotime(sce)[colData(sce)$sample == "Mock", 14])

#### Add Trajectories to Seurat Object
combined.all$pt1 <- sce$slingPseudotime_1
combined.all$pt2 <- sce$slingPseudotime_2
combined.all$pt3 <- sce$slingPseudotime_3
combined.all$pt4 <- sce$slingPseudotime_4
combined.all$pt5 <- sce$slingPseudotime_5
combined.all$pt6 <- sce$slingPseudotime_6
combined.all$pt7 <- sce$slingPseudotime_7
combined.all$pt8 <- sce$slingPseudotime_8
combined.all$pt9 <- sce$slingPseudotime_9
combined.all$pt10 <- sce$slingPseudotime_10
combined.all$pt11 <- sce$slingPseudotime_11
combined.all$pt12 <- sce$slingPseudotime_12
combined.all$pt13 <- sce$slingPseudotime_13
combined.all$pt14 <- sce$slingPseudotime_14

combined.all <- SetIdent(combined.all, value = "main")
# PLot cells within different trajectories
FeaturePlot(combined.all, c("pt1", "pt2", "pt3", "pt4", "pt5", "pt6"), label = T, label.size = 2.5)
FeaturePlot(combined.all, c("pt7", "pt8", "pt9", "pt10", "pt11", "pt12", "pt13", "pt14"), label = T, label.size = 2.5)

########### Plot Density
layout(matrix(c(1, 1, 2, 3), 2))
# integrated
plot(reducedDim(SlingshotDataSet(sce)), col = cell_colors_clust, pch = 10, cex = 0.4)
legend("topright", legend = vars$Var1 , col = cols, pch = 16, bty = "n", cex = 2.7 / 3)
lines(SlingshotDataSet(sce), lwd = 2, type = 'lineages', col = 'black', show.constraints = T)

par(mar = c(4.5, 4, 1, 1))
plot(reducedDims(sce)$UMAP[shuffle, ], asp = 1, pch = 16, xlab = "UMAP-1", ylab = "UMAP-2",
     col = hcl.colors(100, alpha = .5)[cut(sce$slingPseudotime_6, breaks = 100)][shuffle])
lines(SlingshotDataSet(sce))
plot(reducedDims(sce)$UMAP[shuffle, ], asp = 1, pch = 16, xlab = "UMAP-1", ylab = "UMAP-2",
     col = hcl.colors(100, alpha = .5)[cut(sce$slingPseudotime_7, breaks = 100)][shuffle])
lines(SlingshotDataSet(sce))
plot(reducedDims(sce)$UMAP[shuffle, ], asp = 1, pch = 16, xlab = "UMAP-1", ylab = "UMAP-2",
     col = hcl.colors(100, alpha = .5)[cut(sce$slingPseudotime_12, breaks = 100)][shuffle])
lines(SlingshotDataSet(sce))
plot(reducedDims(sce)$UMAP[shuffle, ], asp = 1, pch = 16, xlab = "UMAP-1", ylab = "UMAP-2",
     col = hcl.colors(100, alpha = .5)[cut(sce$slingPseudotime_14, breaks = 100)][shuffle])
lines(SlingshotDataSet(sce))
plot(reducedDims(sce)$UMAP[shuffle, ], asp = 1, pch = 16, xlab = "UMAP-1", ylab = "UMAP-2",
     col = hcl.colors(100, alpha = .5)[cut(sce$slingPseudotime_9, breaks = 100)][shuffle])
lines(SlingshotDataSet(sce))

ds <- list(Mock = density(na.omit(slingPseudotime(sce)[colData(sce)$sample == "Mock", 9])),
           HIV_Infected = density(na.omit(slingPseudotime(sce)[colData(sce)$sample == "Infected", 9])))
xlim <- range(c(ds$Mock$x, ds$HIV_Infected$x))
ylim <- range(c(ds$Mock$y, ds$HIV_Infected$y))
plot(xlim, ylim, col = "white", xlab = "Pseudotime", ylab = "")
polygon(c(min(ds$Mock$x),ds$Mock$x,max(ds$Mock$x)),
        c(0,ds$Mock$y,0), col = rgb(0,0,0,.5))
polygon(c(min(ds$HIV_Infected$x),ds$Mock$x,max(ds$Mock$x)),
        c(0,ds$HIV_Infected$y,0), col = alpha(brewer.pal(4,'Set1')[1], alpha = .5))
legend("topright", legend = c("Mock", "Infected"), 
       fill = alpha(c(1, brewer.pal(3, "Set1")[1]), alpha = .5), bty = "n")

par(mar = c(5, 4, 4, 2) + .1)


### Plot all Curves
sce1 <- SlingshotDataSet(sce)
nc <- 3
pt <- slingPseudotime(sce1)
nms <- colnames(pt)
nr <- ceiling(length(nms)/nc)
pal <- viridis(100, end = 0.95)
par(mfrow = c(nr, nc))
for (i in nms) {
  colors <- pal[cut(pt[,i], breaks = 100)]
  plot(reducedDim(sce1), col = colors, pch = 16, cex = 0.5, main = i)
  lines(sce1, lwd = 2, col = 'black', type = 'lineages') 
}

# FINAL VISUALIZATIONS ############################################################
# Dimensional Reduction UMAP
combined.all <- SetIdent(combined.all, value = combined.all$main)
colors = c('#87CEFA', '#20B2AA', '#00FFFF', '#DB7093',  '#FFB6C1', '#1E90FF',  '#9ACD32',  '#A0522D',
           '#CC9966',  '#0000CD', '#98FB98', '#DC143C', '#E0E0E0',  '#FF69B4', '#FF7F50', '#7B68EE',
           '#6B8E23',  '#BC8F8F', '#8A2BE2',  '#FFFF00', '#FFA500',  '#BA55D3', '#B0C4DE')
Seurat::DimPlot(combined.all, label = F, pt.size = 1.5, split.by = "three", cols = colors) +  guides(color = guide_legend(override.aes = list(size=5), ncol=1))
Seurat::DimPlot(combined.all, label = F, pt.size = 1.5, split.by = "sample", cols = colors) +  guides(color = guide_legend(override.aes = list(size=4), ncol=1))

p <- Seurat::DimPlot(combined.all, label = F, pt.size = 2.4, label.size = 4.5, repel = T) +  guides(color = guide_legend(override.aes = list(size=4), ncol=1))
LabelClusters(p, id = "ident",  fontface = "bold", color = "black", size = 4.5)


################################### Dot Plot
Seurat::DotPlot(combined.all, assay = "RNA", features = c('CD38','CD34','AVP', 'PROM1', 'MPO','CEBPD', 'HOPX',
                                  'DNTT', 'CYGB', 'VPREB1', 'VPREB3', 'CD79A', 'MKI67', 
                                  'RAG1','HVCN1', 'ACY3', 'IRF8','BATF3','LILRA4',
                                  'FCER1A', 'ELANE', 'S100A8', 'S100A9','VCAN', 
                                  'GATA1', 'GATA2','KLF1','APOC1','TPSB2','EPX'), cols = c("blue", "red"), dot.scale = 8) + RotatedAxis()
###################################### Volcano Plot
Idents(combined.all) <- "main"
table(combined.all$main)
MPP <- WhichCells(combined.all, idents = c("MPP1", "MPP2", "MPP3", "MPP4"))
dc <- WhichCells(combined.all, idents = c("CDP", "pre-DCs", "pre-DC Cycle"))
mono <- WhichCells(combined.all, idents = c("MP", "pre-Monocytes", 'Mye-DC'))
clpb <- WhichCells(combined.all, idents = c("LMPP", "CLP.1", "CLP.2", "CLP cycle", "pre-Pro-B cycle", "pre-Pro-B.1", "pre-Pro-B.2", "Pro-B"))

combined.all$lin_main <- combined.all$main
combined.all@meta.data[MPP, "lin_main"] <- "MPP"
combined.all@meta.data[dc, "lin_main"] <- "DC-Lineage"
combined.all@meta.data[mono, "lin_main"] <- "Monocyte-Lineage"
combined.all@meta.data[clpb, "lin_main"] <- "B-Cell_Lineage"
Idents(combined.all) <- "lin_main"
combined.all$celltype.stim.combined <- paste(Idents(combined.all), combined.all$sample, sep = "_")
Idents(combined.all) <- "celltype.stim.combined"
table(combined.all@meta.data$celltype.stim.combined)

combined.all <- SetIdent(combined.all, value = "three")


### PROPORTIONS PLOT ############################################################
# show different proportions of cells in Control & Infected for each Cluster
prop = prop.table(table(Idents(combined.all),combined.all$sample),1) %>% data.frame()
names(prop) = c('cluster','group','prop')
ranked = as.character(arrange(prop[prop$group=='Not',],-prop)$cluster)
prop$cluster = factor(prop$cluster,levels=ranked)
counts = table(Idents(combined.all)) %>% data.frame()
names(counts) = c('cluster','counts')
counts$cluster = factor(counts$cluster,levels=ranked)

prop_plot = ggplot(prop,aes(fill=group,x=cluster,label=counts$counts))+
  geom_bar(data=prop[prop$group=='Infected',],aes(y=prop),stat='identity',color='white')+
  geom_bar(data=prop[prop$group=='Not',],aes(y=-prop),stat='identity',color='white')+
  scale_x_discrete(expand=c(0,0.6))+
  scale_y_continuous(expand=c(0,0.01),labels=function(x) abs(x*100))+
  scale_fill_manual(values=c('red3','royalblue4'),labels=c('Infected','Control'))+
  labs(y='Proportion (%)')+
  theme(panel.background=element_rect(fill = "white"),
        panel.grid=element_blank(),
        axis.line=element_line(size=1),
        axis.ticks=element_line(size=1),
        axis.ticks.length=unit(1.5,'mm'),
        axis.text.y=element_text(size=16,hjust=0.5,color='black', face = "bold"),
        axis.text.x=element_text(size=18,hjust=1,color='black',angle=45,vjust=1, face = "bold"),
        legend.position='top',
        legend.text=element_text(size=18,color='black', face = "bold"),
        legend.title=element_blank(),
        axis.title.y= element_text(size=16,hjust=0.2,color='black',margin=margin(r=10), face = "bold"),
        axis.title.x= element_blank()) # + coord_flip()

################################## violin plot #####################################
colors = c('#E0E0E0', '#CC9966', '#DC143C', '#FF7F50', '#FFA500', '#FFFF00', '#9ACD32', '#6B8E23',
           '#98FB98', '#20B2AA', '#00FFFF','#1E90FF', '#87CEFA', '#0000CD', '#8A2BE2', '#7B68EE', '#BA55D3',
           '#DB7093', '#FF69B4', '#FFB6C1', '#FAEBD7', '#A0522D', '#BC8F8F', '#B0C4DE')


p = VlnPlot(combined.all,c('CD34','AVP', 'PROM1', 'MPO','CEBPD', 'DNTT',
                           'CYGB', 'VPREB1', 'VPREB3', 'CD79A', 'MKI67', 
                           'VPREB3','RAG1', 'ACY3', 'IRF8','BATF3','LILRA4', 'FCER1A','ELANE', 'S100A8', 'S100A9','VCAN',
                           'GATA1', 'GATA2','APOC1','EPX'),
            pt.size=0,stack=T,flip=T,fill.by='ident',assay='RNA')+
  NoLegend()+
  scale_x_discrete(position='top')+
  scale_fill_manual(values=colors)+
  theme(panel.background=element_rect(color='black',size=0.8),
        panel.spacing=unit(0,'lines'),
        strip.text.y.right=element_text(size=17,hjust=0),
        axis.text.x=element_text(size=16,angle=45,hjust=0, face= "bold"),
        axis.text.y=element_blank(),
        axis.line=element_blank(),
        axis.ticks.x=element_line(size=1.2),
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        plot.title=element_blank())+
  coord_cartesian(clip='off')



################################################################################
############################ UCELL GENE SIGNATURE SCORING ######################
# Score Cells as Type 1 or Type 2 Interferon Responding based on expression of signature Gene List
library(UCell)
library(dplyr)
combined.all <- SetIdent(combined.all, value = "sample")
ctrl <- subset(combined.all, ident = "Not")
infected <- subset(combined.all, ident = "Infected")

# Find IFN-GAMMA & IFN-BETA responding Cells
expr.matrix <- combined.all[["RNA"]]@counts
expr.matrix.ctrl <- ctrl[["RNA"]]@counts
expr.matrix.inf <- infected[["RNA"]]@counts
ifnB <- c("IFIT1", "MS4A", "PHF11", "MMP13", "KDR", "IFIT2", "ISG20", "HERC6", "MOV10", "GADD45G", 
          "ITPR1", "IFIT3", "ABCG1", "TLR3", "OASL2P", "SLC28A2", "IFIH1","IFI44", "USP18", "MNDA",
          "IL1RN", "OAS1", "IRF7", "NT5C3A", "THBS1", "SP140", "SLFN5", "RSAD2", "CMPK2")
ifnG <- c("CIITA", "CD74", "SOCS1", "RNASE6", "TAGAP", "CXCL9", "IRF1", "H2AB1", "HCAR2", "IRF8", 
          "BATF2", "ARHGEF3", "PRR5L", "CDKN1A", "LRRC8C", "UPP1", "NAAA", "CD83", "GBP7")
gene.sets <- list(Beta_Signature = ifnB, Gamma_Signature = ifnG)
scores <- ScoreSignatures_UCell(expr.matrix, features = gene.sets)
scores.ctrl <- ScoreSignatures_UCell(expr.matrix.ctrl, features = gene.sets)
scores.inf <- ScoreSignatures_UCell(expr.matrix.inf, features = gene.sets)
head(scores)

scores_df <- data.frame(scores)
scores_df_ctrl <- data.frame(scores.ctrl)
scores_df_inf <- data.frame(scores.inf)

sum(scores_df$Beta_Signature_UCell > 0.3)

# No Cells have a score higher than 0.5 for either gene set
# 174 Cells have a score > 0.3 for Beta and only 1 cell > 0.3 for Gamma

########## ALL TOGETHER
beta <- rep("Beta", length(scores_df$Beta_Signature_UCell))
gamma <- rep("Gamma", length(scores_df$Gamma_Signature_UCell))
scoring <- c(scores_df$Beta_Signature_UCell, scores_df$Gamma_Signature_UCell)
samp <- c(beta, gamma)
cells <- row.names(scores_df)
df <- data.frame(score = scoring, sample = samp, cell = c(cells, cells))

meta.data <- data.frame(combined.all@meta.data)
meta.data$cell <- row.names(meta.data) 
head(meta.data)
meta.data <- meta.data %>% select(cell, lin_main, celltype.stim.combined, main, celltype.stim, sample)
head(df)

df <- left_join(df, meta.data, by = "cell")
head(df)
df$SAMP <- paste(df$sample.x, df$sample.y, sep = "_")
########## CONTROL & INFECTED SEPARATE 
beta_ctrl <- rep("Beta_CTRL", length(scores_df_ctrl$Beta_Signature_UCell))
gamma_ctrl <- rep("Gamma_CTRL", length(scores_df_ctrl$Gamma_Signature_UCell))
scoring_ctrl <- c(scores_df_ctrl$Beta_Signature_UCell, scores_df_ctrl$Gamma_Signature_UCell)
beta_inf <- rep("Beta_INF", length(scores_df_inf$Beta_Signature_UCell))
gamma_inf <- rep("Gamma_INF", length(scores_df_inf$Gamma_Signature_UCell))
scoring_inf <- c(scores_df_inf$Beta_Signature_UCell, scores_df_inf$Gamma_Signature_UCell)
samp <- c(beta_ctrl, gamma_ctrl, beta_inf, gamma_inf)
scoring <- c(scoring_ctrl, scoring_inf)
cells_ctrl <- row.names(scores_df_ctrl)
cells_inf <- row.names(scores_df_inf)
df <- data.frame(score = scoring, sample = samp, cell = c(cells_ctrl, cells_inf))

meta.data <- data.frame(combined.all@meta.data)
meta.data$cell <- row.names(meta.data) 
head(meta.data)
meta.data <- meta.data %>% select(cell, lin_main, celltype.stim.combined, main, celltype.stim)
head(df)

df <- left_join(df, meta.data, by = "cell")
head(df)

# Violin PLot
ggplot(df, aes(x=SAMP, y=score)) + geom_violin(draw_quantiles = c(0.85))  + xlab("Score") + ylab("Sample") + theme(axis.text.y=element_text(size=12,color='black', face = "bold"),
                                                                                     axis.text.x=element_text(size=12,color='black', face = "bold"),
                                                                                     axis.title.y= element_text(size=16,color='black',face = "bold"),
                                                                                     axis.title.x= element_text(size=16,color='black',face = "bold"))  #+ geom_jitter(height = 0, width = 0.05)

# Get 99% Quantile
quantile(scores_df$Beta_Signature_UCell, probs = 0.85)  # 0.1741839
quantile(scores_df$Gamma_Signature_UCell, probs = 0.85) # 0.1478965

quantile(scores_df_ctrl$Beta_Signature_UCell, probs = 0.25)  # 0.0249
quantile(scores_df_ctrl$Gamma_Signature_UCell, probs = 0.25) # 0.0635
quantile(scores_df_inf$Beta_Signature_UCell, probs = 0.25)   # 0.0589
quantile(scores_df_inf$Gamma_Signature_UCell, probs = 0.25)  # 0.0673

quantile(scores_df_ctrl$Beta_Signature_UCell, probs = 0.5)  # 0.0408
quantile(scores_df_ctrl$Gamma_Signature_UCell, probs = 0.5) # 0.0931
quantile(scores_df_inf$Beta_Signature_UCell, probs = 0.5)   # 0.0978
quantile(scores_df_inf$Gamma_Signature_UCell, probs = 0.5)  # 0.0922

quantile(scores_df_ctrl$Beta_Signature_UCell, probs = 0.75)  # 0.0635
quantile(scores_df_ctrl$Gamma_Signature_UCell, probs = 0.75) # 0.1213
quantile(scores_df_inf$Beta_Signature_UCell, probs = 0.75)   # 0.1516
quantile(scores_df_inf$Gamma_Signature_UCell, probs = 0.75)  # 0.1136

quantile(scores_df_ctrl$Beta_Signature_UCell, probs = 0.90)  # 0.0904
quantile(scores_df_ctrl$Gamma_Signature_UCell, probs = 0.90) # 0.1530
quantile(scores_df_inf$Beta_Signature_UCell, probs = 0.90)   # 0.2211
quantile(scores_df_inf$Gamma_Signature_UCell, probs = 0.90)  # 0.1423

# How many cells score over that 99% quantile expression
sum(scores_df$Beta_Signature_UCell > 0.0355)  # 98 Cells
sum(scores_df$Gamma_Signature_UCell > 0.0656) # 99 Cells

sigB <- which(scores_df$Beta_Signature_UCell > 0.1741839) 
sigG <- which(scores_df$Gamma_Signature_UCell > 0.1478965)
betaa <- row.names(scores_df)[sigB]
gammaa <- row.names(scores_df)[sigG]

sigBctrl <- which(scores_df_ctrl$Beta_Signature_UCell > 0.1434385) 
sigGctrl <- which(scores_df_ctrl$Gamma_Signature_UCell > 0.1356491)
sigBinf <- which(scores_df_inf$Beta_Signature_UCell > 0.1434385) 
sigGinf <- which(scores_df_inf$Gamma_Signature_UCell > 0.1356491)

betaCTRL <- row.names(scores_df_ctrl)[sigBctrl]
gammaCTRL <- row.names(scores_df_ctrl)[sigGctrl]
betaINF <- row.names(scores_df_inf)[sigBinf]
gammaINF <- row.names(scores_df_inf)[sigGinf]

# Combined
both <- intersect(betaa, gammaa)
Bonly <- setdiff(betaa, gammaa)
Gonly <- setdiff(gammaa, betaa)

combined.all$Ucell <- "neither"
combined.all@meta.data[both, "Ucell"] <- "Both"
combined.all@meta.data[Bonly, "Ucell"] <- "IFN-Beta"
combined.all@meta.data[Gonly, "Ucell"] <- "IFN-Gamma"

# Control
both <- intersect(betaCTRL, gammaCTRL)
Bonly <- setdiff(betaCTRL, gammaCTRL)
Gonly <- setdiff(gammaCTRL, betaCTRL)

ctrl$Ucell <- "neither"
ctrl@meta.data[both, "Ucell"] <- "Both"
ctrl@meta.data[Bonly, "Ucell"] <- "IFN-Beta"
ctrl@meta.data[Gonly, "Ucell"] <- "IFN-Gamma"

# INfected
both <- intersect(betaINF, gammaINF)
Bonly <- setdiff(betaINF, gammaINF)
Gonly <- setdiff(gammaINF, betaINF)

infected$Ucell <- "neither"
infected@meta.data[both, "Ucell"] <- "Both"
infected@meta.data[Bonly, "Ucell"] <- "IFN-Beta"
infected@meta.data[Gonly, "Ucell"] <- "IFN-Gamma"

### Dimplot Combined
combined.all <- SetIdent(combined.all, value = "Ucell")
ctrl <- SetIdent(ctrl, value = "Ucell")
infected <- SetIdent(infected, value = "Ucell")
DimPlot(combined.all, label = F, pt.size = 1, split.by = sample, order = c( "IFN-Beta", "IFN-Gamma", "Both", "neither"), cols = c("grey", "red", "blue", "green"))

combined.all <- SetIdent(combined.all, value = "main")
DimPlot(combined.all, label = F, pt.size = 1) +  guides(color = guide_legend(override.aes = list(size=5), ncol=1))

table(combined.all$main, combined.all$Ucell)

#### Control
ctrl <- SetIdent(ctrl, value = "Ucell")
DimPlot(ctrl, label = F, pt.size = 1, order = c( "IFN-Beta", "IFN-Gamma", "Both", "neither"), cols = c("grey", "red", "blue", "green"))

ctrl <- SetIdent(ctrl, value = "main")
DimPlot(ctrl, label = F, pt.size = 1) +  guides(color = guide_legend(override.aes = list(size=5), ncol=1))

table(ctrl$main, ctrl$Ucell)

#### Infected
infected <- SetIdent(infected, value = "Ucell")

DimPlot(infected, label = F, pt.size = 1, order = c( "IFN-Beta", "IFN-Gamma", "Both", "neither"), cols = c("grey", "red", "blue", "green"))

infected <- SetIdent(infected, value = "main")
DimPlot(infected, label = F, pt.size = 1) +  guides(color = guide_legend(override.aes = list(size=5), ncol=1))

table(infected$main, infected$Ucell)

####################### Frequency Plot #########################################
table(combined.all$main, combined.all$Ucell)
prop.table(table(combined.all$main, combined.all$Ucell))

library(ggsci)
library(ggbreak)
df <- as.data.frame(table(combined.all$main, combined.all$Ucell))
p <- ggplot(df) +
  aes(x = Var1, fill = Var2, weight = Freq, by = Var2) +
  geom_bar(position = "fill") +xlab("Cluster") + ylab("Proportion") 
p +scale_fill_lancet() + theme(  axis.text.y=element_text(size=16,color='black', face = "bold"),
                                 axis.text.x=element_text(size=14,color='black', face = "bold", hjust = 1, angle = 45),
                                 legend.position='right',
                                 legend.text=element_text(size=14,color='black', face = "bold"),
                                 legend.title=element_blank(),
                                 axis.title.y=element_text(size=20,color='black',face = "bold",angle = 90), axis.title.x=element_text(size=20,color='black',face = "bold"))

ggplot(df) +aes(x = Var1, weight = Freq, fill=Var2) + geom_bar() +
  theme(strip.text = element_text(size = 10, face = "bold")) +
  xlab("Cluster") + ylab("proportion") + ylab("Cell Count") +scale_fill_lancet() + theme(  axis.text.y=element_text(size=16,color='black', face = "bold"),
                                                                                           axis.text.x=element_text(size=14,color='black', face = "bold", hjust = 1, angle = 45),
                                                                                           legend.position='right',
                                                                                           legend.text=element_text(size=14,color='black', face = "bold"),
                                                                                           legend.title=element_blank(),
                                                                                           axis.title.y=element_text(size=20,color='black',face = "bold",angle = 90), axis.title.x=element_text(size=20,color='black',face = "bold"))

# + scale_y_break(c(500,1000), scales=2) 

#### PROP PLOT
### PROPORTIONS - does not work
prop = prop.table(table(combined.all$main,combined.all$Ucell),1) %>% data.frame()
names(prop) = c('cluster','group','prop')
ranked = as.character(arrange(prop[prop$group=='Not',],-prop)$cluster)
prop$cluster = factor(prop$cluster,levels=ranked)
counts = table(Idents(combined.all)) %>% data.frame()
names(counts) = c('cluster','counts')
counts$cluster = factor(counts$cluster,levels=ranked)

prop_plot = ggplot(prop,aes(fill=group,x=cluster,label=counts$counts))+
  geom_bar(data=prop[prop$group=='Both',],aes(y=prop),stat='identity',color='white')+
  geom_bar(data=prop[prop$group=='IFN-Beta',],aes(y=-prop),stat='identity',color='white')+
  geom_bar(data=prop[prop$group=='IFN-Gamma',],aes(y=prop),stat='identity',color='white')+
  geom_bar(data=prop[prop$group=='neither',],aes(y=-prop),stat='identity',color='white')+
  scale_x_discrete(expand=c(0,0.6))+
  scale_y_continuous(expand=c(0,0.01),labels=function(x) abs(x*100))+
  scale_fill_manual(values=c('red3','royalblue4', 'green', 'purple'),labels=c('Both',  'IFN-Beta', 'IFN-Gamma',   'neither'), breaks = c('Both',  'IFN-Beta', 'IFN-Gamma',   'neither'))+
  labs(y='Proportion (%)')+
  theme(panel.background=element_rect(fill = "white"),
        panel.grid=element_blank(),
        axis.line=element_line(size=1),
        axis.ticks=element_line(size=1),
        axis.ticks.length=unit(1.5,'mm'),
        axis.text.y=element_text(size=16,hjust=0.5,color='black', face = "bold"),
        axis.text.x=element_text(size=18,hjust=1,color='black',angle=45,vjust=1, face = "bold"),
        legend.position='top',
        legend.text=element_text(size=18,color='black', face = "bold"),
        legend.title=element_blank(),
        axis.title.y= element_text(size=16,hjust=0.2,color='black',margin=margin(r=10), face = "bold"),
        axis.title.x= element_blank()) # + coord_flip()

######################### Final DimPlot ########################################
combined.all <- SetIdent(combined.all, value = combined.all$main)
colors = c('#87CEFA', '#20B2AA', '#00FFFF', '#DB7093',  '#FFB6C1', '#1E90FF',  '#9ACD32',  '#A0522D',
           '#CC9966',  '#0000CD', '#98FB98', '#DC143C', '#E0E0E0',  '#FF69B4', '#FF7F50', '#7B68EE',
           '#6B8E23',  '#BC8F8F', '#8A2BE2',  '#FFFF00', '#FFA500',  '#BA55D3', '#B0C4DE')
Seurat::DimPlot(combined.all, label = F, pt.size = 1.5, split.by = "three", cols = colors) +  guides(color = guide_legend(override.aes = list(size=5), ncol=1))
Seurat::DimPlot(combined.all, label = F, pt.size = 1.5, split.by = "sample", cols = colors) +  guides(color = guide_legend(override.aes = list(size=4), ncol=1))

p <- Seurat::DimPlot(combined.all, label = F, pt.size = 2.4, label.size = 4.5, repel = T) +  guides(color = guide_legend(override.aes = list(size=4), ncol=1))
LabelClusters(p, id = "ident",  fontface = "bold", color = "black", size = 4.5)
p
############################### Add numbers to DimPlots ########################
hsc <- WhichCells(combined.all, idents = "HSC")
mpp1 <- WhichCells(combined.all, idents = "MPP1")
mpp2 <- WhichCells(combined.all, idents = "MPP2")
mpp3 <- WhichCells(combined.all, idents = "MPP3")
mpp4 <- WhichCells(combined.all, idents = "MPP4")
lmpp <- WhichCells(combined.all, idents = "LMPP")
clp1 <- WhichCells(combined.all, idents = "CLP.1")
clp2 <- WhichCells(combined.all, idents = "CLP.2")
clpcyc <- WhichCells(combined.all, idents = "CLP cycle")
preprob1 <- WhichCells(combined.all, idents = "pre-Pro-B.1")
preprob2 <- WhichCells(combined.all, idents = "pre-Pro-B.2")
preprobcyc <- WhichCells(combined.all, idents = "pre-Pro-B cycle")
prob <- WhichCells(combined.all, idents = "Pro-B")
cdp <- WhichCells(combined.all, idents = "CDP")
myedc <- WhichCells(combined.all, idents = "Mye-DC")
predc <- WhichCells(combined.all, idents = "pre-DCs")
predcyc <- WhichCells(combined.all, idents = "pre-DC Cycle")
gmp <- WhichCells(combined.all, idents = "GMP")
mp <- WhichCells(combined.all, idents = "MP")
premono <- WhichCells(combined.all, idents = "pre-Monocytes")
memp <- WhichCells(combined.all, idents = "MEMP")
unident <- WhichCells(combined.all, idents = "unidentified")


combined.all@meta.data[hsc, "umapd"] <- 1
combined.all@meta.data[mpp1, "umapd"] <- 2
combined.all@meta.data[mpp2, "umapd"] <- 3
combined.all@meta.data[mpp3, "umapd"] <- 4
combined.all@meta.data[mpp4, "umapd"] <- 5
combined.all@meta.data[lmpp, "umapd"] <- 6
combined.all@meta.data[clp1, "umapd"] <- 8
combined.all@meta.data[clp2, "umapd"] <- 9
combined.all@meta.data[clpcyc, "umapd"] <- 10
combined.all@meta.data[preprob1, "umapd"] <- 18
combined.all@meta.data[preprob2, "umapd"] <- 19
combined.all@meta.data[preprobcyc, "umapd"] <- 20
combined.all@meta.data[prob, "umapd"] <- 21
combined.all@meta.data[cdp, "umapd"] <- 7
combined.all@meta.data[myedc, "umapd"] <- 17
combined.all@meta.data[predc, "umapd"] <- 12
combined.all@meta.data[predcyc, "umapd"] <- 11
combined.all@meta.data[gmp, "umapd"] <- 14
combined.all@meta.data[mp, "umapd"] <- 15
combined.all@meta.data[premono, "umapd"] <- 16
combined.all@meta.data[memp, "umapd"] <- 13
combined.all@meta.data[unident, "umapd"] <- 22

combined.all@meta.data[hsc, "umapds"] <- "1 : HSC"
combined.all@meta.data[mpp1, "umapds"] <- "2 : MPP1"
combined.all@meta.data[mpp2, "umapds"] <- "3 : MPP2"
combined.all@meta.data[mpp3, "umapds"] <- "4 : MPP3"
combined.all@meta.data[mpp4, "umapds"] <- "5 : MPP4"
combined.all@meta.data[lmpp, "umapds"] <- "6 : LMPP"
combined.all@meta.data[clp1, "umapds"] <- "8 : CLP.1"
combined.all@meta.data[clp2, "umapds"] <- "9 : CLP.2"
combined.all@meta.data[clpcyc, "umapds"] <- "10 : CLP Cycle"
combined.all@meta.data[preprob1, "umapds"] <- "18 : Pre-Pro-B.1"
combined.all@meta.data[preprob2, "umapds"] <- "19 : Pre-Pro-B.2"
combined.all@meta.data[preprobcyc, "umapds"] <- "20 : Pre-Pro-B Cycle"
combined.all@meta.data[prob, "umapds"] <- "21 : Pro-B"
combined.all@meta.data[cdp, "umapds"] <- "7 : CDP"
combined.all@meta.data[myedc, "umapds"] <- "17 : Mye-DC"
combined.all@meta.data[predc, "umapds"] <- "12 : pre-DC"
combined.all@meta.data[predcyc, "umapds"] <- "11 : pre-DC Cycle"
combined.all@meta.data[gmp, "umapds"] <- "14 : GMP"
combined.all@meta.data[mp, "umapds"] <- "15 : MP"
combined.all@meta.data[premono, "umapds"] <- "16 : pre-Monocytes"
combined.all@meta.data[memp, "umapds"] <- "13 : MEMP"
combined.all@meta.data[unident, "umapds"] <- "22 : unidentified"

combined.all$umapds <- factor(combined.all$umapds, levels = c("1 : HSC", "2 : MPP1", "3 : MPP2", "4 : MPP3", "5 : MPP4", "6 : LMPP", "7 : CDP",
                                                              "8 : CLP.1", "9 : CLP.2", "10 : CLP Cycle", "11 : pre-DC Cycle", "12 : pre-DC", 
                                                              "13 : MEMP", "14 : GMP", "15 : MP", "16 : pre-Monocytes", "17 : Mye-DC", "18 : Pre-Pro-B.1",
                                                              "19 : Pre-Pro-B.2", "20 : Pre-Pro-B Cycle", "21 : Pro-B", "22 : unidentified"))

combined.all <- SetIdent(combined.all, value = combined.all$umapds)
colors = c('#87CEFA', '#20B2AA', '#00FFFF', '#DB7093',  '#FFB6C1', '#1E90FF',  '#9ACD32',  '#A0522D',
           '#CC9966',  '#0000CD', '#98FB98', '#DC143C', '#E0E0E0',  '#FF69B4', '#FF7F50', '#7B68EE',
           '#6B8E23',  '#BC8F8F', '#8A2BE2',  '#FFFF00', '#FFA500',  '#BA55D3', '#B0C4DE')
colors = c(colors[13], colors[9], colors[12], colors[21], colors[15], colors[20], colors[10], colors[7], colors[17], colors[11],
           colors[22], colors[16], colors[18], colors[4], colors[14], colors[5], colors[19], colors[2], colors[3], colors[6], colors[1], colors[8])

Seurat::DimPlot(combined.all, label = F, pt.size = 1.5, split.by = "three", cols = colors) +  guides(color = guide_legend(override.aes = list(size=5), ncol=1))
Seurat::DimPlot(combined.all, label = F, pt.size = 1.5, split.by = "sample", cols = colors) +  guides(color = guide_legend(override.aes = list(size=4), ncol=1))

p <- Seurat::DimPlot(combined.all, label = F, pt.size = 2, label.size = 4.5, repel = F, cols = colors) +  guides(color = guide_legend(override.aes = list(size=4), ncol=1))
LabelClusters(p, id = "ident",  fontface = "bold", color = "black", size = 5, repel = F) +
  theme(axis.text.y=element_text(size=12,color='black', face = "bold"),
        axis.text.x=element_text(size=12,color='black', face = "bold"),
        legend.text=element_text(size=16,color='black', face = "bold"),
        axis.title.y= element_text(size=16,color='black',face = "bold"),
        axis.title.x= element_text(size=16,color='black',face = "bold"))

Seurat::DimPlot(combined.all, label = F, pt.size = 2, label.size = 4.5, repel = F, cols = colors, split.by = "sample") +  guides(color = guide_legend(override.aes = list(size=4), ncol=1))+
  theme(axis.text.y=element_text(size=12,color='black', face = "bold"),
        axis.text.x=element_text(size=12,color='black', face = "bold"),
        legend.text=element_text(size=16,color='black', face = "bold"),
        axis.title.y= element_text(size=16,color='black',face = "bold"),
        axis.title.x= element_text(size=16,color='black',face = "bold"))


###### DISTINGUISH IFN-ALPHA + IFN-GAMMA ##########################################
# load libraries
library(dplyr)
library(ggplot2)
library(annotate)
library(gage)
library(gageData)
library("AnnotationDbi")
library("org.Hs.eg.db")
columns(org.Hs.eg.db)

## MSigDB - Molecular Signatures Database
msigdb <- readList("/storage1/fs1/leyao.wang/Active/saran/RNA_Sara/h.all.v2022.1.Hs.entrez.gmt") 

# get gene lists
ifna <- msigdb$HALLMARK_INTERFERON_ALPHA_RESPONSE
ifng <- msigdb$HALLMARK_INTERFERON_GAMMA_RESPONSE

intersect(ifna,ifng)
ifna_only <- setdiff(ifna, ifng)
ifng_only <- setdiff(ifng,ifna)

df_ifna <- data.frame(IFNA = ifna_only)
df_ifng <- data.frame(IFNG = ifng_only)

# Add Annotations
df_ifna$symbol = mapIds(org.Hs.eg.db,
                      keys=df_ifna$IFNA, 
                      column="SYMBOL",
                      keytype="ENTREZID",
                      multiVals="first")
df_ifng$symbol = mapIds(org.Hs.eg.db,
                        keys=df_ifng$IFNG, 
                        column="SYMBOL",
                        keytype="ENTREZID",
                        multiVals="first")

# U-CELL CELL SCORING
library(UCell)
library(dplyr)

expr.matrix <- combined.all[["RNA"]]@counts
ifnA <- df_ifna$symbol
ifnG <- df_ifng$symbol
gene.sets <- list(Alpha_Signature = ifnA, Gamma_Signature = ifnG)
scores <- ScoreSignatures_UCell(expr.matrix, features = gene.sets)
head(scores)

scores_df <- data.frame(scores)
sum(scores_df$Alpha_Signature_UCell > 0.2)
sum(scores_df$Gamma_Signature_UCell > 0.2)
# No Cells have a score higher than 0.5 for either gene set
# 174 Cells have a score > 0.3 for Beta and only 1 cell > 0.3 for Gamma

########## ALL TOGETHER
alpha <- rep("Alpha", length(scores_df$Alpha_Signature_UCell))
gamma <- rep("Gamma", length(scores_df$Gamma_Signature_UCell))
scoring <- c(scores_df$Alpha_Signature_UCell, scores_df$Gamma_Signature_UCell)
samp <- c(alpha, gamma)
cells <- row.names(scores_df)
df <- data.frame(score = scoring, sample = samp, cell = c(cells, cells))

meta.data <- data.frame(combined.all@meta.data)
meta.data$cell <- row.names(meta.data) 
head(meta.data)
meta.data <- meta.data %>% select(cell, lin_main, celltype.stim.combined, main, celltype.stim, sample)
head(df)

df <- left_join(df, meta.data, by = "cell")
head(df)
df$SAMP <- paste(df$sample.x, df$sample.y, sep = "_")

# Violin PLot
ggplot(df, aes(x=SAMP, y=score)) + geom_violin(draw_quantiles = c(0.90))  +
  xlab("Score") + ylab("Sample") + 
  theme(axis.text.y=element_text(size=12,color='black', face = "bold"),
        axis.title.y= element_text(size=16,color='black',face = "bold"),
        axis.title.x= element_text(size=16,color='black',face = "bold"))  #+ geom_jitter(height = 0, width = 0.05)



# Get 99% Quantile
quantile(scores_df$Alpha_Signature_UCell, probs = 0.90)  # 0.1763917
quantile(scores_df$Gamma_Signature_UCell, probs = 0.90) # 0.1638814

# How many cells score over that 99% quantile expression
sum(scores_df$Alpha_Signature_UCell > 0.1763917)  # 98 Cells
sum(scores_df$Gamma_Signature_UCell > 0.1638814) # 99 Cells

sigB <- which(scores_df$Alpha_Signature_UCell > 0.1763917) 
sigG <- which(scores_df$Gamma_Signature_UCell > 0.1638814)
alphaa <- row.names(scores_df)[sigB]
gammaa <- row.names(scores_df)[sigG]


# Combined
both <- intersect(alphaa, gammaa)
Bonly <- setdiff(alphaa, gammaa)
Gonly <- setdiff(gammaa, alphaa)

combined.all$Ucell <- "neither"
combined.all@meta.data[both, "Ucell"] <- "Both"
combined.all@meta.data[Bonly, "Ucell"] <- "IFN-Alpha"
combined.all@meta.data[Gonly, "Ucell"] <- "IFN-Gamma"


### Dimplot
combined.all <- SetIdent(combined.all, value = "Ucell")
DimPlot(combined.all, label = F, pt.size = 1, order = c( "IFN-Alpha", "IFN-Gamma", "Both", "neither"), cols = c("grey", "red", "blue", "green"))

combined.all <- SetIdent(combined.all, value = "main")
DimPlot(combined.all, label = F, pt.size = 1) +  guides(color = guide_legend(override.aes = list(size=5), ncol=1))

table(combined.all$sample, combined.all$Ucell)

