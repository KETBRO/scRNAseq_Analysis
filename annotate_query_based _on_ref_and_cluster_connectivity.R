library(Seurat)
seurat_ACO1 <- readRDS("seurat_ACO1.rds")
seurat_ref <- readRDS("ref_seurat_obj.rds")

# library(patchwork)
# plot1 <- UMAPPlot(seurat_ref, group.by="branch")
# plot2 <- UMAPPlot(seurat_ref, group.by="celltype")
# plot3 <- FeaturePlot(seurat_ref, c("SOX2","DCX","FOXG1"),
#                      ncol=3, pt.size = 0.1)
# ((plot1 / plot2) | plot3) + plot_layout(width = c(1,3))

#Transcriptome similarity on cell cluster level
#calculate the average transcriptome profiles for
# every annotated cell type in the reference data set 
# and every cell cluster in the query data set
options(Seurat.object.assay.version = "v3")
options(Seurat.object.assay.version = "v5")

seurat_ACO1[["percent.mt"]] <- PercentageFeatureSet(seurat_ACO1, pattern = "^MT[-\\.]")
seurat_ACO1 <- subset(seurat_ACO1, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 5)
seurat_ACO1 <- NormalizeData(seurat_ACO1)

# avg_expr_ref <- sapply(sort(unique(seurat_ref$celltype)), function(ct) rowMeans(seurat_ref@assays$RNA@data[,which(seurat_ref$celltype == ct)] ))
# avg_expr_ACO1 <- sapply(levels(seurat_ACO1@active.ident), function(ct) rowMeans(seurat_ACO1@assays$RNA@data[,which(seurat_ACO1@active.ident == ct)]))


# 1. Scaling
seurat_ACO1 <- ScaleData(seurat_ACO1)
seurat_ACO1 <- ScaleData(seurat_ACO1, vars.to.regress = c("nFeature_RNA", "percent.mt"))

# 2. Dimensionality Reduction
seurat_ACO1 <- FindVariableFeatures(seurat_ACO1, nfeatures = 3000)
seurat_ACO1 <- RunPCA(seurat_ACO1, npcs = 50)
seurat_ACO1 <- RunUMAP(seurat_ACO1, dims = 1:20)

# 3. Clustering
seurat_ACO1 <- FindNeighbors(seurat_ACO1, dims = 1:20)
seurat_ACO1 <- FindClusters(seurat_ACO1, resolution = 1) # Adjust resolution for more/fewer clusters

# 4. Cluster Annotation
library(dplyr)
cl_markers <- FindAllMarkers(seurat_ACO1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = log(1.2))
top_markers <- cl_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

# 5. Rename Identities (customize based on your markers!)
new_ident <- setNames(c("Dorsal telen. NPC",
                        "Midbrain-hindbrain boundary neuron",
                        "Dorsal telen. neuron",
                        "Dien. and midbrain excitatory neuron",
                        "MGE-like neuron","G2M dorsal telen. NPC",
                        "Dorsal telen. IP","Dien. and midbrain NPC",
                        "Dien. and midbrain IP and excitatory early neuron",
                        "G2M Dien. and midbrain NPC",
                        "G2M dorsal telen. NPC",
                        "Dien. and midbrain inhibitory neuron",
                        "Dien. and midbrain IP and early inhibitory neuron",
                        "Ventral telen. neuron",
                        "Unknown 1",
                        "Unknown 2"),
                      levels(seurat_ACO1))
seurat_ACO1 <- RenameIdents(seurat_ACO1, new_ident)

avg_expr_ref <- sapply(sort(unique(seurat_ref$celltype)), function(ct) {rowMeans(GetAssayData(seurat_ref, assay = "RNA", slot = "data")[, which(seurat_ref$celltype == ct)])})
avg_expr_ACO1 <- sapply(levels(seurat_ACO1@active.ident), function(ct) {rowMeans(GetAssayData(seurat_ACO1, assay = "RNA", slot = "data")[, which(seurat_ACO1@active.ident == ct)])})

genes2cor <- intersect(VariableFeatures(seurat_ref), rownames(seurat_ACO1))
corr2ref_cl <- cor(avg_expr_ACO1[genes2cor,], avg_expr_ref[genes2cor,], method="spearman")


library(gplots)
heatmap.2(corr2ref_cl, scale="none", trace="none", key=F, keysize=0.5, margins=c(15,17),
          labRow = colnames(avg_expr_ACO1), labCol = colnames(avg_expr_ref), cexRow=0.8, cexCol=0.8,
          col=colorRampPalette(rev(c("#b2182b","#d6604d","#f4a582","#fddbc7","#f7f7f7","#d1e5f0","#92c5de","#4393c3","#2166ac")))(30))


ranked_expr_ref <- apply(avg_expr_ref[genes2cor,],2,rank)
library(presto)
ranked_expr_ACO1 <- rank_matrix(avg_expr_ACO1@assays$RNA@data[genes2cor,])$X_ranked

rank_matrix <- function(mat, assay = "RNA", layer = "data") {
  # If input is a Seurat object
  if (inherits(mat, "Seurat")) {
    mat <- SeuratObject::LayerData(mat, assay = assay, layer = layer)
  }
  
  # Check if we got a valid matrix
  if (!(is.matrix(mat) || is.data.frame(mat) || inherits(mat, "dgCMatrix"))) {
    stop("Input must be a matrix, data.frame, dgCMatrix, or Seurat object")
  }
  
  if (is.matrix(mat) || is.data.frame(mat)) {
    ranked_mat <- apply(mat, 2, rank, ties.method = "average")
  } else if (inherits(mat, "dgCMatrix")) {
    df_mat <- Matrix::summary(mat)
    dfs_mat <- split(df_mat, df_mat$j)
    df_mat_ranked <- do.call(rbind, lapply(dfs_mat, function(df) {
      num_zeros <- nrow(mat) - nrow(df)
      ranks_nonzero <- rank(df$x, ties.method = "average")
      df$x <- ranks_nonzero + num_zeros - (1 + num_zeros)/2
      return(df)
    }))
    ranked_mat <- Matrix::sparseMatrix(
      i = df_mat_ranked$i, 
      j = df_mat_ranked$j, 
      x = df_mat_ranked$x, 
      dims = dim(mat), 
      dimnames = dimnames(mat)
    )
  }
  return(ranked_mat)
}

# Usage:
ranked_expr_ds1 <- rank_matrix(seurat_ACO1[genes2cor, ], assay = "RNA", layer = "data")

library(qlcMatrix)
corr2ref_cell <- corSparse(ranked_expr_ds1, ranked_expr_ref)
ct_maxcor <- colnames(avg_expr_ref)[apply(corr2ref_cell, 1, which.max)]
seurat_ACO1$celltype_maxcor <- ct_maxcor

plot1 <- UMAPPlot(seurat_ACO1, label=T)
plot2 <- UMAPPlot(seurat_ACO1, group.by="celltype_maxcor", label=T)
plot1 | plot2

#--------------Seurat-based label transfer------------------
anchors <- FindTransferAnchors(reference = seurat_ref, query = seurat_ACO1, dims = 1:30, npcs = 30)
predictions <- TransferData(anchorset = anchors, refdata = seurat_ref$celltype, dims = 1:30)
seurat_ACO1$celltype_transfer <- predictions$predicted.id

seurat_ref <- UpdateSeuratObject(seurat_ref)
seurat_ACO1 <- UpdateSeuratObject(seurat_ACO1)

plot1 <- UMAPPlot(seurat_ACO1, label=T)
plot2 <- UMAPPlot(seurat_ACO1, group.by="celltype_transfer", label=T)
plot1 | plot2

#---------------Cluster connectivity analysis with PAGA----------------------------------------
library(Matrix)
install.packages("devtools")
devtools::install_github(repo = "mojaveazure/loomR", ref = "develop")
library(loomR)
cell_attrs <- list(pca = Embeddings(seurat_ACO1,"pca")[,1:20],
                   umap = Embeddings(seurat_ACO1,"umap"),
                   celltype = seurat_ACO1@active.ident)
loom <- loomR::create("loom_obj.loom",
                      data = seurat_ACO1[['RNA']]$data,
                      layers = list(counts = seurat[['RNA']]$counts),
                      cell.attrs = cell_attrs)
loom$close_all()
saveRDS(seurat_ACO1, file="seurat_ACO1.rds")
seurat_ACO1.rds
seurat_ACO1 <- readRDS("seurat_ACO1.rds")
library(anndata)
library(Matrix)
adata <- AnnData(X = t(seurat_ACO1[['RNA']]$data),
                 obs = data.frame(celltype = seurat_ACO1@active.ident, row.names = colnames(seurat_ACO1)),
                 var = seurat_ACO1[['RNA']]@meta.data,
                 layers = list(counts = t(seurat_ACO1[['RNA']]$counts)),
                 obsm = list(pca = Embeddings(seurat_ACO1,"pca")[,1:20],
                             umap = Embeddings(seurat_ACO1,"umap"))
)
adata$write_h5ad("anndata_obj.h5ad")

install.packages("reticulate")
library("reticulate")
py_install("scanpy", pip=T)
# options(reticulate.conda_binary = "C:/Users/Krishnendu/anaconda3/Scripts/conda.exe")
# conda_path <- reticulate::conda_binary()
# conda_create("scanpy_env")  # Create a new environment
# use_condaenv("scanpy_env")  # Activate it
# 
# # Install compatible versions
# py_install(c("numpy==2.2.0", "numba", "scanpy"), pip = TRUE)

# Test
sc <- import("scanpy")


adata_ACO1 <- sc$read("anndata_obj.h5ad") 

sc$pp$neighbors(adata_ACO1, n_neighbors=20L, use_rep='pca')
reticulate::py_install("python-igraph", pip = TRUE)
sc$tl$paga(adata_ACO1, groups='celltype')
adata_ACO1$write_h5ad("anndata_obj.h5ad")
plt <- import("matplotlib")
plt$use("Agg", force = TRUE)
sc$pl$paga(adata_ACO1,
           color='celltype',
           fontsize=7,
           frameon=FALSE,
           save="ACO1_paga.png")






