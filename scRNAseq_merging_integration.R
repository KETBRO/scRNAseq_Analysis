library(Seurat)
library(dplyr)
library(patchwork)


seurat_DS1 <- readRDS("seurat_ACO3.rds")
seurat_DS2 <- readRDS("seurat_ACO1.rds")


seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^MT[-\\.]")
seurat <- subset(seurat, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 5)
seurat <- NormalizeData(seurat)

seurat <- seurat%>%
  FindVariableFeatures(nfeatures = 3000) %>%
  ScaleData() %>%
  RunPCA(npcs = 50) %>%
  RunUMAP(dims = 1:20)

plot1 <- DimPlot(seurat, group.by="orig.ident")
plot2 <- FeaturePlot(seurat, c("FOXG1","EMX1","DLX2","LHX9"), ncol=2, pt.size = 0.1)
plot1 + plot2 + plot_layout(widths = c(1.5, 2))

seurat_DS1 <- NormalizeData(seurat_DS1) %>% FindVariableFeatures(nfeatures = 3000)
seurat_DS2 <- NormalizeData(seurat_DS2) %>% FindVariableFeatures(nfeatures = 3000)

seurat_objs <- list(DS1 = seurat_DS1, DS2 = seurat_DS2)
anchors <- FindIntegrationAnchors(object.list = seurat_objs, dims = 1:30)
seurat <- IntegrateData(anchors, dims = 1:30)

seurat <- ScaleData(seurat)
seurat <- RunPCA(seurat, npcs = 50)
seurat <- RunUMAP(seurat, dims = 1:20)
seurat <- FindNeighbors(seurat, dims = 1:20) %>% FindClusters(resolution = 0.6)

plot1 <- UMAPPlot(seurat, group.by="orig.ident")
plot2 <- UMAPPlot(seurat, label = T)
plot3 <- FeaturePlot(seurat, c("FOXG1","EMX1","DLX2","LHX9"), ncol=2, pt.size = 0.1)
((plot1 / plot2) | plot3) + plot_layout(width = c(1,2))

saveRDS(seurat, file="integrated_seurat.rds")

DefaultAssay(seurat) <- "RNA"
plot1 <- UMAPPlot(seurat, group.by="orig.ident")
plot2 <- UMAPPlot(seurat, label = T)
plot3 <- FeaturePlot(seurat, c("FOXG1","EMX1","DLX2","LHX9"), ncol=2, pt.size = 0.1)
((plot1 / plot2) | plot3) + plot_layout(width = c(1,2))


seurat_harmony <- merge(seurat_DS1, seurat_DS2) %>%
  FindVariableFeatures(nfeatures = 3000) %>%
  ScaleData() %>%
  RunPCA(npcs = 50)
library(harmony)
seurat_harmony <- RunHarmony(seurat_harmony, group.by.vars = "orig.ident", dims.use = 1:20, max.iter.harmony = 50)
seurat_harmony <- RunUMAP(seurat_harmony, reduction = "harmony", dims = 1:20)
seurat_harmony <- FindNeighbors(seurat_harmony, reduction = "harmony", dims = 1:20) %>% FindClusters(resolution = 0.6)

saveRDS(seurat_harmony, file="integrated_harmony.rds")

plot1 <- UMAPPlot(seurat_harmony, group.by="orig.ident")
plot2 <- UMAPPlot(seurat_harmony, label = T)
plot3 <- FeaturePlot(seurat_harmony, c("FOXG1","EMX1","DLX2","LHX9"), ncol=2, pt.size = 0.1)
((plot1 / plot2) | plot3) + plot_layout(width = c(1,2))

