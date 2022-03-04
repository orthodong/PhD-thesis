
#install.packages("devtools")

library(devtools)
#install_github("immunogenomics/harmony")
library(Rcpp)
#library(harmony)
library(dplyr)
library(Seurat)
library(patchwork)
library(celldex)
library(SingleR)
library(SingleCellExperiment)
setwd("C:/Users/cdong/Documents/Din-WT-analysis")
setwd("H:/Changchun/WT-analysis")
getwd()

# Load the datasets
ctrl1d.data <- Read10X("1dctrl")
cut1d.data <- Read10X("1dcut")
ctrl3d.data <- Read10X("3dctrl")
cut3d.data <- Read10X("3dcut")
WT.data <- Read10X("cd1") 

Ctrl1 <- CreateSeuratObject(counts = ctrl1d.data, project = "Ctrl1",names.field = 1, min.cells = 3, min.features = 200)
Ctrl3<- CreateSeuratObject(counts = ctrl3d.data, project = "Ctrl3",names.field = 1, min.cells = 3, min.features = 200)
Ctrl2<- CreateSeuratObject(counts = cut1d.data, project = "Ctrl2",names.field = 1, min.cells = 3, min.features = 200)
Ctrl4<- CreateSeuratObject(counts = cut3d.data, project = "Ctrl4",names.field = 1, min.cells = 3, min.features = 200)
WT<- CreateSeuratObject(counts = WT.data, project = "Ctrl5",min.cells = 3, min.features = 200)
WT
Ctrl1
Ctrl2
Ctrl3
Ctrl4

#add "stim" to meta.data 
WT[["replicate"]]<-c("wt-5")
Ctrl1[["replicate"]]<-c("wt-1")
Ctrl2[["replicate"]]<-c("wt-2")
Ctrl3[["replicate"]]<-c("wt-3")
Ctrl4[["replicate"]]<-c("wt-4")

din.list<-list(WT,Ctrl1,Ctrl2,Ctrl3,Ctrl4)

din.list[["percent.mt"]]<- PercentageFeatureSet(din.list, pattern = "^mt-")
din.list <- lapply(X =din.list, FUN = function(x) {
  x[["percent.mt"]]<- PercentageFeatureSet(x, pattern = "^mt-")
  })
WT[["percent.mt"]]<- PercentageFeatureSet(WT, pattern = "^mt-")
Ctrl1[["percent.mt"]]<- PercentageFeatureSet(Ctrl1, pattern = "^mt-")
Ctrl2[["percent.mt"]]<- PercentageFeatureSet(Ctrl2, pattern = "^mt-")
Ctrl3[["percent.mt"]]<- PercentageFeatureSet(Ctrl3, pattern = "^mt-")
Ctrl4[["percent.mt"]]<- PercentageFeatureSet(Ctrl4, pattern = "^mt-")
wt.list<-list(WT,Ctrl1,Ctrl2,Ctrl3,Ctrl4)


VlnPlot(WT, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(Ctrl1, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3)
VlnPlot(Ctrl2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(Ctrl3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(Ctrl4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

VlnPlot(Ctrl2, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
VlnPlot(din.list,features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)

wt.list <- lapply(X =wt.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features <- SelectIntegrationFeatures(object.list = wt.list)
wt.anchors <- FindIntegrationAnchors(object.list = wt.list, anchor.features = features)
# this command creates an 'integrated' data assay
wt.combined <- IntegrateData(anchorset = wt.anchors)
# specify that we will perform downstream analysis on the corrected data note that the original
# unmodified data still resides in the 'RNA' assay
DefaultAssay(wt.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
wt.combined <- ScaleData(wt.combined, verbose = T)
wt.combined <- RunPCA(wt.combined, npcs = 30, verbose = FALSE)
wt.combined <- RunUMAP(wt.combined, reduction = "pca", dims = 1:5)
wt.combined <- FindNeighbors(wt.combined, reduction = "pca", dims = 1:5)
wt.combined <- FindClusters(wt.combined, resolution = 0.5)
# Visualization
p1 <- DimPlot(wt.combined, reduction = "umap", group.by ="replicate")
p2 <- DimPlot(wt.combined, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2
p1
p2
###########################################################################################
#combine  4 4w samples together and intend to compare 2w and 4w

wt4w<-CreateSeuratObject(counts = cbind(ctrl1d.data,ctrl3d.data,cut1d.data,cut3d.data), project = "wt4w",names.field = 1, min.cells = 3, min.features = 200)
str(wt4w)
wt4w
wt4w<- CreateSeuratObject(counts = cbind(ctrl1d.data,ctrl3d.data,cut1d.data,cut3d.data), project = "wt4w", min.cells = 5,min.features = 200) %>%
  Seurat::NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData(verbose = FALSE) %>% 
  RunPCA(pc.genes = wt4w@var.genes, npcs = 20, verbose = FALSE)
#Make sure that the dataset ID is in the object's metadata. Here, we define datasets with the variable stim
wt4w@meta.data$stim <- c(rep("ctrl1d", ncol(ctrl1d.data)), rep("ctrl3d", ncol(ctrl3d.data)),rep("cut1d", ncol(cut1d.data)),rep("cut3d", ncol(cut3d.data)))

VlnPlot(object =wt.combined, features = "PC_1", group.by = "replicate",  pt.size = .1)
wt4w[["Time"]]<-c("4w")
wt2w<-WT
wt2w[["Time"]]<-c("2w")
wt.list<-list(wt2w,wt4w)

wt.list <- lapply(X =wt.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features <- SelectIntegrationFeatures(object.list = wt.list)
wt.anchors <- FindIntegrationAnchors(object.list = wt.list, anchor.features = features)
# this command creates an 'integrated' data assay
wt.combined <- IntegrateData(anchorset = wt.anchors)
# specify that we will perform downstream analysis on the corrected data note that the original
# unmodified data still resides in the 'RNA' assay
DefaultAssay(wt.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
wt.combined <- ScaleData(wt.combined, verbose = T)
wt.combined <- RunPCA(wt.combined, npcs = 30, verbose = FALSE)
wt.combined <- RunUMAP(wt.combined, reduction = "pca", dims = 1:5)
wt.combined <- FindNeighbors(wt.combined, reduction = "pca", dims = 1:5)
wt.combined <- FindClusters(wt.combined, resolution = 0.5)
# Visualization
p1 <- DimPlot(wt.combined, reduction = "umap", group.by ="Time")
p2 <- DimPlot(wt.combined, reduction = "umap", label = TRUE,pt.size = 1)
p1 + p2
p1
p2
DimPlot(wt.combined, reduction = "umap",shape.by ="Time",pt.size = 2,split.by = "Time")
wt.sum1<-table(wt.combined$Time,wt.combined$seurat_clusters)
write.csv(wt.sum1,"wt cluster information.csv")
ElbowPlot(wt.combined)
saveRDS(wt.combined, file = "wt-combined0418.rds")


wt.markers <- FindAllMarkers(wt.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
wt.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
top10 <- wt.markers %>% group_by(cluster) %>% top_n(n =5, wt= avg_log2FC)
DoHeatmap(wt.combined, features = top10$gene) + NoLegend()
DotPlot(wt.combined, features = top10$gene) + NoLegend()

top20 <- wt.markers %>% group_by(cluster) %>% top_n(n =20, wt= avg_log2FC)
table(top20)
write.csv(top20,"wt cluster top20.csv")

library(celldex)
library(SingleR)
library(SingleCellExperiment)

mouseref <- MouseRNAseqData()

wt.qury<-wt.combined
wt.qury<-as.SingleCellExperiment(wt.qury)
library(SingleR)
predictions <- SingleR(test=wt.qury, assay.type.test=1, 
                       ref=mouseref, labels=mouseref$label.main)

#
wt.combined[["lables"]] <- predictions$labels
wt.annotation.sum<-table(wt.combined$lables)
wt.annotation.sum<-table(wt.combined$lables,wt.combined$seurat_clusters)
write.csv(wt.annotation.sum,"wt-annotation sum.csv")

DimPlot(wt.combined, reduction = "umap")
DimPlot(wt.combined, reduction = "umap", group.by = "lables", label = TRUE, repel = TRUE)
DimPlot(wt.combined, reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE,pt.size = 1.5)
# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster.markers <- FindMarkers(wt.combined, ident.1 = 4, ident.2 = c(7), min.pct = 0.25)
cluster_top10 <- cluster.markers %>% top_n(n =50, wt= avg_log2FC)
write.csv(cluster.markers, "cluster4 vs 7.csv")
#

#####give names to clusters
new.cluster.ids <- c("Dental pulp-1","Dental pulp-2","Dental pulp-3","Macrophages-1",
                     "Ameloblasts-1","Ameloblasts-2","B cells","Ameloblasts-3","Granulocytes","Erythrocytes-1",
                     "Endothelial cells-1","Odontoblasts","Eryghrocytes-2","Endothelial cells-2","Macrophages-2","Dental pulp-4",
                     "T cell/Monocytes","Endothelial cells-3","Eryghrocytes-3")
names(new.cluster.ids) <- levels(wt.combined)
wt.combined<- RenameIdents(wt.combined, new.cluster.ids)
DimPlot(wt.combined, reduction = "umap", label = TRUE, pt.size = 1.5) + NoLegend()
DimPlot(wt.combined, reduction = "umap", label = TRUE, pt.size = 1.5) +ggsci::scale_color_igv()
wt.combined[['label']]<-wt.combined@active.ident
saveRDS(wt.combined,'wt.combined.02092022.rds')


# you can plot raw counts as well
VlnPlot(wt.combined, features = c("Ambn", "Enam"), slot = "counts", log = TRUE)+
VlnPlot(wt.combined, features = c("Ambn", "Enam"), slot = "counts", log = TRUE)

library(reshape2)
vln.df<-as.data.frame(wt.combined[["RNA"]]@data[dimnames()])
vln.df$gene<-rownames(vln.df)
vln.df=melt(vln,df.id="gene")
colnames(vln.df)[c(2,3)]=c("CB","exp")
genes<-c("Ambn", "Enam")
VlnPlot(wt.combined, features = genes, slot = "counts", log = TRUE)+facet_grid(genes)

#stacking vollin plots
# devtools::install_github("road2stat/ggsci")
library(ggsci)

DefaultAssay(wt.combined) <- "RNA"
markers<-c("Bglap2",'Msx2',"Krt17","Msx2","Sox9","Ifitm1","Sfrp2","Igfbp3","Wnt10a","Gli1","Dspp","Mmp20","Krt14","Enam","Mki67")
markers<-CaseMatch(markers,rownames(wt.combined))
VlnPlot(wt.combined,features = markers, pt.size = 2,stack = T)+ggsci::scale_fill_jco()
VlnPlot(wt.combined,features = markers, pt.size = 2,stack = T)+NoLegend()

##########arrange orders displaying
wt.combined[["cluster_annotation"]] <- wt.combined@active.ident
group_list<-factor(wt.combined@meta.data$cluster_annotation,levels =c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18),
                   labels=c("Dental pulp-1","Dental pulp-2","Dental pulp-3",
                           "Dental pulp-4","Odontoblasts","Ameloblasts-1","Ameloblasts-2","Ameloblasts-3","B cells",
                           "Granulocytes","Erythrocytes-1","Eryghrocytes-2","Endothelial cells-1","Endothelial cells-2",
                           "Endothelial cells-3","Macrophages-1","Macrophages-2","Eryghrocytes-3","T cell/Monocytes"))
wt.combined@meta.data$cluster_annotation<-group_list
VlnPlot(wt.combined,features = markers, pt.size = 0,stack =T)+ggsci::scale_fill_jco()
VlnPlot(wt.combined,features = markers, pt.size = 0)
wt.combined@meta.data$lables<-as.character(wt.combined@meta.data$lables)
Idents(wt.combined)<-wt.combined$seurat_clusters
#not successfull above
#########
FeaturePlot(wt.combined, features = c("4930453N24Rik"))


wt_pulp_ep<-subset(wt.combined,idents=c("Dental pulp-1","Dental pulp-2","Dental pulp-3","Dental pulp-4","Odontoblasts","Ameloblasts-1","Ameloblasts-2","Ameloblasts-3"))
wt_pulp_ep
DimPlot(wt_pulp_ep,pt.size = 1.5)

# wt_pulp_ep <- ScaleData(wt_pulp_ep , verbose = T)
# wt_pulp_ep <- RunPCA(wt_pulp_ep , npcs = 30, verbose = FALSE)
wt_pulp_ep <- RunUMAP(wt_pulp_ep , reduction = "pca", dims = 1:5)
wt_pulp_ep <- FindNeighbors(wt_pulp_ep , reduction = "pca", dims = 1:5)
wt_pulp_ep<- FindClusters(wt_pulp_ep , resolution = 0.5)
DimPlot(wt_pulp_ep,pt.size = 1.5)+ggsci::scale_color_igv()
DimPlot(wt_pulp_ep,pt.size = 1.5,group.by = 'label',label = T)+ggsci::scale_color_igv()
saveRDS(wt_pulp_ep,'wt_pulp_ep.rds')


wt.markers <- FindAllMarkers(wt_pulp_ep, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
wt.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
top10 <- wt.markers %>% group_by(cluster) %>% top_n(n =5, wt= avg_log2FC)
DoHeatmap(wt_pulp_ep, features = top10$gene) + NoLegend()
top20 <- wt.markers %>% group_by(cluster) %>% top_n(n =20, wt= avg_log2FC)
table(top20)
write.csv(top20,"wt cluster top20.csv")


#BiocManager::install("monocle")
library(monocle)
#construct cds by DIY, don't need this, not run
mono_matrix<-as(as.matrix(GetAssayData(wt_pulp_ep,slot = "counts")),"sparseMatrix")
feature_ann<-data.frame(gene_id=rownames(mono_matrix),gene_short_name=rownames(mono_matrix))
rownames(feature_ann)<-rownames(mono_matrix)
mono_fd<-new("AnnotatedDataFrame",data=feature_ann)
sample_ann<-wt_pulp_ep@meta.data
mono_pd<-new("AnnotatedDataFrame",data=sample_ann)
#use this directly
wt_cds<-as.CellDataSet(wt_pulp_ep,assay = "RNA")
wt_cds<-estimateSizeFactors(wt_cds)
wt_cds<-estimateDispersions(wt_cds)
disp_table<-dispersionTable(wt_cds)
unsup_clustering_genes<-subset(disp_table,mean_expression>=0.1)
wt_cds<-setOrderingFilter(wt_cds,unsup_clustering_genes$gene_id)

wt_cds<-reduceDimension(wt_cds,max_components = 2,method="DDRTree")

wt_cds<-orderCells(wt_cds)
head
saveRDS(wt_cds, file = "wt-cds-after DDRTree.rds")

plot_cell_trajectory(wt_cds,cell_size=1)
plot_cell_trajectory(wt_cds,cell_size=1,color_by = "Pseudotime",lables="cluster_annotation")
plot_cell_trajectory(wt_cds,cell_size=1,color_by = "seurat_clusters",lables="cluster_annotation")
plot_cell_trajectory(wt_cds,cell_size=1,color_by = "lables",lables="cluster_annotation")


my_genes <- row.names(subset(fData(wt_cds),
                             gene_short_name %in% c("Sox9", "Dspp","Msx2",'Sfrp2','Enam','Igfbp3','Igfbp5','Mki67','4930453N24Rik')))
cds_subset <- wt_cds[my_genes,]
plot_genes_in_pseudotime(cds_subset, color_by = "State")




wt_expressed_genes <-  row.names(subset(fData(wt_cds),
                                          wt_cds$nFeature_RNA >= 10))
wt_filtered <- wt_cds[wt_expressed_genes,]
my_genes <- row.names(subset(fData(wt_filtered),
                             gene_short_name %in% c("Dspp","Gli1","Enam","Amelx","Sox9","4930453N24Rik")))
my_genes <- row.names(subset(fData(wt_filtered),
                             gene_short_name %in% c("4930453N24Rik")))
cds_subset <- wt_filtered[my_genes,]
plot_genes_in_pseudotime(cds_subset, color_by = "State",cell_size = 2)




marker_genes <- row.names(subset(fData(wt_cds),
                                 gene_short_name %in% c("Dspp","Gli1","Enam","Amelx","Sox9","4930453N24Rik")))
marker_genes <- row.names(subset(fData(wt_cds),))
diff_test_res <- differentialGeneTest(wt_cds[marker_genes,],
                                      fullModelFormulaStr = "~sm.ns(Pseudotime)")
sig_gene_names <- row.names(subset(diff_test_res, qval < 0.1))
plot_pseudotime_heatmap(wt_cds[sig_gene_names,],
                        num_clusters = 3,
                        cores = 1,
                        show_rownames = F)
p<-plot_pseudotime_heatmap(wt_cds[sig_gene_names,],
                           num_clusters = 15,
                           cores = 1,
                           show_rownames = T,
                           return_heatmap = T)
p
#stop############################
str(p)
p$tree_row
clusters<-cutree(p$tree_row,k=15)
clustering<-data.frame(clusters)
clustering[,1]<-as.character(clustering[,1])
colnames(clustering)<-"Gene_Clusters"
table(clustering)
write.csv(clustering,"clustering1.csv")
#stop here

write.csv(clustering,"clustering.csv")

#############################dental pulp cluster
wt_pulp<-subset(wt.combined,idents=c("Dental pulp-1","Dental pulp-2","Dental pulp-3","Dental pulp-4","Odontoblasts"))
wt_pulp
DimPlot(wt_pulp,pt.size = 1.5)

wt_pulp <- ScaleData(wt_pulp , verbose = T)
wt_pulp<- RunPCA(wt_pulp, npcs = 30, verbose = FALSE)
wt_pulp<- RunUMAP(wt_pulp , reduction = "pca", dims = 1:5)
wt_pulp<- FindNeighbors(wt_pulp, reduction = "pca", dims = 1:5)
wt_pulp<- FindClusters(wt_pulp, resolution = 0.5)
DimPlot(wt_pulp,pt.size = 1.5)

wt.markers <- FindAllMarkers(wt_pulp, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
wt.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
top10 <- wt.markers %>% group_by(cluster) %>% top_n(n =5, wt= avg_log2FC)
DoHeatmap(wt_pulp, features = top10$gene) + NoLegend()
top20 <- wt.markers %>% group_by(cluster) %>% top_n(n =20, wt= avg_log2FC)
table(top20)
write.csv(top20,"wt cluster top20.csv")


library(monocle)
#construct cds by DIY, don't need this, not run
mono_matrix<-as(as.matrix(GetAssayData(wt_pulp,slot = "counts")),"sparseMatrix")
feature_ann<-data.frame(gene_id=rownames(mono_matrix),gene_short_name=rownames(mono_matrix))
rownames(feature_ann)<-rownames(mono_matrix)
mono_fd<-new("AnnotatedDataFrame",data=feature_ann)
sample_ann<-wt_pulp@meta.data
mono_pd<-new("AnnotatedDataFrame",data=sample_ann)
#use this directly
wt_pulp_cds<-as.CellDataSet(wt_pulp)
wt_pulp_cds<-estimateSizeFactors(wt_pulp_cds)
wt_pulp_cds<-estimateDispersions(wt_pulp_cds)
disp_table<-dispersionTable(wt_pulp_cds)
unsup_clustering_genes<-subset(disp_table,mean_expression>=0.1)
wt_pulp_cds<-setOrderingFilter(wt_pulp_cds,unsup_clustering_genes$gene_id)

wt_pulp_cds<-reduceDimension(wt_pulp_cds,max_components = 2,method="DDRTree")
wt_pulp_cds<-orderCells(wt_pulp_cds)
head
plot_cell_trajectory(wt_pulp_cds,cell_size=1)
plot_cell_trajectory(wt_pulp_cds,cell_size=1,color_by = "Pseudotime",lables="cluster_annotation")
plot_cell_trajectory(wt_pulp_cds,cell_size=1,color_by = "seurat_clusters",lables="cluster_annotation")


wt_expressed_genes <-  row.names(subset(fData(wt_pulp_cds),
                                        wt_pulp_cds$nFeature_RNA >= 10))
wt_filtered <- wt_pulp_cds[wt_expressed_genes,]
my_genes <- row.names(subset(fData(wt_filtered),
                             gene_short_name %in% c("Dspp","Gli1","Enam","Amelx","Sox9","4930453N24Rik")))
my_genes <- row.names(subset(fData(wt_filtered),
                             gene_short_name %in% c("Wnt7a")))
cds_subset <- wt_filtered[my_genes,]
plot_genes_in_pseudotime(cds_subset, color_by = "State",cell_size = 2)




marker_genes <- row.names(subset(fData(wt_pulp_cds),
                                 gene_short_name %in% c("Dspp","Gli1","Enam","Amelx","Sox9","4930453N24Rik")))
marker_genes <- row.names(subset(fData(wt_pulp_cds),))
diff_test_res <- differentialGeneTest(wt_pulp_cds[marker_genes,],
                                      fullModelFormulaStr = "~sm.ns(Pseudotime)")
sig_gene_names <- row.names(subset(diff_test_res, qval < 0.1))
plot_pseudotime_heatmap(wt_pulp_cds[sig_gene_names,],
                        num_clusters = 3,
                        cores = 1,
                        show_rownames = F)
p<-plot_pseudotime_heatmap(wt_pulp_cds[sig_gene_names,],
                           num_clusters = 15,
                           cores = 1,
                           show_rownames = T,
                           return_heatmap = T)
p$tree_row
clusters<-cutree(p$tree_row,k=15)
clustering<-data.frame(clusters)
clustering[,1]<-as.character(clustering[,1])
colnames(clustering)<-"Gene_Clusters"
table(clustering)
write.csv(clustering,"pulp_clustering.csv")





## 02/09/2022

library(Nebulosa)
DefaultAssay(wt_pulp_ep)
DefaultAssay(wt_pulp_ep)<-"RNA"
plot_density(wt_pulp_ep,'4930453N24Rik')
plot_density(wt_pulp_ep,c('4930453N24Rik','Gli1'),joint = TRUE)
plot_density(wt_pulp_ep,c('4930453N24Rik','Axin2'),joint = TRUE)
