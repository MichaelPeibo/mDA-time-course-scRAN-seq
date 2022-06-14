

### conda activate urd ###

library(Seurat)
library(dplyr)
library(future)
options(future.globals.maxSize = 4000 * 1024^2)

#### neuron only comparison ####
setwd('/home/xupb/scRNA_data/cell2016')
load('manno_seu.Robj')
manno_seu=UpdateSeuratObject(manno_seu)
manno_seu@meta.data$dataset <- 'cell2016'
table(manno_seu@meta.data$celltype)
Idents(manno_seu)=manno_seu@meta.data$celltype

manno_neu=subset(manno_seu, idents=c('hDA0','hDA1','hDA2','hGaba','hNbGaba','hNbM','hNbML1',
												'hNbML5','hOMTN','hRN','hSert'))


load('/home/xupb/scRNA_data/time_course/manuscript/seurat_pipe_results/woday51/merged_seurat_pipe.Robj')
mda=mda.merge
mda@meta.data$dataset <- 'this_project'
meta=read.csv('/scanpy/woday51_scanpy_meta.csv',row.names=1)
Idents(mda)=meta$cluster
mda.neu=subset(mda,idents=c('11-N_DA','8-N_Glut','18-N_GABA',
                                    '14-N_DA_Neuroblast','17-N_Sero','16-N_Interneuron_Neuroblast',
                                    '12-N_Motor','15-N_Glut_Neuroblast'))
table(Idents(mda.neu))
ob.list <- list(manno_neu, mda.neu)
mda.features <- SelectIntegrationFeatures(object.list = ob.list, nfeatures = 2000)

mda.anchors <- FindIntegrationAnchors(object.list = ob.list, 
    anchor.features = mda.features, verbose = FALSE)

mda.integrated <- IntegrateData(anchorset = mda.anchors,
    verbose = FALSE)

DefaultAssay(mda.integrated) <- "integrated"
mda.integrated <- ScaleData(mda.integrated, verbose = FALSE)
mda.integrated <- RunPCA(mda.integrated, verbose = FALSE,npcs = 30)

mda.integrated <- FindNeighbors(object = mda.integrated, reduction = "pca", dims = 1:30)
mda.integrated <- FindClusters(mda.integrated, resolution =0.4)

table(Idents(object = mda.integrated),mda.integrated@meta.data$dataset)

save(mda.integrated,file='integrated_manno_mda_neuron_only_woday51.Robj')
load('integrated_manno_mda_neuron_only_woday51.Robj')

### try metaneighbour ###

var.gene=VariableFeatures(object = mda.integrated)

combined_mat=as.matrix(mda.integrated@assays$integrated@scale.data)
mda.neu@meta.data$celltype=Idents(mda.neu)
Study_ID = rep(c('1', '2'), c(ncol(manno_neu@assays$RNA@data), ncol(mda.neu@assays$RNA@data)))
Celltype = c(as.character(manno_neu@meta.data$celltype),as.character(mda.neu@meta.data$celltype))


### conda activate urd ###
library(MetaNeighbor) ### packageVersion('MetaNeighbor') = 1.6.0
library(SummarizedExperiment)
dat=SummarizedExperiment(assays=list(counts=combined_mat))

celltype_NV=MetaNeighbor::MetaNeighborUS(var_genes = var.gene,
dat = dat,
study_id = Study_ID,
cell_type = Celltype,fast_version=T)

library(gplots)
library(RColorBrewer)
cols = rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdYlBu"))(100))
breaks = seq(0, 1, length=101)

pdf('heatmap_ neuron only manno our woday51 .pdf')
gplots::heatmap.2(celltype_NV,
margins=c(15,15),
keysize=1,
key.xlab="AUROC",
key.title="NULL",
trace = "none",
density.info =c("none"),
col = cols,
RowSideColors=c(rep("darkgreen", 11), rep('deeppink',8)),
breaks = breaks,
dendrogram="row",     # only draw a row dendrogram
offsetRow=0.1,
offsetCol=0.1,
cexRow = 1.5,
cexCol = 1.5)

par(lend = 1)      # square line ends for the color legend
legend("top",      # location of the legend on the heatmap plot
    legend = c("1-La Manno el al.2016", "2-This project"), # category labels
    col = c("darkgreen", "deeppink"),  # color key
    lty= 1,             # line style
    lwd = 10            # line width
)
dev.off()

library(pheatmap)

ann_row=data.frame(dataset=c(rep("cell2016", 11), rep('This_study',8)))
rownames(ann_row)=rownames(celltype_NV)

pdf('pheatmap_ neuron only manno our.pdf',width=10)
pheatmap(celltype_NV, 
         color = cols, cutree_rows=5,cutree_cols=5,
         breaks = breaks, cellwidth = 12, cellheight = 12,
         cluster_rows = TRUE,cluster_cols = TRUE, treeheight_col = 0,
         clustering_method = "ward.D2",
         annotation_row = ann_row,
         annotation_colors = list(dataset = c(cell2016 = "darkgreen", This_study = "deeppink"))
         )
dev.off()

#### progenitor only comparison ####
manno_prog=subset(manno_seu, idents=c('hDA0','hDA1','hDA2','hGaba','hNbGaba','hNbM','hNbML1',
												'hNbML5','hOMTN','hRN','hSert'),invert=T)

mda.prog=subset(mda,idents=c('11-N_DA','8-N_Glut','18-N_GABA',
                                    '14-N_DA_Neuroblast','17-N_Sero','16-N_Interneuron_Neuroblast',
                                    '12-N_Motor','15-N_Glut_Neuroblast'),invert=T)

ob.list <- list(manno_prog, mda.prog)
mda.features <- SelectIntegrationFeatures(object.list = ob.list, nfeatures = 2000)

mda.anchors <- FindIntegrationAnchors(object.list = ob.list, 
    anchor.features = mda.features, verbose = FALSE)

mda.integrated <- IntegrateData(anchorset = mda.anchors,
    verbose = FALSE)

DefaultAssay(mda.integrated) <- "integrated"
mda.integrated <- ScaleData(mda.integrated, verbose = FALSE)
mda.integrated <- RunPCA(mda.integrated, verbose = FALSE,npcs = 30)

mda.integrated <- FindNeighbors(object = mda.integrated, reduction = "pca", dims = 1:30)
mda.integrated <- FindClusters(mda.integrated, resolution =0.4)


save(mda.integrated,file='integrated_manno_mda_prog_only_woday51.Robj')
load('integrated_manno_mda_prog_only_woday51.Robj')


var.gene=VariableFeatures(object = mda.integrated)

combined_mat=as.matrix(mda.integrated@assays$integrated@scale.data)
#combined_mat=as.matrix(mda.integrated@assays$integrated@data)

mda.prog@meta.data$celltype=Idents(mda.prog)
Study_ID = rep(c('1', '2'), c(ncol(manno_prog@assays$RNA@data), ncol(mda.prog@assays$RNA@data)))
Celltype = c(as.character(manno_prog@meta.data$celltype),as.character(mda.prog@meta.data$celltype))

library(MetaNeighbor)
library(SummarizedExperiment)
dat=SummarizedExperiment(assays=list(counts=combined_mat))

celltype_NV=MetaNeighbor::MetaNeighborUS(var_genes = var.gene,
dat = dat,
study_id = Study_ID,
cell_type = Celltype,fast_version=T)

library(gplots)
library(RColorBrewer)
cols = rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdYlBu"))(100))
breaks = seq(0, 1, length=101)

ann_row=data.frame(dataset=c(rep("cell2016", 15), rep('This_study',11)))
rownames(ann_row)=rownames(celltype_NV)

pdf('pheatmap_ progenitor only manno our.pdf',width=10)
pheatmap(celltype_NV, 
         color = cols, cutree_rows=5,cutree_cols=5,
         breaks = breaks, cellwidth = 12, cellheight = 12,
         cluster_rows = TRUE,cluster_cols = TRUE, treeheight_col = 0,
         clustering_method = "ward.D2",
         annotation_row = ann_row,
         annotation_colors = list(dataset = c(cell2016 = "darkgreen", This_study = "deeppink"))
         )
dev.off()

pdf('heatmap progenitor only manno our scaledata woday51.pdf',height=10,width=10)
gplots::heatmap.2(celltype_NV,
margins=c(15,15),
keysize=1,
key.xlab="AUROC",
key.title="NULL",
trace = "none",
density.info =c("none"),
col = cols,
RowSideColors=c(rep("darkgreen", 15), rep('deeppink',11)),
breaks = breaks,
dendrogram="row",     # only draw a row dendrogram
offsetRow=0.1,
offsetCol=0.1,
cexRow = 0.8,
cexCol = 0.8)

par(lend = 1)      # square line ends for the color legend
legend("top",      # location of the legend on the heatmap plot
    legend = c("1-La Manno el al.2016", "2-This project"), # category labels
    col = c("darkgreen", "deeppink"),  # color key
    lty= 1,             # line style
    lwd = 10            # line width
)
dev.off()

