
library(Seurat)
library(tidyverse)
source('/pipe_function/seurat_pipe_cluster.R')

### step 1 ###
# for seperate time points datasets---------------------------------------------

path = "/all_filtered_feature_bc_matrix"

matrix.dir <- dir(path)
for(i in 1:length(matrix.dir)){
    data.dir=paste(path,'/',matrix.dir[i],sep='')
	mda=seuratv3_pipe(data.dir,project_name=matrix.dir[i])
    seuratv3_to_loom(mda)
}

### step 2 ###
##############################
### processed in scanpy ###
##############################

### step 3 ###



### step 4 ###

# for merged time points datasets---------------------------------------------
seuratv3_create_rename_pipe <- function(data.dir="",project_name) {
    project_name <- Read10X(data.dir = data.dir)
    print(paste('Creating ',project_name,' seurat object......',sep=''))
    mda <- CreateSeuratObject(counts = mda, project = project_name, min.cells = 3, min.features = 200)
    mda@meta.data$stage <- project_name
    mda <- RenameCells(mda, add.cell.id = project_name)
    return(mda)
}

path = "/all_filtered_feature_bc_matrix"
matrix.dir <- dir(path)
samples_list <- vector("list", length = length(matrix.dir))

for(i in 1:length(matrix.dir)){
    data.dir=paste(path,'/',matrix.dir[i],sep='')
    samples_list[[i]]=seuratv3_create_rename_pipe(data.dir,matrix.dir[i])
}

mda.merge <- merge(x = samples_list[[6]], y = samples_list[[1]], project = "mda_time_course")
mda.merge <- merge(x = mda.merge, y = samples_list[[2]], project = "mda_time_course")
mda.merge <- merge(x = mda.merge, y = samples_list[[3]], project = "mda_time_course")
mda.merge <- merge(x = mda.merge, y = samples_list[[4]], project = "mda_time_course")
mda.merge <- merge(x = mda.merge, y = samples_list[[5]], project = "mda_time_course")
mda.merge[["percent.mt"]] <- PercentageFeatureSet(object = mda.merge, pattern = "^MT-")
mda.merge[["percent.rp"]] <- PercentageFeatureSet(object = mda.merge, pattern = "^RPS|^RPL")
save(mda.merge,file='mda.merge.day8.14.21.28.35.51.rawstart.Robj')

setwd('/seurat_pipe_results')
load('mda.merge.day8.14.21.28.35.51.rawstart.Robj')

#------------------------------------------------------------------------------
day8_scanpy_meta=read.csv('/scanpy/day8_scanpy_meta.csv')
day8_celllist_tofilter=as_tibble(day8_scanpy_meta) %>% filter(leiden %in% c('3','4'))#c3 lowq, c4 apop

day14_scanpy_meta=read.csv('/scanpy/day14_scanpy_meta.csv')
day14_celllist_tofilter=as_tibble(day14_scanpy_meta) %>% 
						filter(leiden %in% c('4','6','8'))#c4 c6 lowq, c8 stress
#mutate_at('leiden',factor)
day21_scanpy_meta=read.csv('/scanpy/day21_scanpy_meta.csv')
day21_celllist_tofilter=as_tibble(day21_scanpy_meta) %>% 
						filter(leiden %in% c('4','8','9'))#c4 hypoxia, c8 lowq stress, c9 unkown, translation
#mutate_at('leiden',factor)

day28_scanpy_meta=read.csv('/scanpy/day28_scanpy_meta.csv')
day28_celllist_tofilter=as_tibble(day28_scanpy_meta) %>% 
						filter(leiden %in% c('2','5','12','13'))#c2 c5 c12 hypoxia, c13 stress

day35_scanpy_meta=read.csv('/scanpy/day35_scanpy_meta.csv')
day35_celllist_tofilter=as_tibble(day35_scanpy_meta) %>% 
						filter(leiden %in% c('7','11','13'))#c7 lowq, c11 highribo no marker, c13 lowq,

# day51_scanpy_meta=read.csv('/home/xupb/scRNA_data/time_course/manuscript/scanpy/day51_scanpy_meta.csv')
# day51_celllist_tofilter=as_tibble(day51_scanpy_meta) %>% 
# 						filter(leiden %in% c('6','7','8','9'))#
cells_to_filter=c(as.character(pull(day8_celllist_tofilter,1)),
				as.character(pull(day14_celllist_tofilter,1)),
				as.character(pull(day21_celllist_tofilter,1)),
				as.character(pull(day28_celllist_tofilter,1)),
				as.character(pull(day35_celllist_tofilter,1)))
cells_to_keep=setdiff(colnames(mda.merge),cells_to_filter)
length(cells_to_keep)

mda.merge <- subset(x = mda.merge, subset=stage == 'day51', invert=T)
mda.merge <- subset(x = mda.merge, subset=nFeature_RNA > 1000 & nFeature_RNA < 6000 & nCount_RNA>1000 &
                  nCount_RNA<40000 & percent.mt < 5)
mda.merge
mda.merge=subset(mda.merge,cells=cells_to_keep)

setwd('/seurat_pipe_results')
seuratv3_merge_pipe <- function(mda) {
	mda <- NormalizeData(object = mda, normalization.method = "LogNormalize", scale.factor = 40000)
	mda <- FindVariableFeatures(object = mda, selection.method = "vst", nfeatures = 2000)
	
	cc.genes <- readLines(con = "/seurat_pipe_results/regev_lab_cell_cycle_genes.txt")
	s.genes <- cc.genes[1:43]
	g2m.genes <- cc.genes[44:97]
	mda <- CellCycleScoring(object = mda, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
	mda@meta.data$CC.Difference <- mda@meta.data$S.Score - mda@meta.data$G2M.Score
	print(head(x = mda@meta.data, 5))
	
	#all.genes=rownames(mda) #features=all.genes,
	mda <- ScaleData(object = mda, 
	                vars.to.regress = c("nCount_RNA",'nFeature_RNA','CC.Difference','percent.mt'),
	                verbose = TRUE)
	var.genes=VariableFeatures(object = mda)
	print('Filtering TOP2A highly correlated genes from HVGs......')
	common_var_genes=hvg_filter_TOP2A_highly_corr_genes(mda,var.genes=var.genes)
	VariableFeatures(object = mda)=common_var_genes
	mda <- RunPCA(object =  mda, features=common_var_genes,npcs = 100, verbose = FALSE)
	var.genes.filename='merged_var_genes.csv'
	write.csv(VariableFeatures(object = mda),file=var.genes.filename)
	mda <- RunUMAP(object = mda, dims = 1:30)
	mda <- RunTSNE(object = mda, dims = 1:30)
	mda <- FindNeighbors(object = mda, dims = 1:30)
	mda <- FindClusters(object = mda, resolution = 0.6)
	print('Saving file in Robj format......')
	filename = "merged_seurat_pipe.Robj"
	mda.merge=mda
	save(mda.merge,file=filename)
	return(mda.merge)
}
mda.merge=seuratv3_merge_pipe(mda.merge)
seuratv3_to_loom(mda.merge)

library(harmony)
load('/merged_seurat_pipe.Robj')
mda=mda.merge
mda <- RunHarmony(mda, group.by.vars = "orig.ident",dims.use=1:70)
mda <- RunUMAP(mda, reduction = "harmony", dims = 1:25)
mda <- FindNeighbors(mda, reduction = "harmony", dims = 1:25) 
mda <- FindClusters(object = mda, resolution = 0.2)
table(Idents(object = mda),mda@meta.data$orig.ident)


harmony_embed=Embeddings(object = mda[["harmony"]])
write.csv(harmony_embed,file='harmony_embed.csv')

load('/seurat_pipe_results/merged_seurat_pipe.Robj')
meta=read.csv('/scanpy/merged_obj_meta.csv',row.names=1)
Idents(mda.merge)=meta$cluster
all.mda.leiden.markers<- FindAllMarkers(object = mda.merge, only.pos = TRUE, min.pct = 0.25,
                                           logfc.threshold = 0.25)

write.csv(all.mda.leiden.markers,file = 'all.mda.leiden.markers.csv')

top10 <-  all.mda.leiden.markers%>% group_by(cluster) %>% top_n(10, avg_logFC)
pdf('heatmap_all.mda.leiden.markers.pdf',height = 20,width =35)
DoHeatmap(mda.merge, features = top10$gene) + NoLegend()
dev.off()
### step 5 ###
##############################
### processed in scanpy ###
##############################


#### more ####
# calculate markers of all clusters
setwd('/seurat_pipe_results')
load('/seurat_pipe_results/merged_seurat_pipe.Robj')
mda.merge

meta=read.csv('/scanpy/woday51_scanpy_meta.csv',header=T,row.names=1)
Idents(mda.merge)=meta$leiden_number_new
table(Idents(mda.merge))
all.mda.scanpy.markers.leiden.number<- FindAllMarkers(object = mda.merge, only.pos = TRUE, min.pct = 0.25,
                                           logfc.threshold = 0.25)
write.csv(all.mda.scanpy.markers.leiden.number,file = 'all.mda.scanpy.markers.leiden.number.csv')

top30 <-  all.mda.scanpy.markers.leiden.number%>% group_by(cluster) %>% top_n(30, avg_logFC)
print(top30,n=800)


meta=read.csv('/scanpy/woday51_scanpy_meta.csv',header=T,row.names=1)
Idents(mda.merge)=meta$cluster
table(Idents(mda.merge))
all.mda.scanpy.markers<- FindAllMarkers(object = mda.merge, only.pos = TRUE, min.pct = 0.25,
                                           logfc.threshold = 0.25)

write.csv(all.mda.scanpy.markers,file = 'all.mda.scanpy.markers.csv')
top50 <-  all.mda.scanpy.markers%>% group_by(cluster) %>% top_n(50, avg_logFC)
print(top50,n=1000)



#### more ####
# calculate markers
load('/seurat_pipe_results/mda_day21_scaled.Robj')
mda_day21
day21_scanpy_meta=read.csv('/scanpy/day21_scanpy_meta.csv')
cells_to_filter=as_tibble(day21_scanpy_meta) %>% 
						filter(leiden %in% c('4','8','9'))#c4 hypoxia, c8 lowq stress, c9 unkown, translation

cells_to_keep=setdiff(colnames(mda_day21),cells_to_filter)
length(cells_to_keep)

meta=read.csv('/scanpy/day21_meta.csv')
Idents(mda_day21)=meta$Day21_union
table(Idents(mda_day21))

day21_all_marker=FindAllMarkers(mda_day21,only.pos = TRUE, min.pct = 0.25,
                                           logfc.threshold = 0.25)
write.csv(day21_all_marker,file = 'day21_all_marker.csv')

day21_dap_marker=FindMarkers(mda_day21,ident.1=c('P_MesenFP_LMX1A_Early'),only.pos = TRUE, min.pct = 0.25,
                                           logfc.threshold = 0.25)
write.csv(day21_dap_marker,file = 'day21_dap_marker.csv')


load('/seurat_pipe_results/mda_day28_scaled.Robj')
mda_day28

meta=read.csv('/scanpy/day28_meta.csv')
Idents(mda_day28)=meta$Day28_union
table(Idents(mda_day28))

day28_dap_marker=FindMarkers(mda_day28,ident.1=c('P_MesenFP_LMX1A_Late'),only.pos = TRUE, min.pct = 0.25,
                                           logfc.threshold = 0.25)
write.csv(day28_dap_marker,file = 'day28_dap_marker.csv')
