library(Seurat)
library(dplyr)
library(loomR)

hvg_filter_TOP2A_highly_corr_genes = function(seurat_obj,var.genes){
    matrix<-seurat_obj@assays$RNA@data
    matrix_mod<-as.matrix(matrix)
    cycle_genes='TOP2A'
    gene<-as.numeric(matrix_mod[cycle_genes,])
    correlations<-apply(matrix_mod,1,function(x){cor(gene,x)})
    low_cellcyle_corr=which(correlations<0.15)
    common_var_genes=intersect(var.genes,names(low_cellcyle_corr))
    cat('There is',length(common_var_genes),'filtered hvg genes', '\n')
    rm(matrix,matrix_mod,cycle_genes,gene,correlations,low_cellcyle_corr)
    return(common_var_genes)
}

seuratv3_pipe <- function(data.dir="",project_name) {
    mda <- Read10X(data.dir = data.dir)
    print(paste('Creating ',project_name,' seurat object......',sep=''))
    mda <- CreateSeuratObject(counts = mda, project = project_name, min.cells = 3, min.features = 200)
    mda[["percent.mt"]] <- PercentageFeatureSet(object = mda, pattern = "^MT-")
    mda[["percent.rp"]] <- PercentageFeatureSet(object = mda, pattern = "^RPS|^RPL")
  
    mda <- subset(x = mda, subset=nFeature_RNA > 1000 & nFeature_RNA < 6000 & nCount_RNA>1000 &
                  nCount_RNA<40000 & percent.mt < 5)
    print(paste('Subsetting ',project_name,' seurat object......',sep=''))
    print(mda)
	mda <- NormalizeData(object = mda, normalization.method = "LogNormalize", scale.factor = 40000)
	mda <- FindVariableFeatures(object = mda, selection.method = "vst", nfeatures = 2000)
	
	cc.genes <- readLines(con = "/home/xupb/regev_lab_cell_cycle_genes.txt")
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
	mda <- RunPCA(object =  mda, features=common_var_genes,npcs = 60, verbose = FALSE)
	var.genes.filename=paste(project_name,'_var_genes.csv',sep='')
	write.csv(VariableFeatures(object = mda),file=var.genes.filename)
	print('Saving file in Robj format......')
	mda <- RunUMAP(object = mda, dims = 1:30)
	mda <- RunTSNE(object = mda, dims = 1:30)
	mda <- FindNeighbors(object = mda, dims = 1:30)
	mda <- FindClusters(object = mda, resolution = 0.6)
	filename = paste(project_name,"_seurat_pipe.Robj",sep="")
	project_name=mda
	save(project_name,file=filename)
	return(project_name)
}

#convert to loom file-----------------------------------------
seuratv3_to_loom=function(seurat_obj){
seurat_obj@graphs <- list()
filename = paste(seurat_obj@project.name,"_seurat2loom.loom",sep="")
seurat_obj=as.loom(x = seurat_obj,filename=filename)
seurat_obj$close_all()
}

#-------------------------------------------------------------

