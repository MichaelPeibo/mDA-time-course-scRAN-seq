

library(URD)
library(Seurat)
setwd('/seurat_pipe_results')
source('/pipe_function/seurat_pipe_cluster.R')

load('/seurat_pipe_results/merged_seurat_pipe.Robj')
mda.merge

scanpy_meta=read.csv('/subset/subset_da_lin_clusterinfo_woday51.csv',
	header=T,row.names=1)
var_meta=read.csv('/subset/subset_da_vars_woday51.csv',
	header=T,row.names=1)
mda.sub=subset(mda.merge,cells=rownames(scanpy_meta),features=rownames(var_meta))


mda.sub <- AddMetaData(
  object = mda.sub,
  metadata = as.character(scanpy_meta$cluster),
  col.name = 'cluster'
)

seuratv3_sub_pipe <- function(mda) {
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
	mda <- RunPCA(object =  mda, features=common_var_genes,npcs = 30, verbose = FALSE)
	var.genes.filename='sub_var_genes.csv'
	write.csv(VariableFeatures(object = mda),file=var.genes.filename)
	print('Saving file in Robj format......')
	mda <- RunUMAP(object = mda, dims = 1:30)
	mda <- RunTSNE(object = mda, dims = 1:30)
	mda <- FindNeighbors(object = mda, dims = 1:30)
	filename = "sub_seurat_pipe.Robj"
	save(mda,file=filename)
	return(mda)
}

mda.sub=seuratv3_sub_pipe(mda.sub)

source('/home/xupb/scRNA_data/time_course/day8_14_21_28_35/trimmed_final/urd/urd_seurat3_convert_function.R')
urd=seuratToURD2(mda.sub)
tsne_embed=Embeddings(object = mda.sub, reduction = "tsne")
colnames(tsne_embed)=c('tSNE1','tSNE2')
urd@tsne.y=as.data.frame(tsne_embed)

urd <- calcDM(urd, knn = 200, sigma='local')
saveRDS(urd@dm, file="urd/knn200_sigmalocal_dasub_rootcluster_woday51.rds")
root.cells<- cellsInCluster(urd, "cluster", '2-P_MesenFP_D14')
flood.result <- floodPseudotime(urd, root.cells=root.cells, n=100, minimum.cells.flooded=2, verbose=T)
saveRDS(flood.result, file="urd/flood-knn200_sigmalocal_dasub_rootcluster_woday51.rds")

#dm<- readRDS("urd/dm/knn200_sigmalocal_dasub_rootcluster.rds")
# Add it to the URD object
#urd <- importDM(urd, dm)

urd <- floodPseudotimeProcess(urd, flood.result, floods.name="pseudotime", 
                                  max.frac.NA=0.4, pseudotime.fun=mean, stability.div=20)
saveRDS(urd,file='urd/urd.floodpseudotime.rds')

		# pdf('density urd by cluster .pdf')
		# plotDists(urd, "pseudotime", "cluster", plot.title="Pseudotime by cluster")
		# dev.off()
		# scanpy_meta=read.csv('/home/xupb/scRNA_data/time_course/manuscript/scanpy/write/subset_da_lin_clusterinfo_new.csv',
		# 	header=T,row.names=1)

		# # load('sub_seurat_pipe.Robj')
		# # head(mda@meta.data$cluster)
		# # df1=scanpy_meta['cluster']
		# # df2=as.data.frame(urd@pseudotime)
		# # df=cbind(df1,df2)

		# # p=ggplot(mapping = aes(x = df$pseudotime, fill = factor(df$cluster)))+
		# #   geom_density(alpha = 0.5)
		# #   #scale_fill_manual(values = cluster_colours, labels = cluster_labels, name = "") +
		# #   #theme(legend.position = c(0.0, 0.8))

		# library(ggpubr)
		# p=ggdensity(df, x = "pseudotime", fill = "cluster",
		#    add = "mean", rug = F, xscale='log2')

		# pdf('test density.pdf',width=15)
		# p
		# dev.off()

#urd<- readRDS("urd/urd.floodpseudotime.rds")
head(urd@pseudotime)
write.csv(as.data.frame(urd@pseudotime),file='urd/urd_pseudotime.csv')


# Consider all genes expressed in 1% of cells
frac.exp <- rowSums(urd@logupx.data > 0) / ncol(urd@logupx.data)
expressed.genes <- names(frac.exp)[which(frac.exp > 0.01)]
length(expressed.genes)
interestingGenes <- c("TH", "PTPRO", "LMX1A",'CLSTN2','CD83','CORIN','PITX3','NR4A2')
# any missing?
interestingGenes[which(!interestingGenes %in% expressed.genes)]


# Calculate spline fit
expressed.spline.5cell <- geneSmoothFit(urd, method="spline", pseudotime="pseudotime", 
cells = colnames(urd@count.data), genes = expressed.genes, moving.window = 1, 
cells.per.window = 5, spar=0.875) #.875 previously
saveRDS(expressed.spline.5cell,file='urd/expressed.spline.5cell.count.0.875%.rds')


# Which genes change in their actual mean expression value by at least 0.5?
spline.change.real <- apply(expressed.spline.5cell$mean.smooth, 1, function(x) diff(range(x)))
genes.change.real <- names(which(spline.change.real >= 0.6)) # 0.5, 0.6 is good 
interestingGenes <- c("TH", "PTPRO", "LMX1A",'CLSTN2','CD83','CORIN','NEUROD1','NEUROG2','PITX3','NR4A2')
# any missing?
interestingGenes[which(!interestingGenes %in% genes.change.real)]


# Which genes are well fit by the spline curves? (Noise is usually poorly fit by a curve)
spline.fit.5cell <- apply(expressed.spline.5cell$scaled.smooth - expressed.spline.5cell$scaled.expression.red, 1, function(i) sum(i^2))
spline.fit.norm.5cell <- spline.fit.5cell / ncol(expressed.spline.5cell$scaled.expression.red)
genes.wellfit.5cell <- names(which(spline.fit.norm.5cell <= 0.045)) ## 0.045

interestingGenes[which(!interestingGenes %in% genes.wellfit.5cell)]

# Which genes change in their scaled log2 mean expression value sufficiently?
# At least 30%, and requiring more change (up to 40%) as the data is less well fit
# by its spline curve.
# TO GET genes such as neurod1 and neurog2, I set threshold 0.15
change.scale.5cell <- apply(expressed.spline.5cell$scaled.smooth, 1, function(x) diff(range(x)))
genes.scale.5cell <- names(which(change.scale.5cell >= spline.fit.norm.5cell * .10/0.045 + 0.3))
interestingGenes[which(!interestingGenes %in% genes.scale.5cell)]


# Ensure that genes are fit by the spline curve significantly better than a flat line of slope 0.
# (Weighted by distance to next point in pseudotime to compensate for point density)
w <- ncol(expressed.spline.5cell$scaled.smooth) - 1
weight <- diff(as.numeric(colnames(expressed.spline.5cell$scaled.smooth))) * 1000
spline.fit.weighted <- apply(expressed.spline.5cell$scaled.smooth[,1:w] - expressed.spline.5cell$scaled.expression.red[,1:w], 1, function(i) sum(weight * i^2))
spline.flat.fit.weighted <- apply(expressed.spline.5cell$scaled.expression.red[,1:w], 1, function(x) sum(weight * (x - mean(x))^2))
spline.fit.ratio <- log2(spline.flat.fit.weighted / spline.fit.weighted)
spline.fit.betterthanflat <- names(which(spline.fit.ratio > 0.25))#0.25
interestingGenes[which(!interestingGenes %in% spline.fit.betterthanflat)]

# Take the intersection of those genes and use them in the heatmap & analysis
varying.genes.5cell <- intersect(intersect(intersect(genes.change.real, genes.scale.5cell), genes.wellfit.5cell), spline.fit.betterthanflat)
length(varying.genes.5cell)

interestingGenes[which(!interestingGenes %in% varying.genes.5cell)]

## in Rstudio
setwd("D:/Data/manuscripts/mda_manu/figure2/figs/woday51/")
expressed.spline.5cell=readRDS('expressed.spline.5cell.count.0.875%.rds')
library(ggplot2)
library(ggsci)
library(ComplexHeatmap)

plotSmoothFit_v1 <- function(smoothed.fit, genes){
	library(ggsci)
	expression <- smoothed.fit$scaled.expression[genes,,drop=F]
 	expression$Gene <- rownames(expression)
 	expression.melt <- reshape2::melt(expression, id.vars="Gene", variable.name="Pseudotime", value.name="Expression")
 	expression.melt$Pseudotime <- as.numeric(as.character(expression.melt$Pseudotime))
 
 g=ggplot(data=expression.melt,
       aes_string(x="Pseudotime", y="Expression",group='Gene', color='Gene')) + 
  theme_classic() + 
  #geom_point(data=expression.melt, alpha=0.2) + 
  geom_smooth(data=expression.melt,size=1, aes(color=Gene),span=0.5)+
  #scale_color_jco() 
  #scale_color_d3()
  #scale_color_rickandmorty()
  scale_color_igv() + 
  #scale_color_jama() + 
  theme(legend.text=element_text(size=14),axis.title=element_text(size=14)) # for last one with legend, comment out this line axis.title=element_blank()
return(g)
}

genes1=c('ASCL1','PTPRO','CD83')
genes2=c('PITX3','NR4A2','TH')
genes3=c('CNPY1','ETV4','PAX8')
genes4=c('CLSTN2','ADAMTS9','SLC25A39')
genes5=c('SLC39A8','WNT5A','RSPO2')


p1=plotSmoothFit_v1(expressed.spline.5cell, genes = genes1)
p2=plotSmoothFit_v1(expressed.spline.5cell, genes = genes2)
p3=plotSmoothFit_v1(expressed.spline.5cell, genes = genes3)
p4=plotSmoothFit_v1(expressed.spline.5cell, genes = genes4)
p5=plotSmoothFit_v1(expressed.spline.5cell, genes = genes5)

library(cowplot)
all <- plot_grid(p1, p2, p3, p4,p5,ncol=1, align = "v")
save_plot("cow_all_v1.pdf", all,base_width=4,base_height=9)

# Heatmap Basics
cols <- scales::gradient_n_pal(RColorBrewer::brewer.pal(9, "YlOrRd"))(seq(0,1,length.out = 50))

# Reduce the data in the splines, such that each block is minimum 0.01
# pseudotime; in this case, each block will be >= 5 cells, >= pseudotime 0.01
s <- expressed.spline.5cell
colnames(s$mean.expression) <- as.character(round(as.numeric(colnames(s$mean.expression)), digits=2))
s$mean.expression <- matrixReduce(s$mean.expression)
colnames(s$mean.smooth) <- as.character(round(as.numeric(colnames(s$mean.smooth)), digits=2))
s$mean.smooth <- matrixReduce(s$mean.smooth)

# Re-scale spline curves min/max for the heatmap
ss <- sweep(s$mean.smooth, 1, apply(s$mean.smooth, 1, min), "-")
ss <- sweep(ss, 1, apply(ss, 1, max), "/")

# Hierarchical cluster based on smoothed expression
h.ss <- hclust(dist(as.matrix(ss[varying.genes.5cell,])), method = "complete")
h.ss.d <- as.dendrogram(h.ss)
k=5 # Cluster number chosen by eye.
h.ss.clust <- cutree(h.ss, k=k)

# Get cluster order as it will be in the heatmap
clust.order <- unique(h.ss.clust[h.ss$order])
h.ss.clust.ord <- plyr::mapvalues(from=clust.order, to=1:k, x = h.ss.clust)

# Generate cluster color vector
cluster.h <- seq(0,1,length.out=k+1)
cluster.s <- rep(c(1, 0.75), length=k)
cluster.v <- rep(c(1, 0.75), length=k)
cluster.colors <- hsv(h=cluster.h[1:k], s=cluster.s, v=cluster.v)
h.ss.clust.col <- plyr::mapvalues(from=as.character(1:k), to=cluster.colors, x=h.ss.clust.ord)

pdf('urd/varing genes heatmap k5 meth complete.pdf',width=17, height=30)
gplots::heatmap.2(
  x = as.matrix(ss[varying.genes.5cell[h.ss$order],]), Rowv=F, RowSideColors=h.ss.clust.col[h.ss$order], 
  Colv=F, dendrogram="none", col=cols, trace="none", density.info="none", 
  key=F, cexCol=0.8, cexRow=0.08, margins = c(8,8), lwid=c(0.3,4), lhei=c(0.3, 4), labCol=NA)
dev.off()

	# mat=as.matrix(ss[varying.genes.5cell[h.ss$order],])
	# right_mark_gene <- c("DDC","SYT1")
	# right_gene_pos <- which(rownames(mat) %in% right_mark_gene)

	# right_row_anno <-  rowAnnotation(mark_gene = anno_mark(at = right_gene_pos,
	#                                      labels = right_mark_gene))

	# pdf('ComplexHeatmap cascade.pdf')
	# Heatmap(mat,
	#         cluster_rows = F,
	#         cluster_columns = F,
	#         show_column_names = FALSE,
	#         show_row_names = FALSE,
	#         #column_split = 2,
	#         #column_title =c('Others','LMX1A+EN1+'),
	#         #top_annotation=top_annotation,
	#         right_annotation = right_row_anno,
	#         heatmap_legend_param = list(
	#          title = "Z-Score",
	#          title_position = "leftcenter-rot"
	#         ))
	# dev.off()


				# h.ss.clust.col=readRDS('h.ss.clust.col.rds')
				# clust_order=readRDS('clsut_order.rds')
				# h.ss.clust.col[clust_order]
				# pdf('urd/phenograph test.pdf',width=17, height=22)
				# gplots::heatmap.2(
				#   x = as.matrix(ss[varying.genes.5cell[clust_order],]), Rowv=F, RowSideColors=as.character(h.ss.clust.col[clust_order]), 
				#   Colv=F, dendrogram="none", col=cols, trace="none", density.info="none", 
				#   key=F, cexCol=0.8, cexRow=0.08, margins = c(8,8), lwid=c(0.3,4), lhei=c(0.3, 4), labCol=NA)
				# dev.off()

# Generate the actual heatmap and save as a PDF.
so <- expressed.spline.5cell
# Aggregate by mean across clusters
so$scaled.expression <- stats::aggregate(s$scaled.expression[varying.genes.5cell,
], by = list(h.ss.clust.ord[varying.genes.5cell]), FUN = mean)
# Rename rows to add a leading 0 to 1-9 so ggplot2 will sort correctly
rownames(so$scaled.expression) <- sprintf("%02d", as.numeric(rownames(so$scaled.expression)))
so$mean.expression <- stats::aggregate(s$mean.expression[varying.genes.5cell, ],
by = list(h.ss.clust.ord[varying.genes.5cell]), FUN = mean)
rownames(so$mean.expression) <- sprintf("%02d", as.numeric(rownames(so$mean.expression)))
so$scaled.smooth <- stats::aggregate(s$scaled.smooth[varying.genes.5cell, ], by = list(h.ss.clust.ord[varying.genes.5cell]),
FUN = mean)
rownames(so$scaled.smooth) <- sprintf("%02d", as.numeric(rownames(so$scaled.smooth)))
so$mean.smooth <- stats::aggregate(s$mean.smooth[varying.genes.5cell, ], by = list(h.ss.clust.ord[varying.genes.5cell]),
FUN = mean)
rownames(so$mean.smooth) <- sprintf("%02d", as.numeric(rownames(so$mean.smooth)))
so$mean.expression.red <- stats::aggregate(s$mean.expression.red[varying.genes.5cell,
], by = list(h.ss.clust.ord[varying.genes.5cell]), FUN = mean)
rownames(so$mean.expression.red) <- sprintf("%02d", as.numeric(rownames(so$mean.expression.red)))
so$scaled.expression.red <- stats::aggregate(s$scaled.expression.red[varying.genes.5cell,
], by = list(h.ss.clust.ord[varying.genes.5cell]), FUN = mean)
rownames(so$scaled.expression.red) <- sprintf("%02d", as.numeric(rownames(so$scaled.expression.red)))
# Plot expression profiles of each cluster
p=plotSmoothFit(so, sort(rownames(so$mean.expression)), scaled = T, plot.data = T,
multiplot = T) + ggtitle("") + ylim(0,1)
p=p+ theme_bw()+ theme(panel.border = element_blank(), panel.grid.major = element_blank(),
panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

pdf("aggregates meth complete k5.pdf",width = 15, height = 15)
p
dev.off()

tf=read.csv('/home/xupb/scRNA_data/time_course/manuscript/seurat_pipe_results/urd/human_tf.csv',row.names=1)
tf.genes <- tf$gene_symbol
# Get ectoderm varying genes that are also TFs
varying.tfs <- intersect(tf.genes, varying.genes.5cell)

# Filter out certain genes to remove them from heatmaps for presentational purposes
filter.heatmap.genes <- function(genes) {
  mt.genes <- grep("^MT-", ignore.case=T, genes, value=T)
  ribo.genes <- grep("^RPL|^RPS", ignore.case=T, genes, value=T)
  return(setdiff(genes, c(mt.genes, ribo.genes)))
}
varying.tfs.rest=filter.heatmap.genes(varying.tfs)
# Figure out how they should be ordered on the plot
tf.order <- t(apply(ss[varying.tfs,] > 0.75, 1, function(x) {
  y <- which(x)
  return(c(min(y), max(y)))
}))
tfs.ordered <- varying.tfs[order(tf.order[,1], tf.order[,2], decreasing=c(F, F), method="radix")]

# Make a heatmap of just the TFs
pdf("varying TFs.pdf", width=8.5, height=10)
gplots::heatmap.2(as.matrix(ss[tfs.ordered,]),Rowv=F, Colv=F, dendrogram="none", col=cols, 
	trace="none", density.info="none", key=F, cexCol=0.8, cexRow=0.2, margins = c(5,13), 
	lwid=c(0.3,4), lhei=c(0.35,4), labCol=NA)
dev.off()


