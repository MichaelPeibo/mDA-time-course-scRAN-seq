
library(dplyr)
library(msigdf)
library(Seurat)
library(ggplot2)
library(ggsci)
library(readxl)
library(viridis)
library(ComplexHeatmap)


expressed.spline.5cell=readRDS('/home/xupb/scRNA_data/time_course/manuscript/seurat_pipe_results/woday51/urd/expressed.spline.5cell.count.0.875%.rds')
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


lr=read_excel('/home/xupb/scRNA_data/time_course/manuscript/pseudotime/Curated and putative ligand-receptor pairs in human.xlsx',
  sheet = 'All.Pairs')

reac_wnt = msigdf.human %>% filter(geneset=="REACTOME_SIGNALING_BY_WNT")
reac_shh = msigdf.human %>% filter(geneset=="REACTOME_SIGNALING_BY_HEDGEHOG")
reac_bmp = msigdf.human %>% filter(geneset=="REACTOME_SIGNALING_BY_BMP")
reac_tgfb = msigdf.human %>% filter(geneset=="REACTOME_SIGNALING_BY_TGF_BETA_FAMILY_MEMBERS")
reac_fgfr = msigdf.human %>% filter(geneset=="REACTOME_SIGNALING_BY_FGFR")
reac_fgfr1 = msigdf.human %>% filter(geneset=="REACTOME_SIGNALING_BY_FGFR1")
pid_fgf = msigdf.human %>% filter(geneset=="PID_FGF_PATHWAY")

# varying ligands
reac_wnt_varying_lig=intersect(intersect(reac_wnt$symbol,varying.genes.5cell),lr$Ligand.ApprovedSymbol)
reac_shh_varying_lig=intersect(intersect(reac_shh$symbol,varying.genes.5cell),lr$Ligand.ApprovedSymbol)
# reac_bmp_varying_lig=intersect(intersect(reac_bmp$symbol,varying.genes.5cell),lr$Ligand.ApprovedSymbol) # none
# reac_tgfb_varying_lig=intersect(intersect(reac_tgfb$symbol,varying.genes.5cell),lr$Ligand.ApprovedSymbol)# none
reac_fgfr_varying_lig=intersect(intersect(reac_fgfr$symbol,varying.genes.5cell),lr$Ligand.ApprovedSymbol)
reac_fgfr1_varying_lig=intersect(intersect(reac_fgfr1$symbol,varying.genes.5cell),lr$Ligand.ApprovedSymbol) #same as fgfr
pid_fgf_varying_lig=intersect(intersect(pid_fgf$symbol,varying.genes.5cell),lr$Ligand.ApprovedSymbol) #same as fgfr
fgf_lig_union=union(reac_fgfr1_varying_lig,pid_fgf_varying_lig)

# varying receptors
reac_wnt_varying_rec=intersect(intersect(reac_wnt$symbol,varying.genes.5cell),lr$Receptor.ApprovedSymbol)
reac_shh_varying_rec=intersect(intersect(reac_shh$symbol,varying.genes.5cell),lr$Receptor.ApprovedSymbol)
# reac_bmp_varying_rec=intersect(intersect(reac_bmp$symbol,varying.genes.5cell),lr$Receptor.ApprovedSymbol) # none
# reac_tgfb_varying_rec=intersect(intersect(reac_tgfb$symbol,varying.genes.5cell),lr$Receptor.ApprovedSymbol) # none
# reac_fgfr_varying_rec=intersect(intersect(reac_fgfr$symbol,varying.genes.5cell),lr$Receptor.ApprovedSymbol) # covered in pid
# reac_fgfr1_varying_rec=intersect(intersect(reac_fgfr1$symbol,varying.genes.5cell),lr$Receptor.ApprovedSymbol) #same as fgfr
pid_fgf_varying_rec=intersect(intersect(pid_fgf$symbol,varying.genes.5cell),lr$Receptor.ApprovedSymbol)

pathway.genes=unique(c(reac_wnt_varying_lig,reac_shh_varying_lig,fgf_lig_union,
          reac_wnt_varying_rec,reac_shh_varying_rec,pid_fgf_varying_rec))

s <- expressed.spline.5cell
colnames(s$mean.expression) <- as.character(round(as.numeric(colnames(s$mean.expression)), digits=2))
s$mean.expression <- URD::matrixReduce(s$mean.expression)
colnames(s$mean.smooth) <- as.character(round(as.numeric(colnames(s$mean.smooth)), digits=2))
s$mean.smooth <-URD:: matrixReduce(s$mean.smooth)

# Re-scale spline curves min/max for the heatmap
ss <- sweep(s$mean.smooth, 1, apply(s$mean.smooth, 1, min), "-")
ss <- sweep(ss, 1, apply(ss, 1, max), "/")
h.ss <- hclust(dist(as.matrix(ss[pathway.genes,])), method = "complete")
h.ss.d <- as.dendrogram(h.ss)
k=4 # Cluster number chosen by eye.
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
# Heatmap Basics
cols <- scales::gradient_n_pal(RColorBrewer::brewer.pal(9, "YlOrRd"))(seq(0,1,length.out = 50))

pdf('pathway lig rec heatmap.pdf',width=17, height=30)
gplots::heatmap.2(
  x = as.matrix(ss[pathway.genes,]), Rowv=F, RowSideColors=h.ss.clust.col[h.ss$order], 
  Colv=F, dendrogram="none", col=cols, trace="none", density.info="none", 
  key=F, cexCol=0.8, cexRow=2, margins = c(8,8), lwid=c(0.3,4), lhei=c(0.3, 4), labCol=NA)
dev.off()


df = data.frame(Pathway = c(rep("WNT", 3),rep("SHH", 2),rep("FGF", 1),
            rep("WNT", 1),rep("SHH", 1),rep("FGF", 3)),
        Ligand_Receptor= c(rep('Ligand',6),rep('Receptor',5))
        )
row_ann_col = list(Pathway = c("WNT" = "red",
             "SHH" = "blue",
             "FGF" = "pink"),
      Ligand_Receptor = c("Ligand" = "orange",
          "Receptor" = "darkgreen"))
ha = rowAnnotation(df = df,width = unit(1, "cm"),col=row_ann_col)
#levels(df$Pathway)=c("WNT_Ligand","SHH_Ligand","FGF_Ligand" ,"WNT_Receptor","SHH_Receptor" , "FGF_Receptor")

pdf('complexheatmap pathway lig rec test.pdf',width=8.5,height=10)
Heatmap(as.matrix(ss[pathway.genes,]),name = "Scaled Expression", 
  show_column_names=F,
  show_row_names=T,
  cluster_rows =hclust(dist(as.matrix(ss[pathway.genes,])),method='ward.D2'),
  #row_order  = h.ss$order,
    cluster_columns = F,
    right_annotation = ha ,
    col = cols  
  )
dev.off()


























Pseudotime=readRDS('/home/xupb/scRNA_data/time_course/manuscript/pseudotime/da_lin_pseudotime.rds')
colnames(Pseudotime)='Pseudotime'

load('/home/xupb/scRNA_data/time_course/manuscript/seurat_pipe_results/sub_seurat_pipe.Robj')
mda@meta.data$Pseudotime=Pseudotime

shh_ptch1=c('SHH','PTCH1')
wnt5a_lrp5=c('WNT5A','LRP5')
wnt4_fzd6=c('WNT4','FZD6')
WNT4=c('WNT4')

mda <- AddModuleScore(mda,  features = list(shh_ptch1), name='SHH_PTCH1')
mda <- AddModuleScore(mda,  features = list(wnt5a_lrp5), name='WNT5A_LRP5')
mda <- AddModuleScore(mda,  features = list(wnt4_fzd6), name='WNT4_FZD6')
mda <- AddModuleScore(mda,  features = list(WNT4), name='WNT4')

plotLR <- function(Pseudotime, seurat_object, lr_name){
  meta=seurat_object@meta.data
  data <- data.frame(Pseudotime = Pseudotime,
                     LR_score = meta[,lr_name])
  data_plot <- ggplot(data, aes(x=Pseudotime, y=LR_score)) +
    geom_point(size=0.3, shape=20,aes(color=Pseudotime)) +#, aes(color=Sample)
    geom_smooth(color="black") +
    scale_color_viridis(option = "D") + 
    #scale_color_manual(values=color_ramp) +
    scale_x_continuous(expand=c(0,0), breaks=c(0.1, 0.3, 0.5)) +
    scale_y_continuous(expand=c(0,0)) +
    ylab("LR Pair Score") +
    theme_classic() +
    theme(legend.position="none",
          axis.text = element_text(size=17, color="black"),
          axis.title = element_blank(),
          axis.line = element_line(size=0.5),
          axis.ticks = element_line(size=0.5))
  
  data_plot
}
p=plotLR(Pseudotime,mda,'WNT41')

pdf('WNT4 score.pdf')
p
dev.off()



# ann_colors = list( binary_proj=RColorBrewer::brewer.pal(n=8, name = 'Dark2'))
# pheno=data.frame(pheno[,3])
# colnames(pheno)=c('binary_proj')

library(wesanderson)
ann_colors = list(
    Cluster=c('#ffff00', '#1ce6ff', '#ff34ff', '#ff4a46', '#008941', '#006fa6', '#a30059'),
    Layer=c("red", "blue", "green"),
    #binary_proj = c(wes_palette("GrandBudapest1"),wes_palette("GrandBudapest2"))
    binary_proj=RColorBrewer::brewer.pal(n=8, name='Set1')
    #binary_proj=ggsci::pal_ucscgb("default")(8)
    )

names(ann_colors$Layer)=levels(pheno$Layer)
names(ann_colors$Cluster)=levels(pheno$Cluster)
names(ann_colors$binary_proj)=levels(pheno$binary_proj)

library(ComplexHeatmap)
col<- circlize::colorRamp2(c(-1,0,1,2,3),viridis::inferno(5))

pdf('complexheatmap barcodes binary projection inferno.pdf',width=15,height=10)
Heatmap(deg.mat,name = "Scaled Expression", 
	show_column_names=F,
	cluster_rows =hclust(dist(deg.mat)),
	#row_order=as.character(unique(deg_genes)[drop=T]), 
    cluster_columns = hclust(dist(t(deg.mat))),
	top_annotation = HeatmapAnnotation(df=pheno,col=ann_colors),
	col = col
	#viridis::inferno(20)
	#col = circlize::colorRamp2(c(-3, 0, 3), c("green", "white", "red"))
	)
dev.off()


