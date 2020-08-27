########################################################################
# Seurat_compare_PBMC.r
#
# This contains some typical Seurat clustering to compare to our own 
# clustering workflow for 10 PBMC data. See Seurat Tutorial at:
# https://satijalab.org/seurat/v3.1/pbmc3k_tutorial.html
#
#
# Copyright (C) 2020  Mahmoud M Ibrahim
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see http://www.gnu.org/licenses.
#
# Contact: mmibrahim@pm.me (Mahmoud M Ibrahim)
########################################################################

rm(list = ls())

##load packages
#plotting
suppressPackageStartupMessages(library(gplots))
suppressPackageStartupMessages(library(pheatmap))
#scRNA specific
suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(scran))
#tSNE
suppressPackageStartupMessages(library(umap))
#for graph clustering
suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(RANN))
suppressPackageStartupMessages(library(genesorteR))
library(mclust)

library(dplyr)
library(Seurat)

dist2d <- function(a,b,c) {
 v1 <- b - c
 v2 <- a - b
 m <- cbind(v1,v2)
 d <- abs(det(m))/sqrt(sum(v1*v1))
 return(d)
} 
makeTransparent = function(..., alpha=0.5) {
alpha = floor(255*alpha)  
newColor = col2rgb(col=unlist(list(...)), alpha=FALSE)
.makeTransparent = function(col, alpha) {
rgb(red=col[1], green=col[2], blue=col[3], alpha=alpha, maxColorValue=255)
}
newColor = apply(newColor, 2, .makeTransparent, alpha=alpha)
return(newColor)
}

col_palette = c("#d64e3c","#7dd554","#8641c6","#cfca45","#7778cb","#59803d","#d04d9c","#73d6a8","#492f60","#ccc497","#7f343b","#72acc0","#b97d40","#c796b5","#45483a","purple","green","yellow","red","blue","black")
col_palette_trans = makeTransparent(col_palette, alpha = 0.4)


##Seurat analysis (following exactly Seurat tutroial on this data, see link above)
pbmc.data = Read10X(data.dir = "/home/mibrahim/Dropbox/neverDelete/corrNet_tutorial/filtered_gene_bc_matrices/hg19")
pbmc = CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc[["percent.mt"]] = PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc = subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc = NormalizeData(pbmc)
pbmc = FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes = rownames(pbmc)
pbmc = ScaleData(pbmc, features = all.genes)
pbmc = RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc = FindNeighbors(pbmc, dims = 1:10)
pbmc = FindClusters(pbmc, resolution = 0.5)
pbmc = RunUMAP(pbmc, dims = 1:10)



#own analysis
l = pbmc@assays$RNA@counts #we start with the same cells used in Seurat analysis
l = SingleCellExperiment(assays=list(counts=l))
l = calculateQCMetrics(l)
l = scran::computeSumFactors(l, sizes = seq(10, 200, 10), positive = TRUE)
l = scater::normalize(l)
l = scater::calculateQCMetrics(l)


#initialize clustering
markerFile = "immune_markers_Villani2017.txt" #prior markers, see assets/public folder
ct = read.table(markerFile)[[1]]
markers = unique(as.vector(ct))
mcut = which(rownames(l) %in% markers)
l_v = l[mcut,] #get reduced dims
q = t(logcounts(l_v))
q = t(apply(q,1, function(x) x-mean(x)))
redDim = irlba::irlba(q, 10) #we use the same n dims as Seurat
temp = t(apply(redDim$u, 1, function(x) x / (sqrt(sum(x^2))))) #scale
allInd = 1:ncol(l) #build cell-cell similarity graph
kk = floor(sqrt(ncol(l)) * 1) * ncol(l)
set.seed(111)
kn = nn2(data = temp, k = ceiling((kk*1) / ncol(l_v)), searchtype = "priority", treetype = "bd")
kn = data.frame(node2 = as.vector(kn$nn.idx[,-1]), node1 = rep(kn$nn.idx[,1], ncol(kn$nn.idx) - 1), sim = 1 / (1+as.vector(kn$nn.dists[,-1])))
cut = min(head(sort(kn$sim, decreasing = TRUE), n = kk))
kn$sim2 = kn$sim
kn$sim2[kn$sim2 < cut] = 0
kn$sim2[kn$sim2 != 0] = 1
edge_list = cbind(kn$node1[which(kn$sim2 == 1)], kn$node2[which(kn$sim2 == 1)])
g = igraph::graph.edgelist(edge_list, directed = FALSE)
g = as.undirected(g, mode = "collapse")
g = simplify(g)
set.seed(111)
km = igraph::cluster_infomap(g) #cluster using infomap
class_info = km$membership
l$orig_class_info = class_info



#normalize (here we start iterative clustering)
l = scran::computeSumFactors(l, sizes = seq(10, 200, 10), positive = TRUE, clusters = class_info)
l = scater::normalize(l)
l = scater::calculateQCMetrics(l)
sg = sortGenes(counts(l), class_info, binarizeMethod = "naive")
var_genes = unique(unlist(plotTopMarkerHeat(sg, top_n = 500, outs = T, plotheat=F))) #1914 genes


#those lists will hold intermediate results from iterative clustering
varGenes = list()
meanClassAUC = list()
clusterLabels = list()
umapCoords = list()
graph = list()
rand = list()
source("/home/mibrahim/Dropbox/neverDelete/corrNet_tutorial/iterativeClustering.r")
pbmc@meta.data$cluster_labels = clusterLabels[[length(clusterLabels)]]


#merge clusters (same procedure we have in paper)
sg_j = sortGenes(counts(l), clusterLabels[[length(clusterLabels)]], binarizeMethod = "naive")
class_info = clusterLabels[[length(clusterLabels)]]
marks = plotTopMarkerHeat(sg_j, top_n=100, plotheat=FALSE, outs = TRUE)
xx = matrix(0, nrow = length(marks), ncol = length(marks))
for (i in 1:length(marks)) {
	for (j in 1:length(marks)) {
		xx[i,j] = length(intersect(marks[[i]], marks[[j]]))
	}
}
rownames(xx) = names(marks)
colnames(xx) = names(marks)
mergeMap = apply(xx, 2, function(x) names(which(x >= 80)))
mergethis = which(sapply(mergeMap,length) > 1)
if (length(mergethis) > 0) {
	mergeMap = mergeMap[mergethis]
	class_info_merging = class_info
	for (i in 1:length(mergeMap)) {
		class_info_merging[which(class_info_merging %in% mergeMap[[i]])] = paste0("merge", i)
	}
} else {
	class_info_merging = class_info
}
pbmc@meta.data$cluster_labels_merged = class_info_merging
l = scran::computeSumFactors(l, sizes = seq(10, 200, 10), positive = TRUE, clusters = class_info_merging)
l = scater::normalize(l)
l = scater::calculateQCMetrics(l)


#make a umap like in the paper (however we use 10 neighbours like in Seurat)
sg_j = sortGenes(logcounts(l), class_info_merging, cores = 16, binarizeMethod = "naive")
var_genes = unique(unlist(plotTopMarkerHeat(sg, top_n = 500, outs = T, plotheat=F)))
l_v = l[which(rownames(l) %in% var_genes),]
q = t(logcounts(l_v))
q = t(apply(q,1, function(x) x-mean(x)))
sv = irlba::irlba(q, 10) #in paper clustering this would be the result from fastMNN for batch correction. We use 10 dimensions because that's what Seurat does here
temp = t(apply(sv$u, 1, function(x) x / (sqrt(sum(x^2))))) #scale
write.csv(temp, file = "/home/mibrahim/Dropbox/neverDelete/corrNet_tutorial/svd_dims.txt", row.names = FALSE)
mapCoords = read.table("/home/mibrahim/Dropbox/neverDelete/corrNet_tutorial/umap_coords.txt", header = FALSE, sep=",") #run in python 


#compare embeddings on our clusters
colors_info = rep("", length = ncol(l))
for (i in 1:length(unique(class_info_merging))) {
	wtemp = which((class_info_merging) == (sort(unique(class_info_merging)))[i])
	colors_info[wtemp] = col_palette_trans[i]
}
colors_info_s = rep("", length = ncol(l))
for (i in 1:length(unique(Idents(pbmc)))) {
	wtemp = which((Idents(pbmc)) == (sort(unique(Idents(pbmc))))[i])
	colors_info_s[wtemp] = col_palette_trans[i]
}
pdf("/home/mibrahim/Dropbox/neverDelete/corrNet_tutorial/umap_compare.pdf", width = 10, height = 10)
par(mfrow = c(2,2))
plot(mapCoords[,1], mapCoords[,2], col = colors_info, pch = 20, xlab = "UMAP-1", ylab = "UMAP-2", main = "Our Embedding - Our Clusters")
legend(min(mapCoords[,1]), max(mapCoords[,2]) + 1, legend = abbreviate(sort(unique(class_info_merging))), pch=21, col=col_palette, pt.bg=col_palette_trans, pt.cex=2, cex=.8, bty="n", ncol = ceiling(length(unique(class_info_merging)) / 4))
plot(pbmc@reductions$umap@cell.embeddings, col = colors_info, pch = 20, xlab = "UMAP-1", ylab = "UMAP-2", main = "Seurat Embedding - Our Clusters")
plot(mapCoords[,1], mapCoords[,2], col = colors_info_s, pch = 20, xlab = "UMAP-1", ylab = "UMAP-2", main = "Our Embedding - Seurat Clusters")
legend(min(mapCoords[,1]), max(mapCoords[,2]) + 1, legend = abbreviate(sort(unique(Idents(pbmc)))), pch=21, col=col_palette, pt.bg=col_palette_trans, pt.cex=2, cex=.8, bty="n", ncol = ceiling(length(unique(Idents(pbmc))) / 4))
plot(pbmc@reductions$umap@cell.embeddings, col = colors_info_s, pch = 20, xlab = "UMAP-1", ylab = "UMAP-2", main = "Seurat Embedding - Seurat Clusters")
dev.off()



#comparison of Seurat clusters vs own clustering
tt = table(class_info_merging, Idents(pbmc))
tt = apply(tt, 2, function(x) x / sum(x))
tt = tt * 100
cluster_name = c("Naive T", "Memory T", "CD14+ Mono", "B", "CD8+ T", "FCGR3A+ Mono", "NK", "DC", "Platelet") #Seurat cluster names from the Seurat tutorial
pdf("/home/mibrahim/Dropbox/neverDelete/corrNet_tutorial/compare_sj.pdf")
heatmap.2(tt, trace = "none", dendrogram="none", col = colorRampPalette(c("white","dodgerblue4"))(n=100), margins = c(10,10), cellnote=round(tt,1), notecol="black", key = FALSE, labCol = cluster_name)
dev.off()
