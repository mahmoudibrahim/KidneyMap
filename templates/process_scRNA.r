########################################################################
# process_scRNA.r
#
# Copyright (C) 2019-2020  Mahmoud M Ibrahim
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
options(scipen=1000000)

#load packages
#plotting
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(vioplot))
suppressPackageStartupMessages(library(gplots))
#scRNA specific
suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(scran))
suppressPackageStartupMessages(library(genesorteR))
#IRLB (fast approximate SVD)
suppressPackageStartupMessages(library(irlba))
#for graph clustering
suppressPackageStartupMessages(library(RANN))
suppressPackageStartupMessages(library(igraph))
#general
suppressPackageStartupMessages(library(mclust))
suppressPackageStartupMessages(library(umap))

#custom functions
source("https://raw.githubusercontent.com/mahmoudibrahim/KidneyMap/master/.customlib.r")




species = "human" #"human" or "mouse"


if (species == "human") {
	rRNA = "https://raw.githubusercontent.com/mahmoudibrahim/KidneyMap/master/assets/public/gencode/rRNA_gencode29_geneNames.txt"
	mtRNA = "https://raw.githubusercontent.com/mahmoudibrahim/KidneyMap/master/assets/public/gencode/mtRNA_gencode29_geneNames.txt"
} else if (species == "mouse") {
	rRNA = "https://raw.githubusercontent.com/mahmoudibrahim/KidneyMap/master/assets/public/gencode/rrnanames.gencode.vM20.txt"
	mtRNA = "https://raw.githubusercontent.com/mahmoudibrahim/KidneyMap/master/assets/public/gencode/chrMnames.gencode.vM20.txt"
}


direc = "/directory/to/alevin/output"
filename = "/prefix/for/saving/plots_and_files"



###########################################################################
########################START ANALYSIS#####################################
###########################################################################

#load data
l = ReadAlevin(direc)
l = ceiling(as(l, "dgCMatrix"))



#removing rRNA and mtRNA genes from main matrix
rRNA_genes = read.table(rRNA)[[1]]
tempGenes = unlist(lapply(strsplit(as.character(rownames(l)), split = ";", fixed = T), function(x) as.character(x[1])))
rRNA_genes = which(tempGenes %in% rRNA_genes) 
mtRNA_genes = read.table(mtRNA)[[1]]
tempGenes = unlist(lapply(strsplit(as.character(rownames(l)), split = ";", fixed = T), function(x) as.character(x[1])))
mtRNA_genes = which(tempGenes %in% mtRNA_genes)
egfp = which(tempGenes == "eGFP") #if there is eGFP trans-gene

if (species == "human") {
	par_y = which(sapply(rownames(l), function(x) grepl("_PAR_Y",x)) == TRUE) #remove PAR_Y genes, identical between X and Y chromosomes
	removethis = unique(sort(c(rRNA_genes, par_y, mtRNA_genes)))
} else {
	removethis = unique(sort(c(rRNA_genes, mtRNA_genes, egfp))) #note egfp!
}

#store mtRNA and rRNA in a separate matrix
if (length(rRNA_genes) > 0) {
	rRNA_genes = l[rRNA_genes,]
}
if (length(mtRNA_genes) > 0) {
	mtRNA_genes = l[mtRNA_genes,] 
}
egfp = l[which(tempGenes == "eGFP"),] #for egfp
if (length(removethis) > 0) {	#now clean main matrix
	l = l[-removethis,] 
}





#create SingleCellExperiment objects
l = SingleCellExperiment(assays=list(counts=l))
lx = l #save the full matrix before any filtering
l = calculateQCMetrics(l)

if(any(l$total_counts == 0)) { #some cells might only have rRNA or mtRNA
	ww = which(l$total_counts == 0)
	l = l[,-ww]
	l = calculateQCMetrics(l)
	mtRNA_genes = mtRNA_genes[,-ww]
	rRNA_genes = rRNA_genes[,-ww]
}

mtRNA_genes = SingleCellExperiment(assays=list(counts=mtRNA_genes))
mtRNA_genes = calculateQCMetrics(mtRNA_genes)
rRNA_genes = SingleCellExperiment(assays=list(counts=rRNA_genes))
rRNA_genes = calculateQCMetrics(rRNA_genes)




#some cell quality plots
zzz = pdf(paste0(filename, "_DepthPerCell.pdf"), height = 9)
par(mfrow = c(2,2))
hist(l$log10_total_counts, breaks = "scott", xlab = "log10[Library Depth]", main = "Depth per Cell", ylab = "# Cells", col = "#D55E00")
plot(log10(1:(ncol(l))), sort(l$total_counts, decreasing = TRUE), col = "#99999980", pch = 20, ylab = "Library Depth", xlab = "Decreasing Cell Read Depth (# of Cells)", log = "y", xaxt = "n", main = paste0(ncol(l), " cells"))
minor.ticks.axis(1,9,mn=0,mx=6)
hist(l$total_features_by_counts, breaks = "scott", xlab = "# Expressed Genes", main = "Captured Genes per Cell", ylab = "# Cells", col = "#D55E00")
plot(log10(l$total_counts+mtRNA_genes$total_counts), (mtRNA_genes$total_counts / (l$total_counts+mtRNA_genes$total_counts))*100, col = "#99999950", pch = 20, ylab = "% mtRNA Reads", xlab = "log10[Library Depth]")
zzz = dev.off()
zzz = pdf(paste0(filename, "_DepthPerCell_rRNA.pdf"), height = 9)
par(mfrow = c(2,2))
hist(l$log10_total_counts, breaks = "scott", xlab = "log10[Library Depth]", main = "Depth per Cell", ylab = "# Cells", col = "#D55E00")
plot(log10(rRNA_genes$total_counts), (mtRNA_genes$total_counts / (l$total_counts+mtRNA_genes$total_counts))*100, col = "#99999950", pch = 20, ylab = "% mtRNA Reads", xlab = "log10[rRNA Library Depth]")
plot(log10(l$total_counts+rRNA_genes$total_counts), (rRNA_genes$total_counts / (l$total_counts+rRNA_genes$total_counts))*100, col = "#99999950", pch = 20, ylab = "% rRNA Reads", xlab = "log10[Library Depth]")
zzz = dev.off()




######cell and gene filtering####### Note: each one of these filtering steps is optional. Always check diagnostic plots and decide. It's good to include borderline case cells than to remove them

#(1) filter cells based on expression, only if the distribution of log10 total counts looks bimodal
k_libSize = Mclust(log10(ceiling(l$total_counts)), G = 2, modelNames="E")
k_libSize_sel = which.max(k_libSize$parameters$mean)
sele = which(k_libSize$classification == k_libSize_sel)
l = l[,sele]
mtRNA_genes = mtRNA_genes[,sele]
rRNA_genes = rRNA_genes[,sele]

geneZero.drop = which(Matrix::rowSums(counts(l)) == 0) #if some genes now are not expressed in any cells
if (length(geneZero.drop) > 0) {
	l = l[-geneZero.drop,]
}
l = calculateQCMetrics(l)
mtRNA_genes = calculateQCMetrics(mtRNA_genes)
rRNA_genes = calculateQCMetrics(rRNA_genes)
a = plotExprsFreqVsMean(l) #calculate gene stats



#(2) mitochondrial filtering, the idea is to select an automated reasonable cutoff given the data at hand, based on overall RNA content and Mt content (which are often opposite, high Mt corresponds to low RNA content)
mtPer = (mtRNA_genes$total_counts / (l$total_counts+mtRNA_genes$total_counts))*100
mm = Mclust(cbind(l$log10_total_counts, mtPer), 2, modelNames="EII")
choose_mt = which.min(mm$parameters$mean[2,])
choose_mt = which(mm$classification == choose_mt)
mtRNA_genes = mtRNA_genes[, choose_mt]
rRNA_genes = rRNA_genes[, choose_mt]
zzz = pdf(paste0(filename, "_MtFiltering.pdf"), height = 9)
plot(mm, what = "classification", xlab = "log10 Total Counts", ylab = "Mt Percent")
dev.off()


#(3) filtering on RNA content and features detected (logic: counts should increase with features, outliers (often very few) be gone!)
mm = Mclust(cbind(l$log10_total_counts, l$log10_total_features_by_counts), 2, modelNames=c("VEV","VEE"))
choose_count = which.max(mm$parameters$pro)
choose_count = which(mm$classification == choose_count)
l = l[,choose_count]
mtRNA_genes = mtRNA_genes[, choose_count]
rRNA_genes = rRNA_genes[, choose_count]
zzz = pdf(paste0(filename, "_countANDfeatureFiltering.pdf"), height = 9)
plot(mm, what = "classification", xlab = "log10 Total Counts", ylab = "log10 Number of Genes")
dev.off()

#(4) filtering on feature capture and Ribosomal protein content. Some cells have relatively few expressed genes. Percent of counts in the top 500 features should not be too high. Here also the cutoff is adaptive to the dataset at hand. Ribosomal proteins and pseudogenes tend to be hallmarks of low quality clusters. We check here for cells that have very high content of these genes. Also here the cutoff is adaptive to the dataset at hand.
mad_percent = median(l$pct_counts_in_top_500_features) + (median(abs(l$pct_counts_in_top_500_features - median(l$pct_counts_in_top_500_features))) * 5)
if (mad_percent > 100) {mad_percent = 100}
choose_percent = which(l$pct_counts_in_top_500_features < mad_percent)
zzz = pdf(paste0(filename, "_topFeatureFiltering.pdf"), height = 9)
plot(l$pct_counts_in_top_500_features, l$log10_total_counts, pch = 19, col = "#99999930", main = round(mad_percent, 2), xlab = "% of Counts in Top 500 Genes", ylab = "log10 Total Counts")
dev.off()

rbRNA = grep("MRP|RPL|RPS", rownames(l), value=T) 
rrr = Matrix::colSums(counts(l[which(rownames(l) %in% rbRNA),]))
norrr = Matrix::colSums(counts(l[-(which(rownames(l) %in% rbRNA)),]))
rrr = rrr*100 / (rrr + norrr)
mad_percent2 = median(rrr) + (median(abs(rrr - median(rrr))) * 5)
if (mad_percent2 > 50) {mad_percent2 = 50}
choose_percent2 = which(rrr < mad_percent2)

zzz = pdf(paste0(filename, "_ribosomalProteinFiltering.pdf"), height = 9)
plot(rrr, l$log10_total_counts, pch = 19, col = "#99999930", main = round(mad_percent2, 2), xlab = "% of Ribosomal Protein Genes", ylab = "log10 Total Counts")
dev.off()
zzz = pdf(paste0(filename, "_ribosomalProtein+topFeatureFiltering.pdf"), height = 9)
plot(rrr, l$pct_counts_in_top_500_features, col = "#99999950", pch = 19, xlab = "Percent of Ribosomal Protein Genes", ylab = "% of Counts in Top 500 Genes")
lines(0:100, rep(mad_percent, 101), lwd = 2, lty = 2, col = "red")
lines(rep(mad_percent2, 101), 0:100, lwd = 2, lty = 2, col = "red")
dev.off()

choose_percent = intersect(choose_percent, choose_percent2) #now remove cells
l = l[,choose_percent]
mtRNA_genes = mtRNA_genes[, choose_percent]
rRNA_genes = rRNA_genes[, choose_percent]
geneZero.drop = which(Matrix::rowSums(counts(l)) == 0) #check if some genes should be dropped
if (length(geneZero.drop) > 0) {
	l = l[-geneZero.drop,]
}
l = calculateQCMetrics(l)
mtRNA_genes = calculateQCMetrics(mtRNA_genes)
rRNA_genes = calculateQCMetrics(rRNA_genes)
a = plotExprsFreqVsMean(l)


#gene filtering (only for this script, full matrix is exported for further analysis)
cutGenes = log10(floor(ncol(l) * 0.005)) #on average at least 1UMI per 0.5% of cells
choose_genes = which((log10(Matrix::rowSums(counts(l)))) >= cutGenes)
zzz = pdf(paste0(filename, "_geneFiltering.pdf"), height = 9)
plot(log10(Matrix::rowSums(counts(l))), a$data$Y, pch = 20, col = "#99999950", xlab = "log10[Total Counts]", ylab = "% of Cells Detected", cex = 0.8)
lines(rep(cutGenes,101), 0:100, col = "red", lty = 2, lwd = 2)
dev.off()
l = l[choose_genes,]
l = calculateQCMetrics(l)
a = plotExprsFreqVsMean(l)
aveRpG = ( (l$total_counts) / (l$total_features_by_counts) )




#some cell quality plots (after filtering)
zzz = pdf(paste0(filename, "_DepthPerCell_afterFiltering.pdf"), height = 9)
par(mfrow = c(2,2))
hist(l$log10_total_counts, breaks = "scott", xlab = "log10[Library Depth]", main = "Depth per Cell", ylab = "# Cells", col = "#D55E00")
plot(log10(1:(ncol(l))), sort(l$total_counts, decreasing = TRUE), col = "#99999980", pch = 20, ylab = "Library Depth", xlab = "Decreasing Cell Read Depth (# of Cells)", log = "y", xaxt = "n", main = paste0(ncol(l), " cells"))
minor.ticks.axis(1,9,mn=0,mx=6)
hist(l$total_features_by_counts, breaks = "scott", xlab = "# Expressed Genes", main = "Captured Genes per Cell", ylab = "# Cells", col = "#D55E00")
plot(log10(l$total_counts+mtRNA_genes$total_counts), (mtRNA_genes$total_counts / (l$total_counts+mtRNA_genes$total_counts))*100, col = "#99999950", pch = 20, ylab = "% mtRNA Reads", xlab = "log10[Library Depth]")
zzz = dev.off()
zzz = pdf(paste0(filename, "_DepthPerCell_rRNA_afterFiltering.pdf"), height = 9)
par(mfrow = c(2,2))
hist(l$log10_total_counts, breaks = "scott", xlab = "log10[Library Depth]", main = "Depth per Cell", ylab = "# Cells", col = "#D55E00")
plot(log10(rRNA_genes$total_counts), (mtRNA_genes$total_counts / (l$total_counts+mtRNA_genes$total_counts))*100, col = "#99999950", pch = 20, ylab = "% mtRNA Reads", xlab = "log10[rRNA Library Depth]")
plot(log10(l$total_counts+rRNA_genes$total_counts), (rRNA_genes$total_counts / (l$total_counts+rRNA_genes$total_counts))*100, col = "#99999950", pch = 20, ylab = "% rRNA Reads", xlab = "log10[Library Depth]")
dev.off()
zzz = pdf(paste0(filename, "_TotalGeneCounts_afterFiltering.pdf"))
hist(log10(Matrix::rowSums(counts(l))), breaks = "scott", xlab = "log10[Total Counts]", main = "Total Counts per Gene", ylab = "# Genes", col = "#D55E00", xlim = c(0,5), xaxp = c(0,5,10))
zzz = dev.off()
zzz = pdf(paste0(filename, "_pct_top_features_afterFiltering.pdf"))
plot(l$pct_counts_in_top_50_features, l$pct_counts_in_top_500_features, pch = 16, ylim = c(40,100), yaxp = c(40,100,6), col = "#99999940")
zzz = dev.off()
zzz = pdf(paste0(filename, "_RpG.pdf"))
hist(aveRpG, breaks = "scott", main = "Average Transcripts per Gene per Cell", xlab = "Ave. Transcripts per Gene per Cell", ylab = "# of Genes", col = "#D55E00", xlim = c(0,10))
zzz = dev.off()

######cell and gene filtering (done!)####### 




#normalize expression
l = scran::computeSumFactors(l, sizes = seq(10, 200, 10), positive = TRUE)
ws = which(sizeFactors(l) <= 0)
if(length(ws) > 0) {
	l= l[,-ws]
	mtRNA_genes = mtRNA_genes[,-ws]
	rRNA_genes = rRNA_genes[,-ws]
	l = calculateQCMetrics(l)
	mtRNA_genes = calculateQCMetrics(mtRNA_genes)
	rRNA_genes = calculateQCMetrics(rRNA_genes)
	message("Warning: some cells had negative size factors and were excluded!")
}
zzz = pdf(paste0(filename, "_sizeFactors.pdf"))
plot(sizeFactors(l), l$total_counts/1e3, xlab = "Size factors", ylab = "Library Depth (thousands)", col = "#99999950", pch = 20)
zzz = dev.off()
l = scater::normalize(l)
l = scater::calculateQCMetrics(l)
a = plotExprsFreqVsMean(l)





#reduce dims
markerFile = "https://raw.githubusercontent.com/mahmoudibrahim/KidneyMap/master/assets/public/all_markers_MMI_Aug2019.txt" #a priori highly variable genes
ct = read.table(markerFile, sep = "\t", header = T)
markers = unique(as.vector(ct$Gene))
gene_symbols = toupper(unlist(lapply(strsplit(as.character(rownames(l)), split = ";", fixed = T), function(x) as.character(x[2]))))
mcut = which(gene_symbols %in% markers)
l_v = l[mcut,]

#SVD
q = t(logcounts(l_v))
q = t(apply(q,1, function(x) x-mean(x)))
ss = svd(q)
d = ss$d
d2 =  d^2 / (sum(d^2))

#take the knee
if (length(d) < 1000) {
	xx = 1:length(d)
	yy = d
} else {
	xx = 1:1000
	yy = d[1:1000]
}
p1 = c(xx[1],yy[1])
p2 = c(xx[length(xx)], yy[length(yy)])
dw = which.max(sapply(1:length(yy), function(x) dist2d(c(xx[x], yy[x]), p1, p2)))
zzz = pdf(paste0(filename, "_pcComp.pdf"))
plot(1:100, d[1:100], xlim = c(1,100), xaxp = c(1,100,99), pch = 19, col = "#99999949")
lines(rep(dw, 1000), seq(0,10000,length.out = 1000), col = "red", lty = 2)
zzz = dev.off()
dw = 1:dw



#UMAP just to tell if some specific group of cells is low quality (high rRNA, low eGFP...etc.)
kk = floor(sqrt(nrow(ss$u)) * 1) * nrow(ss$u)
temp = t(apply(ss$u[,dw], 1, function(x) x / (sqrt(sum(x^2)))))
load("https://github.com/mahmoudibrahim/KidneyMap/blob/master/templates/umap.defaults.RData?raw=true")
set.seed(111)
umap.defaults$n_neighbors = kk/(nrow(ss$u))
rr = umap(temp, config = umap.defaults)
rr = rr$layout

zzz = pdf(paste0(filename, "_tSNE_totalCounts.pdf"))
plot(rr, col = color.gradient(scale(l$total_counts), colors = c("gray",rev(c("#d7191c","#fdae61","#abdda4")))), pch = 20, xlab = "tSNE-1", ylab = "tSNE-2", frame.plot = FALSE, ylim = c(min(rr[,2])-10,max(rr[,2])+10))
zzz = dev.off()
zzz = pdf(paste0(filename, "_tSNE_MtPercent.pdf"))
plot(rr, col = color.gradient(scale((mtRNA_genes$total_counts / (l$total_counts+mtRNA_genes$total_counts))), colors = c("gray",rev(c("#d7191c","#fdae61","#abdda4")))), pch = 20, xlab = "tSNE-1", ylab = "tSNE-2", frame.plot = FALSE, ylim = c(min(rr[,2])-10,max(rr[,2])+10))
zzz = dev.off()
zzz = pdf(paste0(filename, "_tSNE_rRNAPercent.pdf"))
plot(rr, col = color.gradient(scale((rRNA_genes$total_counts / (l$total_counts+rRNA_genes$total_counts))), colors = c("gray",rev(c("#d7191c","#fdae61","#abdda4")))), pch = 20, xlab = "tSNE-1", ylab = "tSNE-2", frame.plot = FALSE, ylim = c(min(rr[,2])-10,max(rr[,2])+10))
zzz = dev.off()





#clustering
temp = t(apply(ss$u[,dw], 1, function(x) x / (sqrt(sum(x^2)))))
allInd = 1:nrow(ss$u)
kk = floor(sqrt(nrow(ss$u)) * 1) * nrow(ss$u)

set.seed(111)
kn = nn2(data = temp, k = ceiling((kk*1) / ncol(l_v)), searchtype = "priority", treetype = "bd")
kn = data.frame(node2 = as.vector(kn$nn.idx[,-1]), node1 = rep(kn$nn.idx[,1], ncol(kn$nn.idx) - 1), sim = 1 / (1+as.vector(kn$nn.dists[,-1])))
cut = min(head(sort(kn$sim, decreasing = TRUE), n = kk))
kn$sim2 = kn$sim
kn$sim2[kn$sim2 < cut] = 0
kn$sim2[kn$sim2 != 0] = 1
edge_list = cbind(kn$node1[which(kn$sim2 == 1)], kn$node2[which(kn$sim2 == 1)])

g = igraph::graph.edgelist(edge_list, directed = FALSE)
g = simplify(g)

g = as.undirected(g, mode = "collapse")
g = simplify(g)
set.seed(111)
km = igraph::cluster_infomap(g)
class_info = km$membership
colors_info = rep("", length = ncol(l_v))
for (i in 1:length(unique(class_info))) {
	wtemp = which((class_info) == (sort(unique(class_info)))[i])
	colors_info[wtemp] = col_palette_trans[i]
}
zzz = pdf(paste0(filename, "_UMAP.pdf"))
plot(rr, col = colors_info, pch = 20, xlab = "UMAP-1", ylab = "UMAP-2", frame.plot = FALSE, ylim = c(min(rr[,2])-10,max(rr[,2])+10))
legend(min(rr[,1]), max(rr[,2]) + 10, legend = 1:length(unique(class_info)), pch=21, col=col_palette, pt.bg=col_palette_trans, pt.cex=2, cex=.8, bty="n", ncol = ceiling(length(unique(class_info)) / 2))
zzz = dev.off()



###remove singular cell clusters and remove bad clusters?
if (any(table(class_info) == 1)) {
	www = which(table(class_info) == 1)
	www = which(class_info %in% www)
	l = l[,-www]
	class_info = class_info[-www]
	rr = rr[-www,]
	colors_info = colors_info[-www]
	message(paste0(length(www), " cell(s) excluded because they formed singular clusters."))
}


pdf(paste0(filename, "_beforeFiltering_totalCounts.pdf"))
boxplot(l$log10_total_counts ~ class_info)
dev.off()
sg = sortGenes(counts(l), class_info, binarizeMethod = "naive")
write.files(sg, prefix = paste0(filename,"_beforeFiltering_gs"))
marks = unique(unlist(plotTopMarkerHeat(sg, top_n=500, outs = TRUE, plotheat=FALSE)))

##if there are low quality clusters, based on a combination of no diff exp genes, marker genes are all ribo-proteins or house keeping genes, very low RNA content. Example:
www = which(class_info %in% c("8","15")) #remove clusters
l = l[,-www]
class_info = class_info[-www]
rr = rr[-www,]
colors_info = colors_info[-www]
##########################################






#Supervised assignment of clusters to epithelial, endothelial, immune, mesenchymal and neuronal (rather ad hoc, some manual checking is still done and some assignments are corrected)
gene_symbols = toupper(unlist(lapply(strsplit(as.character(rownames(l)), split = ";", fixed = T), function(x) as.character(x[2]))))
ll = l
rownames(ll) = gene_symbols
sg = sortGenes(logcounts(ll), class_info, binarizeMethod = "naive")
write.files(sg, prefix = paste0(filename,"_afterFiltering_gs"))
zzz = pdf(paste0(filename, "_clusterSimilarity.pdf"))
plotCorrelationHeat(sg, getMarkers(sg, quant=0.95)$markers)
dev.off()
toplist = plotTopMarkerHeat(sg, top_n = 50, plotheat=FALSE, outs=TRUE)
markerFile = "/home/mibrahim/Dropbox/dev/kidneyMap/public/assets/public/all_markers_MMI_Aug2019_Mar2020.txt"

ct = read.table(markerFile, sep = "\t", header = T)
ct = ct[which(ct$Source != "Clark2019"),]
cell_types = names(table(ct$Cell_Type2))
cell_types = cell_types[c(-1,-5)]
inter_mat = matrix(0, ncol = ncol(sg$specScore), nrow = length(cell_types))
for (i in 1:length(cell_types)) {
	for (j in 1:ncol(sg$specScore)) {
		pm1 = unique(toupper(ct$Gene[which(ct$Cell_Type2 == cell_types[i])]))
		pm2 = toplist[[j]]
		inter_mat[i,j] = length(intersect(pm1, pm2))
	}
}

rownames(inter_mat) = c(cell_types)
norma = table(as.character(ct$Cell_Type2[which(ct$Cell_Type2 %in% cell_types)])) / rowSums(inter_mat)
enri = apply(inter_mat,2,function(x)x/norma)
enri = apply(enri, 2, function(x) x/sum(x))
enri = t(scale(t(enri)))
rownames(enri) = rownames(inter_mat)
colnames(enri) = names(sg$classProb)
enri = enri * 100
enri[which(is.nan(enri))] = 0

zzz = pdf(paste0(filename, "_japanium_cellTypeAssignment.pdf"), height = 15)
heatmap.2(t(enri), trace = "none", scale = "none", col = colorRampPalette(c("blue","white","red"))(n=100), margins = c(16,4), notecol="black", notecex=1, density.info = "none", cellnote=round(t(enri),1))
dev.off()



whichCat = rownames(enri)[apply(enri,2,which.max)]
#here do manual correction if needed. For example, whichCat[1] = "SMC"
classes = as.integer(as.factor(class_info))
map = data.frame(oldO = 1:length(unique(classes)), newO = whichCat)
classes = map$newO[match(classes, map$oldO)]
colors_info2 = rep("", length = ncol(l))
for (i in 1:length(unique(classes))) {
	wtemp = which((classes) == (sort(unique(classes)))[i])
	colors_info2[wtemp] = col_palette_trans[i]
}
zzz = pdf(paste0(filename, "_UMAP_cellCat.pdf")) #plot as a UMAP to check
plot(rr, col = colors_info2, pch = 20, xlab = "UMAP-1", ylab = "UMAP-2", frame.plot = FALSE, ylim = c(min(rr[,2])-10,max(rr[,2])+10))
legend(min(rr[,1]), max(rr[,2]) + 10, legend = levels(classes), pch=21, col=col_palette, pt.bg=col_palette_trans, pt.cex=2, cex=.8, bty="n", ncol = ceiling(length(unique(classes)) / 2))
zzz = dev.off()

#main major type markers
hend_markers = toupper(c("Pdgfrb", "Ptprc", "Pecam1", "Epcam", "Nrxn3", "MME"))
names1 = toupper(unlist(lapply(strsplit(as.character(rownames(l)), split = ";", fixed = T), function(x) as.character(x[2]))))
for (i in 1:length(hend_markers)) {
	if(any(toupper(names1) == hend_markers[i])) {
		zzz = pdf(paste0(filename, "_markers_box_", hend_markers[i], ".pdf"))
		exp = logcounts(l)
		gn = which(toupper(names1) == hend_markers[i])
		plot(as.numeric(as.vector(class_info))+rnorm(ncol(exp), 0, 0.13), exp[gn,], col = makeTransparent(colors_info, alpha = 0.6), pch = 20, ylab = "Log Expression", xlab = "Cell Cluster", cex = 0.5, main = hend_markers[i])
		boxplot(exp[gn,] ~ class_info, col = makeTransparent("white", alpha = 0), border = makeTransparent("black", alpha = 0.5), range = 0.1, outline = F, add = T, notch = T, lwd = 3)
		dev.off()
	}
}





###doublet score
set.seed(111)
dbs = log10(doubletCells(l, subset.row = which(rownames(l) %in% marks), d = length(dw)) + 1)
zzz = pdf(paste0(filename, "_doubletScore_.pdf"))
boxplot(dbs ~ class_info, col = makeTransparent("white", alpha = 0), border = makeTransparent("black", alpha = 0.5), range = 0, outline = F, add = F, notch = T, lwd = 3)
dev.off()
#############################################



#messages
l = calculateQCMetrics(l)
a = plotExprsFreqVsMean(l) #re-calculate gene stats
aveRpG = ( (l$total_counts) / (l$total_features_by_counts) )

cat(paste0("There are ", ncol(l), " cells.\n"), file = paste0(filename, "_japanium.info.txt"))
cat(paste0("There are ", nrow(l), " expressed genes, with an average of ", round(mean(l$total_features_by_counts)), " expressed genes per cell.\n"), file = paste0(filename, "_japanium.info.txt"), append = T)
cat(paste0("There are on average ", round(mean(l$total_counts)), " UMIs per cell.\n"), file = paste0(filename, "_japanium.info.txt"), append = T)
cat(paste0("In each cell, there are on average ", round(mean(aveRpG),2), " transcripts per gene.\n"), file = paste0(filename, "_japanium.info.txt"), append = T)
cat(paste0("Timestamp: ", Sys.time()), file = paste0(filename, "_japanium.info.txt"), append = T)


#write data
lx = lx[,which(colnames(lx) %in% colnames(l))]
save(lx = lx, file = paste0(filename, "-japanium.filtered_matrix.RData"))
write(ws, file = paste0(filename, "_noNormaCells.txt"), ncolumns = 1)
write(colnames(l), file = paste0(filename, "_cellNames.txt"), ncolumns = 1)
write(class_info, file = paste0(filename, "_cellClasses.txt"), ncolumns = 1)
write(as.character(classes), file = paste0(filename, "_cellCats.txt"), ncolumns = 1)


##Fin
