########################################################################
# clusterCells.r
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



cellType = "SMC" #SMC(=mesenchymal) or Hematopoietic or Neuronal or Epithelial or Endothelial
direc = "/path/to/folder/for/output"
outName = paste0(cellType, "_combined")
outName = paste0(direc, "/", outName)
setwd(direc)


#character vector including paths to all RData output objects of individual libraries to contain the cell category of interest (SMC, Hematopoietic..etc.), example:
filenames = c("CK84_CD10minus_alevin-japanium.filtered_matrix.RData", "CK85_CD10minus_alevin-japanium.filtered_matrix.RData")

#those lists will hold intermediate results from iterative clustering
varGenes = list()
meanClassAUC = list()
clusterLabels = list()
umapCoords = list()
graph = list()
rand = list()
mat = list()
file_id = list()
cellCluster = list()

#load the RData objects
for (i in 1:length(filenames)) {
	load(filenames[i])
	mat[[i]] = lx
	noNorma = paste0(strsplit(as.character(filenames[i]), split = "-japanium", fixed = T)[[1]][1], "_noNormaCells.txt")
	if (file.info(noNorma) > 1) {
		noNorma = read.table(noNorma)[[1]]
		mat[[i]] = mat[[i]][,-noNorma]
	}
}


#filter for desired cell type
for (i in 1:length(filenames)) {
	cellCat = as.character(read.table(paste0(strsplit(as.character(filenames[i]), split = "-japanium", fixed = T)[[1]][1], "_cellCats.txt"))[[1]])
	ww = which(cellCat == cellType)
	
	mat[[i]] = mat[[i]][,ww]


	file_id[[i]] = strsplit(as.character(filenames[i]), split = "/", fixed = T)
	file_id[[i]] = unlist(file_id[[i]])[length(unlist(file_id[[i]]))]
	file_id[[i]] = rep(unlist(file_id[[i]]), ncol(mat[[i]]))

	cellCluster[[i]] = as.character(read.table(paste0(strsplit(as.character(filenames[i]), split = "-japanium", fixed = T)[[1]][1], "_cellClasses.txt"))[[1]])
	cellCluster[[i]] = cellCluster[[i]][ww]
	cellCluster[[i]] = paste0("File", i, "-", cellCluster[[i]])
}


#make one object
l = mat[[1]]
for (i in 2:length(filenames)) {
	l = cbind(l, mat[[i]])
}
l$colData = data.frame(fileID = unlist(file_id), origCellCluster = unlist(cellCluster)) #this doesn't really define colData() object as SingleCellExperiment understand
l = calculateQCMetrics(l)
lx = l #save the full matrix for later




#some cell quality plots
zzz = pdf(paste0(outName, "_DepthPerCell.pdf"), height = 9)
par(mfrow = c(2,2))
hist(l$log10_total_counts, breaks = "scott", xlab = "log10[Library Depth]", main = "Depth per Cell", ylab = "# Cells", col = "#D55E00")
plot(log10(1:(ncol(l))), sort(l$total_counts, decreasing = TRUE), col = "#99999980", pch = 20, ylab = "Library Depth", xlab = "Decreasing Cell Read Depth (# of Cells)", log = "y", xaxt = "n", main = paste0(ncol(l), " cells"))
minor.ticks.axis(1,9,mn=0,mx=6)
hist(l$total_features_by_counts, breaks = "scott", xlab = "# Expressed Genes", main = "Captured Genes per Cell", ylab = "# Cells", col = "#D55E00")
zzz = dev.off()
zzz = pdf(paste0(outName, "_TotalGeneCounts.pdf"))
hist(log10(Matrix::rowSums(counts(l))), breaks = "scott", xlab = "log10[Total Counts]", main = "Total Counts per Gene", ylab = "# Genes", col = "#D55E00", xlim = c(0,5), xaxp = c(0,5,10))
zzz = dev.off()
zzz = pdf(paste0(outName, "_pct_top_features.pdf"))
plot(l$pct_counts_in_top_50_features, l$pct_counts_in_top_500_features, pch = 16, ylim = c(40,100), yaxp = c(40,100,6), col = "#99999940", xlab = "% of Counts in Top 50 Genes", ylab = "% of Counts in Top 500 Genes")
zzz = dev.off()

a = plotExprsFreqVsMean(l)





#initial gene filtering
geneZero.drop = which(Matrix::rowSums(counts(l)) == 0)
if (length(geneZero.drop) > 0) {
	l = l[-geneZero.drop,]
}
l = calculateQCMetrics(l)
a = plotExprsFreqVsMean(l) #re-calculate gene stats
cutGenes = log10(floor(ncol(l) * 0.005)) #on average at least 1UMI in 0.5% of cells 
choose_genes = which((log10(Matrix::rowSums(counts(l)))) >= cutGenes)
zzz = pdf(paste0(outName, "_geneFiltering.pdf"), height = 9)
plot(log10(Matrix::rowSums(counts(l))), a$data$Y, pch = 20, col = "#99999950", xlab = "log10[Total Counts]", ylab = "% of Cells Detected", cex = 0.8)
lines(rep(cutGenes,101), 0:100, col = "red", lty = 2, lwd = 2)
dev.off()
l = l[choose_genes,]
l = calculateQCMetrics(l)
a = plotExprsFreqVsMean(l)

aveRpG = ( (l$total_counts) / (l$total_features_by_counts) )
zzz = pdf(paste0(outName, "_RpG.pdf"))
hist(aveRpG, breaks = "scott", main = "Average Transcripts per Gene per Cell", xlab = "Ave. Transcripts per Gene per Cell", ylab = "# of Genes", col = "#D55E00", xlim = c(0,5))
zzz = dev.off()




#some cell quality plots (after fitlering)
zzz = pdf(paste0(outName, "_DepthPerCell_afterFiltering.pdf"), height = 9)
par(mfrow = c(2,2))
hist(l$log10_total_counts, breaks = "scott", xlab = "log10[Library Depth]", main = "Depth per Cell", ylab = "# Cells", col = "#D55E00")
plot(log10(1:(ncol(l))), sort(l$total_counts, decreasing = TRUE), col = "#99999980", pch = 20, ylab = "Library Depth", xlab = "Decreasing Cell Read Depth (# of Cells)", log = "y", xaxt = "n", main = paste0(ncol(l), " cells"))
minor.ticks.axis(1,9,mn=0,mx=6)
hist(l$total_features_by_counts, breaks = "scott", xlab = "# Expressed Genes", main = "Captured Genes per Cell", ylab = "# Cells", col = "#D55E00")
zzz = dev.off()
zzz = pdf(paste0(outName, "_TotalGeneCounts_afterFiltering.pdf"))
hist(log10(Matrix::rowSums(counts(l))), breaks = "scott", xlab = "log10[Total Counts]", main = "Total Counts per Gene", ylab = "# Genes", col = "#D55E00", xlim = c(0,5), xaxp = c(0,5,10))
zzz = dev.off()
zzz = pdf(paste0(outName, "_pct_top_features_afterFiltering.pdf"))
plot(l$pct_counts_in_top_50_features, l$pct_counts_in_top_500_features, pch = 16, ylim = c(40,100), yaxp = c(40,100,6), col = "#99999940", xlab = "% of Counts in Top 50 Genes", ylab = "% of Counts in Top 500 Genes")
zzz = dev.off()
a = plotExprsFreqVsMean(l) #calculate gene stats





#normalize expression
l = scran::computeSumFactors(l, sizes = seq(10, 200, 20), positive = TRUE)
ws = which(sizeFactors(l) <= 0)
if(length(ws) > 0) {
	message("Warning: some cells had negative size factors!")
	sizeFactors(l)[ws] = min(sizeFactors(l)[sizeFactors(l) > 0])
}
zzz = pdf(paste0(outName, "_japanium_scran_sizeFactors.pdf"))
plot(sizeFactors(l), l$total_counts/1e3, xlab = "Size factors", ylab = "Library Depth (thousands)", col = "#99999950", pch = 20)
zzz = dev.off()
l = scater::normalize(l)
l = scater::calculateQCMetrics(l)
a = plotExprsFreqVsMean(l)





#(supervised) gene selection, markers of clusters from the process_scRNA.r clustering
sg = sortGenes(logcounts(l), l$colData$origCellCluster, binarizeMethod = "naive", cores = 16)
sg_files = sortGenes(logcounts(l), l$colData$fileID, binarizeMethod = "naive", cores = 16)
var_genes = unique(unlist(plotTopMarkerHeat(sg, top_n = 500, outs = T, plotheat=F))) #the top 500 genes from each cluster (there will many overlaps as clusters will be shared between libraries)







#iterative clustering (you actually need to edit this according to how many libraries there are, and run from a local copy)
source("https://raw.githubusercontent.com/mahmoudibrahim/KidneyMap/master/templates/iterativeClustering.r")





#genesorteR (some marker plots and cluster quality checks
sg = sortGenes(logcounts(l), class_info, binarizeMethod = "naive", cores = 16)
zzz = pdf(paste0(outName, "_japanium_genesorteR_clusterSimilarity.pdf"))
corr_out = plotCorrelationHeat(sg, getMarkers(sg, quant=0.95)$markers, outs = TRUE)
dev.off()

pdf(paste0(outName,"_japanium_genesorteR_top10.pdf"), height = 12)
plotTopMarkerHeat(sg, top_n = 10, gaps = FALSE, averageCells=10^6, newOrder = corr_out$pheat$tree_col$order)
dev.off()

zzz = pdf(paste0(outName, "_japanium_genesorteR_clusterQuality.pdf"), height = 10)
barplot(getClassAUC(sg), ylim = c(0,1), yaxp = c(0,1,20))
zzz = dev.off()

zzz = pdf(paste0(outName, "_japanium_genesorteR_markerSets.pdf"), height = 10)
mmm = getMarkers(sg, quant = 0.95)
plotMarkerHeat(sg$inputMat, sg$inputClass, mmm$markers, clusterGenes=T, averageCells=10, newOrder = corr_out$pheat$tree_col$order)
dev.off()

zzz = pdf(paste0(outName, "_japanium_genesorteR_permuteMarkerSets.pdf"), height = 10)
pp = getPValues(sg, numPerm=20, cores=16)
pMark = names(which(apply(pp$adjpval, 1, function(x) any(x < 0.05))))
plotMarkerHeat(logcounts(l_bc), sg$inputClass, pMark, clusterGenes=T, averageCells=50, newOrder = corr_out$pheat$tree_col$order)
dev.off()

write.files(sg, prefix = paste0(outName, "_japanium_genesorteR"), markers = pMark)




#single plots of "important" genes
hend_markers = toupper(c("Pdgfrb", "Ptprc", "Pecam1", "Epcam", "Nrxn3", "MME","COL1A1","COL15A1","VIM","ACTA2","MYH11"))
names1 = toupper(unlist(lapply(strsplit(as.character(rownames(l)), split = ";", fixed = T), function(x) as.character(x[2]))))
for (i in 1:length(hend_markers)) {
	if(any(toupper(names1) == hend_markers[i])) {
		zzz = pdf(paste0(outName, "_markers_box_", hend_markers[i], ".pdf"))
		exp = logcounts(l)
		gn = which(toupper(names1) == hend_markers[i])
		plot(as.numeric(as.vector(class_info))+rnorm(ncol(exp), 0, 0.13), exp[gn,], col = makeTransparent(colors_info, alpha = 0.6), pch = 20, ylab = "Log Expression", xlab = "Cell Cluster", cex = 0.5, main = hend_markers[i])
		boxplot(exp[gn,] ~ class_info, col = makeTransparent("white", alpha = 0), border = makeTransparent("black", alpha = 0.5), range = 0.1, outline = F, add = T, notch = T, lwd = 3)
		dev.off()
	}
}





###OPTIONAL
#need to remove doublets (based on expression of major genes above) and clusters with no differential expression based on gene sorter (optional if these things exist). Example: from CD10- mesenchymal cell integration
wsd = which(class_info %in% c(2,7)) #cluster 2 is low quality cells (no diff. exp genes), cluster 7 is potential endothel./smc (mesenchymal doublet based on pecam1 and mesenchymal genes expression)
l = l[,-wsd]
class_info = class_info[-wsd]

#we will now start clustering again (if we removed clusters, otherwise we're done and skip to output)
varGenes = list()
meanClassAUC = list()
clusterLabels = list()
umapCoords = list()
graph = list()
rand = list()
mat = list()
file_id = list()
cellCluster = list()

#determine variable genes
sg = sortGenes(logcounts(l), class_info, cores = 16, binarizeMethod = "naive")
zzz = pdf(paste0(outName, "_japanium_genesorteR_clusterQuality-take", scorpions, ".pdf"), height = 10)
barplot(getClassAUC(sg), ylim = c(0,1), yaxp = c(0,1,20))
zzz = dev.off()
var_genes = unique(unlist(plotTopMarkerHeat(sg, top_n = 500, outs = T, plotheat=F)))

#normalize expression
l = scran::computeSumFactors(l, sizes = seq(10, 200, 20), positive = TRUE)
ws = which(sizeFactors(l) <= 0)
if(length(ws) > 0) {
	message("Warning: some cells had negative size factors!")
	sizeFactors(l)[ws] = min(sizeFactors(l)[sizeFactors(l) > 0])
}
zzz = pdf(paste0(outName, "_japanium_scran_sizeFactors_take", scorpions, ".pdf"))
plot(sizeFactors(l), l$total_counts/1e3, xlab = "Size factors", ylab = "Library Depth (thousands)", col = "#99999950", pch = 20)
zzz = dev.off()
l = scater::normalize(l)
l = scater::calculateQCMetrics(l)
a = plotExprsFreqVsMean(l)






#iterative clustering (you actually need to edit this according to how many libraries there are, and run from a local copy)
source("https://raw.githubusercontent.com/mahmoudibrahim/KidneyMap/master/templates/iterativeClustering.r")





#genesorteR (some marker plots and cluster quality checks
sg = sortGenes(logcounts(l), class_info, binarizeMethod = "naive", cores = 16)
zzz = pdf(paste0(outName, "_japanium_genesorteR_clusterSimilarity.pdf"))
corr_out = plotCorrelationHeat(sg, getMarkers(sg, quant=0.95)$markers, outs = TRUE)
dev.off()

pdf(paste0(outName,"_japanium_genesorteR_top10.pdf"), height = 12)
plotTopMarkerHeat(sg, top_n = 10, gaps = FALSE, averageCells=10^6, newOrder = corr_out$pheat$tree_col$order)
dev.off()

zzz = pdf(paste0(outName, "_japanium_genesorteR_clusterQuality.pdf"), height = 10)
barplot(getClassAUC(sg), ylim = c(0,1), yaxp = c(0,1,20))
zzz = dev.off()

zzz = pdf(paste0(outName, "_japanium_genesorteR_markerSets.pdf"), height = 10)
mmm = getMarkers(sg, quant = 0.95)
plotMarkerHeat(sg$inputMat, sg$inputClass, mmm$markers, clusterGenes=T, averageCells=10, newOrder = corr_out$pheat$tree_col$order)
dev.off()

zzz = pdf(paste0(outName, "_japanium_genesorteR_permuteMarkerSets.pdf"), height = 10)
pp = getPValues(sg, numPerm=20, cores=16)
pMark = names(which(apply(pp$adjpval, 1, function(x) any(x < 0.05))))
plotMarkerHeat(logcounts(l_bc), sg$inputClass, pMark, clusterGenes=T, averageCells=50, newOrder = corr_out$pheat$tree_col$order)
dev.off()

write.files(sg, prefix = paste0(outName, "_japanium_genesorteR"), markers = pMark)




#single plots of "important" genes
hend_markers = toupper(c("Pdgfrb", "Ptprc", "Pecam1", "Epcam", "Nrxn3", "MME","COL1A1","COL15A1","VIM","ACTA2","MYH11"))
names1 = toupper(unlist(lapply(strsplit(as.character(rownames(l)), split = ";", fixed = T), function(x) as.character(x[2]))))
for (i in 1:length(hend_markers)) {
	if(any(toupper(names1) == hend_markers[i])) {
		zzz = pdf(paste0(outName, "_markers_box_", hend_markers[i], ".pdf"))
		exp = logcounts(l)
		gn = which(toupper(names1) == hend_markers[i])
		plot(as.numeric(as.vector(class_info))+rnorm(ncol(exp), 0, 0.13), exp[gn,], col = makeTransparent(colors_info, alpha = 0.6), pch = 20, ylab = "Log Expression", xlab = "Cell Cluster", cex = 0.5, main = hend_markers[i])
		boxplot(exp[gn,] ~ class_info, col = makeTransparent("white", alpha = 0), border = makeTransparent("black", alpha = 0.5), range = 0.1, outline = F, add = T, notch = T, lwd = 3)
		dev.off()
	}
}





##messages and output objects that are then used in integrated maps
l = calculateQCMetrics(l)
a = plotExprsFreqVsMean(l)
aveRpG = ( (l$total_counts) / (l$total_features_by_counts) )

cat(paste0("There are ", ncol(l), " cells.\n"), file = paste0(outName, "_japanium.info.txt"))
cat(paste0("There are ", nrow(l), " expressed genes, with an average of ", round(mean(l$total_features_by_counts)), " expressed genes per cell.\n"), file = paste0(outName, "_japanium.info.txt"), append = T)
cat(paste0("There are on average ", round(mean(l$total_counts)), " UMIs per cell.\n"), file = paste0(outName, "_japanium.info.txt"), append = T)
cat(paste0("In each cell, there are on average ", round(mean(aveRpG),2), " transcripts per gene.\n"), file = paste0(outName, "_japanium.info.txt"), append = T)
cat(paste0("Timestamp: ", Sys.time()), file = paste0(outName, "_japanium.info.txt"), append = T)


##write data
lx = lx[,which(colnames(lx) %in% colnames(l))]
lx$colData$clusterLabels = class_info
save(lx = lx, file = paste0(outName, "-japanium.filtered_matrix.RData"))
write(wsd, file = paste0(outName, "_noNormaCells.txt"), ncolumns = 1)
write(colnames(l), file = paste0(outName, "_cellNames.txt"), ncolumns = 1)
write(paste(colnames(l), l$colData$fileID, l$colData$origCellCluster, class_info, sep = "\t"), file = paste0(outName, "_japanium_infoClusters", ".txt"), ncolumns = 1)


save(varGenes, meanClassAUC, umapCoords, graph, rand, file = paste0(outName, "_iterInfo.RData"))


##Fin
