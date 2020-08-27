########################################################################
# human_PDGFRBpositive.r
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

#custom functions
source("https://raw.githubusercontent.com/mahmoudibrahim/KidneyMap/master/.customlib.r")


#define color scheme
level1_col = c("#4363d8","#42006a","#800000","#9A6324","#808000")
level1_col_trans = makeTransparent(level1_col, alpha = 0.4)
level2_col_one = (colorRampPalette(c("#4363d8","white"))(n=6))[-6]
level2_col_one_trans = makeTransparent(level2_col_one, alpha = 0.9)
level2_col_two = (colorRampPalette(c("#42006a","white"))(n=6))[-6]
level2_col_two_trans = makeTransparent(level2_col_two, alpha = 0.9)
level2_col_three = (colorRampPalette(c("#800000","white"))(n=6))[-6]
level2_col_three_trans = makeTransparent(level2_col_three, alpha = 0.9)
level2_col_four = (colorRampPalette(c("#9A6324","white"))(n=6))[-6]
level2_col_four_trans = makeTransparent(level2_col_four, alpha = 0.9)
level2_col_five = c("#808000", "#d8d35f","#e7e16d")
level2_col_five_trans = makeTransparent(level2_col_five, alpha = 0.9)
kf_col = c("#999999","#fc8d62")
kf_col_trans = makeTransparent(kf_col, alpha = 0.4)



prefix = "human_PDGFRBpositive_"



###########################################################################
########################START ANALYSIS#####################################
###########################################################################

#load data
load("human_PDGFRBpositive.RData") #R object containing aggregated data, saved as SingleCellExperiment. Variable name: l
lx = l #save full matrix for later



#gene filtering
cutGenes = log10(floor(ncol(l) * 0.001)) 
choose_genes = which((log10(Matrix::rowSums(counts(l)))) >= cutGenes)
l = l[choose_genes, ]


#merge highly similar clusters (no merging happens here)
cutoff= 80 #if 80 genes common in the top 100 markers, merge the clusters
sg = sortGenes(counts(l), colData(l)$class_info, binarizeMethod="naive", cores = 16)
marks = plotTopMarkerHeat(sg, top_n=100, plotheat=FALSE, outs = TRUE)

xx = matrix(0, nrow = length(marks), ncol = length(marks))
for (i in 1:length(marks)) {
	for (j in 1:length(marks)) {
		xx[i,j] = length(intersect(marks[[i]], marks[[j]]))
	}
}
rownames(xx) = names(marks)
colnames(xx) = names(marks)
mergeMap = apply(xx, 2, function(x) names(which(x >= cutoff)))
mergethis = which(sapply(mergeMap,length) > 1)
mergeMap = mergeMap[mergethis]
class_info_merging = colData(l)$class_info
if (length(mergethis) > 0) {
	mergeMap = mergeMap[mergethis]
	class_info_merging = colData(l)$class_info
	for (i in 1:length(mergeMap)) {
		class_info_merging[which(class_info_merging %in% mergeMap[[i]])] = paste0("merge", i)
	}
} else {
	class_info_merging = colData(l)$class_info
}




#remove clusters with no differential gene expression
sg = sortGenes(counts(l), class_info_merging, binarizeMethod="naive", cores = 16)
removethis = names(which(apply(pp$adjpval,2,function(x) length(x[x<0.05])) == 0))

#find pdgfrb low clusters and remove them
pdgfrb_percent = sg$condGeneProb[which(rownames(sg$condGeneProb) == "ENSG00000113721.13;PDGFRB"),-(which(colnames(sg$condGeneProb) %in% removethis))] * 100
pdfgrb_cut = median(pdgfrb_percent) - (median(abs(pdgfrb_percent - median(pdgfrb_percent))) * 1)
pdfgrb_cut = names(which(pdgfrb_percent < pdfgrb_cut))

removethis = which(class_info_merging %in% removethis)
removethis = unique(c(removethis, which(class_info_merging %in% pdfgrb_cut)))




#reform matrix and refilter genes
class_info_merging_remove = class_info_merging[-removethis]
l_remove = lx
l_remove = l_remove[,-removethis]
l_remove = calculateQCMetrics(l_remove)
a = plotExprsFreqVsMean(l_remove)
colData(l_remove)$final_classes = class_info_merging_remove
cutGenes = log10(floor(ncol(l_remove) * 0.001)) 
choose_genes = which((log10(Matrix::rowSums(counts(l_remove)))) >= cutGenes)
l_remove = l_remove[choose_genes, ]





#make table of healthy vs disease
pdf(paste0(prefix, "percent_healthy_disease.pdf"))
barplot((table(l_remove$kidney_function) / ncol(l_remove)) * 100, col = rev(kf_col), ylim = c(0,70), yaxp = c(0,70,14), ylab = "% of Cells")
dev.off()

KFlow = rep(0,length(unique(class_info_merging_remove)))
for (i in 1:length(unique(class_info_merging_remove))) {
	where = which(class_info_merging_remove == unique(class_info_merging_remove)[i])
	tab = table(l_remove$kidney_function[where]) / length(where)
	if (length(tab) == 2) {
		KFlow[i] = tab[2]
	} else { #edit this otherwise
		KFlow[i] = tab
	}
}
names(KFlow) = unique(class_info_merging_remove)





#load cluster annotation info
annot = read.table("https://raw.githubusercontent.com/mahmoudibrahim/KidneyMap/master/assets/clusterInfo/human_PDGFRBpositive.txt", sep = "\t")
annot = annot[match(names(KFlow), annot$V1),]
annot$V7 = KFlow
annot$V8 = log(KFlow) - log((table(l_remove$kidney_function) / ncol(l_remove))[2])
nums = table(class_info_merging_remove)
annot = annot[match(names(nums), annot$V1),]
annot$V9 = as.numeric(nums)
annot$V10 = as.numeric(log10(nums))
colnames(annot) = c("ID", "Label", "Annotation 1", "Annotation 2", "Annotation 3", "Percent Low KF Cells", "Log Enrichment of Low KF Cells", "Number of Cells", "log10 Number of Cells")
##write cluster annotation information
write.csv(annot, paste0(prefix, "cluster_annotation_info_withAddData.csv"), row.names = FALSE)
#annotation levels
colData(l_remove)$level1 = l_remove$class_info
colData(l_remove)$level2 = l_remove$class_info
colData(l_remove)$level3 = l_remove$class_info
for (i in 1:length(unique(l_remove$final_classes))) {
	where = which(l_remove$final_classes == unique(l_remove$final_classes)[i])
	l_remove$level3[where] = as.character(annot[[5]][which(annot[[1]] == unique(l_remove$final_classes)[i])])
	l_remove$level2[where] = as.character(annot[[4]][which(annot[[1]] == unique(l_remove$final_classes)[i])])
	l_remove$level1[where] = as.character(annot[[3]][which(annot[[1]] == unique(l_remove$final_classes)[i])])
}





#get markers for different levels
sg_level1 = sortGenes(counts(l_remove), l_remove$level1, binarizeMethod = "naive", cores = 16)
write.files(sg_level1, prefix = paste0(prefix, "genesorteRouts_level1"))

sg_level2 = sortGenes(counts(l_remove), l_remove$level2, binarizeMethod = "naive", cores = 16)
write.files(sg_level2, prefix = paste0(prefix, "genesorteRouts_level2"))

sg_level3 = sortGenes(counts(l_remove), l_remove$level3, binarizeMethod = "naive", cores = 16)
write.files(sg_level3, prefix = paste0(prefix, "genesorteRouts_level3"))





#add colors
colors_info_level1 = rep("", length = ncol(l_remove))
for (i in 1:length(unique(l_remove$level1))) {
	wtemp = which(l_remove$level1 == (sort(unique(l_remove$level1)))[i])
	colors_info_level1[wtemp] = level1_col_trans[i]
}
colors_cluster_net = c("","")
colors_info_level2 = rep("", length = ncol(l_remove))
ww = which(l_remove$level1 == "Endothelial")
for (i in 1:length(unique(l_remove[,ww]$level2))) {
	wtemp = which(l_remove$level2 == (sort(unique(l_remove[,ww]$level2)))[i])
	colors_info_level2[wtemp] = level2_col_one_trans[i]
	colors_cluster_net = rbind(colors_cluster_net,c((sort(unique(l_remove[,ww]$level2)))[i],level2_col_one_trans[i]))
}
ww = which(l_remove$level1 == "Epithelial")
for (i in 1:length(unique(l_remove[,ww]$level2))) {
	wtemp = which(l_remove$level2 == (sort(unique(l_remove[,ww]$level2)))[i])
	colors_info_level2[wtemp] = level2_col_two_trans[i]
	colors_cluster_net = rbind(colors_cluster_net,c((sort(unique(l_remove[,ww]$level2)))[i],level2_col_two_trans[i]))
}
ww = which(l_remove$level1 == "Immune")
for (i in 1:length(unique(l_remove[,ww]$level2))) {
	wtemp = which(l_remove$level2 == (sort(unique(l_remove[,ww]$level2)))[i])
	colors_info_level2[wtemp] = level2_col_three_trans[i]
	colors_cluster_net = rbind(colors_cluster_net,c((sort(unique(l_remove[,ww]$level2)))[i],level2_col_three_trans[i]))
}
ww = which(l_remove$level1 == "Mesenchymal")
for (i in 1:length(unique(l_remove[,ww]$level2))) {
	wtemp = which(l_remove$level2 == (sort(unique(l_remove[,ww]$level2)))[i])
	colors_info_level2[wtemp] = level2_col_four_trans[i]
	colors_cluster_net = rbind(colors_cluster_net,c((sort(unique(l_remove[,ww]$level2)))[i],level2_col_four_trans[i]))
}
ww = which(l_remove$level1 == "Neuronal")
for (i in 1:length(unique(l_remove[,ww]$level2))) {
	wtemp = which(l_remove$level2 == (sort(unique(l_remove[,ww]$level2)))[i])
	colors_info_level2[wtemp] = level2_col_five_trans[i]
	colors_cluster_net = rbind(colors_cluster_net,c((sort(unique(l_remove[,ww]$level2)))[i],level2_col_five_trans[i]))
}
colors_cluster_net = colors_cluster_net[-1,]

colors_info_kf = rep("", length = ncol(l_remove))
for (i in 1:length(unique(l_remove$kidney_function))) {
	wtemp = which(l_remove$kidney_function == (sort(unique(l_remove$kidney_function)))[i])
	colors_info_kf[wtemp] = rev(kf_col_trans)[i]
}
colors_file = rep("", length = ncol(l_remove))
for (i in 1:length(unique(l_remove$fileID))) {
	wtemp = which((l_remove$fileID) == (sort(unique(l_remove$fileID)))[i])
	colors_file[wtemp] = col_palette_trans[i]
}

library(RColorBrewer)
nb.cols = length(unique(l_remove$level3))
mycolors = colorRampPalette(brewer.pal(30, "Set2"))(nb.cols)
colors_info_level3 = rep("", length = ncol(l_remove))
for (i in 1:length(unique(l_remove$level3))) {
	wtemp = which((l_remove$level3) == (sort(unique(l_remove$level3)))[i])
	colors_info_level3[wtemp] = makeTransparent(mycolors[i])
}



#dimension reduction
sg = sortGenes(counts(l_remove), class_info_merging_remove, binarizeMethod = "naive", cores = 16)
pp = getPValues(sg, numPerm = 20, cores = 16)
l_remove = scran::computeSumFactors(l_remove, sizes = seq(10, 200, 20), clusters = class_info_merging_remove, positive = TRUE)
l_remove = scater::normalize(l_remove)
l_remove = scater::calculateQCMetrics(l_remove)
a = plotExprsFreqVsMean(l_remove)
var_genes = names(which(apply(pp$adjpval, 1, function(x) any(x < 0.05))))

mnn = batchelor::fastMNN(l_remove, batch = l_remove$fileID, d = 30, auto.order = TRUE, subset.row = var_genes, BSPARAM = BiocSingular::IrlbaParam(tol=1e-8))
mnn_u = mnn@reducedDims$corrected
temp = t(apply(mnn_u, 1, function(x) x / (sqrt(sum(x^2))))) #norm based scaling

write.csv(temp, file = paste0(prefix, "dims.txt"), row.names = FALSE, col.names = FALSE) #file saved @https://raw.githubusercontent.com/mahmoudibrahim/KidneyMap/master/assets/reducedDims/SVD/human_PDGFRBpositive_dims.txt




#UMAP (produced in python, just for visualization, independent of clustering)
mapCoords = read.table("https://raw.githubusercontent.com/mahmoudibrahim/KidneyMap/master/assets/reducedDims/UMAP/human_PDGFRBpositive_umapCoords.csv", header = FALSE, sep=",")

pdf(paste0(prefix, "UMAP_level1_level3TEXT.pdf"))
plot(mapCoords[,1], mapCoords[,2], pch = 20, col = colors_info_level1, cex = 0.3)
allClasses = unique(l_remove$level3)
for (i in 1:length(allClasses)) {
	ww = which(l_remove$level3 == allClasses[i])
	xx = median(mapCoords[ww,1])
	yy = median(mapCoords[ww,2])
	text(xx,yy,allClasses[i], cex = 0.4, adj = 0.5)
}
dev.off()
pdf(paste0(prefix, "UMAP_KF.pdf"))
plot(mapCoords[,1], mapCoords[,2], pch = 20, col = colors_info_kf, cex = 0.3, xlim = c(-7,15))
dev.off()
pdf(paste0(prefix, "UMAP_fileID.pdf"))
plot(mapCoords[,1], mapCoords[,2], pch = 20, col = colors_file, cex = 0.3)
legend(min(mapCoords[,1]) + 1, min(mapCoords[,2]) + 15, legend = abbreviate(unique(l_remove$fileID), minlength=6), pch=21, col = col_palette, pt.bg = col_palette_trans, pt.cex=2, cex=.8, bty="n", ncol = ceiling(length(unique(l_remove$fileID)) / 10))
dev.off()
pdf(paste0(prefix, "UMAP_clusterColors.pdf"))
plot(mapCoords, col = colors_info_level3, pch = 20, cex = 0.3)
legend(max(mapCoords[,1]) - 10, max(mapCoords[,2]), legend = abbreviate(sort(unique(l_remove$level3)), minlength=6), pch=21, col = mycolors, pt.bg = makeTransparent(mycolors), pt.cex=2, cex=.8, bty="n", ncol = ceiling(length(unique(l_remove$level3)) / 20))
dev.off()








#core matrisome score
l_remove$level3KF = paste(l_remove$level3, l_remove$kidney_function)
l_remove$level2KF = paste(l_remove$level2, l_remove$kidney_function)
l_remove$level1KF = paste(l_remove$level1, l_remove$kidney_function)

matrisome_set = read.table("https://raw.githubusercontent.com/mahmoudibrahim/KidneyMap/master/assets/public/ecm_genes_human.txt", sep = "\t", header = TRUE)
collagens = matrisome_set$Gene.Symbol[which(matrisome_set$Division %in% c("Core matrisome"))]
reads_single_phase = logcounts(l_remove)
rownames(reads_single_phase) = toupper(unlist(lapply(strsplit(as.character(rownames(reads_single_phase)), split = ";", fixed = T), function(x) as.character(x[2]))))
reads_single_phase = as.matrix(reads_single_phase[rownames(reads_single_phase) %in% toupper(collagens) ,])
combined_matrix = rbind(reads_single_phase,average=apply(reads_single_phase,2,mean))
cor_matrix = cor(t(combined_matrix))
cor_vector = cor_matrix[,dim(cor_matrix)[1]]
reads_single_phase_restricted = reads_single_phase[rownames(reads_single_phase) %in% names(cor_vector[cor_vector >= 0.1]),]
aaa = apply(reads_single_phase_restricted,2,mean)


pdf(paste0(prefix, "_core_matrisome_score_level1.pdf"))
size = scale(sqrt(table(l_remove$level1)),center = FALSE) * 0.8
vioplot(aaa ~ l_remove$level1, col = level1_col, ylim = c(0,2.2), plotCentre="line", wex = size, lwd = 3, border = level1_col)
dev.off()

pdf(paste0(prefix, "_core_matrisome_score_level2_mesenchymal.pdf"))
size = scale(sqrt(table(l_remove$level2[which(l_remove$level1 == "Mesenchymal")])),center = FALSE)
vioplot(aaa[which(l_remove$level1 == "Mesenchymal")] ~ l_remove$level2[which(l_remove$level1 == "Mesenchymal")], col = level1_col[4], ylim = c(0,2.2), plotCentre="line", wex = size, lwd = 3, border = level1_col[4], names = abbreviate(rownames(size)))
dev.off()

pdf(paste0(prefix, "_core_matrisome_score_level2_epithelial.pdf"))
size = scale(sqrt(table(l_remove$level2[which(l_remove$level1 == "Epithelial")])),center = FALSE)
vioplot(aaa[which((l_remove$level1 == "Epithelial"))] ~ l_remove$level2[which((l_remove$level1 == "Epithelial"))], col = level1_col[2], ylim = c(0,2.2), plotCentre="line", wex = size, lwd = 3, border = level1_col[2])
dev.off()

zzz = pdf(paste0(prefix, "umap_ecm_continuous.pdf"))
mapCoords = read.table("https://raw.githubusercontent.com/mahmoudibrahim/KidneyMap/master/assets/reducedDims/UMAP/human_PDGFRBpositive_umapCoords.csv", header = FALSE, sep=",")
plot(mapCoords, col = color.gradient(scale(aaa), colors = c("#ECECEC50",rev(c("#d7191c","#abdda4")))), pch = 20, cex = 1)
legend.col(col = colorRampPalette(c("#ECECEC50",rev(c("#d7191c","#abdda4"))))(100), lev = seq(min(aaa), max(aaa), length.out=100))
zzz = dev.off()


#marker heatmaps
level3_order = c(1,2,9,10,24,26,13,3,22,23,12,4,15,11,14,25,21,7,8,5,6,17,18,19,20,16)
pdf(paste0(prefix, "level3_top_auto_markers_5.pdf"), height = 10)
plotTopMarkerHeat(sg_level3, averageCells=10^2, newOrder=level3_order, colors = colorRampPalette(c("cornflowerblue", "black", "gold"))(n=100), gaps = TRUE, top_n=5)
dev.off()



#doublet scores
pdf(paste0(prefix, "doubletScore_box.pdf"), width = 40, height = 10)
exp = l_remove$doublet_score
facts = factor(as.factor(l_remove$level3),levels(as.factor(l_remove$level3))[level3_order])
plot(as.numeric(facts)+rnorm(length(exp), 0, 0.13), exp, pch = 20, xlab = "", ylab = "", cex = 0.3, main = "Doublet Score", col = "#99999920", xaxt = "n")
axis(1, at=sort(unique(as.numeric(facts))), labels=abbreviate(levels(facts)))
vioplot(as.vector(exp) ~ facts, col = makeTransparent("white", alpha = 0), border = makeTransparent("black", alpha = 0.3), range = 0.1, outline = F, add = T, notch = T, lwd = 4,  plotCentre="line", xlab = "Cell Clusters", ylab = "Doublet Score")
dev.off()





#############################################################
####################focus on mesenchymal cells###############
#############################################################

ww = c(which(l_remove$level2 %in% c("Fibroblasts", "Myofibroblasts","Pericytes")))
l_v = l_remove[,ww]
l_v = scran::computeSumFactors(l_v, sizes = seq(10, 200, 20), clusters = l_v$level3, positive = TRUE)
l_v = scater::normalize(l_v)

sg_mesen = sortGenes(counts(l_v), l_v$level3, binarizeMethod="naive", cores = 16)
pp_mesen = getPValues(sg_mesen, numPerm = 20, cores = 16)
marker_mesen = names(which(apply(pp_mesen$adjpval, 1, function(x) any(x < 0.05))))
sg_mesen_l2 = sortGenes(logcounts(l_v), l_v$level2, binarizeMethod="naive", cores = 16)
pp_mesen_l2 = getPValues(sg_mesen_l2, numPerm = 20, cores = 16)
marker_mesen_l2 = names(which(apply(pp_mesen_l2$adjpval, 1, function(x) any(x < 0.05))))


#dimension reduction for mesenchymal (fibroblast-related) cells
mnn = batchelor::fastMNN(l_v, batch = l_v$fileID, d = 30, auto.order = TRUE, subset.row = marker_mesen, BSPARAM = BiocSingular::IrlbaParam(tol=1e-8))
mnn_u = mnn@reducedDims$corrected
temp = t(apply(mnn_u, 1, function(x) x / (sqrt(sum(x^2)))))
write.csv(temp, file = paste0(prefix, "dims_pseudotime.txt"), row.names = FALSE, col.names = FALSE)





#marker heatmap fibro-zoom
sg_lv_level3 = sortGenes(logcounts(l_v), l_v$level3, binarizeMethod="median")
pdf(paste0(prefix, "fibrozoom_marker_heat.pdf"))
plotTopMarkerHeat(sg_lv_level3, top_n=5, averageCells=10^2, newOrder=c(10,1,2,3,4,8,9,6,7,5), colors = colorRampPalette(c("cornflowerblue", "black", "gold"))(n=100))
dev.off()


#matrisome score
matrisome_set = read.table("https://raw.githubusercontent.com/mahmoudibrahim/KidneyMap/master/assets/public/ecm_genes_human.txt", sep = "\t", header = TRUE)
collagens = matrisome_set$Gene.Symbol[which(matrisome_set$Division %in% c("Core matrisome"))]
reads_single_phase = logcounts(l_remove)
rownames(reads_single_phase) = toupper(unlist(lapply(strsplit(as.character(rownames(reads_single_phase)), split = ";", fixed = T), function(x) as.character(x[2]))))
reads_single_phase = as.matrix(reads_single_phase[rownames(reads_single_phase) %in% toupper(collagens) ,])
combined_matrix = rbind(reads_single_phase,average=apply(reads_single_phase,2,mean))
cor_matrix = cor(t(combined_matrix))
cor_vector = cor_matrix[,dim(cor_matrix)[1]]
reads_single_phase_restricted = reads_single_phase[rownames(reads_single_phase) %in% names(cor_vector[cor_vector >= 0.1]),]
aaa = apply(reads_single_phase_restricted,2,mean)
ecm = aaa[ww]


#UMAP-based
mapCoords = read.table("https://raw.githubusercontent.com/mahmoudibrahim/KidneyMap/master/assets/reducedDims/UMAP/human_PDGFRBpositive_mesenchymal_umapCoords.csv", header = FALSE, sep=",") #generated in python

#set the colors
colors_info_mesenl3 = rep("", length = ncol(l_v))
for (i in 1:length(unique(l_v$level3))) {
	wtemp = which(l_v$level3 == (sort(unique(l_v$level3)))[i])
	colors_info_mesenl3[wtemp] = col_palette_trans[i]
}
colors_info_mesenl2 = rep("", length = ncol(l_v))
for (i in 1:length(unique(l_v$level2))) {
	wtemp = which(l_v$level2 == (sort(unique(l_v$level2)))[i])
	colors_info_mesenl2[wtemp] = col_palette_short_trans[c(2,4,1)][i]
}
colors_info_patient = rep("", length = ncol(l_v))
for (i in 1:length(unique(l_v$fileID))) {
	wtemp = which((l_v$fileID) == (sort(unique(l_v$fileID)))[i])
	colors_info_patient[wtemp] = col_palette_trans[i]
}


pdf(paste0(prefix, "pseudotime_umap_level3.pdf"))
plot(mapCoords, col = colors_info_mesenl3, pch = 20)
legend(min(mapCoords[,1]), max(mapCoords[,2]), legend = sort(unique(l_v$level3)), pch=21, col = col_palette_trans, pt.bg = col_palette, pt.cex=2, cex=.8, bty="n", ncol = ceiling(length(unique(l_v$level3)) / 20))
dev.off()


pdf(paste0(prefix, "pseudotime_specific_umap_level2.pdf"))
plot(mapCoords, col = colors_info_mesenl2, pch = 20)
legend(min(mapCoords[,1]), max(mapCoords[,2]), legend = sort(unique(l_v$level2)), pch=21, col = col_palette_short_trans[c(2,4,1)], pt.bg = col_palette_short[c(2,4,1)], pt.cex=2, cex=.8, bty="n", ncol = ceiling(length(unique(l_v$level2)) / 4))
dev.off()

pdf(paste0(prefix, "pseudotime_specific_umap_ECMscore.pdf"))
plot(mapCoords, col = color.gradient(scale(aaa[ww]), colors = c("#ECECEC50",rev(c("#d7191c","#abdda4")))), pch = 20)
dev.off()


hend_markers = toupper(c("Pdgfrb", "Col1a1", "Col1a2", "Vim","Dcn","notch3","scara5","rspo3","pcolce2","postn","pdgfra","ogn","CCL19","adh1b","nkd2","col15a1", "CLMP", "ENTPD3", "SLC26A7","RGS5","NDUFA4L2","HIGD1B","GJA4","COX4I2","MCAM","CPE","STEAP4","TM4SF1","GJC1","ADIRF","NOTCH3","IMPA2","PMEPA1","KLHL23","ADGRF5","PLAC9","ISYNA1","PTP4A3","ETV1","NKD2","IGFBP5","COL1A1","TF","IGF1R","IGF1R"))
names1 = toupper(unlist(lapply(strsplit(as.character(rownames(l_remove)), split = ";", fixed = T), function(x) as.character(x[2]))))
for (i in 1:length(markersp)) {
	if(any(names1 == markersp[i])) {
		aaax = as.vector(logcounts(l_remove)[which(names1 == markersp[i]),])
		aaax = aaax[ww]
		zzz = pdf(paste0(prefix, "pseudotime_specific_umap_", markersp[i], ".pdf"))
		plot(mapCoords, col = color.gradient(scale(aaax), colors = c("#ECECEC50",rev(c("#d7191c","#abdda4")))), pch = 20, cex = 1)
		legend.col(col = colorRampPalette(c("#ECECEC50",rev(c("#d7191c","#abdda4"))))(100), lev = seq(min(aaax), max(aaax), length.out=100))
		zzz = dev.off()
	}
}



library(slingshot)
sl1 = slingshot(mapCoords, l_v$level3, start.clus = "Myofibroblasts 1")
st1 = rowMeans(slingPseudotime(sl1), na.rm = TRUE)
pdf(paste0(prefix, "pseudotime_umap_pseudotime.pdf"))
col_pst = color.gradient(scale(-st1), colors = c("#ECECEC50",rev(c("#d7191c","#abdda4"))))
plot(mapCoords, col = col_pst, pch = 20)
dev.off()
pdf(paste0(prefix, "pseudotime_umap_level3_withLineageTree.pdf"))
plot(mapCoords, col = colors_info_mesenl3, pch = 20)
legend(min(mapCoords[,1]), min(mapCoords[,2]) + 6, legend = sort(unique(l_v$level3)), pch=21, col = col_palette_trans, pt.bg = col_palette, pt.cex=2, cex=.8, bty="n", ncol = ceiling(length(unique(l_v$level3)) / 20))
lines(sl1, lwd=5, col="#00000095", type = "lineages")
dev.off()




#diffusion map based
library(destiny) 
set.seed(111)
df = DiffusionMap(mnn_u, k = 80, n_eigs = 3)
write.csv(df@eigenvectors, file = paste0(prefix, "dims_diffusionmap.txt"), row.names = FALSE, col.names = FALSE) #file now saved @https://raw.githubusercontent.com/mahmoudibrahim/KidneyMap/master/assets/reducedDims/DM/human_PDGFRBpositive_dm.txt

#get cell ordering
df2 = as.matrix(read.csv("https://raw.githubusercontent.com/mahmoudibrahim/KidneyMap/master/assets/reducedDims/DM/human_PDGFRBpositive_dm.txt"))
sl2 = slingshot(df2, l_v$level3, start.clus = "Myofibroblasts 1")
st2 = rowMeans(slingPseudotime(sl2), na.rm = TRUE)
col_pst = color.gradient(scale(-st2), colors = c("#ECECEC50",rev(c("#d7191c","#abdda4"))))
curv = getCurves(sl2)


library(plot3D)
pdf(paste0(prefix, "pseudotime_specific_dm_level3.pdf"))
scatter3D(df2[,1], df2[,2], df2[,3], col = colors_info_mesenl3, pch = 20, cex = 0.5, colvar = NULL, phi = 120, theta = 230, bty = "n")
dev.off()
pdf(paste0(prefix, "pseudotime_specific_dm_level3_bty.pdf"))
scatter3D(df2[,1], df2[,2], df2[,3], col = colors_info_mesenl3, pch = 20, cex = 0.5, colvar = NULL, phi = 120, theta = 230) #90,310
dev.off()
pdf(paste0(prefix, "pseudotime_specific_dm_level2.pdf"))
scatter3D(df2[,1], df2[,2], df2[,3], col = colors_info_mesenl2, pch = 20, cex = 0.5, colvar = NULL, bty = "n", phi = 120, theta = 230)
dev.off()
pdf(paste0(prefix, "pseudotime_specific_dm_pseudotime.pdf"))
scatter3D(df2[,1], df2[,2], df2[,3], col = col_pst, pch = 19, cex = 0.8, colvar = NULL, bty = "n", phi = 120, theta = 230)
dev.off()
pdf(paste0(prefix, "pseudotime_specific_dm_ecm.pdf"))
scatter3D(df2[,1], df2[,2], df2[,3], col = color.gradient(scale(aaa[ww]), colors = c("#ECECEC50",rev(c("#d7191c","#abdda4")))), pch = 19, cex = 1, colvar = NULL, bty = "n", phi = 120, theta = 230)
dev.off()
pdf(paste0(prefix, "pseudotime_specific_dm_kidneyFunction.pdf"))
scatter3D(df2[,1], df2[,2], df2[,3], col = colors_info_kf[ww], pch = 19, cex = 0.8, colvar = NULL, bty = "n", phi = 120, theta = 230)
dev.off()
tiff(paste0(prefix,"pseudotime_specific_dm_lines.tiff"), res=300, height = 12, width = 12, units = "in")
points3D(df2[,1], df2[,2], df2[,3], col = "#99999901", pch = 20, cex = 0.1, colvar = NULL, bty = "n", phi = 120, theta = 230, alpha = 0)
lines3D(curv@curves$curve1$s[curv@curves$curve1$ord,1], curv@curves$curve1$s[curv@curves$curve1$ord,2], curv@curves$curve1$s[curv@curves$curve1$ord,3],  phi = 120, theta = 230, colvar = NULL, col = "#DE5246", lwd = 4, add = TRUE, alpha = 0.5)
lines3D(curv@curves$curve2$s[curv@curves$curve2$ord,1], curv@curves$curve2$s[curv@curves$curve2$ord,2], curv@curves$curve2$s[curv@curves$curve2$ord,3],  phi = 120, theta = 230, colvar = NULL, col = "#DE5246", lwd = 4, add = TRUE, alpha = 0.5)
lines3D(curv@curves$curve3$s[curv@curves$curve3$ord,1], curv@curves$curve3$s[curv@curves$curve3$ord,2], curv@curves$curve3$s[curv@curves$curve3$ord,3],  phi = 120, theta = 230, colvar = NULL, col = "#DE5246", lwd = 4, add = TRUE, alpha = 0.5)
dev.off()


#plot some markers
markersp = toupper(c("Pdgfrb", "Col1a1", "Col1a2", "Vim","Dcn","notch3","scara5","rspo3","pcolce2","postn","pdgfra","ogn","CCL19","adh1b","nkd2","col15a1", "CLMP", "ENTPD3", "SLC26A7","col14a1"))
names1 = toupper(unlist(lapply(strsplit(as.character(rownames(l_remove)), split = ";", fixed = T), function(x) as.character(x[2]))))
for (i in 1:length(markersp)) {
	if(any(names1 == markersp[i])) {
		aaax = as.vector(logcounts(l_remove)[which(names1 == markersp[i]),])
		aaax = aaax[ww]
		zzz = pdf(paste0(prefix, "pseudotime_specific_dm_", markersp[i], ".pdf"))
		scatter3D(df2[,1], df2[,2], -df2[,3], col = color.gradient(scale(aaax), colors = c("#ECECEC50",rev(c("#d7191c","#abdda4")))), pch = 19, cex = 1, colvar = NULL, bty = "n", phi = 0, theta = 310, alpha = 1)
		zzz = dev.off()
	}
}




#get genes that are correlated with cell ordering for each lineage
st = slingPseudotime(sl2)
time_test = list()
for (i in 1:3) {
	negSel = which(!is.na(st[,i]))
	set.seed(111)
	time_test1 = apply(logcounts(l_v[which(rownames(l_v) %in% marker_mesen),negSel]), 1, function(x) if (length(x[x!=0]) > 5) {cor.test(x, st[negSel,i], method = "spearman", exact = FALSE)$p.value} else {1})
	time_test[[i]] = which(p.adjust(time_test1, method = "BH") < 0.001)
}


#plot gene heatmap and gene clustering for lineage 1
tiff(paste0(prefix,"pseudotime_heat_L1.tiff"), res=300, height = 12, width = 12, units = "in")
negSel = which(!is.na(st[,1]))
set.seed(111)
kk = Mclust(-st[negSel,1], 5, modelNames="E")
roses = logcounts(l_v[,negSel])
l1 = plotMarkerHeat(roses[,order(st[negSel,1], decreasing = TRUE)], kk$classification[order(st[negSel,1], decreasing = TRUE)], names(time_test[[1]]), averageCells=10^1, clusterGenes = TRUE, clusterGenesK = 7, outs = TRUE, plotheat = FALSE) #cluster the genes
gene_order = l1$gene_class_info #reorder gene clusters (manually for now)
gene_order[l1$gene_class_info == 1] = 7
gene_order[l1$gene_class_info == 2] = 2
gene_order[l1$gene_class_info == 3] = 4
gene_order[l1$gene_class_info == 4] = 1
gene_order[l1$gene_class_info == 5] = 3
gene_order[l1$gene_class_info == 6] = 6
gene_order[l1$gene_class_info == 7] = 5
l1 = plotMarkerHeat(roses[,order(st[negSel,1], decreasing = TRUE)], kk$classification[order(st[negSel,1], decreasing = TRUE)], names(sort(gene_order)), averageCells=10^1, clusterGenes = FALSE, gaps = FALSE, colors = colorRampPalette(c("cornflowerblue","black","gold"))(n=100)) #make the heatmap with the new gene order
go1 = gene_order
dev.off()


#plot gene heatmap and gene clustering for lineage 2
tiff(paste0(prefix,"pseudotime_heat_L2.tiff"), res=300, height = 12, width = 12, units = "in")
negSel = which(!is.na(st[,2]))
set.seed(111)
kk = Mclust(-st[negSel,2], 6, modelNames="E")
roses = logcounts(l_v[,negSel])
l2 = plotMarkerHeat(roses[,order(st[negSel,2], decreasing = TRUE)], kk$classification[order(st[negSel,2], decreasing = TRUE)], names(time_test[[2]]), averageCells=10^1, clusterGenes = TRUE, clusterGenesK = 7, outs = TRUE, plotheat = FALSE)
gene_order = l2$gene_class_info
gene_order[l2$gene_class_info == 1] = 3
gene_order[l2$gene_class_info == 2] = 4
gene_order[l2$gene_class_info == 3] = 2
gene_order[l2$gene_class_info == 4] = 6
gene_order[l2$gene_class_info == 5] = 1
gene_order[l2$gene_class_info == 6] = 7
gene_order[l2$gene_class_info == 7] = 5
l2 = plotMarkerHeat(roses[,order(st[negSel,2], decreasing = TRUE)], kk$classification[order(st[negSel,2], decreasing = TRUE)], names(sort(gene_order)), averageCells=10^1, clusterGenes = FALSE, gaps = FALSE, colors = colorRampPalette(c("cornflowerblue","black","gold"))(n=100))
go2 = gene_order
dev.off()


#plot gene heatmap and gene clustering for lineage 3
tiff(paste0(prefix,"pseudotime_heat_L3.tiff"), res=300, height = 12, width = 12, units = "in")
negSel = which(!is.na(st[,3]))
set.seed(111)
kk = Mclust(-st[negSel,3], 5, modelNames="E")
roses = logcounts(l_v[,negSel])
l3 = plotMarkerHeat(roses[,order(st[negSel,3], decreasing = TRUE)], kk$classification[order(st[negSel,3], decreasing = TRUE)], names(time_test[[3]]), averageCells=10^1, clusterGenes = TRUE, clusterGenesK = 7, outs = TRUE, plotheat = FALSE)
gene_order = l3$gene_class_info
gene_order[l3$gene_class_info == 1] = 3
gene_order[l3$gene_class_info == 2] = 6
gene_order[l3$gene_class_info == 3] = 4
gene_order[l3$gene_class_info == 4] = 5
gene_order[l3$gene_class_info == 5] = 7
gene_order[l3$gene_class_info == 6] = 1
gene_order[l3$gene_class_info == 7] = 2
l3 = plotMarkerHeat(roses[,order(st[negSel,3], decreasing = TRUE)], kk$classification[order(st[negSel,3], decreasing = TRUE)], names(sort(gene_order)), averageCells=10^1, clusterGenes = FALSE, gaps = FALSE, colors = colorRampPalette(c("cornflowerblue","black","gold"))(n=100))
go3 = gene_order
dev.off()

#export gene clustering information
write.csv(data.frame(gene = names(go1), cluster = go1), file = paste0(prefix,"gene_cluster_info_L1.txt"), quote = FALSE, row.names = FALSE) 
write.csv(data.frame(gene = names(go2), cluster = go2), file = paste0(prefix,"gene_cluster_info_L2.txt"), quote = FALSE, row.names = FALSE)
write.csv(data.frame(gene = names(go3), cluster = go3), file = paste0(prefix,"gene_cluster_info_L3.txt"), quote = FALSE, row.names = FALSE)



#cell cycle & GO along cell ordering
library(tibble)
library(clusterProfiler)
G1_S = read.table("https://raw.githubusercontent.com/mahmoudibrahim/KidneyMap/master/assets/public/cell_cycle/G1S_human")[[1]]
G2_M = read.table("https://raw.githubusercontent.com/mahmoudibrahim/KidneyMap/master/assets/public/cell_cycle/G2M_human")[[1]]
G0 = read.table("https://raw.githubusercontent.com/mahmoudibrahim/KidneyMap/master/assets/public/cell_cycle/G0_human")[[1]]
cycle_list = list(G1_S, G2_M, G0)
flexible_normalization <- function(data_in,by_row=TRUE){
  if(by_row){
    row_mean <- apply(data_in,1,mean)
    row_sd   <- apply(data_in,1,sd)
    output <- data_in
    for(i in 1:dim(data_in)[1]){
      output[i,] <- (data_in[i,] - row_mean[i])/row_sd[i]
    }
  }
  if(!by_row){
    col_mean <- apply(data_in,2,mean)
    col_sd   <- apply(data_in,2,sd)
    output <- data_in
    for(i in 1:dim(data_in)[2]){
      output[,i] <- (data_in[,i] - col_mean[i])/col_sd[i]
    }
  }
  output
}



#Lineage 1 (cell cycle)
negSel = which(!is.na(st[,1]))
set.seed(111)
kk = Mclust(-st[negSel,1], 5, modelNames="E")
roses = logcounts(l_v[,negSel])
roses = roses[,order(st[negSel,1], decreasing = TRUE)]
names1 = toupper(unlist(lapply(strsplit(as.character(rownames(roses)), split = ";", fixed = T), function(x) as.character(x[2]))))
rownames(roses) = names1
reads_single = as.matrix(roses)
ans =
lapply(cycle_list,function(xx){
  reads_single_phase <- reads_single[rownames(reads_single) %in% toupper(unlist(xx)) ,]
  combined_matrix <- rbind(reads_single_phase,average=apply(reads_single_phase,2,mean))
  cor_matrix <- cor(t(combined_matrix))
  cor_vector <- cor_matrix[,dim(cor_matrix)[1]]
  reads_single_phase_restricted <- reads_single_phase[rownames(reads_single_phase) %in% names(cor_vector[cor_vector >= 0.1]),]
  apply(reads_single_phase_restricted,2,mean)
})
ans = matrix(unlist(ans), ncol = 3, byrow = FALSE)
colnames(ans) = c("G1_S", "G2_M", "G0")
ans_normed <- flexible_normalization(ans,by_row=FALSE)
ans_normed_normed <- flexible_normalization(ans_normed,by_row=TRUE)
cell_phase = apply(ans_normed_normed,1,function(x) colnames(ans_normed_normed)[which.max(x)])
n1 = seq(1,length(cell_phase), by = 2000)
n2 = n1 + 2000
n2[length(n2)] = length(cell_phase)
n3 = n2
for (i in 1:length(n2)) {
	n3[(n1[i]):(n2[i])] = i
}
n3[n3 == 12] = 11 #avoid having only a few cells at the end of the chain
ballet = lapply(1:length(n1), function(x) table(cell_phase[(n1[x]):(n2[x])]) / length((n1[x]):(n2[x])))
ballet = (matrix(unlist(ballet), byrow = TRUE, ncol = 3)) * 100
colnames(ballet) = c("G0", "G1_S", "G2_M")
pdf(paste0(prefix, "L1_cellcyle.pdf"))
plot(n1, ballet[,1], type = "l", ylim = c(15,60), col = makeTransparent(col_palette_short[1], 0.9), lwd = 5, yaxp = c(15,60,9), xlab = "Cells ordered by pseudotime (Fib4 to MF1)", ylab = "% of Cells")
lines(n1, ballet[,2], type = "l", col = makeTransparent(col_palette_short[2], 0.9), lwd = 5)
lines(n1, ballet[,3], type = "l", col = makeTransparent(col_palette_short[4], 0.9), lwd = 5)
legend("topright", legend = colnames(ballet), pch=15, col = col_palette_short[c(1,2,4)], bty = "n", cex = 1.5)
dev.off()
#gene ontology
sg_n3  = sortGenes(roses, n3, binarizeMethod = "naive")
n3_l = plotTopMarkerHeat(sg_n3, averageCells=10^2, top_n=100, plotheat=FALSE, outs = TRUE)
m_t2g =  as.tibble(read.gmt("/home/mibrahim/Dropbox/dev/kidneyMap/public/assets/public/c2.cp.pid.v7.0.symbols.gmt"))
en = mclapply(n3_l, function(x) enricher(x, TERM2GENE=m_t2g, minGSSize = 10, maxGSSize = 200, pvalueCutoff = 0.05), mc.cores = 12)
en2 = list()
for (i in 1:length(en)) {
	en2[[i]] =  en[[i]]@result$ID[1:5]
}
en2 = unique(unlist(en2))
en3 = list()
for (i in 1:length(en)) {
	en3[[i]] =  cbind(i, en[[i]]@result$ID[which(en[[i]]@result$ID %in% en2)], -log10(as.numeric(en[[i]]@result$qvalue[which(en[[i]]@result$ID %in% en2)])))
}
en3 = data.frame(do.call(rbind, en3))
en3 = data.frame(en3$i, en3$V2, en3$V3)
colnames(en3) = c("Cluster", "Term", "PValue")
en3 = reshape(en3, idvar = "Cluster", timevar = "Term", direction = "wide")
en4 = data.frame(en3[2:length(en3)])
en4 = as.matrix(en4)
en4[which(is.na(en4))] = "0"
en4 = apply(en4, 2, as.numeric)
en4[which(en4 > -log10(0.0001))] = -log10(0.0001)
rownames(en4) = en3[[1]]
pdf(paste0(prefix, "PID_Lineage1.pdf"), width = 10, height = 10)
heatmap.2(t(en4), trace = "none", scale = "none", col = colorRampPalette(c("white","dodgerblue4"))(n=100), margins = c(20,25), density = "none", Colv = NA, key = FALSE, labCol = names(n3_l))
dev.off()
sg_n3  = sortGenes(roses, n3, binarizeMethod = "naive")
n3_l = plotTopMarkerHeat(sg_n3, averageCells=10^2, top_n=100, plotheat=FALSE, outs = TRUE)
m_t2g =  as.tibble(read.gmt("/home/mibrahim/Dropbox/dev/kidneyMap/public/assets/public/c2.cp.kegg.v7.0.symbols.gmt"))
en = mclapply(n3_l, function(x) enricher(x, TERM2GENE=m_t2g, minGSSize = 10, maxGSSize = 200, pvalueCutoff = 0.05), mc.cores = 12)
en2 = list()
for (i in 1:length(en)) {
	en2[[i]] =  en[[i]]@result$ID[1:5]
}
en2 = unique(unlist(en2))
en3 = list()
for (i in 1:length(en)) {
	en3[[i]] =  cbind(i, en[[i]]@result$ID[which(en[[i]]@result$ID %in% en2)], -log10(as.numeric(en[[i]]@result$qvalue[which(en[[i]]@result$ID %in% en2)])))
}
en3 = data.frame(do.call(rbind, en3))
en3 = data.frame(en3$i, en3$V2, en3$V3)
colnames(en3) = c("Cluster", "Term", "PValue")
en3 = reshape(en3, idvar = "Cluster", timevar = "Term", direction = "wide")
en4 = data.frame(en3[2:length(en3)])
en4 = as.matrix(en4)
en4[which(is.na(en4))] = "0"
en4 = apply(en4, 2, as.numeric)
en4[which(en4 > -log10(0.0001))] = -log10(0.0001)
rownames(en4) = en3[[1]]
pdf(paste0(prefix, "KEGG_Lineage1.pdf"), width = 10, height = 10)
heatmap.2(t(en4), trace = "none", scale = "none", col = colorRampPalette(c("white","dodgerblue4"))(n=100), margins = c(20,25), density = "none", Colv = NA, key = FALSE, labCol = names(n3_l))
dev.off()





#Lineage 2 (cell cycle)
negSel = which(!is.na(st[,2]))
set.seed(111)
kk = Mclust(-st[negSel,2], 5, modelNames="E")
roses = logcounts(l_v[,negSel])
roses = roses[,order(st[negSel,2], decreasing = TRUE)]
names1 = toupper(unlist(lapply(strsplit(as.character(rownames(roses)), split = ";", fixed = T), function(x) as.character(x[2]))))
rownames(roses) = names1
reads_single = as.matrix(roses)
ans =
lapply(cycle_list,function(xx){
  reads_single_phase <- reads_single[rownames(reads_single) %in% toupper(unlist(xx)) ,]
  combined_matrix <- rbind(reads_single_phase,average=apply(reads_single_phase,2,mean))
  cor_matrix <- cor(t(combined_matrix))
  cor_vector <- cor_matrix[,dim(cor_matrix)[1]]
  reads_single_phase_restricted <- reads_single_phase[rownames(reads_single_phase) %in% names(cor_vector[cor_vector >= 0.1]),]
  apply(reads_single_phase_restricted,2,mean)
})
ans = matrix(unlist(ans), ncol = 3, byrow = FALSE)
colnames(ans) = c("G1_S", "G2_M", "G0")
ans_normed <- flexible_normalization(ans,by_row=FALSE)
ans_normed_normed <- flexible_normalization(ans_normed,by_row=TRUE)
cell_phase = apply(ans_normed_normed,1,function(x) colnames(ans_normed_normed)[which.max(x)])
n1 = seq(1,length(cell_phase), by = 2000)
n2 = n1 + 2000
n2[length(n2)] = length(cell_phase)
n3 = n2
for (i in 1:length(n2)) {
	n3[(n1[i]):(n2[i])] = i
}
ballet = lapply(1:length(n1), function(x) table(cell_phase[(n1[x]):(n2[x])]) / length((n1[x]):(n2[x])))
ballet = (matrix(unlist(ballet), byrow = TRUE, ncol = 3)) * 100
colnames(ballet) = c("G0", "G1_S", "G2_M")
pdf(paste0(prefix, "L2_cellcyle.pdf"))
plot(n1, ballet[,1], type = "l", col = makeTransparent(col_palette_short[1], 0.9), lwd = 5, xlab = "Cells ordered by pseudotime (Pericytes to MF1)", ylab = "% of Cells", ylim = c(20,50))
lines(n1, ballet[,2], type = "l", col = makeTransparent(col_palette_short[2], 0.9), lwd = 5)
lines(n1, ballet[,3], type = "l", col = makeTransparent(col_palette_short[4], 0.9), lwd = 5)
legend("topright", legend = colnames(ballet), pch=15, col = col_palette_short[c(1,2,4)], bty = "n", cex = 1.5)
dev.off()
#gene ontology
sg_n3  = sortGenes(roses, n3, binarizeMethod = "naive")
n3_l = plotTopMarkerHeat(sg_n3, averageCells=10^2, top_n=100, plotheat=FALSE, outs = TRUE)
m_t2g =  as.tibble(read.gmt("/home/mibrahim/Dropbox/dev/kidneyMap/public/assets/public/c2.cp.pid.v7.0.symbols.gmt"))
en = mclapply(n3_l, function(x) enricher(x, TERM2GENE=m_t2g, minGSSize = 10, maxGSSize = 200, pvalueCutoff = 0.05), mc.cores = 12)
en2 = list()
for (i in 1:length(en)) {
	en2[[i]] =  en[[i]]@result$ID[1:5]
}
en2 = unique(unlist(en2))
en3 = list()
for (i in 1:length(en)) {
	en3[[i]] =  cbind(i, en[[i]]@result$ID[which(en[[i]]@result$ID %in% en2)], -log10(as.numeric(en[[i]]@result$qvalue[which(en[[i]]@result$ID %in% en2)])))
}
en3 = data.frame(do.call(rbind, en3))
en3 = data.frame(en3$i, en3$V2, en3$V3)
colnames(en3) = c("Cluster", "Term", "PValue")
en3 = reshape(en3, idvar = "Cluster", timevar = "Term", direction = "wide")
en4 = data.frame(en3[2:length(en3)])
en4 = as.matrix(en4)
en4[which(is.na(en4))] = "0"
en4 = apply(en4, 2, as.numeric)
en4[which(en4 > -log10(0.0001))] = -log10(0.0001)
rownames(en4) = en3[[1]]
pdf(paste0(prefix, "PID_Lineage2.pdf"), width = 10, height = 10)
heatmap.2(t(en4), trace = "none", scale = "none", col = colorRampPalette(c("white","dodgerblue4"))(n=100), margins = c(20,25), density = "none", Colv = NA, key = FALSE, labCol = names(n3_l))
dev.off()
sg_n3  = sortGenes(roses, n3, binarizeMethod = "naive")
n3_l = plotTopMarkerHeat(sg_n3, averageCells=10^2, top_n=100, plotheat=FALSE, outs = TRUE)
m_t2g =  as.tibble(read.gmt("/home/mibrahim/Dropbox/dev/kidneyMap/public/assets/public/c2.cp.kegg.v7.0.symbols.gmt"))
en = mclapply(n3_l, function(x) enricher(x, TERM2GENE=m_t2g, minGSSize = 10, maxGSSize = 200, pvalueCutoff = 0.05), mc.cores = 12)
en2 = list()
for (i in 1:length(en)) {
	en2[[i]] =  en[[i]]@result$ID[1:5]
}
en2 = unique(unlist(en2))
en3 = list()
for (i in 1:length(en)) {
	en3[[i]] =  cbind(i, en[[i]]@result$ID[which(en[[i]]@result$ID %in% en2)], -log10(as.numeric(en[[i]]@result$qvalue[which(en[[i]]@result$ID %in% en2)])))
}
en3 = data.frame(do.call(rbind, en3))
en3 = data.frame(en3$i, en3$V2, en3$V3)
colnames(en3) = c("Cluster", "Term", "PValue")
en3 = reshape(en3, idvar = "Cluster", timevar = "Term", direction = "wide")
en4 = data.frame(en3[2:length(en3)])
en4 = as.matrix(en4)
en4[which(is.na(en4))] = "0"
en4 = apply(en4, 2, as.numeric)
en4[which(en4 > -log10(0.0001))] = -log10(0.0001)
rownames(en4) = en3[[1]]
pdf(paste0(prefix, "KEGG_Lineage2.pdf"), width = 10, height = 10)
heatmap.2(t(en4), trace = "none", scale = "none", col = colorRampPalette(c("white","dodgerblue4"))(n=100), margins = c(20,25), density = "none", Colv = NA, key = FALSE, labCol = names(n3_l))
dev.off()






#Lineage 3 (cell cycle)
negSel = which(!is.na(st[,3]))
set.seed(111)
kk = Mclust(-st[negSel,3], 5, modelNames="E")
roses = logcounts(l_v[,negSel])
roses = roses[,order(st[negSel,3], decreasing = TRUE)]
names1 = toupper(unlist(lapply(strsplit(as.character(rownames(roses)), split = ";", fixed = T), function(x) as.character(x[2]))))
rownames(roses) = names1
reads_single = as.matrix(roses)
ans =
lapply(cycle_list,function(xx){
  reads_single_phase <- reads_single[rownames(reads_single) %in% toupper(unlist(xx)) ,]
  combined_matrix <- rbind(reads_single_phase,average=apply(reads_single_phase,2,mean))
  cor_matrix <- cor(t(combined_matrix))
  cor_vector <- cor_matrix[,dim(cor_matrix)[1]]
  reads_single_phase_restricted <- reads_single_phase[rownames(reads_single_phase) %in% names(cor_vector[cor_vector >= 0.1]),]
  apply(reads_single_phase_restricted,2,mean)
})
ans = matrix(unlist(ans), ncol = 3, byrow = FALSE)
colnames(ans) = c("G1_S", "G2_M", "G0")
ans_normed <- flexible_normalization(ans,by_row=FALSE)
ans_normed_normed <- flexible_normalization(ans_normed,by_row=TRUE)
cell_phase = apply(ans_normed_normed,1,function(x) colnames(ans_normed_normed)[which.max(x)])
n1 = seq(1,length(cell_phase), by = 2000)
n2 = n1 + 2000
n2[length(n2)] = length(cell_phase)
n3 = n2
for (i in 1:length(n2)) {
	n3[(n1[i]):(n2[i])] = i
}
ballet = lapply(1:length(n1), function(x) table(cell_phase[(n1[x]):(n2[x])]) / length((n1[x]):(n2[x])))
ballet = (matrix(unlist(ballet), byrow = TRUE, ncol = 3)) * 100
colnames(ballet) = c("G0", "G1_S", "G2_M")
pdf(paste0(prefix, "L3_cellcyle.pdf"))
plot(n1, ballet[,1], type = "l", col = makeTransparent(col_palette_short[1], 0.9), lwd = 5, xlab = "Cells ordered by pseudotime (Fib2 to MF1)", ylab = "% of Cells", ylim = c(25,45))
lines(n1, ballet[,2], type = "l", col = makeTransparent(col_palette_short[2], 0.9), lwd = 5)
lines(n1, ballet[,3], type = "l", col = makeTransparent(col_palette_short[4], 0.9), lwd = 5)
legend("bottomright", legend = colnames(ballet), pch=15, col = col_palette_short[c(1,2,4)], bty = "n", cex = 1.5)
dev.off()
#gene ontology
sg_n3  = sortGenes(roses, n3, binarizeMethod = "naive")
n3_l = plotTopMarkerHeat(sg_n3, averageCells=10^2, top_n=100, plotheat=FALSE, outs = TRUE)
m_t2g =  as.tibble(read.gmt("/home/mibrahim/Dropbox/dev/kidneyMap/public/assets/public/c2.cp.pid.v7.0.symbols.gmt"))
en = mclapply(n3_l, function(x) enricher(x, TERM2GENE=m_t2g, minGSSize = 10, maxGSSize = 200, pvalueCutoff = 0.05), mc.cores = 12)
en2 = list()
for (i in 1:length(en)) {
	en2[[i]] =  en[[i]]@result$ID[1:5]
}
en2 = unique(unlist(en2))
en3 = list()
for (i in 1:length(en)) {
	en3[[i]] =  cbind(i, en[[i]]@result$ID[which(en[[i]]@result$ID %in% en2)], -log10(as.numeric(en[[i]]@result$qvalue[which(en[[i]]@result$ID %in% en2)])))
}
en3 = data.frame(do.call(rbind, en3))
en3 = data.frame(en3$i, en3$V2, en3$V3)
colnames(en3) = c("Cluster", "Term", "PValue")
en3 = reshape(en3, idvar = "Cluster", timevar = "Term", direction = "wide")
en4 = data.frame(en3[2:length(en3)])
en4 = as.matrix(en4)
en4[which(is.na(en4))] = "0"
en4 = apply(en4, 2, as.numeric)
en4[which(en4 > -log10(0.0001))] = -log10(0.0001)
rownames(en4) = en3[[1]]
pdf(paste0(prefix, "PID_Lineage3.pdf"), width = 10, height = 10)
heatmap.2(t(en4), trace = "none", scale = "none", col = colorRampPalette(c("white","dodgerblue4"))(n=100), margins = c(20,25), density = "none", Colv = NA, key = FALSE, labCol = names(n3_l))
dev.off()
sg_n3  = sortGenes(roses, n3, binarizeMethod = "naive")
n3_l = plotTopMarkerHeat(sg_n3, averageCells=10^2, top_n=100, plotheat=FALSE, outs = TRUE)
m_t2g =  as.tibble(read.gmt("/home/mibrahim/Dropbox/dev/kidneyMap/public/assets/public/c2.cp.kegg.v7.0.symbols.gmt"))
en = mclapply(n3_l, function(x) enricher(x, TERM2GENE=m_t2g, minGSSize = 10, maxGSSize = 200, pvalueCutoff = 0.05), mc.cores = 12)
en2 = list()
for (i in 1:length(en)) {
	en2[[i]] =  en[[i]]@result$ID[1:5]
}
en2 = unique(unlist(en2))
en3 = list()
for (i in 1:length(en)) {
	en3[[i]] =  cbind(i, en[[i]]@result$ID[which(en[[i]]@result$ID %in% en2)], -log10(as.numeric(en[[i]]@result$qvalue[which(en[[i]]@result$ID %in% en2)])))
}
en3 = data.frame(do.call(rbind, en3))
en3 = data.frame(en3$i, en3$V2, en3$V3)
colnames(en3) = c("Cluster", "Term", "PValue")
en3 = reshape(en3, idvar = "Cluster", timevar = "Term", direction = "wide")
en4 = data.frame(en3[2:length(en3)])
en4 = as.matrix(en4)
en4[which(is.na(en4))] = "0"
en4 = apply(en4, 2, as.numeric)
en4[which(en4 > -log10(0.0001))] = -log10(0.0001)
rownames(en4) = en3[[1]]
pdf(paste0(prefix, "KEGG_Lineage3.pdf"), width = 10, height = 10)
heatmap.2(t(en4), trace = "none", scale = "none", col = colorRampPalette(c("white","dodgerblue4"))(n=100), margins = c(20,25), density = "none", Colv = NA, key = FALSE, labCol = names(n3_l))
dev.off()



#AP1 analysis
matrisome_set = read.table("https://raw.githubusercontent.com/mahmoudibrahim/KidneyMap/master/assets/public/ecm_genes_human.txt", sep = "\t", header = TRUE)
collagens = matrisome_set$Gene.Symbol[which(matrisome_set$Category %in% c("Collagens"))]
reads_single_phase = logcounts(l_remove)
rownames(reads_single_phase) = toupper(unlist(lapply(strsplit(as.character(rownames(reads_single_phase)), split = ";", fixed = T), function(x) as.character(x[2]))))
reads_single_phase = as.matrix(reads_single_phase[rownames(reads_single_phase) %in% toupper(collagens) ,])
combined_matrix = rbind(reads_single_phase,average=apply(reads_single_phase,2,mean))
cor_matrix = cor(t(combined_matrix))
cor_vector = cor_matrix[,dim(cor_matrix)[1]]
reads_single_phase_restricted = reads_single_phase[rownames(reads_single_phase) %in% names(cor_vector[cor_vector >= 0.1]),]
aaa = apply(reads_single_phase_restricted,2,mean)
ecm = aaa[ww]


library(ggpubr)
ap1 = grep("^JUN|^FOS", dorothea_hs$tf, value=T)
collagens = ap1
reads_single_phase = logcounts(l_remove)
rownames(reads_single_phase) = toupper(unlist(lapply(strsplit(as.character(rownames(reads_single_phase)), split = ";", fixed = T), function(x) as.character(x[2]))))
reads_single_phase = as.matrix(reads_single_phase[rownames(reads_single_phase) %in% toupper(collagens) ,])
combined_matrix = rbind(reads_single_phase,average=apply(reads_single_phase,2,mean))
cor_matrix = cor(t(combined_matrix))
cor_vector = cor_matrix[,dim(cor_matrix)[1]]
reads_single_phase_restricted = reads_single_phase[rownames(reads_single_phase) %in% names(cor_vector[cor_vector >= 0.1]),]
aaa = apply(reads_single_phase_restricted,2,mean)
aaa = aaa[ww]
pdf(paste0(prefix, "ap1_score_collagen_score.pdf"), width = 14)
www = which(l_v$level2 == "Pericytes")
df = data.frame(y=aaa[www],x=ecm[www]) #ecm is 
o1 <- ggplot(df, aes(x, y)) +
  geom_point(alpha = 0.15, col = colors_info_mesenl3[www]) +
  geom_density_2d(col = "#00000080") + 
   lims(x= c(0.1,2.2), y = c(0.1,5)) + 
	theme_bw()
www = which(l_v$level2 == "Fibroblasts")
df = data.frame(y=aaa[www],x=ecm[www])
o2 <- ggplot(df, aes(x, y)) +
  geom_point(alpha = 0.15, col = colors_info_mesenl3[www]) +
  geom_density_2d(col = "#00000080") + 
   lims(x= c(0.1,2.2), y = c(0.1,5)) + 
	theme_bw()
www = which(l_v$level2 == "Myofibroblasts")
df = data.frame(y=aaa[www],x=ecm[www])
o3 <- ggplot(df, aes(x, y)) +
  geom_point(alpha = 0.15, col = colors_info_mesenl3[www]) +
  geom_density_2d(col = "#00000080") + 
   lims(x= c(0.1,2.2), y = c(0.1,5)) + 
	theme_bw()
ggarrange(o1,o2,o3,ncol = 3, nrow = 1, labels = c("Peri","Fib","MF"))
dev.off()



ap1 = grep("^JUN|^FOS", dorothea_hs$tf, value=T)
ap1 = which(dorothea_hs$tf %in% ap1)
ap1 = dorothea_hs[ap1,]
ap1 = ap1[which(ap1$confidence == "A"),]
ap1 = ap1[which(ap1$mor == 1),]
ap1 = sort(unique(ap1$target))
collagens = ap1
reads_single_phase = logcounts(l_remove)
rownames(reads_single_phase) = toupper(unlist(lapply(strsplit(as.character(rownames(reads_single_phase)), split = ";", fixed = T), function(x) as.character(x[2]))))
reads_single_phase = as.matrix(reads_single_phase[rownames(reads_single_phase) %in% toupper(collagens) ,])
combined_matrix = rbind(reads_single_phase,average=apply(reads_single_phase,2,mean))
cor_matrix = cor(t(combined_matrix))
cor_vector = cor_matrix[,dim(cor_matrix)[1]]
reads_single_phase_restricted = reads_single_phase[rownames(reads_single_phase) %in% names(cor_vector[cor_vector >= 0.1]),]
aaa = apply(reads_single_phase_restricted,2,mean)
aaa = aaa[ww]
pdf(paste0(prefix, "ap1_score_collagen_score_doro.pdf"), width = 14)
www = which(l_v$level2 == "Pericytes")
df = data.frame(y=aaa[www],x=ecm[www])
o1 <- ggplot(df, aes(x, y)) +
  geom_point(alpha = 0.15, col = colors_info_mesenl3[www]) +
  geom_density_2d(col = "#00000080") + 
   lims(x= c(0.1,2.2), y = c(0.1,1.2)) + 
	theme_bw()
www = which(l_v$level2 == "Fibroblasts")
df = data.frame(y=aaa[www],x=ecm[www])
o2 <- ggplot(df, aes(x, y)) +
  geom_point(alpha = 0.15, col = colors_info_mesenl3[www]) +
  geom_density_2d(col = "#00000080") + 
   lims(x= c(0.1,2.2), y = c(0.1,1.2)) + 
	theme_bw()
www = which(l_v$level2 == "Myofibroblasts")
df = data.frame(y=aaa[www],x=ecm[www])
o3 <- ggplot(df, aes(x, y)) +
  geom_point(alpha = 0.15, col = colors_info_mesenl3[www]) +
  geom_density_2d(col = "#00000080") + 
   lims(x= c(0.1,2.2), y = c(0.1,1.2)) + 
	theme_bw()
ggarrange(o1,o2,o3,ncol = 3, nrow = 1, labels = c("Peri","Fib","MF"))
dev.off()




#Nkd2 analysis
ll = logcounts(l_v)
ll = apply(ll, 1, function(x) x/(sqrt(sum(x^2))))
ll = t(ll)
rownames(ll) = toupper(unlist(lapply(strsplit(as.character(rownames(ll)), split = ";", fixed = T), function(x) as.character(x[2]))))
ll = ll[!duplicated(rownames(ll)),]

nk = as.vector(ll[which(rownames(ll) == "NKD2"),])
corr = apply(ll, 1, function(x) cor(x, nk))

top = head(names(sort(corr, decreasing = TRUE)), n = 100)
bottom = head(names(sort(corr, decreasing = FALSE)), n = 100)
write(c(top, bottom), file = paste0(prefix, "_NKD2_geneList.csv"), ncolumns = 1)


ng = list("correlated" = top, "anti-correlated" = bottom)
msig = msigdbr::msigdbr(species = "Homo sapiens", category = "C2")
msig = data.frame(msig$gs_name[which(msig$gs_subcat == "CP:PID")], msig$gene_symbol[which(msig$gs_subcat == "CP:PID")])
en = mclapply(ng, function(x) enricher(x, TERM2GENE=msig, minGSSize = 10, maxGSSize = 200, pvalueCutoff = 0.05), mc.cores = 2)
en2 = list()
for (i in 1:length(en)) {
	en2[[i]] =  en[[i]]@result$ID[1:5]
}
en2 = unique(unlist(en2))
en3 = list()
for (i in 1:length(en)) {
	en3[[i]] =  cbind(i, en[[i]]@result$ID[which(en[[i]]@result$ID %in% en2)], -log10(as.numeric(en[[i]]@result$qvalue[which(en[[i]]@result$ID %in% en2)])))
}
en3 = data.frame(do.call(rbind, en3))
en3 = data.frame(en3$i, en3$V2, en3$V3)
colnames(en3) = c("Cluster", "Term", "PValue")
en3 = reshape(en3, idvar = "Cluster", timevar = "Term", direction = "wide")
en4 = data.frame(en3[2:length(en3)])
en4 = as.matrix(en4)
en4[which(is.na(en4))] = "0"
en4 = apply(en4, 2, as.numeric)
en4[which(en4 > -log10(0.0001))] = -log10(0.0001)
rownames(en4) = en3[[1]]
pdf(paste0(prefix, "PID_NKD2.pdf"), width = 10, height = 10)
heatmap.2(t(en4), trace = "none", scale = "none", col = colorRampPalette(c("white","dodgerblue4"))(n=100), margins = c(20,25), density = "none", key = FALSE, labCol = names(ng), Colv = NA)
dev.off()



top = head(names(sort(corr, decreasing = TRUE)), n = 100)
bottom = head(names(sort(corr, decreasing = FALSE)), n = 100)
ng = list("correlated" = top, "anti-correlated" = bottom)
msig = msigdbr::msigdbr(species = "Homo sapiens", category = "C2")
msig = data.frame(msig$gs_name[which(msig$gs_subcat == "CP:KEGG")], msig$gene_symbol[which(msig$gs_subcat == "CP:KEGG")])
cp = clusterProfiler::compareCluster(geneCluster = ng, fun = "enricher", TERM2GENE = msig)
en = mclapply(ng, function(x) enricher(x, TERM2GENE=msig, minGSSize = 10, maxGSSize = 200, pvalueCutoff = 0.05), mc.cores = 2)
en2 = list()
for (i in 1:length(en)) {
	en2[[i]] =  en[[i]]@result$ID[1:5]
}
en2 = unique(unlist(en2))
en3 = list()
for (i in 1:length(en)) {
	en3[[i]] =  cbind(i, en[[i]]@result$ID[which(en[[i]]@result$ID %in% en2)], -log10(as.numeric(en[[i]]@result$qvalue[which(en[[i]]@result$ID %in% en2)])))
}
en3 = data.frame(do.call(rbind, en3))
en3 = data.frame(en3$i, en3$V2, en3$V3)
colnames(en3) = c("Cluster", "Term", "PValue")
en3 = reshape(en3, idvar = "Cluster", timevar = "Term", direction = "wide")
en4 = data.frame(en3[2:length(en3)])
en4 = as.matrix(en4)
en4[which(is.na(en4))] = "0"
en4 = apply(en4, 2, as.numeric)
en4[which(en4 > -log10(0.0001))] = -log10(0.0001)
rownames(en4) = en3[[1]]
pdf(paste0(prefix, "KEGG_NKD2.pdf"), width = 10, height = 10)
heatmap.2(t(en4), trace = "none", scale = "none", col = colorRampPalette(c("white","dodgerblue4"))(n=100), margins = c(20,25), density = "none", key = FALSE, labCol = names(ng), cexRow = 0.8, cexCol = 0.8, Colv = NA)
dev.off()


pdf(paste0(prefix, "NKD2_correlated_genes.pdf"), height = 15)
plotMarkerHeat(ll, l_v$level2, unlist(ng), averageCells=10^3, clusterGenes = TRUE, newOrder = c(3,1,2), colors = colorRampPalette(c("cornflowerblue","black","gold"))(n=100))
dev.off()
pdf(paste0(prefix, "NKD2_select_genes.pdf"))
plotMarkerHeat(ll, l_v$level2, c(head(ng[[1]], 11), head(ng[[2]], 10)), averageCells=10^3, clusterGenes = TRUE, newOrder = c(3,1,2), colors = colorRampPalette(c("cornflowerblue","black","gold"))(n=100))
dev.off()


#Nkd2 GRN
library(igraph)
library(ggraph)
library(tibble)
library(clusterProfiler)
col_palette = c("#d64e3c","#8641c6","#7778cb","#59803d","#d04d9c","#73d6a8","#492f60","#ccc497","#7f343b","#72acc0","#b97d40","#c796b5","#45483a","purple","green","yellow")
l= read.table("HumanPDGFRBpositive_Nkd2_grnboost2.csv")
l = l[which(l$V3 > 10),]
gg = graph_from_data_frame(l)
E(gg)$weight = log10(l$V3)
gg2 = as.undirected(gg)
cl = cluster_louvain(gg2)
pdf("NKD2_geneList_net.pdf", width = 10, height = 10)
set.seed(111)
lay = ggraph::create_layout(gg, layout = "fr")
ggraph(lay) + 
  geom_edge_link(aes(alpha = E(gg)$weight / 10), edge_colour = "gray") + 
  geom_node_point(aes(), colour = "#99999910") + 
  geom_node_text(aes(label = V(gg)$name), repel=TRUE, cex = 5, color = col_palette[cl$membership]) +
  theme(panel.background = element_blank())
dev.off()
df = data.frame("gene" = V(gg)$name, "cluster" = cl$membership)
write.csv(df, "HumanPDGFRBpositive_Nkd2_grnboost2_netClusters.csv")
df = list("1" = as.character(df$gene[which(df$cluster == 1)]), "2" = as.character(df$gene[which(df$cluster == 2)]), "3" = as.character(df$gene[which(df$cluster == 3)]), "4" = as.character(df$gene[which(df$cluster == 4)]))
gg_myofib = induced_subgraph(gg, which(cl$membership == 2))
nkd2_edges = apply(as.data.frame(get.edgelist(gg_myofib)), 1, function(x) any(x%in%c("NKD2")))
pdf("NKD2_geneList_net_Myofibroblast.pdf", width = 12, height = 10)
set.seed(111)
lay = ggraph::create_layout(gg_myofib, layout = "fr")
ggraph(lay) + 
  geom_edge_link(aes(col = as.factor(nkd2_edges))) + 
  scale_edge_colour_manual(values = c("#99999920","#d64e3c")) + 
  geom_node_point(aes(), colour = "#99999910") + 
  geom_node_text(aes(label = V(gg_myofib)$name), repel=TRUE, cex = 5, color = "black") +
  theme(panel.background = element_blank())
dev.off()





#save output
writeMM(counts(l_remove), file = paste0(prefix, "UMI_counts.mtx"))
names1 = toupper(unlist(lapply(strsplit(as.character(rownames(l_remove)), split = ";", fixed = T), function(x) as.character(x[2]))))
names2 = toupper(unlist(lapply(strsplit(as.character(rownames(l_remove)), split = ";", fixed = T), function(x) as.character(x[1]))))
names2 = toupper(unlist(lapply(strsplit(as.character(names2), split = ".", fixed = T), function(x) as.character(x[1]))))
rowdat = cbind(names1, names2)
colnames(rowdat) = c("Gene.Symbol","ENSEMBL.ID")
write.table(rowdat, file = paste0(prefix, "UMI_counts_rowData.txt"), sep = "\t", row.names=FALSE, quote = FALSE)
coldat = cbind(l_remove$level1, l_remove$level2, l_remove$level3, l_remove$kidney_function, l_remove$fileID)
colnames(coldat) = c("Annotation.Level.1","Annotation.Level.2","Annotation.Level.3","Kidney.Function","Patient ID")
write.table(coldat, file = paste0(prefix, "UMI_counts_colData.txt"), sep = "\t", row.names=FALSE, quote = FALSE)



##Fin
