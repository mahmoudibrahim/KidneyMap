########################################################################
# human_CD10negative.r
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
level2_col_one = (colorRampPalette(c("#4363d8","white"))(n=7))[-7]
level2_col_one_trans = makeTransparent(level2_col_one, alpha = 0.9)
level2_col_two = (colorRampPalette(c("#42006a","white"))(n=10))[-10]
level2_col_two_trans = makeTransparent(level2_col_two, alpha = 0.9)
level2_col_three = (colorRampPalette(c("#800000","white"))(n=9))[-9]
level2_col_three_trans = makeTransparent(level2_col_three, alpha = 0.9)
level2_col_four = (colorRampPalette(c("#9A6324","white"))(n=6))[-6]
level2_col_four_trans = makeTransparent(level2_col_four, alpha = 0.9)
level2_col_five = c("#808000", "#d8d35f","#e7e16d")
level2_col_five_trans = makeTransparent(level2_col_five, alpha = 0.9)
kf_col = c("#999999","#fc8d62")
kf_col_trans = makeTransparent(kf_col, alpha = 0.4)

prefix = "human_CD10negative_"



###########################################################################
########################START ANALYSIS#####################################
###########################################################################

#load data
load("human_CD10negative.RData") #R object containing aggregated data, saved as SingleCellExperiment. Variable name: l
lx = l #save full matrix for later


#gene filtering
cutGenes = log10(floor(ncol(l) * 0.001))
choose_genes = which((log10(Matrix::rowSums(counts(l)))) >= cutGenes)
l = l[choose_genes, ]


#merge highly similar clusters
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
for (i in 1:length(mergeMap)) {
	class_info_merging[which(class_info_merging %in% mergeMap[[i]])] = paste0("merge", i)
}



#remove clusters with no differential gene expression
sg = sortGenes(counts(l), class_info_merging, binarizeMethod="naive", cores = 16)
removethis = c(names(which(apply(pp$adjpval,2,function(x) length(x[x<0.05])) == 0)), "en7") #en7 is endothelial cell doublets
removethis = which(class_info_merging %in% removethis)




#reform matrix and refilter genes
class_info_merging_remove = class_info_merging[-removethis]
l_remove = lx
l_remove = l_remove[,-removethis]
l_remove = calculateQCMetrics(l_remove)
colData(l_remove)$final_classes = class_info_merging_remove
cutGenes = log10(floor(ncol(l_remove) * 0.001)) #on average at least 1UMI in 0.1% of cells 
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
annot = read.table("https://raw.githubusercontent.com/mahmoudibrahim/KidneyMap/master/assets/clusterInfo/human_CD10negative.txt", sep = "\t")
annot = annot[match(names(KFlow), annot$V1),]
annot$V7 = KFlow
annot$V8 = log(KFlow) - log((table(l_remove$kidney_function) / ncol(l_remove))[2])
nums = table(class_info_merging_remove)
annot = annot[match(names(nums), annot$V1),]
annot$V9 = as.numeric(nums)
annot$V10 = as.numeric(log10(nums))
colnames(annot) = c("ID", "Label", "Annotation 1", "Annotation 2", "Annotation 3", "Percent Low KF Cells", "Log Enrichment of Low KF Cells", "Number of Cells", "log10 Number of Cells")
#write cluster annotation information
write.csv(annot, paste0(prefix, "cluster_annotation_info_withAddData.csv"), row.names = FALSE)
#annotation levels
colData(l_remove)$level1 = l_remove$class_info
colData(l_remove)$level2 = l_remove$class_info
colData(l_remove)$level3 = l_remove$class_info
colData(l_remove)$label = l_remove$class_info
for (i in 1:length(unique(l_remove$final_classes))) {
	where = which(l_remove$final_classes == unique(l_remove$final_classes)[i])
	l_remove$level3[where] = as.character(annot[[5]][which(annot[[1]] == unique(l_remove$final_classes)[i])])
	l_remove$level2[where] = as.character(annot[[4]][which(annot[[1]] == unique(l_remove$final_classes)[i])])
	l_remove$level1[where] = as.character(annot[[3]][which(annot[[1]] == unique(l_remove$final_classes)[i])])
	l_remove$label[where] = as.character(annot[[2]][which(annot[[1]] == unique(l_remove$final_classes)[i])])
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
	colors_info_kf[wtemp] = (rev(kf_col_trans)[i])
}
colors_file = rep("", length = ncol(l_remove))
for (i in 1:length(unique(l_remove$fileID))) {
	wtemp = which((l_remove$fileID) == (sort(unique(l_remove$fileID)))[i])
	colors_file[wtemp] = col_palette_trans[i]
}
colors_patient = rep("", length = ncol(l_remove))
for (i in 1:length(unique(l_remove$patientID))) {
	wtemp = which((l_remove$patientID) == (sort(unique(l_remove$patientID)))[i])
	colors_patient[wtemp] = col_palette_trans[i]
}



#dimension reduction
sg = sortGenes(counts(l_remove), class_info_merging_remove, binarizeMethod = "naive", cores = 16)
pp = getPValues(sg, numPerm = 20, cores = 16)
l_remove = scran::computeSumFactors(l_remove, sizes = seq(10, 200, 20), clusters = class_info_merging_remove, positive = TRUE)
l_remove = scater::normalize(l_remove)
l_remove = scater::calculateQCMetrics(l_remove)
l_v = l_remove[which(rownames(l_remove) %in% names(which(apply(pp$adjpval, 1, function(x) any(x < 0.05))))),]
q = t(logcounts(l_v))
q = t(apply(q,1, function(x) x-mean(x)))
set.seed(111)
sv = irlba(q, 50)
d = sv$d
dw = 21 #based on the knee of the singular values
dw = 1:dw
temp = t(apply(ss$u[,dw], 1, function(x) x / (sqrt(sum(x^2)))))
write.csv(temp, file = paste0(prefix, "dims.txt"), row.names = FALSE, col.names = FALSE) #file saved @https://raw.githubusercontent.com/mahmoudibrahim/KidneyMap/master/assets/reducedDims/SVD/human_CD10negative_dims.txt



#UMAP (produced in python, just for visualization, independent of clustering)
mapCoords = read.table("https://raw.githubusercontent.com/mahmoudibrahim/KidneyMap/master/assets/reducedDims/UMAP/human_CD10negative_umapCoords.csv", header = FALSE, sep=",")
pdf(paste0(prefix, "regularUMAP_level1_level3TEXT.pdf"))
plot(mapCoords[,1], mapCoords[,2], pch = 20, col = colors_info_level1, cex = 0.2, xlim = c(-12,16))
allClasses = unique(l_remove$label)
for (i in 1:length(allClasses)) {
	ww = which(l_remove$label == allClasses[i])
	xx = median(mapCoords[ww,1])
	yy = median(mapCoords[ww,2])
	text(xx,yy,allClasses[i], cex = 0.7, adj = 0.5)
}
dev.off()
pdf(paste0(prefix, "regularUMAP_KF.pdf"))
plot(mapCoords[,1], mapCoords[,2], pch = 20, col = colors_info_kf, cex = 0.2, xlim = c(-12,16))
dev.off()
pdf(paste0(prefix, "clusterPercent_versus_KF.pdf"), height = 15)
howmany = table(l_remove$level3, l_remove$kidney_function)
howmany = t(apply(howmany, 1, function(x) x/sum(x)))
heatmap.2(howmany, trace = "none", col = colorRampPalette(c("white","#0072b2"))(n=100), cellnote=round(howmany*100,2), notecol="black", margins = c(15,15), key = FALSE, dendrogram = "none")
dev.off()

#corrnet visualization
corr = plotCorrelationHeat(sg, displayNumbers=F, markers=getMarkers(sg, quant = 0.99)$markers, outs = TRUE)
classNet = corr$corr
diag(classNet) = 0
classNet[which(corr$corr < 0)] = 0
classNet = graph_from_adjacency_matrix(classNet, weighted = TRUE)
classNet = simplify(classNet, edge.attr.comb = "max")
classNetel = as_edgelist(classNet)
classNetew = E(classNet)$weight
classNet = cbind(classNetel, classNetew)
colnames(classNet) = c("source","target","Weight")
write.csv(classNet, paste0(prefix, "_cluster_annotation_addInfo_forGephi.txt"), row.names = FALSE)
############################################################



###big marker heatmaps [last run Oct18, should be no change]
level3_order = c(25, 9:14, 28:30, 33:34, 53:54, 20:22, 37:40, 7, 43:45, 5, 8, 4, 23, 36, 1, 15:16, 46:52,35 , 55, 24, 2, 3, 31,32, 42, 6, 26:27, 17:19, 41)
sg_temp = sortGenes(logcounts(l_remove), l_remove$level3, binarizeMethod = "naive", cores = 16)


tiff(paste0(prefix,"level3_top_auto_markers_10.tiff"), res=300, height = 25, width = 12, units = "in")
plotTopMarkerHeat(sg_temp, averageCells=10^5, newOrder=level3_order, gaps = FALSE, colors = colorRampPalette( c("cornflowerblue","black","gold"))(n=100), top_n=10)
dev.off()




###core matrisome score
l_remove$level3KF = paste0(l_remove$level3, "_", l_remove$kidney_function)
l_remove$level2KF = paste0(l_remove$level2, "_",  l_remove$kidney_function)
l_remove$level1KF = paste0(l_remove$level1, "_",  l_remove$kidney_function)
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
size = scale(sqrt(table(l_remove$level1KF)),center = FALSE)
vioplot(aaa[which(l_remove$kidney_function == "Healthy")] ~ l_remove$level1[which(l_remove$kidney_function == "Healthy")], col = makeTransparent(level1_col, alpha=0.4),side = "left", ylim = c(0,1.4), plotCentre="line", wex = size[c(2,4,6,8,10)], lwd = 3, border = level1_col)
vioplot(aaa[which(l_remove$kidney_function == "CKD")] ~ l_remove$level1[which(l_remove$kidney_function == "CKD")], col = makeTransparent(level1_col,alpha=0.7), side = "right", add = TRUE, plotCentre="line", wex = size[c(1,3,5,7,9)], lwd = 3, border = level1_col)
dev.off()

pdf(paste0(prefix, "_core_matrisome_score_level2_mesenchymal.pdf"))
size = scale(sqrt(table(l_remove$level2KF[which(l_remove$level1 == "Mesenchymal")])),center = FALSE)
vioplot(aaa[which((l_remove$level1 == "Mesenchymal") & (l_remove$kidney_function == "Healthy"))] ~ l_remove$level2KF[which((l_remove$level1 == "Mesenchymal") & (l_remove$kidney_function == "Healthy"))], col = makeTransparent(level1_col[4], alpha=0.4),side = "left", ylim = c(0,1.4), plotCentre="line", wex = size[c(2,4,6,8,10)], lwd = 3, border = level1_col[4])
vioplot(aaa[which((l_remove$level1 == "Mesenchymal") & (l_remove$kidney_function == "CKD"))] ~ l_remove$level2KF[which((l_remove$level1 == "Mesenchymal") & (l_remove$kidney_function == "CKD"))], col = makeTransparent(level1_col[4], alpha=0.8),side = "right", ylim = c(0,1.4), plotCentre="line", wex = size[c(1,3,5,7,9)], lwd = 3, border = level1_col[4],add = TRUE)
dev.off()

pdf(paste0(prefix, "_core_matrisome_score_level2_epithelial.pdf"))
size = scale(sqrt(table(l_remove$level2KF[which(l_remove$level1 == "Epithelial")])),center = FALSE)
vioplot(aaa[which((l_remove$level1 == "Epithelial") & (l_remove$kidney_function == "Healthy"))] ~ l_remove$level2KF[which((l_remove$level1 == "Epithelial") & (l_remove$kidney_function == "Healthy"))], col = makeTransparent(level1_col[2], alpha=0.4),side = "left", ylim = c(0,1.4), plotCentre="line", wex = size[c(2,4,6,8,10,12,14,16,18)], lwd = 3, border = level1_col[2])
vioplot(aaa[which((l_remove$level1 == "Epithelial") & (l_remove$kidney_function == "CKD"))] ~ l_remove$level2KF[which((l_remove$level1 == "Epithelial") & (l_remove$kidney_function == "CKD"))], col = makeTransparent(level1_col[2], alpha=0.8),side = "right", ylim = c(0,1.4), plotCentre="line", wex = size[c(1,3,5,7,9,11,13,15,17)], lwd = 3, border = level1_col[2],add = TRUE)
dev.off()



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


pdf(paste0(prefix, "collagen_score_level1.pdf"))
size = scale(sqrt(table(l_remove$level1KF)),center = FALSE)
vioplot(aaa[which(l_remove$kidney_function == "Healthy")] ~ l_remove$level1[which(l_remove$kidney_function == "Healthy")], col = makeTransparent(level1_col, alpha=0.4),side = "left", ylim = c(0,1.8), plotCentre="line", wex = size[c(2,4,6,8,10)], lwd = 3, border = level1_col)
vioplot(aaa[which(l_remove$kidney_function == "CKD")] ~ l_remove$level1[which(l_remove$kidney_function == "CKD")], col = makeTransparent(level1_col,alpha=0.7), side = "right", add = TRUE, plotCentre="line", wex = size[c(1,3,5,7,9)], lwd = 3, border = level1_col)
dev.off()

pdf(paste0(prefix, "collagen_score_level2_mesenchymal.pdf"))
size = scale(sqrt(table(l_remove$level2KF[which(l_remove$level1 == "Mesenchymal")])),center = FALSE)
vioplot(aaa[which((l_remove$level1 == "Mesenchymal") & (l_remove$kidney_function == "Healthy"))] ~ l_remove$level2KF[which((l_remove$level1 == "Mesenchymal") & (l_remove$kidney_function == "Healthy"))], col = makeTransparent(level1_col[4], alpha=0.4),side = "left", ylim = c(0,1.8), plotCentre="line", wex = size[c(2,4,6,8,10)], lwd = 3, border = level1_col[4])
vioplot(aaa[which((l_remove$level1 == "Mesenchymal") & (l_remove$kidney_function == "CKD"))] ~ l_remove$level2KF[which((l_remove$level1 == "Mesenchymal") & (l_remove$kidney_function == "CKD"))], col = makeTransparent(level1_col[4], alpha=0.8),side = "right", ylim = c(0,1.8), plotCentre="line", wex = size[c(1,3,5,7,9)], lwd = 3, border = level1_col[4],add = TRUE)
dev.off()

pdf(paste0(prefix, "collagen_level2_epithelial.pdf"))
size = scale(sqrt(table(l_remove$level2KF[which(l_remove$level1 == "Epithelial")])),center = FALSE)
vioplot(aaa[which((l_remove$level1 == "Epithelial") & (l_remove$kidney_function == "Healthy"))] ~ l_remove$level2KF[which((l_remove$level1 == "Epithelial") & (l_remove$kidney_function == "Healthy"))], col = makeTransparent(level1_col[2], alpha=0.4),side = "left", ylim = c(0,1.8), plotCentre="line", wex = size[c(2,4,6,8,10,12,14,16,18)], lwd = 3, border = level1_col[2])
vioplot(aaa[which((l_remove$level1 == "Epithelial") & (l_remove$kidney_function == "CKD"))] ~ l_remove$level2KF[which((l_remove$level1 == "Epithelial") & (l_remove$kidney_function == "CKD"))], col = makeTransparent(level1_col[2], alpha=0.8),side = "right", ylim = c(0,1.8), plotCentre="line", wex = size[c(1,3,5,7,9,11,13,15,17)], lwd = 3, border = level1_col[2],add = TRUE)
dev.off()



matrisome_set = read.table("https://raw.githubusercontent.com/mahmoudibrahim/KidneyMap/master/assets/public/ecm_genes_human.txt", sep = "\t", header = TRUE)
collagens = matrisome_set$Gene.Symbol[which(matrisome_set$Category %in% c("ECM Glycoproteins"))]
reads_single_phase = logcounts(l_remove)
rownames(reads_single_phase) = toupper(unlist(lapply(strsplit(as.character(rownames(reads_single_phase)), split = ";", fixed = T), function(x) as.character(x[2]))))
reads_single_phase = as.matrix(reads_single_phase[rownames(reads_single_phase) %in% toupper(collagens) ,])
combined_matrix = rbind(reads_single_phase,average=apply(reads_single_phase,2,mean))
cor_matrix = cor(t(combined_matrix))
cor_vector = cor_matrix[,dim(cor_matrix)[1]]
reads_single_phase_restricted = reads_single_phase[rownames(reads_single_phase) %in% names(cor_vector[cor_vector >= 0.1]),]
aaa = apply(reads_single_phase_restricted,2,mean)


pdf(paste0(prefix, "glycoproteins_score_level1.pdf"))
size = scale(sqrt(table(l_remove$level1KF)),center = FALSE)
vioplot(aaa[which(l_remove$kidney_function == "Healthy")] ~ l_remove$level1[which(l_remove$kidney_function == "Healthy")], col = makeTransparent(level1_col, alpha=0.4),side = "left", ylim = c(0,1.3), plotCentre="line", wex = size[c(2,4,6,8,10)], lwd = 3, border = level1_col)
vioplot(aaa[which(l_remove$kidney_function == "CKD")] ~ l_remove$level1[which(l_remove$kidney_function == "CKD")], col = makeTransparent(level1_col,alpha=0.7), side = "right", add = TRUE, plotCentre="line", wex = size[c(1,3,5,7,9)], lwd = 3, border = level1_col)
dev.off()


pdf(paste0(prefix, "glycoproteins_score_level2_mesenchymal.pdf"))
size = scale(sqrt(table(l_remove$level2KF[which(l_remove$level1 == "Mesenchymal")])),center = FALSE)
vioplot(aaa[which((l_remove$level1 == "Mesenchymal") & (l_remove$kidney_function == "Healthy"))] ~ l_remove$level2KF[which((l_remove$level1 == "Mesenchymal") & (l_remove$kidney_function == "Healthy"))], col = makeTransparent(level1_col[4], alpha=0.4),side = "left", ylim = c(0,1.3), plotCentre="line", wex = size[c(2,4,6,8,10)], lwd = 3, border = level1_col[4])
vioplot(aaa[which((l_remove$level1 == "Mesenchymal") & (l_remove$kidney_function == "CKD"))] ~ l_remove$level2KF[which((l_remove$level1 == "Mesenchymal") & (l_remove$kidney_function == "CKD"))], col = makeTransparent(level1_col[4], alpha=0.8),side = "right", ylim = c(0,1.3), plotCentre="line", wex = size[c(1,3,5,7,9)], lwd = 3, border = level1_col[4],add = TRUE)
dev.off()


pdf(paste0(prefix, "glycoproteins_level2_epithelial.pdf"))
size = scale(sqrt(table(l_remove$level2KF[which(l_remove$level1 == "Epithelial")])),center = FALSE)
vioplot(aaa[which((l_remove$level1 == "Epithelial") & (l_remove$kidney_function == "Healthy"))] ~ l_remove$level2KF[which((l_remove$level1 == "Epithelial") & (l_remove$kidney_function == "Healthy"))], col = makeTransparent(level1_col[2], alpha=0.4),side = "left", ylim = c(0,1.3), plotCentre="line", wex = size[c(2,4,6,8,10,12,14,16,18)], lwd = 3, border = level1_col[2])
vioplot(aaa[which((l_remove$level1 == "Epithelial") & (l_remove$kidney_function == "CKD"))] ~ l_remove$level2KF[which((l_remove$level1 == "Epithelial") & (l_remove$kidney_function == "CKD"))], col = makeTransparent(level1_col[2], alpha=0.8),side = "right", ylim = c(0,1.3), plotCentre="line", wex = size[c(1,3,5,7,9,11,13,15,17)], lwd = 3, border = level1_col[2],add = TRUE)
dev.off()


matrisome_set = read.table("https://raw.githubusercontent.com/mahmoudibrahim/KidneyMap/master/assets/public/ecm_genes_human.txt", sep = "\t", header = TRUE)
collagens = matrisome_set$Gene.Symbol[which(matrisome_set$Category %in% c("Proteoglycans"))]
reads_single_phase = logcounts(l_remove)
rownames(reads_single_phase) = toupper(unlist(lapply(strsplit(as.character(rownames(reads_single_phase)), split = ";", fixed = T), function(x) as.character(x[2]))))
reads_single_phase = as.matrix(reads_single_phase[rownames(reads_single_phase) %in% toupper(collagens) ,])
combined_matrix = rbind(reads_single_phase,average=apply(reads_single_phase,2,mean))
cor_matrix = cor(t(combined_matrix))
cor_vector = cor_matrix[,dim(cor_matrix)[1]]
reads_single_phase_restricted = reads_single_phase[rownames(reads_single_phase) %in% names(cor_vector[cor_vector >= 0.1]),]
aaa = apply(reads_single_phase_restricted,2,mean)

pdf(paste0(prefix, "proteoglycans_score_level1.pdf"))
size = scale(sqrt(table(l_remove$level1KF)),center = FALSE)
vioplot(aaa[which(l_remove$kidney_function == "Healthy")] ~ l_remove$level1[which(l_remove$kidney_function == "Healthy")], col = makeTransparent(level1_col, alpha=0.4),side = "left", ylim = c(0,2), plotCentre="line", wex = size[c(2,4,6,8,10)], lwd = 3, border = level1_col)
vioplot(aaa[which(l_remove$kidney_function == "CKD")] ~ l_remove$level1[which(l_remove$kidney_function == "CKD")], col = makeTransparent(level1_col,alpha=0.7), side = "right", add = TRUE, plotCentre="line", wex = size[c(1,3,5,7,9)], lwd = 3, border = level1_col)
dev.off()

pdf(paste0(prefix, "proteoglycans_score_level2_mesenchymal.pdf"))
size = scale(sqrt(table(l_remove$level2KF[which(l_remove$level1 == "Mesenchymal")])),center = FALSE)
vioplot(aaa[which((l_remove$level1 == "Mesenchymal") & (l_remove$kidney_function == "Healthy"))] ~ l_remove$level2KF[which((l_remove$level1 == "Mesenchymal") & (l_remove$kidney_function == "Healthy"))], col = makeTransparent(level1_col[4], alpha=0.4),side = "left", ylim = c(0,2), plotCentre="line", wex = size[c(2,4,6,8,10)], lwd = 3, border = level1_col[4])
vioplot(aaa[which((l_remove$level1 == "Mesenchymal") & (l_remove$kidney_function == "CKD"))] ~ l_remove$level2KF[which((l_remove$level1 == "Mesenchymal") & (l_remove$kidney_function == "CKD"))], col = makeTransparent(level1_col[4], alpha=0.8),side = "right", ylim = c(0,2), plotCentre="line", wex = size[c(1,3,5,7,9)], lwd = 3, border = level1_col[4],add = TRUE)
dev.off()

pdf(paste0(prefix, "proteoglycans_level2_epithelial.pdf"))
size = scale(sqrt(table(l_remove$level2KF[which(l_remove$level1 == "Epithelial")])),center = FALSE)
vioplot(aaa[which((l_remove$level1 == "Epithelial") & (l_remove$kidney_function == "Healthy"))] ~ l_remove$level2KF[which((l_remove$level1 == "Epithelial") & (l_remove$kidney_function == "Healthy"))], col = makeTransparent(level1_col[2], alpha=0.4),side = "left", ylim = c(0,2), plotCentre="line", wex = size[c(2,4,6,8,10,12,14,16,18)], lwd = 3, border = level1_col[2])
vioplot(aaa[which((l_remove$level1 == "Epithelial") & (l_remove$kidney_function == "CKD"))] ~ l_remove$level2KF[which((l_remove$level1 == "Epithelial") & (l_remove$kidney_function == "CKD"))], col = makeTransparent(level1_col[2], alpha=0.8),side = "right", ylim = c(0,2), plotCentre="line", wex = size[c(1,3,5,7,9,11,13,15,17)], lwd = 3, border = level1_col[2],add = TRUE)
dev.off()







#diff exp (Tubule)
ww = c(which(l_remove$level2 %in% c("Proximal Tubule", "Injured Proximal Tubule")))
l_v = l_remove[,ww]
pdf(paste0(prefix, "prox_tubule_top5genes.pdf"))
tube_marks = plotTopMarkerHeat(sg_level3, outs = TRUE, plotheat=FALSE, top_n = 5)[c(20:22,37:40)]
tube_marks = c(tube_marks, "ENSG00000026025.15;VIM")
plotMarkerHeat(logcounts(l_v), l_v$level3, markers = unique(unlist(tube_marks)), clusterGenes=TRUE, gaps = TRUE, averageCells=10^2, clusterGenesK=3, colors = colorRampPalette(c("cornflowerblue","black","gold"))(n=100))
dev.off()



#mesenchymal lineage map
ww = c(which(l_remove$level2 %in% c("Fibroblast", "Myofibroblast","Pericytes")))
l_v = l_remove[,ww]
sg_mesen = sortGenes(counts(l_v), l_v$level3, binarizeMethod="naive", cores = 16)
pp_mesen = getPValues(sg_mesen, numPerm = 20, cores = 16)
marker_mesen = names(which(apply(pp_mesen$adjpval, 1, function(x) any(x < 0.05))))
l_v = scran::computeSumFactors(l_v, sizes = seq(10, 200, 20), clusters = l_v$level3, positive = TRUE)
l_v = scater::normalize(l_v)
l_v = scater::calculateQCMetrics(l_v)
a = plotExprsFreqVsMean(l_v)


#MNN batch effect correction
kk = floor(sqrt(ncol(l_v)) * 1) * ncol(l_v)
newOrder = order(table(l_v$fileID), decreasing = TRUE)
var_genes = marker_mesen
mnn = mnnCorrect(as.matrix(logcounts(l_v[,which(l_v$fileID == (unique(l_v$fileID)[1]))])), as.matrix(logcounts(l_v[,which(l_v$fileID == (unique(l_v$fileID)[2]))])), as.matrix(logcounts(l_v[,which(l_v$fileID == (unique(l_v$fileID)[3]))])), as.matrix(logcounts(l_v[,which(l_v$fileID == (unique(l_v$fileID)[4]))])), as.matrix(logcounts(l_v[,which(l_v$fileID == (unique(l_v$fileID)[5]))])), as.matrix(logcounts(l_v[,which(l_v$fileID == (unique(l_v$fileID)[6]))])), as.matrix(logcounts(l_v[,which(l_v$fileID == (unique(l_v$fileID)[7]))])), as.matrix(logcounts(l_v[,which(l_v$fileID == (unique(l_v$fileID)[8]))])), k = floor(sqrt(ncol(l_v)) * 1), subset.row=var_genes, cos.norm.out=FALSE, order = newOrder, cos.norm.in=TRUE)
l_bc = cbind(mnn$corrected[[1]], mnn$corrected[[2]], mnn$corrected[[3]], mnn$corrected[[4]], mnn$corrected[[5]], mnn$corrected[[6]], mnn$corrected[[7]], mnn$corrected[[8]])
l_bc = SingleCellExperiment(assays=list(logcounts=l_bc))
l_bc$colData = l_v$colData

q = t(logcounts(l_bc[which(rownames(l_bc) %in% marker_mesen),]))
q = t(apply(q,1, function(x) x-mean(x)))
ss = svd(q)
d = ss$d
xx = 1:1000
yy = d[1:1000]
p1 = c(xx[1],yy[1])
p2 = c(xx[length(xx)], yy[length(yy)])
dw = which.max(sapply(1:length(yy), function(x) dist2d(c(xx[x], yy[x]), p1, p2)))
dw = 1:dw
temp = t(apply(ss$u[,dw], 1, function(x) x / (sqrt(sum(x^2))))) #norm based normalization
write.csv(temp, file = paste0(prefix, "dims_mesenchymal.txt"), row.names = FALSE, col.names = FALSE) #this file is saved at reducedDims/SVD/

mapCoords = read.table("https://raw.githubusercontent.com/mahmoudibrahim/KidneyMap/master/assets/reducedDims/UMAP/human_CD10negative_mesenchymal_umapCoords.csv", header = FALSE, sep=",") #learned in python
colors_info_mesenl3 = rep("", length = ncol(l_v))
for (i in 1:length(unique(l_v$level3))) {
	wtemp = which(l_v$level3 == (sort(unique(l_v$level3)))[i])
	colors_info_mesenl3[wtemp] = col_palette_trans[i]
}
colors_info_mesenl2 = rep("", length = ncol(l_v))
for (i in 1:length(unique(l_v$level2))) {
	wtemp = which(l_v$level2 == (sort(unique(l_v$level2)))[i])
	colors_info_mesenl2[wtemp] = col_palette_trans[i]
}
colors_info_patient = rep("", length = ncol(l_v))
for (i in 1:length(unique(l_v$patientID))) {
	wtemp = which((l_v$patientID) == (sort(unique(l_v$patientID)))[i])
	colors_info_patient[wtemp] = col_palette_trans[i]
}
pdf(paste0(prefix, "mesenchymal_specific_umap_level3.pdf"))
plot(mapCoords, col = colors_info_mesenl3, pch = 20)
legend(max(mapCoords[,1]) - 5, min(mapCoords[,2]) + 6, legend = sort(unique(l_v$level3)), pch=21, col = col_palette_trans, pt.bg = col_palette, pt.cex=2, cex=.8, bty="n", ncol = ceiling(length(unique(l_v$level3)) / 20))
dev.off()
pdf(paste0(prefix, "mesenchymal_specific_umap_level2.pdf"))
plot(mapCoords, col = colors_info_mesenl2, pch = 20)
legend(max(mapCoords[,1]) - 5, min(mapCoords[,2]) + 6, legend = sort(unique(l_v$level2)), pch=21, col = col_palette_trans, pt.bg = col_palette, pt.cex=2, cex=.8, bty="n", ncol = ceiling(length(unique(l_v$level2)) / 4))
dev.off()
pdf(paste0(prefix, "mesenchymal_specific_umap_patient.pdf"))
plot(mapCoords, col = colors_info_patient, pch = 20)
legend(max(mapCoords[,1]) - 4, min(mapCoords[,2]) + 6, legend = abbreviate(unique(l_v$patientID), minlength=6), pch=21, col = col_palette, pt.bg = col_palette_trans, pt.cex=2, cex=.8, bty="n", ncol = ceiling(length(unique(l_v$patientID)) / 10))
dev.off()
pdf(paste0(prefix, "mesenchymal_specific_umap_ECMscore.pdf"))
plot(mapCoords, col = color.gradient(scale(aaa[ww]), colors = c("#ECECEC50",rev(c("#d7191c","#abdda4")))), pch = 20)
dev.off()

library(slingshot)
sl = slingshot(mapCoords, l_v$level3, end.clus = "Myofibroblast 1")
pdf(paste0(prefix, "mesenchymal_specific_umap_level2_lineageCurves.pdf"))
plot(mapCoords, col = colors_info_mesenl2, pch = 20)
pst = rep(NA, nrow(mapCoords))
lines(sl, lwd=5, col="#00000095", type = "lineages")
dev.off()
pdf(paste0(prefix, "mesenchymal_specific_umap_level3_lineageCurves.pdf"))
plot(mapCoords, col = colors_info_mesenl3, pch = 20)
lines(sl, lwd=6, col="#00000095", type = "lineages")
dev.off()





#diff exp (fibroblast)
ww = c(which(l_remove$level2 %in% c("Fibroblast", "Myofibroblast","Pericytes")))
l_v = l_remove[,ww]
pdf(paste0(prefix, "fibroblast_top5genes_level2.pdf"))
tube_marks = plotTopMarkerHeat(sg_level2, outs = TRUE, plotheat=FALSE, top_n = 5)[c(9,17,20)]
plotMarkerHeat(logcounts(l_v), l_v$level2, markers = unique(unlist(tube_marks)), clusterGenes=TRUE, gaps = TRUE, averageCells=10^2, clusterGenesK=3, colors = colorRampPalette(c("cornflowerblue","black","gold"))(n=100))
dev.off()
pdf(paste0(prefix, "fibroblast_top10genes_level2.pdf"))
tube_marks = plotTopMarkerHeat(sg_level2, outs = TRUE, plotheat=FALSE, top_n = 10)[c(9,17,20)]
plotMarkerHeat(logcounts(l_v), l_v$level2, markers = unique(unlist(tube_marks)), clusterGenes=TRUE, gaps = TRUE, averageCells=10^2, clusterGenesK=3, colors = colorRampPalette(c("cornflowerblue","black","gold"))(n=100))
dev.off()


#diff exp (endothelial)
ww = c(which(l_remove$level1 %in% c("Endothelial")))
l_v = l_remove[,ww]
pdf(paste0(prefix, "endothelial_top5genes_level2.pdf"))
tube_marks = plotTopMarkerHeat(sg_level2, outs = TRUE, plotheat=FALSE, top_n = 5)[c(1,10,14,21,27,29)]
plotMarkerHeat(logcounts(l_v), l_v$level2, markers = unique(unlist(tube_marks)), clusterGenes=TRUE, gaps = TRUE, averageCells=10^2, clusterGenesK=8, colors = colorRampPalette(c("cornflowerblue","black","gold"))(n=100))
dev.off()
pdf(paste0(prefix, "endothelial_top10genes_level2.pdf"))
tube_marks = plotTopMarkerHeat(sg_level2, outs = TRUE, plotheat=FALSE, top_n = 10)[c(1,10,14,21,27,29)]
plotMarkerHeat(logcounts(l_v), l_v$level2, markers = unique(unlist(tube_marks)), clusterGenes=TRUE, gaps = TRUE, averageCells=10^2, clusterGenesK=8, colors = colorRampPalette(c("cornflowerblue","black","gold"))(n=100))
dev.off()



#diff exp (immune)
ww = c(which(l_remove$level1 %in% c("Immune")))
l_v = l_remove[,ww]
pdf(paste0(prefix, "immune_top5genes_level2.pdf"))
tube_marks = plotTopMarkerHeat(sg_level2, outs = TRUE, plotheat=FALSE, top_n = 5)[c(2,3,6,11,16,19,25)]
plotMarkerHeat(logcounts(l_v), l_v$level2, markers = unique(unlist(tube_marks)), clusterGenes=TRUE, gaps = TRUE, averageCells=10^2, clusterGenesK=8, colors = colorRampPalette(c("cornflowerblue","black","gold"))(n=100))
dev.off()
pdf(paste0(prefix, "immune_top10genes_level2.pdf"))
tube_marks = plotTopMarkerHeat(sg_level2, outs = TRUE, plotheat=FALSE, top_n = 10)[c(2,3,6,11,16,19,25)]
plotMarkerHeat(logcounts(l_v), l_v$level2, markers = unique(unlist(tube_marks)), clusterGenes=TRUE, gaps = TRUE, averageCells=10^2, clusterGenesK=8, colors = colorRampPalette(c("cornflowerblue","black","gold"))(n=100))
dev.off()



#diff exp (epithelial)
ww = c(which(l_remove$level1 %in% c("Epithelial")))
www = c(which(l_remove$level2 %in% c("Injured Proximal Tubule"))) #we don't consider injured proximal tubule here
ww = setdiff(ww,www)
l_v = l_remove[,ww]
pdf(paste0(prefix, "epithelial_top5genes_level2.pdf"))
tube_marks = plotTopMarkerHeat(sg_level2, outs = TRUE, plotheat=FALSE, top_n = 5)[sort(c(23,7,26,5,8,4,13,22))]
plotMarkerHeat(logcounts(l_v), l_v$level2, markers = unique(unlist(tube_marks)), clusterGenes=TRUE, gaps = TRUE, averageCells=10^2, clusterGenesK=8, colors = colorRampPalette(c("cornflowerblue","black","gold"))(n=100), newOrder = c(7,3,8,4,2,1,5,6))
dev.off()
pdf(paste0(prefix, "epithelial_top10genes_level2.pdf"))
tube_marks = plotTopMarkerHeat(sg_level2, outs = TRUE, plotheat=FALSE, top_n = 10)[sort(c(23,7,26,5,8,4,13,22))]
plotMarkerHeat(logcounts(l_v), l_v$level2, markers = unique(unlist(tube_marks)), clusterGenes=TRUE, gaps = TRUE, averageCells=10^2, clusterGenesK=8, colors = colorRampPalette(c("cornflowerblue","black","gold"))(n=100), newOrder = c(7,3,8,4,2,1,5,6))
dev.off()





#core matrisome score
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

#mixture model clustering of matrisome score (vs kidney function)
set.seed(111)
mmm = Mclust(aaa,3)
ttt = table(mmm$classification, l_remove$kidney_function)

pdf(paste0(prefix, "collagen_score_kidney_function.pdf"))
barplot(t(apply(ttt,1,function(x) x / sum(x))) * 100, beside=T, col = col_palette_short[1:3], ylim = c(0,100), yaxp = c(0,100,10))
dev.off()

pdf(paste0(prefix, "collagen_score_histogram.pdf"))
h1 = hist(aaa[which(mmm$classification == 1)], breaks = seq(0,1.5,length.out = 60), plot = FALSE)
h2 = hist(aaa[which(mmm$classification == 2)], breaks = seq(0,1.5,length.out = 60), plot = FALSE)
h3 = hist(aaa[which(mmm$classification == 3)], breaks = seq(0,1.5,length.out = 60), plot = FALSE)
plot(h1, col = makeTransparent(col_palette_short[1],alpha = 0.7), border = col_palette_short[1], ylim = c(0,2500), xaxp = c(0,1.4,14), xlab = "ECM Expression Score")
plot(h2, add = TRUE, col = makeTransparent(col_palette_short[2],alpha = 0.7), border = col_palette_short[2], ylim = c(0,2500))
plot(h3, add = TRUE, col = makeTransparent(col_palette_short[3],alpha=0.7), border = col_palette_short[3], ylim = c(0,2500))
dev.off()

zzz = pdf(paste0(prefix, "umap_ecm_continuous.pdf"))
mapCoords = read.table("https://raw.githubusercontent.com/mahmoudibrahim/KidneyMap/master/assets/reducedDims/UMAP/human_CD10negative_umapCoords.csv", header = FALSE, sep=",")
plot(mapCoords, col = color.gradient(scale(aaa), colors = c("gray",rev(c("#d7191c","#abdda4")))), pch = 20, cex = 0.2)
legend.col(col = colorRampPalette(c("gray",rev(c("#d7191c","#abdda4"))))(100), lev = seq(min(aaa), max(aaa), length.out=100))
zzz = dev.off()
pdf(paste0(prefix,"level3_collagens.pdf"))
matrisome_set = read.table("https://raw.githubusercontent.com/mahmoudibrahim/KidneyMap/master/assets/public/ecm_genes_human.txt", sep = "\t", header = TRUE)
collagens = matrisome_set$Gene.Symbol[which(matrisome_set$Category %in% c("Collagens"))]
plotMarkerHeat(reads_single_phase, sg_level3$inputClass, as.vector(collagens), averageCells= 10^6, newOrder=level3_order, gaps = FALSE, clusterGenes = TRUE, clusterGenesK = 5, colors = colorRampPalette(c("cornflowerblue","black","gold"))(n=100))
dev.off()
pdf(paste0(prefix,"level3_proteoglycans.pdf"))
matrisome_set = read.table("https://raw.githubusercontent.com/mahmoudibrahim/KidneyMap/master/assets/public/ecm_genes_human.txt", sep = "\t", header = TRUE)
collagens = matrisome_set$Gene.Symbol[which(matrisome_set$Category %in% c("Proteoglycans"))]
plotMarkerHeat(reads_single_phase, sg_level3$inputClass, as.vector(collagens), averageCells= 10^6, newOrder=level3_order, gaps = FALSE, clusterGenes = TRUE, clusterGenesK = 4, colors = colorRampPalette(c("cornflowerblue","black","gold"))(n=100))
dev.off()
pdf(paste0(prefix,"level3_ecmGlycoproteins.pdf"), height = 10)
matrisome_set = read.table("https://raw.githubusercontent.com/mahmoudibrahim/KidneyMap/master/assets/public/ecm_genes_human.txt", sep = "\t", header = TRUE)
collagens = matrisome_set$Gene.Symbol[which(matrisome_set$Category %in% c("ECM Glycoproteins"))]
plotMarkerHeat(reads_single_phase, sg_level3$inputClass, as.vector(collagens), averageCells= 10^6, newOrder=level3_order, gaps = FALSE, clusterGenes = TRUE, clusterGenesK = 6, colors = colorRampPalette(c("cornflowerblue","black","gold"))(n=100))
dev.off()




#for ISH
threeGene = t(sg_level3$binary[which(rownames(sg_level3$binary) %in% c("ENSG00000134853.11;PDGFRA", "ENSG00000113721.13;PDGFRB", "ENSG00000108821.13;COL1A1")),])

bin = threeGene
bin[bin > 0] = 1
bin = cbind(bin, rep(1, nrow(bin)))
colnames(bin) = c("PDGFRa","PDGFRb","Col1a1","Other Genes")
bin2 = rowSums(bin[,1:2])
bin3 = rep(0, length(bin2))
bin4 = rep(0, length(bin2))
bin5 = rep(0, length(bin2))
bin3[(bin2 == 1) & (bin[,2] == 1)] = 1
bin4[(bin2 == 1) & (bin[,1] == 1)] = 1
bin5[which(rowSums(bin[,1:2]) == 0)] = 1
bin2[bin2 < 2] = 0
bin2[bin2 == 2] = 1
bin2 = cbind(bin2, bin3, bin4, bin5, bin[,3])
colnames(bin2) = c("PDGFRb_and_PDGFRa","PDGFRb_Only","PDGFRa_Only","Other_Genes","Col1a1")


pdf(paste0(prefix, "percent_heat_withDAPI.pdf"))
agg = aggregate(bin2[,1:4], by = list(bin2[,5]), sum)
agg_per = apply(agg[,2:5], 1, function(x) x/sum(x)) * 100
heatmap.2(agg_per, trace = "none", col = colorRampPalette(c("white","dodgerblue4"))(n=100), scale = "none", cellnote = round(agg_per, 2), notecol = "black", notecex = 1.5, dendrogram = "none", Rowv = NA, Colv = NA, breaks = seq(0,100,length.out=101),key = FALSE, margin = c(10,20))
dev.off()

bb = rep("",nrow(bin2))
bb[which(bin2[,1] == 1)] = "d_Both"
bb[which(bin2[,2] == 1)] = "c_b_only"
bb[which(bin2[,3] == 1)] = "b_a_only"
bb[which(bin2[,4] == 1)] = "a_dapi_only"

threeGene = t(logcounts(l_remove)[which(rownames(sg_level3$binary) %in% c("ENSG00000134853.11;PDGFRA", "ENSG00000113721.13;PDGFRB", "ENSG00000108821.13;COL1A1")),])
pdf(paste0(prefix, "Collagen_by_cell.pdf"))
size = scale(sqrt(table(bb)),center = FALSE)
vioplot(as.matrix(threeGene[,3]) ~ bb, col = col_palette_short[c(1,3,2,4)], plotCentre="line", lwd = 3, border = col_palette_short[c(1,3,2,4)], names = c("other_genes","a_only","b_only","both"), ylab = "Log Collagen Expression", main = "All CD10- Cells (n= 21834)")
dev.off()


##Fin
