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
removethis = names(which(apply(pp$adjpval,2,function(x) length(x[x<0.05])) == 0))
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
var_genes = names(which(apply(pp$adjpval, 1, function(x) any(x < 0.05))))
mnn_u = batchelor::fastMNN(l_remove, batch = l_remove$fileID, d = 40, auto.order = TRUE, subset.row = var_genes, BSPARAM = BiocSingular::IrlbaParam(tol=1e-8))
mnn_u = mnn_u@reducedDims$corrected
temp = t(apply(mnn_u, 1, function(x) x / (sqrt(sum(x^2)))))


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
pdf(paste0(prefix, "regularUMAP_fileID.pdf"))
plot(mapCoords[,1], mapCoords[,2], pch = 20, col = colors_file, cex = 0.3)
legend(max(mapCoords[,1]) - 4, min(mapCoords[,2]) + 20, legend = abbreviate(sort(unique(l_remove$fileID)), minlength=6), pch=21, col = col_palette, pt.bg = col_palette_trans, pt.cex=2, cex=.8, bty="n", ncol = ceiling(length(unique(l_remove$fileID)) / 10))
dev.off()
pdf(paste0(prefix, "regularUMAP_patientID.pdf"))
plot(mapCoords[,1], mapCoords[,2], pch = 20, col = colors_patient, cex = 0.3)
legend(min(mapCoords[,1]), min(mapCoords[,2]) + 10, legend = abbreviate(sort(unique(l_remove$patientID)), minlength=6), pch=21, col = col_palette, pt.bg = col_palette_trans, pt.cex=2, cex=.8, bty="n", ncol = ceiling(length(unique(l_remove$patientID)) / 10))
dev.off()


#corrnet visualization
library(igraph)
library(ggraph)
pv = genesorteR::getPValues(sg_level3, numPerm = 20, cores = 12)
hvg = names(which(apply(pv$adjpval, 1, function(x) any(x < 0.05))))
corr = genesorteR::plotCorrelationHeat(sg_level3, markers = hvg, outs = TRUE)
classNet = corr$corr
diag(classNet) = 0
classNet[which(corr$corr < 0.05)] = 0
net = igraph::graph_from_adjacency_matrix(classNet, weighted = TRUE, mode = "undirected")
net = igraph::simplify(net, edge.attr.comb = "max")
Correlation = E(net)$weight
Percent_of_Cells = as.vector(sg_level3$classProb*100)
cluster_name = names(sg_level3$classProb)
tekka = annot[,3][match(names(sg_level3$classProb), annot[,5])]
coloring = level1_col[as.numeric(factor(tekka))]
set.seed(111)
lay = ggraph::create_layout(net, layout = "fr")
pdf(paste0(prefix, "corrnet.pdf"), width = 10, height = 10)
ggraph(lay) + 
  geom_edge_link(aes(alpha = Correlation), edge_colour = "gray") + 
  geom_node_point(aes(size = Percent_of_Cells), colour = coloring) + 
  geom_node_text(aes(label = cluster_name), repel=TRUE, cex = 3) +
  theme(panel.background = element_blank())
dev.off()




###bmarker heatmap
level3_order = c(2,38,31,34,27,5,28,23,24,25,1,11,50,22,43:48,12,49,32,33,8,9,10,29,30,36,13,6,4,26,39:41,7,3,14:21,42,35,37)
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
vioplot(aaa[which(l_remove$kidney_function == "Healthy")] ~ l_remove$level1[which(l_remove$kidney_function == "Healthy")], col = makeTransparent(level1_col, alpha=0.4),side = "left", ylim = c(0,1.6), plotCentre="line", wex = size[c(2,4,6,8,10)], lwd = 3, border = level1_col, xlim = c(-1,7))
vioplot(aaa[which(l_remove$kidney_function == "CKD")] ~ l_remove$level1[which(l_remove$kidney_function == "CKD")], col = makeTransparent(level1_col,alpha=0.7), side = "right", add = TRUE, plotCentre="line", wex = size[c(1,3,5,7,9)], lwd = 3, border = level1_col, xlim = c(-1,7))
dev.off()
pdf(paste0(prefix, "_core_matrisome_score_level2_mesenchymal.pdf"))
size = scale(sqrt(table(l_remove$level2KF[which(l_remove$level1 == "Mesenchymal")])),center = FALSE)
vioplot(aaa[which((l_remove$level1 == "Mesenchymal") & (l_remove$kidney_function == "Healthy"))] ~ l_remove$level2[which((l_remove$level1 == "Mesenchymal") & (l_remove$kidney_function == "Healthy"))], col = makeTransparent(level1_col[4], alpha=0.4),side = "left", ylim = c(0,1.6), plotCentre="line", wex = size[c(2,4,6,8)], lwd = 3, border = level1_col[4], names = abbreviate(levels(factor(l_remove$level2[which(l_remove$level1 == "Mesenchymal")]))))
vioplot(aaa[which((l_remove$level1 == "Mesenchymal") & (l_remove$kidney_function == "CKD"))] ~ l_remove$level2[which((l_remove$level1 == "Mesenchymal") & (l_remove$kidney_function == "CKD"))], col = makeTransparent(level1_col[4], alpha=0.8),side = "right", ylim = c(0,1.6), plotCentre="line", wex = size[c(1,3,5,7)], lwd = 3, border = level1_col[4],add = TRUE)
dev.off()
pdf(paste0(prefix, "_core_matrisome_score_level3_mesenchymal.pdf"), width=10)
size = scale(sqrt(table(l_remove$level3KF[which(l_remove$level1 == "Mesenchymal")])),center = FALSE)
vioplot(aaa[which((l_remove$level1 == "Mesenchymal") & (l_remove$kidney_function == "Healthy"))] ~ l_remove$level3[which((l_remove$level1 == "Mesenchymal") & (l_remove$kidney_function == "Healthy"))], col = makeTransparent(level1_col[4], alpha=0.4),side = "left", ylim = c(0,1.6), plotCentre="line", wex = size[c(2,4,6,8,10,12,14,16)], lwd = 3, border = level1_col[4], names = abbreviate(levels(factor(l_remove$level3[which(l_remove$level1 == "Mesenchymal")]))))
vioplot(aaa[which((l_remove$level1 == "Mesenchymal") & (l_remove$kidney_function == "CKD"))] ~ l_remove$level3[which((l_remove$level1 == "Mesenchymal") & (l_remove$kidney_function == "CKD"))], col = makeTransparent(level1_col[4], alpha=0.8),side = "right", ylim = c(0,1.6), plotCentre="line", wex = size[c(1,3,5,7,9,11,13,15)], lwd = 3, border = level1_col[4],add = TRUE)
dev.off()
pdf(paste0(prefix, "_core_matrisome_score_level2_epithelial.pdf"))
size = scale(sqrt(table(l_remove$level2KF[which(l_remove$level1 == "Epithelial")])),center = FALSE)
vioplot(aaa[which((l_remove$level1 == "Epithelial") & (l_remove$kidney_function == "Healthy"))] ~ l_remove$level2KF[which((l_remove$level1 == "Epithelial") & (l_remove$kidney_function == "Healthy"))], col = makeTransparent(level1_col[2], alpha=0.4),side = "left", ylim = c(0,1.6), plotCentre="line", wex = size[c(2,4,6,8,10,12,14,16,18,20)], lwd = 3, border = level1_col[2])
vioplot(aaa[which((l_remove$level1 == "Epithelial") & (l_remove$kidney_function == "CKD"))] ~ l_remove$level2KF[which((l_remove$level1 == "Epithelial") & (l_remove$kidney_function == "CKD"))], col = makeTransparent(level1_col[2], alpha=0.8),side = "right", ylim = c(0,1.6), plotCentre="line", wex = size[c(1,3,5,7,9,11,13,15,17,19)], lwd = 3, border = level1_col[2],add = TRUE)
dev.off()
#pval
mesen_p = wilcox.test(aaa[which(l_remove$level1 == "Mesenchymal")] ~ l_remove$kidney_function[which(l_remove$level1 == "Mesenchymal")])$p.value
immune_p = wilcox.test(aaa[which(l_remove$level1 == "Immune")] ~ l_remove$kidney_function[which(l_remove$level1 == "Immune")])$p.value
epith_p = wilcox.test(aaa[which(l_remove$level1 == "Epithelial")] ~ l_remove$kidney_function[which(l_remove$level1 == "Epithelial")])$p.value
endoth_p = wilcox.test(aaa[which(l_remove$level1 == "Endothelial")] ~ l_remove$kidney_function[which(l_remove$level1 == "Endothelial")])$p.value
pvals_level1 = p.adjust(c(mesen_p, immune_p, epith_p, endoth_p), method = "bonferroni")
#pval 2
f2_p = wilcox.test(aaa[which(l_remove$level3 == "Fibroblast 2")] ~ l_remove$kidney_function[which(l_remove$level3 == "Fibroblast 2")])$p.value
f4_p = wilcox.test(aaa[which(l_remove$level3 == "Fibroblast 4")] ~ l_remove$kidney_function[which(l_remove$level3 == "Fibroblast 4")])$p.value
f6_p = wilcox.test(aaa[which(l_remove$level3 == "Fibroblast 6")] ~ l_remove$kidney_function[which(l_remove$level3 == "Fibroblast 6")])$p.value
mf1a_p = wilcox.test(aaa[which(l_remove$level3 == "Myofibroblast 1a")] ~ l_remove$kidney_function[which(l_remove$level3 == "Myofibroblast 1a")])$p.value
mf1b_p = wilcox.test(aaa[which(l_remove$level3 == "Myofibroblast 1b")] ~ l_remove$kidney_function[which(l_remove$level3 == "Myofibroblast 1b")])$p.value
peri1_p = wilcox.test(aaa[which(l_remove$level3 == "Pericytes 1")] ~ l_remove$kidney_function[which(l_remove$level3 == "Pericytes 1")])$p.value
peri2_p = wilcox.test(aaa[which(l_remove$level3 == "Pericytes 2")] ~ l_remove$kidney_function[which(l_remove$level3 == "Pericytes 2")])$p.value
smc_p = wilcox.test(aaa[which(l_remove$level3 == "Vascular Smooth Muscle Cells")] ~ l_remove$kidney_function[which(l_remove$level3 == "Vascular Smooth Muscle Cells")])$p.value
pvals_mesen = p.adjust(c(f2_p, f4_p, f6_p, mf1a_p, mf1b_p, peri1_p, peri2_p, smc_p), method = "bonferroni")
names(pvals_mesen) = c("Fib2","Fib4","Fib6","MF1a","MF1b","Peri1","Peri2","SMC")

#collagen score
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
pdf(paste0(prefix, "collagen_level2_epithelial.pdf"), width = 10)
size = scale(sqrt(table(l_remove$level2KF[which(l_remove$level1 == "Epithelial")])),center = FALSE)
vioplot(aaa[which((l_remove$level1 == "Epithelial") & (l_remove$kidney_function == "Healthy"))] ~ l_remove$level2KF[which((l_remove$level1 == "Epithelial") & (l_remove$kidney_function == "Healthy"))], col = makeTransparent(level1_col[2], alpha=0.4),side = "left", ylim = c(0,2), plotCentre="line", wex = size[c(2,4,6,8,10,12,14,16,18,20)], lwd = 3, border = level1_col[2], names = abbreviate(names(size[c(2,4,6,8,10,12,14,16,18,20),])), xlab = "", ylab = "")
vioplot(aaa[which((l_remove$level1 == "Epithelial") & (l_remove$kidney_function == "CKD"))] ~ l_remove$level2KF[which((l_remove$level1 == "Epithelial") & (l_remove$kidney_function == "CKD"))], col = makeTransparent(level1_col[2], alpha=0.8),side = "right", ylim = c(0,2), plotCentre="line", wex = size[c(1,3,5,7,9,11,13,15,17,19)], lwd = 3, border = level1_col[2],add = TRUE)
dev.off()

#glycoprotein score
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
pdf(paste0(prefix, "glycoproteins_level2_epithelial.pdf"), width = 10)
size = scale(sqrt(table(l_remove$level2KF[which(l_remove$level1 == "Epithelial")])),center = FALSE)
vioplot(aaa[which((l_remove$level1 == "Epithelial") & (l_remove$kidney_function == "Healthy"))] ~ l_remove$level2KF[which((l_remove$level1 == "Epithelial") & (l_remove$kidney_function == "Healthy"))], col = makeTransparent(level1_col[2], alpha=0.4),side = "left", ylim = c(0,2), plotCentre="line", wex = size[c(2,4,6,8,10,12,14,16,18,20)], lwd = 3, border = level1_col[2], names = abbreviate(names(size[c(2,4,6,8,10,12,14,16,18,20),])), xlab = "", ylab = "")
vioplot(aaa[which((l_remove$level1 == "Epithelial") & (l_remove$kidney_function == "CKD"))] ~ l_remove$level2KF[which((l_remove$level1 == "Epithelial") & (l_remove$kidney_function == "CKD"))], col = makeTransparent(level1_col[2], alpha=0.8),side = "right", ylim = c(0,2), plotCentre="line", wex = size[c(1,3,5,7,9,11,13,15,17,19)], lwd = 3, border = level1_col[2],add = TRUE)
dev.off()

#proteoglycan score
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
pdf(paste0(prefix, "proteoglycans_level2_epithelial.pdf"), width = 10)
size = scale(sqrt(table(l_remove$level2KF[which(l_remove$level1 == "Epithelial")])),center = FALSE)
vioplot(aaa[which((l_remove$level1 == "Epithelial") & (l_remove$kidney_function == "Healthy"))] ~ l_remove$level2KF[which((l_remove$level1 == "Epithelial") & (l_remove$kidney_function == "Healthy"))], col = makeTransparent(level1_col[2], alpha=0.4),side = "left", ylim = c(0,2), plotCentre="line", wex = size[c(2,4,6,8,10,12,14,16,18,20)], lwd = 3, border = level1_col[2], names = abbreviate(names(size[c(2,4,6,8,10,12,14,16,18,20),])), xlab = "", ylab = "")
vioplot(aaa[which((l_remove$level1 == "Epithelial") & (l_remove$kidney_function == "CKD"))] ~ l_remove$level2KF[which((l_remove$level1 == "Epithelial") & (l_remove$kidney_function == "CKD"))], col = makeTransparent(level1_col[2], alpha=0.8),side = "right", ylim = c(0,2), plotCentre="line", wex = size[c(1,3,5,7,9,11,13,15,17,19)], lwd = 3, border = level1_col[2],add = TRUE)
dev.off()







#diff exp (Tubule)
ww = c(which(l_remove$level3 %in% c("Proximal Tubule", "Injured Proximal tubule")))
l_v = l_remove[,ww]
pdf(paste0(prefix, "prox_tubule_top5genes.pdf"))
tube_marks = plotTopMarkerHeat(sg_level3, outs = TRUE, plotheat=FALSE, top_n = 20)[c(13,36)]
tube_marks = c(tube_marks, "ENSG00000026025.15;VIM")
plotMarkerHeat(logcounts(l_v), l_v$level3, markers = unique(unlist(tube_marks)), clusterGenes=TRUE, gaps = TRUE, averageCells=10^2, clusterGenesK=3, colors = colorRampPalette(c("cornflowerblue","black","gold"))(n=100))
dev.off()
#GO-BP
ng = plotTopMarkerHeat(sg_level3, outs = TRUE, plotheat=FALSE, top_n = 100)[c(13,36)]
ng[[1]] = toupper(unlist(lapply(strsplit(as.character(ng[[1]]), split = ";", fixed = T), function(x) as.character(x[2]))))
ng[[2]] = toupper(unlist(lapply(strsplit(as.character(ng[[2]]), split = ";", fixed = T), function(x) as.character(x[2]))))
library(tibble)
library(clusterProfiler)
m_t2g =  as.tibble(read.gmt("/home/mibrahim/Dropbox/dev/kidneyMap/public/assets/public/c5.bp.v7.1.symbols.gmt"))
en = mclapply(ng, function(x) enricher(x, TERM2GENE=m_t2g, minGSSize = 100, maxGSSize = 500, pvalueCutoff = 0.05), mc.cores = 12)
en2 = list()
for (i in 1:length(en)) {
	en2[[i]] =  en[[i]]@result$ID[1:10]
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
pdf(paste0(prefix, "bp_pt.pdf"), width = 10, height = 10)
heatmap.2(t(en4), trace = "none", scale = "none", col = colorRampPalette(c("white","dodgerblue4"))(n=100), margins = c(20,30), density = "none", Colv = NA, key = FALSE, labCol = names(ng), cexRow =0.5, cexCol = 0.8)
dev.off()




#mesenchymal lineage map
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


#start here
ww = c(which(l_remove$level2 %in% c("Fibroblast", "Myofibroblast","Pericytes")))
l_v = l_remove[,ww]
sg_mesen = sortGenes(counts(l_v), l_v$level3, binarizeMethod="naive", cores = 16)
pp_mesen = getPValues(sg_mesen, numPerm = 20, cores = 16)
marker_mesen = names(which(apply(pp_mesen$adjpval, 1, function(x) any(x < 0.05))))
l_v = scran::computeSumFactors(l_v, sizes = seq(10, 200, 20), clusters = l_v$level3, positive = TRUE)
l_v = scater::normalize(l_v)
l_v = scater::calculateQCMetrics(l_v)
a = plotExprsFreqVsMean(l_v)
var_genes = marker_mesen
mnn_u = batchelor::fastMNN(l_v, batch = l_v$fileID, d = 40, auto.order = TRUE, subset.row = var_genes, BSPARAM = BiocSingular::IrlbaParam(tol=1e-8))
mnn_u = mnn_u@reducedDims$corrected
temp = t(apply(ss$u[,dw], 1, function(x) x / (sqrt(sum(x^2)))))
write.csv(temp, file = paste0(prefix, "dims_mesenchymal.txt"), row.names = FALSE, col.names = FALSE) #this file is saved at reducedDims/SVD/
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


#umap based
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
sl = slingshot(mapCoords, l_v$level3, end.clus = "Myofibroblast 1a", start.clus = "Pericytes 1")
st1 = rowMeans(slingPseudotime(sl), na.rm = TRUE)
pst = st1
col_pst = color.gradient(scale(pst), colors = c("#ECECEC50",rev(c("#d7191c","#abdda4"))))
pdf(paste0(prefix, "mesenchymal_specific_umap_level2_lineageCurves.pdf"))
plot(mapCoords, col = colors_info_mesenl2, pch = 20)
pst = rep(NA, nrow(mapCoords))
lines(sl, lwd=5, col="#00000095", type = "curves")
dev.off()
pdf(paste0(prefix, "mesenchymal_specific_umap_pseudoTime_lineageCurves.pdf"))
plot(mapCoords, col = col_pst, pch = 20)
lines(sl, lwd=6, col="#00000095", type = "curves")
dev.off()
pdf(paste0(prefix, "mesenchymal_specific_umap_pseudoTime.pdf"))
plot(mapCoords, col = col_pst, pch = 20)
dev.off()
pdf(paste0(prefix, "mesenchymal_specific_umap_level3_lineageCurves.pdf"))
plot(mapCoords, col = colors_info_mesenl3, pch = 20)
lines(sl, lwd=6, col="#00000095", type = "curves")
dev.off()


#diffusion map
library(destiny)
set.seed(111)
df = DiffusionMap(mnn_u, k = 25, n_eigs = 3)
write.csv(df@eigenvectors, file = paste0(prefix, "dims_diffusionmap.txt"), row.names = FALSE, col.names = FALSE) #fle saved at reducedDims/DM
sl2 = slingshot(df@eigenvectors, l_v$level3, start.clus = "Myofibroblast 1a")
st2 = rowMeans(slingPseudotime(sl2), na.rm = TRUE)
col_pst = color.gradient(scale(-st2), colors = c("#ECECEC50",rev(c("#d7191c","#abdda4"))))
curv = getCurves(sl2)
library(plot3D)
pdf(paste0(prefix, "pseudotime_specific_dm_level3.pdf"))
scatter3D(df2[,1], df2[,2], df2[,3], col = colors_info_mesenl3, pch = 20, cex = 0.5, colvar = NULL, phi = 30, theta = 230, bty = "n")
dev.off()
pdf(paste0(prefix, "pseudotime_specific_dm_level3_bty.pdf"))
scatter3D(df2[,1], df2[,2], df2[,3], col = colors_info_mesenl3, pch = 20, cex = 0.5, colvar = NULL, phi = 30, theta = 230)
dev.off()
pdf(paste0(prefix, "pseudotime_specific_dm_level2.pdf"))
scatter3D(df2[,1], df2[,2], df2[,3], col = colors_info_mesenl2, pch = 20, cex = 0.5, colvar = NULL, bty = "n", phi = 30, theta = 230)
dev.off()
pdf(paste0(prefix, "pseudotime_specific_dm_pseudotime.pdf"))
scatter3D(df2[,1], df2[,2], df2[,3], col = col_pst, pch = 19, cex = 0.8, colvar = NULL, bty = "n", phi = 30, theta = 230)
dev.off()
pdf(paste0(prefix, "pseudotime_specific_dm_ecm.pdf"))
scatter3D(df2[,1], df2[,2], df2[,3], col = color.gradient(scale(aaa[ww]), colors = c("#ECECEC50",rev(c("#d7191c","#abdda4")))), pch = 19, cex = 1, colvar = NULL, bty = "n", phi = 30, theta = 230)
dev.off()
pdf(paste0(prefix, "pseudotime_specific_dm_kidneyFunction.pdf"))
scatter3D(df2[,1], df2[,2], df2[,3], col = colors_info_kf[ww], pch = 19, cex = 0.8, colvar = NULL, bty = "n", phi = 30, theta = 230)
dev.off()

hend_markers = unique(unlist(plotTopMarkerHeat(sg_mesen, outs = TRUE, plotheat=FALSE))) #plot some markers
for (i in 1:length(hend_markers)) {
	if(any(names1 == hend_markers[i])) {
		aaax = as.vector(logcounts(l_remove)[which(names1 == hend_markers[i]),])
		aaax = aaax[ww]
		zzz = pdf(paste0(prefix, "pseudotime_specific_dm_", hend_markers[i], ".pdf"))
		scatter3D(df2[,1], df2[,2], df2[,3], col = color.gradient(scale(aaax), colors = c("#ECECEC50",rev(c("#d7191c","#abdda4")))), pch = 19, cex = 1, colvar = NULL, bty = "n", phi = 30, theta = 230, alpha = 1)
		
		zzz = dev.off()
	}
}










#diff exp (fibroblast)
ww = c(which(l_remove$level2 %in% c("Fibroblast", "Myofibroblast","Pericytes","Smooth Muscle Cells")))
l_v = l_remove[,ww]
pdf(paste0(prefix, "fibroblast_kindeyFunction_barplot.pdf"))
par(mar=c(5,14,4,2))
barplot(table(l_v$kidney_function, l_v$level2), col = rev(kf_col_trans), horiz = TRUE, las = 2, xaxp = c(0,1200,10), xlab = "# Cells")
dev.off()
pdf(paste0(prefix, "fibroblast_top10genes_level2.pdf"))
tube_marks = plotTopMarkerHeat(sg_level2, outs = TRUE, plotheat=FALSE, top_n = 10)[c(9,17,19)]
plotMarkerHeat(logcounts(l_v), l_v$level2, markers = unique(unlist(tube_marks)), clusterGenes=TRUE, gaps = TRUE, averageCells=10^2, clusterGenesK=3, colors = colorRampPalette(c("cornflowerblue","black","gold"))(n=100))
dev.off()
dev.off()
#GO-BP
ng = plotTopMarkerHeat(sg_level3, outs = TRUE, plotheat=FALSE, top_n = 100)[c(8,9,10,29,30,32,33,49)]
ng[[1]] = toupper(unlist(lapply(strsplit(as.character(ng[[1]]), split = ";", fixed = T), function(x) as.character(x[2]))))
ng[[2]] = toupper(unlist(lapply(strsplit(as.character(ng[[2]]), split = ";", fixed = T), function(x) as.character(x[2]))))
ng[[3]] = toupper(unlist(lapply(strsplit(as.character(ng[[3]]), split = ";", fixed = T), function(x) as.character(x[2]))))
ng[[4]] = toupper(unlist(lapply(strsplit(as.character(ng[[4]]), split = ";", fixed = T), function(x) as.character(x[2]))))
ng[[5]] = toupper(unlist(lapply(strsplit(as.character(ng[[5]]), split = ";", fixed = T), function(x) as.character(x[2]))))
ng[[6]] = toupper(unlist(lapply(strsplit(as.character(ng[[6]]), split = ";", fixed = T), function(x) as.character(x[2]))))
ng[[7]] = toupper(unlist(lapply(strsplit(as.character(ng[[7]]), split = ";", fixed = T), function(x) as.character(x[2]))))
ng[[8]] = toupper(unlist(lapply(strsplit(as.character(ng[[8]]), split = ";", fixed = T), function(x) as.character(x[2]))))
library(tibble)
library(clusterProfiler)
m_t2g =  as.tibble(read.gmt("/home/mibrahim/Dropbox/dev/kidneyMap/public/assets/public/c5.bp.v7.1.symbols.gmt"))
en = mclapply(ng, function(x) enricher(x, TERM2GENE=m_t2g, minGSSize = 100, maxGSSize = 500, pvalueCutoff = 0.05), mc.cores = 12)
en2 = list()
for (i in 1:length(en)) {
	en2[[i]] =  en[[i]]@result$ID[1:10]
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
pdf(paste0(prefix, "bp_mesenchymal.pdf"), width = 10, height = 10)
heatmap.2(t(en4), trace = "none", scale = "none", col = colorRampPalette(c("white","dodgerblue4"))(n=100), margins = c(20,30), density = "none", key = FALSE, labCol = names(ng), cexRow =0.5, cexCol = 0.8)
dev.off()







#diff exp (endothelial)
ww = c(which(l_remove$level1 %in% c("Endothelial")))
l_v = l_remove[,ww]
pdf(paste0(prefix, "endothelial_kindeyFunction_barplot.pdf"))
par(mar=c(5,14,4,2))
barplot(table(l_v$kidney_function, l_v$level2), col = rev(kf_col_trans), horiz = TRUE, las = 2, xaxp = c(0,16000,10), xlab = "# Cells")
dev.off()
pdf(paste0(prefix, "endothelial_top10genes_level2.pdf"))
tube_marks = plotTopMarkerHeat(sg_level2, outs = TRUE, plotheat=FALSE, top_n = 10)[c(1,10,11,14,28,29)]
plotMarkerHeat(logcounts(l_v), l_v$level2, markers = unique(unlist(tube_marks)), clusterGenes=TRUE, gaps = TRUE, averageCells=10^2, clusterGenesK=8, colors = colorRampPalette(c("cornflowerblue","black","gold"))(n=100))
dev.off()


#diff exp (immune)
ww = c(which(l_remove$level1 %in% c("Immune")))
l_v = l_remove[,ww]
pdf(paste0(prefix, "immune_kindeyFunction_barplot.pdf"))
par(mar=c(5,14,4,2))
barplot(table(l_v$kidney_function, l_v$level2), col = rev(kf_col_trans), horiz = TRUE, las = 2, xaxp = c(0,3200,10), xlab = "# Cells")
dev.off()
pdf(paste0(prefix, "immune_top10genes_level2.pdf"))
tube_marks = plotTopMarkerHeat(sg_level2, outs = TRUE, plotheat=FALSE, top_n = 10)[c(2,3,6,15,16,18,20,25)]
plotMarkerHeat(logcounts(l_v), l_v$level2, markers = unique(unlist(tube_marks)), clusterGenes=TRUE, gaps = TRUE, averageCells=10^2, clusterGenesK=9, colors = colorRampPalette(c("cornflowerblue","black","gold"))(n=100))
dev.off()



#diff exp (epithelial)
ww = c(which(l_remove$level1 %in% c("Epithelial")))
www = c(which(l_remove$level2 %in% c("Injured Tubule")))
ww = setdiff(ww,www)
l_v = l_remove[,ww]

pdf(paste0(prefix, "epithelial_kindeyFunction_barplot.pdf"))
par(mar=c(5,14,4,2))
barplot(table(l_v$kidney_function, l_v$level2), col = rev(kf_col_trans), horiz = TRUE, las = 2, xaxp = c(0,2400,10), xlab = "# Cells")
dev.off()
pdf(paste0(prefix, "epithelial_top10genes_level2.pdf"))
tube_marks = plotTopMarkerHeat(sg_level2, outs = TRUE, plotheat=FALSE, top_n = 10)[c(4,5,7,8,13,21,22,26,27)]
plotMarkerHeat(logcounts(l_v), l_v$level2, markers = unique(unlist(tube_marks)), clusterGenes=TRUE, gaps = TRUE, averageCells=10^2, clusterGenesK=8, colors = colorRampPalette(c("cornflowerblue","black","gold"))(n=100))
dev.off()





#mixture model clustering of matrisome score (vs kidney function)
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


set.seed(111)
mmm = Mclust(aaa,3) #cluster
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




#for ISH (collagen by cell)
threeGene = t(sg_level3$binary[which(rownames(sg_level3$binary) %in% c("ENSG00000134853.11;PDGFRA", "ENSG00000113721.13;PDGFRB", "ENSG00000108821.13;COL1A1")),])
bin = threeGene
bin[bin > 0] = 1
bin = cbind(bin, rep(1, nrow(bin)))
colnames(bin) = c("PDGFRa","PDGFRb","Col1a1","Other Genes")
library(UpSetR)
pdf(paste0(prefix, "UpSet.pdf"))
upset(data.frame(as.matrix(bin)), nsets = 100, group.by = "degree", keep.order=TRUE, main.bar.color=col_palette_short[2], sets.bar.color=col_palette_short[1], decreasing = c(TRUE,TRUE))
dev.off()
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
library(gplots)
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
vioplot(as.matrix(threeGene[,3]) ~ bb, col = col_palette_short[c(1,3,2,4)], plotCentre="line", lwd = 3, border = col_palette_short[c(1,3,2,4)], names = c("other_genes","a_only","b_only","both"), ylab = "Log Collagen Expression", main = "All CD10- Cells (n= 51849)")
dev.off()
#tests
ww = which(bb %in% c("a_dapi_only","d_Both"))
tt1 = wilcox.test(threeGene[ww,3] ~ bb[ww])$p.value
ww = which(bb %in% c("b_a_only","d_Both"))
tt2 = wilcox.test(threeGene[ww,3] ~ bb[ww])$p.value
ww = which(bb %in% c("c_b_only","d_Both"))
tt3 = wilcox.test(threeGene[ww,3] ~ bb[ww])$p.value
ww = which(bb %in% c("a_dapi_only","b_a_only"))
tt4 = wilcox.test(threeGene[ww,3] ~ bb[ww])$p.value
ww = which(bb %in% c("c_b_only","b_a_only"))
tt5 = wilcox.test(threeGene[ww,3] ~ bb[ww])$p.value
ww = which(bb %in% c("a_dapi_only","c_b_only"))
tt6 = wilcox.test(threeGene[ww,3] ~ bb[ww])$p.value
pvals = p.adjust(c(tt1,tt2,tt3,tt4,tt5,tt6), method = "bonferroni")
names(pvals) = c("dapi_vs_both","a_vs_both","b_vs_both","dapi_vs_a","b_vs_a","dapi_vs_b")


#for ISH, ECM by Cell
threeGene = t(logcounts(l_remove)[which(rownames(sg_level3$binary) %in% c("ENSG00000134853.11;PDGFRA", "ENSG00000113721.13;PDGFRB", "ENSG00000108821.13;COL1A1")),])
threeGene[,3] = aaa #aaa is ECM score
pdf(paste0(prefix, "ECM_by_cell.pdf"))
size = scale(sqrt(table(bb)),center = FALSE)
vioplot(as.matrix(threeGene[,3]) ~ bb, col = col_palette_short[c(1,3,2,4)], plotCentre="line", lwd = 3, border = col_palette_short[c(1,3,2,4)], names = c("other_genes","a_only","b_only","both"), ylab = "ECM Score", main = "All CD10- Cells (n= 51849)")
dev.off()
ww = which(bb %in% c("a_dapi_only","d_Both"))
tt1 = wilcox.test(threeGene[ww,3] ~ bb[ww])$p.value
ww = which(bb %in% c("b_a_only","d_Both"))
tt2 = wilcox.test(threeGene[ww,3] ~ bb[ww])$p.value
ww = which(bb %in% c("c_b_only","d_Both"))
tt3 = wilcox.test(threeGene[ww,3] ~ bb[ww])$p.value
ww = which(bb %in% c("a_dapi_only","b_a_only"))
tt4 = wilcox.test(threeGene[ww,3] ~ bb[ww])$p.value
ww = which(bb %in% c("c_b_only","b_a_only"))
tt5 = wilcox.test(threeGene[ww,3] ~ bb[ww])$p.value
ww = which(bb %in% c("a_dapi_only","c_b_only"))
tt6 = wilcox.test(threeGene[ww,3] ~ bb[ww])$p.value
pvals = p.adjust(c(tt1,tt2,tt3,tt4,tt5,tt6), method = "bonferroni")
names(pvals) = c("dapi_vs_both","a_vs_both","b_vs_both","dapi_vs_a","b_vs_a","dapi_vs_b")



#write output
writeMM(counts(l_remove), file = paste0(prefix, "UMI_counts.mtx"))
names1 = toupper(unlist(lapply(strsplit(as.character(rownames(l_remove)), split = ";", fixed = T), function(x) as.character(x[2]))))
names2 = toupper(unlist(lapply(strsplit(as.character(rownames(l_remove)), split = ";", fixed = T), function(x) as.character(x[1]))))
names2 = toupper(unlist(lapply(strsplit(as.character(names2), split = ".", fixed = T), function(x) as.character(x[1]))))
rowdat = cbind(names1, names2)
colnames(rowdat) = c("Gene.Symbol","ENSEMBL.ID")
write.table(rowdat, file = paste0(prefix, "UMI_counts_rowData.txt"), sep = "\t", row.names=FALSE, quote = FALSE)

coldat = cbind(l_remove$level1, l_remove$level2, l_remove$level3, l_remove$kidney_function, l_remove$patientID)
colnames(coldat) = c("Annotation.Level.1","Annotation.Level.2","Annotation.Level.3","Kidney.Function","Patient ID")
write.table(coldat, file = paste0(prefix, "UMI_counts_colData.txt"), sep = "\t", row.names=FALSE, quote = FALSE)

##Fin
