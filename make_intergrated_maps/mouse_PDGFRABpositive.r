########################################################################
# mouse_PDGFRABpositive.r
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


prefix = "mouse_PDGFRABpositive_"



###########################################################################
########################START ANALYSIS#####################################
###########################################################################

#load data
load("mouse_PDGFRAB.RData") #R object containing aggregated data, saved as SingleCellExperiment. Variable name: l
lx = l #save full matrix for later


#gene filtering
cutGenes = log10(floor(ncol(l) * 0.001))
choose_genes = which((log10(Matrix::rowSums(counts(l)))) >= cutGenes)
l = l[choose_genes, ]




###first stuff [merging similar clusters]
sg = sortGenes(counts(l), colData(l)$class_info, binarizeMethod="naive", cores = 16)
write.files(sg, prefix = paste0(prefix, "genesorteRouts_noMerging_noRemove"))
marks = plotTopMarkerHeat(sg, top_n=100, plotheat=FALSE, outs = TRUE)
pdf(paste0(prefix, "correlationBetweenClusters_noMerging_noRemove.pdf"), height = 12, width = 12)
corr = plotCorrelationHeat(sg, displayNumbers=F, markers=getMarkers(sg)$markers, outs = TRUE)
dev.off()

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
removethis = which(class_info_merging %in% removethis)



#reform matrix and refilter genes
class_info_merging_remove = class_info_merging[-removethis]
l_remove = lx
l_remove = l_remove[,-removethis]
l_remove = calculateQCMetrics(l_remove)
a = plotExprsFreqVsMean(l_remove)
colData(l_remove)$final_classes = class_info_merging_remove
cutGenes = log10(floor(ncol(l_remove) * 0.001)) #on average at least 1UMI in 0.1% of cells 
choose_genes = which((log10(Matrix::rowSums(counts(l_remove)))) >= cutGenes)
l_remove = l_remove[choose_genes, ]




#make table of healthy vs disease
pdf(paste0(prefix, "percent_healthy_disease.pdf"))
barplot((table(l_remove$kidney_function) / ncol(l_remove)) * 100, col = kf_col, ylim = c(0,90), yaxp = c(0,90,18), ylab = "% of Cells")
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
annot = read.table("https://raw.githubusercontent.com/mahmoudibrahim/KidneyMap/master/assets/clusterInfo/mouse_PDGFRABpositive.txt", sep = "\t")
annot = annot[match(names(KFlow), annot$V1),]
annot$V7 = KFlow
annot$V8 = log(KFlow) - log((table(l_remove$kidney_function) / ncol(l_remove))[2])
nums = table(class_info_merging_remove)
annot = annot[match(names(nums), annot$V1),]
annot$V9 = as.numeric(nums)
annot$V10 = as.numeric(log10(nums))
colnames(annot) = c("ID", "Label", "Annotation 1", "Annotation 2", "Annotation 3", "Percent Low KF Cells", "Log Enrichment of Low KF Cells", "Number of Cells", "log10 Number of Cells")
##write cluster annotation information
write.csv(annot, paste0(prefix, "_cluster_annotation_info_withAddData.csv"), row.names = FALSE)
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





#set colors
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
	colors_info_kf[wtemp] = (kf_col_trans)[i]
}

colors_file = rep("", length = ncol(l_remove))
for (i in 1:length(unique(l_remove$fileID))) {
	wtemp = which((l_remove$fileID) == (sort(unique(l_remove$fileID)))[i])
	colors_file[wtemp] = col_palette_trans[i]
}




#start reduce dims (norm expressin)
sg = sortGenes(counts(l_remove), class_info_merging_remove, binarizeMethod = "naive", cores = 16)
pp = getPValues(sg, numPerm = 20, cores = 16)
l_remove = scran::computeSumFactors(l_remove, sizes = seq(10, 200, 20), clusters = class_info_merging_remove, positive = TRUE)
l_remove = scater::normalize(l_remove)
l_remove = scater::calculateQCMetrics(l_remove)


#MNN is used here for visualization
kk = floor(sqrt(ncol(l_remove)) * 1) * ncol(l_remove)
newOrder = order(table(l_remove$fileID), decreasing = TRUE)
var_genes = names(which(apply(pp$adjpval, 1, function(x) any(x < 0.05))))
mnn = mnnCorrect(as.matrix(logcounts(l_remove[,which(l_remove$fileID == (unique(l_remove$fileID)[1]))])), as.matrix(logcounts(l_remove[,which(l_remove$fileID == (unique(l_remove$fileID)[2]))])), k = floor(sqrt(ncol(l_remove)) * 1), subset.row=var_genes, cos.norm.out=FALSE, order = newOrder, cos.norm.in=TRUE)
l_bc = cbind(mnn$corrected[[1]], mnn$corrected[[2]])
l_bc = SingleCellExperiment(assays=list(logcounts=l_bc))
l_bc$colData = l_remove$colData

#reduce dimensions of MNN-corrected expression
l_v = l_bc[which(rownames(l_bc) %in% names(which(apply(pp$adjpval, 1, function(x) any(x < 0.05))))),]
q = t(logcounts(l_v))
q = t(apply(q,1, function(x) x-mean(x)))
set.seed(111)
sv = irlba(q, 50)
d = sv$d
dw = 21 #based on the knee of singular values
dw = 1:dw
temp = t(apply(ss$u[,dw], 1, function(x) x / (sqrt(sum(x^2)))))
write.csv(temp, file = paste0(prefix, "dims.txt"), row.names = FALSE, col.names = FALSE) #file saved @https://raw.githubusercontent.com/mahmoudibrahim/KidneyMap/master/assets/reducedDims/SVD/mouse_PDGFRAB_dims.txt




#UMAP (produced in python, just for visualization, independent of clustering)
mapCoords = read.table("https://raw.githubusercontent.com/mahmoudibrahim/KidneyMap/master/assets/reducedDims/UMAP/mouse_PDGFRABpositive_umapCoords.csv", header = FALSE, sep=",")
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


#percent of cells per time point
pdf(paste0(prefix, "percentCell_heat.pdf"))
timeheat = t(apply(table(l_remove$level3,l_remove$kidney_function),1,function(x) x/sum(x))) * 100
level3_order = c(1,3,4,c(2,5,6,7,8,9,10)[c(1,4,5,2,6,3,7)]) #reorder the clusters manually
timeheat = timeheat[level3_order,]
heatmap.2(timeheat, trace = "none", scale = "none", col = colorRampPalette(c("white","dodgerblue4"))(n=100), cellnote=round(timeheat,2), notecol="black", Rowv = NA, Colv = NA, margins = c(15,15))
dev.off()




#some markers
markersp = toupper(c("Pdgfrb","Pdgfra","Col1a1","DCN","PTPRC","EPCAM","PECAM1","KRT18","PAX8","RNASE1","CXCR4","PLP1","GFRA3","SOX10","EMCN","EGFL7","hsd11b2","clcnka","atp6v1b1","MMRN1","Prox1","c7","col15a1","nkd2","scara5","ugt2b7","glyat","rergl","mustn1","ackr1","vwf","ramp2","frzb","tagln","krt8","adamts6","igf2","ssuh2","GPR183","CD69","srp1","mitf","fgf7","sall1","pax2","dnmt1","tp53","acta2","vim","ramp2","fn1","tmeff2","irf8","creb5","nrf1","atf3"))
names1 = toupper(unlist(lapply(strsplit(as.character(rownames(l_remove)), split = ";", fixed = T), function(x) as.character(x[2]))))
for (i in 1:length(markersp)) {
	if(any(names1 == markersp[i])) {
		aaax = as.vector(logcounts(l_remove)[which(names1 == markersp[i]),])
		zzz = pdf(paste0(prefix, "umap_", markersp[i], ".pdf"))
		plot(mapCoords, col = color.gradient(scale(aaax), colors = c("#ECECEC50",rev(c("#d7191c","#abdda4")))), pch = 20, cex = 1)
		legend.col(col = colorRampPalette(c("#ECECEC50",rev(c("#d7191c","#abdda4"))))(100), lev = seq(min(aaax), max(aaax), length.out=100))
		zzz = dev.off()
	}
}





#gene ontology
library(tibble)
library(clusterProfiler)
ww = c(which(l_remove$level3 %in% c("Fibroblast 2", "Myofibroblasts 2", "Myofibroblasts 4", "Myofibroblasts 1"))) #four clusters that match ATAC-Seq plots
l_g = l_remove[,ww]
names1 = toupper(unlist(lapply(strsplit(as.character(rownames(l_g)), split = ";", fixed = T), function(x) as.character(x[2]))))
rownames(l_g) = names1
sg_g  = sortGenes(counts(l_g), l_g$level3, binarizeMethod = "naive")
ng = plotTopMarkerHeat(sg_g, averageCells=10^2, top_n=100, plotheat=FALSE, outs = TRUE)
ng = ng[c(1,3,2,4)]

m_t2g =  as.tibble(read.gmt("https://raw.githubusercontent.com/mahmoudibrahim/KidneyMap/master/assets/public/c2.cp.pid.v7.0.symbols.gmt"))
en = mclapply(ng, function(x) enricher(x, TERM2GENE=m_t2g, minGSSize = 10, maxGSSize = 200, pvalueCutoff = 0.05), mc.cores = 12)
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
pdf(paste0(prefix, "PID_toMatchATAC.pdf"), width = 10, height = 10)
heatmap.2(t(en4), trace = "none", scale = "none", col = colorRampPalette(c("white","dodgerblue4"))(n=100), margins = c(20,25), density = "none", Colv = NA, key = FALSE, labCol = names(ng))
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


#violin plots (matrosome score)
pdf(paste0(prefix, "_core_matrisome_score_level1.pdf"))
size = scale(sqrt(table(l_remove$level1)),center = FALSE) * 0.8
vioplot(aaa ~ l_remove$level1, col = level1_col, ylim = c(0,2.2), plotCentre="line", wex = size, lwd = 3, border = level1_col)
dev.off()
pdf(paste0(prefix, "_core_matrisome_score_level2_mesenchymal.pdf"))
size = scale(sqrt(table(l_remove$level2[which(l_remove$level1 == "Mesenchymal")])),center = FALSE)
vioplot(aaa[which(l_remove$level1 == "Mesenchymal")] ~ l_remove$level2[which(l_remove$level1 == "Mesenchymal")], col = level1_col[4], ylim = c(0,2.2), plotCentre="line", wex = size, lwd = 3, border = level1_col[4], names = abbreviate(rownames(size)))
dev.off()
zzz = pdf(paste0(prefix, "umap_ecm_continuous.pdf"))
mapCoords = read.table("https://raw.githubusercontent.com/mahmoudibrahim/KidneyMap/master/assets/reducedDims/UMAP/mouse_PDGFRABpositive_umapCoords.csv", header = FALSE, sep=",")
plot(mapCoords, col = color.gradient(scale(aaa), colors = c("#ECECEC50",rev(c("#d7191c","#abdda4")))), pch = 20, cex = 1)
legend.col(col = colorRampPalette(c("#ECECEC50",rev(c("#d7191c","#abdda4"))))(100), lev = seq(min(aaa), max(aaa), length.out=100))
zzz = dev.off()



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

#violin plots (collagen score)
pdf(paste0(prefix, "_collagens_score_level1.pdf"))
size = scale(sqrt(table(l_remove$level1)),center = FALSE) * 0.8
vioplot(aaa ~ l_remove$level1, col = level1_col, ylim = c(0,4), plotCentre="line", wex = size, lwd = 3, border = level1_col)
dev.off()
pdf(paste0(prefix, "_collagens_score_level2_mesenchymal.pdf"))
size = scale(sqrt(table(l_remove$level2[which(l_remove$level1 == "Mesenchymal")])),center = FALSE)
vioplot(aaa[which(l_remove$level1 == "Mesenchymal")] ~ l_remove$level2[which(l_remove$level1 == "Mesenchymal")], col = level1_col[4], ylim = c(0,4), plotCentre="line", wex = size, lwd = 3, border = level1_col[4], names = abbreviate(rownames(size)))
dev.off()
zzz = pdf(paste0(prefix, "umap_collagens_continuous.pdf"))
mapCoords = read.table("https://raw.githubusercontent.com/mahmoudibrahim/KidneyMap/master/assets/reducedDims/UMAP/mouse_PDGFRABpositive_umapCoords.csv", header = FALSE, sep=",")
plot(mapCoords, col = color.gradient(scale(aaa), colors = c("#ECECEC50",rev(c("#d7191c","#abdda4")))), pch = 20, cex = 1)
legend.col(col = colorRampPalette(c("#ECECEC50",rev(c("#d7191c","#abdda4"))))(100), lev = seq(min(aaa), max(aaa), length.out=100))
zzz = dev.off()
###################################################




#marker heatmaps
names1 = tolower(unlist(lapply(strsplit(as.character(rownames(l_remove)), split = ";", fixed = T), function(x) as.character(x[2]))))
sg_heat = logcounts(l_remove)
rownames(sg_heat) = names1
sg_heat = sortGenes(sg_heat, sg_level3$inputClass, binarizeMethod = "naive", cores = 16)
level3_order = c(1,3,4,c(2,5,6,7,8,9,10)[c(1,4,5,2,6,3,7)])

pdf(paste0(prefix, "level3_top_auto_markers_10.pdf"), height = 15)
plotTopMarkerHeat(sg_heat, averageCells=10^2, newOrder=level3_order, gaps = FALSE, colors = colorRampPalette( c("cornflowerblue","black","gold"))(n=100), top_n=10)
dev.off()

pdf(paste0(prefix, "select_markers_wholeClusterAverage_expression.pdf"))
hend_markers =  tolower(c("Pdgfrb", "Pdgfra", "Itgb1","Vim", "Acta2", "Ly86","Ptprc","Cdh5","Pecam1","Cldn3", "epcam", "meg3", "Scara5", "Col1a1","Col15a1"))
plotMarkerHeat(sg_heat$inputMat, sg_heat$inputClass, hend_markers, clusterGenes = F, averageCells=10^6, newOrder=level3_order, colors = colorRampPalette(c("cornflowerblue","black","gold"))(n=100))
dev.off()






#diffusion map
ww = c(which(l_remove$level2 %in% c("Fibroblasts", "Myofibroblasts"))) #only mesenchymal cells
l_v = l_remove[,ww]
l_v = scran::computeSumFactors(l_v, sizes = seq(10, 200, 20), clusters = l_v$level3, positive = TRUE)
l_v = scater::normalize(l_v)
sg_mesen = sortGenes(counts(l_v), l_v$level3, binarizeMethod="naive", cores = 16)
pp_mesen = getPValues(sg_mesen, numPerm = 20, cores = 16)
marker_mesen = names(which(apply(pp_mesen$adjpval, 1, function(x) any(x < 0.05))))

#reduce dims after MNN
kk = floor(sqrt(ncol(l_v)) * 1) * ncol(l_v)
newOrder = order(table(l_v$fileID), decreasing = TRUE)
mnn = mnnCorrect(as.matrix(logcounts(l_v[,which(l_v$fileID == (unique(l_v$fileID)[1]))])), as.matrix(logcounts(l_v[,which(l_v$fileID == (unique(l_v$fileID)[2]))])), k = floor(sqrt(ncol(l_v)) * 1), subset.row=marker_mesen, cos.norm.out=FALSE, order = newOrder, cos.norm.in=TRUE)
l_bc = cbind(mnn$corrected[[1]], mnn$corrected[[2]])
l_bc = SingleCellExperiment(assays=list(logcounts=l_bc))
colData(l_bc)= colData(l_v)

l_v = l_bc[which(rownames(l_bc) %in% names(which(apply(pp$adjpval, 1, function(x) any(x < 0.05))))),]
q = t(logcounts(l_v[which(rownames(l_v) %in% marker_mesen),]))
q = t(apply(q,1, function(x) x-mean(x)))
svs = svd(q)
d = sv$d
xx = 1:1000
yy = d[1:1000]
p1 = c(xx[1],yy[1])
p2 = c(xx[length(xx)], yy[length(yy)])
dw = which.max(sapply(1:length(yy), function(x) dist2d(c(xx[x], yy[x]), p1, p2)))
dw = 1:dw
temp = t(apply(ss$u[,dw], 1, function(x) x / (sqrt(sum(x^2)))))
write.csv(temp, file = paste0(prefix, "dims_pseudotime.txt"), row.names = FALSE, col.names = FALSE)


matrisome_set = read.table("/home/mibrahim/Dropbox/Kidney_Interstitium_collab_Henderson/Public_data/Matrisome/ecm_genes_human.txt", sep = "\t", header = TRUE)
collagens = matrisome_set$Gene.Symbol[which(matrisome_set$Category %in% c("Collagens"))]
reads_single_phase = logcounts(l_remove)
rownames(reads_single_phase) = toupper(unlist(lapply(strsplit(as.character(rownames(reads_single_phase)), split = ";", fixed = T), function(x) as.character(x[2]))))
reads_single_phase = as.matrix(reads_single_phase[rownames(reads_single_phase) %in% toupper(collagens) ,])
combined_matrix = rbind(reads_single_phase,average=apply(reads_single_phase,2,mean))
cor_matrix = cor(t(combined_matrix))
cor_vector = cor_matrix[,dim(cor_matrix)[1]]
reads_single_phase_restricted = reads_single_phase[rownames(reads_single_phase) %in% names(cor_vector[cor_vector >= 0.1]),]
aaa = apply(reads_single_phase_restricted,2,mean)

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

########diffusion map based
library(destiny) 
set.seed(111)
df = DiffusionMap(ss$u[,dw], k = 40, n_eigs = 3) 
write.csv(df@eigenvectors, file = paste0(prefix, "dims_diffusionmap.txt"), row.names = FALSE, col.names = FALSE) #file now saved @https://raw.githubusercontent.com/mahmoudibrahim/KidneyMap/master/assets/reducedDims/DM/mouse_PDGFRBpositive_dm.txt


#plot 3d
library(plot3D)
df2 = as.matrix(read.csv("https://raw.githubusercontent.com/mahmoudibrahim/KidneyMap/master/assets/reducedDims/DM/mouse_PDGFRBpositive_dm.txt"))
pdf(paste0(prefix, "pseudotime_specific_dm_level3.pdf"))
scatter3D(-df2[,1], df2[,2], df2[,3], col = colors_info_mesenl3, pch = 20, cex = 0.5, colvar = NULL, bty = "n", phi = 00, theta = 310)
legend("topright", legend = sort(unique(l_v$level3)), pch=21, col = col_palette_trans, pt.bg = col_palette, pt.cex=2, cex=.8, bty="n", ncol = ceiling(length(unique(l_v$level3)) / 20))
dev.off()
pdf(paste0(prefix, "pseudotime_specific_dm_level2.pdf"))
scatter3D(-df2[,1], df2[,2], df2[,3], col = colors_info_mesenl2, pch = 20, cex = 0.5, colvar = NULL, bty = "n", phi = 00, theta = 310)
legend("topright", legend = sort(unique(l_v$level3)), pch=21, col = col_palette_trans, pt.bg = col_palette, pt.cex=2, cex=.8, bty="n", ncol = ceiling(length(unique(l_v$level2)) / 20))
dev.off()
pdf(paste0(prefix, "pseudotime_specific_dm_ecm.pdf"))
scatter3D(-df2[,1], df2[,2], df2[,3], col = color.gradient(scale(aaa[ww]), colors = c("#ECECEC50",rev(c("#d7191c","#abdda4")))), pch = 19, cex = 1, colvar = NULL, bty = "n", phi = 00, theta = 310)
dev.off()




#for ISH
threeGene = t(sg_level3$binary[which(rownames(sg_level3$binary) %in% c("ENSMUSG00000029231.15;Pdgfra", "ENSMUSG00000021567.15;Nkd2", "ENSMUSG00000001506.10;Col1a1")),])

bin = threeGene
bin[bin > 0] = 1
bin = cbind(bin, rep(1, nrow(bin)))
bin = bin[,c(1,3,2,4)]
colnames(bin) = c("PDGFRa","Nkd2","Col1a1","Other Genes")
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
colnames(bin2) = c("Nkd2_and_PDGFRa","Nkd2_Only","PDGFRa_Only","Other_Genes","Col1a1")

pdf(paste0(prefix, "percent_heat_withDAPI.pdf"))
agg = aggregate(bin2[,1:4], by = list(bin2[,5]), sum)
agg_per = apply(agg[,2:5], 1, function(x) x/sum(x)) * 100
heatmap.2(agg_per, trace = "none", col = colorRampPalette(c("white","dodgerblue4"))(n=100), scale = "none", cellnote = round(agg_per, 2), notecol = "black", notecex = 1.5, dendrogram = "none", Rowv = NA, Colv = NA, breaks = seq(0,100,length.out=101),key = FALSE, margin = c(10,20))
dev.off()
bb = rep("",nrow(bin2))
bb[which(bin2[,1] == 1)] = "d_Both"
bb[which(bin2[,2] == 1)] = "c_n_only"
bb[which(bin2[,3] == 1)] = "b_a_only"
bb[which(bin2[,4] == 1)] = "a_dapi_only"
threeGene = t(logcounts(l_remove)[which(rownames(sg_level3$binary) %in% c("ENSMUSG00000029231.15;Pdgfra", "ENSMUSG00000021567.15;Nkd2", "ENSMUSG00000001506.10;Col1a1")),])
pdf(paste0(prefix, "Collagen_by_cell.pdf"))
size = scale(sqrt(table(bb)),center = FALSE)
vioplot(as.vector(threeGene[,2]) ~ bb, col = col_palette_short[c(1,3,2,4)], plotCentre="line", wex = size, lwd = 3, border = col_palette_short[c(1,3,2,4)], names = c("other_genes","a_only","n_only","both"), ylab = "Log Collagen Expression", main = "All PDGFRa+/b+ Cells (n= 7245)")
dev.off()






#heatmap plot to compare with ATAC-seq
ww = c(which(l_remove$level3 %in% c("Fibroblast 2", "Myofibroblasts 1", "Myofibroblasts 2", "Myofibroblasts 4")))
l_v = l_remove[,ww]
sg_tf = sortGenes(counts(l_v), l_v$level3, binarizeMethod = "naive")
genes = c("ENSMUSG00000052684.4;Jun", "ENSMUSG00000071076.6;Jund", "ENSMUSG00000052837.6;Junb", "ENSMUSG00000003545.3;Fosb",  "ENSMUSG00000021250.13;Fos", "ENSMUSG00000026628.13;Atf3", "ENSMUSG00000037174.18;Elf2", "ENSMUSG00000036461.16;Elf1", "ENSMUSG00000036602.14;Alx1", "ENSMUSG00000032035.15;Ets1", "ENSMUSG00000022895.16;Ets2", "ENSMUSG00000019947.10;Arid5b", "ENSMUSG00000053007.9;Creb5", "ENSMUSG00000030551.14;Nr2f2", "ENSMUSG00000042406.8;Atf4", "ENSMUSG00000023034.7;Nr4a1",  "ENSMUSG00000034168.7;Irf2bpl", "ENSMUSG00000041515.10;Irf8",  "ENSMUSG00000009079.16;Ewsr1", "ENSMUSG00000055148.7;Klf2", "ENSMUSG00000033863.1;Klf9", "ENSMUSG00000058440.14;Nrf1")
pdf(paste0(prefix, "TFs_expression.pdf"))
plotMarkerHeat(sg_tf$inputMat, sg_tf$inputClass, genes, averageCells=10^4, colors = colorRampPalette(c("cornflowerblue","black","gold"))(n=100), newOrder=c(1,3,2,4), clusterGenes=T, clusterGenesK=4, gaps = FALSE)
dev.off()







##Fin
