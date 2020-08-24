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
pp = getPValues(sg, numPerm = 20, cores = 16)
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
rem = which(apply(pp$adjpval,2,function(x) length(x[x<0.05])) == 0)
removethis = names(rem)
removethis = which(class_info_merging %in% removethis)


#find pdgfrb low clusters
pdgfrb_percent = sg$condGeneProb[which(rownames(sg$condGeneProb) == "ENSMUSG00000024620.11;Pdgfrb"),-rem] * 100
pdfgrb_cut = median(pdgfrb_percent) - (median(abs(pdgfrb_percent - median(pdgfrb_percent))) * 1)
pdfgrb_cut = names(which(pdgfrb_percent < pdfgrb_cut))
pdgfra_percent = sg$condGeneProb[which(rownames(sg$condGeneProb) == "ENSMUSG00000029231.15;Pdgfra"),-rem] * 100
pdfgra_cut = median(pdgfra_percent) - (median(abs(pdgfra_percent - median(pdgfra_percent))) * 1)
pdfgra_cut = names(which(pdgfra_percent < pdfgra_cut))
pdg_cut = intersect(pdfgrb_cut, pdfgra_cut) #low on both PDGFRa and PDGFRb
pdg_cut = which(class_info_merging %in% pdg_cut)
removethis = sort(unique(c(removethis, pdg_cut)))



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
var_genes = names(which(apply(pp$adjpval, 1, function(x) any(x < 0.05))))
mnn = batchelor::fastMNN(l_remove, batch = l_remove$fileID, d = 30, auto.order = TRUE, subset.row = var_genes, BSPARAM = BiocSingular::IrlbaParam(tol=1e-8))
mnn_u = mnn@reducedDims$corrected
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
pdf(paste0(prefix, "superUMAP_level1_level3Colors.pdf"))
plot(mapCoords[,1], mapCoords[,2], col = colors_info_level3, pch = 20, xlab = "UMAP-1", ylab = "UMAP-2")
legend(min(mapCoords[,1]), min(mapCoords[,2]) + 1, legend = abbreviate(sort(unique(l_remove$level3))), pch=21, col=col_palette, pt.bg=col_palette_trans, pt.cex=2, cex=.8, bty="n", ncol = ceiling(length(unique(l_remove$level3)) / 2))
dev.off()
pdf(paste0(prefix, "UMAP_KF.pdf"))
plot(mapCoords[,1], mapCoords[,2], pch = 20, col = colors_info_kf, cex = 0.3, xlim = c(-7,15))
dev.off()
pdf(paste0(prefix, "UMAP_fileID.pdf"))
plot(mapCoords[,1], mapCoords[,2], pch = 20, col = colors_file, cex = 0.3)
legend(min(mapCoords[,1]) + 1, min(mapCoords[,2]) + 15, legend = abbreviate(unique(l_remove$fileID), minlength=6), pch=21, col = col_palette, pt.bg = col_palette_trans, pt.cex=2, cex=.8, bty="n", ncol = ceiling(length(unique(l_remove$fileID)) / 10))
dev.off()


#percent of cells per time point
library(gplots)
pdf(paste0(prefix, "percentCell_heat2.pdf"))
timeheat = (apply(table(l_remove$level3,l_remove$kidney_function),2,function(x) x/sum(x))) * 100
timeheat = timeheat[level3_order = c(1,3,4,10,2,6,7,8,9,5),]
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



ww = c(which(l_remove$level2 %in% c("Fibroblasts", "Myofibroblasts")))
l_v = l_remove[,ww]
l_v = scran::computeSumFactors(l_v, sizes = seq(10, 200, 20), clusters = l_v$level3, positive = TRUE)
l_v = scater::normalize(l_v)
names1 = unlist(lapply(strsplit(as.character(rownames(l_v)), split = ";", fixed = T), function(x) as.character(x[2])))
rownames(l_v) = names1
l_v = l_v[!duplicated(rownames(l_v)),]


colors_info_mesenl3 = rep("", length = ncol(l_v))
for (i in 1:length(unique(l_v$level3))) {
	wtemp = which(l_v$level3 == (sort(unique(l_v$level3)))[i])
	colors_info_mesenl3[wtemp] = col_palette_trans[i]
}



#microarray
load("GSE121190_exprs.RData")
pa = rowMeans(pdgfra_fib[,1:3]) - rowMeans(pdgfra_fib[,4:6])
inter = intersect(names(pa), rownames(l_v))
pa_int = pa[which(names(pa) %in% inter)]
pa_int = pa_int[order(names(pa_int))]
exp_int = logcounts(l_v)[which(rownames(l_v) %in% inter),]
exp_int = exp_int[order(rownames(exp_int)),]
exp_int_fc = exp_int - rowMeans(exp_int)
pdf(paste0(prefix, "microarray_box_pearson.pdf"))
aaa = apply(as.matrix(exp_int_fc), 2, function(x) cor(x, pa_int, method = "pearson"))
facts = factor(as.factor(l_v$level3),levels(as.factor(l_v$level3))[c(1,3,4,5,6,2)])
plot(as.numeric(facts)+rnorm(length(aaa), 0, 0.13), aaa, pch = 20, xlab = "", ylab = "", cex = 0.3, main = "Correlation to UUO PDGFRa+ Fibroblasts", col = colors_info_mesenl3, xaxt = "n")
axis(1, at=sort(unique(as.numeric(facts))), labels=abbreviate(levels(facts)))
vioplot(as.vector(aaa) ~ facts, col = makeTransparent("white", alpha = 0), border = makeTransparent("black", alpha = 0.3), range = 0.1, outline = F, add = T, notch = T, lwd = 4,  plotCentre="line", xlab = "Cell Clusters", ylab = "Correlation to UUO PDGFRa+ Fibroblasts")
dev.off()




#PID
library(clusterProfiler)
library(gplots)
library(msigdbr)
l_g = l_remove
names1 = unlist(lapply(strsplit(as.character(rownames(l_g)), split = ";", fixed = T), function(x) as.character(x[2])))
rownames(l_g) = names1
sg_g  = sortGenes(counts(l_g), l_g$level3, binarizeMethod = "naive")
ng = plotTopMarkerHeat(sg_g, averageCells=10^2, top_n=200, plotheat=FALSE, outs = TRUE)
ng = ng[c(2,6,7,8,9,5)]
msig = msigdbr::msigdbr(species = "Mus musculus", category = "C2")
msig = data.frame(msig$gs_name[which(msig$gs_subcat == "CP:PID")], msig$gene_symbol[which(msig$gs_subcat == "CP:PID")])
en = mclapply(ng, function(x) enricher(x, TERM2GENE=msig, minGSSize = 10, maxGSSize = 200, pvalueCutoff = 0.05), mc.cores = 12)
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
pdf(paste0(prefix, "PID.pdf"), width = 10, height = 10)
heatmap.2(t(en4), trace = "none", scale = "none", col = colorRampPalette(c("white","dodgerblue4"))(n=100), margins = c(20,25), density = "none", key = FALSE, labCol = names(ng))
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





#marker heatmaps
names1 = tolower(unlist(lapply(strsplit(as.character(rownames(l_remove)), split = ";", fixed = T), function(x) as.character(x[2]))))
sg_heat = logcounts(l_remove)
rownames(sg_heat) = names1
sg_heat = sortGenes(sg_heat, sg_level3$inputClass, binarizeMethod = "naive", cores = 16)
level3_order = c(1,3,4,10,2,7,8,9,5,6)
pdf(paste0(prefix, "level3_top_auto_markers_10.pdf"), height = 15)
plotTopMarkerHeat(sg_heat, averageCells=10^2, newOrder=level3_order, gaps = FALSE, colors = colorRampPalette( c("cornflowerblue","black","gold"))(n=100), top_n=10)
dev.off()
pdf(paste0(prefix, "select_markers_wholeClusterAverage.pdf"))
hend_markers =  tolower(c("Pdgfrb", "Pdgfra", "Itgb1","Vim", "Acta2", "Ly86","Ptprc","Cdh5","Pecam1","Cldn3", "epcam", "bmp5",  "meg3", "Scara5", "Col1a1","Col15a1"))
plotBinaryHeat(sg_heat$binary, sg_heat$inputClass, hend_markers, clusterGenes = F, averageCells=10^6, newOrder=level3_order, colors = colorRampPalette(c("black","gold"))(n=100))
dev.off()




#Nkd2-centric
ww = c(which(l_remove$level2 %in% c("Fibroblasts", "Myofibroblasts"))) #only mesenchymal cells
l_v = l_remove[,ww]
l_v = scran::computeSumFactors(l_v, sizes = seq(10, 200, 20), clusters = l_v$level3, positive = TRUE)
l_v = scater::normalize(l_v)
ll = logcounts(l_v)
ll = apply(ll, 1, function(x) x/(sqrt(sum(x^2))))
ll = t(ll)
nk = as.vector(ll[which(rownames(ll) == "Nkd2"),])
corr = apply(ll, 1, function(x) cor(x, nk))

top = head(names(sort(corr, decreasing = TRUE)), n = 100)
bottom = head(names(sort(corr, decreasing = FALSE)), n = 100)
ng = list("correlated" = top, "anti-correlated" = bottom)
msig = msigdbr::msigdbr(species = "Mus musculus", category = "C2")
msig = data.frame(msig$gs_name[which(msig$gs_subcat == "CP:PID")], msig$gene_symbol[which(msig$gs_subcat == "CP:PID")])
en = mclapply(ng, function(x) enricher(x, TERM2GENE=msig, minGSSize = 10, maxGSSize = 200, pvalueCutoff = 0.05), mc.cores = 12)
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
msig = msigdbr::msigdbr(species = "Mus musculus", category = "C2")
msig = data.frame(msig$gs_name[which(msig$gs_subcat == "CP:KEGG")], msig$gene_symbol[which(msig$gs_subcat == "CP:KEGG")])
cp = clusterProfiler::compareCluster(geneCluster = ng, fun = "enricher", TERM2GENE = msig)
en = mclapply(ng, function(x) enricher(x, TERM2GENE=msig, minGSSize = 10, maxGSSize = 200, pvalueCutoff = 0.05), mc.cores = 12)
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
pdf(paste0(prefix, "PID_KEGG.pdf"), width = 10, height = 10)
heatmap.2(t(en4), trace = "none", scale = "none", col = colorRampPalette(c("white","dodgerblue4"))(n=100), margins = c(20,25), density = "none", key = FALSE, labCol = names(ng), cexRow = 0.8, cexCol = 0.8, Colv = NA)
dev.off()





#fibro lineage
ww = c(which(l_remove$level2 %in% c("Fibroblasts", "Myofibroblasts"))) #only mesenchymal cells
l_v = l_remove[,ww]
l_v = scran::computeSumFactors(l_v, sizes = seq(10, 200, 20), clusters = l_v$level3, positive = TRUE)
l_v = scater::normalize(l_v)
sg_mesen = sortGenes(counts(l_v), l_v$level3, binarizeMethod="naive", cores = 16)
pp_mesen = getPValues(sg_mesen, numPerm = 20, cores = 16)
marker_mesen = names(which(apply(pp_mesen$adjpval, 1, function(x) any(x < 0.05))))
mnn = batchelor::fastMNN(l_v, batch = l_v$fileID, d = 30, auto.order = TRUE, subset.row = marker_mesen, BSPARAM = BiocSingular::IrlbaParam(tol=1e-8))
mnn_u = mnn@reducedDims$corrected
temp = t(apply(mnn_u, 1, function(x) x / (sqrt(sum(x^2)))))
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
ecm = aaa[ww]


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

#diffusion map based
library(destiny) 
set.seed(111)
df = DiffusionMap(ss$u[,dw], k = 40, n_eigs = 3) 
write.csv(df@eigenvectors, file = paste0(prefix, "dims_diffusionmap.txt"), row.names = FALSE, col.names = FALSE) #file now saved @https://raw.githubusercontent.com/mahmoudibrahim/KidneyMap/master/assets/reducedDims/DM/mouse_PDGFRBpositive_dm.txt
library(plot3D)
pdf(paste0(prefix, "pseudotime_specific_dm_level3.pdf"))
scatter3D(df2[,1],-df2[,2], df2[,3], col = colors_info_mesenl3, pch = 20, cex = 0.5, colvar = NULL, bty = "n", phi = 60, theta = 210)
legend("topright", legend = sort(unique(l_v$level3)), pch=21, col = col_palette_trans, pt.bg = col_palette, pt.cex=2, cex=.8, bty="n", ncol = ceiling(length(unique(l_v$level3)) / 20))
dev.off()
pdf(paste0(prefix, "pseudotime_specific_dm_level2.pdf"))
scatter3D(df2[,1],-df2[,2], df2[,3], col = colors_info_mesenl3, pch = 20, cex = 0.5, colvar = NULL, bty = "n", phi = 60, theta = 210)
legend("topright", legend = sort(unique(l_v$level3)), pch=21, col = col_palette_trans, pt.bg = col_palette, pt.cex=2, cex=.8, bty="n", ncol = ceiling(length(unique(l_v$level2)) / 20))
dev.off()
pdf(paste0(prefix, "pseudotime_specific_dm_ecm.pdf"))
scatter3D(df2[,1], -df2[,2], df2[,3], col = color.gradient(scale(aaa[ww]), colors = c("#ECECEC50",rev(c("#d7191c","#abdda4")))), pch = 20, cex = 0.8, colvar = NULL, bty = "n", phi = 60, theta = 210)
dev.off()



#umap based
mapCoords = read.table("https://raw.githubusercontent.com/mahmoudibrahim/KidneyMap/master/assets/reducedDims/UMAP/mouse_PDGFRABpositive_mesenchymal_umapCoords.csv", header = FALSE, sep=",")
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
library(slingshot)
sl1 = slingshot(mapCoords, l_v$level3, start.clus = c("Myofibroblasts 1a", "Myofibroblasts 1b"))
st1 = rowMeans(slingPseudotime(sl1), na.rm = TRUE)
pdf(paste0(prefix, "pseudotime_specific_umap_pseudotime.pdf"))
col_pst = color.gradient(scale(-st1), colors = c("#ECECEC50",rev(c("#d7191c","#abdda4"))))
plot(mapCoords, col = col_pst, pch = 20)
dev.off()
pdf(paste0(prefix, "pseudotime_umap_level3_withLineageTree.pdf"))
plot(mapCoords, col = colors_info_mesenl3, pch = 20)
legend(min(mapCoords[,1]), min(mapCoords[,2]) + 6, legend = sort(unique(l_v$level3)), pch=21, col = col_palette_trans, pt.bg = col_palette, pt.cex=2, cex=.8, bty="n", ncol = ceiling(length(unique(l_v$level3)) / 20))
lines(sl1, lwd=5, col="#00000095", type = "lineages")
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
genes = toupper(c("Atf3","Atf4","Runx1","Nr4a","Jun","Fos","Mef2a","Mef2c","Mef2d","Nrf4a1","Nrfa2","Meis1","Sox4","Irf8","Sox4","Gata6","Arid5a","Alx1","Creb5","Dbp","Nrf1","Yy1","Atoh8","Klf2","klf9","elf1","elf2","ets1","ets2","smad1","smad3","smad5","smad7","Ewsr1","max","myc","maf1","mafb","maff","egr1","egr2"))
ww = c(which(l_remove$level3 %in% c("Fibroblasts 2", "Myofibroblasts 1a", "Myofibroblasts 2", "Myofibroblasts 1b","Myofibroblasts 3a","Myofibroblasts 3b")))
l_v = l_remove[,ww]
names1 = toupper(unlist(lapply(strsplit(as.character(rownames(l_v)), split = ";", fixed = T), function(x) as.character(x[2]))))
rownames(l_v) = names1
sg_tf = sortGenes(logcounts(l_v), l_v$level3, binarizeMethod = "naive")
 pdf(paste0(prefix, "TFs_expression.pdf"))
plotMarkerHeat(sg_tf$inputMat, sg_tf$inputClass, genes, averageCells=10^2, colors = colorRampPalette(c("cornflowerblue","black","gold"))(n=100), clusterGenes=T, clusterGenesK=7, gaps = TRUE, newOrder= c(1,4,5,6,2,3))
dev.off()




#saveoutput
writeMM(counts(l_remove), file = paste0(prefix, "UMI_counts.mtx"))
names1 = toupper(unlist(lapply(strsplit(as.character(rownames(l_remove)), split = ";", fixed = T), function(x) as.character(x[2]))))
names2 = toupper(unlist(lapply(strsplit(as.character(rownames(l_remove)), split = ";", fixed = T), function(x) as.character(x[1]))))
names2 = toupper(unlist(lapply(strsplit(as.character(names2), split = ".", fixed = T), function(x) as.character(x[1]))))
rowdat = cbind(names1, names2)
colnames(rowdat) = c("Gene.Symbol","ENSEMBL.ID")
write.table(rowdat, file = paste0(prefix, "UMI_counts_rowData.txt"), sep = "\t", row.names=FALSE, quote = FALSE)
coldat = cbind(l_remove$level1, l_remove$level2, l_remove$level3, l_remove$kidney_function, l_remove$fileID)
colnames(coldat) = c("Annotation.Level.1","Annotation.Level.2","Annotation.Level.3","Kidney.Function","File ID")
write.table(coldat, file = paste0(prefix, "UMI_counts_colData.txt"), sep = "\t", row.names=FALSE, quote = FALSE)




##Fin
