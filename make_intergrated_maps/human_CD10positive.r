########################################################################
# human_CD10positive.r
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




prefix = "human_CD10positive_"



###########################################################################
########################START ANALYSIS#####################################
###########################################################################

#load data
load("human_CD10positive.RData") #R object containing aggregated data, saved as SingleCellExperiment. Variable name: l
lx = l #save full matrix for later


#gene filtering
cutGenes = log10(floor(ncol(l) * 0.001))
choose_genes = which((log10(Matrix::rowSums(counts(l)))) >= cutGenes)
l = l[choose_genes, ]




#first stuff [merging similar clusters]
sg = sortGenes(counts(l), colData(l)$class_info, binarizeMethod="naive", cores = 16)
write.files(sg, prefix = paste0(prefix, "genesorteRouts_noMerging_noRemove"))
marks = plotTopMarkerHeat(sg, top_n=100, plotheat=FALSE, outs = TRUE)
pdf(paste0(prefix, "correlationBetweenClusters_noMerging_noRemove.pdf"), height = 12, width = 12)
corr = plotCorrelationHeat(sg, displayNumbers=F, markers=getMarkers(sg)$markers, outs = TRUE)
dev.off()


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
removethis = c(names(which(apply(pp$adjpval,2,function(x) length(x[x<0.05])) == 0)), "en9") #en9 based on top marker genes (ribo proteins)
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
barplot((table(l_remove$kidney_function) / ncol(l_remove)) * 100, col = rev(kf_col), ylim = c(0,80), yaxp = c(0,80,4), ylab = "% of Cells")
dev.off()

KFlow = rep(0,length(unique(class_info_merging_remove)))
for (i in 1:length(unique(class_info_merging_remove))) {
	where = which(class_info_merging_remove == unique(class_info_merging_remove)[i])
	tab = table(l_remove$kidney_function[where]) / length(where)
	if (length(tab) == 2) {
		KFlow[i] = tab[1]
	} else { #edit this otherwise
		KFlow[i] = tab
	}
}
names(KFlow) = unique(class_info_merging_remove)



#load cluster annotation info
annot = read.table("https://raw.githubusercontent.com/mahmoudibrahim/KidneyMap/master/assets/clusterInfo/human_CD10positive.txt", sep = "\t")
annot = annot[match(names(KFlow), annot$V1),]
annot$V7 = KFlow
annot$V8 = log(KFlow) - log((table(l_remove$kidney_function) / ncol(l_remove))[2])
nums = table(class_info_merging_remove)
annot = annot[match(names(nums), annot$V1),]
annot$V9 = as.numeric(nums)
annot$V10 = as.numeric(log10(nums))
colnames(annot) = c("ID", "Label", "Annotation 1", "Annotation 2", "Annotation 3", "Percent Low KF Cells", "Log Enrichment of Low KF Cells", "Number of Cells", "log10 Number of Cells")
write.csv(annot, paste0(prefix, "_cluster_annotation_info_withAddData.csv"), row.names = FALSE) ##write cluster annotation information
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
sg_level3 = sortGenes(counts(l_remove), l_remove$level3, binarizeMethod = "naive", cores = 16)
write.files(sg_level3, prefix = paste0(prefix, "genesorteRouts_level3"))




#Kidney Function
temp = counts(l_remove)
names1 = toupper(unlist(lapply(strsplit(as.character(rownames(l_remove)), split = ";", fixed = T), function(x) as.character(x[2]))))
rownames(temp) = names1
rRNA = grep("^MRP|^RPL|^RPS", rownames(temp), value=T) #remove ribsomoal proteins and pseudogenes for this comparison
temp = temp[-(which(rownames(temp) %in% rRNA)),]
sg_kf = sortGenes(temp, l_remove$kidney_function, binarizeMethod = "naive")
pp_kf = getPValues(sg_kf, numPerm = 20, cores = 16)
ng = plotTopMarkerHeat(sg_kf, averageCells=10^2, top_n=100, plotheat=FALSE, outs = TRUE)
write.files(sg_kf, prefix = paste0(prefix, "Kidney_Function")) #this is used for PantherDB
#KEGG pathway 
library(tibble)
library(clusterProfiler)
m_t2g =  as.tibble(read.gmt("https://raw.githubusercontent.com/mahmoudibrahim/KidneyMap/master/assets/public/c2.cp.kegg.v7.0.symbols.gmt"))
en = mclapply(ng, function(x) enricher(x, TERM2GENE=m_t2g, minGSSize = 10, maxGSSize = 200, pvalueCutoff = 0.05), mc.cores = 12)
en2 = list()
for (i in 1:length(en)) {
	en2[[i]] =  en[[i]]@result$ID[1:3]
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
pdf(paste0(prefix, "kegg_KidneyFunction.pdf"), width = 10, height = 10)
heatmap.2(t(en4), trace = "none", scale = "none", col = colorRampPalette(c("white","dodgerblue4"))(n=100), margins = c(20,30), density = "none", Colv = NA, key = FALSE, labCol = names(ng), labcex = 0.9)
dev.off()
pdf(paste0(prefix, "kegg_KidneyFunction_bar.pdf"), width = 10, height = 10)
barplot(matrix(as.numeric(t(en3[,c(2,3,4)])),nrow = 2, byrow = TRUE), beside=TRUE, horiz=TRUE, col = rev(c("#99999995","#fc8d6295")), xlim = c(0,10), xlab = "-log10(q-value)")
dev.off()



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



#start reducing dimensions (norm. expression)
sg = sortGenes(counts(l_remove), class_info_merging_remove, binarizeMethod = "naive", cores = 16)
pp = getPValues(sg, numPerm = 20, cores = 16)
l_remove = scran::computeSumFactors(l_remove, sizes = seq(10, 200, 20), clusters = class_info_merging_remove, positive = TRUE)
l_remove = scater::normalize(l_remove)
l_remove = scater::calculateQCMetrics(l_remove)
a = plotExprsFreqVsMean(l_remove)

#MNN is used here for visualization
kk = floor(sqrt(ncol(l_remove)) * 1) * ncol(l_remove)
newOrder = order(table(l_remove$fileID), decreasing = TRUE)
var_genes = names(which(apply(pp$adjpval, 1, function(x) any(x < 0.05))))
mnn = mnnCorrect(as.matrix(logcounts(l_remove[,which(l_remove$fileID == (unique(l_remove$fileID)[1]))])), as.matrix(logcounts(l_remove[,which(l_remove$fileID == (unique(l_remove$fileID)[2]))])), as.matrix(logcounts(l_remove[,which(l_remove$fileID == (unique(l_remove$fileID)[3]))])), as.matrix(logcounts(l_remove[,which(l_remove$fileID == (unique(l_remove$fileID)[4]))])), as.matrix(logcounts(l_remove[,which(l_remove$fileID == (unique(l_remove$fileID)[5]))])), as.matrix(logcounts(l_remove[,which(l_remove$fileID == (unique(l_remove$fileID)[6]))])), k = 139, subset.row=var_genes, cos.norm.out=FALSE, order = newOrder, cos.norm.in=TRUE)
l_bc = cbind(mnn$corrected[[1]], mnn$corrected[[2]], mnn$corrected[[3]], mnn$corrected[[4]], mnn$corrected[[5]], mnn$corrected[[6]])
l_bc = SingleCellExperiment(assays=list(logcounts=l_bc))
l_bc$colData = l_remove$colData
	

#reduce dimensions of MNN-corrected expression
l_v = l_bc[which(rownames(l_bc) %in% names(which(apply(pp$adjpval, 1, function(x) any(x < 0.05))))),]
q = t(logcounts(l_v))
q = t(apply(q,1, function(x) x-mean(x)))
set.seed(111)
sv = irlba(q, 50)
d = sv$d
dw = 13 #based on the knee of the singular values
dw = 1:dw
temp = t(apply(ss$u[,dw], 1, function(x) x / (sqrt(sum(x^2)))))
write.csv(temp, file = paste0(prefix, "dims.txt"), row.names = FALSE, col.names = FALSE) #file saved @https://raw.githubusercontent.com/mahmoudibrahim/KidneyMap/master/assets/reducedDims/SVD/human_CD10positive_dims.txt



#UMAP (produced in python, just for visualization, independent of clustering)
mapCoords = read.table("https://raw.githubusercontent.com/mahmoudibrahim/KidneyMap/master/assets/reducedDims/UMAP/human_CD10positive_umapCoords.csv", header = FALSE, sep=",")
pdf(paste0(prefix, "regularUMAP_level1_level3TEXT.pdf"))
plot(mapCoords[,1], mapCoords[,2], pch = 20, col = colors_info_level1, cex = 0.2)
allClasses = unique(l_remove$label)
for (i in 1:length(allClasses)) {
	ww = which(l_remove$label == allClasses[i])
	xx = median(mapCoords[ww,1])
	yy = median(mapCoords[ww,2])
	text(xx,yy,allClasses[i], cex = 0.7, adj = 0.5)
}
dev.off()
dev.off()
pdf(paste0(prefix, "regularUMAP_KF.pdf"))
plot(mapCoords[,1], mapCoords[,2], pch = 20, col = colors_info_kf, cex = 0.4)
dev.off()
pdf(paste0(prefix, "regularUMAP_patientID.pdf"))
plot(mapCoords[,1], mapCoords[,2], pch = 20, col = colors_patient, cex = 0.3)
legend(max(mapCoords[,1]) - 1.5, min(mapCoords[,2]) + 2, legend = abbreviate(unique(l_remove$patientID), minlength=6), pch=21, col = col_palette, pt.bg = col_palette_trans, pt.cex=2, cex=.8, bty="n", ncol = ceiling(length(unique(l_remove$patientID)) / 10))
dev.off()
pdf(paste0(prefix, "clusterPercent_versus_KF.pdf"), height = 15)
howmany = table(l_remove$level3, l_remove$kidney_function)
howmany = t(apply(howmany, 1, function(x) x/sum(x)))
heatmap.2(howmany, trace = "none", col = colorRampPalette(c("white","#0072b2"))(n=100), cellnote=round(howmany*100,2), notecol="black", margins = c(15,15), key = FALSE, dendrogram = "none")
dev.off()




#general marker heatmaps
level3_order = c(1,2,3,4,5,6)
sg_temp = sortGenes(logcounts(l_remove), l_remove$level3, binarizeMethod = "naive", cores = 16)
tiff(paste0(prefix,"level3_top_auto_markers_5.tiff"), res=300, height = 25, width = 12, units = "in")
plotTopMarkerHeat(sg_temp, averageCells=10^5, newOrder=level3_order, gaps = FALSE, colors = colorRampPalette( c("cornflowerblue","black","gold"))(n=100), top_n=5)
dev.off()
tiff(paste0(prefix,"level3_top_auto_markers_10.tiff"), res=300, height = 25, width = 12, units = "in")
plotTopMarkerHeat(sg_temp, averageCells=10^5, newOrder=level3_order, gaps = FALSE, colors = colorRampPalette( c("cornflowerblue","black","gold"))(n=100), top_n=10)
dev.off()
tiff(paste0(prefix,"level3_top_auto_markers_20.tiff"), res=300, height = 25, width = 12, units = "in")
plotTopMarkerHeat(sg_temp, averageCells=10^5, newOrder=level3_order, gaps = FALSE, colors = colorRampPalette( c("cornflowerblue","black","gold"))(n=100), top_n=20)
dev.off()




#core matrisome score
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
zzz = pdf(paste0(prefix, "umap_ecm_continuous.pdf"))
mapCoords = read.table("https://raw.githubusercontent.com/mahmoudibrahim/KidneyMap/master/assets/reducedDims/UMAP/human_CD10positive_umapCoords.csv", header = FALSE, sep=",")
plot(mapCoords, col = color.gradient(scale(aaa), colors = c("gray",rev(c("#d7191c","#abdda4")))), pch = 20, cex = 0.4)
legend.col(col = colorRampPalette(c("gray",rev(c("#d7191c","#abdda4"))))(100), lev = seq(min(aaa), max(aaa), length.out=100))
zzz = dev.off()



#manual select markers
markersp = toupper(c("Pdgfrb","Pdgfra","Col1a1","DCN","PTPRC","EPCAM","PECAM1","KRT18","PAX8","RNASE1","CXCR4","PLP1","GFRA3","SOX10","EMCN","EGFL7","hsd11b2","clcnka","atp6v1b1","MMRN1","Prox1","c7","col15a1","nkd2","scara5","ugt2b7","glyat","rergl","mustn1","ackr1","vwf","ramp2","frzb","tagln","krt8","adamts6","igf2","ssuh2","GPR183","CD69","srp1","mitf","fgf7","sall1","pax2","dnmt1","tp53","acta2","vim","col1a1","mme","ugt1a5","slc16a7","dpep1","slc3a2","slc22a8","Fth1","ftl","aldob","ugt1a3","mt1h","mt1g","skp1","prdx2","Ech1","FBP1","Prdx6","mt1x","cldn1","slpi","NPL", "Slc4a4", "Slc36a2", "Slc13a1","UNK", "ACSM3", "DCDC2","Slc7a13", "Slc30a8", "RBP4"))
markersp = toupper(c("Slc22a6","Slc5a2","agt"))
mapCoords = read.table("https://raw.githubusercontent.com/mahmoudibrahim/KidneyMap/master/assets/reducedDims/UMAP/human_CD10positive_umapCoords.csv", header = FALSE, sep=",")
names1 = toupper(unlist(lapply(strsplit(as.character(rownames(l_remove)), split = ";", fixed = T), function(x) as.character(x[2]))))
for (i in 1:length(markersp)) {
	if(any(names1 == markersp[i])) {
		aaax = as.vector(logcounts(l_remove)[which(names1 == markersp[i]),])
		zzz = pdf(paste0(prefix, "umap_", markersp[i], ".pdf"))
		plot(mapCoords, col = color.gradient(scale(aaax), colors = c("#ECECEC50",rev(c("#d7191c","#abdda4")))), pch = 20, cex = 0.5)
		legend.col(col = colorRampPalette(c("#ECECEC50",rev(c("#d7191c","#abdda4"))))(100), lev = seq(min(aaax), max(aaax), length.out=100))
		zzz = dev.off()
	}
}






#cell cycle
G1_S = read.table("https://raw.githubusercontent.com/mahmoudibrahim/KidneyMap/master/assets/public/cell_cycle/G1S")[[1]]
G2_M = read.table("https://raw.githubusercontent.com/mahmoudibrahim/KidneyMap/master/assets/public/cell_cycle/G2M")[[1]]
G0 = read.table("https://raw.githubusercontent.com/mahmoudibrahim/KidneyMap/master/assets/public/cell_cycle/G0")[[1]]
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
roses = logcounts(l_remove)
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
cycle_tab = t(table(cell_phase, l_remove$level3))
cycle_tab = t(apply(cycle_tab, 1, function(x) x/sum(x)))
zzz = pdf(paste0(prefix, "cell_cycle.pdf"))
heatmap.2(cycle_tab, trace = "none", scale="none", col = colorRampPalette(c("white","dodgerblue4"))(n=100), Rowv = NA, Colv = NA, margins = c(15,15), key = FALSE, cellnote= round(cycle_tab,2)*100, notecol = "black", notecex = 2)
dev.off()



##Fin
