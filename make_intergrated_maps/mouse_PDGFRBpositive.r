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

#general color palette
col_palette = c("#d64e3c","#7dd554","#8641c6","#cfca45","#7778cb","#59803d","#d04d9c","#73d6a8","#492f60","#ccc497","#7f343b","#72acc0","#b97d40","#c796b5","#45483a","purple","green","yellow")
col_palette_trans = makeTransparent(col_palette, alpha = 0.5)


#define color scheme
level1_col = c("#42006a","#9A6324")
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


###########################################################################
########################START ANALYSIS#####################################
###########################################################################

#load data
load("kidney_mesen_cleansed.RData") ##from Hendersson lab
exp = kidney_mesen_normexprs

#some filtering
ercc = exp[(nrow(exp)-91):(nrow(exp)),]
ercc = ercc[which(rowSums(ercc) > 0),]
eGFP = exp[which(rownames(exp) == "eGFP"),]
Pdgfrb = exp[which(rownames(exp) == "Pdgfrb"),] 
exp = exp[1:(nrow(exp)-93),]

#cell filtering
exp = exp[,-(which(eGFP == 0 & Pdgfrb == 0))] #take out cells that do not express both eGFP and Pdgfrb
ercc = ercc[,-(which(eGFP == 0 & Pdgfrb == 0))]
kidney_mesen_metadata$Clustering = kidney_mesen_metadata$Clustering[-(which(eGFP == 0 & Pdgfrb == 0))] #hend. lab clusters
kidney_mesen_metadata$Timepoint = kidney_mesen_metadata$Timepoint[-(which(eGFP == 0 & Pdgfrb == 0))] #experiment time-point

#gene filtering (1% of cells)
pctExp = apply(exp, 1, function(x) (length(x[x>0]) / length(x)) * 100)
exp = exp[which(pctExp >= 1),]
rRNA = grep("^Mrp|^Rpl|^Rps", rownames(exp), value=T) #remove ribosomal proteins
exp = exp[-(which(rownames(exp) %in% rRNA)),]



#############################
###Clustering, Iteration 1###
#############################

#determine variable genes
fit = scran::trendVar(ercc, parametric = TRUE)
decomp = scran::decomposeVar(exp, fit = fit)
mcut = which(decomp$FDR < 0.01 & decomp$bio > 1)
message(paste0("Number of variable genes is ", length(mcut), "."))

pdf("NormCVvsMean.pdf")
par(mfrow = c(2,2))
plot(decomp$mean, decomp$total, xlab="Mean log-expression", ylab="Total Variance", pch = 16, col = "#99999930")
points(decomp$mean[mcut], decomp$total[mcut], pch = 16, col = "red")
plot(decomp$mean, decomp$bio, xlab="Mean log-expression", ylab="Biological Variance", pch = 16, col = "#99999930")
points(decomp$mean[mcut], decomp$bio[mcut], pch = 16, col = "red")
plot(decomp$bio, decomp$total, xlab="Biological Variance", ylab="Total Variance", pch = 16, col = "#99999930")
points(decomp$bio[mcut], decomp$total[mcut], pch = 16, col = "red")
zzz = dev.off()


#dimension reduction
names_v = rownames(cv[mcut,])
mat_v = as.matrix(exp[which(rownames(exp) %in% names_v),])
q = t(mat_v)
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
zzz = pdf(paste0("pcComp.pdf"))
plot(1:100, d[1:100], xlim = c(1,100), xaxp = c(1,100,99), pch = 19, col = "#99999949")
lines(rep(dw, 1000), seq(0,10000,length.out = 1000), col = "red", lty = 2)
zzz = dev.off()
dw = 1:dw



pdf("SVD_clusterColor.pdf")
plot(ss$u, col = as.factor(kidney_mesen_metadata$Clustering), pch = 16)
dev.off()
pdf("SVD_timeColor.pdf")
plot(ss$u, col = as.factor(kidney_mesen_metadata$Timepoint), pch = 16)
dev.off()
temp = t(apply(ss$u[,dw], 1, function(x) x / (sqrt(sum(x^2)))))
write.csv(temp, file = "mouse_PDGFRBpositive_dims.txt", row.names = FALSE, col.names = FALSE) #file saved @https://raw.githubusercontent.com/mahmoudibrahim/KidneyMap/master/assets/reducedDims/SVD/mouse_PDGFRBpositive_dims.txt


#clustering
allInd = 1:nrow(ss$u)
kk = floor(sqrt(nrow(ss$u)) * 1) * nrow(ss$u)
set.seed(111)
kn = nn2(data = temp, k = ceiling((kk*1) / nrow(ss$u)), searchtype = "priority", treetype = "bd")
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
km = igraph::cluster_louvain(g)
class_info = km$membership
colors_info = rep("", length = ncol(exp))
for (i in 1:length(unique(class_info))) {
	wtemp = which((class_info) == sort(unique(class_info)))[i]
	colors_info[wtemp] = col_palette_trans[i]
}


#UMAP (produced in python, just for visualization, independent of clustering)
rr = read.table("https://raw.githubusercontent.com/mahmoudibrahim/KidneyMap/master/assets/reducedDims/UMAP/mouse_PDGFRBpositive_umapCoords.csv", header = FALSE, sep=",") 
pdf("umap_clusterColor.pdf")
plot(rr, col = as.factor(kidney_mesen_metadata$Clustering), pch = 16)
dev.off()
pdf("umap_timeColor.pdf")
plot(rr, col = as.factor(kidney_mesen_metadata$Timepoint), pch = 16)
dev.off()
pdf(paste0("eGFP", ".pdf"))
aaa = eGFP
plot(rr, col = color.gradient(scale(aaa), colors = c("gray",rev(c("#d7191c","#abdda4")))), pch = 19, xlab = "umap-1", ylab = "umap-2", frame.plot = FALSE)
legend.col(col = colorRampPalette(c("gray",rev(c("#d7191c","#abdda4"))))(100), lev = seq(min(aaa), max(aaa), length.out=100))
dev.off()
zzz = pdf(paste0("japanium_UMAP.pdf"))
plot(rr, col = colors_info, pch = 20, xlab = "UMAP-1", ylab = "UMAP-2")
legend(min(rr[,1]), max(rr[,2]), legend = sort(unique(class_info))), pch=21, col=col_palette, pt.bg=col_palette_trans, pt.cex=2, cex=.8, bty="n", ncol = ceiling(length(unique(class_info)) / 2))
zzz = dev.off()


#get cell cluster gene markers
sg = sortGenes(exp, class_info, binarizeMethod = "adaptiveMedian")
markers = unique(unlist(plotTopMarkerHeat(sg, top_n = 500, outs = T, plotheat=F)))
mcut = which(rownames(exp) %in% markers)



#############################
###Clustering, Iteration 2###
#############################

l_v = exp[mcut,] #highly variable genes are the gene markers from Iteration 1

q = t(l_v)
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
zzz = pdf(paste0("pcComp.pdf"))
plot(1:100, d[1:100], xlim = c(1,100), xaxp = c(1,100,99), pch = 19, col = "#99999949")
lines(rep(dw, 1000), seq(0,10000,length.out = 1000), col = "red", lty = 2)
zzz = dev.off()
dw = 1:dw

temp = t(apply(ss$u[,dw], 1, function(x) x / (sqrt(sum(x^2)))))
write.csv(temp, file = "mouse_PDGFRBpositive_dims_take2.txt", row.names = FALSE, col.names = FALSE) #file saved @https://raw.githubusercontent.com/mahmoudibrahim/KidneyMap/master/assets/reducedDims/SVD/mouse_PDGFRBpositive_dims_take2.txt



#clustering
allInd = 1:nrow(ss$u)
kk = floor(sqrt(nrow(ss$u)) * 1) * nrow(ss$u)

set.seed(111)
kn = nn2(data = temp, k = ceiling((kk*1) / nrow(ss$u)), searchtype = "priority", treetype = "bd")
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
km = igraph::cluster_louvain(g)
class_info = km$membership
colors_info = rep("", length = ncol(exp))
for (i in 1:length(unique(class_info))) {
	wtemp = which((class_info) == (sort(unique(class_info)))[i])
	colors_info[wtemp] = col_palette_trans[i]
}


#UMAP (produced in python, just for visualization, independent of clustering)
rr = read.table("https://raw.githubusercontent.com/mahmoudibrahim/KidneyMap/master/assets/reducedDims/UMAP/mouse_PDGFRBpositive_umapCoords_take2.csv", header = FALSE, sep=",") 
pdf("umap_clusterColor.pdf")
plot(rr, col = as.factor(kidney_mesen_metadata$Clustering), pch = 16)
dev.off()
pdf("umap_timeColor.pdf")
plot(rr, col = as.factor(kidney_mesen_metadata$Timepoint), pch = 16)
dev.off()
pdf(paste0("eGFP", ".pdf"))
aaa = eGFP
plot(rr, col = color.gradient(scale(aaa), colors = c("gray",rev(c("#d7191c","#abdda4")))), pch = 19, xlab = "umap-1", ylab = "umap-2", frame.plot = FALSE)
legend.col(col = colorRampPalette(c("gray",rev(c("#d7191c","#abdda4"))))(100), lev = seq(min(aaa), max(aaa), length.out=100))
dev.off()
zzz = pdf(paste0("japanium_UMAP.pdf"))
plot(rr, col = colors_info, pch = 20, xlab = "UMAP-1", ylab = "UMAP-2")
legend(min(rr[,1]), max(rr[,2]), legend = 1:length(unique(class_info)), pch=21, col=col_palette, pt.bg=col_palette_trans, pt.cex=2, cex=.8, bty="n", ncol = ceiling(length(unique(class_info)) / 2))
zzz = dev.off()


#############################
#####clustering done#########
#############################





#load cluster annotation info
annot = read.table("https://raw.githubusercontent.com/mahmoudibrahim/KidneyMap/master/assets/clusterInfo/mouse_PDGFRBpositive.txt", sep = "\t")
nums = table(class_info)
annot = annot[match(names(nums), annot$V1),]
annot$V6 = as.numeric(nums)
annot$V7 = as.numeric(log10(nums))
perc = t(table(kidney_mesen_metadata$Timepoint, class_info))
perc = (t(apply(perc, 1, function(x) x/sum(x)))) * 100
annot = cbind(annot, perc)

colnames(annot) = c("ID", "Label", "Annotation 1", "Annotation 2", "Annotation 3", "Number of Cells", "log10 Number of Cells", "Percent of Cells (No Injury)", "Percent of Cells (Day 2)", "Percent of Cells (Day 10)")
##write cluster annotation information
write.csv(annot, "cluster_annotation_info_withAddData.csv", row.names = FALSE)

level1 = class_info
level2 = class_info
level3 = class_info
for (i in 1:length(unique(class_info))) {
	where = which(class_info == unique(class_info)[i])
	level3[where] = as.character(annot[[5]][which(annot[[1]] == unique(class_info)[i])])
	level2[where] = as.character(annot[[4]][which(annot[[1]] == unique(class_info)[i])])
	level1[where] = as.character(annot[[3]][which(annot[[1]] == unique(class_info)[i])])
}





#get gene markers for different annotation levels (ranking for each gene in each cell cluster)
sg_level1 = sortGenes(exp, level1, binarizeMethod = "naive", cores = 16)
write.files(sg_level1, prefix = paste0("genesorteRouts_level1"))
sg_level2 = sortGenes(exp, level2, binarizeMethod = "naive", cores = 16)
write.files(sg_level2, prefix = paste0("genesorteRouts_level2"))
sg_level3 = sortGenes(exp, level3, binarizeMethod = "naive", cores = 16)
write.files(sg_level3, prefix = paste0("genesorteRouts_level3"))



#define cell colors for plotting
colors_info_level1 = rep("", length = ncol(exp))
for (i in 1:length(unique(level1))) {
	wtemp = which(level1 == (sort(unique(level1)))[i])
	colors_info_level1[wtemp] = level1_col_trans[i]
}
colors_info_level2 = rep("", length = ncol(exp))
ww = which(level1 == "Epithelial")
for (i in 1:length(unique(level2[ww]))) {
	wtemp = which(level2 == (sort(unique(level2[ww])))[i])
	colors_info_level2[wtemp] = level2_col_two_trans[i]
}
ww = which(level1 == "Mesenchymal")
for (i in 1:length(unique(level2[ww]))) {
	wtemp = which(level2 == (sort(unique(level2[ww])))[i])
	colors_info_level2[wtemp] = level2_col_four_trans[i]
}


#make publication umap
pdf(paste0("superUMAP_level1_level3TEXT.pdf"))
plot(rr[,1], rr[,2], pch = 19, col = colors_info_level1, cex = 0.6)
allClasses = unique(level3)
for (i in 1:length(allClasses)) {
	ww = which(level3 == allClasses[i])
	xx = median(rr[ww,1])
	yy = median(rr[ww,2])
	text(xx,yy,allClasses[i], cex = 0.4, adj = 0.5)
}
dev.off()
pdf(paste0("superUMAP_level1_level2TEXT.pdf"))
plot(rr[,1], rr[,2], pch = 19, col = colors_info_level1, cex = 0.6)
allClasses = unique(level2)
for (i in 1:length(allClasses)) {
	ww = which(level2 == allClasses[i])
	xx = median(rr[ww,1])
	yy = median(rr[ww,2])
	text(xx,yy,allClasses[i], cex = 0.65, adj = 0.5)
}
dev.off()



#make publication heatmaps
level3_order = c(1,6,5,9,4,2,8,7,3,10)
sg_temp = sortGenes(exp, class_info, binarizeMethod = "naive", cores = 16)
pdf(paste0("level3_top_auto_markers_10.pdf"), height = 15)
plotTopMarkerHeat(sg_temp, averageCells=10^5, newOrder=level3_order, gaps = FALSE, colors = colorRampPalette( c("cornflowerblue","black","gold"))(n=100), top_n=10)
dev.off()

pdf("select_markers_wholeClusterAverage.pdf")
hend_markers =  c("Pdgfrb", "Itgb1","Vim", "Acta2", "Rgs5", "Myh11", "Cnn1", "Ephx3", "Pax8", "Cldn1","Itih4", "Efhd1","Gja5", "Ren1", "Esr1", "Pdgfra", "Col1a1","Col15a1")
plotBinaryHeat(sg_temp$binary, sg_temp$inputClass, hend_markers, clusterGenes = F, averageCells=10^6, newOrder=level3_order, colors = colorRampPalette(c("black","gold"))(n=100))
dev.off()

timeline = table(kidney_mesen_metadata$Timepoint, class_info)
rownames(timeline) = c("d0","d2","d10")
timeline2 = apply(timeline, 2, function(x) x/sum(x))
timeline3 = t(apply(timeline, 1, function(x) x/sum(x)))
pdf("timeHeat_mahmoudCluster_level3.pdf")
heatmap.2((t(timeline3)*100)[level3_order,], trace = "none", col = colorRampPalette(c("white","dodgerblue4"))(n=100), scale = "none", cellnote = (round(t(timeline3)*100, 2))[level3_order,], notecol = "black", notecex = 1.5, dendrogram = "none", Rowv = NA, Colv = NA, breaks = seq(0,100,length.out=101),key = FALSE)
dev.off()

timeline = table(kidney_mesen_metadata$Timepoint, level2)
rownames(timeline) = c("d0","d2","d10")
timeline3 = t(apply(timeline, 1, function(x) x/sum(x)))
level2_order = c(4,6,2,5,3,1)
pdf("timeHeat_mahmoudCluster_level2.pdf")
heatmap.2((t(timeline3)*100)[level2_order,], trace = "none", col = colorRampPalette(c("white","dodgerblue4"))(n=100), scale = "none", cellnote = (round(t(timeline3)*100, 2))[level2_order,], notecol = "black", notecex = 1.5, dendrogram = "none", Rowv = NA, Colv = NA, breaks = seq(0,100,length.out=101),key = FALSE, margin = c(10,20))
dev.off()


#make marker UMAPs for some genes
hend_markers = c("Pdgfrb", "Pdgfra", "Acta2", "Cnn1", "Col1a1", "Col1a2", "Bcam", "Cldn1", "Cnn1", "Vim", "Cd34", "Ephx3", "Dcn", "Ncam1", "Miki67", "Col3a1","Mcam","Runx1","Ren1","Nkd2","Col15a1","Scara5","Rspo3","Postn","Col14a1","Ogn","Comp","Crlf1","Spon2")
names1 = rownames(exp)
for (i in 1:length(hend_markers)) {
	if(any(names1 == hend_markers[i])) {
		aaa = as.vector(exp[which(names1 == hend_markers[i]),])
		pdf(paste0("hendersson_markers_", hend_markers[i], ".pdf"))
		plot(rr, col = color.gradient(scale(aaa), colors = c("#ECECEC50",rev(c("#d7191c","#abdda4")))), pch = 19, xlab = "umap-1", ylab = "umap-2", frame.plot = FALSE, ylim = c(min(rr[,2]),max(rr[,2])))
		legend.col(col = colorRampPalette(c("#ECECEC50",rev(c("#d7191c","#abdda4"))))(100), lev = seq(min(aaa), max(aaa), length.out=100))
		dev.off()
	}
}





#matrisome gene set analysis
matrisome_set = read.table("https://raw.githubusercontent.com/mahmoudibrahim/KidneyMap/master/assets/public/ecm_genes_mouse.txt", sep = "\t", header = TRUE)

#collagens
collagens = matrisome_set$Gene.Symbol[which(matrisome_set$Category %in% c("Collagens"))]
reads_single_phase = exp
reads_single_phase = as.matrix(reads_single_phase[rownames(reads_single_phase) %in% (collagens) ,])
combined_matrix = rbind(reads_single_phase,average=apply(reads_single_phase,2,mean))
cor_matrix = cor(t(combined_matrix))
cor_vector = cor_matrix[,dim(cor_matrix)[1]]
reads_single_phase_restricted = reads_single_phase[rownames(reads_single_phase) %in% names(cor_vector[cor_vector >= 0.1]),]
aaa = apply(reads_single_phase_restricted,2,mean)
pdf("collagens_score.pdf")
plot(rr, col = color.gradient(scale(aaa), colors = c("#ECECEC50",rev(c("#d7191c","#abdda4")))), pch = 19, xlab = "umap-1", ylab = "umap-2", frame.plot = FALSE, ylim = c(min(rr[,2]),max(rr[,2])))
legend.col(col = colorRampPalette(c("#ECECEC50",rev(c("#d7191c","#abdda4"))))(100), lev = seq(min(aaa), max(aaa), length.out=100))
dev.off()
pdf("collagens_score_violinplot.pdf")
level2[which(level2 == "(Myo)fibroblast")] = "Myofibroblast"
facts = factor(as.factor(level2),levels(as.factor(level2))[c(4,6,1,5,2,3)])
plot(as.numeric(facts)+rnorm(length(aaa), 0, 0.1), aaa, pch = 19, xlab = "", ylab = "", cex = 0.8, main = "Collagens Score", col = colors_info_level1, xaxt = "n")
axis(1, at=sort(unique(as.numeric(facts))), labels=abbreviate(levels(facts)))
vioplot(as.vector(aaa) ~ as.factor(facts), col = makeTransparent("white", alpha = 0), border = makeTransparent("black", alpha = 0.3), range = 0.1, outline = F, add = T, notch = T, lwd = 4,  plotCentre="line", xlab = "Cell Clusters", ylab = "Log Expression")
dev.off()


#glycoproteins
matrisome_set = read.table("/home/mibrahim/Dropbox/Kidney_Interstitium_collab_Henderson/Public_data/Matrisome/ecm_genes_mouse.txt", sep = "\t", header = TRUE)
collagens = matrisome_set$Gene.Symbol[which(matrisome_set$Category %in% c("ECM Glycoproteins"))]
reads_single_phase = exp
reads_single_phase = as.matrix(reads_single_phase[rownames(reads_single_phase) %in% (collagens) ,])
combined_matrix = rbind(reads_single_phase,average=apply(reads_single_phase,2,mean))
cor_matrix = cor(t(combined_matrix))
cor_vector = cor_matrix[,dim(cor_matrix)[1]]
reads_single_phase_restricted = reads_single_phase[rownames(reads_single_phase) %in% names(cor_vector[cor_vector >= 0.1]),]
aaa = apply(reads_single_phase_restricted,2,mean)
pdf("glycoproteins_score.pdf")
plot(rr, col = color.gradient(scale(aaa), colors = c("#ECECEC50",rev(c("#d7191c","#abdda4")))), pch = 19, xlab = "umap-1", ylab = "umap-2", frame.plot = FALSE, ylim = c(min(rr[,2]),max(rr[,2])))
legend.col(col = colorRampPalette(c("#ECECEC50",rev(c("#d7191c","#abdda4"))))(100), lev = seq(min(aaa), max(aaa), length.out=100))
dev.off()
pdf("glycoproteins_score_violinplot.pdf")
level2[which(level2 == "(Myo)fibroblast")] = "Myofibroblast"
facts = factor(as.factor(level2),levels(as.factor(level2))[c(4,6,1,5,2,3)])
plot(as.numeric(facts)+rnorm(length(aaa), 0, 0.1), aaa, pch = 19, xlab = "", ylab = "", cex = 0.8, main = "Glycoproteins Score", col = colors_info_level1, xaxt = "n")
axis(1, at=sort(unique(as.numeric(facts))), labels=abbreviate(levels(facts)))
vioplot(as.vector(aaa) ~ as.factor(facts), col = makeTransparent("white", alpha = 0), border = makeTransparent("black", alpha = 0.3), range = 0.1, outline = F, add = T, notch = T, lwd = 4,  plotCentre="line", xlab = "Cell Clusters", ylab = "Log Expression")
dev.off()


#proteoglycans
matrisome_set = read.table("/home/mibrahim/Dropbox/Kidney_Interstitium_collab_Henderson/Public_data/Matrisome/ecm_genes_mouse.txt", sep = "\t", header = TRUE)
collagens = matrisome_set$Gene.Symbol[which(matrisome_set$Category %in% c("Proteoglycans"))]
reads_single_phase = exp
reads_single_phase = as.matrix(reads_single_phase[rownames(reads_single_phase) %in% (collagens) ,])
combined_matrix = rbind(reads_single_phase,average=apply(reads_single_phase,2,mean))
cor_matrix = cor(t(combined_matrix))
cor_vector = cor_matrix[,dim(cor_matrix)[1]]
reads_single_phase_restricted = reads_single_phase[rownames(reads_single_phase) %in% names(cor_vector[cor_vector >= 0.1]),]
aaa = apply(reads_single_phase_restricted,2,mean)
pdf("proteoglycans_score.pdf")
plot(rr, col = color.gradient(scale(aaa), colors = c("#ECECEC50",rev(c("#d7191c","#abdda4")))), pch = 19, xlab = "umap-1", ylab = "umap-2", frame.plot = FALSE, ylim = c(min(rr[,2]),max(rr[,2])))
legend.col(col = colorRampPalette(c("#ECECEC50",rev(c("#d7191c","#abdda4"))))(100), lev = seq(min(aaa), max(aaa), length.out=100))
dev.off()
pdf("proteoglycans_score_violinplot.pdf")
level2[which(level2 == "(Myo)fibroblast")] = "Myofibroblast"
facts = factor(as.factor(level2),levels(as.factor(level2))[c(4,6,1,5,2,3)])
plot(as.numeric(facts)+rnorm(length(aaa), 0, 0.1), aaa, pch = 19, xlab = "", ylab = "", cex = 0.8, main = "Proteoglycans Score", col = colors_info_level1, xaxt = "n")
axis(1, at=sort(unique(as.numeric(facts))), labels=abbreviate(levels(facts)))
vioplot(as.vector(aaa) ~ as.factor(facts), col = makeTransparent("white", alpha = 0), border = makeTransparent("black", alpha = 0.3), range = 0.1, outline = F, add = T, notch = T, lwd = 4,  plotCentre="line", xlab = "Cell Clusters", ylab = "Log Expression")
dev.off()


#all ECM genes
matrisome_set = read.table("/home/mibrahim/Dropbox/Kidney_Interstitium_collab_Henderson/Public_data/Matrisome/ecm_genes_mouse.txt", sep = "\t", header = TRUE)
collagens = matrisome_set$Gene.Symbol[which(matrisome_set$Division %in% c("Core matrisome"))]
reads_single_phase = exp
reads_single_phase = as.matrix(reads_single_phase[rownames(reads_single_phase) %in% (collagens) ,])
combined_matrix = rbind(reads_single_phase,average=apply(reads_single_phase,2,mean))
cor_matrix = cor(t(combined_matrix))
cor_vector = cor_matrix[,dim(cor_matrix)[1]]
reads_single_phase_restricted = reads_single_phase[rownames(reads_single_phase) %in% names(cor_vector[cor_vector >= 0.1]),]
aaa = apply(reads_single_phase_restricted,2,mean)
pdf("ecm_score.pdf")
plot(rr, col = color.gradient(scale(aaa), colors = c("#ECECEC50",rev(c("#d7191c","#abdda4")))), pch = 19, xlab = "umap-1", ylab = "umap-2", frame.plot = FALSE, ylim = c(min(rr[,2]),max(rr[,2])))
legend.col(col = colorRampPalette(c("#ECECEC50",rev(c("#d7191c","#abdda4"))))(100), lev = seq(min(aaa), max(aaa), length.out=100))
dev.off()
pdf("ECM_score_violinplot.pdf")
level2[which(level2 == "(Myo)fibroblast")] = "Myofibroblast"
facts = factor(as.factor(level2),levels(as.factor(level2))[c(4,6,1,5,2,3)])
plot(as.numeric(facts)+rnorm(length(aaa), 0, 0.1), aaa, pch = 19, xlab = "", ylab = "", cex = 0.8, main = "ECM Score", col = colors_info_level1, xaxt = "n")
axis(1, at=sort(unique(as.numeric(facts))), labels=abbreviate(levels(facts)))
vioplot(as.vector(aaa) ~ as.factor(facts), col = makeTransparent("white", alpha = 0), border = makeTransparent("black", alpha = 0.3), range = 0.1, outline = F, add = T, notch = T, lwd = 4,  plotCentre="line", xlab = "Cell Clusters", ylab = "Log Expression")
dev.off()




#cell cycle analysis
https://raw.githubusercontent.com/mahmoudibrahim/KidneyMap/master/assets/public/cell_cycle/G0
G1_S = read.table("https://raw.githubusercontent.com/mahmoudibrahim/KidneyMap/master/assets/public/cell_cycle/G1S")[[1]]
G2_M = read.table("https://raw.githubusercontent.com/mahmoudibrahim/KidneyMap/master/assets/public/cell_cycle/G2M")[[1]]
G0 = read.table("https://raw.githubusercontent.com/mahmoudibrahim/KidneyMap/master/assets/public/cell_cycle/G0")[[1]]


cycle_list = list(G1_S, G2_M, G0)
reads_single = exp
rownames(reads_single) = toupper(rownames(reads_single))
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
ans_normed <- flexible_normalization(ans,by_row=FALSE)
ans_normed_normed <- flexible_normalization(ans_normed,by_row=TRUE)
cell_phase <- apply(ans_normed_normed,1,function(x) colnames(ans_normed_normed)[which.max(x)])

cycle_tab = t(table(cell_phase, level2))[level2_order,]
pdf("cell_cycle_level2.pdf", width = 10)
barplot((apply(cycle_tab, 1, function(x) x/sum(x)) * 100), beside=T, col = c("#999999","#56B4E9","#0072B2"), ylim = c(0,60), yaxp = c(0,60,12), names = abbreviate(rownames(cycle_tab)), ylab = "% of Cells")
legend("topright", legend = colnames(cycle_tab), pch=15, col = c("#999999","#56B4E9","#0072B2"), pt.bg =  c("#999999","#56B4E9","#0072B2"), bty = "n", cex = 1.5)
dev.off()

cycle_tab = t(table(cell_phase, level3))[level3_order,]
pdf("cell_cycle_level3.pdf", width = 10)
barplot((apply(cycle_tab, 1, function(x) x/sum(x)) * 100), beside=T, col = c("#999999","#56B4E9","#0072B2"), ylim = c(0,60), yaxp = c(0,60,12), names = abbreviate(rownames(cycle_tab)), ylab = "% of Cells")
legend("topright", legend = colnames(cycle_tab), pch=15, col = c("#999999","#56B4E9","#0072B2"), pt.bg =  c("#999999","#56B4E9","#0072B2"), bty = "n", cex = 1.5)
dev.off()


##Fin
