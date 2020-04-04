########################################################################
# Seurat_compare_PDGFRab_mouse.r
#
# This contains some typical Seurat clustering to compare to our own 
# clustering for the PDGFRab mouse data.
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
options(scipen=1000000)


#plotting
makeTransparent = function(..., alpha=0.5) {
	alpha = floor(255*alpha)  
	newColor = col2rgb(col=unlist(list(...)), alpha=FALSE)
	.makeTransparent = function(col, alpha) {
		rgb(red=col[1], green=col[2], blue=col[3], alpha=alpha, maxColorValue=255)
	}
	newColor = apply(newColor, 2, .makeTransparent, alpha=alpha)
	return(newColor)
}
col_palette = c("#d64e3c","#7dd554","#8641c6","#cfca45","#7778cb","#59803d","#d04d9c","#73d6a8","#492f60","#ccc497","#7f343b","#72acc0","#b97d40","#c796b5","#45483a","purple","green","yellow","cornflowerblue")
col_palette_trans = makeTransparent(col_palette, alpha = 0.4)




#prefix and load data
prefix = "Seurat_demo_" #prefix for output
load("pdgfrabMap_japanium.filtered_matrix.RData") ## this contains data and our analysis


#Seurat processing and clustering, MNN for batch correction
library(Seurat)
library(SeuratWrappers)
library(gplots)
seu = Seurat::CreateSeuratObject(counts = counts(l_remove), project = "pdgfr_ab", min.cells = 3, min.features = 200)
seu@meta.data$experiment = l_remove$fileID
seu = Seurat::NormalizeData(seu)
seu = Seurat::FindVariableFeatures(seu)
seu = RunFastMNN(object.list = SplitObject(seu, split.by = "experiment"))
seu = FindNeighbors(seu, reduction = "mnn", dims = 1:30)
seu = FindClusters(seu, resolution = 0.6)

class_info = Idents(seu)
colors_info = rep("", length = ncol(counts(l_remove)))
for (i in 1:length(unique(Idents(seu)))) {
	wtemp = which((Idents(seu)) == (unique(sort(Idents(seu)))[i]))
	colors_info[wtemp] = col_palette_trans[i]
}


#compare 
tt = table(l_remove$level3, Idents(seu))
tt = apply(tt, 2, function(x) x / sum(x))
tt = tt * 100
pdf(paste0(prefix, "_compare_seurat_level3.pdf"), width = 8)
heatmap.2(tt, trace = "none", dendrogram="none", col = colorRampPalette(c("white","dodgerblue4"))(n=100), margins = c(5,12), cellnote=round(tt,1), notecol="black", key = FALSE)
dev.off()
tt = table(l_remove$level2, Idents(seu))
tt = apply(tt, 2, function(x) x / sum(x))
tt = tt * 100
pdf(paste0(prefix, "_compare_seurat_level2.pdf"), width = 8)
heatmap.2(tt, trace = "none", dendrogram="none", col = colorRampPalette(c("white","dodgerblue4"))(n=100), margins = c(5,15), cellnote=round(tt,1), notecol="black", key = FALSE)
dev.off()


class_info = l_remove$level3
colors_info2 = rep("", length = ncol(Idents(seu)))
for (i in 1:length(unique(class_info))) {
	wtemp = which((class_info) == (sort(unique(class_info)))[i])
	colors_info2[wtemp] = col_palette_trans[i]
}

class_info = l_remove$level2
colors_info3 = rep("", length = ncol(Idents(seu)))
for (i in 1:length(unique(class_info))) {
	wtemp = which((class_info) == (sort(unique(class_info)))[i])
	colors_info3[wtemp] = col_palette_trans[i]
}


class_info = l_remove$fileID
colors_info4 = rep("", length = ncol(Idents(seu)))
for (i in 1:length(unique(class_info))) {
	wtemp = which((class_info) == (sort(unique(class_info)))[i])
	colors_info4[wtemp] = col_palette_trans[i]
}


#plot UMAP for Seurat
mapCoords = read.table("", header = FALSE, sep=",")##add link to UMAP coords
pdf(paste0(prefix, "_compare_seurat_seuratcolors.pdf"))
plot(mapCoords, col = colors_info, pch = 19, cex = 0.6)
legend(min(mapCoords[,1]), min(mapCoords[,2]) + 1, legend = abbreviate(sort(unique(Idents(seu)))), pch=21, col=col_palette, pt.bg=col_palette_trans, pt.cex=2, cex=.8, bty="n", ncol = ceiling(length(unique(Idents(seu))) / 2))
dev.off()
