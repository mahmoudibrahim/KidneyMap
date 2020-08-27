#Nov 13: pseudotime
#Nov 22: add cell cycle
#Dec 09: remove pseudotime, reoorder clusters
rm(list = ls())

##load packages
#plotting
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(vioplot))
#scRNA specific
suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(scran))
#IRLB (fast SVD)
suppressPackageStartupMessages(library(irlba))
#for graph clustering
suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(RANN))
#markers
suppressPackageStartupMessages(library(genesorteR))
suppressPackageStartupMessages(library(mclust))
#match cell types
suppressPackageStartupMessages(library(CHETAH))
#biomart
suppressPackageStartupMessages(library(biomaRt))


options(warn = -1, scipen = 1000)


convertHumanGeneList <- function(x){ 
	human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
	mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
	genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
	return(genesV2)
}





#load human reference data
load("/data/kramann/neil/human/pdgfrb/combined/v2/combined/pdgfrbMap_japanium.filtered_matrix.RData")
human = l_remove
names1 = unlist(lapply(strsplit(as.character(rownames(l_remove)), split = ";", fixed = T), function(x) as.character(x[2])))
names2 = convertHumanGeneList(names1)
inter = intersect(names1, names2$HGNC.symbol)
rownames(human) = names1
human = human[which(rownames(human) %in% inter),]
names2 = names2[which(names2$HGNC.symbol %in% inter),]

human = human[!(duplicated(rownames(human))),]
names2 = names2[!(duplicated(names2$HGNC.symbol)),]
human = human[order(rownames(human)),]
names2 = names2[order(names2$HGNC.symbol),]
rownames(human) = names2$MGI.symbol
human = human[!(duplicated(rownames(human))),] #final matrix with mouse gene names and no duplicates

#load query mouse data
load("/data/kramann/neil/mouse/pdgfrab/combined/map/pdgfrabMap_take2japanium.filtered_matrix.RData")
mouse = l_remove
names1 = unlist(lapply(strsplit(as.character(rownames(l_remove)), split = ";", fixed = T), function(x) as.character(x[2])))
rownames(mouse) = names1
mouse = mouse[!duplicated(rownames(mouse)),] #final mouse matrix 
reducedDims(mouse)$tSNE = as.matrix(read.table("/data/kramann/neil/mouse/pdgfrab/combined/map/pdgfrabMap_take2umapCoordsP_default.csv", sep = ","))

cheat = CHETAHclassifier(mouse, human, input_c = "logcounts", ref_c = "logcounts", ref_ct = "level3", thresh = 0.1)
pdf(paste0("cheatah_tree.pdf"), width = 20)
PlotCHETAH(input = cheat, tree = TRUE, interm = TRUE)
dev.off()
tab = table(cheat$level3, cheat$celltype_CHETAH)
tab[which(tab < 10)] = 0

tab = t(apply(tab, 1, function(x) x/sum(x)) * 100)
tab = tab[,-(which(colSums(tab) == 0))]

pdf(paste0("compare_mouse_human_heat.pdf"), width = 8)
pheatmap(t(tab), col = colorRampPalette(c("white","dodgerblue4"))(n=100), display_numbers=TRUE, number_color="black")
dev.off()
