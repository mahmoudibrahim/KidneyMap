########################################################################
# ecm_score_demo.r
#
# Data describes bulk kidney RNA-seq from Diabetic Nephropathy Patients
# and Control (healthy part) of tumor nephrectomies. 
# Data source: Fan et al. Diabetes, 2019. 
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE142025 (GSE142025_RAW.tar) 
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

#load packages
#plotting
suppressPackageStartupMessages(library(vioplot))
makeTransparent = function(..., alpha=0.5) {
	alpha = floor(255*alpha)  
	newColor = col2rgb(col=unlist(list(...)), alpha=FALSE)
	.makeTransparent = function(col, alpha) {
		rgb(red=col[1], green=col[2], blue=col[3], alpha=alpha, maxColorValue=255)
	}
	newColor = apply(newColor, 2, .makeTransparent, alpha=alpha)
	return(newColor)
}
col_palette_short = c("#999999","#56B4E9","#009E73","#0072B2")
col_palette_short_trans = makeTransparent(col_palette_short, alpha = 0.3)


#get data
files = list.files(path=".", pattern="*.txt", full.names=TRUE, recursive=FALSE)
dat = list()
for (f in 1:length(files)) {
	dat[[f]] = read.table(files[[f]], header = TRUE)[[2]]
}
dat = do.call(cbind, dat)
rownames(dat) = read.table(files[[1]], header = TRUE)[[1]]
colnames(dat) = files
cond = c(rep("Advanced_DN", 21), rep("Early_DN", 6), rep("Ctrl", 9))


#matrisome gene set analysis
matrisome_set = read.table("https://raw.githubusercontent.com/mahmoudibrahim/KidneyMap/master/assets/public/ecm_genes_human.txt", sep = "\t", header = TRUE)



#collagens
collagens = matrisome_set$Gene.Symbol[which(matrisome_set$Category %in% c("Collagens"))]
reads_single_phase = dat
reads_single_phase = as.matrix(reads_single_phase[rownames(reads_single_phase) %in% (collagens) ,])
combined_matrix = rbind(reads_single_phase,average=apply(reads_single_phase,2,mean))
cor_matrix = cor(t(combined_matrix))
cor_vector = cor_matrix[,dim(cor_matrix)[1]]
reads_single_phase_restricted = reads_single_phase[rownames(reads_single_phase) %in% names(cor_vector[cor_vector >= 0.1]),]
col_score = apply(reads_single_phase_restricted,2,mean)


#glycoproteins
collagens = matrisome_set$Gene.Symbol[which(matrisome_set$Category %in% c("ECM Glycoproteins"))]
reads_single_phase = dat
reads_single_phase = as.matrix(reads_single_phase[rownames(reads_single_phase) %in% (collagens) ,])
combined_matrix = rbind(reads_single_phase,average=apply(reads_single_phase,2,mean))
cor_matrix = cor(t(combined_matrix))
cor_vector = cor_matrix[,dim(cor_matrix)[1]]
reads_single_phase_restricted = reads_single_phase[rownames(reads_single_phase) %in% names(cor_vector[cor_vector >= 0.1]),]
gp_score = apply(reads_single_phase_restricted,2,mean)


#proteoglycans
collagens = matrisome_set$Gene.Symbol[which(matrisome_set$Category %in% c("Proteoglycans"))]
reads_single_phase = dat
reads_single_phase = as.matrix(reads_single_phase[rownames(reads_single_phase) %in% (collagens) ,])
combined_matrix = rbind(reads_single_phase,average=apply(reads_single_phase,2,mean))
cor_matrix = cor(t(combined_matrix))
cor_vector = cor_matrix[,dim(cor_matrix)[1]]
reads_single_phase_restricted = reads_single_phase[rownames(reads_single_phase) %in% names(cor_vector[cor_vector >= 0.1]),]
pg_score = apply(reads_single_phase_restricted,2,mean)


#ECM score
collagens = matrisome_set$Gene.Symbol[which(matrisome_set$Division %in% c("Core matrisome"))]
reads_single_phase = dat
reads_single_phase = as.matrix(reads_single_phase[rownames(reads_single_phase) %in% (collagens) ,])
combined_matrix = rbind(reads_single_phase,average=apply(reads_single_phase,2,mean))
cor_matrix = cor(t(combined_matrix))
cor_vector = cor_matrix[,dim(cor_matrix)[1]]
reads_single_phase_restricted = reads_single_phase[rownames(reads_single_phase) %in% names(cor_vector[cor_vector >= 0.1]),]
ecm_score = apply(reads_single_phase_restricted,2,mean)



pdf("ecm_score.pdf", width = 8)
facts = c(rep(1, 21), rep(2, 6), rep(3, 9))
par(mfrow = c(2,2))

plot(facts + rnorm(length(facts),0,0.08), ecm_score, pch = 19, col = makeTransparent(col_palette_short[facts], 0.8), cex = 1.5,  main = paste0("ECM Score"), xaxt = "n", xlab = "Condition", ylab = "ECM Score", xlim = c(0.5,3.5))
axis(1, at=1:3, c("Advanced_DN", "Early_DN", "Ctrl"))
vioplot(ecm_score ~ facts, plotCentre="line", col = col_palette_short_trans, border = col_palette_short, lwd = 3, add = TRUE)

plot(facts + rnorm(length(facts),0,0.08), col_score, pch = 19, col = makeTransparent(col_palette_short[facts], 0.8), cex = 1.5,  main = paste0("Collagen Score"), xaxt = "n", xlab = "Condition", ylab = "Collagen Score", xlim = c(0.5,3.5))
axis(1, at=1:3, c("Advanced_DN", "Early_DN", "Ctrl"))
vioplot(col_score ~ facts, plotCentre="line", col = col_palette_short_trans, border = col_palette_short, lwd = 3, add = TRUE)

plot(facts + rnorm(length(facts),0,0.08), pg_score, pch = 19, col = makeTransparent(col_palette_short[facts], 0.8), cex = 1.5,  main = paste0("Proteoglycan Score"), xaxt = "n", xlab = "Condition", ylab = "Proteoglycan Score", xlim = c(0.5,3.5))
axis(1, at=1:3, c("Advanced_DN", "Early_DN", "Ctrl"))
vioplot(pg_score ~ facts, plotCentre="line", col = col_palette_short_trans, border = col_palette_short, lwd = 3, add = TRUE)

plot(facts + rnorm(length(facts),0,0.08), gp_score, pch = 19, col = makeTransparent(col_palette_short[facts], 0.8), cex = 1.5,  main = paste0("Glycoprotein Score"), xaxt = "n", xlab = "Condition", ylab = "Glycoprotein Score", xlim = c(0.5,3.5))
axis(1, at=1:3, c("Advanced_DN", "Early_DN", "Ctrl"))
vioplot(gp_score ~ facts, plotCentre="line", col = col_palette_short_trans, border = col_palette_short, lwd = 3, add = TRUE)

dev.off()
