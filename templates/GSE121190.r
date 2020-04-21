########################################################################
# GSE121190.r
#
# process expression microarray PDGFRa+ data from Higashi et al. 2019
# PMID: 30545984
#
# This script is generated automatically through the GEO website, up to 
# "summarize gene info"
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


library(Biobase)
library(GEOquery)

setwd("where do you want to save the output")
gset <- getGEO("GSE121190", GSEMatrix =TRUE, AnnotGPL=FALSE) #pulls data from GEO (needs internet)
if (length(gset) > 1) idx <- grep("GPL11180", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group names for all samples
gsms <- "XXXXXXXXX000111"
sml <- c()
for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }

# eliminate samples marked as "X"
sel <- which(sml != "X")
sml <- sml[sel]
gset <- gset[ ,sel]

# log2 transform
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
          (qx[6]-qx[1] > 50 && qx[2] > 0) ||
          (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { ex[which(ex <= 0)] <- NaN
  exprs(gset) <- log2(ex) }


#summarize gene info
geneInfo = data.frame(ID = as.character(fit$genes$ID), GS = as.character(fit$genes$Gene.Symbol), coef = exprs(gset))
geneInfo = geneInfo[-which(geneInfo$GS == ""),]
tt = aggregate(. ~ GS, data = geneInfo, max, drop = TRUE)
exp = as.matrix(tt[,3:8])
rownames(exp) = tt$GS
colnames(exp) = c("uuo1","uuo2","uuo3","sham1","sham2","sham3")
pdgfra_fib = exp

save(pdgfra_fib, file = "GSE121190_exprs.RData")
