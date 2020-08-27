
	for (scorpions in 1:10) {
	
		message(scorpions)
	
		#batch correct
		mnn = batchelor::fastMNN(l, batch = l$colData$fileID, d = 30, auto.order = TRUE, subset.row = var_genes, BSPARAM = BiocSingular::IrlbaParam(tol=1e-8))
		mnn_u = mnn@reducedDims$corrected
		gc()

		for (i in 2:30) {
			zzz = pdf(paste0(outName, "_PCA", i, ".pdf"))
			plot(mnn_u[,c(1,i)], col = color.gradient(scale(l$total_counts), colors = c("gray",rev(c("#d7191c","#fdae61","#abdda4")))), pch = 20)
			zzz = dev.off()
		}
		colors_file = rep("", length = ncol(l))
	
		for (i in 1:length(unique(l$colData$file))) {
			wtemp = which((l$colData$file) == (sort(unique(l$colData$file)))[i])
			colors_file[wtemp] = col_palette_trans[i]
		}
		for (i in 2:30) {
			zzz = pdf(paste0(outName, "_PCA", i, "_fileID.pdf"))
			plot(mnn_u[,c(1,i)], col = colors_file, pch = 16)
			zzz = dev.off()
		}
		
		
		###UMAP
		kk = floor(sqrt(ncol(l)) * 1)
		if (kk > 100) {
			kk = 100
		}
		kk = kk * ncol(l)
		
		temp = t(apply(mnn_u, 1, function(x) x / (sqrt(sum(x^2))))) #norm based scaling
		set.seed(111)
		rr = uwot::umap(temp, n_neighbors = kk / (nrow(mnn_u)), min_dist = 2, verbose = FALSE, fast_sgd = TRUE, spread = 5, n_threads=16)
	
		zzz = pdf(paste0(outName, "_japanium_UMAP_mnn_totalCounts_take", scorpions, ".pdf"))
		plot(rr, col = color.gradient(scale(l$total_counts), colors = c("gray",rev(c("#d7191c","#fdae61","#abdda4")))), pch = 20, xlab = "UMAP-1", ylab = "UMAP-2", frame.plot = FALSE, ylim = c(min(rr[,2])-10,max(rr[,2])+10))
		zzz = dev.off()


		colors_file = rep("", length = ncol(l))
		for (i in 1:length(unique(l$colData$fileID))) {
			wtemp = which((l$colData$fileID) == (sort(unique(l$colData$fileID)))[i])
			colors_file[wtemp] = col_palette_trans[i]
		}
		zzz = pdf(paste0(outName, "_japanium_UMAP_mnn_fileID_take", scorpions, ".pdf"))
		plot(rr, col = colors_file, pch = 20, xlab = "UMAP-1", ylab = "UMAP-2", frame.plot = FALSE, ylim = c(min(rr[,2])-10,max(rr[,2])+10))
		legend(min(rr[,1]), max(rr[,2]) + 10, legend = sort(unique(l$colData$fileID)), pch=21, col=col_palette, pt.bg=col_palette_trans, pt.cex=2, cex=.8, bty="n")
		zzz = dev.off()


		colors_file = rep("", length = ncol(l))
		for (i in 1:length(unique(l$colData$origCellCluster))) {
			wtemp = which((l$colData$origCellCluster) == (sort(unique(l$colData$origCellCluster)))[i])
			colors_file[wtemp] = col_palette_trans[i]
		}
		zzz = pdf(paste0(outName, "_japanium_UMAP_mnn_originalCluster_take", scorpions, ".pdf"))
		plot(rr, col = colors_file, pch = 20, xlab = "UMAP-1", ylab = "UMAP-2", frame.plot = FALSE, ylim = c(min(rr[,2])-10,max(rr[,2])+10))
		legend(min(rr[,1]), max(rr[,2]) + 10, legend = sort(unique(l$colData$origCellCluster)), pch=21, col=col_palette, pt.bg=col_palette_trans, pt.cex=2, cex=.8, bty="n")
		zzz = dev.off()
		umapCoords[[scorpions]] = rr

		
		###Graph Clustering (choose resolution)
		if (nrow(mnn_u) > 1000) {
		
		seqs = c(0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.09, seq(0.1,1,length.out=10))
		seqsres = rep(0, length(seqs))
		for (ii in 1:length(seqs)) {
			resolution = seqs[ii]
			kk = floor(sqrt(nrow(mnn_u) * resolution)) * nrow(mnn_u)
			if (floor(sqrt(nrow(mnn_u) * resolution)) < 10) {
				kk = 10 * nrow(mnn_u)
			}
			kn = nn2(data = temp, k = ceiling((kk*1) / nrow(mnn_u)), searchtype = "priority", treetype = "bd")
			kn = data.frame(node2 = as.vector(kn$nn.idx[,-1]), node1 = rep(kn$nn.idx[,1], ncol(kn$nn.idx) - 1), sim = 1 / (1+as.vector(kn$nn.dists[,-1])))
			cut = min(head(sort(kn$sim, decreasing = TRUE), n = kk))
			kn$sim2 = kn$sim
			kn$sim2[kn$sim2 < cut] = 0
			kn$sim2[kn$sim2 != 0] = 1
			edge_list = cbind(kn$node1[which(kn$sim2 == 1)], kn$node2[which(kn$sim2 == 1)])


			g = igraph::graph.edgelist(edge_list, directed = FALSE)
			g = simplify(g)
			set.seed(111)
			km = igraph::cluster_louvain(g)
			class_info = km$membership
			colors_info = rep("", length = nrow(mnn_u))
			for (i in 1:length(unique(class_info))) {
				wtemp = which((class_info) == (sort(unique(class_info)))[i])
				colors_info[wtemp] = col_palette_trans[i]
			}

			zzz = pdf(paste0(outName, "_japanium_UMAP_mnn-take", scorpions, "-", resolution, ".pdf"))
			plot(rr, col = colors_info, pch = 20, xlab = "UMAP-1", ylab = "UMAP-2", frame.plot = FALSE, ylim = c(min(rr[,2])-10,max(rr[,2])+10))
			legend(min(rr[,1]), max(rr[,2]) + 10, legend = 1:length(unique(class_info)), pch=21, col=col_palette, pt.bg=col_palette_trans, pt.cex=2, cex=.8, bty="n", ncol = ceiling(length(unique(class_info)) / 2))
			zzz = dev.off()

	
	
			###merge singular clusters
			if (any(table(class_info) == 1)) {
				www = which(table(class_info) == 1)
				www = which(class_info %in% www)
				class_info[www] = 1 #merge with the largest cluster
				message(paste0(length(www), " cells were merged because they formed singular clusters."))
			}
		
			###determine variable genes
			sg = sortGenes(logcounts(l), class_info, cores = 1, binarizeMethod = "naive")
			var_genes = unique(unlist(plotTopMarkerHeat(sg, top_n = 500, outs = T, plotheat=F)))
			seqsres[ii] = mean(getClassAUC(sg, plotCurves=F))

		}###
		gc()

		#now final clustering
		resolution = seqs[which.max(seqsres)]
		} else {resolution = 1} ##cell number check
		message(resolution)

		kk = floor(sqrt(nrow(mnn_u) * resolution)) * nrow(mnn_u)
		if (floor(sqrt(nrow(mnn_u) * resolution)) < 10) {
			kk = 10 * nrow(mnn_u)
		}
		kn = nn2(data = temp, k = ceiling((kk*1) / nrow(mnn_u)), searchtype = "priority", treetype = "bd")
		kn = data.frame(node2 = as.vector(kn$nn.idx[,-1]), node1 = rep(kn$nn.idx[,1], ncol(kn$nn.idx) - 1), sim = 1 / (1+as.vector(kn$nn.dists[,-1])))
		cut = min(head(sort(kn$sim, decreasing = TRUE), n = kk))
		kn$sim2 = kn$sim
		kn$sim2[kn$sim2 < cut] = 0
		kn$sim2[kn$sim2 != 0] = 1
		edge_list = cbind(kn$node1[which(kn$sim2 == 1)], kn$node2[which(kn$sim2 == 1)])

		g = igraph::graph.edgelist(edge_list, directed = FALSE)
		g = simplify(g)
		set.seed(111)
		km = igraph::cluster_louvain(g)
		class_info = km$membership
		colors_info = rep("", length = nrow(mnn_u))
		for (i in 1:length(unique(class_info))) {
			wtemp = which((class_info) == (sort(unique(class_info)))[i])
			colors_info[wtemp] = col_palette_trans[i]
		}
		clusterLabels[[scorpions]] = class_info
		zzz = pdf(paste0(outName, "_japanium_UMAP_mnn-take", scorpions, ".pdf"))
		plot(rr, col = colors_info, pch = 20, xlab = "UMAP-1", ylab = "UMAP-2", frame.plot = FALSE, ylim = c(min(rr[,2])-10,max(rr[,2])+10))
		legend(min(rr[,1]), max(rr[,2]) + 10, legend = 1:length(unique(class_info)), pch=21, col=col_palette, pt.bg=col_palette_trans, pt.cex=2, cex=.8, bty="n", ncol = ceiling(length(unique(class_info)) / 2))
		zzz = dev.off()

		###merge singular clusters
		if (any(table(class_info) == 1)) {
			www = which(table(class_info) == 1)
			www = which(class_info %in% www)
			class_info[www] = 1 #merge with the largest cluster
			message(paste0(length(www), " cells were merged because they formed singular clusters."))
		}
		gc()


		###determine variable genes
		sg = sortGenes(logcounts(l), class_info, cores = 16, binarizeMethod = "naive")
		zzz = pdf(paste0(outName, "_japanium_genesorteR_clusterQuality-take", scorpions, ".pdf"), height = 10)
		barplot(getClassAUC(sg), ylim = c(0,1), yaxp = c(0,1,20))
		zzz = dev.off()
		var_genes = unique(unlist(plotTopMarkerHeat(sg, top_n = 500, outs = T, plotheat=F)))
		varGenes[[scorpions]] = var_genes
		meanClassAUC[[scorpions]] = mean(getClassAUC(sg, plotCurves = FALSE))


		###normalize expression
		minko = min(table(clusterLabels[[scorpions]]))
		if (minko > 10) {
			minko = 10
		}
		l = scran::computeSumFactors(l, sizes = seq(minko, 200, 20), clusters = clusterLabels[[scorpions]], positive = TRUE)
		ws = which(sizeFactors(l) <= 0)
		if(length(ws) > 0) {
			message("Warning: some cells had negative size factors!")
			sizeFactors(l)[ws] = min(sizeFactors(l)[sizeFactors(l) > 0])
		}
		zzz = pdf(paste0(outName, "_japanium_scran_sizeFactors_take", scorpions, ".pdf"))
		plot(sizeFactors(l), l$total_counts/1e3, xlab = "Size factors", ylab = "Library Depth (thousands)", col = "#99999950", pch = 20)
		zzz = dev.off()
		l = scater::normalize(l)
		l = scater::calculateQCMetrics(l)
		a = plotExprsFreqVsMean(l)
	
	
		###get stats
		if (scorpions == 1) {
			rand[[1]] = adjustedRandIndex(l$colData$origCellCluster, clusterLabels[[1]])
		} else {
			rand[[scorpions]] = adjustedRandIndex(clusterLabels[[scorpions-1]],clusterLabels[[scorpions]])

			llDiff = (rand[[scorpions]] - rand[[scorpions-1]]) / rand[[scorpions]]
			message(llDiff)
			if (llDiff < 0.1) {
				break;
			}
		}
	}#scorpions
