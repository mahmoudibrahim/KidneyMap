
	for (scorpions in 1:10) {
	
		message(scorpions)
	
	
		###batch correct
		kk = floor(sqrt(ncol(l)) * 1) * ncol(l)
		newOrder = order(table(l$colData$fileID), decreasing = TRUE)
	
		mnn = mnnCorrect(as.matrix(logcounts(l[,which(l$colData$file == (unique(l$colData$file)[1]))])), ..., k = kk/(ncol(l)), subset.row=var_genes, cos.norm.out=FALSE, order = newOrder, cos.norm.in=TRUE) #here one has to list as many matrices as there are samples

		l_bc = cbind(mnn$corrected[[1]], mnn$corrected[[2]], ...) #here one has to list as many matrices as there are samples
		l_bc = SingleCellExperiment(assays=list(logcounts=l_bc))
		l_bc$colData = l$colData



		###SVD dimension reduction
		l_v = l_bc[which(rownames(l_bc) %in% var_genes),]
		q = t(logcounts(l_v))
		q = t(apply(q,1, function(x) x-mean(x)))
		sv = svd(q)
		d = sv$d
		d2 =  d^2 / (sum(d^2))

		###take the knee
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
		dw = 1:dw


		###UMAP (optional, just to have this but this is kind of slow.)
		temp = t(apply(ss$u[,dw], 1, function(x) x / (sqrt(sum(x^2))))) #norm based normalization
		load("umap.defaults.RData")
		umap.defaults$n_neighbors = kk/(nrow(ss$u))
		set.seed(111)
		rr = umap(temp, config = umap.defaults)
		rr = rr$layout
		umapCoords[[scorpions]] = rr

		
		###Graph Clustering (choose resolution)
		if (nrow(sv$u) > 1000) { #cell number check
		seqs = seq(0.1,1,length.out=10)
		seqsres = rep(0, length(seqs))
		for (ii in 1:length(seqs)) {
			resolution = seqs[ii]
			kk = floor(sqrt(nrow(ss$u)) * resolution) * nrow(ss$u)
			
			kn = nn2(data = temp, k = ceiling((kk*1) / ncol(l_v)), searchtype = "priority", treetype = "bd")
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
			colors_info = rep("", length = ncol(l_v))
			for (i in 1:length(unique(class_info))) {
				wtemp = which((class_info) == (sort(unique(class_info)))[i])
				colors_info[wtemp] = col_palette_trans[i]
			}

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


		#now final clustering
		resolution = seqs[which.max(seqsres)]
		} else {resolution = 1} ##cell number check
		message(resolution)
		kk = floor(sqrt(nrow(ss$u)) * resolution) * nrow(ss$u)
		kn = nn2(data = temp, k = ceiling((kk*1) / ncol(l_v)), searchtype = "priority", treetype = "bd")
		kn = data.frame(node2 = as.vector(kn$nn.idx[,-1]), node1 = rep(kn$nn.idx[,1], ncol(kn$nn.idx) - 1), sim = 1 / (1+as.vector(kn$nn.dists[,-1])))
		cut = min(head(sort(kn$sim, decreasing = TRUE), n = kk))
		kn$sim2 = kn$sim
		kn$sim2[kn$sim2 < cut] = 0
		kn$sim2[kn$sim2 != 0] = 1
		edge_list = cbind(kn$node1[which(kn$sim2 == 1)], kn$node2[which(kn$sim2 == 1)])

		g = igraph::graph.edgelist(edge_list, directed = FALSE)
		g = simplify(g)
		graph[[scorpions]] = g
		set.seed(111)
		km = igraph::cluster_louvain(g)
		class_info = km$membership
		colors_info = rep("", length = ncol(l_v))
		for (i in 1:length(unique(class_info))) {
			wtemp = which((class_info) == (sort(unique(class_info)))[i])
			colors_info[wtemp] = col_palette_trans[i]
		}
		clusterLabels[[scorpions]] = class_info

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
		varGenes[[scorpions]] = var_genes
		meanClassAUC[[scorpions]] = mean(getClassAUC(sg, plotCurves = FALSE))


		###normalize expression
		minko = min(table(clusterLabels[[scorpions]])) #if there are clusters that are less than 10 cells
		if (minko > 10) {
			minko = 10
		}
		l = scran::computeSumFactors(l, sizes = seq(minko, 200, 20), clusters = clusterLabels[[scorpions]], positive = TRUE)
		ws = which(sizeFactors(l) <= 0)
		if(length(ws) > 0) {
			message("Warning: some cells had negative size factors!")
			sizeFactors(l)[ws] = min(sizeFactors(l)[sizeFactors(l) > 0])
		}
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
			if (llDiff < 0.1) { #when to break the loop
				break;
			}
		}
	}#scorpions
