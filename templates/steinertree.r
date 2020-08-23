#from Afshin Sadeghi and Holger Froehlich, "Steiner tree methods for optimal sub-network identification: an empirical study", BMC Bioinformatics 2013 14:144 doi:10.1186/1471-2105-14-144.
#https://github.com/krashkov/SteinerNet 

# All shortest paths between terminals (ASP)
asp_steiner <- function (optimize, terminals, glist, color) {
        g <- glist[[1]]
    
        paths <- lapply(terminals, function (x) get.all.shortest.paths(g, x, terminals)$res)
        #nodes <- unique(unlist(paths))
        nodes <- unique(names(unlist(paths)))
        
        if (optimize) {
                steinert <- minimum.spanning.tree(induced_subgraph(graph = g, vids = nodes))
                a   <- V(steinert)$color
                b   <- igraph::degree(steinert, v = V(steinert), mode = c("all"))
                a1  <- match(a, "yellow")
                b1  <- match(b, "1")
                opt <- sapply(1:length(a1), function (r) a1[r] * b1[r])
                new_g <- delete.vertices(steinert, grep(1, opt))
                steinert <- new_g
        } else
                steinert <- induced_subgraph(graph = g, vids = nodes)
        
        glst <- c()
        
        if (color) {
                V(g)[setdiff(x = nodes, y = terminals)]$color <- "green"
                glst[[length(glst) + 1]] <- g
        }
    
        glst[[length(glst) + 1]] <- steinert
    
        return(glst)
}



# Randomized all shortest paths approximation (RSP)
appr_steiner <- function (repeattimes, optimize, terminals, glist, color) {
        set <- c()
        g <- glist[[1]]
        
        # Start with the sub-graph G* consisting of all nodes and edges appearing on shortest paths between terminals
        paths <- lapply(terminals, function (x) get.all.shortest.paths(g, x, terminals)$res)
        
        r  <- 1:length(paths)
        t1 <- lapply(r, function (r) length(paths[[r]]))
        
        distances <- lapply(r, function (r) lapply(1:t1[[r]], function(x, y) length(paths[[y]][[x]]), y = r))
        neighbour_distance <- max(unlist(distances))
        
        # Note, graph has to have name attribute, because we assign names of vertices to
        # path variable. It is much more convenient to work with names, not with ids.
        #paths <- unique(unlist(paths))
        paths <- unique(names(unlist(paths)))
        #set   <- V(g)[paths]
        set   <- V(g)[paths]$name
        size  <- length(E(minimum.spanning.tree(induced_subgraph(g, union(terminals, set)))))
        
        i <- 1
        while (i <= repeattimes) {
        	#seed_list <- unlist(neighborhood(graph = g, order = neighbour_distance, nodes = terminals, mode = "all"))
                seed_list <- names(unlist(neighborhood(graph = g, order = neighbour_distance, nodes = terminals, mode = "all")))
                seed_list <- seed_list[!(seed_list %in% terminals)]
                seed <- sample(seed_list, 1)
            
                paths2 <- get.all.shortest.paths(g, seed, terminals)
                paths2 <- paths2$res

                #seedpaths <- unique(unlist(paths2))
                seedpaths <- unique(names(unlist(paths2)))
                        
                set2  <- union(set, V(g)[seedpaths]$name)
                size2 <- length(E(minimum.spanning.tree(induced_subgraph(g, union(terminals, set2)))))
                        
                if (size2 < size) {
                        size <- size2
                        set  <- set2
                }
                        
                seed  <- sample(set, 1, prob = NULL)
                set2  <- V(g)[setdiff(set, seed)]$name
                size2 <- length(E(minimum.spanning.tree(induced_subgraph(g, union(terminals, set2)))))
            
                if (size2  < size && is.connected(minimum.spanning.tree(induced_subgraph(g, union(terminals, set2))))) {
                        size <- size2
                        set  <- set2
                }
                        
                i <- i + 1
        }
        
        # Perform "optimization": find minimum spanning tree and remove nodes of degree 1
        if (optimize) {
                steinert <- minimum.spanning.tree(induced_subgraph(g, union(terminals, set)))
                a     <- V(steinert)$color
                b     <- igraph::degree(steinert, v = V(steinert), mode = c("all"))
                a1    <- match(a, "yellow")
                b1    <- match(b, "1")
                opt   <- sapply(1:length(a1), function(r) a1[r] * b1[r])
                new_g <- delete.vertices(steinert, grep(1, opt))
                steinert <- new_g
        } else
                steinert <- induced_subgraph(g, union(terminals, set))
    
        glst <- c()
        if (color) {
                V(g)[setdiff(set, terminals)]$color <- "green"
                glst[[length(glst) + 1]]  <- g
        }
    
        glst[[length(glst) + 1]] <- steinert
    
        return(glst)
}



# Shortest Path Based Approximation (SP)
steinertree2 <- function (optimize, terminals, glist, color) {
        g <- glist[[1]]
        
        # Pick a terminal randomly and Form a subtree (sub-graph G')
        prob     <- sample(1:length(terminals), 1)
        subtree  <- terminals[[prob]]
        nsubtree <- setdiff(terminals, subtree)
        
        # Proceed until all terminals not in G'
        while ( !all(is.element(terminals, intersect(subtree, terminals))) ) {
                # Compute shortest paths and their lengths between each node in subtree (G') and the remaining nodes
                paths <- lapply(subtree, function (x) get.all.shortest.paths(g, x, nsubtree))
                
                r <- 1:length(paths)
                t <- sapply(r, function (r) sapply(paths[[r]]$res, length))
        
                # Compute a minimum for each set of lengths from each node to other nodes
                if (class(t) == "list" || class(t) == "integer") {
                        r  <- 1:length(t)
                        t2 <- sapply(r, function (r) min(t[[r]]))
                }
                if (class(t) == "matrix") {
                        r  <- 1:dim(t)[2]
                        t2 <- sapply(r, function (r) min(t[, r]))
                }
        
                # Find a path with minimum among minimum length
                t3 <- which(t2 == min(t2))
                
                # Note, graph has to have name attribute, because in found variable we assign names
                # of vertices. It is much more convenient to work with names, not with ids.
                if (length(paths) > 1) {
                        if (class(t) == "list" || class(t) == "integer")
                                t4 <- which(t[[t3[1]]] == min(t[[t3[1]]]))
            
                        if (class(t) == "matrix")
                                t4 <- which( t[ , t3[1]] == min(t[ , t3[1]]) )
                        
                        #found <- unlist(paths[[t3[1]]][t4][1]$res)
                        found <- names(unlist(paths[[t3[1]]][t4][1]$res))
                } else {
                	#found <- unlist(paths[[1]][t3][1]$res)
                	found <- names(unlist(paths[[1]][t3][1]$res))
                }
                
                # Add all vertices from all shortest paths to subtree
                #subtree  <- union(subtree, V(g)[unique(found)])
                subtree  <- union(subtree, V(g)[unique(found)]$name)
                #nsubtree <- setdiff(nsubtree, V(g)[unique(found)])
                nsubtree <- setdiff(nsubtree, V(g)[unique(found)]$name)
        }
        
        # Perform "optimization": find minimum spanning tree and remove nodes of degree 1
        if (optimize) {
                steinert <- minimum.spanning.tree(induced_subgraph(g, subtree))
                a   <- V(steinert)$color
                b   <- igraph::degree(steinert, v = V(steinert), mode = c("all"))
                a1  <- match(a, "yellow")
                b1  <- match(b, "1")
                opt <- sapply(1:length(a1), function (r) a1[r] * b1[r])
                new_g <- delete.vertices(steinert, grep(1, opt))
                steinert <- new_g
        } else
                steinert <- induced_subgraph(g, subtree)
    
        glst <- c()
        if (color) {
                V(g)[subtree]$color   <- "green"
                V(g)[terminals]$color <- "red"
                
                glst[[length(glst) + 1]] <- g
        }
    
        glst[[length(glst) + 1]] <- steinert
    
        return(glst)
}



# Minimum spanning tree based approximation (Kruskal's minimum spanning tree algorithm)
steinertree3 <- function (optimize, terminals, glist, color) {
        makesubtrees <- function (x) {
                if ( !is.na(any(match(t3, x))) )
                	#return(union(subtrees[[x]],
                	#	     found[[grep(1, match(t3, x))]][[1]]))
                        return(union(subtrees[[x]],
                                     names(found[[grep(1, match(t3, x))]][[1]])))
                else return(subtrees[[x]])
        }
        
        subtreenum <- c()
        x <- c()
        g <- glist[[1]]
    
        # Make a Streiner Tree from every terminal
        r <- 1:length(terminals)
        subtrees  <- lapply(r, function (r) terminals[[r]])
        terminals <- subtrees
        nsubtrees <- lapply(r, function (r) setdiff(terminals, subtrees[r]))
    
        # Proceed until all terminals won't be added to a subtree
        while (length(subtrees) > 1) {
        	# Find shortest paths between different Steiner Trees and compute their lengths
                r     <- 1:length(subtrees)
                #paths <- lapply(r, function (r) lapply(subtrees[[r]],
                #				       function (x, y) get.all.shortest.paths(g, x, y)$res,
                #				       y = nsubtrees[[r]]))
                paths <- lapply(r, function (r) lapply(subtrees[[r]],
                                                       function (x, y) get.all.shortest.paths(g, x, y)$res,
                                                       y = unlist(nsubtrees[[r]])))
                
                r <- 1:length(paths)
                t <- sapply(r, function (r) sapply(paths[[r]][[1]], length))
                
                # Compute a minimum for each set of lengths from each Steiner tree to other trees
                if (class(t) == "list" | class(t) == "integer") {
                        r  <- 1:length(t)
                        t2 <- sapply(r, function (x) min(t[[x]]))
                }
                if (class(t) == "matrix") {
                        r  <- 1:dim(t)[2]
                        t2 <- sapply(r, function (r) min(t[, r]))
                }
                
                # Find a minimum among minimum length and paths corresponding to it
                t3    <- which(t2 == min(t2))
                t3len <- 1:length(t3)
                
                if (length(paths) > 1) {
                        if (class(t) == "list" || class(t) == "integer" )
                                t4 <- lapply(t3len, function (x) which(t[[t3[x]]] == min(t[[t3[x]]])))
                        if (class(t) == "matrix")
                                t4 <- lapply(t3len, function (x) which((t[ , t3[x]]) == min(t[ , t3[x]])))
                        
                        found <- lapply( t3len, function (x) paths[t3[x]][[1]][[1]][t4[[x]][1]] )
                } else {
                        intersect(subtrees[[x]], V(g)[unlist(terminals)])
                        print("Error")
                }
        
                # Merge subgraphs and paths
                subtrees <- lapply(1:length(subtrees), function (x) makesubtrees(x))
        
                # Delete repeated subtrees (presume that length is more than 1)
                i <- 1
                j <- 2
                while (i <= (length(subtrees) - 1)) {
                        j <- i + 1
                        while (j <= length(subtrees)) {
                                if (length(intersect(subtrees[[i]], subtrees[[j]])) > 0) {
                                        subtrees[[i]] <- union(subtrees[[i]], subtrees[[j]])
                                        subtrees <- subtrees[-j]
                                        j <- j - 1
                                }
                                j <- j + 1
                        }
                        i <- i + 1
                }
                nsubtrees <- lapply(1:length(subtrees), function (x) setdiff(terminals, subtrees[[x]]))
        }
        
        # Perform "optimization": find minimum spanning tree and remove nodes of degree 1
        if (optimize) {
                steinert <- minimum.spanning.tree(induced_subgraph(g, subtrees[[1]]))
                a   <- V(steinert)$color
                b   <- igraph::degree(steinert, v = V(steinert), mode = c("all"))
                a1  <- match(a, "yellow")
                b1  <- match(b, "1")
                opt <- sapply(1:length(a1), function (r) a1[r] * b1[r] )
                new_g <- delete.vertices(steinert, grep(1, opt))
                steinert <- new_g
        } else
                steinert <- induced_subgraph(g, subtrees[[1]])
    
        glst <- c()
        if (color) {
                V(g)[subtrees[[1]]]$color     <- "green"
                V(g)[unlist(terminals)]$color <- "red"
                #V(g)[terminals]$color <- "red"
                
                glst[[length(glst) + 1]]  <- g
        }
    
        glst[[length(glst) + 1]] <- steinert
    
        return(glst)
}



# Sub-graph of merged steiner trees (SPM or STM)
steinertree8 <- function (optimize, terminals, glist, color) {
        g <- glist[[1]]
        
        queue         <- c()
        results_queue <- c()
        edgeslist     <- c()
    
        prob <- sample(1:length(terminals), 1)
    
        subtree    <- terminals[[prob]]
        nsubtree   <- setdiff(terminals, subtree)
        startpoint <- subtree
    
        paths <- get.all.shortest.paths(g, subtree, nsubtree)
        paths <- paths$res
    
        t  <- sapply(paths, length)
        t2 <- which(t == min(t))
    
        # Put in queue paths with minimal lengths
        for (i in 1:length(t2))
        	#queue[length(queue) + 1] <- paths[t2[i]]
                queue[[length(queue) + 1]] <- names(unlist(paths[t2[i]]))
    
        index <- length(t2)
        while (index > 0) {
                edgeslist <- queue[1]
                queue[1]  <- NULL
                index     <- index - 1
        
                if (length(intersect(unlist(terminals), unlist(edgeslist))) == length(terminals)) {
                #if (length(intersect(unlist(terminals), names(unlist(edgeslist)))) == length(terminals)) {
                        graph_is_new <- TRUE
            
                        if (length(results_queue) == 0)
                                results_queue[length(results_queue) + 1] <- edgeslist
            
                        for (count_path in 1:length(results_queue)) {
                                t1 <- unlist(edgeslist[[1]])
                                t2 <- unlist(results_queue[[count_path]])
                
                                if (length(union(t1, t2)) == length(t1))
                                        if (all(union(t1, t2) %in% t2))
                                                graph_is_new <- FALSE
                        }
                        
                        if (graph_is_new == TRUE)
                                results_queue[length(results_queue) + 1] <- edgeslist
                } else {
                	subtree  <- intersect(unlist(terminals), unlist(edgeslist))
                        #subtree  <- intersect(unlist(terminals), names(unlist(edgeslist)))
                        nsubtree <- setdiff(terminals, subtree)
                        
                        paths    <- get.all.shortest.paths(g, subtree[length(subtree)], nsubtree)
                        paths    <- paths$res
                        
                        t  <- sapply(paths, length)
                        t2 <- which(t == min(t))

                        for (i in 1:length(t2))
                        	#queue[[index + i]] <- union(unlist(edgeslist), unlist(paths[t2[i]]))
                                queue[[index + i]] <- union(unlist(edgeslist), names(unlist(paths[t2[i]])))
                        
                        index <- index + length(t2)
                }
        }
    
        paths <- results_queue
        t     <- sapply(paths, length)
        t2    <- which(t == min(t))
        queue <- paths[t2]
    
        steinert_list <- c()
        glst <- c()
    
        for (i in 1:length(t2)) {
                steinert = minimum.spanning.tree(induced_subgraph(g, queue[[i]]))
        
                if (optimize) {
                        a   <- V(steinert)$color
                        b   <- igraph::degree(steinert, v = V(steinert), mode = c("all"))
                        a1  <- match(a, "yellow")
                        b1  <- match(b, "1")
                        opt <- sapply(1:length(a1), function (r) a1[r] * b1[r])
                        new_g <- delete.vertices(steinert, grep(1, opt))
                        steinert <- new_g
                }
                
                if (color)
                        V(g)[queue[[i]]]$color <- "green"
        
                steinert_list[[length(steinert_list) + 1]] <- steinert
        }
    
        if (color) {
        	#V(g)[as.numeric(terminals)]$color <- "red"
                V(g)[terminals]$color <- "red"
                
                glst[[length(glst) + 1]] <- g
                glst[[length(glst) + 1]] <- steinert_list
        } else
                glst <- steinert_list
    
        return (glst)
}


# Exact algorithm
steinerexact <- function (terminals, glist, color) {
        rwhile <- function (lim) {
                if (get("runloop", envir = en)) {
                        r <- length(V(g)) - lim
                        
                        allcom <- combn(t[1:length(V(g))], r)
                        allmst <- lapply(1:dim(allcom)[2],
                                         function (x)  minimum.spanning.tree(induced_subgraph(g, allcom[ , x])))
                        assign("allmst", allmst, envir = en)
            
                        edgmst <- lapply(1:dim(allcom)[2],
                                         function (x) get.edgelist(allmst[[x]], names = TRUE))
                        assign("edgmst", edgmst, envir = en)
            
                        # Check connectivity
                        connectedlist <- lapply(1:dim(allcom)[2], function (x) is.connected(allmst[[x]]))
                        # Check terminals availability
                        withterminals <- lapply(1:dim(allcom)[2], function (x) all(is.element(terminals, V(allmst[[x]])$name)))
            
                        # Both previous conditions
                        smst <- lapply(1:dim(allcom)[2], function (x) connectedlist[[x]] && withterminals[[x]])
            
                        assign("runloop",   !is.element(TRUE, unlist(smst)),  envir = en)
                        assign("sol_place", get("sol_place", envir = en) + 1, envir = en)
                }
                return(smst)
        }
    
        g <- glist[[1]]
        t <- V(g)$name
    
        lim <- length(V(g)) - length(terminals)
    
        en <- new.env(hash = TRUE, parent = emptyenv(), size = NA)
        assign("runloop",  TRUE, envir = en)
        assign("sol_place",   0, envir = en)
        smst <- c()

        res <- lim:1
        sol <- sapply(res, function (x) rwhile(x))
    
        sol_place <- get("sol_place", envir = en)
        allmst    <- get("allmst",    envir = en)
        edgmst    <- get("edgmst",    envir = en)
    
        # Size of trees
        iter <- length(sol[[sol_place]])
        size <- lapply(1:iter, function (x) length(edgmst[[x]]) / 2)
        midresult <- lapply(1:iter, function (x) size[[x]] * as.integer(sol[[sol_place]][[x]]))
    
        min_len <- min(unlist(midresult)[unlist(midresult) > 0])
        poslist <- which(unlist(midresult) == min_len)
        stgraphlist <- allmst[poslist]
    
        stgraphlist2 <- c()
    
        if (color) {
                green_guys <- lapply(stgraphlist, function (x) V(x)$name)
                green_guys <- unique(unlist(green_guys))
        
                V(g)[green_guys]$color <- "green"
                #V(g)[as.numeric(terminals)]$color <- "red"
                V(g)[terminals]$color  <- "red"
        
                stgraphlist2[[length(stgraphlist2) + 1]] <- g
                stgraphlist2[[length(stgraphlist2) + 1]] <- stgraphlist
                stgraphlist <- stgraphlist2
        }
        
        return(stgraphlist)
}



merge_steiner <- function (treelist) {
    
        merged <- treelist[[1]]
    
        if (length(treelist) > 1) {
                for (i in 2:length(treelist))
                        merged <- union(merged, treelist[[i]])
        } else
                print("Nothing to merge. Only one solution was found")
    
        glist      <- c()
        glist[[1]] <- merged
    
        return(glist)
}



check_input <- function (type, terminals, glist) {
    
        g <- glist[[1]]
        g <- as.undirected(g)
    
        # Checking terminals
    
        if (is.null(terminals) || is.na(terminals) || length(terminals) == 0)
                stop("Error: Terminals not found")
    
        # Checking graph
    
        if (is.null(g))
                stop("Error: The graph object is Null.")
    
        if (length(V(g)) == 0 )
                stop("Error: The graph doesn't contain vertices.")
    
        if (is.null(V(g)$name)) {
                # creating name attribute
                V(g)$name <- as.character(1:length(V(g)))
                attr_flag <- FALSE
        } else {
                # creating new name and realname attributes
                V(g)$realname <- V(g)$name
                V(g)$name     <- as.character(1:length(V(g)))
                attr_flag <- TRUE
        }
    
        # Mathcing names of vertices and terminals, if possible
    
        if (class(terminals) == "character") {
                # terminals contain realname of vertices
                if (sum(terminals %in% V(g)$realname) != length(terminals)) {
                        stop("Error: vertices names do not contain terminal names")
                } else {
                        # Convert realnames of terminals to names (character id's)
                        terminals <- V(g)$name[match(terminals, V(g)$realname)]
                }
        } else if (class(terminals) == "numeric" | class(terminals) == "integer") {
                # terminals contains id's of vertices
                terminals <- V(g)$name[terminals]
        } else
                print("Error: invalid type of terminals")
    
        V(g)$color            <- "yellow"
        V(g)[terminals]$color <- "red"
    
        # Checking type
    
        if ( !(type == "SPM" | type == "EXA" | type == "SP" | type == "RSP" | type == "KB" | type == "ASP") )
                stop("Error: the input type is not correct. Choose one from SPM, EXA, SP, RSP or KB.")
        
        varlist      <- c()
        varlist[[1]] <- g
        varlist[[2]] <- terminals
        varlist[[3]] <- attr_flag
        
        return(varlist)
}



restore_name_attribute <- function (attr_flag, type, result, color) {
	
	if (color) {
		if (attr_flag) {
			V(result[[1]])$name <- V(result[[1]])$realname
			result[[1]] <- delete_vertex_attr(result[[1]], 'realname')
		}
	}
	
	if (type == "EXA" | type == "SPM") {
		if (attr_flag) {
			numSteiner <- length(result[[length(result)]])
			
			for (i in 1:numSteiner) {
				V(result[[length(result)]][[i]])$name <- V(result[[length(result)]][[i]])$realname
				result[[length(result)]][[i]] <- delete_vertex_attr(result[[length(result)]][[i]], 'realname')
			}
		}
        } else {
		if (attr_flag) {
			V(result[[length(result)]])$name <- V(result[[length(result)]])$realname
			result[[length(result)]] <- delete_vertex_attr(result[[length(result)]], 'realname')
		}
	}
	
	return(result)
}



####--------------------------------------- Documentation ---------------------------------------####
#' Find Steiner Tree
#' 
#' @description A set of functions for finding Steiner Tree. Includes both exact and heuristic approaches.
#' 
#' @usage steinertree(type, repeattimes = 70, optimize = TRUE, terminals,
#'             graph, color = TRUE, merge = FALSE)
#' 
#' @param type a character scalar, which indicates type of algorithms to perform. Can be
#'             "EXA", "SP", "KB", "RSP", "SPM" or "ASP".
#' @param repeattimes a numeric scalar to specify "RSP" algorithm; number of times the optimization procedure is repeated.
#' @param optimize a logical scalar to specify all algorithms except "EXA"; if TRUE, an optimization of the resultant
#'                 steiner tree is performed, otherwise nothing is done.
#' @param terminals a numeric vector (ids of terminals are passed) or character vector (vertices must have 'name' attribute).
#' @param graph an igraph graph; should be undirected, otherwise it is converted to undirected.
#' @param color a logical scalar; whether to return an original graph with terminals colored in red and
#'              steiner nodes colored in green. Note, if several trees will be found, steiner nodes from all trees
#'              are colored in green.
#' @param merge a logical scalar to specify "EXA" and "SPM" algorithms; if several trees will be found, whether to return
#'              a list with trees or merge them
#' 
#' @return (color = FALSE) Returns a list first element of which is a steiner tree (or a graph of merged trees).
#'         If several steiner trees are found, return a list, each element of which is a steiner tree.
#'          
#'         (color = TRUE) Returns a list, first element of which is a colored original graph and second element is
#'         a steiner tree (or a graph of merged trees) or list of steiner trees.
#'         
#' @details If input graph doesn't have 'name' attribute, one is created. In this case it will contain character ids of vertices.
#'          Also before execution all vertices will be colored in yellow and terminals will be colored in red.
#' 
#' @seealso \code{\link{generate_st_samples}}
#' 
#' @examples
#' steinertree(type = "RSP", optimize = FALSE,
#'             terminals = c(1, 3),
#'             graph = graph("Cubical"),
#'             color = TRUE, merge = FALSE)
#' 
#' @references 1. Path heuristic and Original path heuristic ,Section 4.1.3 of the book "The Steiner tree Problem",
#'                Petter,L,Hammer
#'                
#'             2. "An approximate solution for the Steiner problem in graphs", H Takahashi, A Matsuyama
#'             
#'             3. F K. Hwang, D S. Richards and P Winter, "The steiner tree Problem", Kruskal-Based Heuristic
#'                Section 4.1.4, ISBN: 978-0-444-89098-6
#'                
#'             4. Afshin Sadeghi and Holger Froehlich, "Steiner tree methods for optimal sub-network
#'                identification: an empirical study", BMC Bioinformatics 2013 14:144
#'                
#'             5. F K. Hwang, D S. Richards and P Winter, "The steiner tree Problem", Kruskal-Based Heuristic Section
#'                4.1.4, The Optimal solution for steiner trees on networks, ISBN: 978-0-444-89098-6.
#'             
#' @export
####------------------------------------- End Documentation -------------------------------------####
steinertree <- function (type, repeattimes = 70, optimize = TRUE, terminals, graph, color = TRUE, merge = FALSE) {
    
        glist      <- c()
        glist[[1]] <- graph
    
        varlist <- check_input(type = type, terminals = terminals, glist = glist)
    
        glist[[1]] <- varlist[[1]]
        terminals  <- varlist[[2]]
        attr_flag  <- varlist[[3]]
    
        if (type == "SP")
                result <- steinertree2(optimize = optimize, terminals = terminals, glist = glist, color = color)
    
        if (type == "KB")
                result <- steinertree3(optimize = optimize, terminals = terminals, glist = glist, color = color)
    
        if (type == "RSP")
                result <- appr_steiner(repeattimes = repeattimes, optimize = optimize, terminals = terminals,
                                       glist = glist, color = color)
    
        if (type == "EXA")
                result <- steinerexact(terminals = terminals, glist = glist, color = color)
    
        if (type == "SPM")
                result <- steinertree8(optimize = optimize, terminals = terminals, glist = glist, color = color)
    
        if (type == "ASP")
                result <- asp_steiner(optimize = optimize, terminals = terminals, glist = glist, color = color)
        
        result <- restore_name_attribute(attr_flag, type, result, color)
        
        if (merge & (type == "EXA" | type == "SPM")) {
                if (color) {
                        result[[2]] <- merge_steiner(treelist = result[[2]])
                } else {
                        result <- merge_steiner(treelist = result)
                }
        }
        
        return(result)
}

