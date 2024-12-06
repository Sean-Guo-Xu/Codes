

select_root_cell <- function(cds, root_state=NULL, reverse=FALSE){
  if (is.null(root_state) == FALSE) {
    if (is.null(pData(cds)$State)){
      stop("Error: State has not yet been set. Please call orderCells() without specifying root_state, then try this call again.")
    }
    # FIXME: Need to gaurd against the case when the supplied root state isn't actually a terminal state in the tree.
    root_cell_candidates <- subset(pData(cds), State == root_state)
    if (nrow(root_cell_candidates) == 0){
      stop(paste("Error: no cells for State =", root_state))
    }
    
    # build a local MST to find a good root cell for this state
    dp <- as.matrix(dist(t(reducedDimS(cds)[,row.names(root_cell_candidates)])))
    gp <- graph.adjacency(dp, mode = "undirected", weighted = TRUE)
    dp_mst <- minimum.spanning.tree(gp)
    
    # Make sure to use the real MST here
    tip_leaves <- names(which(degree(minSpanningTree(cds)) == 1))
    #root_cell_candidates <- root_cell_candidates[row.names(root_cell_candidates) %in% tip_leaves,]
    #sg <- make_ego_graph(dp_mst, nodes=row.names(root_cell_candidates))[[1]]
    
    # if(length(intersect(tip_leaves, row.names(root_cell_candidates))) == 0){
    #   stop(paste("Error: no valid root cells for State =", root_state))
    # }
    
    diameter <- get_diameter(dp_mst)
    
    if (length(diameter) == 0){
      stop(paste("Error: no valid root cells for State =", root_state))
    }
    
    #root_cell = names(diameter)[tip_leaves %in% names(diameter)]
    root_cell_candidates <- root_cell_candidates[names(diameter),]
    if (is.null(cds@auxOrderingData[[cds@dim_reduce_type]]$root_cell) == FALSE &&
        pData(cds)[cds@auxOrderingData[[cds@dim_reduce_type]]$root_cell,]$State == root_state){
      root_cell <- row.names(root_cell_candidates)[which(root_cell_candidates$Pseudotime == min(root_cell_candidates$Pseudotime))]
    }else{
      root_cell <- row.names(root_cell_candidates)[which(root_cell_candidates$Pseudotime == max(root_cell_candidates$Pseudotime))]
    }
    if (length(root_cell) > 1)
      root_cell <- root_cell[1]
    
    # If we used DDRTree, we need to go from this actual cell to the nearst
    # point on the principal graph
    if (cds@dim_reduce_type == "DDRTree"){
      #root_cell_idx <- which(V(minSpanningTree(cds))$name == root_cell, arr.ind=T)
      graph_point_for_root_cell <- cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex[root_cell,]
      root_cell = V(minSpanningTree(cds))[graph_point_for_root_cell]$name
    }
    
  }else{
    if (is.null(minSpanningTree(cds))){
      stop("Error: no spanning tree found for CellDataSet object. Please call reduceDimension before calling orderCells()")
    }
    diameter <- get.diameter(minSpanningTree(cds))
    if (is.null(reverse) == FALSE && reverse == TRUE){
      root_cell = names(diameter[length(diameter)])
    } else {
      root_cell = names(diameter[1])
    }
  }
  return(root_cell)
}

extract_ddrtree_ordering <- function(cds, root_cell, verbose=T)
{
  
  dp <- cellPairwiseDistances(cds)
  dp_mst <- minSpanningTree(cds)
  
  curr_state <- 1
  
  res <- list(subtree = dp_mst, root = root_cell)
  
  states = rep(1, ncol(dp))
  names(states) <- V(dp_mst)$name
  
  pseudotimes = rep(0, ncol(dp))
  names(pseudotimes) <- V(dp_mst)$name
  
  parents = rep(NA, ncol(dp))
  names(parents) <- V(dp_mst)$name
  
  mst_traversal <- graph.dfs(dp_mst,
                             root=root_cell,
                             neimode = "all",
                             unreachable=FALSE,
                             father=TRUE)
  mst_traversal$father <- as.numeric(mst_traversal$father)
  curr_state <- 1
  
  for (i in 1:length(mst_traversal$order)){
    curr_node = mst_traversal$order[i]
    curr_node_name = V(dp_mst)[curr_node]$name
    
    if (is.na(mst_traversal$father[curr_node]) == FALSE){
      parent_node = mst_traversal$father[curr_node]
      parent_node_name = V(dp_mst)[parent_node]$name
      parent_node_pseudotime = pseudotimes[parent_node_name]
      parent_node_state = states[parent_node_name]
      curr_node_pseudotime = parent_node_pseudotime + dp[curr_node_name, parent_node_name]
      if (degree(dp_mst, v=parent_node_name) > 2){
        curr_state <- curr_state + 1
      }
    }else{
      parent_node = NA
      parent_node_name = NA
      curr_node_pseudotime = 0
    }
    
    curr_node_state = curr_state
    pseudotimes[curr_node_name] <- curr_node_pseudotime
    states[curr_node_name] <- curr_node_state
    parents[curr_node_name] <- parent_node_name
  }
  
  ordering_df <- data.frame(sample_name = names(states),
                            cell_state = factor(states),
                            pseudo_time = as.vector(pseudotimes),
                            parent = parents)
  row.names(ordering_df) <- ordering_df$sample_name
  # ordering_df <- plyr::arrange(ordering_df, pseudo_time)
  return(ordering_df)
}

project2MST <- function(cds, Projection_Method){
  dp_mst <- minSpanningTree(cds)
  Z <- reducedDimS(cds)
  Y <- reducedDimK(cds)
  
  cds <- findNearestPointOnMST(cds)
  closest_vertex <- cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex
  
  #closest_vertex <- as.vector(closest_vertex)
  closest_vertex_names <- colnames(Y)[closest_vertex]
  closest_vertex_df <- as.matrix(closest_vertex)
  row.names(closest_vertex_df) <- row.names(closest_vertex)
  #closest_vertex_names <- as.vector(closest_vertex)
  
  tip_leaves <- names(which(degree(dp_mst) == 1))
  
  if(!is.function(Projection_Method)) {
    P <- Y[, closest_vertex]
  }else{
    P <- matrix(rep(0, length(Z)), nrow = nrow(Z)) #Y
    for(i in 1:length(closest_vertex)) {
      neighbors <- names(V(dp_mst) [ suppressWarnings(nei(closest_vertex_names[i], mode="all")) ])
      projection <- NULL
      distance <- NULL
      Z_i <- Z[, i]
      
      for(neighbor in neighbors) {
        if(closest_vertex_names[i] %in% tip_leaves) {
          tmp <- projPointOnLine(Z_i, Y[, c(closest_vertex_names[i], neighbor)]) #projPointOnLine: always perform orthogonal projection to the line
        } else {
          tmp <- Projection_Method(Z_i, Y[, c(closest_vertex_names[i], neighbor)])
        }
        projection <- rbind(projection, tmp)
        distance <- c(distance, dist(rbind(Z_i, tmp)))
      }
      #   if(class(projection) != 'matrix')
      projection <- as.matrix(projection)
      P[, i] <- projection[which(distance == min(distance))[1], ] #use only the first index to avoid assignment error
    }
  }
  # tip_leaves <- names(which(degree(minSpanningTree(cds)) == 1))
  
  colnames(P) <- colnames(Z)
  mat1 <- apply(t(P), 1, crossprod)
  mat1 <- matrix(mat1, nrow=ncol(P), ncol=ncol(P))
  mat3 <- tcrossprod(t(P))
  mat3 <- 2*mat3
  mat1 = mat1+t(mat1)
  gc()
  mat4 = mat1-mat3
  rm(mat1,mat3)
  mat4[mat4<0] = 0
  diag(mat4) <- 0
  mat4 <- sqrt(mat4)
  dp = mat4
  rm(mat4)
  min_dist = min(dp[dp!=0])
  #dp[dp == 0] <- min_dist
  gc()
  dp <- dp + min_dist #to avoid exact Pseudotime for a lot cells
  diag(dp) <- 0
  
  cellPairwiseDistances(cds) <- dp
  gc()
  gp <- graph.adjacency(dp, mode = "undirected", weighted = TRUE)
  rm(dp)
  gc()
  dp_mst <- minimum.spanning.tree(gp)
  
  cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_tree <- dp_mst
  cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_dist <- P #dp, P projection point not output
  cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex <- closest_vertex_df #as.matrix(closest_vertex)
  
  cds
}

findNearestPointOnMST <- function(cds){
  dp_mst <- minSpanningTree(cds)
  Z <- reducedDimS(cds)
  Y <- reducedDimK(cds)
  
  tip_leaves <- names(which(degree(dp_mst) == 1))
  
  distances_Z_to_Y <- proxy::dist(t(Z), t(Y))
  closest_vertex <- apply(distances_Z_to_Y, 1, function(z) { which ( z == min(z) )[1] } )
  #closest_vertex <- which(distance_to_closest == min(distance_to_closest))
  
  #closest_vertex <- as.vector(closest_vertex)
  closest_vertex_names <- colnames(Y)[closest_vertex]
  closest_vertex_df <- as.matrix(closest_vertex) #index on Z
  row.names(closest_vertex_df) <- names(closest_vertex) #original cell names for projection
  
  cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex <- closest_vertex_df #as.matrix(closest_vertex)
  cds
}

project_point_to_line_segment <- function(p, df){
  # returns q the closest point to p on the line segment from A to B
  A <- df[, 1]
  B <- df[, 2]
  # vector from A to B
  AB <- (B-A)
  # squared distance from A to B
  AB_squared = sum(AB^2)
  if(AB_squared == 0) {
    # A and B are the same point
    q <- A
  }
  else {
    # vector from A to p
    Ap <- (p-A)
    # from http://stackoverflow.com/questions/849211/
    # Consider the line extending the segment, parameterized as A + t (B - A)
    # We find projection of point p onto the line.
    # It falls where t = [(p-A) . (B-A)] / |B-A|^2
    # t <- max(0, min(1, sum(Ap * AB) / AB_squared))
    t <- sum(Ap * AB) / AB_squared
    
    if (t < 0.0) {
      # "Before" A on the line, just return A
      q <- A
    }
    else if (t > 1.0) {
      # "After" B on the line, just return B
      q <- B
    }
    else {
      # projection lines "inbetween" A and B on the line
      q <- A + t * AB#
    }
  }
  return(q)
}

projPointOnLine <- function(point, line) {
  ap <- point - line[, 1]
  ab <- line[, 2] - line[, 1]
  
  res <- line[, 1] + c((ap %*% ab) / (ab %*% ab)) * ab
  return(res)
}

ordercell<- function(cds){
  library(igraph)
  root_state = NULL
  num_paths = NULL
  reverse = NULL
  root_cell <- select_root_cell(cds, root_state, reverse)
  cds@auxOrderingData <- new.env(hash = TRUE)
  cc_ordering <- extract_ddrtree_ordering(cds, root_cell)
  pData(cds)$Pseudotime <- cc_ordering[row.names(pData(cds)),]$pseudo_time
  K_old <- reducedDimK(cds)
  old_dp <- cellPairwiseDistances(cds)
  old_mst <- minSpanningTree(cds)
  old_A <- reducedDimA(cds)
  old_W <- reducedDimW(cds)
  dp_mst <- minSpanningTree(cds)
  Z <- reducedDimS(cds)
  Y <- reducedDimK(cds)
  Projection_Method = project_point_to_line_segment
  cds <- findNearestPointOnMST(cds)
  closest_vertex <- cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex
  
  #closest_vertex <- as.vector(closest_vertex)
  closest_vertex_names <- colnames(Y)[closest_vertex]
  closest_vertex_df <- as.matrix(closest_vertex)
  row.names(closest_vertex_df) <- row.names(closest_vertex)
  #closest_vertex_names <- as.vector(closest_vertex)
  
  tip_leaves <- names(which(degree(dp_mst) == 1))
  
  if(!is.function(Projection_Method)) {
    P <- Y[, closest_vertex]
  }else{
    P <- matrix(rep(0, length(Z)), nrow = nrow(Z)) #Y
    for(i in 1:length(closest_vertex)) {
      neighbors <- names(V(dp_mst) [ suppressWarnings(nei(closest_vertex_names[i], mode="all")) ])
      projection <- NULL
      distance <- NULL
      Z_i <- Z[, i]
      
      for(neighbor in neighbors) {
        if(closest_vertex_names[i] %in% tip_leaves) {
          tmp <- projPointOnLine(Z_i, Y[, c(closest_vertex_names[i], neighbor)]) #projPointOnLine: always perform orthogonal projection to the line
        } else {
          tmp <- Projection_Method(Z_i, Y[, c(closest_vertex_names[i], neighbor)])
        }
        projection <- rbind(projection, tmp)
        distance <- c(distance, dist(rbind(Z_i, tmp)))
      }
      #   if(class(projection) != 'matrix')
      projection <- as.matrix(projection)
      P[, i] <- projection[which(distance == min(distance))[1], ] #use only the first index to avoid assignment error
    }
  }
  # tip_leaves <- names(which(degree(minSpanningTree(cds)) == 1))
  
  colnames(P) <- colnames(Z)
  
  cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_tree <- dp_mst
  cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_dist <- P #dp, P projection point not output
  cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex <- closest_vertex_df #as.matrix(closest_vertex)
  minSpanningTree(cds) <- cds@auxOrderingData[[cds@dim_reduce_type]]$pr_graph_cell_proj_tree
  root_cell_idx <- which(V(old_mst)$name == root_cell, arr.ind = T)
  cells_mapped_to_graph_root <- which(cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex == root_cell_idx)
  if (length(cells_mapped_to_graph_root) == 0) {
    cells_mapped_to_graph_root <- root_cell_idx
  }
  cells_mapped_to_graph_root <- V(minSpanningTree(cds))[root_cell_idx]$name
  tip_leaves <- names(which(degree(minSpanningTree(cds)) ==  1))
  root_cell <- cells_mapped_to_graph_root[cells_mapped_to_graph_root %in% tip_leaves][1]
  if (is.na(root_cell)) {
    root_cell <- select_root_cell(cds, root_state, reverse)
  }
  cds@auxOrderingData[[cds@dim_reduce_type]]$root_cell <- root_cell
  cc_ordering_new_pseudotime <- extract_ddrtree_ordering(cds, root_cell)
  pData(cds)$Pseudotime <- cc_ordering_new_pseudotime[row.names(pData(cds)),  ]$pseudo_time
  if (is.null(root_state) == TRUE) {
    closest_vertex <- cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex
    pData(cds)$State <- cc_ordering[closest_vertex[, 1], ]$cell_state
  }
  mst_branch_nodes <- V(minSpanningTree(cds))[which(degree(minSpanningTree(cds)) > 2)]$name
  
  cds@auxOrderingData[[cds@dim_reduce_type]]$branch_points <- mst_branch_nodes
  cds
  
}