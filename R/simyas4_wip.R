


numFeatures_chems = c(10,0,2,0,10)
numFeatures_genes = c(5,0,2,0,0)
numFeatures_bact = c(0,0,2,0,5)


# library(igraph)

#------------------------------------------------------------------
# Helper: generate a simple example component graph (if not supplied)
#------------------------------------------------------------------
generate_component_graph <- function(num_components = 6, edge_prob = 0.3, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  g <- erdos.renyi.game(
    n        = num_components,
    p.or.m   = edge_prob,
    type     = "gnp",
    directed = FALSE
  )
  V(g)$name <- as.character(seq_len(vcount(g)))
  g
}

#------------------------------------------------------------------
# Helper: generate block-specific loading matrix from a subgraph
# Each feature either:
#   - sits on an edge (loads on both incident components), or
#   - sits on a single node (pure component feature).
#------------------------------------------------------------------
generate_block_loadings_from_topology <- function(
    graph_block,
    comps_active,
    num_features,
    distribution   = rnorm,
    dist_params    = list(mean = 1, sd = 0.1),
    pure_frac      = 0.3  # fraction of features that are "pure" (single-component)
) {
  k_active <- length(comps_active)
  if (k_active == 0L || num_features == 0L) {
    return(matrix(0, nrow = num_features, ncol = 0))
  }
  
  # number of pure vs edge-based features
  num_pure  <- round(num_features * pure_frac)
  num_edgef <- max(num_features - num_pure, 0L)
  
  # start with all zeros
  P_sub <- matrix(0, nrow = num_features, ncol = k_active)
  colnames(P_sub) <- as.character(comps_active)
  
  # --- Edge-based features (link components) ---
  edges_df <- as.data.frame(as_edgelist(graph_block))
  n_edges  <- nrow(edges_df)
  
  if (n_edges > 0L && num_edgef > 0L) {
    features_per_edge <- pmax(1L, floor(num_edgef / n_edges))
    edge_feature_idx  <- seq_len(num_edgef)
    
    # read edge "pattern" attribute (if present); default to "same"
    edge_pattern <- edge_attr(graph_block, "pattern")
    if (is.null(edge_pattern)) {
      edge_pattern <- rep("same", n_edges)
    }
    if (length(edge_pattern) != n_edges) {
      stop("Length of E(graph_block)$pattern must equal number of edges.")
    }
    
    for (i in seq_len(n_edges)) {
      edge <- edges_df[i, ]
      # map igraph vertex names (character) to column positions
      comp1 <- as.character(edge$V1)
      comp2 <- as.character(edge$V2)
      col1  <- match(comp1, colnames(P_sub))
      col2  <- match(comp2, colnames(P_sub))
      
      if (is.na(col1) || is.na(col2)) next
      
      idx_start <- (i - 1L) * features_per_edge + 1L
      idx_end   <- min(i * features_per_edge, num_edgef)
      if (idx_start > idx_end) next
      
      idx <- edge_feature_idx[idx_start:idx_end]
      n_edgefeat <- length(idx)
      
      # base random draw
      s_base <- do.call(
        distribution,
        c(list(n = n_edgefeat), dist_params)
      )
      
      # defaults: same sign, same magnitudes
      load1 <- s_base
      load2 <- s_base
      
      pat <- edge_pattern[i]
      if (is.na(pat)) pat <- "same"
      
      ## -------- pattern logic --------
      if (pat == "same") {
        # load1, load2 already same
        
      } else if (pat == "opposite") {
        # perfect negative correlation
        load2 <- -s_base
        
      } else if (pat == "dom12") {
        # comp1 dominant, comp2 weaker
        load2 <- 0.3 * s_base
        
      } else if (pat == "dom21") {
        # comp2 dominant, comp1 weaker
        load1 <- 0.3 * s_base
        
      } else if (pat == "one1") {
        # only comp1 nonzero
        load2 <- 0
        
      } else if (pat == "one2") {
        # only comp2 nonzero
        load1 <- 0
        
      } else if (pat == "corr_pos" || pat == "corr_neg") {
        # noisy (positively / negatively) correlated pair
        rho <- if (pat == "corr_pos") 0.7 else -0.7
        noise <- rnorm(n_edgefeat, mean = 0, sd = sd(s_base))
        load1 <- s_base
        load2 <- rho * s_base + sqrt(1 - rho^2) * noise
        
      } else if (pat == "circle") {
        # points on a circle (or near it) in 2D
        t <- seq(0, 2 * pi, length.out = n_edgefeat)
        r <- 1  # could also use r <- abs(s_base) / max(abs(s_base))
        load1 <- r * cos(t)
        load2 <- r * sin(t)
        
      } else if (pat == "spiral") {
        # spiral in 2D: radius grows with t
        t <- seq(0, 2 * pi, length.out = n_edgefeat)
        r <- seq(0.1, 1, length.out = n_edgefeat)
        load1 <- r * cos(t)
        load2 <- r * sin(t)
        
      } else {
        warning("Unknown edge pattern '", pat, "', using 'same'.")
      }
      
      P_sub[idx, col1] <- load1
      P_sub[idx, col2] <- load2
    }
  }
  
  # --- Pure features (single-component) ---
  if (num_pure > 0L) {
    pure_idx <- (num_edgef + 1L):num_features
    # distribute equally across active components
    features_per_comp <- pmax(1L, floor(length(pure_idx) / k_active))
    for (j in seq_len(k_active)) {
      idx_start <- (j - 1L) * features_per_comp + 1L
      idx_end   <- min(j * features_per_comp, length(pure_idx))
      if (idx_start > idx_end) next
      
      idx <- pure_idx[idx_start:idx_end]
      
      strength_vals <- do.call(
        distribution,
        c(list(n = length(idx)), dist_params)
      )
      
      P_sub[idx, j] <- strength_vals
    }
  }
  
  P_sub
}

build_CLD_loadings_from_graph <- function(
    graph,
    num_features_list,    # named vector: c(chems = 100, genes = 100, bact = 100)
    distribution   = rnorm,
    dist_params    = list(mean = 1, sd = 0.1),
    pure_frac      = 0.3,
    edge_fun       = NULL  # << NEW
) {
  K <- vcount(graph)
  block_names <- names(num_features_list)
  
  P_list <- vector("list", length(block_names))
  names(P_list) <- block_names
  
  for (b in block_names) {
    # components active in this block: those whose V(graph)$blocks contains b
    comps_active <- which(vapply(
      V(graph)$blocks,
      function(bl) b %in% bl,
      logical(1)
    ))
    
    if (length(comps_active) == 0L) {
      P_block <- matrix(
        0, nrow = num_features_list[[b]], ncol = K,
        dimnames = list(NULL, paste0("Comp_", seq_len(K)))
      )
      prefix <- substr(b, 1, 1)
      rownames(P_block) <- paste0(toupper(prefix), "_", seq_len(nrow(P_block)))
      P_list[[b]] <- P_block
      next
    }
    
    # Induced subgraph on active components
    g_sub <- induced_subgraph(graph, vids = comps_active)
    
    P_sub <- generate_block_loadings_from_topology(
      graph_block  = g_sub,
      comps_active = comps_active,
      num_features = num_features_list[[b]],
      distribution = distribution,
      dist_params  = dist_params,
      pure_frac    = pure_frac
    )
    
    # embed into full K-component matrix
    P_block <- matrix(
      0, nrow = num_features_list[[b]], ncol = K,
      dimnames = list(NULL, paste0("Comp_", seq_len(K)))
    )
    P_block[, comps_active] <- P_sub
    
    prefix <- substr(b, 1, 1)
    rownames(P_block) <- paste0(toupper(prefix), "_", seq_len(nrow(P_block)))
    
    P_list[[b]] <- P_block
  }
  
  P_list
}







dataSimWrapper = function(Yweights,
                          Ynoise_sd,
                          componentWeights = list(c(1,1,1,1), c(1,1,1,1), c(1,1,1)),
                          numFeatures_chems = c(10,0,2,0,10),
                          numFeatures_genes = c(5,0,2,0,0),
                          numFeatures_bact  = c(0,0,2,0,5),
                          chems_tot = 100,
                          bact_tot  = 100,
                          genes_tot = 100,
                          error_weights = c(1,1,1),
                          noise_sd = 0.75,
                          seed = 123,
                          polydemo = FALSE,
                          component_graph = NULL,      # optional CLD topology
                          topology_pure_frac = 0.3,    # control pure vs edge features
                          edge_fun = NULL              # << NEW
) {
  set.seed(seed)
  
  # Create study design
  litter = generateDesignFactor(numFactors=3, lengthFactors=c(18,18,18), totalLength=108)
  treatment = generateDesignFactor(numFactors=3, lengthFactors=c(6,6,6), totalLength=108)
  weather = generateDesignFactor(numFactors=2, lengthFactors=c(54,54), totalLength=108)
  individual = generateDesignFactor(numFactors=2, lengthFactors=c(2,4), totalLength=108)
  
  design = cbind(litter, treatment, weather, individual)
  
  
  
  orth_scores <- svd(scale(matrix(rnorm(300*108), ncol = 300)), 6)$u
  
  A <- orth_scores[,c(1:2)]
  B <- orth_scores[,c(3:4)]
  C <- orth_scores[,c(5), drop = F]
  D <- orth_scores[,c(6), drop = F]
  
  if(polydemo == T){
    
    n <- 6
    m <- 108
    
    # Create a sequence of m points over some [-1, 1]
    x <- seq(-1, 1, length.out = m)
    
    POLYS <- poly(x, degree = n, raw = FALSE)
    
    A <- POLYS[,c(1:2)]
    B <- POLYS[,c(3:4)]
    C <- POLYS[,c(5), drop = F]
    D <- POLYS[,c(6), drop = F]
    
  }
  
  
  t_ca = A[,1] #- mean(A[,1])
  t_cb = B[,2] #- mean(B[,2])
  t_lc1 = D[,1] #- mean(D[,1])
  t_lc2 = B[,1] #- mean(B[,1])
  t_2d = C[,1] #- mean(C[,1])
  t_3d = A[,2] #- mean(A[,2])
  
  
  num_features_list <- c(
    chems = chems_tot,
    genes = genes_tot,
    bact  = bact_tot
  )
  
  # If no graph provided, make a simple random one with 6 components 
  if (is.null(component_graph)) {
    component_graph <- generate_component_graph(
      num_components = 6,
      edge_prob      = 0.3,
      seed           = seed
    )
  }
  
  P_all <- build_CLD_loadings_from_graph(
    graph             = component_graph,
    num_features_list = num_features_list,
    distribution      = rnorm,
    dist_params       = list(mean = 1, sd = 0.1),
    pure_frac         = topology_pure_frac,
    edge_fun          = edge_fun          
  )
  
  P_chems <- P_all$chems
  P_genes <- P_all$genes
  P_bact  <- P_all$bact
  
  
  
  
  normalise_F <- function(mat) {
    n <- norm(mat, type = "F")
    if (n == 0) {
      matrix(0, nrow(mat), ncol(mat))
    } else {
      mat / n
    }
  }
  # browser()
  # Chems (X1)
  # All score vectors together (6 components)
  scores_mat <- cbind(t_ca, t_cb, t_lc1, t_lc2, t_2d, t_3d) # n x 6
  
  
  ## ---------------- NEW: topology-driven ---------------- ##
  ## We reuse the previous mapping of weights to scores:
  ##   chems: 1=t_ca,2=t_cb,3=t_lc1,4=t_lc2,6=t_3d
  ##   genes: 1=t_ca,2=t_cb,3=t_lc1,5=t_2d
  ##   bact : 1=t_ca,2=t_cb,3=t_lc1,4=t_lc2
  
  # safety: number of columns in P_* must match number of score components (K)
  K <- ncol(scores_mat)
  
  g <- component_graph  # igraph with V(g)$blocks
  
  build_block_X <- function(block_name, P_block, w_block) {
    # which components are present in this block?
    
    # browser()
    comps_active <- which(vapply(
      V(g)$blocks,
      function(bl) block_name %in% bl,
      logical(1)
    ))
    
    if (length(comps_active) == 0L) {
      return(matrix(0, nrow = nrow(scores_mat), ncol = nrow(P_block)))
    }
    
    # one weight per active component
    # stopifnot(length(w_block) == length(comps_active))
    
    X_block <- matrix(0, nrow = nrow(scores_mat), ncol = nrow(P_block))
    
    for (j in seq_along(comps_active)) {
      k <- comps_active[j]       # component index
      if (k > K) next            # safety
      X_block <- X_block +
        w_block[j] * normalise_F(scores_mat[, k, drop = FALSE] %*% t(P_block[, k]))
    }
    
    X_block
  }
  
  # weights per block: must be length = #active comps in that block
  #need to link to topology - WARNING -
  w_chems <- componentWeights[[1]]
  w_genes <- componentWeights[[2]]
  w_bact  <- componentWeights[[3]]
  
  X1 <- build_block_X("chems", P_chems, w_chems)
  X2 <- build_block_X("genes", P_genes, w_genes)
  X3 <- build_block_X("bact",  P_bact,  w_bact)
  
  X1clean <- X1
  X2clean <- X2
  X3clean <- X3
  
  
  
  
  E1 <- matrix(rnorm(dim(X1)[1] * dim(X1)[2], mean = 0, sd = noise_sd), nrow = nrow(X1))
  E1 <- error_weights[1] * (E1/norm(E1, type = "F"))
  
  E2 <- matrix(rnorm(dim(X2)[1] * dim(X2)[2], mean = 0, sd = noise_sd), nrow = nrow(X2))
  E2 <- error_weights[2] * (E2/norm(E2, type = "F"))
  
  E3 <- matrix(rnorm(dim(X3)[1] * dim(X3)[2], mean = 0, sd = noise_sd), nrow = nrow(X3))
  E3 <- error_weights[3] * (E3/norm(E3, type = "F"))
  
  
  # Add noise to the data
  X1 = X1 + E1
  X2 = X2 + E2
  X3 = X3 + E3
  
  X <- list(X1,X2,X3)
  X <- lapply(X,scale, scale = F)
  
  Data <- X
  names(Data) <- c("chems","genes","bact")
  
  # With topology, the CLD structure is encoded in P_* already.
  # We just stack them; number of columns = number of score components (K).
  l_c <- P_chems
  l_g <- P_genes
  l_b <- P_bact
  
  
  loadings <- rbind(l_c, l_g, l_b)
  
  
  
  scores = cbind(t_ca, t_cb, t_lc1, t_lc2, t_2d, t_3d)
  
  # browser()
  y = generateY(cbind(t_ca, t_cb, t_lc1, t_lc2), Yweights, Ynoise_sd)
  
  
  
  # browser()
  
  return(list(Data, loadings, y, Yweights, scores, Dataclean = cbind(X1clean,X2clean,X3clean),component_graph))   #, B_wt
}

generateY = function(scoreMatrix, weights, noise_sd=1){
  y = weights %*% t(scoreMatrix)
  y = y + rnorm(length(y), sd=noise_sd)
  return(t(y))
}

gram_schmidt <- function(V) {
  U <- V
  for (i in 2:nrow(V)) {
    for (j in 1:(i - 1)) {
      U[i, ] <- U[i, ] - (sum(U[i, ] * U[j, ]) / sum(U[j, ]^2)) * U[j, ]
    }
    U[i, ] <- U[i, ] / sqrt(sum(U[i, ]^2))
  }
  return(U)
}
