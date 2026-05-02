
set.seed(31)          #21


library(tidygraph)
library(igraph)
library(ggraph)
library(ggplot2)
library(ggforce)


make_indexed_feature_block_map <- function(Data, y_rows = 0L, prefix = "X"){
  stopifnot(is.list(Data), !is.null(names(Data)))
  
  x_map <- unlist(lapply(seq_along(Data), function(i) {
    p <- ncol(Data[[i]])
    ids <- paste0(prefix, i, "_", seq_len(p))
    stats::setNames(rep(names(Data)[i], p), ids)
  }), use.names = TRUE)
  
  if (y_rows > 0L) {
    y_map <- stats::setNames(rep("Y", y_rows), paste0("Y_", seq_len(y_rows)))
    x_map <- c(x_map, y_map)
  }
  
  x_map
}

make_feature_block_map <- function(Data, indexed_ids = FALSE, prefix = "X"){
  stopifnot(is.list(Data), !is.null(names(Data)))
  
  if (!indexed_ids) {
    all_feats <- unlist(lapply(names(Data), function(b) {
      feats <- colnames(Data[[b]])
      stats::setNames(rep(b, length(feats)), feats)
    }), use.names = TRUE)
    
    # if a feature appears in multiple blocks, keep the first
    return(all_feats[!duplicated(names(all_feats))])
  }
  
  out <- unlist(lapply(seq_along(Data), function(i) {
    b <- names(Data)[i]
    p <- ncol(Data[[i]])
    ids <- paste0(prefix, i, "_", seq_len(p))
    stats::setNames(rep(b, p), ids)
  }), use.names = TRUE)
  
  out
}


compress_entities_from_Btrue <- function(B_true, feature_block = NULL, Data = NULL,
                                         use_sign = FALSE, prefix = "X"){
  if (is.null(colnames(B_true))) colnames(B_true) <- as.character(seq_len(ncol(B_true)))
  feat_ids <- rownames(B_true)
  
  # get block for each feature
  if (is.null(feature_block)) {
    if (is.null(Data)) {
      stop("Provide either feature_block (named vector) or Data (named list) to infer blocks.")
    }
    
    p_x <- sum(vapply(Data, ncol, integer(1)))
    n_all <- nrow(B_true)
    
    if (n_all < p_x) {
      stop("B_true has fewer rows than the total number of columns in Data.")
    }
    
    n_y <- n_all - p_x
    
    # first try the original-name logic
    feature_block_named <- make_feature_block_map(Data, indexed_ids = FALSE)
    
    can_match_names <- !is.null(feat_ids) &&
      length(feat_ids) == n_all &&
      all(!is.na(feature_block_named[feat_ids[seq_len(p_x)]]))
    
    if (can_match_names) {
      feature_block <- feature_block_named
      
      if (n_y > 0L) {
        y_ids <- feat_ids[p_x + seq_len(n_y)]
        feature_block <- c(feature_block, stats::setNames(rep("Y", n_y), y_ids))
      }
    } else {
      # fallback: assign ordered generic IDs from Data
      feature_block <- make_feature_block_map(Data, indexed_ids = TRUE, prefix = prefix)
      y_ids <- if (n_y > 0L) paste0("Y_", seq_len(n_y)) else character(0)
      
      rownames(B_true) <- c(names(feature_block), y_ids)
      feat_ids <- rownames(B_true)
      
      if (n_y > 0L) {
        feature_block <- c(feature_block, stats::setNames(rep("Y", n_y), y_ids))
      }
    }
  }
  
  feat_ids <- rownames(B_true)
  if (is.null(feat_ids)) stop("B_true must have rownames = feature IDs.")
  
  block_vec <- unname(feature_block[feat_ids])
  if (any(is.na(block_vec))) {
    missing_ids <- feat_ids[is.na(block_vec)]
    stop(paste0("Missing block assignment for ", length(missing_ids),
                " features (e.g. ", missing_ids[1], ")."))
  }
  
  # signature per feature
  if (!use_sign) {
    sig <- apply(B_true != 0, 1, function(v){
      idx <- which(v)
      if (length(idx) == 0) return("none")
      paste0("c", idx, collapse = "_")
    })
  } else {
    sig <- apply(B_true, 1, function(v){
      idx <- which(v != 0)
      if (length(idx) == 0) return("none")
      paste0("c", idx, ":", ifelse(v[idx] > 0, "+", "-"), collapse = "_")
    })
  }
  
  group_key <- paste(block_vec, sig, sep = " | ")
  
  nodes_ent <- data.frame(
    entity = group_key,
    block = block_vec,
    signature = sig,
    stringsAsFactors = FALSE
  )
  nodes_ent <- unique(nodes_ent)
  
  n_feat <- as.data.frame(table(group_key), stringsAsFactors = FALSE)
  names(n_feat) <- c("entity", "n_features")
  nodes_ent <- merge(nodes_ent, n_feat, by = "entity", all.x = TRUE)
  
  feature_to_entity <- data.frame(
    feature = feat_ids,
    entity = group_key,
    block = block_vec,
    signature = sig,
    stringsAsFactors = FALSE
  )
  
  comps <- seq_len(ncol(B_true))
  hyperedges <- setNames(vector("list", length(comps)), colnames(B_true))
  
  mem_counts <- vector("list", length(comps))
  names(mem_counts) <- colnames(B_true)
  
  for (j in comps) {
    feat_in <- feat_ids[B_true[, j] != 0]
    ent_in <- feature_to_entity$entity[match(feat_in, feature_to_entity$feature)]
    ent_in <- ent_in[!is.na(ent_in)]
    hyperedges[[j]] <- unique(ent_in)
    
    tab <- as.data.frame(table(ent_in), stringsAsFactors = FALSE)
    names(tab) <- c("entity", "n_features_in_hyperedge")
    tab$hyperedge <- colnames(B_true)[j]
    mem_counts[[j]] <- tab
  }
  
  mem_counts <- do.call(rbind, mem_counts)
  
  list(
    B_true = B_true,
    feature_block = feature_block,
    hyperedges = hyperedges,
    entity_nodes = nodes_ent,
    feature_to_entity = feature_to_entity,
    membership_counts = mem_counts
  )
}






make_B_plot <- function(B_true,
                        Data = NULL,
                        feature_block = NULL,
                        tau_rel = 0.15,
                        min_components = 1L,
                        keep_y = TRUE,
                        drop_rows = TRUE,
                        prefix = "X") {
  
  B_plot <- as.matrix(B_true)
  
  if (is.null(colnames(B_plot))) {
    colnames(B_plot) <- paste0("c", seq_len(ncol(B_plot)))
  }
  
  # build stable feature map BEFORE filtering
  if (is.null(feature_block)) {
    if (is.null(Data)) stop("Provide either Data or feature_block.")
    
    n_x <- sum(vapply(Data, ncol, integer(1)))
    n_all <- nrow(B_plot)
    
    if (n_all < n_x) {
      stop("B_true has fewer rows than the total number of columns in Data.")
    }
    
    n_y <- n_all - n_x
    feature_block <- make_indexed_feature_block_map(Data, y_rows = n_y, prefix = prefix)
  }
  
  if (length(feature_block) != nrow(B_plot)) {
    stop("length(feature_block) must equal nrow(B_true).")
  }
  if (is.null(names(feature_block))) {
    stop("feature_block must be a named vector.")
  }
  
  # force stable rownames from the feature map
  rownames(B_plot) <- names(feature_block)
  
  y_rows <- unname(feature_block[rownames(B_plot)]) == "Y"
  x_idx <- which(!y_rows)
  y_idx <- which(y_rows)
  
  # sparsify predictor rows component-wise using relative threshold only
  for (j in seq_len(ncol(B_plot))) {
    abs_vals <- abs(B_plot[x_idx, j])
    mx <- max(abs_vals, na.rm = TRUE)
    
    if (is.finite(mx) && mx > 0) {
      thr <- tau_rel * mx
      B_plot[x_idx[abs_vals < thr], j] <- 0
    } else {
      B_plot[x_idx, j] <- 0
    }
  }
  
  # remove rows that are no longer active
  row_nnz <- rowSums(B_plot != 0)
  keep_rows <- row_nnz >= min_components
  
  if (keep_y && length(y_idx) > 0) {
    keep_rows[y_idx] <- TRUE
  }
  
  if (drop_rows) {
    B_plot <- B_plot[keep_rows, , drop = FALSE]
  }
  
  # carry the matching block map along with the filtered matrix
  attr(B_plot, "feature_block") <- feature_block[rownames(B_plot)]
  
  B_plot
}






# --- clique expansion helper 
hyper_to_entity_edges <- function(hyperedges) {
  edge_rows <- do.call(
    rbind,
    lapply(names(hyperedges), function(h) {
      members <- unique(hyperedges[[h]])
      if (length(members) < 2L) return(NULL)
      cmb <- t(combn(members, 2))
      data.frame(from = cmb[, 1], to = cmb[, 2], hyperedge = h, stringsAsFactors = FALSE)
    })
  )
  if (is.null(edge_rows) || nrow(edge_rows) == 0L) {
    return(data.frame(from = character(), to = character(), weight = integer()))
  }
  edges_w <- aggregate(hyperedge ~ from + to, edge_rows, length)
  names(edges_w)[names(edges_w) == "hyperedge"] <- "weight"
  edges_w
}

make_shades <- function(base_col, n) {
  if (n <= 1) return(base_col)
  grDevices::colorRampPalette(c(base_col, "white"))(n)  # lighter variants
}


###########




