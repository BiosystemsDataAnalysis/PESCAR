generateDesignFactor = function(numFactors, lengthFactors, totalLength){
  numReps = totalLength / sum(lengthFactors)
  
  setup = vector()
  for(i in 1:numFactors){
    setup = c(setup, rep(i, lengthFactors[i]))
  }
  
  result = rep(setup, numReps)
  return(result)
}

generateScoresFromDesign = function(design, formula){
  df = prepASCAdata(design)
  asca_model = asca(formula = formula, data=df)
  return(asca_model)
}

prepASCAdata = function(design, data=NA, names=NA){
  numFactors = ncol(design)
  
  # If no data is supplied, make a random dataset
  if(!is.matrix(data)){
    n = nrow(design)
    m = 100 # m sets the number of columns for data generation
    data = matrix(rnorm(n*m), nrow=n)
  }
  df = data.frame("data" = I(data))
  
  # Add all the factors to the dataframe
  for(i in 1:numFactors){
    if(is.character(names)){
      df[names[i]] = as.factor(design[,i])
    }
    else{
      df[LETTERS[i]] = as.factor(design[,i])
    }
  }
  
  return(df)
}

createBinaryLoading = function(nonzeroIndices, length=100){
  result = rep(0, length)
  result[nonzeroIndices] = 1
  return(result)
}

createNoise = function(n, m, mean=0, sd=1, centered=FALSE, sumsqr=0){
  E = matrix(rnorm(n*m, mean=mean, sd=sd), nrow=n, ncol=m)
  
  # Center the noise if desired
  if(centered == TRUE){
    E = sweep(E, 2, colMeans(E), FUN="-")
  }
  
  # Scale the block to a specific sum of squares if desired
  if(sumsqr != 0){
    v = sum(E^2)
    E = E / sqrt(v) * sqrt(sumsqr)
  }
  
  return(E)
}




####PLOTTING


# Scale theme text + any explicit geom_text/label sizes inside the plot object
bump_plot_text <- function(p, base_size = 18, mult = 1.4) {
  p <- p + theme(
    text        = element_text(size = base_size),
    axis.title  = element_text(size = base_size),
    axis.text   = element_text(size = round(base_size * 0.9)),
    legend.title= element_text(size = base_size),
    legend.text = element_text(size = round(base_size * 0.9))
  )
  
  for (i in seq_along(p$layers)) {
    geom_i <- p$layers[[i]]$geom
    is_text_geom <- inherits(geom_i, "GeomText") ||
      inherits(geom_i, "GeomLabel") ||
      inherits(geom_i, "GeomTextRepel") ||
      inherits(geom_i, "GeomLabelRepel")
    
    if (is_text_geom) {
      if (!is.null(p$layers[[i]]$aes_params$size)) {
        p$layers[[i]]$aes_params$size <- p$layers[[i]]$aes_params$size * mult
      }
      if (!is.null(p$layers[[i]]$geom_params$size)) {
        p$layers[[i]]$geom_params$size <- p$layers[[i]]$geom_params$size * mult
      }
    }
  }
  p
}

# helper to add a panel label (allow fontsize control)
label_panel <- function(grob, lab, fontsize = 28) {
  gridExtra::arrangeGrob(
    grob,
    top = grid::textGrob(lab, x = grid::unit(0, "npc"), just = "left",
                         gp = grid::gpar(fontface = "bold", fontsize = fontsize))
  )
}

