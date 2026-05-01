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
