# # library(PERMANOVA)
# # library(MetStaT)
# # library(gASCA)
# 
# run_PE_ASCA = function(datasets, datatypes, design, alphaEstimation_numComponentsRange, alphaEstimation_CVfolds, alpha_tol=1e-4, nTries=100, penalty="L2", fun_concave="gdp", tol=1e-6, maxit=1000, gamma=1, numComponents=1, rand_start=0, thr_path=0, quiet=1, pvalueThreshold=0.05){
#   
#   # Run PESCA over all data blocks to find common and distinct variation
#   pesca_output = runPESCA(datasets, datatypes, alphaEstimation_numComponentsRange, alphaEstimation_CVfolds, alpha_tol, nTries, penalty, fun_concave, tol, maxit, gamma, numComponents, rand_start, thr_path, quiet)
#   
#   # Reinflate the PESCA output back into datasets
#   preparedPESCAoutput = createCommonDistinctDatasets(datasets, pesca_output, design)
#   commonBlocks = preparedPESCAoutput$common
#   distinctBlocks = preparedPESCAoutput$distinct
#   localcommonBlocks = preparedPESCAoutput$localcommon
#   
#   # Run ASCA on common datasets
#   commonModels = list()
#   if(length(commonBlocks) > 0){
#     for(i in 1:length(commonBlocks)){
#       X = preparedPESCAoutput$common[[i]]
#       importantFactors = findImportantFactorsPERMANOVA(X, design)
#       
#       # Option A: run ASCA through MetStaT implementation
#       commonModels[[i]] = ASCA.Calculate(X, design, paste(which(importantFactors<=pvalueThreshold), collapse=","))
#       
#       # Option B: run ASCA through gASCA implementation
#       #f = buildFormula(importantFactors, pvalueThreshold, colnames(design))
#       #print(f)
#       #commonModels[[i]] = ASCA_decompose(apply(design,2,as.factor), X, f)
#       
#       # Option C: run ASCA through multiblock implementation
#       # This approach does not work with the current architecture of the script
#       #f = buildFormula2(importantFactors, pvalueThreshold, colnames(design))
#       #df = prepASCAdata(design, data=X, names=colnames(design))
#       #print(f)
#       #h = as.formula(f)
#       #print(h)
#       #commonModels[[i]] = asca(h, df) 
#     }
#   }
#   
#   # Run ASCA on local common datasets
#   localcommonModels = list()
#   if(length(localcommonBlocks) > 0){
#     for(i in 1:length(localcommonBlocks)){
#       X = preparedPESCAoutput$localcommon[[i]]
#       importantFactors = findImportantFactorsPERMANOVA(X, design)
#       
#       # Option A: run ASCA through MetStaT implementation
#       localcommonModels[[i]] = ASCA.Calculate(X, design, paste(which(importantFactors<=pvalueThreshold), collapse=","))
#       
#       # Option B: run ASCA through gASCA implementation
#       #f = buildFormula(importantFactors, pvalueThreshold, colnames(design))
#       #print(f)
#       #localcommonModels[[i]] = ASCA_decompose(apply(design,2,as.factor), X, f)
#       
#       # Option C: run ASCA through multiblock implementation
#       # This approach does not work with the current architecture of the script
#       #f = buildFormula2(importantFactors, pvalueThreshold, colnames(design))
#       #df = prepASCAdata(design, data=X, names=colnames(design))
#       #print(f)
#       #localcommonModels[[i]] = asca(as.formula(f), df) 
#     }
#   }
#   
#   # Run ASCA on distinct datasets
#   distinctModels = list()
#   if(length(distinctBlocks) > 0){
#     for(i in 1:length(distinctBlocks)){
#       X = preparedPESCAoutput$distinct[[i]]
#       importantFactors = findImportantFactorsPERMANOVA(X, design)
#       
#       # Option A: run ASCA through MetStaT implementation
#       distinctModels[[i]] = ASCA.Calculate(X, design, paste(which(importantFactors<=pvalueThreshold), collapse=","))
#       
#       # Option B: run ASCA through gASCA implementation
#       #f = buildFormula(importantFactors, pvalueThreshold, colnames(design))
#       #distinctModels[[i]] = ASCA_decompose(apply(design,2,as.factor), X, f)
#       
#       # Option C: run ASCA through multiblock implementation
#       # This approach does not work with the current architecture of the script
#       #f = buildFormula2(importantFactors, pvalueThreshold, colnames(design))
#       #df = prepASCAdata(design, data=X, names=colnames(design))
#       #print(f)
#       #localcommonModels[[i]] = asca(as.formula(f), df) 
#     }
#   }
#   
#   return(list("common"=commonModels, "localcommon"=localcommonModels, "distinct"=distinctModels, "pesca_output"=pesca_output, "pesca_blocks"=preparedPESCAoutput))
# }
# 
# runPESCA = function(datasets, datatypes, alphaEstimation_numComponentsRange, alphaEstimation_CVfolds, alpha_tol=1e-4, nTries=100, penalty="L2", fun_concave="gdp", tol=1e-6, maxit=1000, gamma=1, numComponents=1, rand_start=0, thr_path=0, quiet=1){
#   
#   # Estimate alphas
#   alpha_est_output = pesca_estimate_alphas(datasets, alphaEstimation_numComponentsRange, alphaEstimation_CVfolds, tol=alpha_tol, quiet=quiet)
#   alphas = alpha_est_output[[1]]
#   
#   # parameters of a pESCA with concave L2norm penalty model
#   opts = list()
#   opts$gamma      = gamma             # hyper-parameter for the used penalty
#   opts$rand_start = rand_start        # initialization method: 1=random,0=SCA
#   opts$tol_obj    = tol               # stopping criterion
#   opts$maxit      = maxit             # maximum number of iterations
#   opts$alphas     = alphas            # alpha estimation from pesca_estimate_alphas
#   opts$R          = numComponents     # components used
#   opts$thr_path   = thr_path          # generate thresholding path or not
#   opts$quiet      = quiet             # don't show progress while running algorithm
#   
#   # Estimate lambdas
#   lambda_est_output = pesca_estimate_lambdas(datasets, datatypes, alphas, opts)
#   lambdas = lambda_est_output[[1]]
#   
#   # Create PESCA model
#   pesca_model = pESCA(dataSets=datasets, dataTypes=datatypes, lambdas, penalty=penalty, fun_concave=fun_concave, opts=opts)
#   
#   # Return all outputs
#   return(list("model"=pesca_model, "alpha"=alpha_est_output, "lambda"=lambda_est_output))
# }
# 
pesca_estimate_alphas = function(datasets, numComponentsRange, CVfolds, tol=1e-6, quiet=1){
  alphas_mean = rep(NA, length(datasets))
  alphas_std = rep(NA, length(datasets))
  R_selected = list()
  cvErrors = list()

  opts = list()
  opts$tol_obj = tol
  opts$quiet = quiet

  for (i in 1:length(datasets)){
    alpha_est = alpha_estimation(datasets[[i]], K = CVfolds, Rs = numComponentsRange, opts=opts) #5:20 NEED TO THINK ABOUT THIS CHANGE K
    alphas_mean[i] = alpha_est$alphas_mean
    alphas_std[i] = alpha_est$alphas_std
    R_selected[[i]] = alpha_est$R_CV
    cvErrors[[i]] = alpha_est$cvErrors
  }
  names(alphas_mean) = paste0('alpha_mean_', 1:length(alphas_mean))
  names(alphas_std) = paste0('alpha_std_', 1:length(alphas_std))

  return(list(alphas_mean, alphas_std, R_selected, cvErrors))
}
# 
# pesca_estimate_lambdas = function(datasets, datatypes, alphas, opts, nTries=100, penalty="L2", fun_concave="gdp"){
#   
#   # Initialize a sequence of lambdas to try
#   lambdas_CV = log10_seq(from=1, to=500, length.out=nTries)
#   
#   # Run cross-validation
#   result_CV = pESCA_CV(datasets, datatypes, lambdas_CV, penalty, fun_concave, opts=opts)
#   cvErrors_mat = result_CV$cvErrors_mat
#   inits = result_CV$inits
#   outs = result_CV$outs
#   
#   # Find best lambda estimation
#   index = which.min(result_CV$cvErrors_mat[,1])
#   lambdas = rep(lambdas_CV[index], length(datasets))
#   
#   return(list(lambdas, lambdas_CV, cvErrors_mat, inits, outs))
# }
# 
# findCommonLocalDistinctVariation = function(datasets, pesca_output, minVariation=1){
#   numBlocks = length(datasets)
#   
#   varExps = pesca_output$model$varExpPCs
#   varExps = varExps[1:(nrow(varExps)-1),] # remove total variation row
#   
#   mappingComponentsToBlocks = varExps > minVariation
#   explainedNumBlocksPerComponent = colSums(mappingComponentsToBlocks)
#   
#   common_components = which(explainedNumBlocksPerComponent == numBlocks)
#   localcommon_components = which(explainedNumBlocksPerComponent > 1 & explainedNumBlocksPerComponent < numBlocks)
#   distinct_components = which(explainedNumBlocksPerComponent == 1)
#   
#   return(list("mapping"=mappingComponentsToBlocks,"common"=common_components,"distinct"=distinct_components,"localcommon"=localcommon_components))
# }
# 
# 
# findCommonLocalDistinctVariation_dev = function(datasets, pesca_output, minVariation=1){
#   numBlocks = length(datasets)
#   
#   varExps = pesca_output$outcome$varExpPCs
#   varExps = varExps[1:(nrow(varExps)-1),] # remove total variation row
#   
#   mappingComponentsToBlocks = varExps > minVariation
#   explainedNumBlocksPerComponent = colSums(mappingComponentsToBlocks)
#   
#   common_components = which(explainedNumBlocksPerComponent == numBlocks)
#   localcommon_components = which(explainedNumBlocksPerComponent > 1 & explainedNumBlocksPerComponent < numBlocks)
#   distinct_components = which(explainedNumBlocksPerComponent == 1)
#   
#   return(list("mapping"=mappingComponentsToBlocks,"common"=common_components,"distinct"=distinct_components,"localcommon"=localcommon_components))
# }

#cutting out lines from sim
generate_options <- function(gamma, alpha_est) {
  opts <- list(
    gamma = gamma,
    rand_start = 0,
    tol_obj = 1e-6,
    maxit = 1000,
    alphas = alpha_est,
    R = 6,
    thr_path = 0,
    quiet = 1
  )
  return(opts)
}

reset_A0_B0 <- function(opts) {
  opts$A0 <- NULL
  opts$B0 <- NULL
  return(opts)
}

pesca_FS <- function(PESCAmodel, Data, y, spikes, type = c("yasblock","PCR")){
  #used to test performance. 
  #for actual use - CHANGE -
  
  stp <- t(PESCAmodel$outcome$A) %*% y
  
  if(type == "PCR"){
    scaled_lods <- PESCAmodel$outcome$B[] %*% stp
    rownames(scaled_lods) <- colnames(do.call(cbind,Data))[]
    
    
  }else{
    scaled_lods <- PESCAmodel$outcome$B[-nrow(PESCAmodel$outcome$B),] %*% stp
    rownames(scaled_lods) <- colnames(do.call(cbind,Data))[-nrow(PESCAmodel$outcome$B)]
    
  }
  
  
  
  cands <- scaled_lods[order(abs(scaled_lods[,1]), decreasing = T),]
  
  
  RP <- prod(which(names(cands) %in% spikes))^(1/length(spikes))
  
  return(list(RP, scaled_lods))
}




pesca_score_comp <- function(input,output){
  
  #input = simulation
  #output = PESCARmod
  
  input <- input[[5]]
  output <- output$outcome
  
  order_A <- t(output$A) %*% input
  inds <- max.col(abs(order_A))    # - WARNING - ; needs confirmation that for each 1:R the values are only used once - not always true
  
  ordered <- output$A[,order(inds)]
  
  
  image(t(input) %*% ordered)
  
  return(output$varExpPCs[,order(inds)])
  
}


# pesca_score_comp(simulation[[5]],PESCARmod$outcome$A)


RV_coef <- function(X,Y){
  
  num <- sum(diag(X %*% t(X) %*% Y %*% t(Y)))
  denom <- sqrt(sum(diag((X %*% t(X))^2)) * sum(diag((Y %*% t(Y))^2)))
  
  return(num/denom)
  
}



spikes_from_sim <- function(simulation){
  #get spikes from sim
  
  yWeights <- simulation[[4]][c(2,1,3)]      # - WARNING -
  loadings <- simulation[[2]]
  
  yWeights <- c(yWeights,rep(0,ncol(loadings) - length(yWeights)))
  
  
  cands <- loadings %*% yWeights
  
  
  cands <- cands[order(abs(cands[,1]), decreasing = T),]
  
  spikes <- cands[which(cands != 0)]
  return(spikes)
}


F_scale <- function(x){
  y <- x/norm(x, type = "F")
}


# createCommonDistinctDatasets = function(datasets, pesca_output, design){
#   scores = pesca_output$model$A
#   loadings = pesca_output$model$B
#   
#   variation = findCommonLocalDistinctVariation(datasetsPareto, pesca_output)
#   common_PCs = variation$common
#   distinct_PCs = variation$distinct
#   localcommon_PCs = variation$localcommon
#   
#   output = list()
#   output$variation = variation
#   output$common = list()
#   output$distinct = list()
#   output$local_common = list()
#   
#   numFeatures = ncol(datasets[[1]])
#   
#   for(i in 1:ncol(variation$mapping)){
#     if(sum(variation$mapping[,i]) == length(datasets)){
#       compType = "common"
#     }
#     else if (sum(variation$mapping[,i]) == 1){
#       compType = "distinct"
#     }
#     else{
#       compType = "localcommon"
#     }
#     
#     dataBlocks = which(variation$mapping[,i])
#     for(j in 1:length(dataBlocks)){
#       dataBlockNum = dataBlocks[j]
#       indexStart = numFeatures*(dataBlockNum-1) + 1
#       indexEnd = numFeatures*dataBlockNum
#       output[[compType]][[dataBlockNum]] = scores[,i] %*% t(loadings[indexStart:indexEnd, i])
#     }
#   }
#   
#   return(output)
# }
# 
# findImportantFactorsPERMANOVA = function(data, design){
#   result = 1:ncol(design)
#   data = IniTransform(data)
#   D = DistContinuous(data)
#   
#   for(i in 1:ncol(design)){
#     result[i] = PERMANOVA(D, as.factor(design[,i]))$pvalue
#   }
#   
#   return(result)
# }
# 
# buildFormula = function(PERMANOVAoutput, pvalueThreshold, factorNames){
#   # Just joins the terms together into a string
#   result = ""
#   terms = factorNames[PERMANOVAoutput <= pvalueThreshold]
#   for(i in 1:length(terms)){
#     if(i >= 2){
#       result = paste0(result, " + ")
#     }
#     result = paste0(result, terms[i])
#   }
#   #result = reformulate(termlabels=terms, response="data", env=env)
#   return(result)
# }
# 
# buildFormula2 = function(PERMANOVAoutput, pvalueThreshold, factorNames){
#   # Builds an actual formula
#   result = "data ~ "
#   terms = factorNames[PERMANOVAoutput <= pvalueThreshold]
#   for(i in 1:length(terms)){
#     if(i >= 2){
#       result = paste0(result, " + ")
#     }
#     result = paste0(result, terms[i])
#   }
#   #result = reformulate(termlabels=terms, response="data", env=env)
#   return(result)
# }