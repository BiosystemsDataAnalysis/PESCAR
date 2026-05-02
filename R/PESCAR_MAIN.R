
pESCA_CV_DEV_lambda_test <-function(dataSets, dataTypes, y,
                                    lambdas_CV=NULL, lambdas_CVy=NULL, penalty='L2',
                                    fun_concave='gdp', alpha, opts=list(), ORTH_A = TRUE){
  # check if the inputs satisfy the requirements
  stopifnot(class(dataSets) == "list")
  stopifnot(class(penalty) == "character")
  stopifnot(class(fun_concave) == "character")
  if(length(dataTypes)==1){dataTypes <- unlist(strsplit(dataTypes, split=""))}
  if(exists('quiet', where=opts)){quiet <- opts$quiet} else{quiet<-0};
  if(exists('thr_path', where=opts)){thr_path <- opts$thr_path} else{thr_path<-0};
  
  # number of datasets, size of each dataset
  nTries <- length(lambdas_CV)
  nTriesy <- length(lambdas_CVy)
  nDataSets <- length(dataSets) # number of data sets
  n <- rep(0,nDataSets)  # number of samples
  d <- rep(0,nDataSets)  # numbers of variables in different data sets
  for(i in 1:nDataSets){
    n[i] <- dim(dataSets[[i]])[1]
    d[i] <- dim(dataSets[[i]])[2]
  }
  if(length(unique(as.factor(n)))!=1)
    stop("multiple datasets have unequal sample size")
  n <- n[1]
  sumd <- sum(d) # total number of variables
  
  # default dispersion parameters alphas
  if(exists('alphas', where=opts)){alphas<-opts$alphas} else{alphas<-rep(1,nDataSets)};
  
  # create zero matrix to hold results
  cvErrors_mat <- matrix(data=0, # +1 is used for the sum of all the Xi
                         nrow=nTries*nTriesy, ncol=nDataSets+2)
  
  # model selection process
  opts_inner <- opts
  
  
  # split data sets into training set and test set
  splitedData <- dataSplit(dataSets=dataSets,       #CHANGE
                           dataTypes=dataTypes,
                           y = y,
                           ratio_mis=0.1)
  trainSets <- splitedData$trainSets
  
  
  # save the parameters during the model selection
  inits <- as.list(1:nTries)
  if(thr_path == 1){outs <- as.list(1:nTries)}
  
  l <- 1
  
  TrainModel <- list()
  
  # model selection
  for(j in 1:nTries){
    for(k in 1:nTriesy){
      
      lambda <- lambdas_CV[j]
      lambda_y <- lambdas_CVy[k]
      lambdas <- lambda*rep(1,nDataSets)
      
      # browser()
      
      
      trainModel <- pESCA_Ryas(dataSets = trainSets,
                               dataTypes = dataTypes,
                               y = splitedData$y_train,
                               lambdas = lambdas,
                               lambdas_y = lambda_y,
                               penalty=penalty,
                               fun_concave=fun_concave,
                               alpha = alpha,
                               opts=opts_inner)
      
      if((trainModel$iter <= 2) & (quiet==0)){
        print("less than 3 iteration is used.")
      }
      
      TrainModel[[l]] <- trainModel
      
      
      # - WARNING -
      
      # # warm start
      # mu <- trainModel$outcome$mu; A <- trainModel$outcome$A; B <- trainModel$outcome$B
      # opts_inner$mu0 <- mu
      # opts_inner$A0 <- A
      # opts_inner$B0 <- B
      # opts_inner$R <- dim(B)[2]
      
      inits[[l]] <- opts_inner
      if(thr_path == 1){outs[[l]] <- trainModel$Sigmas}
      
      yti <- splitedData$y_index_test
      
      # - BUG HERE? -
      
      # ytri <- splitedData$
      # compute the test error - for coefs;
      # coefs <- solve(t(A[,]) %*% A[,]) %*% t(A[,]) %*% y #splitedData$y_train
      
      
      # ThetaHat <- ones(n) %*% mu + A %*% t(B)
      # ThetaHat <- ones(n) %*% mu + A %*% matrix(diag(c(coefs)), ncol = length(coefs)) %*% t(B)
      
      # testError_vec <- alpha * cvError_comput(splitedData,dataTypes,alphas,ThetaHat,d)
      
      #add y rmse here and alpha weighting?
      # Prediction of y_hat
      # W_k <- ginv(ThetaHat) %*% A
      # y_hat <- predict_y_new(X_new = ThetaHat, W_k = W_k, regression_coefficients = coefs)
      # browser()
      
      # Calculate RMSE between y and y_hat
      # y_contribution <- (1 - alpha) * calculate_rmse(y[yti], y_hat[yti])
      
      
      # browser()
      # cvErrors_tmp <- c(sum(testError_vec), testError_vec, y_contribution)
      
      
      # cvErrors_mat[l,] <- cvErrors_tmp
      
      l <- l + 1
    }
    
  }
  
  # colnames(cvErrors_mat) <- c(paste0(rep("X_"), c("full", as.character(1:nDataSets))),"y")
  
  result_CV <- list()
  # result_CV$cvErrors_mat <- cvErrors_mat
  # result_CV$inits <- inits
  if(thr_path == 1){result_CV$outs <- outs}
  result_CV$TrainModel <- TrainModel
  
  return(result_CV)
}


pESCA_CV_DEV <-function(dataSets, dataTypes, y,
                        lambdas_CV=NULL, lambdas_CVy=NULL, penalty='L2',
                        fun_concave='gdp', alpha, opts=list(),ORTH_A = TRUE, rstartseed = 1){
  # check if the inputs satisfy the requirements
  stopifnot(class(dataSets) == "list")
  stopifnot(class(penalty) == "character")
  stopifnot(class(fun_concave) == "character")
  if(length(dataTypes)==1){dataTypes <- unlist(strsplit(dataTypes, split=""))}
  if(exists('quiet', where=opts)){quiet <- opts$quiet} else{quiet<-0};
  if(exists('thr_path', where=opts)){thr_path <- opts$thr_path} else{thr_path<-0};
  
  # number of datasets, size of each dataset
  nTries <- length(lambdas_CV)
  nTriesy <- length(lambdas_CVy)
  nDataSets <- length(dataSets) # number of data sets
  n <- rep(0,nDataSets)  # number of samples
  d <- rep(0,nDataSets)  # numbers of variables in different data sets
  for(i in 1:nDataSets){
    n[i] <- dim(dataSets[[i]])[1]
    d[i] <- dim(dataSets[[i]])[2]
  }
  if(length(unique(as.factor(n)))!=1)
    stop("multiple datasets have unequal sample size")
  n <- n[1]
  sumd <- sum(d) # total number of variables
  
  # default dispersion parameters alphas
  if(exists('alphas', where=opts)){alphas<-opts$alphas} else{alphas<-rep(1,nDataSets)};
  
  
  
  
  # create zero matrix to hold results
  cvErrors_mat <- matrix(data=0, # +1 is used for the sum of all the Xi
                         nrow=nTries*nTriesy, ncol=nDataSets+2)
  
  # model selection process
  opts_inner <- opts
  
  
  # split data sets into training set and test set
  splitedData <- dataSplit(dataSets=dataSets,       #CHANGE
                           dataTypes=dataTypes,
                           y = y,
                           ratio_mis=0.1)
  trainSets <- splitedData$trainSets
  
  
  # save the parameters during the model selection
  inits <- as.list(1:nTries)
  if(thr_path == 1){outs <- as.list(1:nTries)}
  
  l <- 1
  
  TrainModel <- list()
  
  # model selection
  for(j in 1:nTries){
    for(k in 1:nTriesy){
      
      lambda <- lambdas_CV[j]
      lambda_y <- lambdas_CVy[k]
      # using the training set to construct a PESCAR model
      lambdas <- lambda*rep(1,nDataSets)
      
      # browser()
      #generalise this for B
      if(is.list(rstartseed)) {opts_inner$A0 <- rstartseed[[1]]; 
      opts_inner$B0 <- rstartseed[[2]];
      opts_inner$mu0 <- rstartseed[[3]]}
      
      
      trainModel <- pESCA_Ryas(dataSets = trainSets,
                               dataTypes = dataTypes,
                               y = splitedData$y_train,
                               lambdas = lambdas,
                               lambdas_y = lambda_y,
                               penalty=penalty,
                               fun_concave=fun_concave,
                               alpha = alpha,
                               opts=opts_inner,
                               ORTH_A = ORTH_A,
                               rstartseed = rstartseed)
      
      if((trainModel$iter <= 2) & (quiet==0)){
        print("less than 3 iteration is used.")
      }
      
      # TrainModel[[l]] <- trainModel
      TrainModel[[l]] <- list(outcome = trainModel$outcome)
      
      # - WARNING -
      
      # warm start
      mu <- trainModel$outcome$mu; A <- trainModel$outcome$A; B <- trainModel$outcome$B
      # opts_inner$mu0 <- mu
      # opts_inner$A0 <- A
      # #
      # opts_inner$B0 <- B
      # opts_inner$R <- dim(B)[2]
      
      inits[[l]] <- opts_inner
      if(thr_path == 1){outs[[l]] <- trainModel$Sigmas}
      
      yti <- splitedData$y_index_test
      
      # - BUG HERE? -
      
      # ytri <- splitedData$
      # compute the test error - for coefs;
      coefs <- solve(t(A[,]) %*% A[,]) %*% t(A[,]) %*% y 
      
      
      ThetaHat <- RpESCA::ones(n) %*% mu + A %*% t(B)
      # ThetaHat <- ones(n) %*% mu + A %*% matrix(diag(c(coefs)), ncol = length(coefs)) %*% t(B)
      
      testError_vec <- alpha * cvError_comput(splitedData,dataTypes,alphas,ThetaHat,d)
      
      #add y rmse here and alpha weighting?
      # Prediction of y_hat
      W_k <- ginv(ThetaHat) %*% A
      y_hat <- predict_y_new(X_new = ThetaHat, W_k = W_k, regression_coefficients = coefs)
      # browser()
      
      # Calculate RMSE between y and y_hat
      y_contribution <- (1 - alpha) * calculate_rmse(y[yti], y_hat[yti])
      
      
      # browser()
      cvErrors_tmp <- c(sum(testError_vec), testError_vec, y_contribution)
      
      
      cvErrors_mat[l,] <- cvErrors_tmp
      
      l <- l + 1
    }
    
  }
  
  colnames(cvErrors_mat) <- c(paste0(rep("X_"), c("full", as.character(1:nDataSets))),"y")
  
  result_CV <- list()
  result_CV$cvErrors_mat <- cvErrors_mat
  result_CV$inits <- inits
  if(thr_path == 1){result_CV$outs <- outs}
  result_CV$TrainModel <- TrainModel
  
  return(result_CV)
}




initialisePESCA = function(dataSets, y, dataTypes, lambdas, lambdas_y, penalty="L2", opts=opts, response = NULL, alpha = alpha, rstartseed = rstartseed){
  #' Initialises many aspects to run PESCA
  #'
  #' @param dataSets A list of datasets.
  #' @param dataTypes Data type per dataset, string of letters: "G" = gaussian, "P" = poisson, "B" = binary.
  #' @param lambdas A list of fitted lambdas, one for each dataset.
  #' @param penalty Penalty type for the PESCA function. Default is L2.
  #'
  #' @returns A list of the following: A, B, mu, Theta, X, W, loss, rho, and some useful variables.
  # browser()
  numDatasets = length(dataSets)
  numSamplesPerDataset = unlist(lapply(dataSets, nrow))
  numFeaturesPerDataset = unlist(lapply(dataSets, ncol))
  totalNumFeatures = sum(numFeaturesPerDataset)
  numComponents = opts$R
  alphas = opts$alphas
  rand_start = opts$rand_start
  rhos = rep(NA, numDatasets) # form rhos, the Lipschitz constant for each data type
  
  # form full data set X, X = [X1,...Xl,...XL]
  # form full weighting matrix W, W = [W1,...Wl,...WL]
  
  # browser()
  X = data.matrix(do.call(cbind.data.frame, dataSets))
  W = 1 * (!is.na(X))
  
  W[is.na(X)] <- 0      #redundant?
  X[is.na(X)] <- 0
  
  
  #implement Wy here
  Wy <- 1 * (!is.na(y))
  Wy[is.na(y)] <- 0
  
  y[is.na(y)] <- 0
  
  # browser()
  # initialisation
  # if(exists("A0", where=opts)){
  #   mu = t(opts$mu0)
  #   A = opts$A0
  #   B = opts$B0
  # } else{
  
  if(rand_start == 1){ # use random initialization
    mu = matrix(data=0, nrow=1, ncol=totalNumFeatures)
    set.seed(rstartseed)
    A = matrix(rnorm(numSamplesPerDataset[1]*numComponents), nrow=numSamplesPerDataset[1], ncol=numComponents)
    B = matrix(rnorm(totalNumFeatures*numComponents), nrow=totalNumFeatures, ncol=numComponents)
  }
  else if(rand_start == 0 & !("P" %in% dataTypes)){ # use SCA model as initialization, Poisson distribution is not used
    mu = matrix(colMeans(X), 1, totalNumFeatures)
    
    # browser()      ###fails with CV possibly due to NAs -> also need to not split block if only 1 feature... how to handle this?
    X_svd = RSpectra::svds(scale(X,center=TRUE,scale=FALSE), numComponents, nu=numComponents, nv=numComponents)
    A = X_svd$u
    B = X_svd$v %*% diag(X_svd$d[1:numComponents])
  }
  else if(rand_start == "PLS" & !("P" %in% dataTypes)){ # use SCA model as initialization, Poisson distribution is not used
    mu = matrix(colMeans(X), 1, totalNumFeatures)
    
    # browser()      ###fails with CV possibly due to NAs -> also need to not split block if only 1 feature... how to handle this?
    # X_svd = RSpectra::svds(scale(X,center=TRUE,scale=FALSE), numComponents, nu=numComponents, nv=numComponents)
    
    X_svd <- simpls.fit(X,y,ncomp = numComponents)
    
    A = X_svd$scores
    B = X_svd$loadings   #could also use projection matrix...#X_svd$v %*% diag(X_svd$d[1:numComponents])
    # Y_lods <- X_svd$Yloadings
    
    # check for approximately orthogonal A matrix.
    
    
  }
  # }
  
  # else if(rand_start == "PLS_block" & !("P" %in% dataTypes)){ # use SCA model as initialization, Poisson distribution is not used
  #   mu = matrix(colMeans(X), 1, totalNumFeatures)
  #   
  #   # browser()      ###fails with CV possibly due to NAs -> also need to not split block if only 1 feature... how to handle this?
  #   # X_svd = RSpectra::svds(scale(X,center=TRUE,scale=FALSE), numComponents, nu=numComponents, nv=numComponents)
  #   
  #   Xs <- splitDataset(X, numFeaturesPerDataset)
  #   
  #   X_svd <- simpls.fit(X,y,ncomp = numComponents)
  #   
  #   A = X_svd$scores
  #   B = X_svd$loadings   #could also use projection matrix...#X_svd$v %*% diag(X_svd$d[1:numComponents])
  #   # Y_lods <- X_svd$Yloadings
  # }
  
  else if(rand_start == 0){ # Poisson distribution version
    X_tmp = X
    for(i in 1:numDatasets){ #log transformation applied to Poisson data
      dataType = dataTypes[i]
      if(dataType == 'P'){
        columns_Xi = grabFeatureIndices(i,numFeaturesPerDataset)
        X_tmp[,columns_Xi] = log(X[,columns_Xi] + 1)
      }
    }
    mu = matrix(colMeans(X_tmp), 1, totalNumFeatures)
    X_svd = RSpectra::svds(scale(X_tmp,center=TRUE,scale=FALSE), numComponents, nu=numComponents, nv=numComponents)     # - CHANGE - ??
    A = X_svd$u
    B = X_svd$v %*% diag(X_svd$d[1:numComponents])   
  }
  
  
  else if(rand_start == "FIXED"){
    mu = t(opts$mu0)
    A = opts$A0
    B = opts$B0
  }
  
  else if(rand_start == "rand_orth"){ # use random orthogonal rank R init!
    mu = matrix(data=0, nrow=1, ncol=totalNumFeatures)
    set.seed(rstartseed)
    
    A = svd(matrix(rnorm(numSamplesPerDataset[1]*numComponents), nrow=numSamplesPerDataset[1], ncol=numComponents))$u[, 1:numComponents, drop = FALSE]
    B = matrix(rnorm(totalNumFeatures*numComponents), nrow=totalNumFeatures, ncol=numComponents)
  }
  
  #FAILS CV not first iteration - FIXED - (?) - CHECK -
  if(is.null(dim(mu)) | ncol(mu) == 1){
    mu <- t(mu)
  }
  
  # browser()
  
  # browser()
  Theta = RpESCA::ones(numSamplesPerDataset[1]) %*% mu + tcrossprod(A,B)
  
  W_k <- ginv(Theta) %*% A  #
  
  
  
  
  # Split up result to individual datasets again
  Xs = splitDataset(X, numFeaturesPerDataset)
  Ws = splitDataset(W, numFeaturesPerDataset)
  Thetas = splitDataset(Theta, numFeaturesPerDataset)
  Bs = splitDataset(as.matrix(B), numFeaturesPerDataset, direction="rows")
  
  # Prepare rhos
  for(i in 1:numDatasets){
    dataType = dataTypes[i]
    if(dataType == 'G') rhos[i] = 1
    if(dataType == 'B') rhos[i] = 0.25
    if(dataType == 'P') rhos[i] = max(exp(Thetas[[i]]))
  }
  
  
  # Calculate loss of initial estimation
  if(is.null(response)){
    # lossResult = calculatePESCAloss(dataTypes, Xs, Ws, Thetas, alphas, lambdas, Bs, numComponents, penalty)
    lossResult = calculatePESCAloss(dataTypes = dataTypes,
                                    Xs = Xs,
                                    Ws = Ws,
                                    Wy = Wy,
                                    Thetas = Thetas,
                                    alphas = alphas,
                                    lambdas = lambdas,
                                    lambdas_y = lambdas_y,
                                    Bs = Bs,
                                    numComponents = numComponents,
                                    penalty = penalty,
                                    A_k = A,
                                    y = y,
                                    alpha = alpha,
                                    opts = opts)
    
    
  }else{
    lossResult = calculatePESCARloss(dataTypes, Xs, Ws, Thetas, alphas, lambdas, Bs, numComponents, penalty,  fun_concave = "gdp", response, A_new = A, Wt = Wt)   #CHANGE hardcoded
    # calculatePESCARloss = function(dataTypes, Xs, Ws, Thetas, alphas, lambdas, Bs, numComponents, penalty="L2", fun_concave="gdp", response, A){
    
    
  }
  
  return(list("A" = A,                        # Scores of the SCA model
              "B" = B,                        # Loadings of the SCA model
              "mu" = mu,                      # Means of the SCA model
              "Theta" = Theta,                # Reconstructed data of the SCA model
              "X" = X,                        # Input datasets following X=[X1,...,XL]
              "Xs" = Xs,                      # Input datasets in list format
              "W" = W,                        # Weight matrix of X with element=1 is the value is non-zero
              "Ws" = Ws,                      # Weight matrix of each dataset in list format
              "Wy" = Wy,
              "y" = y,
              "Thetas" = Thetas,              # Reconstructed datasets of the SCA model in list format
              "Bs" = Bs,                      # Loadings of the SCA model for each dataset, in list format
              "loss" = lossResult,            # Loss of iteration 0
              "rhos" = rhos,
              "numDatasets" = numDatasets,
              "numSamplesPerDataset" = numSamplesPerDataset,
              "numFeaturesPerDataset" = numFeaturesPerDataset,
              "totalNumFeatures" = totalNumFeatures,
              "W_k" = W_k
  )) #,
}


# calculateFMS <- function(modelFac, realFac, numComponents) {
#   
#   # Project modelFac onto realFac to make sure components are not split
#   # modelFac = modelFac %*% pracma::pinv(t(modelFac) %*% modelFac) %*% t(modelFac) %*% realFac
#   
#   # Create a similarity matrix for the pairwise comparison
#   similarity_matrix <- matrix(0, nrow = ncol(modelFac), ncol = ncol(modelFac))
#   
#   # Compute pairwise cosine similarity across all modes
#   for (k in 1:ncol(modelFac)) {
#     for (l in 1:ncol(realFac)) {
#       vect1 <- as.matrix(modelFac[, k])
#       vect2 <- as.matrix(realFac[, l])
#       
#       # Cosine similarity
#       similarity_matrix[k, l] <- abs(t(vect1) %*% vect2) / (norm(vect1, "F") * norm(vect2, "F"))
#     }
#   }
#   
#   # Use Hungarian algorithm to find the best matching
#   assignment <- solve_LSAP(similarity_matrix, maximum = TRUE)
#   
#   # Permute columns in the model to fit the input components
#   # modelFac_perm <- modelFac[,assignment]
#   
#   # Calculate FMS based on the best matching
#   FMS <- sum(similarity_matrix[cbind(seq_along(assignment), assignment)])
#   
#   # Average over the number of components
#   FMS <- FMS / numComponents
#   return(FMS)
# }


calculateFMS = function(modelFac, realFac){
  similarityMatrix = matrix(0, nrow = ncol(realFac), ncol = ncol(modelFac))
  
  for (k in 1:ncol(realFac)) {
    for (l in 1:ncol(modelFac)) {
      vect1 = as.matrix(modelFac[, l])
      vect2 = as.matrix(realFac[, k])
      
      # Cosine similarity
      similarityMatrix[k, l] = abs(t(vect2) %*% vect1) / (norm(vect1, "F") * norm(vect2, "F"))
    }
  }
  
  return(similarityMatrix)
}



calculateLSI = function(lambda_hat, lambda_true){
  
  similarityMatrix = matrix(NA, nrow=ncol(lambda_true), ncol=ncol(lambda_hat))
  #lambda_hat = sweep(lambda_hat, 1, norms, FUN="*") # correct for block scaling
  
  for(k in 1:ncol(lambda_true)){
    for(l in 1:ncol(lambda_hat)){
      real = lambda_true[,k]
      hat = lambda_hat[,l]
      similarityMatrix[k,l] = 1 / (1 + sum((real - hat)^2))
    }
  }
  
  return(similarityMatrix)
}

# ######y as block
# pESCA_Ryas <- function(dataSets, y, dataTypes, lambdas, lambdas_y ,penalty="L2", fun_concave="gdp", opts=list(), alpha, ORTH_A = TRUE, rstartseed){
#   # y <- dataSets[[length(dataSets)]]
#   # Run zero-th iteration
#   alphas = opts$alphas
#   numComponents = opts$R
#   gamma <- opts$gamma
#   # dataSets <- dataSets[1:(length(dataSets) -1 )]
#   # browser()
#   
# 
#   initialisation = initialisePESCA(dataSets = dataSets,       # - WARNING -
#                                    y = y,
#                                    dataTypes = dataTypes,
#                                    lambdas = lambdas,
#                                    lambdas_y = lambdas_y,
#                                    penalty = penalty,
#                                    opts = opts,
#                                    alpha = alpha,
#                                    rstartseed = rstartseed)
#   
# 
#   X = initialisation$X
#   W = initialisation$W
#   y <- initialisation$y
#   Wy = initialisation$Wy
#   numSamples = initialisation$numSamplesPerDataset[1]
#   numDatasets = initialisation$numDatasets
#   numFeaturesPerDataset = initialisation$numFeaturesPerDataset
#   totalNumFeatures = initialisation$totalNumFeatures
#   rhos = initialisation$rhos
#   maxit = opts$maxit
#   tol_obj = opts$tol_obj
#   
#   # get update B function based on the penalty chosen.
#   update_B_fun = get(paste0("update_B_",penalty))
#   
#   # browser()
# 
#   #rcpp version of matrix multiplication of JHk and B (?)
#   mat_vec_matC <- function(X, d, Y) {
#     .Call('_RpESCA_mat_vec_matC', PACKAGE = 'RpESCA', X, d, Y)
#   }
#   
# 
#   # Diagnostics
#   result = list()
#   result$init = initialisation
#   result$A = list()
#   result$B = list()
#   result$mu = list()
#   result$Theta = list()
#   result$Sigmas = list()
#   result$loss = list()
#   result$lossChange = initialisation$loss$loss
#   result$outcome = list()
#   result$JHk = list()
#   
#   # Prepare first iteration
#   # Note: these are the only things that change from iteration to iteration.
#   A_old = initialisation$A
#   B_old = initialisation$B
#   Bs_old = initialisation$Bs
#   Theta_old  = initialisation$Theta
#   Thetas_old = initialisation$Thetas
#   Sigmas_old = initialisation$loss$Sigmas
#   loss_old = initialisation$loss$loss
#   
#   W_old <- initialisation$W_k
#   
#   # b4B <- list()
#   
# 
#   for (k in 1:maxit){
#     
#     # majorisation step for p_ESCA model
#     #--- form Hk ---
#     #--- update mu ---
#     #--- form JHk ---
#     JHk = matrix(data=NA, numSamples, totalNumFeatures)
#     mu = matrix(data=NA, 1, totalNumFeatures)
#     cs = rep(NA, totalNumFeatures) # scaling factors
#     
#     for(i in 1:numDatasets){
#       dataType = dataTypes[i]
#       Xl = initialisation$Xs[[i]]
#       Wl = initialisation$Ws[[i]]
#       Theta_l = Thetas_old[[i]]
#       alpha_l = alphas[i]
#       B_l = Bs_old[[i]]
#       rho_l = rhos[i]
#       
#       # specify the gradient of the log-partition function
#       log_partition_g <- get(paste0("log_part_",dataType,"_g"))
#       
#       # form Hk_i  CHANGE THIS FOR SUPERVISED
#       Hk_l = Theta_l - (1/rho_l) * (Wl * (log_partition_g(Theta_l) - Xl))
#       
#       if(is.null(dim(Hk_l))){
#         
#         Hk_l <- matrix(Hk_l)
#       }
#       
#       # browser()
#       # update mu_i
#       mu_l = colMeans(Hk_l)
#       columns_Xi = grabFeatureIndices(i, numFeaturesPerDataset)
#       mu[1, columns_Xi] = mu_l
#       
#       # form JHk_i
#       JHk[, columns_Xi] = scale(Hk_l, center=TRUE, scale=FALSE)
#       
#       # form scaling factors for scaled_JHk_i, scaled_Bk_i
#       cs[columns_Xi]= rep(sqrt(rho_l/alpha_l), numFeaturesPerDataset[i])
#       
#     }
#     
# 
#     # update A
#     # browser()     # use <<- to get chosen object into global env!
#     A_tmp <- mat_vec_matC(JHk, cs^2, B_old)     #same as JHk %*% diag(cs^2) %*% B_old     # 25_02_27 changed back from W_old to B_old
#                                                 #original used B_old instead of W_old here...~
#     
#     if(ORTH_A == TRUE){
#       
#       A_svd <- svd(A_tmp, nu=numComponents, nv=numComponents)
#       
#       A_new <- tcrossprod(A_svd$u, A_svd$v)
#       
#     }
#     
#     else{
#       
#       
#       A_new <- A_tmp
#       print("ORTH_A = F")
#       }   ## need to normalise columns?
# 
#     # update B
#     # browser()
#     #CHANGE - missing Wy related update and na removal
#     B_new = update_B_fun(JHk, A_new, B_old, Sigmas_old, numFeaturesPerDataset, fun_concave, alphas, rhos, lambdas, gamma, y, lambdas_y,alpha)
#     
#     # Recalculation of Theta
#     Theta_new = RpESCA::ones(numSamples) %*% mu + tcrossprod(A_new,B_new)
#     
#     W_new <- ginv(Theta_new) %*% A_new  #  - WARNING -
# 
# 
#     # Split up Theta and B
#     Thetas_new = splitDataset(Theta_new, numFeaturesPerDataset)
#     Bs_new = splitDataset(B_new, numFeaturesPerDataset, direction="rows")
# 
#     
#     # browser()
#     # Calculate loss of this iteration
#     lossResult = calculatePESCAloss(dataTypes = dataTypes,
#                                     Xs = initialisation$Xs,
#                                     Ws = initialisation$Ws,
#                                     Wy = initialisation$Wy,
#                                     Thetas = Thetas_new,
#                                     alphas = alphas,
#                                     lambdas = lambdas,
#                                     lambdas_y = lambdas_y,
#                                     Bs = Bs_new,
#                                     numComponents = numComponents,
#                                     penalty = penalty,
#                                     A_k = A_new,
#                                     y = y,
#                                     alpha = alpha,
#                                     opts = opts)
#     
#     
#     
#     loss_new = lossResult$loss
#     
#     # print("loss" = loss_new)
#     
#     # remove the all zeros columns to simplify the computation and save memory
#     Sigmas_new = lossResult$Sigmas
#     
#     # if(thr_path == 0){
#     #   nonZeros_index <- (colMeans(Sigmas_new) > 0)
#     #   if(sum(nonZeros_index) > 3){
#     #     A_new <- A_new[,nonZeros_index]
#     #     B_new <- B_new[,nonZeros_index]
#     #     Sigmas_new <- Sigmas_new[,nonZeros_index]
#     #   }
#     # }
#     result$JHk[[k]] = JHk
#     
#     # browser()
#     # Check for convergence
#     lossChange = (loss_old-loss_new)/abs(loss_old) # relative change of loss function
#     
#      if((k>1) & (lossChange < tol_obj)) break
#     
#     # Store new version of the variables for next iteration
#     A_old = A_new
#     B_old = B_new
#     Bs_old = Bs_new
#     Theta_old = Theta_new
#     Thetas_old = Thetas_new
#     Sigmas_old = Sigmas_new
#     loss_old = loss_new
#     
#     W_old <- W_new
#     
#     
#     # Diagnostics
#     result$A[[k]] = A_new
#     result$B[[k]] = B_new
#     result$mu[[k]] = mu
#     result$Theta[[k]] = Theta_new
#     result$Sigmas[[k]] = Sigmas_new
#     result$loss[[k]] = lossResult
#     result$lossChange = c(result$lossChange, lossChange)
#     
#     result$iter <- k
#     
#     # result$b4B[[k]] <- b4B
#   }
#   
#   
#   
#   result$outcome$A = A_old
#   result$outcome$B = B_old
#   result$outcome$mu = mu
#   result$outcome$loss = loss_old
#   result$outcome$Theta = Theta_old
#   result$outcome$k = k-1
#   
#   # Calculate W
#   W_k <- ginv(result$outcome$Theta) %*% result$outcome$A  #
#   # calculate size of residuals of T = XW + E from the above line
#   # Calculate the residual matrix E
#   result$outcome$A_E <- result$outcome$A - result$outcome$Theta %*% W_k
#   
#   # Calculate the size of E (Frobenius norm)
#   # size_E <- norm(E, type = "F")
#   
#   
#   # look into predicting y with only one block with this same model.
#   result$outcome$W_k <- W_k
#   
#   # Calculate variances explained
#   varExpResult = calcVarExp(dataTypes, alphas, X, W, JHk, A_new, B_new, mu, numDatasets, numComponents, numSamples, numFeaturesPerDataset)
#   result$outcome$varExp = varExpResult[[1]]
#   result$outcome$varExpPCs = varExpResult[[2]]
#   
#   return(result)
# }


######y as block
pESCA_Ryas <- function(dataSets, y, dataTypes, lambdas, lambdas_y ,penalty="L2", fun_concave="gdp", opts=list(), alpha, ORTH_A = TRUE, rstartseed){
  # y <- dataSets[[length(dataSets)]]
  # Run zero-th iteration
  alphas = opts$alphas
  numComponents = opts$R
  gamma <- opts$gamma
  # dataSets <- dataSets[1:(length(dataSets) -1 )]
  # browser()
  
  
  initialisation = initialisePESCA(dataSets = dataSets,       # - WARNING -
                                   y = y,
                                   dataTypes = dataTypes,
                                   lambdas = lambdas,
                                   lambdas_y = lambdas_y,
                                   penalty = penalty,
                                   opts = opts,
                                   alpha = alpha,
                                   rstartseed = rstartseed)
  
  
  X = initialisation$X
  W = initialisation$W
  y <- initialisation$y
  Wy = initialisation$Wy
  numSamples = initialisation$numSamplesPerDataset[1]
  numDatasets = initialisation$numDatasets
  numFeaturesPerDataset = initialisation$numFeaturesPerDataset
  totalNumFeatures = initialisation$totalNumFeatures
  rhos = initialisation$rhos
  maxit = opts$maxit
  tol_obj = opts$tol_obj
  
  # get update B function based on the penalty chosen.
  update_B_fun = get(paste0("update_B_",penalty))
  
  # browser()
  
  #rcpp version of matrix multiplication of JHk and B (?)
  mat_vec_matC <- function(X, d, Y) {
    .Call('_RpESCA_mat_vec_matC', PACKAGE = 'RpESCA', X, d, Y)
  }
  
  
  # Diagnostics
  result = list()
  result$init = initialisation
  result$A = list()
  result$B = list()
  result$mu = list()
  result$Theta = list()
  result$Sigmas = list()
  result$loss = list()
  result$lossChange = initialisation$loss$loss
  result$outcome = list()
  result$JHk = list()
  
  # Prepare first iteration
  # Note: these are the only things that change from iteration to iteration.
  A_old = initialisation$A
  B_old = initialisation$B
  Bs_old = initialisation$Bs
  Theta_old  = initialisation$Theta
  Thetas_old = initialisation$Thetas
  Sigmas_old = initialisation$loss$Sigmas
  loss_old = initialisation$loss$loss
  
  W_old <- initialisation$W_k
  
  # b4B <- list()
  
  for (k in 1:maxit){
    
    # majorisation step for p_ESCA model
    #--- form Hk ---
    #--- update mu ---
    #--- form JHk ---
    JHk = matrix(data=NA, numSamples, totalNumFeatures)
    mu = matrix(data=NA, 1, totalNumFeatures)
    cs = rep(NA, totalNumFeatures) # scaling factors
    
    for(i in 1:numDatasets){
      dataType = dataTypes[i]
      Xl = initialisation$Xs[[i]]
      Wl = initialisation$Ws[[i]]
      Theta_l = Thetas_old[[i]]
      alpha_l = alphas[i]
      B_l = Bs_old[[i]]
      rho_l = rhos[i]
      
      # specify the gradient of the log-partition function
      # log_partition_g <- get(paste0("log_part_",dataType,"_g"))
      log_partition_g <- getFromNamespace(paste0("log_part_", dataType, "_g"), "RpESCA")
      
      # form Hk_i  CHANGE THIS FOR SUPERVISED
      Hk_l = Theta_l - (1/rho_l) * (Wl * (log_partition_g(Theta_l) - Xl))
      
      if(is.null(dim(Hk_l))){
        
        Hk_l <- matrix(Hk_l)
      }
      
      # browser()
      # update mu_i
      mu_l = colMeans(Hk_l)
      columns_Xi = grabFeatureIndices(i, numFeaturesPerDataset)
      mu[1, columns_Xi] = mu_l
      
      # form JHk_i
      JHk[, columns_Xi] = scale(Hk_l, center=TRUE, scale=FALSE)
      
      # form scaling factors for scaled_JHk_i, scaled_Bk_i
      cs[columns_Xi]= rep(sqrt(rho_l/alpha_l), numFeaturesPerDataset[i])
      
    }
    
    
    
    
    ##################################################
    # browser()
    
    coefs_ <- t(A_old) %*% y
    y_hat_ <- A_old %*% coefs_
    
    
    # coefs <- MASS::ginv(A_old) %*% as.matrix(y)
    
    # browser()
    coefs <- solve(t(A_old) %*% A_old) %*% t(A_old) %*% y
    
    coefs <- t(A_old) %*% y
    
    predict_y_new <- function(X_new, W_k, regression_coefficients) {
      X_new <- as.matrix(X_new)
      
      A_new <- X_new %*% W_k
      
      # Predict y_new using the model projected scores and regression coefficients
      y_new <- A_new %*% regression_coefficients
      
      return(y_new)
    }
    Xs <- initialisation$Xs
    W_k <- ginv(do.call(cbind, Xs)) %*% A_old
    y_hat <- predict_y_new(X_new = do.call(cbind, Xs), W_k = W_k, regression_coefficients = coefs)
    
    # cor(y_hat_ %*% t(coefs_), y_hat_ %*% t(y) %*% A_old)    #all.equal(y_hat_ %*% t(coefs_), y_hat_ %*% t(y) %*% A_old)    #TRUE
    # all.equal(y_hat_ %*% t(coefs_), y_hat_ %*% t(y) %*% A_old)
    # 
    # cor(y_hat %*% t(coefs), y_hat %*% t(y) %*% A_old)       #all.equal(y_hat %*% t(coefs), y_hat %*% t(y) %*% A_old)    #TRUE
    # all.equal(y_hat %*% t(coefs), y_hat %*% t(y) %*% A_old)
    # # yhat*Beta’ = yhat*y’*A_old.
    # ##################################################
    
    
    
    
    # update A
    # browser()     # use <<- to get chosen object into global env!
    # A_tmp <- mat_vec_matC(JHk, cs^2, B_old)     #same as JHk %*% diag(cs^2) %*% B_old     # 25_02_27 changed back from W_old to B_old
    #ignore ~ original used B_old instead of W_old here... ~
    
    
    A_tmp_x <- alpha * mat_vec_matC(JHk, cs^2, B_old)     #same as JHk %*% diag(cs^2) %*% B_old     # 25_02_27 changed back from W_old to B_old
    
    A_tmp_y <- (1-alpha) * y %*% t(coefs)#y_hat %*% t(y) %*% A_old
    
    
    A_tmp <- A_tmp_x + A_tmp_y
    
    
    
    if(ORTH_A == TRUE){
      
      A_svd <- svd(A_tmp, nu=numComponents, nv=numComponents)
      
      A_new <- tcrossprod(A_svd$u, A_svd$v)
      
    }
    
    else{
      
      
      A_new <- A_tmp
      print("ORTH_A = F")
    }   ## need to normalise columns?
    
    
    ## ~ TEMPORARY ADDITIONAL LOSS CALC ~
    Bs_old = splitDataset(B_old, numFeaturesPerDataset, direction="rows")
    
    
    Theta_old = RpESCA::ones(numSamples) %*% mu + tcrossprod(A_new,B_old)
    Thetas_old = splitDataset(Theta_old, numFeaturesPerDataset)
    
    lossResult1 = calculatePESCAloss(dataTypes = dataTypes,
                                     Xs = initialisation$Xs,
                                     Ws = initialisation$Ws,
                                     Wy = initialisation$Wy,
                                     Thetas = Thetas_old,
                                     alphas = alphas,
                                     lambdas = lambdas,
                                     lambdas_y = lambdas_y,
                                     Bs = Bs_old,
                                     numComponents = numComponents,
                                     penalty = penalty,
                                     A_k = A_new,
                                     y = y,
                                     alpha = alpha,
                                     opts = opts)
    
    lossl <- lossResult1
    
    
    
    
    
    
    
    # update B
    # browser()
    #CHANGE - missing Wy related update and na removal
    B_new = update_B_fun(JHk, A_new, B_old, Sigmas_old, numFeaturesPerDataset, fun_concave, alphas, rhos, lambdas, gamma, y, lambdas_y,alpha)
    
    # Recalculation of Theta
    Theta_new = RpESCA::ones(numSamples) %*% mu + tcrossprod(A_new,B_new)
    
    
    # browser()
    # W_new <- ginv(Theta_new) %*% A_new  #  - WARNING -
    
    
    
    # Split up Theta and B
    Thetas_new = splitDataset(Theta_new, numFeaturesPerDataset)
    Bs_new = splitDataset(B_new, numFeaturesPerDataset, direction="rows")
    
    
    # browser()
    # Calculate loss of this iteration
    lossResult = calculatePESCAloss(dataTypes = dataTypes,
                                    Xs = initialisation$Xs,
                                    Ws = initialisation$Ws,
                                    Wy = initialisation$Wy,
                                    Thetas = Thetas_new,
                                    alphas = alphas,
                                    lambdas = lambdas,
                                    lambdas_y = lambdas_y,
                                    Bs = Bs_new,
                                    numComponents = numComponents,
                                    penalty = penalty,
                                    A_k = A_new,
                                    y = y,
                                    alpha = alpha,
                                    opts = opts)
    
    
    
    loss_new = lossResult$loss
    
    # print("loss" = loss_new)
    
    # remove the all zeros columns to simplify the computation and save memory
    Sigmas_new = lossResult$Sigmas
    
    
    result$JHk[[k]] = JHk
    
    # browser()
    # Check for convergence
    lossChange = (loss_old-loss_new)/abs(loss_old) # relative change of loss function
    
    if((k>1) & (lossChange < tol_obj)) break
    
    
    # vx <- varimax(B_new)$rotmat  # R x R, orthogonal
    # A_rot <- A_new %*% vx
    # B_rot <- B_new %*% vx
    # 
    # 
    
    # H <- crossprod(B_new)
    # Off <- H - diag(diag(H))
    # B_new <- B_new - (4 * 0.001) * (B_new %*% Off)    # eta * lambda_lo <- 0.01
    # 
    
    # Store new version of the variables for next iteration
    A_old = A_new
    B_old = B_new
    Bs_old = Bs_new
    Theta_old = Theta_new
    Thetas_old = Thetas_new
    Sigmas_old = Sigmas_new
    loss_old = loss_new
    
    #W_old <- W_new
    
    save_history <- if (exists("save_history", where = opts)) isTRUE(opts$save_history) else TRUE
    
    if (save_history) {
      # Diagnostics
      result$A[[k]] = A_new
      result$B[[k]] = B_new
      result$mu[[k]] = mu
      result$Theta[[k]] = Theta_new
      result$Sigmas[[k]] = Sigmas_new
      result$loss[[k]] = lossResult
      result$lossl[[k]] = lossl
    }
    result$lossChange = c(result$lossChange, lossChange)
    
    result$iter <- k
    
    # result$b4B[[k]] <- b4B
  }
  
  
  
  result$outcome$A = A_new
  result$outcome$B = B_new
  result$outcome$mu = mu
  result$outcome$loss = loss_new
  result$outcome$Theta = Theta_new
  result$outcome$k = k-1
  
  # Calculate W
  W_k <- ginv(result$outcome$Theta) %*% result$outcome$A  #
  # calculate size of residuals of T = XW + E from the above line
  # Calculate the residual matrix E
  result$outcome$A_E <- result$outcome$A - result$outcome$Theta %*% W_k
  
  # Calculate the size of E (Frobenius norm)
  # size_E <- norm(E, type = "F")
  
  
  # look into predicting y with only one block with this same model.
  result$outcome$W_k <- W_k
  
  # Calculate variances explained
  # browser()
  varExpResult = calcVarExp(dataTypes, alphas, X, W, JHk, A_new, B_new, mu, numDatasets, numComponents, numSamples, numFeaturesPerDataset)
  
  # browser()
  #add rot here
  
  # WARNING removed for real data
  # vx <- varimax(B_new)$rotmat  # R x R, orthogonal
  # A_rot <- A_new %*% vx
  # B_rot <- B_new %*% vx
  # 
  # varExpResultROT = calcVarExp(dataTypes, alphas, X, W, JHk, A_rot, B_rot, mu, numDatasets, numComponents, numSamples, numFeaturesPerDataset)
  # 
  
  
  result$outcome$varExp = varExpResult[[1]]
  result$outcome$varExpPCs = varExpResult[[2]]
  
  # #add rot here
  # result$outcome$ROT$A <- A_rot
  # result$outcome$ROT$B <- B_rot
  # result$outcome$ROT$varExp <- varExpResultROT[[1]]
  # result$outcome$ROT$varExpPCs <- varExpResultROT[[2]]
  # 
  return(result)
}





# B_new = update_B_fun(JHk, A_new, B_old, Sigmas_old, numFeaturesPerDataset, fun_concave, alphas, rhos, lambdas, gamma, y, lambdas_y,alpha)


update_B_L2 <- function(JHk, A, B0, Sigmas0, d, fun_concave, alphas, rhos, lambdas, gamma, y, lambdas_y, alpha) {
  sumd <- sum(d)
  nDataSets <- length(d)
  n <- dim(A)[1]
  R <- dim(A)[2]
  hfun_sg <- get(paste0(fun_concave, "_sg"))  # super gradient of penalty function
  B <- matrix(NA, sumd, R)
  
  for (i in 1:nDataSets) {
    columns_Xi <- index_Xi(i, d)
    JHk_i <- JHk[, columns_Xi]
    JHkitA <- crossprod(JHk_i, A)
    alpha_i <- alphas[i]
    rho_i <- rhos[i]
    
    B0_l <- B0[columns_Xi,]
    # browser()
    
    weight_i <- sqrt(d[i])  # weight for L2 norm
    lambda_i <- lambdas[i] * weight_i * alpha_i / rho_i
    
    
    # # Solve for regression coefficients (beta)
    # coefs <- solve(t(A) %*% A) %*% t(A) %*% y
    # 
    
    if (alpha < 1) {
      coefs <- solve(t(A) %*% A) %*% t(A) %*% y
    } else {
      coefs <- rep(0, R)                # makes group_lasso_penalty zero
    }
    
    
    for (r in 1:R) {
      # Calculate the weights of the penalty according to previous sigma0_ir
      sigma0_ir <- Sigmas0[i, r]
      omega_ir <- hfun_sg(sigma0_ir, gamma = gamma, lambda = 1)  # weights
      
      # Group lasso penalty: combining coefs (beta) and A
      # group_lasso_penalty <- lambdas_y * abs(coefs[r]) * norm(A[, r], "2")
      # group_lasso_penalty <- lambdas_y * abs(coefs[r]) * norm(B0_l[,r], "2")        # - WARNING - 
      
      group_lasso_penalty <- lambdas_y * abs(coefs[r]) * norm(B0_l[,r], "2")^2    ## for update function
      
      
      #25/04/04
      # group_lasso_penalty <- lambdas_y * abs(coefs[r]) #* norm(B0[,r], "2")^2    ## non-CLD
      
      
      # group_lasso_penalty <- lambdas_y * abs(coefs[r]) * norm(B0_l[,r], "2")^2 
      # y_penalty <- (1 - alpha) * (y_contribution) + group_lasso_penalty
      
      
      # combined_penalty <- alpha * lambda_i * omega_ir + y_penalty
      
      # browser()
      # Combined penalty (with lambda_i for X and group lasso for y)
      combined_penalty <- alpha * lambda_i * omega_ir + ((1-alpha) * group_lasso_penalty)
      # combined_penalty <- lambda_i * omega_ir + group_lasso_penalty
      
      
      # Apply proximal operator with the combined penalty
      JHkitA_r <- JHkitA[, r]
      JHkitA_r_norm <- norm(JHkitA_r, "2")
      
      B_ir <- max(0, 1 - (combined_penalty / JHkitA_r_norm)) * JHkitA_r
      B[columns_Xi, r] <- B_ir
    }
  }
  
  return(B)
}





calculatePESCAloss <- function(dataTypes, Xs, Ws, Wy, Thetas, alphas, lambdas, lambdas_y, Bs, numComponents, penalty="L2", fun_concave="gdp", A_k, y, alpha, opts = opts) {
  
  #sample y here for y_v^hat rmse
  
  
  numDatasets <- length(Xs)
  f_obj_per_dataset <- rep(0, numDatasets)
  g_obj_per_dataset <- rep(0, numDatasets)
  Sigmas <- matrix(data=0, numDatasets, numComponents)
  penalty_fun <- get(paste0("penalty_concave_", penalty))
  
  y_penalty <- 0
  
  for (i in 1:numDatasets) {
    dataType <- dataTypes[i]
    Xl <- as.matrix(Xs[[i]])
    Wl <- as.matrix(Ws[[i]])
    Theta_l <- as.matrix(Thetas[[i]])
    alpha_l <- alphas[[i]]
    lambda_l <- lambdas[[i]]
    
    Bl <- if (is.null(dim(Bs[[i]]))) t(as.matrix(Bs[[i]])) else as.matrix(Bs[[i]])
    
    # Loss function for the required data type
    # log_partition <- get(paste0("log_part_", dataType))
    log_partition <- getFromNamespace(paste0("log_part_", dataType), "RpESCA")
    
    # browser()
    # Calculation of first loss term (X-related)
    f_result <- (1/alpha_l) * (trace_fast(Wl, log_partition(Theta_l)) - trace_fast(Theta_l, Xl))
    
    
    # f_result <- (1/alpha_l) * (
    #   trace_fast(Wl, log_partition(Theta_l))
    #   - trace_fast(Wl, Theta_l * Xl)
    # )
    # 
    f_obj_per_dataset[i] <- f_result
    # browser()
    # Calculation of second loss term (penalty on B)
    
    
    g_penalty <- penalty_fun(Bl, fun_concave, opts$gamma, numComponents)
    
    
    g_result <- lambda_l * g_penalty$out
    g_obj_per_dataset[i] <- g_result
    
    
    Sigmas[i,] <- g_penalty$sigma
    
  }
  
  
  # Solve for regression coefficients (beta)
  coefs <- solve(t(A_k) %*% A_k) %*% t(A_k) %*% y
  
  
  
  # Prediction of y_hat
  W_k <- ginv(do.call(cbind, Xs)) %*% A_k
  y_hat <- predict_y_new(X_new = do.call(cbind, Xs), W_k = W_k, regression_coefficients = coefs)
  
  
  if (alpha < 1) {
    # Calculate RMSE between y and y_hat
    y_contribution <- calculate_rmse(y, y_hat)     ## - CHANGE - add dispersion estimate weighting as with f_obj....
    
    
    
    # browser()
    
    # # old 24_11_02
    # group_lasso_penalty <- lambdas_y * sum(sapply(1:numComponents, function(r) {
    #   sqrt(coefs[r]^2 + sum(sapply(1:numDatasets, function(l) {
    #     norm(Bs[[l]][, r], "2")^2
    #   })))
    # }))
    # 
    
    
    ## much longer iterations slightly less overlapping loadings ....
    # group_lasso_penalty <- lambdas_y * abs(coefs[r]) * norm(B0_l[,r], "2")^2    ## for update function above
    # group_lasso_penalty <- lambdas_y * sum(sapply(1:numComponents, function(r) {
    #   abs(coefs[r]) * sum(sapply(1:numDatasets, function(l) {
    #     norm(Bs[[l]][, r], "2")^2      #CLD
    #   }))
    # }))
    
    
    #   # - 25_03_17 -
    B00 <- do.call(rbind,Bs)
    # group_lasso_penalty <- lambdas_y * sum(sapply(1:numComponents, function(r) {
    #   sqrt(coefs[r]^2 + norm(B00[, r], "2")^2)
    # }))
    
    group_lasso_penalty <- lambdas_y * sum(sapply(1:numComponents, function(r) {
      abs(coefs[r]) * norm(B00[, r], "2")^2
    }))
    
    # B00 <- do.call(rbind,Bs)
    # group_lasso_penalty <- lambdas_y * sum(abs(coefs))
    
    #   
    #   diag(c(coefs)) %*% t(B00)
    
    # browser()
    
    # Final y-related penalty: RMSE + group lasso penalty
    y_penalty <- (1 - alpha) * (y_contribution + group_lasso_penalty)#/norm(data.matrix(y), type = "F")
    
  } else {
    y_contribution <- 0
    group_lasso_penalty <- 0
    y_penalty <- 0
  }
  
  
  # Sum X-related terms (f_obj and g_obj)
  f_obj <- sum(f_obj_per_dataset)
  g_obj <- sum(g_obj_per_dataset)
  X_penalty <- alpha * (f_obj + g_obj)#/norm(do.call(rbind,Xs),type = "F")
  
  # Final loss calculation with weighted X and y objectives
  loss <- X_penalty +  y_penalty   #(1 - alpha) *
  
  
  
  return(list(
    "loss" = loss,
    "f_obj" = f_obj,
    "g_obj" = g_obj,
    "f_obj_per_dataset" = f_obj_per_dataset,
    "g_obj_per_dataset" = g_obj_per_dataset,
    "Sigmas" = Sigmas,
    "y_RMSE" = y_contribution,
    "group_lasso" = group_lasso_penalty,
    "y_obj_full" = y_penalty,
    "coefs" = coefs,
    "y_hat" = y_hat
  ))
}

predict_y_new <- function(X_new, W_k, regression_coefficients) {
  X_new <- as.matrix(X_new)
  
  A_new <- X_new %*% W_k
  
  # Predict y_new using the model projected scores and regression coefficients
  y_new <- A_new %*% regression_coefficients
  
  return(y_new)
}




# works but not implemented currently - tbd with new cross-block prediction model
predict_y_partial <- function(X_new, Full_Features , W_k, regression_coefficients) {
  
  # get available features in X_new
  available_features <- which(Full_Features %in% colnames(x_new))
  
  # Subset Theta_new_partial and W_k to include only the available blocks
  Theta_subset <- X_new[, , drop = FALSE]
  W_k_subset <- W_k[available_features, , drop = FALSE]
  
  # Compute the latent representation A_new using the available blocks
  A_new <- Theta_subset %*% W_k_subset
  
  # Predict y_new
  y_new <- A_new %*% regression_coefficients
  
  return(y_new)
}


calculate_rmse <- function(y, y_hat) {
  residuals <- y - y_hat
  
  # Calculate the mean of the squared residuals
  mse <- mean(residuals^2)
  
  # Calculate the root of the mean squared error (RMSE)
  rmse <- sqrt(mse)
  
  return(rmse)
}



dataSplit <- function(dataSets, dataTypes, y, ratio_mis = 0.1) {
  # number of data sets, size of each data set
  nDataSets <- length(dataSets)  # number of data sets
  n <- rep(0, nDataSets)  # number of samples
  d <- rep(0, nDataSets)  # numbers of variables in different data sets
  for (i in 1:nDataSets) {
    n[i] <- dim(dataSets[[i]])[1]
    d[i] <- dim(dataSets[[i]])[2]
  }
  n <- n[1]
  
  # split data sets into training set and test set
  trainSets <- as.list(1:nDataSets)  # training set
  testSets <- as.list(1:nDataSets)  # test set
  indexSets <- as.list(1:nDataSets)  # index of the test set
  for (i in 1:nDataSets) {
    # index out the i-th data set
    Xi <- dataSets[[i]]
    dataType_Xi <- dataTypes[i]
    
    # generate the index of the test set
    full_ind_vec <- 1:(n * d[i])
    
    # if it is binary data, using hierachical sampling
    if (dataType_Xi == "B") {
      ones_ind_vec <- full_ind_vec[Xi == 1]
      zeros_ind_vec <- full_ind_vec[Xi == 0]
      index_Xi_ones <- sample(ones_ind_vec, round(ratio_mis * length(ones_ind_vec)))
      index_Xi_zeros <- sample(zeros_ind_vec, round(ratio_mis * length(zeros_ind_vec)))
      
      # test the sampled samples
      if (!(all(Xi[index_Xi_ones] == 1)) | !(all(Xi[index_Xi_zeros] == 0))) 
        message("the hierachical sampling does not work")
      
      index_Xi_test <- c(index_Xi_ones, index_Xi_zeros)
    } else {
      non_NaN_mat <- 1 - is.na(Xi)
      non_NaN_ind_vec <- full_ind_vec[non_NaN_mat > 0]
      index_Xi_test <- sample(non_NaN_ind_vec, round(ratio_mis * length(non_NaN_ind_vec)))
    }
    
    # generate the train set
    Xi_train <- Xi
    Xi_train[index_Xi_test] <- NA
    trainSets[[i]] <- Xi_train
    
    # generate the test set
    Xi_test <- Xi[index_Xi_test]
    testSets[[i]] <- Xi_test
    indexSets[[i]] <- index_Xi_test
  }
  
  # apply the same level of masking to the response vector y
  n_y <- length(y)  # Number of samples in y
  y_index_full <- 1:n_y  # Create index for y
  
  # select y values to mask
  y_index_test <- sample(y_index_full, round(ratio_mis * n_y))
  
  # generate the train and test versions of y
  y_train <- y
  y_train[y_index_test] <- NA  # Mask the test set indices in y
  
  # The test set for y will contain the original values at the masked indices
  y_test <- y[y_index_test]
  
  # return
  result <- list()
  result$trainSets <- trainSets
  result$testSets <- testSets
  result$indexSets <- indexSets
  result$y_train <- y_train  # Masked y for training
  result$y_test <- y_test    # Test set for y
  result$y_index_test <- y_index_test  # Indices of masked y
  return(result)
}




penalty_concave_L2 <- function(B_i, fun_concave, gamma, R) {
  weight_i <- sqrt(dim(B_i)[1])  # weight when L2 norm is used
  hfun <- get(fun_concave)  # name of concave function
  
  out <- 0
  sigmas <- matrix(data = 0, 1, R)
  for (r in 1:R) {
    sigma_ir <- norm(B_i[, r], "2")  # sigma_{lr} = ||b_{lr}||_2
    sigmas[1, r] <- sigma_ir
    # browser()
    out <- out + hfun(sigma_ir, gamma = gamma, lambda = 1)
  }
  out <- weight_i * out
  
  result <- list()
  result$sigmas <- sigmas
  result$out <- out
  
  return(result)
}


#' Compute CV errors
#'
#' This function will compute CV errors for a specific model
#'
#' @param splitedData output of function \code{dataSplit}
#' @param dataTypes the data types for each data set
#' @param alphas dispersion parameters for each data set
#' @param ThetaHat estimated Theta
#' @param d a numeric vector contains the number of variables of data sets
#'
#' @return This function returns a vector contains CV errors
#'
#' @examples
#' \dontrun{cvError_comput(splitedData,dataTypes,alphas,ThetaHat,d)}
cvError_comput <- function(splitedData, dataTypes, alphas, ThetaHat, d) {
  
  nDataSets <- length(d)
  testSets <- splitedData$testSets
  indexSets <- splitedData$indexSets
  
  testError_vec <- rep(0, nDataSets)
  
  for (i in 1:nDataSets) {
    # index out ThetaHat_Xi
    columns_Xi <- index_Xi(i, d)
    ThetaHat_Xi <- ThetaHat[, columns_Xi]
    
    # compute the CV error
    index_Xi_test <- indexSets[[i]]
    
    
    Xi_test <- testSets[[i]]
    dataType_Xi <- dataTypes[i]
    if (dataType_Xi == "G") {
      testError_Xi <- (1/alphas[i]) * 0.5 * norm(Xi_test - ThetaHat_Xi[index_Xi_test], "2")^2
    } else if (dataType_Xi == "B") {
      # testError_Xi <- (1/alphas[i]) * obj_logistic(Xi_test, ThetaHat_Xi[index_Xi_test])
      
      # browser()
      testError_Xi <- (1 / alphas[i]) * getFromNamespace("obj_logistic", "RpESCA")(
        Xi_test,
        ThetaHat_Xi[index_Xi_test]
      )
    }
    testError_vec[i] <- testError_Xi
  }
  
  # return
  return(testError_vec)
}




pesca_FS <- function(PESCAmodel, Data, y, spikes = NULL, type = c("yasblock","PCR")){
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
  
  if(!is.null(spikes)){
    RP <- prod(which(names(cands) %in% spikes))^(1/length(spikes))
  }else(RP <- "NA")
  return(list(RP, scaled_lods, cands))
}



pfs <- function(A, B, Data, y, spikes = NULL, type = c("yasblock","PCR")){
  
  stp <- t(A) %*% y
  
  if(type == "PCR"){
    scaled_lods <- B %*% stp
    rownames(scaled_lods) <- colnames(do.call(cbind,Data))[]
    
    
  }else{
    scaled_lods <- B %*% stp
    rownames(scaled_lods) <- colnames(do.call(cbind,Data))[-nrow(B)]
    
  }
  
  
  
  cands <- scaled_lods[order(abs(scaled_lods[,1]), decreasing = T),]
  
  if(!is.null(spikes)){
    RP <- prod(which(names(cands) %in% spikes))^(1/length(spikes))
    
  }
}

RV_2 <- function(X,Y){
  
  cov_x <- X %*% t(X)
  cov_y <- Y %*% t(Y)
  
  cov_x_wig <- c(cov_x - diag(diag(cov_x)))
  cov_y_wig <- c(cov_y - diag(diag(cov_y)))
  
  num <-  t(cov_x_wig) %*% cov_y_wig
  denom <- sqrt(t(t(cov_x_wig) %*% cov_x_wig) %*% (t(cov_y_wig) %*% cov_y_wig))
  
  num/denom
  
}



calculatePESCAloss <- function(dataTypes, Xs, Ws, Wy, Thetas, alphas, lambdas, lambdas_y,
                               Bs, numComponents, penalty = "L2", fun_concave = "gdp",
                               A_k, y, alpha, opts = opts) {
  
  numDatasets <- length(Xs)
  f_obj_per_dataset <- rep(0, numDatasets)
  g_obj_per_dataset <- rep(0, numDatasets)
  Sigmas <- matrix(data = 0, numDatasets, numComponents)
  penalty_fun <- get(paste0("penalty_concave_", penalty))
  
  y_penalty <- 0
  
  for (i in 1:numDatasets) {
    dataType <- dataTypes[i]
    Xl <- as.matrix(Xs[[i]])
    Wl <- as.matrix(Ws[[i]])
    Theta_l <- as.matrix(Thetas[[i]])
    alpha_l <- alphas[[i]]
    lambda_l <- lambdas[[i]]
    
    Bl <- if (is.null(dim(Bs[[i]]))) t(as.matrix(Bs[[i]])) else as.matrix(Bs[[i]])
    
    # log_partition <- get(paste0("log_part_", dataType))
    log_partition <- getFromNamespace(paste0("log_part_", dataType), "RpESCA")
    
    # X-part 
    f_result <- (1/alpha_l) * (trace_fast(Wl, log_partition(Theta_l)) - trace_fast(Theta_l, Xl))
    f_obj_per_dataset[i] <- f_result
    
    # concave penalty majorization per dataset 
    g_penalty <- penalty_fun(Bl, fun_concave, opts$gamma, numComponents)
    g_obj_per_dataset[i] <- lambdas[[i]] * g_penalty$out
    Sigmas[i, ] <- g_penalty$sigma
  }
  
  # browser()
  # OLS beta on current scores 
  # coefs <- solve(t(A_k) %*% A_k) %*% t(A_k) %*% y
  # coefs <- MASS::ginv(A_k) %*% as.matrix(y)
  coefs <- t(A_k) %*% y
  
  # Predict y_hat 
  W_k <- ginv(do.call(cbind, Xs)) %*% A_k
  y_hat <- predict_y_new(X_new = do.call(cbind, Xs), W_k = W_k, regression_coefficients = coefs)
  
  # browser()
  
  if (alpha < 1) {
    # RMSE term 
    y_contribution <- calculate_rmse(y, y_hat)
    
    # Concatenate loadings across blocks: columns correspond to components r
    B00 <- do.call(rbind, Bs)
    
    # Quadratic coupling on all X: sum_r |beta_r| * ||b_:,r||_2^2
    group_lasso_penalty <- lambdas_y * sum(sapply(1:numComponents, function(r) {
      abs(coefs[r]) * sum(B00[, r]^2)  # ||b_:,r||_2^2
    }))
    
    y_penalty <- (1 - alpha) * (y_contribution + group_lasso_penalty)
  } else {
    y_contribution <- 0
    group_lasso_penalty <- 0
    y_penalty <- 0
  }
  
  f_obj <- sum(f_obj_per_dataset)
  g_obj <- sum(g_obj_per_dataset)
  X_penalty <- alpha * (f_obj + g_obj)
  
  loss <- X_penalty + y_penalty
  
  return(list(
    "loss" = loss,
    "f_obj" = f_obj,
    "g_obj" = g_obj,
    "f_obj_per_dataset" = f_obj_per_dataset,
    "g_obj_per_dataset" = g_obj_per_dataset,
    "Sigmas" = Sigmas,
    "y_RMSE" = y_contribution,
    "group_lasso" = group_lasso_penalty,
    "y_obj_full" = y_penalty,
    "X_obj_full" = X_penalty,
    "coefs" = coefs
  ))
}


update_B_L2 <- function(JHk, A, B0, Sigmas0, d, fun_concave, alphas, rhos, lambdas, gamma,
                        y, lambdas_y, alpha) {
  sumd <- sum(d)
  nDataSets <- length(d)
  R <- ncol(A)
  hfun_sg <- get(paste0(fun_concave, "_sg"))  # supergradient of concave penalty
  B <- matrix(NA, sumd, R)
  eps <- 1e-12
  
  # current beta (OLS on A; same rule you use in the loss)
  coefs <- if (alpha < 1) {
    solve(t(A) %*% A) %*% t(A) %*% y
  } else {
    matrix(0, nrow = R, ncol = 1)
  }
  
  for (i in 1:nDataSets) {
    columns_Xi <- index_Xi(i, d)
    JHk_i <- JHk[, columns_Xi]
    JHkitA <- crossprod(JHk_i, A)
    alpha_i <- alphas[i]
    rho_i <- rhos[i]
    
    # dataset scaling as in your code
    weight_i <- sqrt(d[i])
    lambda_i <- lambdas[i] * weight_i * alpha_i / rho_i
    
    B0_l <- B0[columns_Xi, , drop = FALSE]
    
    for (r in 1:R) {
      # GDP/PESCA weight from previous sigma
      sigma0_ir <- Sigmas0[i, r]
      omega_ir <- hfun_sg(sigma0_ir, gamma = gamma, lambda = 1)
      
      # ---- Quadratic concatenated coupling -> shared ridge for component r ----
      # c_r is identical across all blocks l for this r
      c_r <- (1 - alpha) * lambdas_y * abs(coefs[r])
      
      # PESCA data term and its norm
      JHkitA_r <- JHkitA[, r]
      JHkitA_r_norm <- max(norm(JHkitA_r, "2"), eps)
      
      # PESCA threshold (group-lasso), unchanged
      tau <- alpha * lambda_i * omega_ir
      
      if (JHkitA_r_norm <= tau) {
        B_ir <- rep(0, length(JHkitA_r))
      } else {
        # scaled soft-threshold: divide by (1 + 2*c_r)
        B_ir <- (1 / (1 + 2 * c_r)) * (1 - tau / JHkitA_r_norm) * JHkitA_r
      }
      
      B[columns_Xi, r] <- B_ir
    }
  }
  
  return(B)
}




var_exp_lr <- function(X_blocks, B) {
  L <- length(X_blocks)
  d <- vapply(X_blocks, ncol, 0L)
  R <- ncol(B)
  # block row indices for B
  idx <- split(seq_len(sum(d)), rep(seq_len(L), d))
  # denominators: ||X_l - 1*mu_l||_F^2
  
  
  denom <- sapply(X_blocks, function(X) {
    Xc <- sweep(X, 2, colMeans(X), "-")
    sum(Xc^2)
  })
  
  ve <- matrix(NA_real_, L, R)
  for (l in seq_len(L)) {
    Bl <- B[idx[[l]], , drop = FALSE]
    num <- colSums(Bl^2)
    ve[l, ] <- num / denom[l]
  }
  rownames(ve) <- paste0("block", seq_len(L))
  colnames(ve) <- paste0("PC", seq_len(R))
  
  attr(ve, "block_total") <- rowSums(ve)
  attr(ve, "full_per_component") <- colSums(B^2) / sum(denom)
  attr(ve, "full_total") <- sum(colSums(B^2)) / sum(denom)
  ve
}




match_components_and_LSI <- function(sigma_hat, sigma_true, FMS, tau = 0.6) {
  stopifnot(is.matrix(sigma_hat), is.matrix(sigma_true), is.matrix(FMS))
  # 1) Binarise sigma_hat
  sigma_hat_bin <- binarise_cols_relmax(sigma_hat, tau = tau)
  
  # 2) LSI matrix: rows = TRUE components, cols = HAT components
  LSI <- calculateLSI(sigma_hat_bin, sigma_true)
  
  # 3) Combined similarity (keep it simple: mean of the two)
  Sim <- (FMS + LSI) / 2  # rows: TRUE, cols: HAT
  
  # 4) Pad to square and solve LSAP (Hungarian)
  nT <- nrow(Sim); nH <- ncol(Sim)
  N  <- max(nT, nH)
  Saug <- matrix(0, N, N)
  Saug[1:nT, 1:nH] <- Sim
  
  assign <- clue::solve_LSAP(Saug, maximum = TRUE)  # length N; for each TRUE row -> chosen HAT col
  
  # 5) Extract mappings with NA for dummy matches
  map_true_to_hat <- assign[seq_len(nT)]
  map_true_to_hat[map_true_to_hat > nH] <- NA_integer_
  
  map_hat_to_true <- match(seq_len(nH), assign)  # inverse
  map_hat_to_true[is.na(map_hat_to_true) | map_hat_to_true > nT] <- NA_integer_
  
  # 6) Scores on matched pairs only
  paired_true <- which(!is.na(map_true_to_hat))
  if (length(paired_true) == 0L) {
    return(list(
      LSI_global        = 1 / (1 + 0),  # degenerate (no pairs); treat as perfect? or 0?
      FMS_matched_mean  = NA_real_,
      Sim_matched_mean  = NA_real_,
      map_true_to_hat   = map_true_to_hat,
      map_hat_to_true   = map_hat_to_true
    ))
  }
  idx <- cbind(paired_true, map_true_to_hat[paired_true])  # [TRUE, HAT]
  
  FMS_matched_mean <- mean(FMS[idx])
  Sim_matched_mean <- mean(Sim[idx])
  
  # 7) Global LSI via aligned binary matrices (only on matched subset)
  sigma_true_sub <- sigma_true[, paired_true, drop = FALSE]
  sigma_hat_sub  <- sigma_hat_bin[, map_true_to_hat[paired_true], drop = FALSE]
  
  SSR <- sum((sigma_true_sub - sigma_hat_sub)^2)
  LSI_global <- 1 / (1 + SSR)
  
  list(
    LSI_global       = LSI_global,
    FMS_matched_mean = FMS_matched_mean,
    Sim_matched_mean = Sim_matched_mean,
    map_true_to_hat  = map_true_to_hat,  # length nT; NA if TRUE matched to dummy
    map_hat_to_true  = map_hat_to_true   # length nH; NA if HAT not selected
  )
}

binarise_cols_relmax <- function(M, tau = 0.6) {
  stopifnot(is.matrix(M), tau >= 0, tau <= 1)
  cm <- apply(M, 2, max, na.rm = TRUE)
  thr <- tau * cm
  # if a column max is 0 or NA, make that column all zeros
  thr[!is.finite(thr) | cm <= 0] <- Inf
  sig_hat <- (sweep(M, 2, thr, `>=`)) * 1L
}


pESCA_CV_DEV <-function(dataSets, dataTypes, y,
                        lambdas_CV=NULL, lambdas_CVy=NULL, penalty='L2',
                        fun_concave='gdp', alpha, opts=list(),ORTH_A = TRUE, rstartseed = 1){
  # check if the inputs satisfy the requirements
  stopifnot(class(dataSets) == "list")
  stopifnot(class(penalty) == "character")
  stopifnot(class(fun_concave) == "character")
  if(length(dataTypes)==1){dataTypes <- unlist(strsplit(dataTypes, split=""))}
  if(exists('quiet', where=opts)){quiet <- opts$quiet} else{quiet<-0};
  if(exists('thr_path', where=opts)){thr_path <- opts$thr_path} else{thr_path<-0};
  
  # number of datasets, size of each dataset
  nTries <- length(lambdas_CV)
  nTriesy <- length(lambdas_CVy)
  nDataSets <- length(dataSets) # number of data sets
  n <- rep(0,nDataSets)  # number of samples
  d <- rep(0,nDataSets)  # numbers of variables in different data sets
  for(i in 1:nDataSets){
    n[i] <- dim(dataSets[[i]])[1]
    d[i] <- dim(dataSets[[i]])[2]
  }
  if(length(unique(as.factor(n)))!=1)
    stop("multiple datasets have unequal sample size")
  n <- n[1]
  sumd <- sum(d) # total number of variables
  
  # default dispersion parameters alphas
  if(exists('alphas', where=opts)){alphas<-opts$alphas} else{alphas<-rep(1,nDataSets)};
  
  
  
  
  # create zero matrix to hold results
  cvErrors_mat <- matrix(data=0, # +1 is used for the sum of all the Xi
                         nrow=nTries*nTriesy, ncol=nDataSets+2)
  
  # model selection process
  opts_inner <- opts
  
  
  # split data sets into training set and test set
  splitedData <- PESCAR:::dataSplit(dataSets=dataSets,       #CHANGE
                           dataTypes=dataTypes,
                           y = y,
                           ratio_mis=0.1)
  trainSets <- splitedData$trainSets
  
  
  # save the parameters during the model selection
  inits <- as.list(1:nTries)
  if(thr_path == 1){outs <- as.list(1:nTries)}
  
  l <- 1
  
  TrainModel <- list()
  
  # model selection
  for(j in 1:nTries){
    for(k in 1:nTriesy){
      
      lambda <- lambdas_CV[j]
      lambda_y <- lambdas_CVy[k]
      # using the training set to construct a PESCAR model
      lambdas <- lambda*rep(1,nDataSets)
      
      # browser()
      #generalise this for B
      if(is.list(rstartseed)) {opts_inner$A0 <- rstartseed[[1]]; 
      opts_inner$B0 <- rstartseed[[2]];
      opts_inner$mu0 <- rstartseed[[3]]}
      
      
      trainModel <- pESCA_Ryas(dataSets = trainSets,
                               dataTypes = dataTypes,
                               y = splitedData$y_train,
                               lambdas = lambdas,
                               lambdas_y = lambda_y,
                               penalty=penalty,
                               fun_concave=fun_concave,
                               alpha = alpha,
                               opts=opts_inner,
                               ORTH_A = ORTH_A,
                               rstartseed = rstartseed)
      
      if((trainModel$iter <= 2) & (quiet==0)){
        print("less than 3 iteration is used.")
      }
      
      TrainModel[[l]] <- trainModel
      # TrainModel[[l]] <- list(outcome = trainModel$outcome)
      
      # - WARNING -
      
      # warm start
      mu <- trainModel$outcome$mu; A <- trainModel$outcome$A; B <- trainModel$outcome$B
      # opts_inner$mu0 <- mu
      # opts_inner$A0 <- A
      # #
      # opts_inner$B0 <- B
      # opts_inner$R <- dim(B)[2]
      
      inits[[l]] <- opts_inner
      if(thr_path == 1){outs[[l]] <- trainModel$Sigmas}
      
      yti <- splitedData$y_index_test
      
      # - BUG HERE? -
      
      # ytri <- splitedData$
      # compute the test error - for coefs;
      coefs <- solve(t(A[,]) %*% A[,]) %*% t(A[,]) %*% y 
      
      
      ThetaHat <- RpESCA::ones(n) %*% mu + A %*% t(B)
      # ThetaHat <- ones(n) %*% mu + A %*% matrix(diag(c(coefs)), ncol = length(coefs)) %*% t(B)
      
      testError_vec <- alpha * cvError_comput(splitedData,dataTypes,alphas,ThetaHat,d)
      
      #add y rmse here and alpha weighting?
      # Prediction of y_hat
      W_k <- ginv(ThetaHat) %*% A
      y_hat <- predict_y_new(X_new = ThetaHat, W_k = W_k, regression_coefficients = coefs)
      # browser()
      
      # Calculate RMSE between y and y_hat
      y_contribution <- (1 - alpha) * calculate_rmse(y[yti], y_hat[yti])
      
      
      # browser()
      cvErrors_tmp <- c(sum(testError_vec), testError_vec, y_contribution)
      
      
      cvErrors_mat[l,] <- cvErrors_tmp
      
      l <- l + 1
    }
    
  }
  
  colnames(cvErrors_mat) <- c(paste0(rep("X_"), c("full", as.character(1:nDataSets))),"y")
  
  result_CV <- list()
  result_CV$cvErrors_mat <- cvErrors_mat
  result_CV$inits <- inits
  if(thr_path == 1){result_CV$outs <- outs}
  result_CV$TrainModel <- TrainModel
  
  return(result_CV)
}
