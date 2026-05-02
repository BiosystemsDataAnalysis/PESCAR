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
# pesca_estimate_alphas = function(datasets, numComponentsRange, CVfolds, tol=1e-6, quiet=1){
#   alphas_mean = rep(NA, length(datasets))
#   alphas_std = rep(NA, length(datasets))
#   R_selected = list()
#   cvErrors = list()
#   
#   opts = list()
#   opts$tol_obj = tol
#   opts$quiet = quiet
#   
#   for (i in 1:length(datasets)){
#     alpha_est = alpha_estimation(datasets[[i]], K = CVfolds, Rs = numComponentsRange, opts=opts) #5:20 NEED TO THINK ABOUT THIS CHANGE K
#     alphas_mean[i] = alpha_est$alphas_mean
#     alphas_std[i] = alpha_est$alphas_std
#     R_selected[[i]] = alpha_est$R_CV
#     cvErrors[[i]] = alpha_est$cvErrors
#   }
#   names(alphas_mean) = paste0('alpha_mean_', 1:length(alphas_mean))
#   names(alphas_std) = paste0('alpha_std_', 1:length(alphas_std))
#   
#   return(list(alphas_mean, alphas_std, R_selected, cvErrors))
# }
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

# tweaked for yasblock
varExp_Gaussian <- function(X, mu, A, B, Q) {
  # parameter used
  
  if(is.null(dim(X))){
    X <- matrix(X)
  }
  m = dim(X)[1]
  
  # compute the loglikelihood of mle and null model
  # browser()
  X_centered <- X - RpESCA::ones(m) %*% mu
  QX <- Q * X_centered
  
  # likelihood of the null model
  l_null <- norm(as.matrix(QX), "F")^2  # null model
  
  
  # likelihood of the full model
  E_hat <- X_centered - tcrossprod(A, B)
  QE_hat <- Q * E_hat
  l_model <- norm(as.matrix(QE_hat), "F")^2  # full model
  
  # compute the least squares of an individual PC
  R <- dim(B)[2]
  l_PCs <- rep(0, R)
  
  for (r in 1:R) {
    Ar <- A[, r]
    Br <- B[, r]
    QPCr <- Q * (tcrossprod(Ar, Br))
    l_PCs[r] <- l_null - 2 * (crossprod((as.matrix(QX) %*% Br), Ar)) + crossprod(Ar, (QPCr %*% Br))
  }
  
  # compute variation explained by each PC
  varExp_PCs <- (1 - l_PCs/l_null) * 100
  
  # total variation explained
  varExp_total <- (1 - l_model/l_null) * 100
  
  # return the results
  out <- list()
  out$varExp_total <- varExp_total
  out$varExp_PCs <- varExp_PCs
  return(out)
}

# Not yet modified
index_Xi <- function(i, ds) {
  if (i == 1) {
    columns_Xi <- 1:ds[1]
  } else {
    columns_Xi <- (sum(ds[1:(i - 1)]) + 1):sum(ds[1:i])
  }
  columns_Xi
}

# Not yet modified
penalty_concave_L2 <- function(B_i, fun_concave, gamma, R) {
  weight_i <- sqrt(dim(B_i)[1])  # weight when L2 norm is used
  hfun <- get(fun_concave)  # name of concave function
  
  out <- 0
  sigmas <- matrix(data = 0, 1, R)
  for (r in 1:R) {
    sigma_ir <- norm(B_i[, r], "2")  # sigma_{lr} = ||b_{lr}||_2
    sigmas[1, r] <- sigma_ir
    
    out <- out + hfun(sigma_ir, gamma = gamma, lambda = 1)
  }
  out <- weight_i * out
  
  result <- list()
  result$sigmas <- sigmas
  result$out <- out
  
  return(result)
}



# Not yet modified
log_part_G <- function(Theta) {
  
  0.5 * (Theta^2)
}

# Not yet modified
trace_fast <- function(X, Y) {
  # fast trace function if n>p, trace(X,Y) = trace(X'Y); if n<p, trace(X,Y) = trace(YX');
  
  stopifnot(all(dim(X) == dim(Y)))
  fast_traceC <- function(X, Y) {
    .Call('_RpESCA_fast_traceC', PACKAGE = 'RpESCA', X, Y)
  }
  
  
  result <- fast_traceC(X, Y)
  return(result)
}

# Not yet modified
gdp <- function(x, gamma, lambda) {
  # gdp penalty for non-negative values x>=0
  
  # the domain of gdp() should be in [0,+inf)
  stopifnot(x >= 0)
  
  y <- lambda * log(1 + x/gamma)
  y
}

# Not yet modified
gdp_sg <- function(x, gamma, lambda) {
  # supergradient of gdp penalty x>=0
  
  # the domain of gdp_sg() should be in [0,+inf)
  stopifnot(x >= 0)
  
  y <- lambda/(gamma + x)
  y
}

# Not yet modified
log_part_G_g <- function(Theta) {
  Theta
}



# GRvdP
splitDataset = function(df, numFeaturesPerDataset, direction="columns"){
  
  numDatasets = length(numFeaturesPerDataset)
  result = list()
  
  for(i in 1:numDatasets){
    selection = grabFeatureIndices(i, numFeaturesPerDataset)
    
    if(direction == "columns") result[[i]] = df[,selection]
    if(direction == "rows") result[[i]] = df[selection,]
  }
  
  return(result)
}

# GRvdP
# Modified from index_Xi
grabFeatureIndices = function(i, numFeaturesPerDataset) {
  #' Find feature indices of dataset i in appended dataset X = X1, X2, ..., Xn
  #' 
  #' @param i Number of the dataset
  #' @param numFeaturesPerDataset Vector containing ncol per dataset
  #' 
  #' @returns result Feature indices of dataset i in appended dataset X.
  
  if (i == 1) {
    result = 1:numFeaturesPerDataset[1]
  } 
  else{
    result <- (sum(numFeaturesPerDataset[1:(i - 1)]) + 1):sum(numFeaturesPerDataset[1:i])
  }
  
  return(result)
}

# GRvdP
calcVarExp = function(dataTypes, alphas, X, W, JHk, A, B, mu, numDatasets, numComponents, numSamples, numFeaturesPerDataset){
  
  totalNumFeatures = sum(numFeaturesPerDataset)
  
  # variation explained ratios
  varExpTotals = rep(NA, numDatasets + 1)  # +1 is for the full data set
  varExpPCs    = matrix(NA, numDatasets + 1, numComponents) # +1 is for the full data set
  
  X_full = matrix(NA, nrow=numSamples, ncol=totalNumFeatures) # combine the quantitative data and the pseudo data of nonGaussian data
  X_full = data.frame(X_full)
  for(i in 1:numDatasets){
    columns_Xi = grabFeatureIndices(i, numFeaturesPerDataset)
    Xi    <- X[,columns_Xi]
    mu_Xi <- mu[1,columns_Xi]
    if(!(dataTypes[i]=='G')){
      Xi <- JHk[,columns_Xi] + RpESCA::ones(numSamples)%*%mu_Xi
    }
    X_full[,columns_Xi] <- Xi
    
    B_Xi  <- B[columns_Xi,]
    W_Xi  <- W[,columns_Xi]
    # browser()
    
    if(is.null(dim(B_Xi))){
      
      B_Xi <- t(matrix(B_Xi))
      
    }
    
    varExp_tmp <- varExp_Gaussian(Xi,mu_Xi,A,B_Xi,W_Xi)
    varExpTotals[i] <- varExp_tmp$varExp_total
    varExpPCs[i,] <- varExp_tmp$varExp_PCs
  }
  
  # variance explained ratio of the full data set
  if(length(unique(alphas))==1){
    varExp_tmp <- varExp_Gaussian(X_full,as.vector(mu),A,B,W)
  } else{
    weighted_X  <- data.frame(matrix(NA, numSamples, totalNumFeatures))
    weighted_mu <- matrix(NA, 1, totalNumFeatures)
    weighted_B  <- matrix(NA, totalNumFeatures, numComponents)
    for (i in 1:numDatasets){
      columns_Xi = grabFeatureIndices(i, numFeaturesPerDataset)
      Xi <- X_full[,columns_Xi]
      mu_Xi <- mu[1,columns_Xi]
      B_Xi  <- B[columns_Xi,]
      alpha_i <- alphas[i]
      
      weighted_X[,columns_Xi] <- (1/sqrt(alpha_i))*Xi
      weighted_mu[1,columns_Xi] <- (1/sqrt(alpha_i))*mu_Xi
      weighted_B[columns_Xi,]  <- (1/sqrt(alpha_i))*B_Xi
    }
    
    varExp_tmp <- varExp_Gaussian(weighted_X,as.vector(weighted_mu),A,weighted_B,W)
  }
  varExpTotals[numDatasets+1] <- varExp_tmp$varExp_total
  varExpPCs[numDatasets+1,] <- varExp_tmp$varExp_PCs
  dataSets_names <- paste0(rep("X_"), c(as.character(1:numDatasets), "full"))
  PCs_names <- paste0("PC", c(as.character(1:numComponents)))
  names(varExpTotals) <- dataSets_names
  rownames(varExpPCs) <- dataSets_names
  colnames(varExpPCs) <- PCs_names
  
  # extract the structure index
  S <- matrix(data=0, numDatasets, numComponents)
  S[varExpPCs[1:numDatasets,] > 0] <- 1
  
  # save the results
  return(list(varExpTotals, varExpPCs))
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


#' pESCA model selection based on cross validation error
#'
#' This function implements a missing value based CV model selection
#' approach for the pESCA model on mutliple data sets with the same data type.
#' The details can be found in  \url{https://arxiv.org/abs/1902.06241}.
#'
#' @param lambdas_CV a vector cotains a sequence of values of lambda
#'
#' @return This function returns a list contains the results of a pESCA mdoel. \itemize{
#' \item cvErrors_mat: a matrix contains the CV errors for the full data set and each
#' single data set;
#' \item inits: a list contains the initilizations of all the constructed models;
#' \item outs: a list contains the outputs of all the constructed models;
#' }
#'
#' @examples
#' \dontrun{
#' result_CV <- pESCA_CV(dataSets, dataTypes,
#'                             lambdas_CV, penalty='L2', fun_concave='gdp', opts=opts)
#' }
#'
#' @export
pESCA_CV <-function(dataSets, dataTypes,
                    lambdas_CV=NULL, penalty='L2', fun_concave='gdp', opts=list()){
  # check if the inputs satisfy the requirements
  stopifnot(class(dataSets) == "list")
  stopifnot(class(penalty) == "character")
  stopifnot(class(fun_concave) == "character")
  if(length(dataTypes)==1){dataTypes <- unlist(strsplit(dataTypes, split=""))}
  if(exists('quiet', where=opts)){quiet <- opts$quiet} else{quiet<-0};
  if(exists('thr_path', where=opts)){thr_path <- opts$thr_path} else{thr_path<-0};
  
  # number of data sets, size of each data set
  nTries <- length(lambdas_CV)
  nDataSets <- length(dataSets) # number of data sets
  n <- rep(0,nDataSets)  # number of samples
  d <- rep(0,nDataSets)  # numbers of variables in different data sets
  for(i in 1:nDataSets){
    n[i] <- dim(dataSets[[i]])[1]
    d[i] <- dim(dataSets[[i]])[2]
  }
  if(length(unique(as.factor(n)))!=1)
    stop("multiple data sets have unequal sample size")
  n <- n[1]
  sumd <- sum(d) # total number of variables
  
  # default dispersion parameters alphas
  if(exists('alphas', where=opts)){alphas<-opts$alphas} else{alphas<-rep(1,nDataSets)};
  
  # create zero matrix to hold results
  cvErrors_mat <- matrix(data=0, # +1 is used for the sum of all the Xi
                         nrow=nTries, ncol=nDataSets+1)
  
  # model selection process
  opts_inner <- opts
  
  # split data sets into training set and test set
  splitedData <- dataSplit(dataSets=dataSets,       #CHANGE
                           dataTypes=dataTypes,
                           ratio_mis=0.1)
  trainSets <- splitedData$trainSets
  
  # save the parameters during the model selection
  inits <- as.list(1:nTries)
  if(thr_path == 1){outs <- as.list(1:nTries)}
  
  # model selection
  for(j in 1:nTries){
    lambda <- lambdas_CV[j]
    
    # using the training set to construct a ESCA model
    lambdas <- lambda*rep(1,nDataSets)
    
    trainModel <- pESCA(dataSets = trainSets,                #CHANGE to pESCA_roel
                        dataTypes = dataTypes,
                        lambdas = lambdas,
                        penalty = penalty,
                        fun_concave= fun_concave,
                        opts=opts_inner)
    if((trainModel$iter <= 2) & (quiet==0)){
      print("less than 3 iteration is used.")
    }
    
    # warm start
    mu <- trainModel$mu; A <- trainModel$A; B <- trainModel$B
    opts_inner$mu0 <- mu
    opts_inner$A0 <- A
    opts_inner$B0 <- B
    opts_inner$R <- dim(B)[2]         ##Addition
    
    inits[[j]] <- opts_inner
    if(thr_path == 1){outs[[j]] <- trainModel$Sigmas}
    
    # compute the test error
    ThetaHat <- RpESCA::ones(n) %*% t(mu) + A %*% t(B)
    testError_vec <- cvError_comput(splitedData,dataTypes,alphas,ThetaHat,d)
    
    cvErrors_tmp <- c(sum(testError_vec), testError_vec)
    cvErrors_mat[j,] <- cvErrors_tmp
  }
  
  colnames(cvErrors_mat) <- paste0(rep("X_"), c("full", as.character(1:nDataSets)))
  
  result_CV <- list()
  result_CV$cvErrors_mat <- cvErrors_mat
  result_CV$inits <- inits
  result_CV$trainSets <- trainSets
  if(thr_path == 1){result_CV$outs <- outs}
  
  return(result_CV)
}




