runPESCA = function(datasets, datatypes, alphaEstimation_numComponentsRange, alphaEstimation_CVfolds, alpha_tol=1e-4, nTries=100, penalty="L2", fun_concave="gdp", tol=1e-6, maxit=1000, gamma=1, numComponents=1, rand_start=0, thr_path=0, quiet=1){
  
  # Estimate alphas
  alpha_est_output = pesca_estimate_alphas(datasets, alphaEstimation_numComponentsRange, alphaEstimation_CVfolds, tol=alpha_tol, quiet=quiet)
  alphas = alpha_est_output[[1]]
  
  # parameters of a pESCA with concave L2norm penalty model
  opts = list()
  opts$gamma      = gamma             # hyper-parameter for the used penalty
  opts$rand_start = rand_start        # initialization method: 1=random,0=SCA
  opts$tol_obj    = tol               # stopping criterion
  opts$maxit      = maxit             # maximum number of iterations
  opts$alphas     = alphas            # alpha estimation from pesca_estimate_alphas
  opts$R          = numComponents     # components used
  opts$thr_path   = thr_path          # generate thresholding path or not
  opts$quiet      = quiet             # don't show progress while running algorithm
  
  # Estimate lambdas
  lambda_est_output = pesca_estimate_lambdas(datasets, datatypes, alphas, opts)
  lambdas = lambda_est_output[[1]]
  
  # Create PESCA model
  pesca_model = pESCA(dataSets=datasets, dataTypes=datatypes, lambdas, penalty=penalty, fun_concave=fun_concave, opts=opts)
  
  # Return all outputs
  return(list("model"=pesca_model, "alpha"=alpha_est_output, "lambda"=lambda_est_output))
}

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

pesca_estimate_lambdas = function(datasets, datatypes, alphas, opts, nTries=100, penalty="L2", fun_concave="gdp"){
  
  # Initialize a sequence of lambdas to try
  lambdas_CV = log10_seq(from=1, to=500, length.out=nTries)
  
  # Run cross-validation
  result_CV = pESCA_CV(datasets, datatypes, lambdas_CV, penalty, fun_concave, opts=opts)
  cvErrors_mat = result_CV$cvErrors_mat
  inits = result_CV$inits
  outs = result_CV$outs
  
  # Find best lambda estimation
  index = which.min(result_CV$cvErrors_mat[,1])
  lambdas = rep(lambdas_CV[index], length(datasets))
  
  return(list(lambdas, lambdas_CV, cvErrors_mat, inits, outs))
}

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


# update_B_L2 <- function(JHk, A, B0, Sigmas0, d, fun_concave, alphas, rhos, lambdas, gamma) {
#   sumd <- sum(d)
#   nDataSets <- length(d)
#   n <- dim(A)[1]
#   R <- dim(A)[2]
#   hfun_sg <- get(paste0(fun_concave, "_sg"))  # super gradient of penalty function
#   B <- matrix(NA, sumd, R)
#   for (i in 1:nDataSets) {
#     columns_Xi <- index_Xi(i, d)
#     JHk_i <- JHk[, columns_Xi]
#     JHkitA <- crossprod(JHk_i, A)    #replace A with weights?  crossprod(JHk_i, A) == t(JHk_i) %*% A
#     #   could replace with JHk_i %*% B0[1:100,]
#     alpha_i <- alphas[i]
#     rho_i <- rhos[i]
#     weight_i <- sqrt(d[i])  # weight for L2 norm
#     lambda_i <- lambdas[i] * weight_i * alpha_i/rho_i
#     
#     for (r in 1:R) {
#       # form weights of the penalty according to previous sigma0_ir
#       sigma0_ir <- Sigmas0[i, r]
#       omega_ir <- hfun_sg(sigma0_ir, gamma = gamma, lambda = 1)  # weights
#       
#       # proximal operator of L2 norm
#       JHkitA_r <- JHkitA[, r]
#       lambda_ir <- lambda_i * omega_ir
#       JHkitA_r_norm <- norm(JHkitA_r, "2")
#       
#       B_ir <- max(0, 1 - (lambda_ir/JHkitA_r_norm)) * JHkitA_r
#       B[columns_Xi, r] <- B_ir
#     }
#   }
#   
#   return(B)
# }
# 
# 
# # FW WIP  #can we use this as is but for weights?
# #change lambda estimation to work with y?..
# 
# #CHANGE -> need to account for y in z with indices
# # and possibly in penalty
# update_B_L2_dev <- function(JHk, W0, B0, Sigmas0, d, fun_concave, alphas, rhos, lambdas, gamma) {
#   sumd <- sum(d)
#   nDataSets <- length(d)
#   n <- dim(X)[1]
#   R <- dim(W0)[2]
#   hfun_sg <- get(paste0(fun_concave, "_sg"))  # super gradient of penalty function
#   B <- matrix(NA, sumd, R)
#   for (i in 1:nDataSets) {
#     columns_Xi <- index_Xi(i, d)
#     JHk_i <- JHk[, columns_Xi]
#     
#     A <- JHk_i %*% W0[columns_Xi, ]
#     JHkitA <- crossprod(JHk_i, A)    
#     
#     alpha_i <- alphas[i]
#     rho_i <- rhos[i]
#     weight_i <- sqrt(d[i])  # weight for L2 norm
#     lambda_i <- lambdas[i] * weight_i * alpha_i/rho_i
#     
#     for (r in 1:R) {
#       # form weights of the penalty according to previous sigma0_ir
#       sigma0_ir <- Sigmas0[i, r]
#       omega_ir <- hfun_sg(sigma0_ir, gamma = gamma, lambda = 1)  # weights
#       
#       # proximal operator of L2 norm
#       JHkitA_r <- JHkitA[, r]
#       lambda_ir <- lambda_i * omega_ir
#       JHkitA_r_norm <- norm(JHkitA_r, "2")
#       
#       B_ir <- max(0, 1 - (lambda_ir/JHkitA_r_norm)) * JHkitA_r
#       B[columns_Xi, r] <- B_ir
#     }
#   }
#   
#   return(B)
# }

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

# pESCA_Roel = function(dataSets, dataTypes, lambdas, penalty="L2", fun_concave="gdp", opts=list()){
#   
#   # Run zero-th iteration
#   alphas = opts$alphas
#   numComponents = opts$R
#   initialisation = initialisePESCA(dataSets, dataTypes, lambdas, penalty, opts)
#   
#   X = initialisation$X
#   W = initialisation$W
#   numSamples = initialisation$numSamplesPerDataset[1]
#   numDatasets = initialisation$numDatasets
#   numFeaturesPerDataset = initialisation$numFeaturesPerDataset
#   totalNumFeatures = initialisation$totalNumFeatures
#   rhos = initialisation$rhos
#   maxit = opts$maxit
#   tol_obj = opts$tol_obj
#   
#   # Grab the relevant update B function based on the penalty chosen.
#   update_B_fun = get(paste0("update_B_",penalty))
#   
#   #rcpp version of matrix multiplication of JHk and B (?)
#   mat_vec_matC <- function(X, d, Y) {
#     .Call('_RpESCA_mat_vec_matC', PACKAGE = 'RpESCA', X, d, Y)
#   }
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
#   # b4B <- list()
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
#       # update mu_i
#       mu_l = colMeans(Hk_l)
#       columns_Xi = grabFeatureIndices(i, numFeaturesPerDataset)
#       # browser()
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
#     # update A
#     A_tmp = mat_vec_matC(JHk, cs^2, B_old)
#     A_svd = svd(A_tmp, nu=numComponents, nv=numComponents)
#     A_new = tcrossprod(A_svd$u, A_svd$v)
#     
#     # update B
#     
#     #to check if inputs here are the same as original algorithm
#     #replace totalNumFeatures with a vector of features per block - numFeaturesPerDataset 
#     
#     b4B <- list(JHk, A_new, B_old, Sigmas_old, numFeaturesPerDataset, fun_concave, alphas, rhos, lambdas, gamma)
#     
#     B_new = update_B_fun(JHk, A_new, B_old, Sigmas_old, numFeaturesPerDataset, fun_concave, alphas, rhos, lambdas, gamma)
#     
#     # Recalculation of Theta
#     Theta_new = RpESCA::ones(numSamples) %*% mu + tcrossprod(A_new,B_new)
#     
#     # Split up Theta and B
#     Thetas_new = splitDataset(Theta_new, numFeaturesPerDataset)
#     Bs_new = splitDataset(B_new, numFeaturesPerDataset, direction="rows")
#     
#     # Calculate loss of this iteration
#     lossResult = calculatePESCAloss(dataTypes, initialisation$Xs, initialisation$Ws, Thetas_new, alphas, lambdas, Bs_new, numComponents, penalty)
#     loss_new = lossResult$loss
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
#     
#     # Check for convergence
#     lossChange = (loss_old-loss_new)/abs(loss_old) # relative change of loss function
#     if((k>1) & ((lossChange < tol_obj))) break
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
#     # Diagnostics
#     result$A[[k]] = A_new
#     result$B[[k]] = B_new
#     result$mu[[k]] = mu
#     result$Theta[[k]] = Theta_new
#     result$Sigmas[[k]] = Sigmas_new
#     result$loss[[k]] = lossResult
#     result$lossChange = c(result$lossChange, lossChange)
#     
#     result$b4B[[k]] <- b4B
#   }
#   
#   result$outcome$A = A_old
#   result$outcome$B = B_old
#   result$outcome$mu = mu
#   result$outcome$loss = loss_old
#   result$outcome$Theta = Theta_old
#   result$outcome$k = k-1
#   
#   # Calculate variances explained
#   varExpResult = calcVarExp(dataTypes, alphas, X, W, JHk, A_old, B_old, mu, numDatasets, numComponents, numSamples, numFeaturesPerDataset)
#   result$outcome$varExp = varExpResult[[1]]
#   result$outcome$varExpPCs = varExpResult[[2]]
#   
#   return(result)
# }
#' 
#' # GRvdP, FW;slightly tweaked.....    CV issue here 
#' initialisePESCA = function(dataSets, dataTypes, lambdas, penalty="L2", opts=list(), response = NULL){
#'   #' Initialises many aspects to run PESCA
#'   #' 
#'   #' @param dataSets A list of datasets.
#'   #' @param dataTypes Data type per dataset, string of letters: "G" = gaussian, "P" = poisson, "B" = binary.
#'   #' @param lambdas A list of fitted lambdas, one for each dataset.
#'   #' @param penalty Penalty type for the PESCA function. Default is L2.
#'   #'
#'   #' @returns A list of the following: A, B, mu, Theta, X, W, loss, rho, and some useful variables.
#'   
#'   numDatasets = length(dataSets)
#'   numSamplesPerDataset = unlist(lapply(dataSets, nrow))
#'   numFeaturesPerDataset = unlist(lapply(dataSets, ncol))
#'   totalNumFeatures = sum(numFeaturesPerDataset)
#'   numComponents = opts$R
#'   alphas = opts$alphas
#'   rand_start = opts$rand_start
#'   rhos = rep(NA, numDatasets) # form rhos, the Lipschitz constant for each data type
#'   
#'   # form full data set X, X = [X1,...Xl,...XL]
#'   # form full weighting matrix W, W = [W1,...Wl,...WL]
#'   
#'   # browser()
#'   X = data.matrix(do.call(cbind.data.frame, dataSets))
#'   W = 1 * (!is.na(X))
#'   
#'   W[is.na(X)] <- 0
#'   X[is.na(X)] <- 0
#'   
#'   
#'   
#'   # initialization
#'   if(exists("A0", where=opts)){
#'     mu = t(opts$mu0)
#'     A = opts$A0
#'     B = opts$B0
#'   }
#'   else if(rand_start == 1){ # use random initialization
#'     mu = matrix(data=0, nrow=1, ncol=totalNumFeatures)
#'     set.seed(1)
#'     A = matrix(rnorm(numSamplesPerDataset[1]*numComponents), nrow=numSamplesPerDataset[1], ncol=numComponents)
#'     B = matrix(rnorm(totalNumFeatures*numComponents), nrow=totalNumFeatures, ncol=numComponents)
#'   } 
#'   else if(rand_start == 0 & !("P" %in% dataTypes)){ # use SCA model as initialization, Poisson distribution is not used
#'     mu = matrix(colMeans(X), 1, totalNumFeatures)
#'     
#'     # browser()      ###fails with CV possibly due to NAs -> also need to not split block if only 1 feature... how to handle this?
#'     X_svd = RSpectra::svds(scale(X,center=TRUE,scale=FALSE), numComponents, nu=numComponents, nv=numComponents)
#'     A = X_svd$u
#'     B = X_svd$v %*% diag(X_svd$d[1:numComponents])
#'   } 
#'   else if(rand_start == 0){ # Poisson distribution version
#'     X_tmp = X
#'     for(i in 1:numDatasets){ #log transformation applied to Poisson data
#'       dataType = dataTypes[i]
#'       if(dataType == 'P'){
#'         columns_Xi = grabFeatureIndices(i,numFeaturesPerDataset)
#'         X_tmp[,columns_Xi] = log(X[,columns_Xi] + 1)
#'       }
#'     }
#'     mu = matrix(colMeans(X_tmp), 1, totalNumFeatures)
#'     X_svd = RSpectra::svds(scale(X_tmp,center=TRUE,scale=FALSE), numComponents, nu=numComponents, nv=numComponents)
#'     A = X_svd$u
#'     B = X_svd$v %*% diag(X_svd$d[1:numComponents])   ##why put size on the loadings not scores?
#'   }
#'   
#'   #FAILS CV not first iteration - FIXED - (?) - CHECK -
#'   if(is.null(dim(mu)) | ncol(mu) == 1){
#'     mu <- t(mu)
#'   }
#'   
#'   # browser()
#'   Theta = RpESCA::ones(numSamplesPerDataset[1]) %*% mu + tcrossprod(A,B)
#'   
#'   # Split up result to individual datasets again
#'   Xs = splitDataset(X, numFeaturesPerDataset)
#'   Ws = splitDataset(W, numFeaturesPerDataset)
#'   Thetas = splitDataset(Theta, numFeaturesPerDataset)
#'   Bs = splitDataset(as.matrix(B), numFeaturesPerDataset, direction="rows")
#'   
#'   # Prepare rhos
#'   for(i in 1:numDatasets){
#'     dataType = dataTypes[i]
#'     if(dataType == 'G') rhos[i] = 1
#'     if(dataType == 'B') rhos[i] = 0.25
#'     if(dataType == 'P') rhos[i] = max(exp(Thetas[[i]]))
#'   }
#'   
#'   # Calculate loss of initial estimation
#'   if(is.null(response)){
#'     # lossResult = calculatePESCAloss(dataTypes, Xs, Ws, Thetas, alphas, lambdas, Bs, numComponents, penalty)
#'     lossResult = calculatePESCAloss(dataTypes = dataTypes, 
#'                                     Xs = Xs, 
#'                                     Ws = Ws, 
#'                                     Thetas = Thetas, 
#'                                     alphas = alphas, 
#'                                     lambdas = lambdas, 
#'                                     Bs = Bs, 
#'                                     numComponents = numComponents, 
#'                                     penalty = penalty)
#'     
#'     
#'   }else{
#'     lossResult = calculatePESCARloss(dataTypes, Xs, Ws, Thetas, alphas, lambdas, Bs, numComponents, penalty,  fun_concave = "gdp", response, A_new = A, Wt = Wt)   #CHANGE hardcoded
#'     # calculatePESCARloss = function(dataTypes, Xs, Ws, Thetas, alphas, lambdas, Bs, numComponents, penalty="L2", fun_concave="gdp", response, A){
#'     
#'     
#'   }
#'   
#'   return(list("A" = A,                        # Scores of the SCA model
#'               "B" = B,                        # Loadings of the SCA model
#'               "mu" = mu,                      # Means of the SCA model
#'               "Theta" = Theta,                # Reconstructed data of the SCA model
#'               "X" = X,                        # Input datasets following X=[X1,...,XL]
#'               "Xs" = Xs,                      # Input datasets in list format
#'               "W" = W,                        # Weight matrix of X with element=1 is the value is non-zero
#'               "Ws" = Ws,                      # Weight matrix of each dataset in list format
#'               "Thetas" = Thetas,              # Reconstructed datasets of the SCA model in list format
#'               "Bs" = Bs,                      # Loadings of the SCA model for each dataset, in list format
#'               "loss" = lossResult,            # Loss of iteration 0
#'               "rhos" = rhos,
#'               "numDatasets" = numDatasets,
#'               "numSamplesPerDataset" = numSamplesPerDataset,
#'               "numFeaturesPerDataset" = numFeaturesPerDataset,
#'               "totalNumFeatures" = totalNumFeatures))             
#' }
# 
# # GRvdP
# calculatePESCAloss = function(dataTypes, Xs, Ws, Thetas, alphas, lambdas, Bs, numComponents, penalty="L2", fun_concave="gdp"){
#   
#   numDatasets = length(Xs)
#   f_obj_per_dataset = rep(0, numDatasets)
#   g_obj_per_dataset = rep(0, numDatasets)
#   Sigmas = matrix(data=0, numDatasets, numComponents)
#   penalty_fun = get(paste0("penalty_concave_",penalty))
#   
#   for(i in 1:numDatasets){
#     dataType = dataTypes[i]
#     Xl = as.matrix(Xs[[i]])
#     Wl = as.matrix(Ws[[i]])
#     Theta_l = as.matrix(Thetas[[i]])
#     alpha_l = alphas[[i]]
#     lambda_l = lambdas[[i]]
#     
#     if(is.null(dim(Bs[[i]]))){
#       
#       Bl <- t(as.matrix(Bs[[i]]))
#     }else{Bl = Bs[[i]]}
#     
#     
#     # browser()
#     # Loss function of the required data type
#     log_partition = get(paste0("log_part_", dataType))
#     
#     # browser()
#     # Calculation of first loss term
#     f_result = (1/alpha_l) * (trace_fast(Wl, log_partition(Theta_l)) - trace_fast(Theta_l, Xl))
#     f_obj_per_dataset[i] = f_result
#     
#     # Calculation of second loss term
#     g_penalty = penalty_fun(Bl, fun_concave, gamma, numComponents)   #causes error when opts$R = 6 (correct number of simulated components)
#     g_result = lambda_l * g_penalty$out
#     g_obj_per_dataset[i] = g_result
#     Sigmas[i,] = g_penalty$sigma
#   }
#   
#   f_obj = sum(f_obj_per_dataset)
#   g_obj = sum(g_obj_per_dataset)
#   loss = f_obj + g_obj
#   
#   return(list("loss" = loss,        # Loss result = f_obj + g_obj
#               "f_obj" = f_obj,      # First loss term (makes the residuals as small as possible)
#               "g_obj" = g_obj,      # Second loss term (makes the loading matrix sparse)
#               "f_obj_per_dataset" = f_obj_per_dataset,
#               "g_obj_per_dataset" = g_obj_per_dataset,
#               "Sigmas" = Sigmas))   # Sigmas
# }

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



#CHANGE

#' 
#' #from Yipeng - added as R not updating throughout CV to reflect dropping 0 columns
#' 
#' 
#' #' Split multiple data sets into training and test sets
#' #'
#' #' This function will split multiple data sets into training and test
#' #' sets. Nonmissing elements are randomly selected as the test sets.
#' #' Then the selected elements are taken as missing, and regarded as
#' #' training sets. The details can be found in \url{https://arxiv.org/abs/1902.06241}.
#' #'
#' #' @inheritParams pESCA_CV
#' #' @param ratio_mis how many percent of test set could be? default: 0.1
#' #'
#' #' @return This function returns a list contains \itemize{
#' #' \item trainSets: a list contains the training sets;
#' #' \item testSets: a list contains the test sets;
#' #' \item indexSets: a list contains the index sets.
#' #' }
#' #'
#' #' @examples
#' #' \dontrun{dataSplit(dataSets,dataTypes,ratio_mis=0.1)}
#' dataSplit <- function(dataSets, dataTypes, ratio_mis = 0.1) {
#'   # number of data sets, size of each data set
#'   nDataSets <- length(dataSets)  # number of data sets
#'   n <- rep(0, nDataSets)  # number of samples
#'   d <- rep(0, nDataSets)  # numbers of variables in different data sets
#'   for (i in 1:nDataSets) {
#'     n[i] <- dim(dataSets[[i]])[1]
#'     d[i] <- dim(dataSets[[i]])[2]
#'   }
#'   n <- n[1]
#'   
#'   # split data sets into training set and test set
#'   trainSets <- as.list(1:nDataSets)  # training set
#'   testSets <- as.list(1:nDataSets)  # test set
#'   indexSets <- as.list(1:nDataSets)  # index of the test set
#'   for (i in 1:nDataSets) {
#'     # index out the i-th data set
#'     Xi <- dataSets[[i]]
#'     dataType_Xi <- dataTypes[i]
#'     
#'     # generate the index of the test set
#'     full_ind_vec <- 1:(n * d[i])
#'     
#'     # if it is binary data, using hierachical sampling
#'     if (dataType_Xi == "B") {
#'       ones_ind_vec <- full_ind_vec[Xi == 1]
#'       zeros_ind_vec <- full_ind_vec[Xi == 0]
#'       index_Xi_ones <- sample(ones_ind_vec, round(ratio_mis * length(ones_ind_vec)))
#'       index_Xi_zeros <- sample(zeros_ind_vec, round(ratio_mis * length(zeros_ind_vec)))
#'       
#'       # test the sampled samples
#'       if (!(all(Xi[index_Xi_ones] == 1)) | !(all(Xi[index_Xi_zeros] == 0))) 
#'         message("the hierachical sampling does not work")
#'       
#'       index_Xi_test <- c(index_Xi_ones, index_Xi_zeros)
#'     } else {
#'       non_NaN_mat <- 1 - is.na(Xi)
#'       non_NaN_ind_vec <- full_ind_vec[non_NaN_mat > 0]
#'       index_Xi_test <- sample(non_NaN_ind_vec, round(ratio_mis * length(non_NaN_ind_vec)))
#'     }
#'     
#'     # generate the train set
#'     Xi_train <- Xi
#'     Xi_train[index_Xi_test] <- NA
#'     trainSets[[i]] <- Xi_train
#'     
#'     # generate the test set
#'     Xi_test <- Xi[index_Xi_test]
#'     testSets[[i]] <- Xi_test
#'     indexSets[[i]] <- index_Xi_test
#'   }
#'   
#'   # return
#'   result <- list()
#'   result$trainSets <- trainSets
#'   result$testSets <- testSets
#'   result$indexSets <- indexSets
#'   return(result)
#' }

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
#' @inheritParams pESCA
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




# 
# 
# pESCA_CV_DEV <-function(dataSets, dataTypes,
#                         lambdas_CV=NULL, penalty='L2', fun_concave='gdp', opts=list()){
#   # check if the inputs satisfy the requirements
#   stopifnot(class(dataSets) == "list")
#   stopifnot(class(penalty) == "character")
#   stopifnot(class(fun_concave) == "character")
#   if(length(dataTypes)==1){dataTypes <- unlist(strsplit(dataTypes, split=""))}
#   if(exists('quiet', where=opts)){quiet <- opts$quiet} else{quiet<-0};
#   if(exists('thr_path', where=opts)){thr_path <- opts$thr_path} else{thr_path<-0};
#   
#   # number of data sets, size of each data set
#   nTries <- length(lambdas_CV)
#   nDataSets <- length(dataSets) # number of data sets
#   n <- rep(0,nDataSets)  # number of samples
#   d <- rep(0,nDataSets)  # numbers of variables in different data sets
#   for(i in 1:nDataSets){
#     n[i] <- dim(dataSets[[i]])[1]
#     d[i] <- dim(dataSets[[i]])[2]
#   }
#   if(length(unique(as.factor(n)))!=1)
#     stop("multiple data sets have unequal sample size")
#   n <- n[1]
#   sumd <- sum(d) # total number of variables
#   
#   # default dispersion parameters alphas
#   if(exists('alphas', where=opts)){alphas<-opts$alphas} else{alphas<-rep(1,nDataSets)};
#   
#   # create zero matrix to hold results
#   cvErrors_mat <- matrix(data=0, # +1 is used for the sum of all the Xi
#                          nrow=nTries, ncol=nDataSets+1)
#   
#   # model selection process
#   opts_inner <- opts
#   
#   
#   # browser()
#   # split data sets into training set and test set
#   splitedData <- dataSplit(dataSets=dataSets,       #CHANGE
#                            dataTypes=dataTypes,
#                            ratio_mis=0.1)
#   trainSets <- splitedData$trainSets
#   
#   # browser()
#   
#   # save the parameters during the model selection
#   inits <- as.list(1:nTries)
#   if(thr_path == 1){outs <- as.list(1:nTries)}
#   
#   # model selection
#   for(j in 1:nTries){
#     lambda <- lambdas_CV[j]
#     
#     # using the training set to construct a ESCA model
#     lambdas <- lambda*rep(1,nDataSets)
#     
#     # trainModel <- pESCA(dataSets = trainSets,                #CHANGE to pESCA_ryas
#     #                     dataTypes = dataTypes,
#     #                     lambdas = lambdas,
#     #                     penalty = penalty,
#     #                     fun_concave= fun_concave,
#     #                     opts=opts_inner)
#     
#     # browser()
#     
#     trainModel <- pESCA_Ryas(dataSets = trainSets, 
#                              dataTypes = dataTypes, 
#                              lambdas = lambdas,  
#                              penalty=penalty, 
#                              fun_concave=fun_concave, 
#                              opts=opts_inner)
#     
#     if((trainModel$iter <= 2) & (quiet==0)){
#       print("less than 3 iteration is used.")
#     }
#     
#     
#     # browser()
#     
#     # warm start
#     mu <- trainModel$outcome$mu; A <- trainModel$outcome$A; B <- trainModel$outcome$B
#     opts_inner$mu0 <- mu
#     opts_inner$A0 <- A
#     opts_inner$B0 <- B
#     opts_inner$R <- dim(B)[2]         ##Addition
#     
#     inits[[j]] <- opts_inner
#     if(thr_path == 1){outs[[j]] <- trainModel$Sigmas}
#     
#     # testn <- n   #n cannot be used with browser()
#     # 
#     # browser()
#     
#     # compute the test error
#     ThetaHat <- RpESCA::ones(n) %*% mu + A %*% t(B)     #24_02_27 removed "t("mu")"
#     testError_vec <- cvError_comput(splitedData,dataTypes,alphas,ThetaHat,d)
#     
#     cvErrors_tmp <- c(sum(testError_vec), testError_vec)
#     cvErrors_mat[j,] <- cvErrors_tmp
#   }
#   
#   colnames(cvErrors_mat) <- paste0(rep("X_"), c("full", as.character(1:nDataSets)))
#   
#   result_CV <- list()
#   result_CV$cvErrors_mat <- cvErrors_mat
#   result_CV$inits <- inits
#   if(thr_path == 1){result_CV$outs <- outs}
#   
#   return(result_CV)
# }
# 




# 
# 
# 
# ##PESCA regression
# PESCA_R = function(dataSets, dataTypes, lambdas, penalty="L2", fun_concave="gdp", opts=list(), response, Wt = Wt){
#   
#   # Run zero-th iteration
#   alphas = opts$alphas
#   numComponents = opts$R
#   
#   # initialisePESCA = function(dataSets, dataTypes, lambdas, penalty="L2", opts=list(), response = NULL){
#   initialisation = initialisePESCA(dataSets, dataTypes, lambdas, penalty, opts, response = response)
#   
#   X = initialisation$X
#   W = initialisation$W
#   numSamples = initialisation$numSamplesPerDataset[1]
#   numDatasets = initialisation$numDatasets
#   numFeaturesPerDataset = initialisation$numFeaturesPerDataset
#   totalNumFeatures = initialisation$totalNumFeatures
#   rhos = initialisation$rhos
#   maxit = opts$maxit
#   tol_obj = opts$tol_obj
#   
#   # Grab the relevant update B function based on the penalty chosen.
#   update_B_fun = get(paste0("update_B_",penalty))
#   
#   #rcpp version of matrix multiplication of JHk and B (?)
#   mat_vec_matC <- function(X, d, Y) {
#     .Call('_RpESCA_mat_vec_matC', PACKAGE = 'RpESCA', X, d, Y)
#   }
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
#   # b4B <- list()
#   
#   # Obj0 <- 0       #CHANGE THIS - add into calculatepescaloss  ? make regression version?
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
#     
#     #for loop mainly for dealing with the distributions of the separate datasets
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
#     # update A
#     A_tmp = mat_vec_matC(JHk, cs^2, B_old)
#     A_svd = svd(A_tmp, nu=numComponents, nv=numComponents)
#     A_new = tcrossprod(A_svd$u, A_svd$v)
#     
#     
#     
#     
#     # update B
#     
#     #to check if inputs here are the same as original algorithm
#     #replace totalNumFeatures with a vector of features per block - numFeaturesPerDataset 
#     
#     b4B <- list(JHk, A_new, B_old, Sigmas_old, numFeaturesPerDataset, fun_concave, alphas, rhos, lambdas, gamma)
#     
#     B_new = update_B_fun(JHk, A_new, B_old, Sigmas_old, numFeaturesPerDataset, fun_concave, alphas, rhos, lambdas, gamma)
#     
#     
#     # W <- updateW_cpp(Z = Z, W = W, X = X, P = P, R = R, lambda1 = as.matrix(lambda1), 
#     #                  lambda2 = lambda2, colsumX2 = as.matrix(colsumX2), colsumP2 = as.matrix(colsumP2), blockindex = blockindex, cd = cd)
#     
#     
#     # Recalculation of Theta
#     Theta_new = RpESCA::ones(numSamples) %*% mu + tcrossprod(A_new,B_new)
#     
#     # Split up Theta and B
#     Thetas_new = splitDataset(Theta_new, numFeaturesPerDataset)
#     Bs_new = splitDataset(B_new, numFeaturesPerDataset, direction="rows")
#     
#     # Calculate loss of this iteration
#     # lossResult = calculatePESCAloss(dataTypes, initialisation$Xs, initialisation$Ws, Thetas_new, alphas, lambdas, Bs_new, numComponents, penalty)
#     #CHANGE
#     
#     print(k)
#     
#     lossResult = calculatePESCARloss(dataTypes, initialisation$Xs, initialisation$Ws, Thetas_new, alphas, lambdas, Bs_new, numComponents, penalty,  fun_concave = fun_concave,response, A_new, Wt = Wt)   
#     # lossResults = calculatePESCARloss(dataTypes, Xs, Ws, Thetas, alphas, lambdas, Bs, numComponents, penalty,  fun_concave = "gdp", response, A)
#     loss_new = lossResult$loss
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
#     
#     
#     # # Check for convergence (with r_obj addition)
#     # lossChange = (Obj0-Obj)/abs(Obj0) # relative change of loss function
#     # if((k>1) & ((lossChange < tol_obj))) break
#     # 
#     
#     # Check for convergence (OLD VERSION)
#     lossChange = (loss_old-loss_new)/abs(loss_old) # relative change of loss function
#     if((k>1) & ((lossChange < tol_obj))) break
#     
#     A_new <- lossResult$A_new
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
#     # Diagnostics
#     result$A[[k]] = A_new
#     result$B[[k]] = B_new
#     result$mu[[k]] = mu
#     result$Theta[[k]] = Theta_new
#     result$Sigmas[[k]] = Sigmas_new
#     result$loss[[k]] = lossResult
#     result$lossChange = c(result$lossChange, lossChange)
#     
#     result$b4B[[k]] <- b4B
#     
#     # if(k==1){norm_obj <- obj}
#     # Obj0    <- Obj
#   }
#   
#   result$outcome$A = A_old
#   result$outcome$B = B_old
#   result$outcome$mu = mu
#   result$outcome$loss = loss_old
#   result$outcome$Theta = Theta_old
#   result$outcome$k = k-1
#   
#   # Calculate variances explained
#   #varExpResult = calcVarExp(dataTypes, alphas, X, W, JHk, A, B, mu, numDatasets, numComponents, numSamples, numFeaturesPerDataset)
#   #result$outcome$varExp = varExpResult[[1]]
#   #result$outcome$varExpPCs = varExpResult[[2]]
#   
#   return(result)
# }

# 
# #CHANGE
# # calculatePESCAloss = function(dataTypes, Xs, Ws, Thetas, alphas, lambdas, Bs, numComponents, penalty="L2", fun_concave="gdp"){
# 
# calculatePESCARloss = function(dataTypes, Xs, Ws, Thetas, alphas, lambdas, Bs, numComponents, penalty="L2", fun_concave="gdp", response, A_new, Wt){
#   
#   numDatasets = length(Xs)
#   f_obj_per_dataset = rep(0, numDatasets)
#   g_obj_per_dataset = rep(0, numDatasets)
#   Sigmas = matrix(data=0, numDatasets, numComponents)
#   penalty_fun = get(paste0("penalty_concave_",penalty))
#   
#   for(i in 1:numDatasets){
#     dataType = dataTypes[i]
#     Xl = as.matrix(Xs[[i]])
#     Wl = Ws[[i]]
#     Theta_l = Thetas[[i]]
#     alpha_l = alphas[[i]]
#     lambda_l = lambdas[[i]]
#     Bl = Bs[[i]]
#     
#     # Loss function of the required data type
#     log_partition = get(paste0("log_part_", dataType))
#     
#     # Calculation of first loss term
#     f_result = (1/alpha_l) * (trace_fast(Wl, log_partition(Theta_l)) - trace_fast(Theta_l, Xl))
#     f_obj_per_dataset[i] = f_result
#     
#     # Calculation of second loss term
#     g_penalty = penalty_fun(Bl, fun_concave, gamma, numComponents)   
#     g_result = lambda_l * g_penalty$out
#     g_obj_per_dataset[i] = g_result
#     Sigmas[i,] = g_penalty$sigma
#   }
#   
#   f_obj = sum(f_obj_per_dataset)
#   g_obj = sum(g_obj_per_dataset)
#   loss = f_obj + g_obj
#   
#   
#   
#   #CHANGE
#   # #regression of scores A onto a response
#   # # make some cross-validation approach similar to https://github.com/enorthrop/sup.r.jive/blob/master/R/sJIVE.R line 809-844?
#   # # how to improve this?
#   
#   
#   regres <- glm(t(response) ~ A_new)
#   
#   #use model deviance as regression objective for now
#   r_obj <- regres$deviance
#   
#   
#   #recalculate deviance manually:
#   
#   
#   
#   #residual sum of squares - equivalent for gaussian data
#   # r_obj <- 
#   
#   
#   # A_new <- A_new %*% diag(regres$coefficients[-c(1)])    #iterative regression coefficient weighted scores CHANGE: tie this to the Wt parameter
#   
#   A_distorted <-  A_new %*% diag(regres$coefficients[-c(1)])      #diag(regres$coefficients[-c(1)]) %*% solve(diag(regres$coefficients[-c(1)]))
#   
#   A_diff <- (A_distorted - A_new)
#   
#   A_new <- A_new + (Wt * A_diff)
#   
#   # A_new <- A_new %*% (diag(regres$coefficients[-c(1)]) %*% solve(diag(regres$coefficients[-c(1)])))
#   
#   
#   
#   # CHANGE
#   # obj <- obj/norm_obj
#   #
#   # weight obj and r_obj by a parameter; Wt in [0,1] and 1-Wt respectively
#   # Obj <- (1-Wt * loss) + ((1-Wt) * r_obj)
#   
#   Obj <- ((1-Wt) * loss) + (Wt * r_obj)
#   
#   
#   return(list("loss" = Obj,        # Loss result = see above CHANGE
#               "f_obj" = f_obj,      # First loss term (makes the residuals as small as possible)
#               "g_obj" = g_obj,      # Second loss term (makes the loading matrix sparse)
#               "r_obj" = r_obj,
#               "f_obj_per_dataset" = f_obj_per_dataset,
#               "g_obj_per_dataset" = g_obj_per_dataset,
#               "Sigmas" = Sigmas,
#               "A_new" = A_new,
#               "glm_model" = regres))   # Sigmas
# }
# 


# 
# ##PESCA regression
# PESCA_R = function(dataSets, dataTypes, lambdas, penalty="L2", fun_concave="gdp", opts=list(), response, Wt = Wt){
#   
#   # Run zero-th iteration
#   alphas = opts$alphas
#   numComponents = opts$R
#   
#   # initialisePESCA = function(dataSets, dataTypes, lambdas, penalty="L2", opts=list(), response = NULL){
#   initialisation = initialisePESCA(dataSets, dataTypes, lambdas, penalty, opts, response = response)
#   
#   X = initialisation$X
#   W = initialisation$W
#   numSamples = initialisation$numSamplesPerDataset[1]
#   numDatasets = initialisation$numDatasets
#   numFeaturesPerDataset = initialisation$numFeaturesPerDataset
#   totalNumFeatures = initialisation$totalNumFeatures
#   rhos = initialisation$rhos
#   maxit = opts$maxit
#   tol_obj = opts$tol_obj
#   
#   # Grab the relevant update B function based on the penalty chosen.
#   update_B_fun = get(paste0("update_B_",penalty))
#   
#   #rcpp version of matrix multiplication of JHk and B (?)
#   mat_vec_matC <- function(X, d, Y) {
#     .Call('_RpESCA_mat_vec_matC', PACKAGE = 'RpESCA', X, d, Y)
#   }
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
#   # b4B <- list()
#   
#   # Obj0 <- 0       #CHANGE THIS - add into calculatepescaloss  ? make regression version?
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
#     
#     #for loop mainly for dealing with the distributions of the separate datasets
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
#     # update A
#     A_tmp = mat_vec_matC(JHk, cs^2, B_old)
#     A_svd = svd(A_tmp, nu=numComponents, nv=numComponents)
#     A_new = tcrossprod(A_svd$u, A_svd$v)
#     
#     
#     
#     
#     # update B
#     
#     #to check if inputs here are the same as original algorithm
#     #replace totalNumFeatures with a vector of features per block - numFeaturesPerDataset 
#     
#     b4B <- list(JHk, A_new, B_old, Sigmas_old, numFeaturesPerDataset, fun_concave, alphas, rhos, lambdas, gamma)
#     
#     B_new = update_B_fun(JHk, A_new, B_old, Sigmas_old, numFeaturesPerDataset, fun_concave, alphas, rhos, lambdas, gamma)
#     
#     
#     # W <- updateW_cpp(Z = Z, W = W, X = X, P = P, R = R, lambda1 = as.matrix(lambda1), 
#     #                  lambda2 = lambda2, colsumX2 = as.matrix(colsumX2), colsumP2 = as.matrix(colsumP2), blockindex = blockindex, cd = cd)
#     
#     
#     # Recalculation of Theta
#     Theta_new = RpESCA::ones(numSamples) %*% mu + tcrossprod(A_new,B_new)
#     
#     # Split up Theta and B
#     Thetas_new = splitDataset(Theta_new, numFeaturesPerDataset)
#     Bs_new = splitDataset(B_new, numFeaturesPerDataset, direction="rows")
#     
#     # Calculate loss of this iteration
#     # lossResult = calculatePESCAloss(dataTypes, initialisation$Xs, initialisation$Ws, Thetas_new, alphas, lambdas, Bs_new, numComponents, penalty)
#     #CHANGE
#     
#     print(k)
#     
#     lossResult = calculatePESCARloss(dataTypes, initialisation$Xs, initialisation$Ws, Thetas_new, alphas, lambdas, Bs_new, numComponents, penalty,  fun_concave = fun_concave,response, A_new, Wt = Wt)   
#     # lossResults = calculatePESCARloss(dataTypes, Xs, Ws, Thetas, alphas, lambdas, Bs, numComponents, penalty,  fun_concave = "gdp", response, A)
#     loss_new = lossResult$loss
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
#     
#     
#     # # Check for convergence (with r_obj addition)
#     # lossChange = (Obj0-Obj)/abs(Obj0) # relative change of loss function
#     # if((k>1) & ((lossChange < tol_obj))) break
#     # 
#     
#     # Check for convergence (OLD VERSION)
#     lossChange = (loss_old-loss_new)/abs(loss_old) # relative change of loss function
#     if((k>1) & ((lossChange < tol_obj))) break
#     
#     A_new <- lossResult$A_new
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
#     # Diagnostics
#     result$A[[k]] = A_new
#     result$B[[k]] = B_new
#     result$mu[[k]] = mu
#     result$Theta[[k]] = Theta_new
#     result$Sigmas[[k]] = Sigmas_new
#     result$loss[[k]] = lossResult
#     result$lossChange = c(result$lossChange, lossChange)
#     
#     result$b4B[[k]] <- b4B
#     
#     # if(k==1){norm_obj <- obj}
#     # Obj0    <- Obj
#   }
#   
#   result$outcome$A = A_old
#   result$outcome$B = B_old
#   result$outcome$mu = mu
#   result$outcome$loss = loss_old
#   result$outcome$Theta = Theta_old
#   result$outcome$k = k-1
#   
#   # Calculate variances explained
#   #varExpResult = calcVarExp(dataTypes, alphas, X, W, JHk, A, B, mu, numDatasets, numComponents, numSamples, numFeaturesPerDataset)
#   #result$outcome$varExp = varExpResult[[1]]
#   #result$outcome$varExpPCs = varExpResult[[2]]
#   
#   return(result)
# }
# 
# 
# #CHANGE
# # calculatePESCAloss = function(dataTypes, Xs, Ws, Thetas, alphas, lambdas, Bs, numComponents, penalty="L2", fun_concave="gdp"){
# 
# calculatePESCARloss = function(dataTypes, Xs, Ws, Thetas, alphas, lambdas, Bs, numComponents, penalty="L2", fun_concave="gdp", response, A_new, Wt){
#   
#   numDatasets = length(Xs)
#   f_obj_per_dataset = rep(0, numDatasets)
#   g_obj_per_dataset = rep(0, numDatasets)
#   Sigmas = matrix(data=0, numDatasets, numComponents)
#   penalty_fun = get(paste0("penalty_concave_",penalty))
#   
#   for(i in 1:numDatasets){
#     dataType = dataTypes[i]
#     Xl = as.matrix(Xs[[i]])
#     Wl = Ws[[i]]
#     Theta_l = Thetas[[i]]
#     alpha_l = alphas[[i]]
#     lambda_l = lambdas[[i]]
#     Bl = Bs[[i]]
#     
#     # Loss function of the required data type
#     log_partition = get(paste0("log_part_", dataType))
#     
#     # Calculation of first loss term
#     f_result = (1/alpha_l) * (trace_fast(Wl, log_partition(Theta_l)) - trace_fast(Theta_l, Xl))
#     f_obj_per_dataset[i] = f_result
#     
#     # Calculation of second loss term
#     g_penalty = penalty_fun(Bl, fun_concave, gamma, numComponents)   
#     g_result = lambda_l * g_penalty$out
#     g_obj_per_dataset[i] = g_result
#     Sigmas[i,] = g_penalty$sigma
#   }
#   
#   f_obj = sum(f_obj_per_dataset)
#   g_obj = sum(g_obj_per_dataset)
#   loss = f_obj + g_obj
#   
#   
#   
#   #CHANGE
#   # #regression of scores A onto a response
#   # # make some cross-validation approach similar to https://github.com/enorthrop/sup.r.jive/blob/master/R/sJIVE.R line 809-844?
#   # # how to improve this?
#   
#   
#   regres <- glm(t(response) ~ A_new)
#   
#   #use model deviance as regression objective for now
#   r_obj <- regres$deviance
#   
#   
#   #recalculate deviance manually:
#   
#   
#   
#   #residual sum of squares - equivalent for gaussian data
#   # r_obj <- 
#   
#   
#   # A_new <- A_new %*% diag(regres$coefficients[-c(1)])    #iterative regression coefficient weighted scores CHANGE: tie this to the Wt parameter
#   
#   A_distorted <-  A_new %*% diag(regres$coefficients[-c(1)])      #diag(regres$coefficients[-c(1)]) %*% solve(diag(regres$coefficients[-c(1)]))
#   
#   A_diff <- (A_distorted - A_new)
#   
#   A_new <- A_new + (Wt * A_diff)
#   
#   # A_new <- A_new %*% (diag(regres$coefficients[-c(1)]) %*% solve(diag(regres$coefficients[-c(1)])))
#   
#   
#   
#   # CHANGE
#   # obj <- obj/norm_obj
#   #
#   # weight obj and r_obj by a parameter; Wt in [0,1] and 1-Wt respectively
#   # Obj <- (1-Wt * loss) + ((1-Wt) * r_obj)
#   
#   Obj <- ((1-Wt) * loss) + (Wt * r_obj)
#   
#   
#   return(list("loss" = Obj,        # Loss result = see above CHANGE
#               "f_obj" = f_obj,      # First loss term (makes the residuals as small as possible)
#               "g_obj" = g_obj,      # Second loss term (makes the loading matrix sparse)
#               "r_obj" = r_obj,
#               "f_obj_per_dataset" = f_obj_per_dataset,
#               "g_obj_per_dataset" = g_obj_per_dataset,
#               "Sigmas" = Sigmas,
#               "A_new" = A_new,
#               "glm_model" = regres))   # Sigmas
# }
# 
# 


# 
# ##PESCA regression
# PESCA_Rv2 = function(dataSets, dataTypes, lambdas, penalty="L2", fun_concave="gdp", opts=list(), response, Wt = Wt){
#   
#   # Run zero-th iteration
#   alphas = opts$alphas
#   numComponents = opts$R
#   
#   # initialisePESCA = function(dataSets, dataTypes, lambdas, penalty="L2", opts=list(), response = NULL){
#   initialisation = initialisePESCA(dataSets, dataTypes, lambdas, penalty, opts, response = response)
#   
#   X = initialisation$X
#   W = initialisation$W
#   numSamples = initialisation$numSamplesPerDataset[1]
#   numDatasets = initialisation$numDatasets
#   numFeaturesPerDataset = initialisation$numFeaturesPerDataset
#   totalNumFeatures = initialisation$totalNumFeatures
#   rhos = initialisation$rhos
#   maxit = opts$maxit
#   tol_obj = opts$tol_obj
#   
#   # Grab the relevant update B function based on the penalty chosen.
#   update_B_fun = get(paste0("update_B_",penalty))
#   
#   #rcpp version of matrix multiplication of JHk and B (?)
#   mat_vec_matC <- function(X, d, Y) {
#     .Call('_RpESCA_mat_vec_matC', PACKAGE = 'RpESCA', X, d, Y)
#   }
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
#   #add initialisation by PCoVR here
#   
#   # b4B <- list()
#   
#   # Obj0 <- 0       #CHANGE THIS - add into calculatepescaloss  ? make regression version?
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
#     
#     #for loop mainly for dealing with the distributions of the separate datasets
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
#     # update A
#     A_tmp = mat_vec_matC(JHk, cs^2, B_old)
#     A_svd = svd(A_tmp, nu=numComponents, nv=numComponents)
#     A_new = tcrossprod(A_svd$u, A_svd$v)
#     
#     
#     
#     
#     # update B
#     
#     #to check if inputs here are the same as original algorithm
#     #replace totalNumFeatures with a vector of features per block - numFeaturesPerDataset 
#     
#     b4B <- list(JHk, A_new, B_old, Sigmas_old, numFeaturesPerDataset, fun_concave, alphas, rhos, lambdas, gamma)
#     
#     B_new = update_B_fun(JHk, A_new, B_old, Sigmas_old, numFeaturesPerDataset, fun_concave, alphas, rhos, lambdas, gamma)
#     
#     
#     # W <- updateW_cpp(Z = Z, W = W, X = X, P = P, R = R, lambda1 = as.matrix(lambda1), 
#     #                  lambda2 = lambda2, colsumX2 = as.matrix(colsumX2), colsumP2 = as.matrix(colsumP2), blockindex = blockindex, cd = cd)
#     
#     
#     # Recalculation of Theta
#     Theta_new = RpESCA::ones(numSamples) %*% mu + tcrossprod(A_new,B_new)
#     
#     # Split up Theta and B
#     Thetas_new = splitDataset(Theta_new, numFeaturesPerDataset)
#     Bs_new = splitDataset(B_new, numFeaturesPerDataset, direction="rows")
#     
#     # Calculate loss of this iteration
#     # lossResult = calculatePESCAloss(dataTypes, initialisation$Xs, initialisation$Ws, Thetas_new, alphas, lambdas, Bs_new, numComponents, penalty)
#     #CHANGE
#     
#     print(k)
#     
#     lossResult = calculatePESCARloss(dataTypes, initialisation$Xs, initialisation$Ws, Thetas_new, alphas, lambdas, Bs_new, numComponents, penalty,  fun_concave = fun_concave,response, A_new, Wt = Wt)   
#     # lossResults = calculatePESCARloss(dataTypes, Xs, Ws, Thetas, alphas, lambdas, Bs, numComponents, penalty,  fun_concave = "gdp", response, A)
#     loss_new = lossResult$loss
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
#     
#     
#     # # Check for convergence (with r_obj addition)
#     # lossChange = (Obj0-Obj)/abs(Obj0) # relative change of loss function
#     # if((k>1) & ((lossChange < tol_obj))) break
#     # 
#     
#     # Check for convergence (OLD VERSION)
#     lossChange = (loss_old-loss_new)/abs(loss_old) # relative change of loss function
#     if((k>1) & ((lossChange < tol_obj))) break
#     
#     A_new <- lossResult$A_new
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
#     # Diagnostics
#     result$A[[k]] = A_new
#     result$B[[k]] = B_new
#     result$mu[[k]] = mu
#     result$Theta[[k]] = Theta_new
#     result$Sigmas[[k]] = Sigmas_new
#     result$loss[[k]] = lossResult
#     result$lossChange = c(result$lossChange, lossChange)
#     
#     result$b4B[[k]] <- b4B
#     
#     # if(k==1){norm_obj <- obj}
#     # Obj0    <- Obj
#   }
#   
#   result$outcome$A = A_old
#   result$outcome$B = B_old
#   result$outcome$mu = mu
#   result$outcome$loss = loss_old
#   result$outcome$Theta = Theta_old
#   result$outcome$k = k-1
#   
#   # Calculate variances explained
#   #varExpResult = calcVarExp(dataTypes, alphas, X, W, JHk, A, B, mu, numDatasets, numComponents, numSamples, numFeaturesPerDataset)
#   #result$outcome$varExp = varExpResult[[1]]
#   #result$outcome$varExpPCs = varExpResult[[2]]
#   
#   return(result)
# }
# 




# 
# 
# ######y as block
# pESCA_Ryas = function(dataSets, dataTypes, lambdas, penalty="L2", fun_concave="gdp", opts=list()){
#   
#   # Run zero-th iteration
#   alphas = opts$alphas
#   numComponents = opts$R
#   
#   # browser()
#   initialisation = initialisePESCA(dataSets = dataSets, 
#                                    dataTypes = dataTypes, 
#                                    lambdas = lambdas, penalty = penalty, 
#                                    opts = opts)
#   
#   X = initialisation$X
#   W = initialisation$W
#   numSamples = initialisation$numSamplesPerDataset[1]
#   numDatasets = initialisation$numDatasets
#   numFeaturesPerDataset = initialisation$numFeaturesPerDataset
#   totalNumFeatures = initialisation$totalNumFeatures
#   rhos = initialisation$rhos
#   maxit = opts$maxit
#   tol_obj = opts$tol_obj
#   
#   # Grab the relevant update B function based on the penalty chosen.
#   update_B_fun = get(paste0("update_B_",penalty))
#   
#   #rcpp version of matrix multiplication of JHk and B (?)
#   mat_vec_matC <- function(X, d, Y) {
#     .Call('_RpESCA_mat_vec_matC', PACKAGE = 'RpESCA', X, d, Y)
#   }
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
#   # b4B <- list()
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
#     # update A
#     A_tmp = mat_vec_matC(JHk, cs^2, B_old)
#     A_svd = svd(A_tmp, nu=numComponents, nv=numComponents)
#     A_new = tcrossprod(A_svd$u, A_svd$v)
#     
#     # update B
#     
#     #to check if inputs here are the same as original algorithm
#     #replace totalNumFeatures with a vector of features per block - numFeaturesPerDataset 
#     
#     b4B <- list(JHk, A_new, B_old, Sigmas_old, numFeaturesPerDataset, fun_concave, alphas, rhos, lambdas, gamma)
#     
#     B_new = update_B_fun(JHk, A_new, B_old, Sigmas_old, numFeaturesPerDataset, fun_concave, alphas, rhos, lambdas, gamma)
#     
#     # Recalculation of Theta
#     Theta_new = RpESCA::ones(numSamples) %*% mu + tcrossprod(A_new,B_new)
#     
#     # Split up Theta and B
#     Thetas_new = splitDataset(Theta_new, numFeaturesPerDataset)
#     Bs_new = splitDataset(B_new, numFeaturesPerDataset, direction="rows")
#     
#     # Calculate loss of this iteration
#     lossResult = calculatePESCAloss(dataTypes = dataTypes, 
#                                     Xs = initialisation$Xs, 
#                                     Ws = initialisation$Ws, 
#                                     Thetas = Thetas_new, 
#                                     alphas = alphas, 
#                                     lambdas = lambdas, 
#                                     Bs = Bs_new, 
#                                     numComponents = numComponents, 
#                                     penalty = penalty)
#     # calculatePESCAloss = function(dataTypes, Xs, Ws, Thetas, alphas, lambdas, Bs, numComponents, penalty="L2", fun_concave="gdp"){
#     
#     
#     
#     loss_new = lossResult$loss
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
#     
#     # Check for convergence
#     lossChange = (loss_old-loss_new)/abs(loss_old) # relative change of loss function
#     if((k>1) & ((lossChange < tol_obj))) break
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
#     result$b4B[[k]] <- b4B
#   }
#   
#   result$outcome$A = A_old
#   result$outcome$B = B_old
#   result$outcome$mu = mu
#   result$outcome$loss = loss_old
#   result$outcome$Theta = Theta_old
#   result$outcome$k = k-1
#   
#   # Calculate variances explained
#   varExpResult = calcVarExp(dataTypes, alphas, X, W, JHk, A_old, B_old, mu, numDatasets, numComponents, numSamples, numFeaturesPerDataset)
#   result$outcome$varExp = varExpResult[[1]]
#   result$outcome$varExpPCs = varExpResult[[2]]
#   
#   return(result)
# }
# 
