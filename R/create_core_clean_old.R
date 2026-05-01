maxr <- function(n){
  rep(1,n)
}

minr <- function(n){
  rep(0,n)
}


#CHANGE THIS - max value shouldn't equal 1
fadeup <- function(n){
  if(n != 0){
    c(1:n)/n
  }
}

#CHANGE THIS - max value shouldn't equal 1
fadedown <- function(n){
  if(n != 0){
    c(n:1)/n
  }
}

generateErrorMatrix = function(num_samples, num_features, effectSize=0.05, mean=0, sd=0, seed=123){
  set.seed(seed)
  E = matrix(rnorm(num_samples*num_features, mean, sd), nrow=num_samples, ncol=num_features)
  E_norm = E / norm(E, type="F")
  E_norm_effect = effectSize * E_norm
  return(E_norm_effect)
}

generateLoadingVector = function(n=c(10,10,10), totalLength=30, options=c("low", "fadeup", "high")){
  P = c()
  
  for(i in 1:length(n)){
    numberOfElements = n[i]
    option = options[i]
    
    if (option == "low"){
      newPart = minr(numberOfElements)
    }
    if (option == "fadeup"){
      newPart = fadeup(numberOfElements)
    }
    if (option == "fadedown"){
      newPart = fadedown(numberOfElements)
    }
    if (option == "high"){
      newPart = maxr(numberOfElements)
    }
    
    P = c(P, newPart)
  }
  
  P = c(P, rep(0, totalLength - length(P)))
  return(P)
}

generateLoadingMatrix = function(numFeatures, totalLength, options){
  loadingMatrix = c()
  numVectors = ncol(options)
  for(i in 1:numVectors){
    v = generateLoadingVector(numFeatures, totalLength, options[,i])
    loadingMatrix = cbind(loadingMatrix,v)
  }
  return(loadingMatrix)
}

Create_Core_PESCA_4 <- function(nreps=3,
                                meta="",
                                a_sigma = c(1.5,0.75,0.6),
                                b_sigma = c(1,0.8),
                                e_sigma = c(1,0.8,0.5),
                                noise_sd = c(1,1,1),
                                EffectSize = c(X_a_ab = 1, e1 = 0.5, e2 = 0.25, e3 = 0.1, mu = 1, Struc_E = 1),
                                genes_tot = 750,
                                chems_tot = 500,
                                bact_tot = 250,
                                numFeatures_chems = c(10, 10, 10, 10, 10),
                                numFeatures_genes = c(5, 5, 5, 5, 5),
                                numFeatures_bact = c(5, 5, 5, 5, 5),
                                SCORE_SEED = 1000,
                                Experiment_responders = 12,
                                ts = 1234,
                                struc_seed = 127,
                                E_seed=2){
  
  if(is.null(meta) | is.null(dim(meta))){
    meta <- cbind.data.frame(rep(c(1,-1), each = 12), rep(c(1:4), each = 3), rep(c(-1:-4,1:4), each = 3))
    colnames(meta) <- c("growth_condition","time","interaction")
  }
  
  targets_a_ab <- rep(0,12)
  # unclear what this contains
  
  
  
  # svded matrix of 4 rows and 3 columns hardcoded
  # 3 columns = 3 components
  set.seed(SCORE_SEED)
  scoreMatrix <- data.frame(svd(matrix(rnorm(12), nrow = 4))$u)
  
  # weird effort to create n replicates
  # also times -1 for some reason
  scoreMatrix <- scoreMatrix[rep(rownames(scoreMatrix), each = nreps),]
  scoreMatrix <- rbind.data.frame(scoreMatrix, scoreMatrix*-1)
  
  P_chems = generateLoadingMatrix(numFeatures=numFeatures_chems,
                                  totalLength=chems_tot,
                                  options=matrix(c("high","fadedown","low","low","low",
                                                   "low","fadeup","high","fadedown","low",
                                                   "low","low","low","fadeup","high"), ncol=3, nrow=5))
  rownames(P_chems) <- paste0("C_",1:nrow(P_chems))
  
  prependMatrix = matrix(rep(0, sum(numFeatures_chems)*3), nrow=sum(numFeatures_chems), ncol=3)
  P_bact = generateLoadingMatrix(numFeatures=numFeatures_bact,
                                 totalLength=bact_tot,
                                 options=matrix(c("low","low","low","fadeup","high",
                                                  "low","fadeup","high","fadedown","low",
                                                  "low","low","low","low","low"), ncol=3, nrow=5))
  # P_bact = rbind(prependMatrix, P_bact)
  rownames(P_bact) <- paste0("B_",1:nrow(P_bact))
  
  prependMatrix2 = matrix(rep(0, sum(numFeatures_bact)*3), nrow=sum(numFeatures_bact), ncol=3)
  P_genes = generateLoadingMatrix(numFeatures=numFeatures_genes,
                                  totalLength=genes_tot,
                                  options=matrix(c("low","low","low","low","low",
                                                   "low","fadeup","high","fadedown","low",
                                                   "high","fadedown","low","low","low"), ncol=3, nrow=5))
  # P_genes = rbind(prependMatrix, prependMatrix2, P_genes)
  rownames(P_genes) <- paste0("G_",1:nrow(P_genes))
  # This solves a PESCA-problem that cannot be dealt with in another way.
  # P_genes = P_genes[nrow(P_genes):1,]
  # P_genes = P_genes[,ncol(P_genes):1]
  
  loadingMatrix <- rbind(P_genes,P_chems,P_bact)
  
  #response
  # hardcoded response
  y <- c(2,1,0.7) %*% t(scoreMatrix)
  y <- y + rnorm(length(y), sd = 0.5)
  
  # P <- cbind(a_ab_p1, a_ab_p2, a_ab_p3, a_ab_p4)
  # a sigma has relative sizes of the components
  # gets put into P
  loadingMatrix_with_size <- loadingMatrix %*% diag(a_sigma)
  
  # X_scoreMatrix gets overwritten here
  X_no_noise <- data.matrix(scoreMatrix) %*% t(loadingMatrix_with_size)
  X_no_noise <- EffectSize["X_a_ab"] * (X_no_noise/norm(X_no_noise, type = "F"))
  
  ##CREATE FEATURE MEAN VALUES and ensure they are above 0
  #set.seed here?
  # add means that are normally distributed
  mu <- EffectSize["mu"] * matrix(rep(rnorm(dim(X_no_noise)[2], mean = mean(X_no_noise)), each = nrow(meta)), nrow = nrow(meta))
  mu <- mu + abs(min(mu))
  X_no_noise_mu <- X_no_noise + mu
  
  ####CREATE E
  #seems redundant but isn't due to prependMatrix requirement
  genes_tot <- length(grep("^G",rownames(loadingMatrix)))
  chems_tot <-  length(grep("^C",rownames(loadingMatrix)))
  bact_tot<-  length(grep("^B",rownames(loadingMatrix)))
  
  ## some random noise - i.e. differences between replicates
  num_samples = nrow(scoreMatrix)
  E1 = generateErrorMatrix(num_samples, genes_tot, mean=0, sd=noise_sd[1], seed=E_seed+20000)
  E2 = generateErrorMatrix(num_samples, chems_tot, mean=0, sd=noise_sd[2], seed=E_seed+40000)
  E3 = generateErrorMatrix(num_samples, bact_tot, mean=0, sd=noise_sd[3], seed=E_seed+60000)
  
  ##ADD ALL TOGETHER TO CREATE "CORE"
  X <- matrix(nrow = nrow(X_no_noise_mu), ncol = ncol(X_no_noise_mu))
  
  X[,1:genes_tot] <- X_no_noise_mu[,1:genes_tot] + E1
  X[,(genes_tot+1):(genes_tot+chems_tot)] <- X_no_noise_mu[,(genes_tot+1):(genes_tot+chems_tot)] + E2
  X[,(genes_tot+chems_tot+1):(genes_tot+chems_tot+bact_tot)] <- X_no_noise_mu[,(genes_tot+chems_tot+1):(genes_tot+chems_tot+bact_tot)] + E3
  
  # make all values positive
  X <- t(((t(X) + abs(colMins(X)))))
  colnames(X) <- paste0("X_",c(1:dim(X)[2]))
  
  #change this to get IDs for the different types of variation
  pathway <- which(rowSums(loadingMatrix) != 0)
  pathway <- paste0("X_",pathway)
  
  DataBlocks = list(X[,1:genes_tot], X[,(genes_tot+1):(genes_tot+chems_tot)], X[,(genes_tot+chems_tot+1):(genes_tot+chems_tot+bact_tot)])
  Residuals = list(E1, E2, E3)
  
  return(list("Data" = DataBlocks,
              "scoreMatrix" = scoreMatrix,
              "loadingMatrix" = loadingMatrix,
              "sigma" = c(a_sigma, b_sigma, e_sigma),
              "residuals" = Residuals,
              "response" = y,
              "pathway_ids" = pathway))  #Effect inputs here?..
}


