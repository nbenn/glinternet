glinternet = function(X, Y, numLevels, lambda=NULL, nLambda=50, lambdaMinRatio=0.01, screenLimit=NULL, numToFind=NULL, family=c("gaussian", "binomial"), tol=1e-5, maxIter=5000, verbose=FALSE, numCores=1, filename=NA){

                                        #get call and family
  thisCall = match.call()
  family = match.arg(family)
  
                                        #make sure inputs are valid
  n = length(Y)
  pCat = sum(numLevels > 1)
  pCont = length(numLevels) - pCat
  stopifnot(n==nrow(X), pCat+pCont==ncol(X), family=="gaussian"||family=="binomial")
  if (family=="binomial" && !all(Y %in% 0:1)) stop("Error:family=binomial but Y not in {0,1}")

                                        #separate into categorical and continuous parts
  if (pCont > 0) {
    #Z = as.matrix(apply(as.matrix(X[, numLevels == 1]), 2, standardize))
    Z = .Call("alloc_z", n, pCont, X)
    for (i in 1:pCont) {
      vec = .Call("extract_col", X, i)
      if (length(unique(vec)) == 1) {
        result = rep(0, length(vec))
      } else {
        result = vec - mean(vec)
        result = result / sqrt(t(result)%*%result)
      }
      Z = .Call("import_col", result, Z, i)
    }
  }
  else Z = NULL
  if (pCat > 0){
    catIndices = which(numLevels > 1)
    levels = numLevels[catIndices]
    Xcat = as.matrix(X[, catIndices])
  }
  else {
    levels = NULL
    Xcat = NULL
  }
  
                                        #compute variable norms
  res = .Call("alloc_res", Y)
  #res = Y - mean(Y)
  #res = .Call("copy_vec", Y - mean(Y), res)
  candidates = get_candidates(Xcat, Z, res, n, pCat, pCont, levels, screenLimit, verbose, numCores=numCores)
 
 
                                        #lambda grid if not user provided
  if (is.null(lambda)) lambda = get_lambda_grid(candidates, nLambda, lambdaMinRatio)
  else {
    stopifnot(min(lambda) > 0)
    if (any(diff(lambda) > 0)) stop("Error: input lambda sequence is not monotone decreasing.")
    lambdaMax = max(get_lambda_grid(candidates, nLambda, lambdaMinRatio))
    nLambda = length(lambda)
    if (nLambda == 1){
      lambda = sort(c(lambda, lambdaMax), decreasing=TRUE)
      nLambda = 2
    }
  }

  if (verbose) {
    cat("lambda sequence:\n", lambda, "\n\n")
  }
  
                                        #initialize storage for results
  fitted = matrix(mean(Y), n, nLambda)
  activeSet = vector("list", nLambda)
  betahat = vector("list", nLambda)
  betahat[[1]] = ifelse(family=="gaussian", mean(Y), -log(1/mean(Y)-1))
  objValue = rep(0, nLambda)
  objValue[1] = ifelse(family=="gaussian", sum(res^2)/(2*n), -mean(Y)*betahat[[1]]+log(1/(1-mean(Y))))

                                        #ever-active set + sequential strong rules + group lasso
  for (i in 2:nLambda){
    if (verbose) {
      cat("-> lambda ", i, ": ", lambda[i], "\n")
      time.lambda <- proc.time()
    }
    activeSet[[i]] = strong_rules(candidates, lambda[i], lambda[i-1])
    betahat[[i]] = initialize_betahat(activeSet[[i]], activeSet[[i-1]], betahat[[i-1]], levels)
    while (TRUE){
      #group lasso on strong set
      if (verbose) time.gl <- proc.time()  
      solution = group_lasso(Xcat, Z, Y, activeSet[[i]], betahat[[i]], levels, lambda[i], family, tol, maxIter, verbose)
      if (verbose) {
        time.gl <- proc.time() - time.gl
        cat("--> time group lasso: ", time.gl[1], "(user)\n",
          rep(" ", 11), time.gl[3], "(elapsed)\n")
      }
      activeSet[[i]] = solution$activeSet
      betahat[[i]] = solution$betahat
      res = .Call("copy_vec", solution$res, res)
      #res = solution$res
      objValue[i] = solution$objValue
      #check kkt conditions on the rest
      if (verbose) time.kkt <- proc.time()  
      check = check_kkt(Xcat, Z, res, n, pCat, pCont, levels, candidates, activeSet[[i]], lambda[i], verbose, numCores)
      if (verbose) {
        time.kkt <- proc.time() - time.kkt
        cat("--> time check kkt: ", time.kkt[1], "(user)\n",
          rep(" ", 10), time.kkt[3], "(elapsed)\n")
      }
      candidates$norms = check$norms
      if (check$flag) break
      betahat[[i]] = initialize_betahat(check$activeSet, activeSet[[i]], betahat[[i]], levels)
      activeSet[[i]] = check$activeSet
    }
    #update the candidate set if necessary
    if (!is.null(screenLimit) && (screenLimit<pCat+pCont) && i<nLambda) {
      if (verbose) time.candidates <- proc.time()  
      candidates = get_candidates(Xcat, Z, res, n, pCat, pCont, levels, screenLimit, activeSet[[i]], candidates$norms, verbose, numCores)
      if (verbose) {
        time.candidates <- proc.time() - time.candidates
        cat("-> time candidates: ", time.candidates[1], "(user)\n",
          rep(" ", 10), time.candidates[3], "(elapsed)\n")
      }
    }
    #get fitted values
    fitted[, i] = Y - res
    #compute total number of interactions found
    if (!is.null(numToFind)){
      numFound = sum(sapply(activeSet[[i]][3:5], function(x) ifelse(is.null(x), 0, nrow(x))))
      if (numFound >= numToFind) break
    }
    if (verbose & !is.na(filename)) {
      saveRDS(activeSet, paste(filename, "-activeSet.rds", sep=""))
    }
    if (verbose) {
      time.lambda <- proc.time() - time.lambda
      cat("-> total time: ", time.lambda[1], "(user)\n",
        rep(" ", 7), time.lambda[3], "(elapsed)\n\n")
    }
  }

  #rescale betahat
  Z = as.matrix(X[, numLevels==1])
  betahatRescaled = lapply(1:i, function(j) rescale_betahat(activeSet[[j]], betahat[[j]], Xcat, Z, levels, n))

  output = list(call=thisCall, fitted=fitted[, 1:i], lambda=lambda[1:i], objValue=objValue, activeSet=activeSet[1:i], betahat=betahatRescaled[1:i], numLevels=numLevels, family=family)
  class(output) = "glinternet"
  return (output)
}
