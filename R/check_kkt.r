check_kkt = function(X, Z, res, n, pCat, pCont, numLevels, candidates, activeSet, lambda, cpuInfo, verbose=FALSE){

  norms = vector("list", 5)
  names(norms) = c("cat", "cont", "catcat", "contcont", "catcont")

  mynode = .Call("get_my_numa_node")

  #compute the norms
  if (pCat > 0){
    norms$cat = compute_norms_cat(X, res, n, pCat, numLevels, cpuInfo$num_used_cpus)
    if (!is.null(candidates$variables$catcat)) norms$catcat = compute_norms_cat_cat(X, res, n, numLevels, candidates$variables$catcat, cpuInfo$num_used_cpus)
  }
  if (pCont > 0){
    norms$cont = compute_norms_cont(Z[[mynode+1]], res, n)
    if (!is.null(candidates$variables$contcont)) norms$contcont = compute_norms_cont_cont(Z, norms$cont, res, n, candidates$variables$contcont, cpuInfo, verbose)
  }
  if (!is.null(candidates$variables$catcont)) norms$catcont = compute_norms_cat_cont(X, Z[[mynode+1]], norms$cat, res, n, numLevels, candidates$variables$catcont, cpuInfo$num_used_cpus)

  #check for nonzero variables
  violators = lapply(1:5, function(x){
    indices = which(norms[[x]] > lambda)
    if (length(indices) > 0) matrix(candidates$variables[[x]][indices, ], nrow=length(indices))
    else NULL})

  #check if any variables should be added to active set
  flag = 1
  for (i in 1:5){
    if (is.null(violators[[i]])) next
    if (!is.null(activeSet[[i]])){
      actives = apply(activeSet[[i]], 1, function(x) paste(x, collapse=":"))
      extras = which(!(apply(violators[[i]], 1, function(x) paste(x, collapse=":")) %in% actives))
    }
    else extras = 1:nrow(violators[[i]])
    if (length(extras) > 0){
      activeSet[[i]] = rbind(activeSet[[i]], violators[[i]][extras, ])
      flag = 0
    }
  }

  list(norms=norms, activeSet=activeSet, flag=flag)
}
    
    
