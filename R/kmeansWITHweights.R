if (1 == 0) { ###################################################################
  source("kmeansWITHweights.R");
} #############################################################################

cat("
################################################################
# kmeansWITHweight.R - k-means algorithm incorporating case weights
################################################################
");

kmeansWITHweights <- function(x, centers, iter.max = 10L, nstart = 1L,
                              algorithm = c("Hartigan-Wong", "Lloyd", "Forgy", "MacQueen"),
                              trace = FALSE, theWeights = NULL) {
  if (is.null(theWeights))
    theWeights = rep(1, dim(x)[1])
  fit0 <- kmeansWITHweightsOne(x, centers, iter.max = iter.max, algorithm = algorithm, trace = trace, theWeights = theWeights);
  cat(".");
  flush.console();
  if (nstart <= 1L)
    return(fit0);
  for (run in 2L:nstart) {
    fit <- kmeansWITHweightsOne(x, centers, iter.max = iter.max, algorithm = algorithm, trace = trace, theWeights = theWeights);
    cat(".");
    flush.console();
    if (fit$tot.withinss < fit0$tot.withinss)
      fit0 = fit;
  }
  return(fit0)
}

kmeansWITHweightsOne <- function(x, centers, iter.max = 10L,
                                 algorithm = c("Hartigan-Wong", "Lloyd", "Forgy", "MacQueen"),
                                 trace = FALSE, theWeights = NULL, max_stagnationL = 20, clMinSize = 5,
                                 clMinShare = 0.1, maxRemShare = 0.5, navigation = FALSE, test = FALSE) {
  # the current implementation is simplified as follows:
  # x - matrix with rows as cases
  # theWeights - weights of cases
  # centers- the number of centers
  # other parameters, currently not used
  # returns: a list with contents ti be defined
  stagnationG = 0;
  stagnationL = 0;

  MU = matrix(nrow = centers, ncol = dim(x)[2]);
  DST2mu = matrix(nrow = dim(x)[1], ncol = centers);
  clMemb = rep(NA, dim(x)[1])
  cldst2 = rep(NA, dim(x)[1])
  k = centers;
  n = dim(x)[1]
  d = dim(x)[2]
  # random initialization 
  #idx.seeds= sample.int(n,k,prob=theWeights)
  #MU=x[idx.seeds,]
  idx.seed = sample.int(n, 1, prob = theWeights)
  MU[1,] = x[idx.seed,]
  for (k1 in 2:k) {
    DST2mu[,] = 0;
    for (j in 1:(k1 - 1))
      for (jd in 1:d) {
        DST2mu[, j] = DST2mu[, j] + (x[, jd] - MU[j, jd])^2
      }
    # assignment to clusters 
    for (i in 1:n) {
      cldst2[i] = min(DST2mu[i, 1:(k1 - 1)])
      #    clMemb[i]=which(DST2mu[i,1:(k1-1)]==cldst2[i])[1]
    }
    idx.seed = sample.int(n, 1, prob = cldst2)
    MU[k1,] = x[idx.seed,]

  } #
  # iterative computations of clusters 
  # for (iteration in 1:iter.max)
  iteration = 1;
  prevRem = -1;
  while (iteration <= iter.max) {
    iteration = iteration + 1;
    # computing distances to seeds 
    DST2mu[,] = 0;
    for (j in 1:k)
      for (jd in 1:d) {
        DST2mu[, j] = DST2mu[, j] + (x[, jd] - MU[j, jd])^2
      }
    # assignment to clusters 
    for (i in 1:n) {
      cldst2[i] = min(DST2mu[i,])
      clMemb[i] = which(DST2mu[i,] == cldst2[i])[1]
    }
    # new cluster centers
    MU[,] = 0
    MUweights = rep(0, k)
    MUsize = rep(0, k)
    for (i in 1:n) {
      cl = clMemb[i]
      MU[cl,] = MU[cl,] + theWeights[i] * x[i,]
      MUweights[cl] = MUweights[cl] + theWeights[i]
      MUsize[cl] = MUsize[cl] + 1
    }
    for (j in 1:k)
      if ( (MUweights[j] > 1e-10
          && MUsize[j] >= clMinSize
          && MUsize[j] >= clMinShare * n / k)
        || maxRemShare * n < prevRem) {
        MU[j,] = MU[j,] / MUweights[j] # new cluster center computed
      } else {
        if (MUsize[j] > 0) {
          theWeights[j == clMemb] = 0;
          cldst2[j == clMemb] = 0
        }
        idx.seed = sample.int(n, 1, prob = cldst2 * theWeights)
        MU[j,] = x[idx.seed,] # a new cluster center sampled
        iteration = iteration - 1 # additional iteration necessary, if a cluster was removed

        #      dat=x[theWeights>0,]
        #      col=clMemb[theWeights>0];
        rem = sum(theWeights == 0)
        if (prevRem == rem) {
          stagnationG = stagnationG + 1;
          stagnationL = stagnationL + 1;
          if (stagnationL > max_stagnationL + k)
            iteration = iteration + 1
        } else {
          stagnationL = 0
        }
        left = sum(theWeights > 0)
        if (test)
          plotEmbed(x, clMemb, paste(rem, "datapoints removed"), theWeights = theWeights, navigation = navigation)
        if (test)
          if (navigation) readline(paste("\n", rem, "datapoints removed,", left, " are left"))
        prevRem = rem
      }
  } #  for (iteration in 1:iter.max)

  # result preparation
  res = list()
  res$centers = MU
  res$cluster = clMemb
  res$size = MUsize
  res$withinss = rep(0, k)
  res$theNewWeights = theWeights
  res$removed = sum(theWeights == 0)
  # Quality computation:
  Q = 0;
  for (i in 1:n) {
    cl = clMemb[i]
    ws = sum((MU[cl,] - x[i,])^2) * theWeights[i]
    Q = Q + ws
    res$withinss[cl] = res$withinss[cl] + ws
  }
  res$tot.withinss = Q
  return(res)
}

kmeansWITHweightsQ <- function(x, cluster, iter.max = 10L,
                               algorithm = c("Hartigan-Wong", "Lloyd", "Forgy", "MacQueen"),
                               trace = FALSE, theWeights = NULL) {
  # the current implementation is simplified as follows:
  # x - matrix with rows as cases
  # theWeights - weights of cases

  # other parameters, currently not used
  # returns: a list with contents ti be defined

  if (is.null(theWeights))
    theWeights = rep(1, dim(x)[1])
  clMemb = cluster
  MU = matrix(nrow = max(cluster), ncol = dim(x)[2]);
  # new cluster centers
  MU[,] = 0
  k = max(cluster)
  n = dim(x)[1]
  MUweights = rep(0, k)
  for (i in 1:n) {
    cl = clMemb[i]
    MU[cl,] = MU[cl,] + theWeights[i] * x[i,]
    MUweights[cl] = MUweights[cl] + theWeights[i]
  }
  for (j in 1:k) {
     MU[j,] = MU[j,] / MUweights[j]
  }

  # result preparation
  res = list()
  res$centers = MU
  res$cluster = clMemb
  # Quality computation:
  Q = 0;
  for (i in 1:n) {
    cl = clMemb[i]
    Q = Q + sum((MU[cl,] - x[i,])^2) * theWeights[i]
    #    cat("",cl," ",i," " , sum((MU[cl,]-x[i,])^2) * theWeights[i]," Q" ,Q,"\n");
  }
  res$tot.withinss = Q
  return(res)
}


kmeansWITHweightsQvD <- function(x, cluster, iter.max = 10L,
                                 algorithm = c("Hartigan-Wong", "Lloyd", "Forgy", "MacQueen"),
                                 trace = FALSE, theWeights = NULL) {
  # the current implementation is simplified as follows:
  # x - matrix with rows as cases
  # theWeights - weights of cases

  # other parameters, currently not used
  # returns: a list with contents ti be defined

  # Quality computation:
  Q = 0;
  k = max(cluster)
  for (j in 1:k) {
    clj = cluster == j
    Qj = 0;
    Vj = 0;
    for (i in which(clj))
      Vj = Vj + theWeights[i]
    for (i in which(clj))
      for (l in which(clj)) {
        Qj = Qj + theWeights[i] * theWeights[l] * sum((x[i,] - x[l,])^2)
      }
    Q = Q + Qj / Vj / 2
  }
  # result preparation
  res = list()
  res$tot.withinss = Q
  return(res)
}


##########################
# AD MAIOREM DEI GLORIAM #
##########################


############################################
############################################
##########                        ##########
########## AD MAIOREM DEI GLORIAM ##########
##########                        ##########
############################################
############################################

