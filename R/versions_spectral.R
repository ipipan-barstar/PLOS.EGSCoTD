if (1 == 0) { ##########################################
  source("versions_spectral.R");
} # 1===0 ##########################################

combNSTART = 1000
#---------------- FUNCTIONS --------------
latexMatrix <- function(theMatrix, theComment, theCaption, theLabel, theWidth = 160, cornertext = "TRUE/PRED") {
  theMatrix0 = theMatrix;
  fm = 0;

  cat("\n\n\\begin{table} %", theComment);
  while (fm < dim(theMatrix0)[2]) {
    theMatrix0 = theMatrix0[, (fm + 1):dim(theMatrix0)[2]]
    theMatrix = theMatrix0;
    cn = colnames(theMatrix);
    rn = rownames(theMatrix);
    mxrn = max(nchar(c(rn, cornertext)));

    cucn = cumsum(nchar(cn));

    #print(mxrn);
    #print(cucn);
    #print(cucn+mxrn>theWidth);
    fm = which(cucn + mxrn > theWidth)[1];
    if (!is.na(fm)) {
      fm = fm - 1
    } else {
      fm = length(cn);
    }

    if (dim(theMatrix)[2] > fm)
      theMatrix = theMatrix[, 1:fm];


    cat("\n\\begin{tabular}");
    cat(paste0("{|r", paste(rep("|r", dim(theMatrix)[2]), collapse = ""), "|}"));
    cat("\n\\hline", cornertext);
    for (k in 1:dim(theMatrix)[2])
      cat("&", EN.make.tex.names(colnames(theMatrix)[k]));
    cat("\\\\\n");
    for (j in 1:dim(theMatrix)[1]) {
      cat("\\hline");
      cat(" ", EN.make.tex.names(rownames(theMatrix)[j]));
      for (k in 1:dim(theMatrix)[2]) {
        cat("&", theMatrix[j, k]);
      }
      cat("\\\\\n");
    }
    cat("\\hline");
    cat("\n\\end{tabular}");
  } # while
  cat(paste0("\n\\caption{", theCaption, "}\n\\label{tab:", theLabel, "}\n"));
  cat("\n\\end{table} %", theComment, "\n\n");
}


MYMAIN <- function(NoM, txt) {
  if (NoM) {
    return("")
  } else {
    return(txt);
  }
}

nidx <- function(vi, vlambda) {
  return((vi - 1) / (length(vlambda) - 1));
}

nidx_inv <- function(nindex, vlambda) {
  # nindex=  ((vi-1)/(length(vlambda)-1));
  # nindex * (length(vlambda)-1) = vi-1
  # nindex * (length(vlambda)-1) +1 = vi
  vi = nindex * (length(vlambda) - 1) + 1
  return(vi);
}

rlamb <- function(vi, vlambda) {
  return(vlambda[vi] / vlambda[1]);
}

alamb <- function(vi, vlambda) {
  return(vlambda[vi] / length(vlambda));
}

eigenValueAtNormIndex <- function(ev, nindex) {
  #nindex ranging from 0 to 1; ev ordered decreasingly
  li = floor((length(ev) - 1) * nindex) + 1;
  ui = ceiling((length(ev) - 1) * nindex) + 1;
  if (li < 1)
    li = 1;
  if (ui < 1)
    ui = 1;
  if (ui == 1)
    return(ev[ui]);
  if (ui == li)
    return(ev[ui]);
  return(ev[li] + (ev[ui] - ev[li]) * (nidx_inv(nindex, ev) - li));
}

nzso <- function(S, ZERO = 0) {# non zero similarity objects
  nzs = c();
  for (j in 1:dim(S)[1]) {
    rs = sum(S[j,]);
    if (rs - S[j, j] > ZERO)
      nzs = c(nzs, j);
  }
  return(nzs);
}


cmpClusterings <- function(clsTrue, clsComputed, theComment = "", theCaption = "", theLabel = "") {
  noCls = length(table(clsTrue))
  combcl = paste(clsTrue, ':', clsComputed);
  combclTab = sort(table(combcl), decreasing = TRUE);
  if (length(combclTab) == noCls)
    en = 0
  if (length(combclTab) == noCls) {
    cat("\nperfect match of clusterings\n")
  } else {
    cat("\n error matrix\n");
    em = table(clsTrue, clsComputed);
    print(em);
    latexMatrix(em, theComment, theCaption, theLabel)

    cat("\ntoo many cluster cells above ", length(combclTab), "> required", noCls, "\n");
    summx = 0;
    for (j in 1:dim(em)[2])
      summx = summx + max(em[, j]);
    cat(
      "number of elements in correct clusters", summx
      , "\nincorrectly clustered: ", sum(em) - summx
      , "\n= errors: ", (sum(em) - summx) / sum(em) * 100, "%"
      , "\n");
    en = round((sum(em) - summx) / sum(em) * 100, digit = 2)
  }
  flush.console();
  return(en);
}

similarityToDiagonal <- function(S) {
  no = dim(S)[1];
  d = rep(0, no);
  for (j in 1:no)
    d[j] = sum(S[j,]);
  return(diag(d));
}

similarityToExdiagDiagonal <- function(S) {
  no = dim(S)[1];
  d = rep(0, no);
  for (j in 1:no) {
    d[j] = sum(S[j,]) - S[j, j];
    if (abs(d[j]) < 1e-20)
      d[j] = 1;
  }
  return(diag(d));
}

similaritySVDapproximation <- function(S, r) {
  for (j in 1:dim(S)[1])
    S[j, j] = 0
  A = svd(S);
  #Scopy=A$u%*% diag(A$d)%*%t(A$v)
  #print(A$d[length(A$d)+(-5:0)])
  rem = sum(A$d < 1e-14)
  A$d[A$d < 1e-14] = 0
  #Sa=A$u[,1:r]%*% diag(A$d[1:r])%*%t(A$v[,1:r])
  Sa = A$u[,] %*% diag(A$d[]) %*% t(A$v[,])
  # X\approx U_r\Sigma_rV_r^T$
  for (j in 1:dim(Sa)[1])
     Sa[j, j] = 0
  for (j in 1:dim(Sa)[1])
    for (k in 1:dim(Sa)[1])
      if (Sa[j, k] < 0) {
        Sa[j, k] = 0;
        Sa[k, j] = 0;
      }
  d = Sa - S;
  #print(Sa)
  dd = sum(abs(d))
  for (j in 1:dim(Sa)[1]) {
    s = sum(Sa[j,])
    if (s < 1e-10)
      print(c(j, s))
  }
  cat("\nSVD difference", dd, "after removing", rem, "cols\n");

  return(Sa);
}

similarityToCombinatorialLaplacian <- function(S) {
  D = similarityToDiagonal(S);
  L = D - S;
  return(L);
}

similarityToNormalizedLaplacian <- function(S) {
  D = similarityToDiagonal(S);
  Dexd = similarityToExdiagDiagonal(S);
  Dmhalf = diag(1 / sqrt(diag(Dexd))) # D^{-1/2}
  L = Dmhalf %*% (D - S) %*% Dmhalf;
  return(L);
}

IDM <- function(n) {
  v = rep(1, n);
  idm = diag(v);
  return(idm);
}

similarityToKamvarLaplacian <- function(S) {
  D = similarityToDiagonal(S);
  dmax = max(D);
  L = (S + dmax * IDM(dim(S)[1]) - D) / dmax;
  return(L);
}

dataToSimilarityGauss <- function(Data) {
  no = dim(Data)[1];
  na = dim(Data)[2];
  sigma2 = 0;
  for (ja in 1:na)
    sigma2 = sigma2 + var(Data[, ja]);
  sigma2 = sigma2 / 16;

  cat("\nsigma2", sigma2);

  S = matrix(nrow = no, ncol = no);
  S[,] = 0;
  for (jr in 1:no)
    for (jc in 1:no)
      if (jr != jc) {
        dst2 = sum(1 / 2 * (Data[jr,] - Data[jc,])^2);
        S[jr, jc] = exp(-dst2 / sigma2);
      }
  return(S);
}

CombinatorialSpectral <- function(SimilarityMatrix, trueNoCl, enforced = NULL, unitrow = FALSE, useEv = NULL) {
  # trueNoCl - number of lowest eigenvectors taken into account
  # enforced - the number of clusters to be outputted
  # unitrow - normalize row vectors to unit length
  # useEv - multiply columns with lambda^(1/2*useEv) prior to normalization
  S = SimilarityMatrix;
  L = similarityToCombinatorialLaplacian(S);
  return(UniversalSpectral(L, trueNoCl, enforced, unitrow = unitrow, useEv = useEv));
}

NormalizedSpectral <- function(SimilarityMatrix, trueNoCl, enforced = NULL, highest = FALSE, unitrow = FALSE, useEv = NULL) {
  # trueNoCl - number of lowest eigenvectors taken into account
  # enforced - the number of clusters to be outputted
  # unitrow - normalize row vectors to unit length
  # useEv - multiply columns with lambda^(1/2*useEv) prior to normalization; no improvement notices
  S = SimilarityMatrix;
  L = similarityToNormalizedLaplacian(S);
  return(UniversalSpectral(L, trueNoCl, enforced, highest = highest, unitrow = unitrow, useEv = useEv));
}

KamvarSpectral <- function(SimilarityMatrix, trueNoCl, enforced = NULL) {
  # trueNoCl - number of lowest eigenvectors taken into account
  # enforced - the number of clusters to be outputted
  S = SimilarityMatrix;
  L = similarityToKamvarLaplacian(S);
  return(UniversalSpectral(L, trueNoCl, enforced, highest = TRUE, unitrow = TRUE));
}

# TEST
# cmpClusterings(ClsTXT,KamvarSpectral(S,trueNoCl)$newcls); 
# cmpClusterings(ClsTXT,KamvarSpectral(S,trueNoCl+1,trueNoCl)$newcls); 
# cmpClusterings(ClsTXT,NormalizedSpectral(S,trueNoCl)$newcls); 
# cmpClusterings(ClsTXT,NormalizedSpectral(S,trueNoCl,unitraw=TRUE)$newcls); 
# cmpClusterings(ClsTXT,NormalizedSpectral(S,trueNoCl+1,trueNoCl,unitrow=TRUE)$newcls); 
# cmpClusterings(ClsTXT,NormalizedSpectral(S,trueNoCl+1,trueNoCl,unitrow=TRUE,useEv=1)$newcls); 
# cmpClusterings(ClsTXT,CombinatorialSpectral(S,trueNoCl+1,trueNoCl,unitrow=TRUE)$newcls); 
# cmpClusterings(ClsTXT,CombinatorialSpectral(S,trueNoCl,unitrow=TRUE)$newcls); 
# cmpClusterings(ClsTXT,CombinatorialSpectral(S,trueNoCl,unitrow=TRUE,useEv=-1)$newcls); 
# cmpClusterings(ClsTXT,CombinatorialSpectral(S,trueNoCl)$newcls); 
# cmpClusterings(ClsTXT,CombinatorialSpectral(S,trueNoCl+1)$newcls); 
# cmpClusterings(ClsTXT,CombinatorialSpectral(S,trueNoCl+1)$newcls); 
# cmpClusterings(ClsTXT,NormalizedSpectral(S,trueNoCl+1,trueNoCl,unitrow=TRUE,highest=TRUE)$newcls);
#r=250; cmpClusterings(ClsTXT,NormalizedSpectral(similaritySVDapproximation(S,r),trueNoCl )$newcls);
#r=295; cmpClusterings(ClsTXT,NormalizedSpectral(similaritySVDapproximation(S,r),trueNoCl+1,trueNoCl,unitrow=TRUE )$newcls);
#r=295; cmpClusterings(ClsTXT,NormalizedSpectral(S,trueNoCl+1,trueNoCl,unitrow=TRUE)$newcls); 
#r=250; Saa=similaritySVDapproximation(S,r)
#  cmpClusterings(ClsTXT,NormalizedSpectral(Saa,trueNoCl+1,trueNoCl,unitrow=TRUE )$newcls);

UniversalSpectral <- function(L, trueNoCl, enforced, highest = FALSE, unitrow = FALSE, useEv = NULL, SETZERO = 1e-14) {
  # trueNoCl - number of lowest eigenvectors taken into account
  # enforced - the number of clusters to be outputted
  # unitrow - normalize row vectors to unit length
  # useEv - multiply columns with lambda^(1/2*useEv) prior to normalization; no improvement notices
  # highest - if true, the eigenvectors of highest eigenvalues are used, otherwise those of the lowest

  sss = eigen(L);

  for (j in 1:dim(L)[1])
    if (abs(sss$values[j]) < SETZERO)
      sss$values[j] = 0
  for (j in 1:dim(L)[1])
    for (k in 1:dim(L)[1])
      if (abs(sss$vectors[j, k]) < SETZERO)
        sss$vectors[j, k] = 0
  #print(sss$values);
  idFiedler = length(sss$values) - 1;
  #print(idFiedler);

  #while(sss$values[idFiedler]<1e-10) idFiedler=idFiedler-1;
  #print(sss$values[idFiedler]);
  newcls = as.numeric(sss$vectors[, idFiedler] > 0) + 1
  if (trueNoCl >= length(sss$values))
    trueNoCl = length(sss$values) - 1;
  if (is.null(enforced))
    enforced = trueNoCl;
  if (enforced >= length(sss$values))
    enforced = floor(length(sss$values) / 2);
  #print(highest);
  X = sss$vectors[, length(sss$values) - (trueNoCl - 1):0]

  if (trueNoCl > 2 || enforced > 2) {
    if (trueNoCl == 2)
      trueNoCl = 3;
    if (highest) {
      X = sss$vectors[, 1:trueNoCl]
    } else {
      X = sss$vectors[, length(sss$values) - (trueNoCl - 1):0]
    #print(head(X[1:6+326,]));
    }
    if (!is.null(useEv)) {
      if (highest) {
        lam = sss$values[1:trueNoCl]
      } else {
        lam = sss$values[length(sss$values) - (trueNoCl - 1):0]
      }
      lam[lam < 1e-6] = 0;
      for (j in 1:length(lam))
        if (lam[j] > 0)
          X[, j] = X[, j] * lam[j]^(0.5 * useEv)
    }
    if (unitrow)
      for (j in 1:dim(X)[1]) {
        lrow = sqrt(sum(X[j,] * X[j,]));
        if (lrow > 1e-6)
          X[j,] = X[j,] / lrow;
      }
    #print(X[1:6+326,]);
    #print(X[,1]);
    fit = kmeans(X, center = enforced, nstart = combNSTART);
    newcls = fit$cluster;
  }
  res = list();
  res$newcls = newcls;
  res$X = X;
  res$embed = sss;
  return(res);
}

CombinatorialSpectralOLD <- function(SimilarityMatrix, trueNoCl, enforced = NULL) {
  # trueNoCl - number of lowest eigenvectors taken into account
  # enforced - the number of clusters to be outputted
  S = SimilarityMatrix;
  L = similarityToCombinatorialLaplacian(S);

  sss = eigen(L);
  #print(sss$values);
  idFiedler = length(sss$values) - 1;
  #print(idFiedler);

  #while(sss$values[idFiedler]<1e-10) idFiedler=idFiedler-1;
  #print(sss$values[idFiedler]);
  newcls = as.numeric(sss$vectors[, idFiedler] > 0) + 1
  if (trueNoCl >= length(sss$values))
    trueNoCl = length(sss$values) - 1;
  if (is.null(enforced))
    enforced = trueNoCl;
  if (enforced >= length(sss$values))
    enforced = floor(length(sss$values) / 2);

  if (trueNoCl > 2 || enforced > 2) {
    if (trueNoCl == 2)
      trueNoCl = 3;
    fit = kmeans(sss$vectors[, length(sss$values) - (trueNoCl - 1):0], center = enforced, nstart = combNSTART);
    newcls = fit$cluster;
  }
  res = list();
  res$newcls = newcls;
  res$embed = sss;
  return(res);
}

NormalizedSpectralOLD <- function(SimilarityMatrix, trueNoCl, enforced = NULL) {
  S = SimilarityMatrix;
  L = similarityToNormalizedLaplacian(S);

  sss = eigen(L);
  idFiedler = length(sss$values) - 1;
  #while(sss$values[idFiedler]<1e-10) idFiedler=idFiedler-1;
  #print(sss$values[idFiedler]);
  newcls = as.numeric(sss$vectors[, idFiedler] > 0) + 1
  if (trueNoCl >= length(sss$values))
    trueNoCl = length(sss$values) - 1;
  if (is.null(enforced))
    enforced = trueNoCl;
  if (enforced >= length(sss$values))
    enforced = floor(length(sss$values)) / 2;
  if (trueNoCl > 2 || enforced > 2) {
    if (trueNoCl == 2)
        trueNoCl = 3;
    fit = kmeans(sss$vectors[, length(sss$values) - (trueNoCl - 1):0], center = enforced, nstart = combNSTART);
    newcls = fit$cluster;
  }
  res = list();
  res$newcls = newcls;
  res$embed = sss;
  return(res);
}

EN.make.names <- function(nml, onlyPLK = FALSE) {
  txt = nml
  if (!onlyPLK)
    txt = make.names(txt);
  return(txt);
}

EN.make.tex.names <- function(nml, onlyPLK = TRUE) {
  txt = nml
  if (!onlyPLK)
    txt = make.names(txt);
  txt = gsub('#', "\\#", txt, fixed = TRUE);
  txt = gsub('_', "\\_", txt, fixed = TRUE);
  return(txt);
}

# S2=similaritySVDapproximation(S,r=250);
# e2=eigen(S2); 


############################################
############################################
##########                        ##########
########## AD MAIOREM DEI GLORIAM ##########
##########                        ##########
############################################
############################################


