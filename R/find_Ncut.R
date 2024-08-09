if (1 == 0) { ###################################################################
  source("find_Ncut.R");
} #############################################################################

while (sink.number() > 0)
  sink();
while (!is.null(dev.list()))
  dev.off();

source("cut_criteria.R");
source("kmeansWITHweights.R");

cat("
################################################################
# find_Ncut.R -  Derivation
# M-based Clustering
################################################################
");


mkSYMprint <- function(M, Mname) {
  cat("\nSymmetrize  Matrix ", Mname, "\n");
  print(dim(M))
  for (j in 1:dim(M)[1])
    for (k in 1:dim(M)[2]) {
      s = (M[j, k] + M[k, j]) / 2;
      M[j, k] = s;
      M[k, j] = s;
    }
  diff = 0;
  for (j in 1:dim(M)[1])
    for (k in 1:dim(M)[2])
      if (abs(M[j, k] - M[k, j]) > diff)
        diff = abs(M[j, k] - M[k, j])
  cat("\n assymetry=", diff, "\n");
  flush.console()
  return(M)
}

source("versions_spectral.R");

cat("\n Saving the data: RES0 and nnn_cats");

targetDir = paste0("../Results/", substr(nnn.type, 1, 3));
target.data = cbind(theGROUP = nnn_cats, RES0)

if (!dir.exists(targetDir)) {
  dir.create(targetDir);
  cat("\nDirectory ", targetDir, "created\n\n");
}

write.csv(target.data, paste0(targetDir, "/", nnn.type, "_DATA_and_cats.csv"));
flush.console();

n = dim(S)[1]

# S diagonal must be 0
for (j in 1:dim(S)[1])
  S[j, j] = 0

D = matrix(nrow = n, ncol = n)
D[,] = 0
for (j in 1:n)
  D[j, j] = sum(S[j,])

# correction for unrelated elements 
toRemove = which(diag(D) < 1e-10)
if (length(toRemove) > 0) {
  toRetain = which(diag(D) >= 1e-10)
  nnn_cats = nnn_cats[toRetain]
  S = S[toRetain, toRetain]
  D = D[toRetain, toRetain]
  RES0 = RES0[toRetain,]
  cat("\n the following documents were removed:");
  print(toRemove);
  n = dim(S)[1]
}

w = diag(D)

Dm1 = diag(diag(D)^-1)
nL = Dm1^0.5 %*% (D - S) %*% Dm1^0.5

IV = 1:n
IV[] = 1
I = diag(IV)

Isq = matrix(nrow = n, ncol = n)
Isq[,] = 1
E = Isq - I
DSQ = Dm1 %*% (D %*% E + E %*% D - 2 * S) %*% Dm1

M = -0.5 * (I - Isq / n) %*% DSQ %*% (I - Isq / n)

M = mkSYMprint(M, "M")


egM = eigen(M)
#print(egM)
LambdaV = egM$values
if (length(LambdaV[LambdaV < 0]) < 2)
  LambdaV[LambdaV < 0] = 0

if (IsPNG)
  png(paste0(targetDir, "/", "EigenvaluesMbased", nnn.type, ".png"));
plot(LambdaV, main = "Eigenvalues M-based", ylab = "eigenvalue")
if (IsPNG)
  dev.off()

mEv = min(LambdaV)
cat("\nMinimal eigenvalue of M", mEv, "\n")

if (mEv < 0) {
  sigma = -2 * mEv + 1e-10

  for (j in 1:dim(DSQ)[1])
    for (k in 1:dim(DSQ)[2])
      if (j != k)
        DSQ[j, k] = DSQ[j, k] + sigma
  M = -0.5 * (I - Isq / n) %*% DSQ %*% (I - Isq / n)
  M = mkSYMprint(M, "M")

  egM = eigen(M)
  frLambdaV = LambdaV
  LambdaV = egM$values
  print(frLambdaV[frLambdaV < 0])
  print(LambdaV[LambdaV < 0])
  mEv = min(LambdaV)
  cat("\nretrial Minimal eigenvalue of M", mEv, "\n")
  if (length(LambdaV[LambdaV < 0]) < 2)
    LambdaV[LambdaV < 0] = 0
}
#LambdaV[LambdaV<0]=0
Lambda = diag(LambdaV)
V = egM$vectors
VT = t(V)

i = 2
l = which(S[i,] > 0.001)[1]
x_i = sqrt(Lambda) %*% VT[, i]
x_l = sqrt(Lambda) %*% VT[, l]
dstil = sum((x_i - x_l)^2)
cat("\nS[", i, ",", l, "]=", S[i, l]);
cat("\nsquared distance in the embedding: ", dstil)
cat("\ntrue squared distance              ", DSQ[i, l])
cat("\nintended squared distance          ", (D[i, i] + D[l, l] - 2 * S[i, l]) / (D[i, i] * D[l, l]))

#supertest 
if (1 == 0)
  for (i in 1:n)
    for (l in 1:n) {
      x_i = sqrt(Lambda) %*% VT[, i]
      x_l = sqrt(Lambda) %*% VT[, l]
      dstil = sum((x_i - x_l)^2)
      DSQil = DSQ[i, l]
      trueV = (D[i, i] + D[l, l] - 2 * S[i, l]) / (D[i, i] * D[l, l])
      if (abs(dstil - DSQil) > 1e-10)
        stop(paste("dstil DSQil", i, l));
      if (i != l)
        if (abs(dstil - trueV) > 1e-10)
          stop(paste("dstil trueV", i, l));
    }  #

x_i = sqrt(Lambda) %*% VT[, i]
XT = sqrt(Lambda) %*% VT
sum(abs(x_i - XT[, i]))
X = t(XT)

trueNoCl = length(table(nnn_cats))
cat("\nTrue clustering distribution\n");
print(table(nnn_cats))


Nbased = NormalizedSpectral(S, trueNoCl, unitrow = FALSE)
#rememer Lbased=CombinatorialSpectral(S,trueNoCl+5,unitrow=FALSE)
clNbased = Nbased$newcls
cat("\nN-based clustering distribution\n");
print(table(clNbased))


L_LambdaV = Nbased$embed$values
if (IsPNG)
  png(paste0(targetDir, "/", "EigenvaluesNbased", nnn.type, ".png"));
plot(L_LambdaV, main = "Eigenvalues N-based (traditional)", col = "green3", ylab = "eigenvalue")
if (IsPNG)
  dev.off()
plot(c(L_LambdaV, LambdaV[1]), main = "Eigenvalues M (black) and N (green) -based (traditional)", col = "green3")
points(LambdaV)


dic = names(table(nnn_cats))
n_cats = rep(NA, length(nnn_cats))
for (j in 1:length(n_cats))
  n_cats[j] = which(nnn_cats[j] == dic)[1]

k = max(n_cats)
cat("\n n=", n, ", k=", k, " n-2k=", n - 2 * k)
cat("\n NCut criterion", NCut(S, n_cats))
cat("\n Qncut criterion", (n - 2 * k + NCut(S, n_cats)))
xxx = kmeansWITHweightsQ(X[,], n_cats, theWeights = diag(D))
cat("\n Qncut criterion from sum to cluster centers", xxx$tot.withinss)
flush.console();
yyy = kmeansWITHweightsQvD(X[,], n_cats, theWeights = diag(D))
cat("\n Qncut criterion from sum to other cluster  elements", yyy$tot.withinss)

readline("\n-----> see the criteria resuls above and the eigenvalue plots: >>>");

# HERE CASE WEIGHTED K-MEANS NEEDED!!!!!!
fitX = kmeansWITHweights(X[, 1:(trueNoCl + 1)], trueNoCl, nstart = 3 * 9, theWeights = diag(D))

clMbased = fitX$cluster
clMembed = X[, 1:trueNoCl]

cat("
cluster sizes", paste(fitX$size, collapse = ", "), "
withinss:", paste(fitX$withinss, collapse = ", "), "
tot.withinss: ", fitX$tot.withinss, "
removed: ", fitX$removed, "
");


readline("see weighted clustering results");

plotEmbed(clMembed, clMbased, "M embedding", theWeights = fitX$theNewWeights)


p = dim(Nbased$embed$vectors)[2]
clNembed = Nbased$embed$vectors[, p - 1:trueNoCl]
plotEmbed(clNembed, clNbased, "N embedding")

plotEmbed(clMembed, nnn_cats, "hashtags in M embedding", theWeights = fitX$theNewWeights)
plotEmbed(clNembed, nnn_cats, "hashtags in N embedding")


cat("\nM-based clustering distribution\n");
print(table(clMbased))


cat("\nM-based clustering distribution\n");
print(table(clMbased))
cat("\nM-based clustering distribution versus true one\n");
print(table(clMbased, nnn_cats))
print(table(clMbased[fitX$theNewWeights > 0], nnn_cats[fitX$theNewWeights > 0]))
cat("\nM-based clustering distribution versus N-based\n");
print(table(clMbased, clNbased))
cat("\nN-based clustering distribution versus true one\n");
print(table(clNbased, nnn_cats))

clrandom = sample.int(trueNoCl, dim(S)[1], replace = TRUE, prob = as.numeric(table(nnn_cats)))


prRCut("clNbased", S, clNbased)
prRCut("clMbased", S, clMbased)
prRCut("true", S, nnn_cats)
prRCut("random", S, clrandom)

prNCut("clNbased", S, clNbased)
prNCut("clMbased", S, clMbased)
prNCut("true", S, nnn_cats)
prNCut("random", S, clrandom)


cat("
Some results to be saved in "
  , paste0(targetDir, "/", "normclusterStats", nnn.type, ".txt"), "
");

sink(paste0(targetDir, "/", "normclusterStats", nnn.type, ".txt"))
cat("\n TRUE NUMBER OF CLUSTERS: ", trueNoCl, "\n");
cat("M-based centers\n");
print(fitX$centers)
cat("N-based centers\n");
print(clNbased)

cat("\n Does N-based GSC imply  our M-based\n");
cmpClusterings(clMbased, clNbased, theComment = " Does N-based GSC imply  our M-based GSC")
cat("\n Does our M-based method imply N-based GSC  \n");
cmpClusterings(clNbased, clMbased)

cat("\n Is true clustering implied with our M-based method \n");
cmpClusterings(nnn_cats, clMbased, theComment = "Is true clustering implied with our M-based method")

cat("\n Is true clustering  implied  with N-based GSC\n");
cmpClusterings(nnn_cats, clNbased, theComment = "Is true clustering implied with N-based GSC")

topwords <- function(RES0, clustering, theWeights) {
  for (j in names(table(clustering)))
    if (sum(clustering == j) > 1) {
      RES0j = RES0[clustering == j & theWeights > 0,]
      wc = c()
      for (k in 1:dim(RES0j)[2]) {
        w = sum(RES0j[, k])
        wc = c(wc, w)
      }
      ix = sort(wc, index.return = TRUE, decreasing = TRUE)
      print(colnames(RES0)[ix$ix[1:10]])
    }
}


topwordsdiffs <- function(RES0, clustering, theWeights) {
  for (j in names(table(clustering)))
    if (sum(clustering == j) > 1) {
      RES0j = RES0[clustering == j & theWeights > 0,]
      wc = c()
      for (k in 1:dim(RES0j)[2]) {
        w = sum(RES0j[, k])
        wc = c(wc, w)
      }
      wcXXX = wc
      wcXXX[] = 0
      wcXXXsigned = wcXXX
      for (j2 in names(table(clustering)))
        if (sum(clustering == j2) > 1 && j != j2) {
          RES0j2 = RES0[clustering == j2 & theWeights > 0,]
          wc2 = c()
          for (k in 1:dim(RES0j)[2]) {
            w = sum(RES0j2[, k])
            wc2 = c(wc2, w)
          }
          wcXXX = wcXXX + (wc - wc2)^2
          wcXXXsigned = wcXXXsigned + (wc - wc2)
        }
      #print(wc)
      ix = sort(wcXXX, index.return = TRUE, decreasing = TRUE)
      #print(ix)
      cat("\n Cluster: ", j, "\n");
      cnn = colnames(RES0)[ix$ix[1:10]]
      signs = wcXXXsigned[ix$ix[1:10]] > 0
      cnn[!signs] = paste("(-)", cnn[!signs])
      cnn[signs] = paste("(++)", cnn[signs])
      print(cnn);
   }
}

cat("\n True clusters top words \n");
topwords(RES0, nnn_cats, 1)


cat("\n N-based clusterings top words \n");
topwords(RES0, clNbased, 1)

cat("\n M-based clusterings top words \n");
topwords(RES0, clMbased, fitX$theNewWeights)
cat("\n M-based clusterings top words differentiated \n");
topwordsdiffs(RES0, clMbased, fitX$theNewWeights)


sink();


##########################
# AD MAIOREM DEI GLORIAM #
##########################
cat("
###########
# THE END #
###########

call now
   source('TWT_read_exp.R');
or
   source('BLK_readX.R');
after changing appropriately the parameters in these files 
");
while (sink.number() > 0)
  sink();


