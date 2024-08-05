if (1 == 0) { ###################################################################
  source("find_Rcut.R");
} #############################################################################

while (sink.number() > 0) sink();
while (!is.null(dev.list())) dev.off();

source("cut_criteria.R");
cat("
################################################################
# find_Rcut.R -  Derivation
# K-based Clustering
################################################################
");
prRCut("true", S, nnn_cats)
flush.console();

mkSYMprint <- function(M, Mname)
{
  cat("\nSymmetrize  Matrix ", Mname, "\n"); print(dim(M))
  for (j in 1:dim(M)[1])
    for (k in 1:dim(M)[2])
    { s = (M[j, k] + M[k, j]) / 2; M[j, k] = s; M[k, j] = s; }
  flush.console()
  return(M)
}

# eof


source("versions_spectral.R");

cat("\n Saving the data: RES0 and nnn_cats");

targetDir = paste0("Wyniki_", substr(nnn.type, 1, 3));
target.data = cbind(theGROUP = nnn_cats, RES0)

if (!dir.exists(targetDir))
{
  dir.create(targetDir);
  cat("\nDirectory ", targetDir, "created\n\n");
}

write.csv(target.data, paste0(targetDir, "/", nnn.type, "_DATA_and_cats.csv"));

flush.console();


n = dim(S)[1]

# S diagonal must be 0
for (j in 1:dim(S)[1]) S[j, j] = 0

D = matrix(nrow = n, ncol = n)
D[,] = 0
for (j in 1:n) D[j, j] = sum(S[j,])

L = D - S


IV = 1:n
IV[] = 1
I = diag(IV)

Isq = matrix(nrow = n, ncol = n)
Isq[,] = 1
DSQ = (Isq - I - S)
K = -0.5 * (I - Isq / n) %*% DSQ %*% (I - Isq / n)
K = mkSYMprint(K, "K")


IsPNG = FALSE


egK = eigen(K)
#print(egK)
LambdaV = egK$values
if (length(LambdaV[LambdaV < 0]) < 2)
  LambdaV[LambdaV < 0] = 0

if (IsPNG) png(paste0(targetDir, "/", "EigenvaluesKbased", nnn.type, ".png"));
plot(LambdaV, main = "Eigenvalues K-based", ylab = "eigenvalue")
if (IsPNG) dev.off()

mEv = min(LambdaV)
cat("\nMinimal eigenvalue of K", mEv, "\n")
#readline("\nlook"); 

if (mEv < 0)
{ sigma = -2 * mEv + 1e-10

  for (j in 1:dim(DSQ)[1])
    for (k in 1:dim(DSQ)[2])
      if (j != k)
        DSQ[j, k] = DSQ[j, k] + sigma
  K = -0.5 * (I - Isq / n) %*% DSQ %*% (I - Isq / n)
  K = mkSYMprint(K, "K")

  egK = eigen(K)
  #print(egK)
  frLambdaV = LambdaV
  LambdaV = egK$values
  print(frLambdaV[frLambdaV < 0])
  print(LambdaV[LambdaV < 0])
  mEv = min(LambdaV)
  cat("\nretrial Minimal eigenvalue of K", mEv, "\n")
  if (length(LambdaV[LambdaV < 0]) < 2)
    LambdaV[LambdaV < 0] = 0
}

#LambdaV[LambdaV<0]=0
Lambda = diag(LambdaV)
V = egK$vectors
VT = t(V)

i = 2
l = which(S[i,] > 0.001)[1]
x_i = sqrt(Lambda) %*% VT[, i]
x_l = sqrt(Lambda) %*% VT[, l]
dstil = sum((x_i - x_l)^2)
cat("\nsquared distance in the embedding: ", dstil)
cat("\ntrue squared distance              ", DSQ[i, l])
cat("\nintended squared distance          ", 1 - S[i, l])

trueNoCl = length(table(nnn_cats))
cat("\nTrue clustering distribution\n");
print(table(nnn_cats))

Lbased = CombinatorialSpectral(S, trueNoCl, unitrow = FALSE)
#Lbased=CombinatorialSpectral(S,trueNoCl+5,unitrow=FALSE)
clLbased = Lbased$newcls
cat("\nL-based clustering distribution\n");
print(table(clLbased))
print(table(clLbased, nnn_cats))


L_LambdaV = Lbased$embed$values
if (IsPNG) png(paste0(targetDir, "/", "EigenvaluesLbased", nnn.type, ".png"));
plot(L_LambdaV, main = "Eigenvalues L-based (traditional)", col = "green3", ylab = "eigenvalue")
if (IsPNG) dev.off()
plot(c(L_LambdaV, LambdaV[1]), main = "Eigenvalues K (black) and L(green) -based (traditional)", col = "green3")
points(LambdaV)
readline("look and see the results of L-based clustering above >>>");

fitX = kmeans(V[, 1:(trueNoCl + 1)], trueNoCl, nstart = 90)
print(fitX$tot.withinss)
clKbased = fitX$cluster
cat("\nK-based clustering distribution\n");
print(table(clKbased))
print(table(clKbased, nnn_cats))
print(table(clKbased, clLbased))


clKembed = V[, 1:trueNoCl]

experiment$figex = "Kemb"
plotEmbed(clKembed, clKbased, "K embedding")


p = dim(Lbased$embed$vectors)[2]
clLembed = Lbased$embed$vectors[, p - 1:trueNoCl]
experiment$figex = "Lemb"
plotEmbed(clLembed, clLbased, "L embedding")


experiment$figex = "KembTrue"
plotEmbed(clKembed, nnn_cats, "hashtags in K embedding")
experiment$figex = "LembTrue"
plotEmbed(clLembed, nnn_cats, "hashtags in L embedding")


clrandom = sample.int(trueNoCl, dim(S)[1], replace = TRUE, prob = as.numeric(table(nnn_cats)))


prRCut("clLbased", S, clLbased)
prRCut("clKbased", S, clKbased)
prRCut("true", S, nnn_cats)
prRCut("random", S, clrandom)

prNCut("clLbased", S, clLbased)
prNCut("clKbased", S, clKbased)
prNCut("true", S, nnn_cats)
prNCut("random", S, clrandom)


sink(paste0(targetDir, "/", "clusterStats", nnn.type, ".txt"))

cat("\n Does L-based GSC imply  our K-based\n");

cmpClusterings(clKbased, clLbased, theComment = " Does L-based GSC imply  our K-based GSC")
cat("\n Does our K-based method imply L-based GSC  \n");
cmpClusterings(clLbased, clKbased)

cat("\n Is true clustering implied with our K-based method \n");
cmpClusterings(nnn_cats, clKbased, theComment = "Is true clustering implied with our K-based method", theLabel = "trueK")

cat("\n Is true clustering  implied  with L-based GSC\n");
cmpClusterings(nnn_cats, clLbased, theComment = "Is true clustering implied with L-based GSC", theLabel = "trueL")

topwords <- function(RES0, clustering)
{
  for (j in names(table(clustering)))
    if (sum(clustering == j) > 1)
    {
      RES0j = RES0[clustering == j,]
      wc = c()
      for (k in 1:dim(RES0j)[2])
      {
        w = sum(RES0j[, k])
        wc = c(wc, w)
      }
      #print(wc)
      ix = sort(wc, index.return = TRUE, decreasing = TRUE)
      #print(ix)
      cat("\n Cluster: ", j, "\n");
      print(colnames(RES0)[ix$ix[1:10]]);
    }
}

# eof


topwordsdiffs <- function(RES0, clustering)
{
  for (j in names(table(clustering)))
    if (sum(clustering == j) > 1)
    {
      RES0j = RES0[clustering == j,]
      wc = c()
      for (k in 1:dim(RES0j)[2])
      {
        w = sum(RES0j[, k])
        wc = c(wc, w)
      }
      wcXXX = wc
      wcXXX[] = 0
      wcXXXsigned = wcXXX
      for (j2 in names(table(clustering)))
        if (sum(clustering == j2) > 1 && j != j2)
        {
          RES0j2 = RES0[clustering == j2,]
          wc2 = c()
          for (k in 1:dim(RES0j)[2])
          {
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

# eof


cat("\n True clusters top words \n");
topwords(RES0, nnn_cats)
cat("\n True clusters top words differentiated \n");
topwordsdiffs(RES0, nnn_cats)

cat("\n True clusters top words \n");
topwords(RES0, nnn_cats)
cat("\n True clusters top words differentiated \n");
topwordsdiffs(RES0, nnn_cats)


cat("\n L-based clusterings top words \n");
topwords(RES0, clLbased)
cat("\n L-based clusterings top words differentiated \n");
topwordsdiffs(RES0, clLbased)

cat("\n K-based clusterings top words \n");
topwords(RES0, clKbased)
cat("\n K-based clusterings top words  differentiated  \n");
topwordsdiffs(RES0, clKbased)

sink();

##########################
# AD MAIOREM DEI GLORIAM #
##########################
cat("
###########
# THE END #
###########

call now
   source('find_Ncut.R');
");
while (sink.number() > 0) sink();


