if (1 == 0) { ###################################################################
  source("S_describeBis.R");
} #############################################################################

while (sink.number() > 0) sink();
while (!is.null(dev.list())) dev.off();

cat("
################################################################
# S_describeBis.R -  describer of  S
################################################################
");


ifelse(!dir.exists(file.path(getwd(), "../Results/figs/")), dir.create(file.path(getwd(), "../Results/figs")), FALSE)

S = S - diag(diag(S))

n = dim(S)[1]
ix = sort(nnn_cats, index.return = TRUE)


So = S[ix$ix, ix$ix]

trv = quantile(S, prob = c((n + 400) / n^2, 1 - 400 / n^2));
tr = trv[2]
tl = trv[1]
#sqrt(dim(S)[1]*dim(S)[1]*(1-0.9999)) 
z = which(So >= tr)
x = (z - 1) %% dim(S)[1] + 1;
y = (z - x) / dim(S)[1];
zl = which(So <= tl)
xl = (zl - 1) %% dim(S)[1] + 1;
yl = (zl - xl) / dim(S)[1];

dg = abs(xl - yl) < 2
xl[dg] = NA
yl[dg] = NA

if (length(xl) > 450)
{
  id = sample.int(length(xl), 400)
  xl = xl[id]
  yl = yl[id]
}


print(sum(is.na(S)));


dst_cats = table(nnn_cats)
dst_cum = cumsum(dst_cats)
print(dst_cats);

experiment$figex = "mxlinks"

lims = c(1, n)

for (isPDF in c(FALSE, TRUE))
{
  if (isPDF) pdf(paste0("../Results/figs/", experiment$name, "_", experiment$figex, ".pdf"))

  plot(x, y, main = paste("nodes connected with weight above", round(tr, digit = 3)), pch = 20, col = "green", xlim = lims, ylim = lims)
  points(xl, yl, col = "grey") # 400 lowest similarities

  for (kk in names(dst_cats))
  {
    points(x[as.character(nnn_cats) == kk], y[as.character(nnn_cats) == kk], col = paste0("blue", min(c(which(kk == names(table(nnn_cats))), 4))), pch = 14 + min(c(which(kk == names(table(nnn_cats))), 4)))

    text(dst_cum[kk] - dst_cats[kk] / 2, dst_cum[kk] - dst_cats[kk] / 2, labels = c(kk))
    lines(
      lims
      , c(dst_cum[kk], dst_cum[kk])
      , col = "green"
    );
    lines(
      c(dst_cum[kk], dst_cum[kk])
      , lims
      , col = "green"
    );

  }
  if (!isPDF) readline("see ?node connection weights>>");
  if (isPDF) dev.off()
}

topy = round(n * 0.05)


DV = rep(0, n) # diagonal sum
DVi = rep(0, n) # diagonal sum inside clusters
DVo = rep(0, n) # diagonal sum outside clusters
aDV = rep(0, n) # diagonal average
aDVi = rep(0, n) # diagonal average inside clusters
aDVo = rep(0, n) # diagonal average outside clusters
DVtop100 = rep(0, n) # diagonal sum top 100 elements
DVlow100 = rep(0, n) # diagonal sum lowest 100 elements

for (j in 1:n)
{ DV[j] = sum(S[j,])
  aDV[j] = DV[j] / n
  kk = nnn_cats[j]
  DVi[j] = sum(S[j, nnn_cats == kk])
  aDVi[j] = DVi[j] / sum(nnn_cats == kk)
  DVo[j] = sum(S[j, nnn_cats != kk])
  aDVo[j] = DVo[j] / sum(nnn_cats != kk)
  x = sort(S[j,])
  DVtop100[j] = sum(x[length(x) - 1:topy + 1])
  DVlow100[j] = sum(x[1:topy])

}
aDVlow100 = DVlow100 / 100
aDVtop100 = DVtop100 / 100


experiment$figex = "siminout"
for (isPDF in c(FALSE, TRUE))
{
  if (isPDF) pdf(paste0("../Results/figs/", experiment$name, "_", experiment$figex, ".pdf"))
  plot(sort(aDV), main = "average similarity sorted - \nblack total, green inside blue outside", ylab = "point average similarity", ylim = c(0, max(c(aDV, aDVi, aDVo))));
  points(sort(aDVi), col = "green3");
  points(sort(aDVo), col = "blue3");

  if (!isPDF) readline("see ?similarities inside versus outside");
  if (isPDF) dev.off()
}


experiment$figex = "siminoutdiff"
for (isPDF in c(FALSE, TRUE))
{
  if (isPDF) pdf(paste0("../Results/figs/", experiment$name, "_", experiment$figex, ".pdf"))
  plot(sort(aDVi - aDVo), main = "average similarity in/out difference", ylab = "point average similarity", ylim = c(0, max(c(aDV, aDVi, aDVo))));

  if (!isPDF) readline("see ?similarities inside  outside differernce>>>");
  if (isPDF) dev.off()
}


experiment$figex = "toplowsim"
for (isPDF in c(FALSE, TRUE))
{
  if (isPDF) pdf(paste0("../Results/figs/", experiment$name, "_", experiment$figex, ".pdf"))
  plot(sort(aDVtop100), main = paste0("top/bottom ", topy, " point similarities"), ylim = c(0, max(aDVtop100)));
  points(sort(aDVlow100), col = "blue3");
  if (!isPDF) readline("see ?top and low similarities>>>");
  if (isPDF) dev.off()
}

experiment$figex = "toplowdiff"
for (isPDF in c(FALSE, TRUE))
{
  if (isPDF) pdf(paste0("../Results/figs/", experiment$name, "_", experiment$figex, ".pdf"))
  plot(sort(aDVtop100 - aDVlow100), main = "average differene top bottom mean", ylab = "point similarity 100 span", ylim = c(0, max(c(aDVtop100 - aDVlow100))));
  if (!isPDF) readline("see ?top low sim difference>>>");
  if (isPDF) dev.off()
}

source("cut_criteria.R");

prRCut("true", S, nnn_cats)
prNCut("true", S, nnn_cats)
trueNoCl = length(table(nnn_cats))
clrandom = sample.int(trueNoCl, dim(S)[1], replace = TRUE, prob = as.numeric(table(nnn_cats)))
prRCut("random", S, clrandom)
prNCut("random", S, clrandom)


clb = readline("\nEnter eventually the minimum average similarity diffrence required (to eliminate ambiguous docs):>>> ");
if (nchar(clb) > 2)
{
  #lb=0.01
  cat("

REDUCTION OF S \n");

  lb = as.numeric(clb);
  print(sum((aDVtop100 - aDVlow100) < lb))
  S = S[(aDVtop100 - aDVlow100) >= lb, (aDVtop100 - aDVlow100) >= lb]
  nnn_cats = nnn_cats[(aDVtop100 - aDVlow100) >= lb];
  RES0 = RES0[(aDVtop100 - aDVlow100) >= lb,]
}

print(table(nnn_cats));

##########################
# AD MAIOREM DEI GLORIAM #
##########################


cat("
###########
# THE END #
###########

call now
   source('find_Rcut.R');
");
while (sink.number() > 0) sink();


