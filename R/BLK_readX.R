if (1 == 0) { ##########################################
  source("BLK_readX.R");
} # 1===0 ##########################################
#Pseudo read of block matrix - artificial data
############ PARAMETER set samples ######################

parsamp = 3

if (parsamp == 1)
{
  group_count = 4;

  ext = 2;
  noDocs = 2000;
  overlap = 0.20;
  minprob = 0.5;
  noiseprob = 5.01;
  repetition = 6
  IsPNG = TRUE;

} # 1


if (parsamp == 2)
{
  group_count = 4;

  ext = 2;
  noDocs = 2000;
  overlap = 0.20;
  minprob = 0.5;
  noiseprob = 5.01;
  repetition = 2
  IsPNG = TRUE;

} # 2


if (parsamp == 3)
{
  group_count = 20;

  ext = 2;
  noDocs = 2000;
  overlap = 0.20;
  minprob = 0.5;
  noiseprob = 5.01;
  repetition = 1
  IsPNG = TRUE;

} # 3


cat("
############################################
Creating an artificial data setr (block-matrix) 
driven by the following parameters
############################################
group_count=,", group_count, " - number of groups/clusters
ext=", ext, "
noDocs=", noDocs, " - number of documents
overlap=", overlap, " - overlap between groups
repetition=", repetition, " - cycle in groups with distinct levels of  connectivity
minprob=", minprob, " - minimum probability
noiseprob=", noiseprob, " - noise  - factor telling how many docvument pairs from any locvation should share words
IsPNG=", IsPNG, " - if a png file will be generated ;
");
flush.console();


#############################################
while (sink.number() > 0) sink();
while (!is.null(dev.list()))  dev.off();


#-------------------------------------

targetDir = paste0("../Results/", "BLK");

if (!dir.exists(targetDir))
{
  dir.create(targetDir);
  cat("\nDirectory ", targetDir, "created\n\n");
}


#-------------------------------------
group_noDocs = ceiling(noDocs / group_count);
noWords = 5 * noDocs # 2*noDocs
group_noWords = floor(noWords / group_count);
docs = matrix(nrow = noDocs, ncol = noWords);
truegroup = rep(NA, noDocs);
docs[,] = 0;

for (j in 1:noDocs)
{
  group = floor((j - 1) / group_noDocs)
  noel = floor(group_noWords * minprob / 3 / 5)
  groupX = (group) %% repetition
  subg = 1:noel %% (group + 1 + ext); # group was applied here
  subg = 1:noel %% (groupX + 1 + ext); # groupX was applied here

  p = rnorm(noel, subg * floor(group_noWords / (group_count + ext)), group_noWords / 3 / 4)


  where = table(floor(abs(p)) %% group_noWords)
  docs[j, as.numeric(names(where))
    + group * group_noWords + 1
  ] =
    runif(length(where), minprob, 1) * where;
  truegroup[j] = group + 1;
  if ((j - 1) / group_noDocs - group < overlap)
    if (group > 0)
    { groupm1 = group - 1;
      where = sample(1:group_noWords + groupm1 * group_noWords, floor(group_noWords * minprob * overlap))
      docs[j, where] = runif(length(where), 0, minprob);
    }
  if ((j - 1) / group_noDocs - group > 1 - overlap)
    if (group < group_count - 1)
    { groupp1 = group + 1;
      where = sample(1:group_noWords + groupp1 * group_noWords, floor(group_noWords * minprob * overlap))
      docs[j, where] = runif(length(where), 0, minprob);
    }

} #

#### adding noise 

noisypoints = noiseprob * noDocs;

for (j in 1:noisypoints)
{
  d = sample(noDocs, 2);
  w = sample(noWords, 1);

  docs[d[1], w] = docs[d[1], w] + runif(1, 0, minprob)
  docs[d[2], w] = docs[d[2], w] + runif(1, 0, minprob)
}


# graph matrix presentation

S = docs %*% t(docs)
D = diag(S)
Ds = sqrt(D);
Ds[Ds == 0] = 1000;
Dv = Ds %*% t(Ds)
S = S / Dv;

for (j in 1:dim(S)[1])S[j, j] = 0;

tr = quantile(S, prob = 0.999);

z = which(S > tr)
x = (z - 1) %% dim(S)[1] + 1;
y = (z - x) / dim(S)[1];
plot(x, y, main = paste("nodes connected with weight above", round(tr, digit = 3)))
if (IsPNG) png(paste0(targetDir, "/", "BlockStructure", round(tr, digit = 3), "_", group_count, "GR.png"));

plot(x, y, main = paste("nodes connected with weight above", round(tr, digit = 3)))
if (IsPNG) dev.off() else
  readline("See the connection matrix outline")
if (IsPNG) png(paste0(targetDir, "/", "S_BlockStructure", round(tr, digit = 3), "_", group_count, "GR.png"));
hist(S)
if (IsPNG) dev.off() else
  readline("See the similarity distribution")
if (IsPNG) png(paste0(targetDir, "/", "posS_BlockStructure", round(tr, digit = 3), "_", group_count, "GR.png"));
hist(S[S > 0])
print(quantile(S[S > 0]))
if (IsPNG) dev.off() else
  readline("See the similarity distribution")


colnames(docs) = paste0("W", 1:dim(docs)[2])

######### THE RESULTS ################### 

# nnn.type
nnn.type = paste0("BLK.", group_count, "_", overlap, "_", minprob);
#  nnn_cats ; 
nnn_cats = truegroup
# RES0; -- 
RES0 = docs;


#############################################

cat('
#############################
         THE END",blk.type,"
#############################
call now 
	source("S_describeBis.R");
');


#############################################

cat('
#############################
         THE END Block Data 
#############################

');

##########################
# AD MAJOREM DEI GLORIAM #
##########################


comb_reslist = NULL;
norm_reslist = NULL;

