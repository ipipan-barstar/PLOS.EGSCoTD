if (1 == 0) { ##########################################
  source("TWT_read_exp.R");
  # reading only those from a predefined list of tags
} # 1===0 ########################################## 
############ PARAMETERS ######################

experiment = list()
experiment$name = "TWT3ht"
experiment$figex = ""


twt.type.list = c("EN");
twt.type.id = 1;  # current selection 1:2

twt.type = twt.type.list[twt.type.id];
twt.ver = "tags";
#IsPNG=FALSE
IsPNG = TRUE

linesToRead = 220000 # the first linesToRead tweets will be read

doccountmax = 2200 # the maximum number of docs to consider


theFolder = "../Data" # place where files are located
dataset.list = c(
  "Tweets.10tags"
);

EN.hashes = c(
  "#puredoctrinesofchrist"
  , "#lolinginlove" ### remved temperarily in some experiments
  , "#anjisalvacion"
  , "#nowplaying" ### remved temperarily  in some experiments
)


#------------------------------

############ END PARAMETERS ######################


#############################################
while (sink.number() > 0) sink();
while (!is.null(dev.list()))  dev.off();

#-------------------------------------
alphabetizeENG <- function(txt)
{
  ltxt = nchar(txt);
  wdl = unlist(strsplit(txt, " "));
  wdl2 = c();
  for (w in wdl)
  {
    wd = unlist(strsplit(w, ""));
    wd2 = wd[(wd >= 'a' & wd <= 'z') | wd == '#'];
    wd3 = c();
    for (l in wd2)
      if (!is.na(l))
      { x = utf8ToInt(l);
        if (l == " ")  wd3 = c(wd3, l) else
          if (!is.na(x)) wd3 = c(wd3, l) else
            wd3 = c(wd3, " ")
      }
    wd4 = paste(wd3, collapse = "");
    wdl2 = c(wdl2, wd4)
  } # for w
  return(paste(wdl2, collapse = " "));
}

# eof


#------------------------------
#-------------------------------------
#------------------------------
cat("
Reading ", dataset.list[twt.type.id], "
");
if (twt.type == "EN")
{
  twEN =
    readLines(paste0(theFolder, "/", dataset.list[twt.type.id])
      , encoding = "UTF-8", n = linesToRead
    );
  cat("
Cleaning ", dataset.list[twt.type.id], "
");
  #readline("tst");
  for (jr in 1:length(twEN))
  {
    txt = twEN[jr];
    q = floor(linesToRead / sqrt(linesToRead) - 0.0001) + 1;
    if (jr %% q == 0)
    { cat("."); flush.console(); }
    txt = alphabetizeENG(txt);


    twEN[jr] = txt;
  } # for jr

  # removing the tags
  twENnt = rep(NA, length(twEN))
  twENt = rep(NA, length(twEN))
  tl = c();
  for (jr in 1:length(twEN))
  {
    tst = twEN[jr];
    tst2 = unlist(strsplit(tst, " "))
    tst3 = tst2[nchar(tst2) >= 1]
    tstnh = tst3[substr(tst3, 1, 1) != '#']
    tsth = tst3[substr(tst3, 1, 1) == '#']
    twENnt[jr] = tolower(paste(tstnh, collapse = ' '))
    tsth = tsth[nchar(tsth) > 2]
    twENt[jr] = tolower(paste(sort(tsth), collapse = ' '))
    tl = c(tl, tsth);
  }

  cat("\n");
  goodhashes = twENt %in% EN.hashes


  flush.console();

  twt_docs = twENnt[goodhashes]
  twt_cats = twENt[goodhashes]
  print(table(twt_cats))

} # tw EN 

#------------------------------
#------------------------------

# computing tweet length

mitw = Inf
matw = 0

for (j in  twt_docs)
{ nc = nchar(j)
  if (mitw > nc) mitw = nc
  if (matw < nc) matw = nc

}


docsetlist = NULL;

# splitting the documents into subsets 

if (length(twt_docs) > doccountmax)
{
  docsetlist = list()
  dsc = ceiling(length(twt_docs) / doccountmax)
  ds_size = ceiling(length(twt_docs) / dsc)
  di = sample(length(twt_docs), length(twt_docs))
  for (j in 1:dsc)
  {
    fr = (j - 1) * ds_size + 1
    to = j * ds_size
    #fr= 1
    #to= ds_size
    if (to > length(twt_docs)) to = length(twt_docs)
    j_twt_docs = twt_docs[di[fr:to]]
    j_twt_cats = twt_cats[di[fr:to]]
    docsetlist[[j]] = list()
    docsetlist[[j]]$twt_docs = j_twt_docs
    docsetlist[[j]]$twt_cats = j_twt_cats

    cat("
--------------------------------
Portion number ", j, "

");
    print(table(j_twt_cats))
    cat("
number of twt_docs ", length(j_twt_docs), "
");

  }
  twt_docs = docsetlist[[1]]$twt_docs
  twt_cats = docsetlist[[1]]$twt_cats
}

#------------------------------
cat('
call now
	source("TWT_doc_matrix.R");
');
#############################################

comb_reslist = NULL;
norm_reslist = NULL;


############################################
############################################
##########                        ##########
########## AD MAIOREM DEI GLORIAM ##########
##########                        ##########
############################################
############################################

