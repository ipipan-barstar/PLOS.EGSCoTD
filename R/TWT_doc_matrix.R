if (1 == 0) { ##########################################
  source("TWT_doc_matrix.R");
} # 1===0 ##########################################
############ PARAMETERS ######################

ENexclusionList = c();

# is the the document matrix to be saved
# SaveDocMatrix=TRUE;
SaveDocMatrix = FALSE;


############ END PARAMETERS ######################


#############################################
while (sink.number() > 0) sink();
while (!is.null(dev.list()))  dev.off();
#library(openxlsx);
#library(tidytext)


cat('
Turning  the TWT set ', twt.type, ' to a document matrix.
Assumption -- the data is already in memory via 	source("TWT_read.R" );
'); flush.console();
#-------------------------------------
#------------------------------
cat("
Processing  ", dataset.list[twt.type.id], "
");
if (twt.type == "EN")
{
  exclusionList = ENexclusionList;
} #


cat("
Original Data Dimensionality -- number of records\n");
print(length(twt_docs));
flush.console();
words = c();
pairedwords = c();

for (j in 1:length(twt_docs))
{ txt = twt_docs[j];
  if (j %% 100 == 0) { cat("."); flush.console(); }
  txt = tolower(txt);
  txt = gsub("?", " ", txt, fixed = TRUE);
  txt = gsub(":", " ", txt, fixed = TRUE);
  txt = gsub("\\.", " ", txt, fixed = TRUE);
  wds = unlist(strsplit(txt, " ", fixed = TRUE));
  wds = wds[nchar(wds) > 0];
  wds = setdiff(wds, exclusionList);
  wds = wds[is.na(as.numeric(wds))];


  wp = c();
  for (js in 2:length(wds))
  { wp = c(wp, paste(wds[js - 1], wds[js]));
  }


  wds = paste0("WD", wds);
  wp = paste0("WP", wp);

  words = union(words, wds);
  pairedwords = union(pairedwords, wp);
}

words = words[nchar(words) > 0];
words = sort(words);
pairedwords = sort(pairedwords);

cat(" 
True Data Dimensionality: records ", length(twt_docs)
  , "words", length(words), "word pairs", length(pairedwords), "\n");
flush.console();

hist_pairedwords = pairedwords;
pairedwords = c();
noDocs = length(twt_docs);
docM = matrix(ncol = length(words) + length(pairedwords), nrow = length(twt_docs));
colnames(docM) = make.names(c(words, pairedwords));

docM[,] = 0;
#readline("tst"); 
for (j in 1:length(twt_docs))
{ txt = twt_docs[j];
  txt = tolower(txt);
  txt = tolower(txt);
  txt = gsub("?", " ", txt, fixed = TRUE);
  txt = gsub(":", " ", txt, fixed = TRUE);
  txt = gsub("\\.", " ", txt, fixed = TRUE);
  wds = unlist(strsplit(txt, " ", fixed = TRUE));
  wds = wds[nchar(wds) > 0];
  wds = wds[is.na(as.numeric(wds))];
  wds = wds[nchar(wds) > 0];
  wds = setdiff(wds, exclusionList);

  wp = c();
  if (length(wds) > 2)
    for (js in 2:length(wds))
    { wp = c(wp, paste(wds[js - 1], wds[js]));
    }

  wds = paste0("WD", wds);
  wp = paste0("WP", wp);
  wds = make.names(wds);
  wp = make.names(wp);
  #words=union(words,wds);
  #pairedwords=union(pairedwords,wp);
  wds = union(wds, c());
  wp = union(wp, c());
  if (length(wds) > 0)
    docM[j, wds] = 1;
  # docM[j,wp]=1;
}

#print(docM[1:10,1:10] );

docM = as.data.frame(docM);
RES = docM;

cat("
Remove columns with one non-zero entry only 

"); flush.console();
toremove = c();
for (jc in colnames(RES))
{
  ile = sum(RES[, jc] > 0);
  if (ile < 2)    toremove = c(toremove, jc);
} # for jc
cat("
Columns with one entry only
");
#print(head(toremove)); 
flush.console();
toremove = c(toremove,
             "WDthe", "WDto", "WDand", "WDof", "WDa", "WDis", "WDin", "WDthat", "WDthis", "WDi",
             "WDversion", "WDkjv", "WDjames",
             "WDfor", "WDare", "WDyou", "WDbe", "WDnot", "WDit", "WDhave", "WDwe", "WDunto"
  , "WDhe", "WDbut", "WDhis", "WDyou"
  , "WDall", "WDon", "WDwith", "WDif", "WDso", "WDthey", "WDwho", "WDmy"
  , "WDstandard", "WDenglish", "WDthey", "WDshall"
  , "WDwhat", "WDby", "WDabout", "WDjust", "WDdo", "WDtheir", "WDyour"
  , "WDhim", "WDwhich", "WDthem", "WDwill", "WDyour"
  , "WDas", "WDcan", "WDor", "WDno", "WDits", "WDhow"
  , "WDye", "WDfrom", "WDthou", "WDat"
  , "WDwas", "WDbecause", "WDup", "WDour", "WDwhen", "WDone", "WDmany", "WDme"
  , "WDdont", "WDbeing", "WDim", "WDan", "WDus", "WDout", "WDthere"
  , "WDnow", "WDlike", "WDmore", "WDthese", "WDthy"
  , "WDhas", "WDhave", "WDgreat", "WDhey", "WDwilling", "WDmake", "WDplease", "WDwant"
  , "WDget", "WDsounds", "WDusing", "WDwont", "WDafter", "WDagain", "WDsaid"
)

colstoleave = setdiff(colnames(RES), toremove)
RES = RES[, colstoleave];
rownames(RES) = c()


cat(" 
Corrected True Data Dimensionality: ", length(twt_docs), length(colnames(RES)), "\n");
flush.console();


res = list();
res$RES = RES;
res$nnn_cats = twt_cats;
res$nnn.type = paste0("TWT.", twt.type);


if (SaveDocMatrix) {

  basepath = paste0("../Results/", substr(res$nnn.type, 1, 3))
  ifelse(!dir.exists(file.path(getwd(), basepath)), dir.create(getwd(), basepath), FALSE)


  wb <- createWorkbook("TWT words")
  sheetname = "WORDS_data";
  addWorksheet(wb, sheetname);
  writeData(wb, sheet = sheetname, RES)
  sheetname = "CLSS_data";
  addWorksheet(wb, sheetname);
  writeData(wb, sheet = sheetname, twt_cats)
  saveWorkbook(wb, file = paste0(basepath, "/", twt.type, "allWordColsClasses.xlsx"), overwrite = TRUE)
  cat("

 saved as", paste0("../Results/TWT/", twt.type, "allWordColsClasses.xlsx")
    , "
");
}

#saveRDS(res, file = paste0("../Results/TWT/",twt.type,"allwordcols.RDS") )

cat("
 
   saved as", paste0("../Results/TWT/", twt.type, "allwordcols.RDS")
  , "


 ");
flush.console();


######### THE RESULTS ################### 

#  nnn_cats ; 
nnn_cats = twt_cats
# nnn.type
nnn.type = paste0("TWT.", twt.type, length(table(nnn_cats))) #,".",twt.ver);
# RES0; -- 
RES0 = RES;
S = as.matrix(RES[,]) %*% t(as.matrix(RES[,]))
#if (COSINEDIST) 
{
  D = diag(S)
  D[D == 0] = 1
  Ds = sqrt(D);
  Dv = Ds %*% t(Ds)
  S = S / Dv;
}


sink("../Results/DATASET_descriptions.txt", append = TRUE)
nrec = length(nnn_cats)
cardcats = table(nnn_cats)
ncl = length(cardcats)
cat("
------------------------------------
", nnn.type, " dataset - ", nrec, "records, ", ncl, "classes,
 named", paste(names(cardcats), collapse = ", "), "
with minimal cardinality  ", min(cardcats),
    "  and  maximal cardinality  ", max(cardcats), "
- choice from manually selected tweet tag list  for tweet lengths without tags 
min. ", mitw, " (which implied  max. length ", matw, ")
")
#print(cardcats); 
sink()


S = S - diag(diag(S))

#############################################

cat('
#############################
         THE END', nnn.type, '
#############################
call now 
    source("S_describeBis.R");
');

############################################
############################################
##########                        ##########
########## AD MAIOREM DEI GLORIAM ##########
##########                        ##########
############################################
############################################


