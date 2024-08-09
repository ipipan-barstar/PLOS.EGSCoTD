
# About
Python and R source code for the paper _"Explainable Graph Spectral Clustering of Text Documents"_
by Bartłomiej Starosta, Mieczysław A. Kłopotek, Sławomir T. Wierzchoń, Dariusz Czerski, Marcin Sydow and Piotr Borkowski
published in PLOSONE.

This code was used to produce the tables presented in the paper.


# Main Files and Folders

## ``Data/Tweets.10tags``

Input data: tweets with exactly 1 of 10 hashtags.

## ``Results``

Output directory ``Results`.
Must exist in parent relatively to R code working folder.

## ``R``

Run the scripts interactively and follow the inctructions as they appear.
Optionally fine-tune parameters in files prior to execution. 
Results are stored in folder ``Results``.

### ``TWT_read_exp.R``

Clustering real world data (``Data/Tweets.10tags``). 

### ``BLK_readX.R``

Clustering artificially generated data.

## ``Python`` 

### ``spherical-kmeans.py``

Clustering with spherical k-means algorithm.

### ``scikt-spectral.py``

Clustering with spectral algorithm.



# Requirements

1. ``Python >= 3.7``

2. ``scikit-learn <= 0.24.2``

3.  ``SphericalKMeans`` from ``soyclustering (https://github.com/lovit/clustering4docs)``