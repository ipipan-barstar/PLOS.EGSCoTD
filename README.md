
# About
Python and R source code for the paper _"Explainable Graph Spectral Clustering of Text Documents"_
by Bartłomiej Starosta, Mieczysław A. Kłopotek, Sławomir T. Wierzchoń, Dariusz Czerski, Marcin Sydow and Piotr Borkowski
published in PLOSONE.

This code was used to produce the tables presented in the paper.


# Files and Folders

## Python 

### ``spherical-kmeans.py``

Clustering with spherical k-means algorithm.

### ``scikt-spectral.py``

Clustering with spectral algorithm.

### ``common.py``

Common python procedures.

## Data

### ``Tweets.10tags``

Input data: tweets with exactly 1 of 10 hashtags.

# Requirements

1. ``Python >= 3.7``

2. ``scikit-learn <= 0.24.2``

3.  ``SphericalKMeans`` from ``soyclustering (https://github.com/lovit/clustering4docs)``