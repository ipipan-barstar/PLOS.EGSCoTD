# -*- encoding: utf-8 -*-

from sklearn.cluster import SpectralClustering
from sklearn.metrics.pairwise import cosine_similarity
from common import *


data_file = "Data/Tweets.10tags"
repeat = 10

# read input and produce list of hashtags and vectors of term counts
hashtags, vectors = read_data(data_file)

# mapping of hashtags to hashtag count, list of hashtags as indices (true clustering), hashtag count (== 10)
tagmap, hashids, cluster_cnt = mktagmap(hashtags)


def repeat_spectral_clustering(rpt, clustering_object, data):
    """
    Repeat spectral clustering, collect results

    Parameters
    --------
    rpt : number of repetitions
    clustering_object : spectral clustering object
    data : input data

    Returns
    -------
    result : list of clustering results

    """

    result = []
    for r in range(rpt):
        clustering = clustering_object.fit_predict(data)
        result.append(clustering)
    return result


for affinity in ('precomputed', 'nearest_neighbors'):
    sc = SpectralClustering(n_clusters=cluster_cnt, affinity=affinity)
    if affinity == 'precomputed':
        matrix = cosine_similarity(vectors)
    else:
        matrix = vectors
    results = repeat_spectral_clustering(repeat, sc, matrix)
    clusters, best, worst, avg = select_best_f(hashids, results)
    print("Affinity: %s, repeat = %d" % (affinity, repeat))
    print("F-score: best=%.6f, avg=%.6f, worst=%.6f" % (best, avg, worst))
    print(evaluate_score(clusters, hashids))

