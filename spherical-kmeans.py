# -*- encoding: utf-8 -*-

from soyclustering import SphericalKMeans
from common import *


data_file = "Data/Tweets.10tags"
repeat = 1

hashtags, vectors = read_data(data_file)
tagmap, hashids, cluster_cnt = mktagmap(hashtags)


skmeans_config = {
    "sc.n": ('similar_cut', None),
    "sc.sc": ('similar_cut', 'sculley'),
    "sc.md": ('similar_cut', 'minimum_df'),
    "k++.n": ('k-means++', None),
    "k++.sc": ('k-means++', 'sculley'),
    "k++.md": ('k-means++', 'minimum_df')
}


def run_skmeans(spher_kmeans, data_matrix, rpt):
    result = []
    for n in range(rpt):
        clus = spher_kmeans.fit_predict(data_matrix)
        result.append(clus)
    return result


for config in skmeans_config:
    init, sparsity, = skmeans_config[config]
    km = SphericalKMeans(n_clusters=cluster_cnt, max_iter=100, init=init, sparsity=sparsity)
    results = run_skmeans(km, vectors, repeat)
    clusters, best, worst, avg = select_best_f(hashids, results)
    print("Sk-means config = %s: (init=%s, sparsity=%s), repeat = %d" % (config, skmeans_config[config][0], skmeans_config[config][1], repeat))
    print("F-score: best=%.6f, avg=%.6f, worst=%.6f" % (best, avg, worst))
    print(evaluate_score(hashids, clusters))


