# -*- encoding: utf-8 -*-

import codecs
import sklearn.metrics as sm
from sklearn.feature_extraction.text import CountVectorizer


def read_data(src_file):
    vectorizer = CountVectorizer(strip_accents='unicode')
    with codecs.open(src_file, mode='r', encoding='utf-8') as input_file:
        input_lines = [line.rstrip() for line in input_file.readlines()]
        hash_split = [line.split("#")[1] for line in input_lines]
        hashtags = [line.split(None)[0] for line in hash_split]
        vectors = vectorizer.fit_transform(input_lines)
        return hashtags, vectors


# Since Python 3.7*, dictionaries are order-preserving,
def mktagmap(hashtags):
    tagmap = {}
    for tag in hashtags:
        if tag in tagmap:
            tagmap[tag] = 1 + tagmap[tag]
        else:
            tagmap[tag] = 1
    taglist = sorted(tagmap, key=tagmap.get)
    cluster_cnt = len(taglist)
    hashids = [taglist.index(tag) for tag in hashtags]
    return tagmap, hashids, cluster_cnt


def evaluate_score(labels_true, labels_pred):
    import pkg_resources
    smv = pkg_resources.get_distribution("scikit-learn").version
    v = int(smv[2:4])
    scores = "" \
             + "adjusted_mutual_info_score:   \t%f\n" % (sm.adjusted_mutual_info_score(labels_true, labels_pred)) \
             + "adjusted_rand_score:          \t%f\n" % (sm.adjusted_rand_score(labels_true, labels_pred)) \
             + "completeness_score:           \t%f\n" % (sm.completeness_score(labels_true, labels_pred)) \
             + "fowlkes_mallows_score:        \t%f\n" % (sm.fowlkes_mallows_score(labels_true, labels_pred)) \
             + "homogeneity_score:            \t%f\n" % (sm.homogeneity_score(labels_true, labels_pred)) \
             + "mutual_info_score:            \t%f\n" % (sm.mutual_info_score(labels_true, labels_pred)) \
             + "normalized_mutual_info_score: \t%f\n" % (sm.normalized_mutual_info_score(labels_true, labels_pred)) \
             + "v_measure_score:              \t%f\n" % (sm.v_measure_score(labels_true, labels_pred))
    if v > 22:
        scores = scores + ("rand_score:                   \t%f\n" % (sm.rand_score(labels_true, labels_pred)))
    return scores


def select_best_f(lab_true, lab_pred_list):
    best = -1.0
    worst = 1.0
    avg = 0.0
    best_cl = None
    for lab_pred in lab_pred_list:
        score = sm.fbeta_score(lab_true, lab_pred, average='macro', beta=2)
        if score > best:
            best = score
            best_cl = lab_pred
        if score < worst:
            worst = score
        avg += score
    return best_cl, best, worst, avg / len(lab_pred_list)
