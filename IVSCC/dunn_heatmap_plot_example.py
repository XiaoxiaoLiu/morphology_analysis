__author__ = 'xiaoxiaol'

from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import squareform
import matplotlib.pyplot as plt
import numpy as np


def clusterDistances(cluster1, cluster2, DM):
    return DM[cluster1][:, cluster2].mean()


def Dunn(clustering, DM):
    clusters = len(clustering)
    maxIntra = max([clusterDistances(c, c, DM) for c in clustering])
    minInter = np.inf
    for cl1 in range(1, clusters):
        for cl2 in range(cl1):
            minInter = min(clusterDistances(clustering[cl1], clustering[cl2], DM), minInter)
    return minInter/maxIntra


def flatclusterDunn(D):
    DM = squareform(D)
    Y = linkage(DM, method="average")
    n = len(Y)+1
    clusterdict = dict([(el, [el]) for el in range(n)])
    bestClustering, bestDunn = None, 0
    for cidx, cluster in enumerate(Y):
        cl1, cl2 = map(int, cluster[:2])
        clusterdict[cidx+n] = clusterdict[cl1] + clusterdict[cl2]
        del clusterdict[cl1]
        del clusterdict[cl2]
        clustering = clusterdict.values()
        if not len(clustering) in [1, n]:
            currentDunn = Dunn(clustering, DM)
            print (clustering, currentDunn)
            if currentDunn > bestDunn:
                bestDunn = currentDunn
                bestClustering = clustering
    print "The Best clustering with a maximal Dunn index is", bestClustering, bestDunn

def heatmap(D):
    temp = D
    D = squareform(D)
    Y = linkage(D, method='average')
    fig = plt.figure(figsize=(8, 8))
    dendogram = fig.add_axes([0.01, 0.1, 0.1, 0.8])
    Z = dendrogram(Y, orientation='right')
    idx = Z['leaves']
    dendogram.set_yticklabels(idx)
    D = D[idx, :][:, idx]
    matrix = fig.add_axes([0.37, 0.1, 0.6, 0.8])
    im = matrix.matshow(D, aspect='auto', origin='lower', cmap=plt.cm.YlGnBu)
    matrix.set_xticklabels([])
    matrix.set_yticklabels([])
    plt.show()
    flatclusterDunn(temp)


def main():
    D = [1, 5, 6, 0.9, 4, 5, 0.1, 1, 4.1, 5.1]
    heatmap(D)

main()
