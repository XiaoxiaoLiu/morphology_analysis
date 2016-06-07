__author__ = 'xiaoxiaol'
from scipy.cluster import hierarchy
import numpy as np
import random
import pylab as plt

def dispersion (data,n_clusters):
    linkage =  hierarchy.linkage(data, method='ward', metric='euclidean')
    assignments = hierarchy.fcluster(linkage, n_clusters, criterion="maxclust")

    if  n_clusters== 1:
        cluster_mean = np.mean(data, axis=0)
        distances_from_mean = np.sum((data - cluster_mean)**2,axis=1)
        dispersion_val = np.log(sum(distances_from_mean))
    else:

        distances_from_mean = range(n_clusters)
        for i in range(n_clusters):
            ids = np.nonzero(assignments == i)[0]  # starting from  0
            df_cluster = data[ids]
            cluster_mean= df_cluster.mean(axis=0)
            distances_from_mean[i] = int()
            for idx in ids:
                distances_from_mean[i] += sum((data[idx] - cluster_mean)**2)
        dispersion_val = np.log(sum(distances_from_mean))

    return dispersion_val

def reference_dispersion(data, num_clusters, num_reference_bootstraps):
    dispersions = [dispersion(generate_uniform_points(data), num_clusters) for i in range(num_reference_bootstraps)]
    mean_dispersion = np.mean(dispersions)
    stddev_dispersion = float(np.std(dispersions)) / np.sqrt(1. + 1. / num_reference_bootstraps)
    return mean_dispersion





def generate_uniform_points(data):
    mins = np.argmin(data, axis=0)
    maxs = np.argmax(data, axis=0)

    num_dimensions = data.shape[1]
    num_datapoints = data.shape[0]

    reference_data_set = np.zeros((num_datapoints,num_dimensions))
    for i in range(num_datapoints):
        for j in range(num_dimensions):
            reference_data_set[i][j] = random.uniform(data[mins[j]][j],data[maxs[j]][j])

    return reference_data_set




def gap_statistic (data, nthCluster, num_reference_bootstraps=10):
    actual_dispersion = dispersion(data, nthCluster)
    ref_dispersion = reference_dispersion(data, nthCluster, num_reference_bootstraps)
    return actual_dispersion, ref_dispersion



def gap_analysis (data, low=5, high =15,num_reference_bootstraps = 10, output_dir="."):

    dispersion_values = np.zeros((high,2))

    for cluster in range(low, high+1):
        dispersion_values_actual,dispersion_values_reference = gap_statistic(data, cluster, num_reference_bootstraps)
        dispersion_values[cluster-1][0] = dispersion_values_actual
        dispersion_values[cluster-1][1] = dispersion_values_reference

    gaps = dispersion_values[:,1] - dispersion_values[:,0]

    print "gpas:",gaps
    print "The estimated number of clusters is ", range(high)[np.argmax(gaps)]+1
    plt.figure()
    plt.plot(range(len(gaps)), gaps,"*-")
    plt.xlabel("cluster number")
    plt.ylabel("gap")
    #plt.show()
    plt.savefig(output_dir+'/gap.png')
    plt.close()