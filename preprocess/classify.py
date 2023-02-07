import sys
import os
if __name__ == "__main__":
        abspath= os.getcwd()
        parent_module= os.path.join(abspath.split("LSTrAP-denovo")[0], "LSTrAP-denovo")
        sys.path.insert(0, parent_module)

import numpy as np        
from scipy.signal import argrelextrema
from sklearn.decomposition import PCA
from sklearn.neighbors import KernelDensity
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from sklearn_extra.cluster import KMedoids
import pandas as pd
from scipy.stats import iqr
from sklearn.metrics.pairwise import euclidean_distances

def kdecutoff(mappingvalues):
    """Used to determine minimum mapping rate cutoff (threshold) for quality control based on their distribution (more specifically, estimated kernel density).
    DEPRECATED"""
    valuearray= np.array(mappingvalues).reshape(-1,1)
    #init kde object using mapping values
    kde = KernelDensity(kernel='gaussian', bandwidth= 5).fit(valuearray)
    #linspace: Returns num evenly spaced samples, calculated over the interval [start, stop].
    intervals= np.linspace(np.amin(valuearray),np.amax(valuearray))
    #get kde probablities at the evenly spaced maprate intervals
    probabilities = kde.score_samples(intervals.reshape(-1,1))
    #find valleys in the kde function
    minima = argrelextrema(probabilities, np.less)[0]
    #return minimum mapping rate if there are no valleys present
    if intervals[minima].size ==0:
        return np.amin(valuearray)
    #else: return value of the first valley as threshold
    else:
        return intervals[minima][-1]
        
def lowerfence_iqr_cutoff(mappingvalues):
    '''replaces kdecutoff (more reliable and is a well established statistical method to look for outliers). Basically returns the lower fence (Q1 - 1.5*IQR) using the interquartile range (IQR)
    if lower fence is lower than 20%, return 10 % instead'''
    return max([np.percentile(mappingvalues, 25) - 1.5*iqr(mappingvalues), 20])
    

def thresholder(maprate_dict, cutoff):
    """wrapper function for kdecutoff(). takes in the dictionary of {accession:PS%,...} in log file.
    upacks the dictionary and return lists of total, failed and passed accession as well"""
    if cutoff==0:
        cutoff= lowerfence_iqr_cutoff(list(maprate_dict.values()))
    total= list(maprate_dict.keys())
    failed= [accession for accession, maprate in maprate_dict.items() if maprate < cutoff or maprate == cutoff]
    passed= [accession for accession, maprate in maprate_dict.items() if maprate > cutoff ]
    return total, failed, passed, cutoff

def kmeans_kwalk(data, k_min, k_max):
    """do kmeans iteration walk"""
    kmeans_kwargs = {"init": "k-medoids++", "method":"pam" ,"max_iter": 2000,"random_state": 42}
    silhouette_coefficients = []
    k_cluster_assignment_dict={}
    k_centroids_dict={}
    for k in range(k_min,k_max):
        kmeans = KMedoids(n_clusters=k, **kmeans_kwargs)
        kmeans.fit(data)
        silhouette_coefficients.append(silhouette_score(data, kmeans.labels_))
        k_cluster_assignment_dict[k]= kmeans.labels_
        k_centroids_dict[k]= kmeans.cluster_centers_
        print(f"K-means iteration at k={k} complete.\n")
    return k_cluster_assignment_dict , silhouette_coefficients, k_centroids_dict
    

def PCA_transformer(Matrix, n_pcs=100):
    """normalize Matrix, transpose into dataframe"""
    #normalize colummns (within samples)
    sample_normalized_matrix = pd.DataFrame(StandardScaler().fit_transform(Matrix).T, 
                                            index=Matrix.columns,
                                            columns=Matrix.index)
    #PCA transformation
    pca = PCA()
    pca_data=pca.fit_transform(sample_normalized_matrix)
    pca_data = pd.DataFrame(pca_data, index= Matrix.columns)
    #keeping ony desired pcs
    pca_data = pca_data[pca_data.columns[:n_pcs]]
    #calculate total explained variance in the pcs that were kept
    pc_variances = np.round(pca.explained_variance_ratio_*100, decimals=1)
    
    return pca_data , pc_variances

def mat_parser(pathtomat, accessions):
    ''' returns matrix for PCA x-formation. accessions= list of accessions that passed qc'''
    Matrix= pd.read_csv( pathtomat,sep="\t").set_index("accession")
    Matrix = Matrix[~Matrix.index.duplicated(keep='first')].transpose()
    Matrix= Matrix[accessions]
    return Matrix
    
def optimal_k_silhouette(k_min, k_max, silhouette_coefficients, k_cluster_assignment_dict, k_centroids_dict):
    """Determine the optimal k based on silhouette score. Returns optimal k and respective list of cluster_labels.
    Optimal K is defined as the lowest k  within 0.01 range of max silhoette score"""
    sc_rounded= [np.round(sc, 2) for sc in silhouette_coefficients] 
    optimal_k = [k for k in range(k_min,k_max)][sc_rounded.index(max(sc_rounded))]
    cluster_assignment = k_cluster_assignment_dict.get(optimal_k)
    centroids = k_centroids_dict.get(optimal_k)
    return optimal_k, cluster_assignment , silhouette_coefficients[sc_rounded.index(max(sc_rounded))] , centroids

def report_cluster_assignment_stats(cluster_assignment_dict):
    n_accession_list= [len(val) for val in cluster_assignment_dict.values()]
    median_stat = np.median(n_accession_list)
    mean_stat = np.mean(n_accession_list)
    min_stat = np.min(n_accession_list)
    max_stat = np.max(n_accession_list)
    return int(median_stat) , int(mean_stat), int(min_stat) , int(max_stat)
 
def generate_master_cluster_assignment_dict(cluster_assignment, centroids, pca_data, mapratedict):
    """Step after optimal_k_silhouette().
    used to generate master cluster assignment dict containing cluster as keys and subdict as values 
    subdict have accession as keys and info tuple as values sorted based on distance to cluster centroids.
    info list : (distance, maprate ) 
    e.g. {"0": {"SRR123456":[dist, maprate], ...}, ...}"""
    temp_dict = {}
    for cluster, accessions in cluster_assignment.items():
        temp_dict[cluster] = {}
        for acc in accessions:
            min_euc_dist= euclidean_distances(np.array(centroids), np.array(pca_data.loc[acc]).reshape(1, -1))[int(cluster)]
            temp_dict[cluster][acc]=float(min_euc_dist)
    master_cluster_assignment_dict = {}
    for cluster, dist_dict in temp_dict.items():
        master_cluster_assignment_dict[cluster]={}
        sortedlist= sorted(list(dist_dict.values()))
        inverted_dist_dict={value:key for key, value in dist_dict.items()}
        for value in sortedlist:
            master_cluster_assignment_dict[cluster][inverted_dist_dict[value]]=[ value, mapratedict[inverted_dist_dict[value]]]
    return master_cluster_assignment_dict
    
def mat_transposer(pathtomat,  accessions ,pathtomat_T):
    '''only include accessions that passed + 
    transpose expression matrix such that columns are RNA-seq accessions, index corresponds to genes'''
    try:
        Matrix = mat_parser(pathtomat, accessions)
        Matrix.to_csv(pathtomat_T ,sep="\t")
        return 1
    except:
        return 0