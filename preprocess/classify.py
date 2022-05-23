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
import pandas as pd

def kdecutoff(mappingvalues):
    """Used to determine minimum mapping rate cutoff (threshold) for quality control based on their distribution (more specifically, estimated kernel density)."""
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

def thresholder(maprate_dict, cutoff):
    """wrapper function for kdecutoff(). takes in the dictionary of {accession:PS%,...} in log file.
    upacks the dictionary and return lists of total, failed and passed accession as well"""
    if cutoff==0:
        cutoff= kdecutoff(list(maprate_dict.values()))
    total= list(maprate_dict.keys())
    failed= [accession for accession, maprate in maprate_dict.items() if maprate < cutoff]
    passed= [accession for accession, maprate in maprate_dict.items() if maprate > cutoff or maprate == cutoff]
    return total, failed, passed, cutoff

def kmeans_kwalk(data, k_min, k_max):
    """do kmeans iteration walk"""
    kmeans_kwargs = {"init": "k-means++", "n_init": 100,"max_iter": 2000,"random_state": 42} #remove random state (seed) after development
    silhouette_coefficients = []
    k_cluster_assignment_dict={}
    for k in range(k_min,k_max):
        kmeans = KMeans(n_clusters=k, **kmeans_kwargs)
        kmeans.fit(data)
        silhouette_coefficients.append(silhouette_score(data, kmeans.labels_))
        k_cluster_assignment_dict[k]= kmeans.labels_
        print(f"K-means iteration at k={k} complete.\n")
    return k_cluster_assignment_dict , silhouette_coefficients
    

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
    Matrix= pd.read_csv( pathtomat,sep="\t").set_index("accession").transpose()
    Matrix= Matrix[accessions]
    return Matrix
    
def optimal_k_silhouette(k_min, k_max, silhouette_coefficients, k_cluster_assignment_dict):
    """Determine the optimal k based on silhouette score. Returns optimal k and respective list of cluster_labels"""
    optimal_k = [k for k in range(k_min,k_max)][silhouette_coefficients.index(max(silhouette_coefficients))]
    cluster_assignment = k_cluster_assignment_dict.get(optimal_k)
    return optimal_k, cluster_assignment , max(silhouette_coefficients)

def report_cluster_assignment_stats(cluster_assignment_dict):
    n_accession_list= [len(val) for val in cluster_assignment_dict.values()]
    median_stat = np.median(n_accession_list)
    mean_stat = np.mean(n_accession_list)
    min_stat = np.max(n_accession_list)
    max_stat = np.min(n_accession_list)
    return int(median_stat) , int(mean_stat), int(min_stat) , int(max_stat)
    