#setting sys.path for importing modules
import os
import sys
if __name__ == "__main__":
        abspath= os.getcwd()
        parent_module= os.path.join(abspath.split("HSS-Trans")[0], "HSS-Trans")
        sys.path.insert(0, parent_module)

import concurrent.futures
import numpy as np
from setup import constants
from preprocess import read_map, classify
from assembly import misc
import pandas as pd

tpm_matpath="/home/pengken001/pipeline_output/TAXID_3649/reference/tpm_mat.tsv"
mapratepath= "/home/pengken001/pipeline_output/TAXID_3649/reference/map_rate.tsv"
infopath="/home/pengken001/pipeline_output/TAXID_3649/reference/info.csv"
referencepath="/home/pengken001/pipeline_output/TAXID_3649/reference/cds_from_genomic.fna"
indexpath="/home/pengken001/pipeline_output/TAXID_3649/reference/cds_from_genomic.index"
kaldir="/home/pengken001/pipeline_output/TAXID_3649/reference/kallisto"
failedaccpath="/home/pengken001/pipeline_output/TAXID_3649/reference/failed.tsv"
retrylimit=2
workers=16
threads=2
fastqdir="/home/pengken001/pipeline_output/TAXID_3649/Step_2/fastq"
kmin=5
kmax=20
pseudoalignment_threshold=0

def ps_job(accession, index):
    fastqpath =os.path.join(fastqdir,accession+".fastq.gz")
    kaloutdir= os.path.join(kaldir, accession)
    result= misc.run_with_retries(retrylimit,
    read_map.launch_kallisto_quant,
    [threads, indexpath , kaloutdir , fastqpath],
    f"{accession}: Kallisto pseudoalignment failed. Retrying...",
    f"{accession}: Kallisto pseudoalignment of accession reads against reference...\n")
    if result == "failed":
        with open(failedaccpath, "a") as f:
            f.write(f"{accession}\tPS_failed\n")
        return f"{accession}: Aborted after {retrylimit} retries."
    map_rate = read_map.write_quant_info(accession, kaloutdir, tpm_matpath)
    with open(mapratepath, "a") as f:
        f.write(f"{accession}\t{map_rate}\n")
    print(f"{accession}: Pseudoalignment completed. {index}")
    return f"{accession}: processed."
    
def parallel_job(workers):
    with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as executor:
        results= [executor.submit(ps_job, accession, index) for index ,accession in enumerate(accessions)]
        for f in concurrent.futures.as_completed(results):
            print(f.result())

if __name__ == "__main__":
    accessions=  [file.split(".fastq")[0] for file in os.listdir(fastqdir) if ".gz" in file]
    print("kallisto index....")
    read_map.launch_kallisto_index(referencepath, indexpath)
    print("Parallel PS...")
    parallel_job(workers)
    maprate_df = pd.read_csv(mapratepath, sep="\t", header=None)
    maprate_dict= maprate_df.set_index(0).to_dict()[1]
    total, failed, passed, cutoff= classify.thresholder({key:val for key, val in maprate_dict.items()}, pseudoalignment_threshold)
    df= pd.DataFrame()
    df["total"]=[total]
    df["failed"]=[failed]
    df["passed"]=[passed]
    df["cutoff"]= [cutoff]
    df.to_csv(infopath)
    print(f"Reducing dimensions of TPM expression matrix ({len(passed)} accessions) using PCA-transformation...\n")
    Matrix= classify.mat_parser(tpm_matpath, passed)
    pca_data , pc_variances = classify.PCA_transformer(Matrix)
    df["pc_variances"]= [pc_variances]
    df["PC1"] = [pc_variances[0]]
    df.to_csv(infopath)
    print(f"PCA-transformation complete with {np.round(sum(pc_variances))}% of variance retained. (PC1= {pc_variances[0]}%)\n")
    
    print(f"Initiating k-means clustering of accesions based on PCA data.\nClustering iterations will walk from k={kmin} to k={kmax-1} to determine optimal number of clusters(k)...\n")
    k_cluster_assignment_dict, silhouette_coefficients = classify.kmeans_kwalk(pca_data, kmin, kmax)
    
    optimal_k, cluster_assignment , sc_max = classify.optimal_k_silhouette(kmin, kmax, silhouette_coefficients, k_cluster_assignment_dict)
    df["optimal_k"]= [optimal_k]
    df.to_csv(infopath)
    cluster_assignment_dict = {}
    for accession , cluster in zip(passed ,cluster_assignment):
        if cluster not in cluster_assignment_dict.keys():
            cluster_assignment_dict[int(cluster)] = [accession]
        else:
            cluster_assignment_dict[int(cluster)]+= [accession]
    silhouette_coefficients_dict = {int(k): sc for k , sc in zip(range(kmin,kmax +1), silhouette_coefficients)}
    df["silhouette_coefficients_dict"] = [silhouette_coefficients_dict]
    df["cluster_assignment_dict"] = [cluster_assignment_dict]
    df.to_csv(infopath)
    median_stat , mean_stat, min_stat , max_stat = classify.report_cluster_assignment_stats(cluster_assignment_dict)
    df["cluster_assignment_stats"] = [[optimal_k, sc_max, median_stat , mean_stat, min_stat , max_stat]]
    df.to_csv(infopath)