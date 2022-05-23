#setting sys.path for importing modules
import os
import sys
if __name__ == "__main__":
        abspath= os.getcwd()
        parent_module= os.path.join(abspath.split("LSTrAP-denovo")[0], "LSTrAP-denovo")
        sys.path.insert(0, parent_module)

import argparse
import concurrent.futures
from time import sleep
from datetime import datetime
from tqdm import tqdm
import numpy as np

from assembly import misc , report
from download import aspera
from setup import constants
from preprocess import read_map, classify

def download_PS_job(accession, index):
    '''Job to validate download path -> download via FTP/ascp -> Peudoalignment(PS) by Kallisto'''
    logfile.load()
    if type(logfile.contents["cluster"]["processed_acc"].get(accession)) == float:
        return f"{accession} already processed."
    #to un-sync workers
    if index < workers:
        sleep((index%workers)*5)
    if accession not in logfile.contents["cluster"]["processed_acc"].keys() or logfile.contents["cluster"]["processed_acc"].get(accession) == "Download failed": #check if accession has been downloaded/processed. Proceed with download if not.
        ascp_fullpath, ftp_fullpath, filesize = aspera.get_download_path(accession)
        fastqpath =os.path.join(C_fastqdir,accession+".fastq.gz")
        #download
        if download_method == "ascp":
            result= misc.run_with_retries(retrylimit,
            aspera.launch_ascp,
            [ascp_fullpath,fastqpath,min([filesizelimit,filesize])],
            f"{accession}: Download failed. Retrying...",
            f"{accession}: Downloading file via ascp...\n")
        elif download_method == "ftp":
            result= misc.run_with_retries(retrylimit,
            aspera.launch_curl,
            [ftp_fullpath,fastqpath,min([filesizelimit,filesize])],
            f"{accession}: Download failed. Retrying...",
            f"{accession}: Downloading file via ftp...\n")
        if result == "failed":
            logfile.load()
            logfile.contents["cluster"]["processed_acc"][accession]= "Download failed"
            logfile.update()
            return f"{accession}: Aborted after {retrylimit} retries."
        logfile.load()
        logfile.contents["cluster"]["processed_acc"][accession]= "Downloaded"
        logfile.update()
    
    if  logfile.contents["cluster"]["processed_acc"].get(accession)== "Downloaded" or logfile.contents["cluster"]["processed_acc"].get(accession) == "Kallisto failed":
        fastqpath =os.path.join(C_fastqdir,accession+".fastq.gz")
        kaloutdir= os.path.join(kaldir, accession)
        result= misc.run_with_retries(retrylimit,
        read_map.launch_kallisto_quant,
        [threads, indexpath , kaloutdir , fastqpath],
        f"{accession}: Kallisto pseudoalignment failed. Retrying...",
        f"{accession}: Kallisto pseudoalignment of accession reads against draft CDS of Preliminary assembly...\n")
        if result == "failed":
            logfile.load()
            logfile.contents["cluster"]["processed_acc"][accession]= "PS failed"
            logfile.update()
            return f"{accession}: Aborted after {retrylimit} retries."
        map_rate = read_map.write_quant_info(accession, kaloutdir, tpm_matpath)
        logfile.load()
        logfile.contents["cluster"]["processed_acc"][accession]= float(map_rate)
        logfile.update()
        print(f"{accession}: Pseudoalignment completed.")
        return f"{accession}: processed."
    
    
    

        
            
        
def parallel_job(workers):
    '''Wrapper to parallelize download and pseudolaignment jobs for each accession'''
    logfile.load()
    with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as executor:
        progress_bar= tqdm(total=len(accessions), desc="Accessions processed", unit="Acsn", leave=True)
        results= [executor.submit(download_PS_job, accession, index) for index, accession in enumerate(accessions)]
        for f in concurrent.futures.as_completed(results):
            if "processed" in f.result() or "Aborted" in f.result():
                progress_bar.update(1)
                progress_bar.set_postfix_str(s=f.result())
                print("\n")
            else:
                print(f.result())
        progress_bar.close()
        logfile.load()







if __name__ == "__main__":
    #retry limit determines the number of retries before aborting
    retrylimit= 2
    #arguments
    parser= argparse.ArgumentParser(description="LSTrAP-denovo.Cluster_accessions: Large-scale download , quality control and clustering of accessions based on transcriptome profiles.\n \
    NOTE: This is step 2 of 3 in the LSTrAP-denovo pipeline. Requires prior run of step 1: Preliminary_assembly.py. \n\
    Refer to https://github.com/pengkenlim/LSTrAP-denovo for more information on pipeline usage and implmentation")
    parser.add_argument("-o", "--output_dir", type=str, metavar= "", required=True,
    help= "Directory for data output. Directory needs to be same as for step 1 (Preliminary_assembly.py).")
    parser.add_argument("-ps", "--pseudoalignment_threshold", type=int ,metavar="", default=0 , choices=range(0, 70+1),
    help = "Specifies pseudoalignment threshold (%%PS) during quality control. Accessions that do not meet this threshold will be discarded and not be clustered. Threshold will be determined automatically based on kernel-density minima of %%PS by default.")
    parser.add_argument("-s","--filesizelimit" , type=int, metavar="", default=500 , choices=range(100, 1500),
    help="Specifies the size limit(mb) of accession read files to partially download. Limit set to 500 (mb) by default. Has a direct impact on the download time and pseudoalignment runtime in this step of the pipeline.\
    User is advised to reduce size limit when downloading and processing >500 accessions")
    parser.add_argument("-t", "--threads", type=int, metavar="", default=16, 
    help = "Total thread pool for workers. Needs to be divisible by number of workers.")
    parser.add_argument("-w", "--workers", type=int, metavar="", default=8, 
    help= "Specify the maximum workers for running multiple download-pseudoalignment jobs in parallel. Set to 8 by default.")
    parser.add_argument("-dm", "--download_method", type=str, metavar="", default="ftp", choices=["ascp","ftp"],
    help = "Method to download accession runs. ftp/ascp.")  
    parser.add_argument("-al", "--accessions_limit", type=int, metavar="", default=500,
    help= "Specifies the upper limit for number of accessions to download and process. Accessions will be selected from a pre-randomised list that was fetched during Preliminary_assembly.py run and stored in in the logs.json file.\
    Default set to 500.")
    parser.add_argument("-kr", "--k_range", type=str, metavar="", default="2:30",
    help = "Specifies the range of k (number of clusters) to iterate through during k-means clustering. Lower and upper limit seperated by colon(:). \
    Set to 2:30 by default. Increasing the upper limit for high-resolution clustering is advised when downloading and processing > 500 accessions or if accessions are expected to have various experimental permutations (genotypes, pertubations, tissue-types) ." )    
    parser.add_argument("-ct", "--consensus_threshold", type=int ,metavar="", default=0 , choices=range(0, 50+1),
    help = "Specifies consensus threshold of preliminary assembly. Default set to 0 where optimal threshold determined automatically in step 1 will be used.")
    parser.add_argument("-con", "--conti", action="store_true",
    help = "Resume incomplete run based on output directory. Only requires -o to run.")
    parser.add_argument("-f", "--force", action="store_true",
    help = "Delete data from previous Cluster_accessions.py run.")
    
    #parse args
    args=parser.parse_args()

    #assign output_dir and conti to variables
    outputdir= args.output_dir
    conti=args.conti
    force=args.force
    #check if outputdir and its subdir from previous step exists. If not, exit.
    P_fastqdir=os.path.join(outputdir,"prelim","fastq")
    P_ssadir=os.path.join(outputdir, "prelim","ssa")
    if not os.path.exists(P_fastqdir) or not os.path.exists(P_ssadir):
        sys.exit("Error in output directory. Either directory or its required sub-directory are not found. Exiting...")
    
    #check if logfile is present in output directory. exit if not there.
    if not os.path.exists(os.path.join(outputdir,"logs.json")):
        sys.exit("Error in output directory. logs.json file not found. Exiting...")
    #create logfile object, exit if step 1 not completed.
    logfile=misc.logfile(os.path.join(outputdir,"logs.json"))
    if  logfile.contents["prelim"]["status"] != "completed":
        sys.exit("Step 1 (Preliminary_assembly.py) was either incomplete or not started at all. \nNOTE: This is step 2 of 3 in the LSTrAP-denovo pipeline.\nExiting...")

    
    #assigning arguments to variables, writing to log OR fetching variables from log
    if conti==False:
        if logfile.contents["cluster"]["run_info"].get("init_time") != None:
            if force != True:
                temp_string = logfile.contents["cluster"]["run_info"].get("init_time")
                sys.exit(f"Previous incomplete run started on {temp_string} has been detected in logfile.\nPlease use --conti to continue incomplete run or --force to clear previous run data. Exiting...")
        if force == True:
            print("\n--force argument has been specified by user. Deleting data from previous run and clearing logs....\n")
            logfile.clear("cluster")
            if os.path.exists(os.path.join(outputdir, "cluster")):
                temp_string= os.path.join(outputdir, "cluster")
                os.system(f"rm -r {temp_string}")

        #assign argument vars
        pseudoalignment_threshold = args.pseudoalignment_threshold
        filesizelimit= args.filesizelimit * 1000000
        threadpool= args.threads
        workers=args.workers
        download_method= args.download_method
        accessions_limit= args.accessions_limit
        k_range= args.k_range
        consensus_threshold = args.consensus_threshold
        #check if threads/worker arguments are logical. Correct if necessary.
        if threadpool % workers != 0:
            print(f"Specified thread pool of {threadpool} is not divisible by number of workers.")
            threadpool= threadpool - (threadpool % workers)
            print(f"Using thread pool of {threadpool} instead.\n")
        threads=int(threadpool/workers)
        #check if consensus_threshold given is logical. Else, change it to 0 (auto).
        if consensus_threshold > logfile.contents["prelim"]["run_var"]["n_accessions"]:
            print("User-defined consensus threshold exceeded range. Run will continue with optimal consensus threshold automatically deternine in step 1.\n")
            consensus_threshold = 0
        
        #write run info (for final report) to logfile
        #if fetch optimal consensus threshold from log file if auto
        if consensus_threshold ==0:
            logfile.contents["cluster"]["run_info"]= {"Consensus_threshold_for_preliminary_assembly": logfile.contents["prelim"]["consensus"]["optimal"],
            "command_issued": " ".join(sys.argv),
            "init_time": datetime.now().strftime("%d/%m/%Y %H:%M:%S")}
        else:
            logfile.contents["cluster"]["run_info"]= {"Consensus_threshold_for_preliminary_assembly": consensus_threshold,
            "command_issued": " ".join(sys.argv),
            "init_time": datetime.now().strftime("%d/%m/%Y %H:%M:%S")}
        
        #write argument vars to logfile
        logfile.contents["cluster"]["run_var"]={
        "pseudoalignment_threshold": pseudoalignment_threshold,
        "filesizelimit": filesizelimit,
        "threadpool": threadpool,
        "workers" : workers,
        "download_method": download_method,
        "accessions_limit": accessions_limit,
        "k_range": k_range,
        "consensus_threshold" : consensus_threshold,
        }
        logfile.contents["cluster"]["processed_acc"]={}
        logfile.update()
        
         
                 
    elif conti==True:
        if logfile.contents["cluster"]["run_info"]["init_time"]==None:
            sys.exit(f"\nNo previous run initiation detected in {outputdir}. Exiting...")
        if logfile.contents["cluster"]["status"]== "completed":
            sys.exit(f"\nPrevious run initiated in {outputdir} has fully completed. There is nothing to run. Use --force to delete all previous run data in order to restart run.")
        #inherit run variables from previous run using logfile contents
        pseudoalignment_threshold, filesizelimit, threadpool, workers, download_method, accessions_limit, k_range, consensus_threshold = logfile.contents["cluster"]["run_var"].values()


        
    #make subdirs if needed
    C_fastqdir= os.path.join(outputdir, "cluster", "fastq") #to store downloaded fastq files
    if not os.path.exists(C_fastqdir):
        os.makedirs(C_fastqdir)   
    
    #get accessions and sci name from logfile
    accessions = logfile.contents["prelim"].get("total_acc")
    _, scientific_name, _, _, _ = logfile.contents["prelim"]["run_info"].values()
    taxid, _, _, _, _, _, _ , _, _, _ ,_, _ = logfile.contents["prelim"]["run_var"].values()
    
    print(f"Organism name: {scientific_name} (NCBI TaxID: {taxid})\n")
    print(f"Total accessions fetched from ENA: {len(accessions)}\n")
    
    print("Transfering accession .fastq files downloaded in step 1 to new directory...")
    for file in os.listdir(P_fastqdir): #xfer old files if not already done so
        if "fastq.gz" in file and not os.path.exists(os.path.join(C_fastqdir,file)): 
            os.system(f"cp {os.path.join(P_fastqdir , file)} {C_fastqdir}")
            if os.path.getsize(os.path.join(C_fastqdir, file)) > filesizelimit:
                os.system(f"truncate -s {filesizelimit}  {os.path.join(C_fastqdir, file)}")
                print(f"{file} transfered and truncated to file size limit of {filesizelimit/1000000} mb")
            else:
                print(f"{file} transfered.")
            logfile.contents["cluster"]["processed_acc"][file.split(".")[0]]="Downloaded"
    logfile.update()
    
    
    
    #directory to hold kallisto intermediate files
    kaldir= os.path.join(outputdir, "cluster", "kallisto")
    if not os.path.exists(kaldir):
        os.makedirs(kaldir)
    
    CT_int = logfile.contents["cluster"]["run_info"]["Consensus_threshold_for_preliminary_assembly"]
    
    fastapath=  logfile.contents["prelim"]["consensus"]["stats"].get(str(CT_int))[-1]
    indexpath= os.path.join(kaldir, fastapath.split("/")[-1] + ".index")
    result= misc.run_with_retries(retrylimit,
            read_map.launch_kallisto_index,
            [fastapath,indexpath],
            f"Kallisto index failed. Retrying...",
            f"\nContructing Kallisto index using Preliminary Assembly (CT={str(CT_int)}) generated from step 1...\n")
    if result == "failed":
        sys.exit(f"Kallisto index step failed after {retrylimit} retries. Exiting...")
    
    
    print(f"Initiating parallel download and pseudoalignment of {min(accessions_limit, len(accessions))} accessions...\n")
    #subseting accessions based on accession limit
    accessions = accessions[0:min(accessions_limit, len(accessions))]
    tpm_matpath= os.path.join(kaldir,"Draft_CDS_exp_mat.tsv")
    
    #initiate parallel DL and PS of accessions
    parallel_job(workers)
    
    print("\nParallel download and pseudoalignment complete.\n")
    
    logfile.load()
    #Quality control of accessions bassed on user-defined/ automatically determined PS threshold
    total, failed, passed, cutoff= classify.thresholder({key:val for key, val in logfile.contents["cluster"].get("processed_acc").items() if val is not str}, pseudoalignment_threshold)
    if pseudoalignment_threshold ==0:
        print(f"A total of {len(total)} accessions has been downloaded and pseudoaligned.\n{len(failed)} accessions failed QC based on auto-determined psedoalignment threshold of {cutoff}%\n")
    elif pseudoalignment_threshold > 0:
        print(f"A total of {len(total)} accessions has been downloaded and pseudoaligned.\n{len(failed)} accessions failed QC based on user-defined psedoalignment threshold of {cutoff}%\n")
    #write to log
    logfile.contents["cluster"]["qc"]["threshold"] = cutoff
    logfile.contents["cluster"]["qc"]["total"] = total
    logfile.contents["cluster"]["qc"]["passed"] = passed
    logfile.contents["cluster"]["qc"]["failed"] = failed
    logfile.update()
    
    print(f"Reducing dimensions of TPM expression matrix ({len(passed)} accessions) using PCA-transformation...\n")
    
    #read TPM expression matrix from file path and return subsetted matrix (without failed accessions)\
    Matrix= classify.mat_parser(tpm_matpath, passed)
    #normalise TPM values within each accession, PCA transfrom
    pca_data , pc_variances = classify.PCA_transformer(Matrix)    
    print(f"PCA-transformation complete with {np.round(sum(pc_variances))}% of variance retained. (PC1= {pc_variances[0]}%)\n")
    
    #extract k-means minimum and maximum values from range k_range variable parsed from arguments
    kmin = int(k_range.split(":")[0])
    kmax= int(k_range.split(":")[1])+1
    print(f"Initiating k-means clustering of accesions based on PCA data.\nClustering iterations will walk from k={kmin} to k={kmax} to determine optimal number of clusters(k)...\n")
    
    #k-means walk proper
    k_cluster_assignment_dict, silhouette_coefficients = classify.kmeans_kwalk(pca_data, kmin, kmax)
    #feed silhouette_coefficients and cluster assignments at different ks into function that determines optimal k 
    optimal_k, cluster_assignment , sc_max = classify.optimal_k_silhouette(kmin, kmax, silhouette_coefficients, k_cluster_assignment_dict)
    
    cluster_assignment_dict = {}
    for accession , cluster in zip(passed ,cluster_assignment):
        if cluster not in cluster_assignment_dict.keys():
            cluster_assignment_dict[int(cluster)] = [accession]
        else:
            cluster_assignment_dict[int(cluster)]+= [accession]
    silhouette_coefficients_dict = {int(k): sc for k , sc in zip(range(kmin,kmax +1), silhouette_coefficients)}
    print(cluster_assignment_dict)
    #write cluster assignments to log
    logfile.contents["cluster"]["kmeans"]["s_coeficient"]= silhouette_coefficients_dict
    logfile.contents["cluster"]["kmeans"]["cluster_assignment_dict"]= cluster_assignment_dict
    
    #get some stats on cluster assignments
    median_stat , mean_stat, min_stat , max_stat = classify.report_cluster_assignment_stats(cluster_assignment_dict)
    logfile.contents["cluster"]["kmeans"]["cluster_assignment_stats"]= [optimal_k, sc_max, median_stat , mean_stat, min_stat , max_stat]
    logfile.update()
    
    #report kmeans stats to user
    print(f"\nOptimal K-means iteration determined to be at k={optimal_k} with a silhouette coefficient of {sc_max}.\n\nAverage cluster size: {mean_stat} accessions\nMedian cluster size: {median_stat} accessions\nSize of largest cluster: {max_stat} accessions\nSize of smallest cluster: {min_stat} accessions\n")
    
    logfile.contents["cluster"]["status"]= "completed"
    logfile.update()
    print("LSTrAP-denovo.Cluster_accessions.py completed.\nGenerating html report...")
    report.generate_from_json_log(logfile.path, os.path.join(outputdir, "LSTrAP-denovo.html"), 2)
    
    
    
    
    
        
            