#setting sys.path for importing modules
import os
import sys
if __name__ == "__main__":
        abspath= os.getcwd()
        parent_module= os.path.join(abspath.split("HSS-Trans")[0], "HSS-Trans")
        sys.path.insert(0, parent_module)

from threadpoolctl import threadpool_limits
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
    #to slow down jobs
    sleep(1)
    logfile.load()
    #to unsync workers
    if index < workers:
        sleep((index%workers)*5)
    if type(logfile.contents["Step_2"]["processed_acc"].get(accession)) == float:
        return accession, index , "Already processed"
    if accession not in logfile.contents["Step_2"]["processed_acc"].keys() or logfile.contents["Step_2"]["processed_acc"].get(accession) == "Download failed" or logfile.contents["Step_2"]["processed_acc"].get(accession) == "Download link not found": #check if accession has been downloaded/processed. Proceed with download if not.
        ascp_fullpath, ftp_fullpath, filesize = aspera.get_download_path_ffq(accession)
        fastqpath =os.path.join(C_fastqdir,accession+".fastq.gz")
        if ascp_fullpath == "NOT_FOUND":
            return accession , index , "Download link not found"
            #logfile.load()
            #logfile.contents["Step_2"]["processed_acc"][accession]= "Download link not found"
            #logfile.update()
            #return f"{accession}: Aborted. Download link not found."
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
            return accession , index , "Download failed"
            #logfile.load()
            #logfile.contents["Step_2"]["processed_acc"][accession]= "Download failed"
            #logfile.update()
            #return f"{accession}: Aborted after {retrylimit} retries."
        #logfile.load()
        logfile.contents["Step_2"]["processed_acc"][accession]= "Downloaded"
        #logfile.update()
    
    if  logfile.contents["Step_2"]["processed_acc"].get(accession)== "Downloaded":
        fastqpath =os.path.join(C_fastqdir,accession+".fastq.gz")
        kaloutdir= os.path.join(kaldir, accession)
        result= misc.run_with_retries(retrylimit,
        read_map.launch_kallisto_quant,
        [threads, indexpath , kaloutdir , fastqpath],
        f"{accession}: Kallisto pseudoalignment failed. Retrying...",
        f"{accession}: Kallisto pseudoalignment of accession reads against draft CDSs...\n")
        if result == "failed":
            #os.system(f"rm {fastqpath}") add in final build
            return accession , index , "PS failed"
            #logfile.load()
            #logfile.contents["Step_2"]["processed_acc"][accession]= "PS failed"
            #logfile.update()
            #return f"{accession}: Aborted after {retrylimit} retries."
        map_rate = read_map.write_quant_info(accession, kaloutdir, tpm_matpath)
        #os.system(f"rm {fastqpath}") add in in final build
        print(f"{accession}: Pseudoalignment completed.")
        return accession , index , float(map_rate)
        #logfile.load()
        #logfile.contents["Step_2"]["processed_acc"][accession]= float(map_rate)
        #logfile.update()
        #return f"{accession}: processed."
    return accession , index , "Unknown exception"
    
    
    

        
            
        
def parallel_job(workers):
    '''Wrapper to parallelize download and pseudolaignment jobs for each accession'''
    logfile.load()
    with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as executor:
        progress_bar= tqdm(total=len(accessions), desc="Accessions processed", unit="Acsn", leave=True)
        results= [executor.submit(download_PS_job, accession, index) for index, accession in enumerate(accessions)]
        for f in concurrent.futures.as_completed(results):
            #print(f.result())
            accession , index , map_rate = f.result()
            if map_rate== "Already processed":
                msg = f"{accession} already processed."
            else:
                logfile.contents["Step_2"]["processed_acc"][accession]=map_rate 
                if map_rate == "Download link not found":
                    msg= f"{accession}: Aborted. Download link not found."
                elif map_rate == "Download failed":
                    msg= f"{accession}: Aborted. Download link not found."
                elif map_rate == "PS failed":
                    msg= f"{accession}: Aborted after {retrylimit} retries."
                elif type(map_rate) == float:
                    msg= f"{accession}: processed."
                else:
                    msg= f"{accession}: ERROR. Unknown exception occurred."
            progress_bar.update(1)
            progress_bar.set_postfix_str(s=msg)
            logfile.update()
            print("\n")
        progress_bar.close()
        logfile.load()
        

def download_job(link, index):
    '''Job to download accessions'''
    #to slow down jobs
    sleep(1)
    logfile.load()
    filename= link.split("/")[-1]
    #to unsync workers
    if index < workers:
        sleep((index%workers)*5)
    if filename in logfile.contents["Step_2"]["selected_accessions"]["download_progress"].keys():
        if logfile.contents["Step_2"]["selected_accessions"]["download_progress"].get(filename) == "Downloaded":
            return filename , index, "Already downloaded"
            #return f"Skipping {filename} as it had already been downloaded."
    fastqpath= os.path.join(F_fastqdir, filename)
    if download_method == "ascp":
        ascp_fullpath = link.replace("ftp://ftp.sra.ebi.ac.uk/", "era-fasp@fasp.sra.ebi.ac.uk:")
        result= misc.run_with_retries(retrylimit,
        aspera.launch_ascp,
        [ascp_fullpath,fastqpath,0],#0 = no limit, full download
        f"{filename}: Download failed. Retrying...",
        f"{filename}: Downloading file via ascp...\n")
    elif download_method == "ftp":
         result= misc.run_with_retries(retrylimit,
            aspera.launch_curl,
            [link,fastqpath,0], #0 = no limit, full download
            f"{filename}: Download failed. Retrying...",
            f"{filename}: Downloading file via ftp...\n")
    if result == "failed":
        #logfile.load()
        #logfile.contents["Step_2"]["selected_accessions"]["download_progress"][filename]= "Download failed"
        #logfile.update()
        return filename , index, "Download failed"
        #return f"{filename}: Download aborted after {retrylimit} retries."
    else:
        #unzip and rename headers (because trinity sometimes have issues with fastq directly downloaded from ENA)
        unzippedfastqpath= fastqpath.split(".gz")[0]
        if "_1" in unzippedfastqpath:
            os.system(f"zcat {fastqpath} | awk \'{{print (NR%4 == 1) ? \"@1_\" ++i \"/1\": $0}}\' > {unzippedfastqpath}")
        elif "_2" in unzippedfastqpath:
            os.system(f"zcat {fastqpath} | awk \'{{print (NR%4 == 1) ? \"@1_\" ++i \"/1\": $0}}\' > {unzippedfastqpath}")
        os.system(f"rm {fastqpath}")
        #logfile.load()
        #logfile.contents["Step_2"]["selected_accessions"]["download_progress"][filename]= "Downloaded"
        #logfile.update()
        print(f"{filename}: Downloaded.")
        return filename , index, "Downloaded"
        #return f"{filename}: Downloaded."
    return filename , index, "Unknown exception"
    
    
def parallel_download(workers):
    '''Wrapper to parallelize download each file'''
    logfile.load()
    with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as executor:
        progress_bar= tqdm(total=len(FTP_links), desc="Fastq Downloaded", unit="fq", leave=True)
        results= [executor.submit(download_job, link, index) for index, link in enumerate(FTP_links)]
        for f in concurrent.futures.as_completed(results):
            
            filename , index, returncode = f.result()
            if returncode == "Already downloaded":
                msg= f"{filename}: Already been downloaded."
            else:
                logfile.contents["Step_2"]["selected_accessions"]["download_progress"][filename]= returncode
                logfile.update()
                if returncode== "Download failed":
                   msg=  f"{filename}: Download aborted after {retrylimit} retries."
                elif returncode == "Downloaded":
                   msg= f"{filename}: Downloaded."
            progress_bar.update(1)
            progress_bar.set_postfix_str(s=msg)
            print("\n")
        progress_bar.close()
        logfile.load()


if __name__ == "__main__":
    #retry limit determines the number of retries before aborting
    retrylimit= 2
    #arguments
    parser= argparse.ArgumentParser(description="HSS-Trans.SelectAccessions.py: Selection of representative accessions for transcriptome assembly.\n \
    NOTE: This is step 2 of 2 in the HSS-Trans pipeline. Requires prior run of step 1: MakeDraftCDS.py. \n\
    Refer to https://github.com/pengkenlim/HSS-Trans for more information on pipeline usage and implmentation")
    parser.add_argument("-o", "--output_dir", type=str, metavar= "", required=True,
    help= "Directory for data output. Directory needs to be same as for step 1 (MakeDraftCDS.py).")
    parser.add_argument("-ps", "--pseudoalignment_threshold", type=int ,metavar="", default=0 , choices=range(0, 70+1),
    help = "Specifies reads pseudoaligned (%%PS) threshold for quality control. Expression data of accessions that do not meet this threshold will be excluded as features for clustering.\
    Set to 0 by default where %%PS threshold will be set to be the lower bound of the %%PS distribution (Q1 - 1.5 *IQR) or 20%%, whichever is higher.")
    parser.add_argument("-s","--filesizelimit" , type=int, metavar="", default=500 , choices=range(100, 1000),
    help="Specifies the size limit(mb) for the partial download of accession read files (gzip compressed) for pseudoalignment. Limit set to 500 (mb) by default. \
    Size maybe decreased to improve the overall runtime (Download and Psudoalignment) and storage used in this step of the pipeline.\
    However, doing so might compromise accurate gene expression quantification as a result of limited sequencing depth.")
    parser.add_argument("-t", "--threads", type=int, metavar="", default=8, 
    help = "Total thread pool for workers. Needs to be divisible by number of workers.\n")
    parser.add_argument("-w", "--workers", type=int, metavar="", default=4, 
    help= "Specify the maximum workers for running multiple download-pseudoalignment jobs in parallel. Set to 4 by default.")
    parser.add_argument("-dm", "--download_method", type=str, metavar="", default="ascp", choices=["ascp","ftp"],
    help = "Method to download accession runs. ftp/ascp. Set to aspera download (ascp) by default.")  
    parser.add_argument("-al", "--accessions_limit", type=int, metavar="", default=500,
    help= "Specifies the upper limit for number of accessions to download and process. Accessions will be selected from a pre-randomised list that was fetched during MakeDraftCDS.py run and stored in in the logs.json file.\
    Default set to 500.")
    parser.add_argument("-kr", "--k_range", type=str, metavar="", default="auto",
    help = "Specifies the range of k (number of clusters) to iterate through during clustering. Lower and upper limit seperated by colon(:). \
    Set to auto(5:20) by default. Optimal k will be chosen from this range based on silhouette coefficient, a clustering performance metric.\
    As such, please set range within expectations based on heterogenity of expression data for that organism.")    
    parser.add_argument("-ct", "--consensus_threshold", type=int ,metavar="", default=0 , choices=range(0, 50+1),
    help = "Specifies consensus threshold of preliminary assembly. Default set to 0 where optimal threshold determined automatically in step 1 will be used.")
    parser.add_argument("-clib", "--cluster_lib_size", type=int ,metavar="", default=2000 , choices=range(500, 10*1000),
    help = "Specifies the minimum library size (mb) for each cluster to guide the selection of representative accessions to download. Total file sizes of read files (gzipped) to download for each cluster is used as an approximation of libary size instead of number of reads.")
    parser.add_argument("-con", "--conti", action="store_true",
    help = "Resume incomplete run based on output directory. Only requires -o to run.")
    parser.add_argument("-f", "--force", action="store_true",
    help = "Delete data from previous SelectAccessions.py run and start a fresh SelectAccessions.py in output directory.")
    
    #parse args
    args=parser.parse_args()

    #assign output_dir and conti to variables
    outputdir= args.output_dir
    conti=args.conti
    force=args.force
    #check if outputdir and its subdir from previous step exists. If not, exit.
    P_fastqdir=os.path.join(outputdir,"Step_1","fastq")
    P_ssadir=os.path.join(outputdir, "Step_1","ssa")
    if not os.path.exists(P_fastqdir) or not os.path.exists(P_ssadir):
        sys.exit("Error in output directory. Either directory or its required sub-directory are not found. Exiting...")
    
    #check if logfile is present in output directory. exit if not there.
    if not os.path.exists(os.path.join(outputdir,"logs.json")):
        sys.exit("Error in output directory. logs.json file not found. Exiting...")
    #create logfile object, exit if step 1 not completed.
    logfile=misc.logfile(os.path.join(outputdir,"logs.json"))
    if  logfile.contents["Step_1"]["status"] != "completed":
        sys.exit("Step 1 (MakeDraftCDS.py) is either incomplete or not started at all. \nNOTE: This is step 2 of 3 in the HSS-Trans pipeline.\nExiting...")

    
    #assigning arguments to variables, writing to log OR fetching variables from log
    if conti==False:
        if logfile.contents["Step_2"]["run_info"].get("init_time") != None:
            if force != True:
                temp_string = logfile.contents["Step_2"]["run_info"].get("init_time")
                sys.exit(f"Previous incomplete run started on {temp_string} has been detected in logfile.\nPlease use --conti to continue incomplete run or --force to clear previous run data. Exiting...")
        if force == True:
            print("\n--force argument has been specified by user. Deleting data from previous run and clearing logs....\n")
            logfile.clear("Step_2")
            if os.path.exists(os.path.join(outputdir, "Step_2")):
                temp_string= os.path.join(outputdir, "Step_2")
                os.system(f"rm -r {temp_string}")

        #assign argument vars
        pseudoalignment_threshold = args.pseudoalignment_threshold
        filesizelimit= args.filesizelimit * 1048576
        threadpool= args.threads
        workers=args.workers
        download_method= args.download_method
        accessions_limit= args.accessions_limit
        k_range= args.k_range
        consensus_threshold = args.consensus_threshold
        cluster_lib_size = args.cluster_lib_size
        #check if threads/worker arguments are logical. Correct if necessary.
        if threadpool % workers != 0:
            print(f"Specified thread pool of {threadpool} is not divisible by number of workers.")
            threadpool= threadpool - (threadpool % workers)
            print(f"Using thread pool of {threadpool} instead.\n")
        threads=int(threadpool/workers)
        #check if consensus_threshold given is logical. Else, change it to 0 (auto).
        if consensus_threshold > logfile.contents["Step_1"]["run_var"]["n_accessions"]:
            print("User-defined consensus threshold exceeded range. Run will continue with optimal consensus threshold automatically deternine in step 1.\n")
            consensus_threshold = 0
        
        #write run info (for final report) to logfile
        #if fetch optimal consensus threshold from log file if auto
        if consensus_threshold ==0:
            logfile.contents["Step_2"]["run_info"]= {"Consensus_threshold_for_preliminary_assembly": logfile.contents["Step_1"]["consensus"]["optimal"],
            "command_issued": " ".join(sys.argv),
            "init_time": datetime.now().strftime("%d/%m/%Y %H:%M:%S")}
        else:
            logfile.contents["Step_2"]["run_info"]= {"Consensus_threshold_for_preliminary_assembly": consensus_threshold,
            "command_issued": " ".join(sys.argv),
            "init_time": datetime.now().strftime("%d/%m/%Y %H:%M:%S")}
        
        #write argument vars to logfile
        logfile.contents["Step_2"]["run_var"]={
        "pseudoalignment_threshold": pseudoalignment_threshold,
        "filesizelimit": filesizelimit,
        "threadpool": threadpool,
        "workers" : workers,
        "download_method": download_method,
        "accessions_limit": accessions_limit,
        "k_range": k_range,
        "consensus_threshold" : consensus_threshold,
        "cluster_lib_size" : cluster_lib_size}
        logfile.contents["Step_2"]["processed_acc"]={}
        logfile.update()
        
         
                 
    elif conti==True:
        if logfile.contents["Step_2"]["run_info"]["init_time"]==None:
            sys.exit(f"\nNo previous run initiation detected in {outputdir}. Exiting...")
        if logfile.contents["Step_2"]["status"]== "completed":
            sys.exit(f"\nPrevious run initiated in {outputdir} has fully completed. There is nothing to run. Use --force to delete all previous run data in order to restart run.")
        #inherit run variables from previous run using logfile contents
        pseudoalignment_threshold, filesizelimit, threadpool, workers, download_method, accessions_limit, k_range, consensus_threshold , cluster_lib_size= logfile.contents["Step_2"]["run_var"].values()
        threads=int(threadpool/workers)
        print("\n--conti argument has been specified by user. Inheriting arguments from previous run ....\n")


        
    #make subdirs if needed
    C_fastqdir= os.path.join(outputdir, "Step_2", "fastq") #to store downloaded fastq files
    if not os.path.exists(C_fastqdir):
        os.makedirs(C_fastqdir)   
    
    #get accessions and sci name from logfile
    accessions = logfile.contents["Step_1"].get("total_acc")
    _, scientific_name, _, _, _ = logfile.contents["Step_1"]["run_info"].values()
    taxid, _, _, _, _, _, _ , _, _, _ ,_, _ = logfile.contents["Step_1"]["run_var"].values()
    
    print(f"Organism name: {scientific_name} (NCBI TaxID: {taxid})\n")
    if len(accessions) == 2000:
        print(f"Total accessions fetched from ENA: {len(accessions)} (capped)\n")
    else:    
        print(f"Total accessions fetched from ENA: {len(accessions)}\n")
    
    print("Transfering accession .fastq files downloaded in step 1 to new directory...")
    for file in os.listdir(P_fastqdir): #xfer old files if not already done so
        if "fastq.gz" in file and not os.path.exists(os.path.join(C_fastqdir,file)): 
            os.system(f"cp {os.path.join(P_fastqdir , file)} {C_fastqdir}")
            if os.path.getsize(os.path.join(C_fastqdir, file)) > filesizelimit:
                os.system(f"truncate -s {filesizelimit}  {os.path.join(C_fastqdir, file)}")
                print(f"{file} transfered and truncated to file size limit of {filesizelimit/1048576} mb")
            else:
                print(f"{file} transfered.")
            logfile.contents["Step_2"]["processed_acc"][file.split(".")[0]]="Downloaded"
    logfile.update()
    
    
    
    #directory to hold kallisto intermediate files
    kaldir= os.path.join(outputdir, "Step_2", "kallisto")
    if not os.path.exists(kaldir):
        os.makedirs(kaldir)
    
    CT_int = logfile.contents["Step_2"]["run_info"]["Consensus_threshold_for_preliminary_assembly"]
    
    fastapath=  logfile.contents["Step_1"]["consensus"]["stats"].get(str(CT_int))[-1]
    indexpath= os.path.join(kaldir, fastapath.split("/")[-1] + ".index")
    if not os.path.exists(indexpath):
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
    total, failed, passed, cutoff= classify.thresholder({key:val for key, val in logfile.contents["Step_2"].get("processed_acc").items() if type(val) is not str}, pseudoalignment_threshold)
    if pseudoalignment_threshold ==0:
        print(f"A total of {len(total)} accessions has been downloaded and pseudoaligned.\n{len(failed)} accessions failed QC based on auto-determined psedoalignment threshold of {cutoff}%\n")
    elif pseudoalignment_threshold > 0:
        print(f"A total of {len(total)} accessions has been downloaded and pseudoaligned.\n{len(failed)} accessions failed QC based on user-defined psedoalignment threshold of {cutoff}%\n")
    #write to log
    logfile.contents["Step_2"]["qc"]["threshold"] = cutoff
    logfile.contents["Step_2"]["qc"]["total"] = total
    logfile.contents["Step_2"]["qc"]["passed"] = passed
    logfile.contents["Step_2"]["qc"]["failed"] = failed
    logfile.update()
    
    print(f"Reducing dimensions of TPM expression matrix ({len(passed)} accessions) using PCA-transformation...\n")
    
    #read TPM expression matrix from file path and return subsetted matrix (without failed accessions)\
    Matrix= classify.mat_parser(tpm_matpath, passed)
    #normalise TPM values within each accession, PCA transfrom. Context manager to limit core usage
    with threadpool_limits(user_api="openmp", limits=threadpool):
        pca_data , pc_variances = classify.PCA_transformer(Matrix)    
    print(f"PCA-transformation complete with {np.round(sum(pc_variances))}% of variance retained. (PC1= {pc_variances[0]}%)\n")
    
    #extract k-means minimum and maximum values from range k_range variable parsed from arguments
    if k_range == "auto":
        kmin = 5
        kmax= 20
        if kmax > len(passed):
            print(f"Number of accessions that passed QC is below upper limit of k range {kmax}.\nRange changed to {kmin}:{len(passed)}.")
            kmax = len(passed)
    else:
        kmin = int(k_range.split(":")[0])
        kmax= int(k_range.split(":")[1])+1 #non-inclusive
        if kmax > len(passed):
            print(f"Number of accessions that passed QC is below upper limit of k range {kmax}.\nRange changed to {kmin}:{len(passed)}.")
            kmax = len(passed)
    if "cluster_assignment_stats" not in logfile.contents["Step_2"]["kmeans"].keys():
        print(f"Initiating k-means clustering of accesions based on PCA data.\nClustering iterations will walk from k={kmin} to k={kmax-1} to determine optimal number of clusters(k)...\n")
        #k-means walk proper under context manager (threadpool_limits) to limit core usage. kmeans package uses all available cores by default.
        with threadpool_limits(user_api="openmp", limits=threadpool):
            k_cluster_assignment_dict, silhouette_coefficients = classify.kmeans_kwalk(pca_data, kmin, kmax)
        #feed silhouette_coefficients and cluster assignments at different ks into function that determines optimal k 
        optimal_k, cluster_assignment , sc_max = classify.optimal_k_silhouette(kmin, kmax, silhouette_coefficients, k_cluster_assignment_dict)
        
        cluster_assignment_dict = {}
        for accession , cluster in zip(passed ,cluster_assignment):
            if str(cluster) not in cluster_assignment_dict.keys():
                cluster_assignment_dict[str(cluster)] = [accession]
            else:
                cluster_assignment_dict[str(cluster)]+= [accession]
        silhouette_coefficients_dict = {int(k): sc for k , sc in zip(range(kmin,kmax +1), silhouette_coefficients)}
        #write cluster assignments to log
        logfile.contents["Step_2"]["kmeans"]["s_coeficient"]= silhouette_coefficients_dict
        logfile.contents["Step_2"]["kmeans"]["cluster_assignment_dict"]= cluster_assignment_dict
        #get some stats on cluster assignments
        median_stat , mean_stat, min_stat , max_stat = classify.report_cluster_assignment_stats(cluster_assignment_dict)
        logfile.contents["Step_2"]["kmeans"]["cluster_assignment_stats"]= [optimal_k, sc_max, median_stat , mean_stat, min_stat , max_stat]
        logfile.update()
    else:
        optimal_k, sc_max, median_stat , mean_stat, min_stat , max_stat = logfile.contents["Step_2"]["kmeans"]["cluster_assignment_stats"]
        cluster_assignment_dict = logfile.contents["Step_2"]["kmeans"]["cluster_assignment_dict"]
        silhouette_coefficients_dict = logfile.contents["Step_2"]["kmeans"]["s_coeficient"]
    print(cluster_assignment_dict)
    #report kmeans stats to user
    print(f"\nOptimal K-means iteration determined to be at k={optimal_k} with a silhouette coefficient of {sc_max}.\n\nAverage cluster size: {mean_stat} accessions\nMedian cluster size: {median_stat} accessions\nSize of largest cluster: {max_stat} accessions\nSize of smallest cluster: {min_stat} accessions\n")
    print(f"Selecting representative accessions from each cluster based on target library size....\nNote: Fetching metadata might take some time.")
    logfile.load()
    clusters = [ int(k) for k in list(cluster_assignment_dict.keys())]
    clusters.sort()
    clusters= [str(k) for k in clusters]
    processed_acc_dict = logfile.contents["Step_2"].get("processed_acc")
    if "FTP_links" not in logfile.contents["Step_2"]["selected_accessions"].keys():
        logfile.contents["Step_2"]["selected_accessions"]["FTP_links"]=[]
        logfile.update()
    logfile.load()
    for cluster in clusters:
        #check if selection is complete for cluster
        if cluster in logfile.contents["Step_2"]["selected_accessions"].keys():
            print(f"Cluster{cluster}: cluster representatives selected.\n")
        else:    
            accessions= cluster_assignment_dict.get(cluster)

            #this dict shall accessions as key and PS% stats as values
            accession_dict = {acc: processed_acc_dict.get(acc) for acc in accessions}

            ranked_ps = list(accession_dict.values())
            ranked_ps.sort()

            accession_dict_sorted={}
            #brute force using nested loop to sort the dictionary based on %PS values
            for ps in ranked_ps:
                for key, val in accession_dict.items():
                    if ps == val:
                        accession_dict_sorted[key]= val
            logfile.contents["Step_2"]["selected_accessions"][cluster]={}
            libsize=0
            for accession in accession_dict_sorted.keys():
                if libsize < (cluster_lib_size * 1048*1048): #minimum library size of 2gb (forward + reverse, gunzipped) per cluster
                    forward, reverse = "ERR", "ERR"
                    ftp_links= aspera.launch_ffq_ftp(accession)
                        #check if paired end
                    if len(ftp_links) > 1:
                        for file in ftp_links:
                            if file["filenumber"] == 1:
                                forward = file
                            elif file["filenumber"] == 2:
                                reverse = file
                        #add read files to total libary size
                        if forward != "ERR" and  reverse != "ERR":
                            libsize+= forward.get("filesize") + reverse.get("filesize")
                            logfile.contents["Step_2"]["selected_accessions"][cluster][accession]={"Library_size": forward.get("filesize") + reverse.get("filesize")}
                            logfile.contents["Step_2"]["selected_accessions"]["FTP_links"] += [forward.get("url") ,reverse.get("url")]
                            logfile.update()
            print(f"Cluster {cluster}: cluster representatives selected.\n")
    
    FTP_links = logfile.contents["Step_2"]["selected_accessions"]["FTP_links"]
    print(f"A total of {int(len(FTP_links)/2)} accessions have been selected.\nInitiating download...\n")
    #prepare to track download progress
    if "download_progress" not in logfile.contents["Step_2"]["selected_accessions"].keys():
        logfile.contents["Step_2"]["selected_accessions"]["download_progress"] = {}
        logfile.update()
    
    F_fastqdir=os.path.join(outputdir,"Step_2","selected_accessions")
    if not os.path.exists(F_fastqdir):
        os.makedirs(F_fastqdir)
     
    parallel_download(workers)
    
    #make samples_file for trinity
    with open(os.path.join(F_fastqdir, "Samples_for_trinity.tsv"), "w") as f:
        for cluster in clusters:
            condition= f"Cluster_{cluster}"
            for index , accession in enumerate(logfile.contents["Step_2"]["selected_accessions"][cluster]):
                replicate=  f"Cluster_{cluster}_rep{index}"
                forwardpath , reversepath = os.path.join(F_fastqdir,accession+"_1.fastq"), os.path.join(F_fastqdir,accession+"_2.fastq")
                f.write(f"{condition}\t{replicate}\t{forwardpath}\t{reversepath}")
                
                
    print("\nTrinity sample file have been created at "+ os.path.join(F_fastqdir, "Samples_for_trinity.tsv.\n"))
    print("Install Trinity (https://github.com/trinityrnaseq/), edit and run Run_Trinity.py to start de novo / genome-guided transcriptome assembly.\n")
        
    
    logfile.contents["Step_2"]["status"]= "completed"
    logfile.update()
    print("HSS-Trans.SelectAccessions.py completed.\nGenerating html report...\n")
    report.generate_from_json_log(logfile.path, os.path.join(outputdir, "HSS-Trans.html"), 2)
    
    
    
    
    
        
            