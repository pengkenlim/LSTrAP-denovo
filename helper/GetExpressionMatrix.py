#setting sys.path for importing modules
import os
import sys
if __name__ == "__main__":
        abspath= os.getcwd()
        parent_module= os.path.join(abspath.split("LSTrAP-denovo")[0], "LSTrAP-denovo")
        sys.path.insert(0, parent_module)

from threadpoolctl import threadpool_limits
import argparse
import concurrent.futures
from time import sleep
from datetime import datetime
from tqdm import tqdm
import numpy as np
import func_timeout
import random


from assembly import misc , report
from download import aspera, ena
from setup import constants
from preprocess import read_map, classify


def download_PS_job(accession, index):
    '''Job to validate download path -> download via FTP/ascp -> Peudoalignment(PS) by Kallisto'''
    try:
    #to slow down jobs
        sleep(0.2)
        logfile.load()
        #to unsync workers
        if index < workers:
            sleep((index%workers)*5)
        if type(logfile.contents["processed_acc"].get(accession)) == float:
            return accession, index , "Already processed"
        if logfile.contents["processed_acc"].get(accession) == "PS failed":
            return accession, index , "PS failed"
        if accession not in logfile.contents["processed_acc"].keys() or logfile.contents["processed_acc"].get(accession) == "Download failed" or logfile.contents["processed_acc"].get(accession) == "Download failed because link not found" or logfile.contents["processed_acc"].get(accession) == "Unknown exception" : #check if accession has been downloaded/processed. Proceed with download if not.
            ascp_fullpath, ftp_fullpath, filesize = aspera.get_download_path_ffq(accession)
            fastqpath =os.path.join(fastqdir,accession+".fastq.gz")
            #remove fastq if exists
            if os.path.exists(fastqpath):
                os.system(f"rm {fastqpath}")
            if ascp_fullpath == "NOT_FOUND":
                return accession , index , "Download failed because link not found"
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
            logfile.contents["processed_acc"][accession]= "Downloaded"

        if  logfile.contents["processed_acc"].get(accession)== "Downloaded" and os.path.exists(fastqpath):
            if os.path.getsize(fastqpath)==0:
                os.system(f"rm {fastqpath}")
                return accession , index , "Download failed"
            fastqpath =os.path.join(fastqdir,accession+".fastq.gz")
            kaloutdir= os.path.join(kaldir, accession)
            result= misc.run_with_retries(retrylimit,
            read_map.launch_kallisto_quant,
            [threads, indexpath , kaloutdir , fastqpath],
            f"{accession}: Kallisto pseudoalignment failed. Retrying...",
            f"{accession}: Kallisto pseudoalignment of accession reads against draft CDSs...\n")
            if result == "failed" or not os.path.exists(kaloutdir):
                os.system(f"rm {fastqpath}") #add in final build
                return accession , index , "PS failed"
            map_rate = read_map.write_quant_info(accession, kaloutdir, tpm_matpath)
            os.system(f"rm {fastqpath}") #add in in final build
            os.system(f"rm -r {kaloutdir}") #add in in final build
            print(f"{accession}: Pseudoalignment completed.")
            return accession , index , float(map_rate)
        return accession , index , "Unknown exception"
    except:
        return accession , index , "Unknown exception"
        
def runjob(f, accession , index ,max_wait ):
    '''Timeout wrapper for task'''
    try:
        return func_timeout.func_timeout(max_wait, f, args=(accession,index))
    except func_timeout.FunctionTimedOut:
        pass
    return accession , index, "Unknown exception"

def parallel_job(workers):
    '''Wrapper to parallelize download and pseudolaignment jobs for each accession'''
    with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as executor:
        progress_bar= tqdm(total=len(accessions), desc="Accessions processed", unit="Acsn", leave=True)
        results= [executor.submit(runjob,download_PS_job, accession, index, 1200) for index, accession in enumerate(accessions)]
        for f in concurrent.futures.as_completed(results): # 20 minutes
            accession , index , map_rate = f.result()
            if map_rate== "Already processed":
                msg = f"{accession} already processed."
            else:
                with open(pathtoprocessed, "a") as f:
                    f.write(f"{accession}\t{map_rate}\n")
                if map_rate == "Download failed because link not found":
                    msg= f"{accession}: Aborted. Download link not found."
                elif map_rate == "Download failed":
                    msg= f"{accession}: Aborted. Download failed."
                elif map_rate == "PS failed":
                    msg= f"{accession}: Aborted after {retrylimit} retries."
                elif type(map_rate) == float:
                    msg= f"{accession}: processed."
                else:
                    msg= f"{accession}: ERROR. Unknown exception occurred."
            progress_bar.update(1)
            progress_bar.set_postfix_str(s=msg)
            
            print("\n")
        progress_bar.close()

if __name__ == "__main__":
    #ARGPARSE
    parser= argparse.ArgumentParser(description="LSTrAP-denovo.GetExpressionMatrix.py: Helper script to generate gene expression matrix via downloading and pseudolaigning of RNA-seq accessions from ENA en masse to assembled CDS.\n\
    Refer to https://github.com/pengkenlim/LSTrAP-denovo for more information on pipeline usage and implementation .\n")
    
    parser.add_argument("-o", "--output_dir", type=str, metavar= "", required=True,
    help= "Directory for data output. Directory needs to be same as for Step 1 and Step 2 of the LSTrAP-denovo pipeline.")
    
    parser.add_argument("-s","--filesizelimit" , type=int, metavar="", default=500,
    help="Specifies the size limit(mb) for the partial download of accession read files (gzip compressed) for pseudoalignment. Limit set to 500 (mb) by default. \
    Size can be decreased to improve the overall runtime (Download and Psudoalignment) and storage used in this step of the pipeline.\
    However, doing so might compromise accurate gene expression quantification as a result of limited sequencing depth.")
    
    parser.add_argument("-t", "--threads", type=int, metavar="", default=8, 
    help = "Total thread pool for workers. Needs to be divisible by number of workers.\n")
    
    parser.add_argument("-w", "--workers", type=int, metavar="", default=4, 
    help= "Specify the maximum workers for running multiple download-pseudoalignment jobs in parallel. Set to 4 by default.")
    
    parser.add_argument("-dm", "--download_method", type=str, metavar="", default="ascp", choices=["ascp","ftp"],
    help = "Method to download accession runs. ftp/ascp. Set to aspera download (ascp) by default.")
    
    parser.add_argument("-al", "--accessions_limit", type=int, metavar="", default=500,
    help= "Specifies the upper limit for number of accessions to download and process. Accessions will be selected from a pre-randomised list that was fetched fetched from ENA stored in in the logs.json file.\
    Default set to 500.")
    
    
    parser.add_argument("-f", "--force", action="store_true",
    help = "Delete data from previous GetExpressionMatrix.py run.")
    
    ME_group_1 = parser.add_mutually_exclusive_group(required=True)
    ME_group_1.add_argument("-i", "--id", type=int, metavar="", 
    help= "NCBI TaxID of organism for fetching SRA run accessions.")
    ME_group_1.add_argument("-con", "--conti", action="store_true",
    help = "Resume incomplete run based on output directory. Only requires -o to run.")
    
    #banner
    misc.print_logo("GetExpressionMatrix.py")
    
    #parse args
    args=parser.parse_args()
    
    #assign output_dir and conti to variables
    retrylimit= 0
    outputdir= args.output_dir
    conti=args.conti
    force=args.force
    kaldir= os.path.join(outputdir, "GetExpressionMatrix", "kallisto")
    fastqdir = os.path.join(outputdir,"GetExpressionMatrix","fastq")
    fastapath= os.path.join(outputdir, "Annotations", "cds_from_transcripts.fasta")
    indexpath= os.path.join(kaldir, fastapath.split("/")[-1] + ".index")
    logfile=misc.logfile_expmat(os.path.join(outputdir,"logs_expmat.json"))
    tpm_matpath= os.path.join(kaldir,"exp_mat_untransposed.tsv")
    tpm_matpath_T= os.path.join(outputdir, "GetExpressionMatrix" ,"exp_mat.tsv")
    

    taxid= args.id
    filesizelimit= args.filesizelimit * 1048576
    threadpool= args.threads
    workers=args.workers
    download_method= args.download_method
    accessions_limit= args.accessions_limit
    
    #check requirements here.
    if not os.path.exists(fastapath):
        sys.exit(f"\n{fastapath} not found. Exiting...")
    #######################
    if conti == False:
        if logfile.contents["run_info"].get("init_time") != None:
            if force != True:
                temp_string = logfile.contents["run_info"].get("init_time")
                sys.exit(f"Previous incomplete run started on {temp_string} has been detected in logfile.\nPlease use --conti to continue incomplete run or --force to clear previous run data. Exiting...")
        if force == True:
            print("\n--force argument has been specified by user. Deleting data from previous run and clearing logs....\n")
            logfile.clear()
            temp_string= os.path.join(outputdir, "GetExpressionMatrix")
            os.system(f"rm -r {temp_string}")
        
        scientific_name= ena.get_sciname(taxid)
        if type(scientific_name) is  not list:
            sys.exit(f"TaxID {taxid} is invalid/not found. Exiting...")
        elif len(scientific_name) > 1:
            sys.exit(f"More than one organism found for TaxID {taxid}. Exiting...")
        scientific_name= scientific_name[0]["scientific_name"]
        print(f"\nFetching RNA-seq accessions of {scientific_name} (NCBI TaxID: {taxid}) from ENA..\n")
        accessions = ena.get_runs(taxid)
        random.shuffle(accessions)
        print(f"Total accessions fetched from ENA: {len(accessions)}\n")
        
        if threadpool % workers != 0:
            print(f"Specified thread pool of {threadpool} is not divisible by number of workers.")
            threadpool= threadpool - (threadpool % workers)
            print(f"Using thread pool of {threadpool} instead.\n")
        threads=int(threadpool/workers)
        
        
        #write argument vars to logfile
        logfile.contents["run_var"]={"taxid": taxid ,"filesizelimit":filesizelimit , "threadpool": threadpool, "workers" : workers, "download_method": download_method,"accessions_limit": accessions_limit}
        logfile.contents["processed_acc"]={}
        logfile.contents["total_acc"] = accessions
        logfile.contents["run_info"]={"taxid":taxid, "sci_name": scientific_name, "n_total_acc": len(accessions), "command_issued": " ".join(sys.argv), "init_time": datetime.now().strftime("%d/%m/%Y %H:%M:%S")}
        logfile.update()
    elif conti == True:
        #if logfile.contents["run_info"]["init_time"]==None:
            #sys.exit(f"\nNo previous run initiation detected in {outputdir}. Exiting...")
        if logfile.contents["status"]== "completed":
            sys.exit(f"\nPrevious run initiated in {outputdir} has fully completed. There is nothing to run. Use --force to delete all previous run data in order to restart run.")
        #inherit run variables from previous run using logfile contents
        taxid, filesizelimit, threadpool, workers, download_method, accessions_limit = logfile.contents["run_var"].values()
        _, scientific_name, _, command_issued, init_time = logfile.contents["run_info"].values()
        threads=int(threadpool/workers)
        accessions = logfile.contents.get("total_acc")
        print(f"\nPrevious incomplete run initiated on {init_time} detected:\n{command_issued}\n\nResuming run...\n")        
        print(f"Organism name: {scientific_name} (NCBI TaxID: {taxid})\n")
        print(f"Total accessions fetched from ENA: {len(accessions)}\n")
    
    #making subdirs if needed
    if not os.path.exists(fastqdir):
        os.makedirs(fastqdir)
    #make file to write processed files
    pathtoprocessed= os.path.join(outputdir, "GetExpressionMatrix", "processed.tsv")
    if not os.path.exists(pathtoprocessed):
        with open(pathtoprocessed, "w") as f:
            f.write("Accession\tMap_rate\n")
    
    #load processed files into log and  update logfile.\
    logfile.update()
    logfile.load()
    with open(pathtoprocessed, "r") as f:
        logfile.contents["processed_acc"] = {chunk.split("\t")[0]: chunk.split("\t")[1] for chunk in f.read().split("\n") if chunk != "Accession\tMap_rate" and chunk != ""}
        logfile.contents["processed_acc"] = {key : value if "failed" in value or "exception" in value else float(value) for key, value in logfile.contents["processed_acc"].items()}
    logfile.update()
    
    #directory to hold kallisto intermediate files
    if not os.path.exists(kaldir):
        os.makedirs(kaldir)
    if not os.path.exists(indexpath):
        print(f'CDS at {fastapath} will be used to build Kallisto index at {indexpath}.')
        print(fastapath)
        result= misc.run_with_retries(retrylimit, read_map.launch_kallisto_index, [fastapath,indexpath],
        f"Kallisto index failed. Retrying...",
        f"\nContructing Kallisto index ...\n")
        if result == "failed":
            sys.exit(f"Kallisto index step failed after {retrylimit} retries. Exiting...")
    
    print(f"Initiating parallel download and pseudoalignment of {min(accessions_limit, len(accessions))} accessions...\n")
    #subseting accessions based on accession limit
    accessions = accessions[0:min(accessions_limit, len(accessions))]
    
    
    #initiate parallel DL and PS of accessions
    parallel_job(workers)
    logfile.load()
    with open(pathtoprocessed, "r") as f:
            logfile.contents["processed_acc"] = {chunk.split("\t")[0]: chunk.split("\t")[1] for chunk in f.read().split("\n") if chunk != "Accession\tMap_rate" and chunk != ""}
            logfile.contents["processed_acc"] = {key : value if "failed" in value or "exception" in value else float(value) for key, value in logfile.contents["processed_acc"].items()}
    logfile.update()
    print("\nChecking logs to re-attempt download of failed accessions....")
    parallel_job(workers)
    logfile.load()
    with open(pathtoprocessed, "r") as f:
        logfile.contents["processed_acc"] = {chunk.split("\t")[0]: chunk.split("\t")[1] for chunk in f.read().split("\n") if chunk != "Accession\tMap_rate" and chunk != ""}
        logfile.contents["processed_acc"] = {key : value if "failed" in value or "exception" in value else float(value) for key, value in logfile.contents["processed_acc"].items()}
    logfile.update()
    #transpose expression matrix
    returncode = classify.mat_transposer(tpm_matpath, tpm_matpath_T)
    if returncode == 0:
        sys.exit("Error transposing expression matrix. Exiting...")
    
    logfile.contents["status"]= "completed"
    logfile.update()
    print(f"Expression Matrix generated in {tpm_matpath_T}")
    print("GetExpressionMatrix.py finished running on ", datetime.now().strftime("%d/%m/%Y %H:%M:%S"))