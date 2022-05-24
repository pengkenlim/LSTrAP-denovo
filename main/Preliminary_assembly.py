#setting sys.path for importing modules
import os
import sys
if __name__ == "__main__":
        abspath= os.getcwd()
        parent_module= os.path.join(abspath.split("LSTrAP-denovo")[0], "LSTrAP-denovo")
        sys.path.insert(0, parent_module)

import argparse
import concurrent.futures
import random
import numpy as np
from time import sleep
from datetime import datetime

from tqdm import tqdm
from download import ena
from download import aspera
from assembly import soapdenovo, misc , consensus, report
from setup import constants
from preprocess import trim




def single_sample_assembly(accession,index):
    '''Job to generate Single-sample assembly. 
    Validate download path -> download via ftp/ascp-> read trimming by fastp -> assembly by soapdenovo-Trans'''
    logfile.load()
    if accession in logfile.contents["prelim"]["processed_acc"].keys():
        print(f"{accession} already processed. Skipping accession...")
        return f"{accession} already processed."
    #to un-sync workers 
    sleep((index%workers)*30)
    ascp_fullpath,ftp_fullpath = logfile.contents["prelim"]["run_var"]["selected_accessions"].get(accession)
    fastqpath=os.path.join(fastqdir,accession+".fastq.gz")
    #use the appropriate download method to download accession
    if download_method == "ascp":
        result= misc.run_with_retries(retrylimit,
        aspera.launch_ascp,
        [ascp_fullpath,fastqpath,filesizelimit],
        f"{accession}: Download failed. Retrying...",
        f"{accession}: Downloading file via ascp...\n")
        
    elif download_method == "ftp":
        result= misc.run_with_retries(retrylimit,
        aspera.launch_curl,
        [ftp_fullpath,fastqpath,filesizelimit],
        f"{accession}: Download failed. Retrying...",
        f"{accession}: Downloading file via ftp...\n")
    if result == "failed":
        logfile.load()
        logfile.contents["prelim"]["processed_acc"][accession]= "Download failed."
        logfile.update()
    #trim and uncompress
    result= misc.run_with_retries(retrylimit,
    trim.launch_fastp,
    [fastqpath, fastqpath.split(".gz")[0],threads],
    f"{accession}: Fastp trimming failed. Retrying...",
    f"{accession}: Trimming file using Fastp...\n")
    if result == "failed":
        logfile.load()
        logfile.contents["prelim"]["processed_acc"][accession]= "Fastp failed."
        logfile.update()
        return f"{accession}: Aborted after {retrylimit} retries."
    
    #make config file for soapdenovotrans to parse
    fastqpath= fastqpath.split(".gz")[0]
    configoutpath = os.path.join(ssadir, accession + "_temp.config")
    soapdenovo.make_config(fastqpath,configoutpath)
    #Single-sample-assembly process
    outputpath_prefix= os.path.join(ssadir, accession)
    results=misc.run_with_retries(retrylimit,
    soapdenovo.launch_soap,
    [configoutpath, kmerlen, outputpath_prefix, threads],
    f"{accession}: Soapdenovo-Trans failed. Retrying...",
    f"{accession}: Assembling transcripts with Soapdenovo-Trans...\n")
    if result == "failed":
        logfile.load()
        logfile.contents["prelim"]["processed_acc"][accession]= "Assembly failed."
        logfile.update()
        return f"{accession}: Aborted after {retrylimit} retries."
    
    #renaming transcript file to keep and deleting others
    os.system(f"mv {outputpath_prefix}_temp.scafSeq {outputpath_prefix}.fasta")
    os.system(f"rm -r {outputpath_prefix}_temp*")

    #remove uncompressed and trimmed fastq file to save space
    os.system(f"rm {fastqpath}")
        
    #extract orf from assembly to get cds.fasta
    results= misc.run_with_retries(retrylimit,
    soapdenovo.extract_orf,
    [outputpath_prefix + ".fasta", outputpath_prefix + "_cds.fasta", orfminlen, startcodon , geneticcode],
    f"{accession}: ORFfinder failed. Retrying...",
    f"{accession}: Extracting CDS with ORFfinder...\n")
    if result == "failed":
        logfile.load()
        logfile.contents["prelim"]["processed_acc"][accession]= "Assembly failed."
        logfile.update()
        return f"{accession}: Aborted after {retrylimit} retries."
    os.system(f"rm {outputpath_prefix}.fasta")
    n_cds, _, _ = misc.get_assembly_stats(outputpath_prefix + "_cds.fasta")
    logfile.load()
    logfile.contents["prelim"]["processed_acc"][accession]= n_cds
    logfile.update()
    print(f"{accession}: Single-sample assembly completed.")
    return f"{accession} processed"    
    
def parallel_ssa(workers):
    ''' Wrapper to parallelize SSA jobs. 
    Includes progress bar visualisation.'''
    logfile.load()
    with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as executor:
                progress_bar= tqdm(total=len(logfile.contents["prelim"]["run_var"]["selected_accessions"]), desc= "Accessions processed", unit="Acsn", leave=True)
                results= [executor.submit(single_sample_assembly, accession, index) for index, accession in enumerate(logfile.contents["prelim"]["run_var"]["selected_accessions"])]
                for f in concurrent.futures.as_completed(results):
                    if "processed" in f.result():
                        progress_bar.update(1)
                        #print("\n")
                        progress_bar.set_postfix_str(s=f.result())
                        print("\n")
                    else:
                        print(f.result())
                progress_bar.close()
                logfile.load()
                
                #conditional to sense when something is really wrong (i.e. every accession fails)
                if len([k for k in logfile.contents["prelim"]["processed_acc"].values() if type(k) is int ]) ==0:
                    sys.exit("Unexpected error occured. Exiting...")

                    
                    
                

def ssa_consensus(assemblydir):
    """ generate consensus assembly from SSAs"""
    logfile.load()
    print("\nConcatenating all Single-sample assemblies...\n")
    #make subdir within assemblydir to store consensus assemblies
    outputdir=os.path.join(assemblydir, "concat")
    if not os.path.exists(outputdir):
        os.makedirs(outputdir)
    #concatenate and rename transcript ids
    concatname="ssa_concat_cds.fasta"
    concatpath=os.path.join(outputdir,concatname)
    consensus.concat_rename_assemblies(assemblydir,concatpath)
    #launch CD-HIT in shell with retry wrapper
    clstr_concatpath= os.path.join(outputdir,"ssa_concat_cds_CT1.fasta")
    result=misc.run_with_retries(retrylimit, 
    consensus.launch_cdhit, 
    [concatpath,0.998,clstr_concatpath, threadpool],
    "CD-HIT-EST failed. Retrying...\n", 
    "Running CD-HIT-EST on concatenated assembly...\n")
    if result == "failed":
        sys.exit(f"CD-HIT-EST aborted after aborted after {retrylimit} retries. Exiting...")
    clstrinfopath = clstr_concatpath + ".clstr" 
    print("Parsing output from CD-HIT-EST and extracting sequences...\n")
    #extract sequence IDs to retain for each consensus threshold
    logfile.contents["prelim"]["consensus"]["stats"]={}
    logfile.update()
    #temporary list to hold n_CDS
    temp_list=[]
    for n_threshold in range(1,len([cds for cds in os.listdir(assemblydir) if "cds.fasta" in cds])+1):
        seq_to_retain= consensus.cluster_seq_extractor(n_threshold,clstrinfopath)
        consensus_ssa_path= os.path.join(outputdir,f"ssa_concat_cds_CT{n_threshold}.fasta")
        consensus.fasta_subset(clstr_concatpath, consensus_ssa_path, seq_to_retain)
        print(f"Preliminary assembly generated at {consensus_ssa_path} using consensus threshold of {n_threshold}.")
        n_cds, avg_cds_len, GC = misc.get_assembly_stats(consensus_ssa_path)
        print(f"No. of CDS: {n_cds}\nAvg. CDS len: {avg_cds_len}\nGC content: {GC}%\n")
        logfile.contents["prelim"]["consensus"]["stats"][n_threshold]= [n_cds, avg_cds_len ,GC, consensus_ssa_path]
        temp_list+=[n_cds]
    logfile.update()
    target_cds= np.median([k for k in logfile.contents["prelim"]["processed_acc"].values() if type(k) is int])
    logfile.contents["prelim"]["consensus"]["optimal"]= int(consensus.CT_from_target_CDS(temp_list,target_cds))
    logfile.update()
    print("Consensus threshold of " + str(logfile.contents["prelim"]["consensus"]["optimal"])+" has been determined to be optimal.\n")


     
    
        
        
        
    


if __name__ == "__main__":
	#retry limit determines the number of retries before aborting
    retrylimit=2
    
    #arguments
    parser= argparse.ArgumentParser(description="LSTrAP-denovo.Preliminary_assembly: Assemble a reduced but high-confidence assembly from public RNA-seq data")
    parser.add_argument("-o", "--output_dir", type=str, metavar= "", required=True,
    help= "Directory for data output.")
    parser.add_argument("-k", "--kmer_len", type=int, metavar="", default=35, choices=range(21, 49+1,2), 
    help = "Specifies K-mer length (odd integer only) for assembly using Soapdenovo-Trans. K-mer length will be set to 35 by default.")
    ##parser.add_argument("-ct", "--consensus_threshold", type=int ,metavar="", default=0 , choices=range(0, 10+1),
    ##help = "Specifies consensus threshold during filtering. Threshold will be determined automatically by default.")
    parser.add_argument("-s","--filesizelimit" , type=int, metavar="", default=1500, 
    help="Specifies the parital download limit/ file size requirement(mb) of accession read files. Limit set to 1500 (mb) by default.")
    parser.add_argument("-t", "--threads", type=int, metavar="", default=4, 
    help = "Total thread pool for workers. Needs to be divisible by number of workers.")
    parser.add_argument("-w", "--workers", type=int, metavar="", default=2, 
    help= "Specifies the maximum workers for running multiple download-assembly jobs in parallel. Set to 2 by default.")
    parser.add_argument("-g", "--gene_code", type=int, metavar="", default=1, choices=range(1, 31), 
    help= "Genetic code (codon table) passed to ORFfinder during ORF extraction. Set to 1 (universal) by default. Refer to https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi for more information.")
    parser.add_argument("-sc", "--start_codon", type=int, metavar="", default=0, choices=range(0, 2+1),
    help= "ORF start codon passed to ORFfinder during ORF extraction. Set to 0 (ATG only) by default. Refer to ORFfinder usage https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/ORFfinder/USAGE.txt for more information")
    parser.add_argument("-ml", "--min_len", type=int, metavar="", default=300, choices=range(30, 500),
    help= "Minimal ORF length (nt) passed to ORFfinder during ORF extraction. Set to 300 by default.")
    parser.add_argument("-dm", "--download_method", type=str, metavar="", default="ascp", choices=["ascp","ftp"],
    help = "Method to download accession runs. ftp/ascp.")
    parser.add_argument("-na", "--n_accessions", type=int, metavar="", default=10, choices=range(10,50+1), 
    help = "Number of single-accession-assemblies to combine in order to generate the preliminary assembly.")
    parser.add_argument("-a", "--accessions", type=str, metavar="",
    help= "User-defined list of SRA run accessions to fetch for preliminary assembly. If insufficient accessions provided, run will be supplemented with other public accessions. E.g.: SRR123456,SRR654321,ERR246810,...")
    
    #mutually excusive args for initial run or resume incomplete run
    ME_group_1 = parser.add_mutually_exclusive_group(required=True)
    ME_group_1.add_argument("-i", "--id", type=int, metavar="", 
    help= "NCBI TaxID of organism for fetching SRA run accessions.")
    ME_group_1.add_argument("-con", "--conti", action="store_true",
    help = "Resume incomplete run based on output directory. Only requires -o to run.")
    
    args=parser.parse_args()
    outputdir= args.output_dir
    conti=args.conti
    #create outputdir , fastqdir and ssadir if not found
    if not os.path.exists(outputdir):
        os.makedirs(outputdir)    
    fastqdir=os.path.join(outputdir,"prelim","fastq")
    ssadir=os.path.join(outputdir, "prelim","ssa")
    if not os.path.exists(fastqdir):
        os.makedirs(fastqdir)
    if not os.path.exists(ssadir):
        os.makedirs(ssadir)
    logfile=misc.logfile(os.path.join(outputdir,"logs.json"))
    
    #assigning arguments to variables, writing to log OR fetching variables from log
    if conti==False:
        taxid= args.id
        selected_accessions= args.accessions
        ##consensus_threshold= args.consensus_threshold
        filesizelimit= args.filesizelimit * 1048576
        threadpool= args.threads
        workers=args.workers
        kmerlen=args.kmer_len
        orfminlen=args.min_len
        startcodon=args.start_codon
        geneticcode=args.gene_code
        download_method= args.download_method
        n_accessions = args.n_accessions
        
        #getting sciname and accessions from ena using taxid
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
        #check if there is a previous run in the same outputdir. if so, exit and print error message
        if logfile.contents["prelim"]["run_info"]["init_time"] is not None:
            if logfile.contents["prelim"]["status"] == "completed":
                sys.exit(f"Previous completed run detected in {outputdir}. Exiting...")
            else:
                sys.exit(f"Previous incomplete run detected at {outputdir}.\nUse either -con to continue previous run or remove output directory to start a fresh run.\nExiting...")
        
        #check if accessions are given. if not, select accessions from total accessions.
        #check file size.
        if selected_accessions is not None:
            selected_accessions = selected_accessions.split(",")
            print("Checking file sizes of accessions provided by user...")
            selected_accessions_dict={}
            for accession in selected_accessions:
                ascp_fullpath, ftp_fullpath, filesize = aspera.get_download_path_ffq(accession)
                if filesize >= filesizelimit:
                    selected_accessions_dict[accession]=(ascp_fullpath,ftp_fullpath)
                else:
                    print(f"{accession} not included due to insufficient file size.")
                if len(selected_accessions_dict)==n_accessions:
                    break
            if len(selected_accessions_dict)<n_accessions:
                print(f"User-provided accessions are insufficient. Will supplement with other accessions...")
                for accession in accessions:
                    ascp_fullpath, ftp_fullpath, filesize = aspera.get_download_path_ffq(accession)
                    if filesize >= filesizelimit:
                        selected_accessions_dict[accession]=(ascp_fullpath,ftp_fullpath)
                    if len(selected_accessions_dict)==n_accessions:
                        break
        else:
            print("Selecting accessions with appropriate file sizes to build preliminary assembly...")
            selected_accessions_dict={}
            for accession in accessions:
                ascp_fullpath, ftp_fullpath, filesize = aspera.get_download_path_ffq(accession)
                if filesize >= filesizelimit:
                    selected_accessions_dict[accession]=(ascp_fullpath,ftp_fullpath)
                if len(selected_accessions_dict)==n_accessions:
                    break
                
        if threadpool % workers != 0:
            print(f"Specified thread pool of {threadpool} is not divisible by number of workers.")
            threadpool= threadpool - (threadpool % workers)
            print(f"Using thread pool of {threadpool} instead.\n")
        threads=int(threadpool/workers)
        
        #write information relavent to fresh run into log file
        logfile.contents["prelim"]["run_var"]={"taxid":taxid,
        "selected_accessions":selected_accessions_dict,
        "outputdir": outputdir,
        ##"consensus_threshold": consensus_threshold,
        "filesizelimit":filesizelimit,
        "threadpool":threadpool,
        "workers":workers,
        "kmerlen": kmerlen,
        "orfminlen": orfminlen,
        "geneticcode": geneticcode,
        "startcodon": startcodon,
        "download_method":download_method,
        "n_accessions": n_accessions}
        
        logfile.contents["prelim"]["run_info"]={"taxid":taxid,
        "sci_name": scientific_name, "n_total_acc": len(accessions), "command_issued": " ".join(sys.argv), "init_time": datetime.now().strftime("%d/%m/%Y %H:%M:%S")}
        logfile.contents["prelim"]["total_acc"]= accessions
        logfile.contents["prelim"]["processed_acc"]={}
        logfile.update()

    elif conti==True:
        #exit if log file contains command args
        if logfile.contents["prelim"]["run_info"]["init_time"]==None:
            sys.exit(f"\nNo previous run initiation detected in {outputdir}. Exiting...")
        if logfile.contents["prelim"]["status"]== "completed":
            sys.exit(f"\nPrevious run initiated in {outputdir} has fully completed. There is nothing to run.")
        #taxid, selected_accessions_dict, outputdir, consensus_threshold, filesizelimit, threadpool, workers, kmerlen , orfminlen, geneticcode, startcodon ,download_method, n_accessions = logfile.contents["prelim"]["run_var"].values()
        taxid, selected_accessions_dict, outputdir, filesizelimit, threadpool, workers, kmerlen , orfminlen, geneticcode, startcodon ,download_method, n_accessions = logfile.contents["prelim"]["run_var"].values()
        _, scientific_name, _, command_issued, init_time = logfile.contents["prelim"]["run_info"].values()
        accessions = logfile.contents["prelim"].get("total_acc")
        print(f"\nPrevious incomplete run initiated on {init_time} detected:\n{command_issued}\n\nResuming run...\n")
        
        if threadpool % workers != 0:
            threadpool= threadpool - (threadpool % workers)
        threads=int(threadpool/workers)
    
    logfile.load()
    #assemble SSAs in parallel
    parallel_ssa(workers)   
    ssa_consensus(ssadir)
    logfile.load()
    logfile.contents["prelim"]["status"]= "completed"
    logfile.update()
    print("LSTrAP-denovo.Preliminary_assembly.py completed.\nGenerating html report...")
    report.generate_from_json_log(logfile.path, os.path.join(outputdir, "LSTrAP-denovo.html"))
    
