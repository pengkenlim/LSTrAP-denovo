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

from tqdm import tqdm
from download import ena
from download import aspera
from assembly import soapdenovo, misc , consensus
from setup import constants
from preprocess import trim
from datetime import datetime


def single_sample_assembly(accession,index):
    '''Job to generate Single-sample assembly. 
    Validate download path -> download via ftp/ascp-> read trimming by fastp -> assembly by soapdenovo-Trans'''
    logfile.load()
    #check if accessions is already processed in the case of a resumed run.
    if accession in logfile.contents["prelim"]["processed_acc"]:
        return f"{accession} processed"
    #to un-sync processes 
    sleep((index%workers)*5)
    #get download path and filesize of accession
    print(f"{accession}: checking file size...")
    ascp_fullpath, ftp_fullpath, filesize = aspera.get_download_path(accession)
    if filesize < filesizelimit:
        logfile.load()
        logfile.contents["prelim"]["processed_acc"][accession]= "File size requirements not met."
        logfile.update()
        return f"{accession}: Aborted. Filesize {np.round(filesize/1000000)} mb is below limit requirements. "
    else:
        #download
        fastqpath=os.path.join(fastqdir,accession+".fastq.gz")
        if download_method == "ascp":
            
            result= misc.run_with_retries(retrylimit, 
            aspera.launch_ascp, 
            [ascp_fullpath,fastqpath,filesizelimit],
            f"{accession}: Download failed. Retrying...", 
            f"{accession}: downloading file via ascp...")
            
        elif download_method == "ftp":
            
            result= misc.run_with_retries(retrylimit,
            aspera.launch_curl,
            [ftp_fullpath,fastqpath,filesizelimit],
            f"{accession}: Download failed. Retrying...",
            f"{accession}: downloading file via ftp...")
        if result == "failed":
            logfile.load()
            logfile.contents["prelim"]["processed_acc"][accession]= "Download failed."
            logfile.update()
            return f"{accession}: aborted after {retrylimit} retries."

        #trim and uncompress
        result= misc.run_with_retries(retrylimit,
        trim.launch_fastp,
        [fastqpath, fastqpath.split(".gz")[0],threads],
        f"{accession}: Fastp trimming failed. Retrying...",
        f"{accession}: trimming file using Fastp...")
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
        f"{accession}: Assembling transcripts with Soapdenovo-Trans...")
        if result == "failed":
            logfile.load()
            logfile.contents["prelim"]["processed_acc"][accession]= "Assembly failed."
            logfile.update()
            return f"{accession}: Aborted after {retrylimit} retries."
        
        #remove uncompressed and trimmed fastq file to save space
        os.system(f"rm {fastqpath}")
        
        #extract orf from assembly to get cds.fasta
        results=misc.run_with_retries(retrylimit,
        soapdenovo.extract_orf,
        [outputpath_prefix + ".fasta", outputpath_prefix + "_cds.fasta", orfminlen, startcodon ,geneticcode],
        f"{accession}: ORFfinder failed. Retrying...",
        f"{accession}: Extracting CDS with ORFfinder...")
        if result == "failed":
            logfile.load()
            logfile.contents["prelim"]["processed_acc"][accession]= "Assembly failed."
            logfile.update()
            return f"{accession}: Aborted after {retrylimit} retries."
        
        os.system(f"rm {outputpath_prefix}.fasta")
        print(f"{accession}: Single-sample assembly completed.")
        return f"{accession} processed"


def parellel_ssa(workers):
    ''' Wrapper to parellelize SSA jobs. 
    Includes progress bar visualisation.'''
    logfile.load()
    progress_bar= tqdm(total=len(logfile.contents["prelim"]["run_var"]["selected_accessions"]), desc= "Accessions processed", unit="Acsn", leave=True)
    with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as executor:
                results= [executor.submit(single_sample_assembly, accession, index) for index, accession in enumerate(selected_accessions)]
                for f in concurrent.futures.as_completed(results):
                    if "processed" in f.result():
                        progress_bar.update(1)
                        progress_bar.set_postfix_str(s=f.result())
                        print("\n")
                    else:
                        print(f.result())
                logfile.load()
                #conditional to sense when something is really wrong (i.e. every accession fails)
                if len(logfile.contents["prelim"]["processed_acc"]) ==0:
                    sys.exit("Unexpected error occured. Exiting...")
                #rerun to get at least 10 successful runs
                while len(logfile.contents["prelim"]["processed_acc"]) <10:
                    logfile.contents["prelim"]["run_var"]["selected_accessions"] = accessions[:len(logfile.contents["prelim"]["run_var"]["selected_accessions"])- len(logfile.contents["prelim"]["processed_acc"]) +10]
                    logfile.update()
                    results= [executor.submit(single_sample_assembly, accession, index) for index, accession in enumerate(selected_accessions)]
                    for f in concurrent.futures.as_completed(results):
                        if "processed" in f.result():
                            progress_bar.update(1)
                            progress_bar.set_postfix_str(s=f.result())
                            print("\n")
                        else:
                            print(f.result())
                    logfile.load()
                    
                    
                

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
    [concatpath,0.995,clstr_concatpath, threadpool],
    "CD-HIT-EST failed. Retrying...\n", 
    "Running CD-HIT-EST on concatenated assembly...\n")
    if result == "failed":
        sys.exit(f"CD-HIT-EST aborted after aborted after {retrylimit} retries. Exiting...")
    clstrinfopath = clstr_concatpath + ".clstr" 
    print("Parsing output from CD-HIT-EST and writing to log...\n")
    #extract sequence IDs to retain for each consensus threshold
    CT_seqid_dict={} #to hold SeqIDs for each threshold
    logfile.contents["prelim"]["consensus"]["CDS"]={}
    logfile.update()
    for n_threshold in range(1,len([cds for cds in os.listdir(assemblydir) if "cds.fasta" in cds])+1):
        seq_to_retain= consensus.cluster_seq_extractor(n_threshold,clstrinfopath)
        CT_seqid_dict[n_threshold]= seq_to_retain
        logfile.contents["prelim"]["consensus"]["CDS"][n_threshold]=len(seq_to_retain)
    logfile.update()
    #for auto determination of consensus_threshold
    if consensus_threshold == 0:
        logfile.contents["prelim"]["consensus"]["stats"]["CT"] = consensus.select_CT(list(logfile.contents["prelim"]["consensus"]["stats"].values()))
        logfile.update()
        print("Consensus threshold of " + str(logfile.contents["prelim"]["consensus"]["stats"]["CT"])+" has been determined automatically. Generating preliminary assembly....\n")
        consensus_ssa_path= os.path.join(outputdir,"ssa_concat_cds_CT"+ str(logfile.contents["prelim"]["consensus"]["stats"]["CT"])+".fasta")
        consensus.fasta_subset(clstr_concatpath, consensus_ssa_path, CT_seqid_dict[logfile.contents["prelim"]["consensus"]["stats"]["CT"]])
    else:
         print(f"Using user-defined consensus threshold of {consensus_threshold} to generate preliminary assembly....\n")
         consensus_ssa_path= os.path.join(outputdir,f"ssa_concat_cds_CT{consensus_threshold}.fasta")
         consensus.fasta_subset(clstr_concatpath, consensus_ssa_path, CT_seqid_dict[f"CT{consensus_threshold}"])
    logfile.contents["prelim"]["consensus"]["stats"]["path"]= consensus_ssa_path
    print("Calculating preliminary assembly statistics...\n")
    n_cds, avg_cds_len, GC = misc.get_assembly_stats(consensus_ssa_path)
    logfile.contents["prelim"]["consensus"]["stats"]["n_CDS"] = n_cds
    logfile.contents["prelim"]["consensus"]["stats"]["CDS_len"] = avg_cds_len
    logfile.contents["prelim"]["consensus"]["stats"]["GC"] = GC
    logfile.update()
     
    
        
        
        
    


if __name__ == "__main__":
	#retry limit determines the number of retries before aborting
    retrylimit=2
    
    #arguments
    parser= argparse.ArgumentParser(description="LSTrAP-denovo.Preliminary_assembly: Assemble a reduced but high-confidence assembly from public RNA-seq data")
    parser.add_argument("-o", "--output_dir", type=str, metavar= "", required=True,
    help= "Directory for data output.")
    parser.add_argument("-k", "--kmer_len", type=int, metavar="", default=35, choices=range(21, 49+1,2), 
    help = "Specifies K-mer length (odd integer only) for assembly using Soapdenovo-Trans. K-mer length will be set to 35 by default.")
    parser.add_argument("-ct", "--consensus_threshold", type=int ,metavar="", default=0 , choices=range(0, 10+1),
    help = "Specifies consensus threshold during filtering. Threshold will be determined automatically by default.")
    parser.add_argument("-s","--filesizelimit" , type=int, metavar="", default=1500, 
    help="Specifies the size limit(mb) of accession read files to partially download. Limit set to 1500 (mb) by default.")
    parser.add_argument("-t", "--threads", type=int, metavar="", default=4, 
    help = "Total thread pool for workers. Needs to be divisible by number of workers.")
    parser.add_argument("-w", "--workers", type=int, metavar="", default=2, 
    help= "Specifies the maximum workers for running multiple download-assembly jobs in parellel. Set to 2 by default.")
    parser.add_argument("-g", "--gene_code", type=int, metavar="", default=1, choices=range(1, 31), 
    help= "Genetic code (codon table) passed to ORFfinder during ORF extraction. Set to 1 (universal) by default. Refer to https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi for more information.")
    parser.add_argument("-sc", "--start_codon", type=int, metavar="", default=0, choices=range(0, 2+1),
    help= "ORF start codon passed to ORFfinder during ORF extraction. Set to 0 (ATG only) by default. Refer to ORFfinder usage https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/ORFfinder/USAGE.txt for more information")
    parser.add_argument("-ml", "--min_len", type=int, metavar="", default=300, choices=range(30, 500),
    help= "Minimal ORF length (nt) passed to ORFfinder during ORF extraction. Set to 300 by default.")
    parser.add_argument("-dm", "--download_method", type=str, metavar="", default="ftp", choices=["ascp","ftp"],
    help = "Method to download accession runs. ftp/ascp.")  
    parser.add_argument("-a", "--accessions", type=str, metavar="",
    help= "User-defined list of SRA run accessions to fetch for preliminary assembly. Requires at least 10 accessions. E.g.: SRR123456,SRR654321,ERR246810,...")
    
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
    fastqdir=os.path.join(outputdir,"fastq")
    ssadir=os.path.join(outputdir, "ssa")
    if not os.path.exists(fastqdir):
        os.makedirs(fastqdir)
    if not os.path.exists(ssadir):
        os.makedirs(ssadir)
    logfile=misc.logfile(os.path.join(outputdir,"logs.json"))
    
    #assigning arguments to variables, writing to log OR fetching variables from log
    if conti==False:
        taxid= args.id
        selected_accessions= args.accessions
        consensus_threshold= args.consensus_threshold
        filesizelimit= args.filesizelimit * 1000000
        threadpool= args.threads
        workers=args.workers
        kmerlen=args.kmer_len
        orfminlen=args.min_len
        startcodon=args.start_codon
        geneticcode=args.gene_code
        download_method= args.download_method
        
        #getting sciname and accessions from ena using taxid
        scientific_name= ena.get_sciname(taxid)
        if type(scientific_name) is  not list:
            sys.exit("TaxID {taxid} is invalid/not found. Exiting...")
        elif len(scientific_name) > 1:
            sys.exit("More than one organism found for TaxID {taxid}. Exiting...")
        scientific_name= scientific_name[0]["scientific_name"]
        print(f"\nFetching RNA-seq accessions of {scientific_name}( NCBI TaxID: {taxid}) from ENA..\n")
        accessions = ena.get_runs(taxid)
        random.shuffle(accessions)
        print(f"Total accessions fetched from ENA: {len(accessions)}\n")
        #check if accessions are given. if not, select accessions from total accessions
        if selected_accessions is not None:
            selected_accessions = selected_accessions.split(",")
            if len(selected_accessions) < 10:
                sys.exit("Not enough accessions provided. Refer to --help for more information.")
            #bring selected accessions to the front of the total accessions list
            accessions = list(set(selected_accessions + accessions))
        else:
            selected_accessions= accessions[:10]
                
        if threadpool % workers != 0:
            print(f"Specified thread pool of {threadpool} is not divisible by number of workers.")
            threadpool= threadpool - (threadpool % workers)
            print(f"Using thread pool of {threadpool} instead.\n")
        threads=int(threadpool/workers)
        
        #clear log file and write information relavent to fresh run
        logfile.clear("prelim")
        logfile.load()
        logfile.contents["prelim"]["run_var"]={"taxid":taxid,
        "selected_accessions":selected_accessions,
        "outputdir": outputdir,
        "consensus_threshold": consensus_threshold,
        "filesizelimit":filesizelimit,
        "threadpool":threadpool,
        "workers":workers,
        "kmerlen": kmerlen,
        "orfminlen": orfminlen,
        "geneticcode": geneticcode,
        "download_method":download_method}
        
        logfile.contents["prelim"]["run_info"]={"taxid":taxid,
        "sci_name": scientific_name, "n_total_acc": len(accessions), "command_issued": " ".join(sys.argv), "init_time": datetime.now().strftime("%d/%m/%Y %H:%M:%S")}
        logfile.contents["prelim"]["total_acc"]= accessions
        logfile.update()

    elif conti==True:
        #exit if log file contains command args
        if logfile.contents["prelim"]["run_info"]["init_time"]==None:
            sys.exit(f"\nNo previous run initiation detected in {outputdir}. Exiting...")
        if logfile.contents["prelim"]["status"]== "completed":
            sys.exit(f"\nPrevious run initiated in {outputdir} has fully completed. There is nothing to run.")
        taxid, selected_accessions, outputdir, consensus_threshold, filesizelimit, threadpool,workers, kmerlen , orfminlen, geneticcode, download_method, = logfile.contents["prelim"]["run_var"].values()
        accessions, scientific_name, _, command_issued, init_time, = logfile.contents["prelim"]["ena"].values()
        print(f"\nPrevious incomplete run: {command_issued} \ninitiated on {init_time} detected.\nResuming run...\n")
    
    logfile.load()
    logfile.contents["prelim"]["processed_acc"]={}
    logfile.update()
    #assemble SSAs in parellel
    parellel_ssa(workers)   
    ssa_consensus(ssadir)
    logfile.load()
    logfile.contents["prelim"]["status"]= "completed"
    logfile.update()
    
