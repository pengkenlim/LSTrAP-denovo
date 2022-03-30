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
from time import sleep

from tqdm import tqdm
from download import ena
from download import aspera
from assembly import soapdenovo, misc
from setup import constants
from preprocess import trim

#job to generate Single-sample assembly
def single_sample_assembly(accession,index):
    global processed_accessions
    if accession in processed_accessions:
        return f"{accession} processed"
    #to un-sync processes 
    sleep((index%workers)*5)
    #get download path and filesize of accession
    print(f"{accession}: checking file size...")
    ascp_fullpath, ftp_fullpath, filesize = aspera.get_download_path(accession)
    if filesize < filesizelimit:
        return f"Accession {accession} does not meet size requirement"
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
            return f"{accession}: aborted after {retrylimit} retries."

        #trim and uncompress
        result= misc.run_with_retries(retrylimit,
        trim.launch_fastp,
        [fastqpath, fastqpath.split(".gz")[0],threads],
        f"{accession}: Fastp trimming failed. Retrying...",
        f"{accession}: trimming file using Fastp...")
        if result == "failed":
            return f"{accession}: aborted after {retrylimit} retries."
        
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
            return f"{accession}: aborted after {retrylimit} retries."
        
        #remove uncompressed and trimmed fastq file to save space
        os.system(f"rm {fastqpath}")
        
        #extract orf from assembly to get cds.fasta
        results=misc.run_with_retries(retrylimit,
        soapdenovo.extract_orf,
        [outputpath_prefix + ".fasta", outputpath_prefix + "_cds.fasta", orfminlen, startcodon ,geneticcode],
        f"{accession}: ORFfinder failed. Retrying...",
        f"{accession}: Extracting CDS with ORFfinder...")
        if result == "failed":
            return f"{accession}: aborted after {retrylimit} retries."
        
        os.system(f"rm {outputpath_prefix}.fasta")
        print(f"{accession}: Single-sample assembly completed.")
        processed_accessions+=accession
        return f"{accession} processed"


def parellel_ssa(workers, selected_accessions):
    progress_bar= tqdm(total=len(selected_accessions), desc= "SSA of selected accessions", unit="Acsn", leave=True)
    with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as executor:
                results= [executor.submit(single_sample_assembly, accession, index) for index, accession in enumerate(accessions)]
                for f in concurrent.futures.as_completed(results):
                    if "processed" in f.result():
                        progress_bar.update(1)
                        progress_bar.set_postfix_str(s=f.result())
                    else:
                        print(f.result())

if __name__ == "__main__":
	#REMOVE THIS AFTER TESTING
    retrylimit=2
    
    #arguments
    parser= argparse.ArgumentParser(description="LSTrAP-denovo.preliminary_assembly: Assemble a reduced but high-confidence assembly from public RNA-seq data")
    parser.add_argument("-o", "--output_dir", type=str, metavar= "", required=True,
    help= "Directory for data output.")
    parser.add_argument("-k", "--kmer_len", type=int, metavar="", default=35, choices=range(21, 49+1,2), 
    help = "Specify K-mer length (odd integer only) for assembly using Soapdenovo-Trans. K-mer length will be set to 35 by default.")
    parser.add_argument("-ct", "--consensus_threshold", type=int ,metavar="", default=0 , choices=range(0, 10+1),
    help = "Specify consensus threshold during filtering. Threshold will be determined automatically by default.")
    parser.add_argument("-s","--filesizelimit" , type=int, metavar="", default=1500, 
    help="Specify the size limit(mb) of accession read files to partially download. Limit set to 1500 by default.")
    parser.add_argument("-t", "--threads", type=int, metavar="", default=4, 
    help = "Total thread pool for workers. Needs to be divisible by number of workers.")
    parser.add_argument("-w", "--workers", type=int, metavar="", default=2, 
    help= "Specify the maximum workers for running multiple download-assembly jobs in parellel. Set to 2 by default.")
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
        print(f"\nFetching RNA-seq accessions of {scientific_name}, NCBI TaxID {taxid} from ENA..\n")
        accessions = ena.get_runs(taxid)
        random.shuffle(accessions)
        print(f"Total accessions fetched from ENA: {len(accessions)}\n")
        #check if accessions are given. if not, select accessions from total accessions
        if selected_accessions is not None:
            selected_accessions = selected_accessions.split(",")
            if len(selected_accessions) < 10:
                sys.exit("Not enough accessions provided. Refer to --help for more information.")
        else:
            selected_accessions= accessions[:10]
                
        if threadpool % workers != 0:
            print(f"Specified thread pool of {threadpool} is not divisible by number of workers.")
            threadpool= threadpool - (threadpool % workers)
            print(f"Using thread pool of {threadpool} instead.\n")
        threads=int(threadpool/workers)
        
        #clear log file and write information relavent to fresh run
        logfile.clear()
        logfile.contents["prelim"]["cmd_args"]={"taxid":taxid,
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
        
        logfile.contents["prelim"]["ena"]={"assessions":accessions,
        "scientific_name": scientific_name}
        logfile.update()

    elif conti==True:
        #exit if log file contains command args
        if logfile.contents["prelim"]["cmd_args"]=={}:
            sys.exit(f"No previous run initiation detected in {outputdir}. Exiting...")

        print(f"Previous run detected. Resuming run...\n")
        taxid, selected_accessions, outputdir, consensus_threshold, filesizelimit, threadpool,workers, kmerlen , orfminlen, geneticcode, download_method, = logfile.contents["prelim"]["cmd_args"].values()
        accessions, scientific_name = logfile.contents["prelim"]["ena"].values()
    
    processed_accessions= logfile.contents["prelim"]["processed"]
    #assemble SSAs in parellel
    parellel_ssa(workers, selected_accessions)
