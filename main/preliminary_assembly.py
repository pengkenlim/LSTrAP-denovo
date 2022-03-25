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
from assembly import soapdenovo
from setup import constants
from preprocess import trim

#job to generate Single-sample assembly
def single_sample_assembly(accession):
    #a very stupid method to make sure that each process is out of sync. will have to fix
    sleep(random.randint(0,10))
    #get download path and filesize of accession
    ascp_fullpath, ftp_fullpath, filesize = aspera.get_download_path(accession)
    if filesize < filesizelimit:
        return f"Accession {accession} does not meet size requirement"
    else:
        #download
        fastqpath=os.path.join(fastqdir,accession+".fastq.gz")
        if download_method == "ascp":
            aspera.launch_ascp(ascp_fullpath, fastqpath, filesizelimit)
        elif download_method == "ftp":
            aspera.launch_curl(ftp_fullpath, fastqpath, filesizelimit)
        #trim and uncompress
        trim.launch_fastp(fastqpath,fastqpath.split(".gz")[0],threads)
        #mnake config file for soapdenovotrans to parse
        fastqpath= fastqpath.split(".gz")[0]
        configoutpath = os.path.join(ssadir, accession + "_temp.config")
        soapdenovo.make_config(fastqpath,configoutpath)
        #start assembly process
        outputpath_prefix= os.path.join(ssadir, accession)
        soapdenovo.launch_soap(configoutpath, kmerlen, outputpath_prefix, threads)
        #remove uncompressed and trimmed fastq file to save space
        os.system(f"rm {fastqpath}")
        #extract orf from assembly to get cds.fasta
        soapdenovo.extract_orf(outputpath_prefix + ".fasta", outputpath_prefix + "_cds.fasta", orfminlen, startcodon ,geneticcode)
        os.system(f"rm {outputpath_prefix}.fasta")
        return accession

def parellel_ssa(workers, accessions):
    #progress_bar= tqdm(total=len(accessions), desc= "Accessions processed", unit="Acsn", leave=True)
    with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as executor:
                results= [executor.submit(single_sample_assembly, accession) for accession in accessions]
                for f in concurrent.futures.as_completed(results):
                    #progress_bar.update(1)
                    #progress_bar.set_postfix_str(s=f.result())
                    print(f.result())

if __name__ == "__main__":
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
    parser.add_argument("-t", "--threads", type=int, metavar="", required=True, 
    help = "Total thread pool for workers. Needs to be divisible by number of workers.")
    parser.add_argument("-w", "--workers", type=int, metavar="", default=2, 
    help= "Specify the maximum workers for running multiple download-assembly jobs in parellel. Set to 2 by default.")
    parser.add_argument("-g", "--gene_code", type=int, metavar="", default=1, choices=range(1, 31), 
    help= "Genetic code (codon table) passed to ORFfinder during ORF extraction. Set to 1 (universal) by default. Refer to https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi for more information.")
    parser.add_argument("-sc", "--start_codon", type=int, metavar="", default=1, choices=range(0, 2+1),
    help= "ORF start codon passed to ORFfinder during ORF extraction. Set to 1 (ATG only) by default. Refer to ORFfinder usage https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/ORFfinder/USAGE.txt for more information")
    parser.add_argument("-ml", "--min_len", type=int, metavar="", default=300, choices=range(30, 500),
    help= "Minimal ORF length (nt) passed to ORFfinder during ORF extraction. Set to 300 by default.")
    parser.add_argument("-dm", "--download_method", type=str, metavar="", required=True, choices=["ascp","ftp"],
    help = "Method to download accession runs. ftp/ascp.")   
    
    #mutually excusive args fetch by taxid or by userdefined accessions
    ME_group_1 = parser.add_mutually_exclusive_group(required=True)
    ME_group_1.add_argument("-i", "--id", type=int, metavar="", 
    help= "NCBI TaxID of organism for fetching SRA run accessions.")
    ME_group_1.add_argument("-a", "--accessions", type=str, metavar="", 
    help= "User-defined list of SRA run accessions to fetch. Requires at least 10 accessions. E.g.: SRR123456,SRR654321,ERR246810,...")
    
    args=parser.parse_args()
    
    #assigning arguments to variables
    taxid= args.id
    accessions= args.accessions
    outputdir= args.output_dir
    consensus_threshold= args.consensus_threshold
    filesizelimit= args.filesizelimit * 1000000
    threadpool= args.threads
    workers=args.workers
    kmerlen=args.kmer_len
    orfminlen=args.min_len
    startcodon=args.start_codon
    geneticcode=args.gene_code
    download_method= args.download_method
    
    #check if threads pool is divisible by number of workers
    if threadpool % workers != 0:
        print(f"Specified thread pool of {threadpool} is not divisible by number of workers.")
        threadpool= threadpool - (threadpool % workers)
        print(f"Using thread pool of {threadpool} instead.")
    threads=int(threadpool/workers)
        
    #create outputdir , fastqdir and ssadir if not found
    if not os.path.exists(outputdir):
        os.makedirs(outputdir)
    
    fastqdir=os.path.join(outputdir,"fastq")
    ssadir=os.path.join(outputdir, "ssa")
    
    if not os.path.exists(fastqdir):
        os.makedirs(fastqdir)
        
    if not os.path.exists(ssadir):
        os.makedirs(ssadir)

	#check if accessions are given
    if accessions is not None:
        accessions = accessions.split(",")
        if len(accessions) < 10:
            sys.exit("Not enough accessions provided. Refer to --help for more information.")
        else:
            print("User defined accessions:")
            for index , accession in enumerate(accessions):
                print(f"{index +1}. ",accession)
            
            #assemble SSAs in parellel
            parellel_ssa(workers, accessions)
                    
                    
	#check if taxid is given
    elif taxid is not None:
        scientific_name= ena.get_sciname(taxid)
        if type(scientific_name) is  not list:
            sys.exit("TaxID {taxid} is invalid/not found.")
        elif len(scientific_name) > 1:
            sys.exit("More than one organism found for TaxID {taxid}.")
        else:
            scientific_name= scientific_name[0]["scientific_name"]
            print(f"Fetching RNA-seq accessions of {scientific_name}, NCBI TaxID {taxid} from ENA..")
            accessions = ena.get_runs(taxid)
            print(f"Total accessions: {len(accessions)}")
        
        random.shuffle(accessions)
        accessions=accessions[:10]
        print("Randomly selected accessions:")
        for index , accession in enumerate(accessions):
            print(f"{index +1}. ",accession)
        #assemble SSAs in parellel
        parellel_ssa(workers, accessions)

