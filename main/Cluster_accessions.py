#setting sys.path for importing modules
import os
import sys
if __name__ == "__main__":
        abspath= os.getcwd()
        parent_module= os.path.join(abspath.split("LSTrAP-denovo")[0], "LSTrAP-denovo")
        sys.path.insert(0, parent_module)

import argparse
import concurrent.futures

from assembly import misc

if __name__ == "__main__":
    
    #arguments
    parser= argparse.ArgumentParser(description="LSTrAP-denovo.Cluster_accessions: Large-scale download , quality control and clustering of accessions based on transcriptome profiles.\n \
    NOTE: This is step 2 of 3 in the LSTrAP-denovo pipeline. Requires prior run of step 1: Preliminary_assembly.py. \n\
    Refer to https://github.com/pengkenlim/LSTrAP-denovo for more information on pipeline usage and implmentation")
    parser.add_argument("-o", "--output_dir", type=str, metavar= "", required=True,
    help= "Directory for data output. Directory needs to be same as for Preliminary_assembly.py run prior to .")
    parser.add_argument("-ps", "--pseudoalignment_threshold", type=int ,metavar="", default=0 , choices=range(0, 70+1),
    help = "Specifies pseudoalignment threshold (%%PS) during quality control. Accessions that do not meet this threshold will be discarded and not be clustered. Threshold will be determined automatically based on kernel-density minima of %%PS by default.")
    parser.add_argument("-s","--filesizelimit" , type=int, metavar="", default=500 , choices=range(100, 1500),
    help="Specifies the size limit(mb) of accession read files to partially download. Limit set to 500 (mb) by default. Has a direct impact on the download time, pseudoalignment runtime in this step of the pipeline.\
    User is advised to reduce size limit when downloading and processing >500 accessions")
    parser.add_argument("-t", "--threads", type=int, metavar="", default=4, 
    help = "Total thread pool for workers. Needs to be divisible by number of workers.")
    parser.add_argument("-w", "--workers", type=int, metavar="", default=4, 
    help= "Specify the maximum workers for running multiple download-pseudoalignment jobs in parellel. Set to 2 by default.")
    parser.add_argument("-dm", "--download_method", type=str, metavar="", default="ftp", choices=["ascp","ftp"],
    help = "Method to download accession runs. ftp/ascp.")  
    parser.add_argument("-al", "--accessions_limit", type=int, metavar="", default=500,
    help= "Specifies the limit for number of accessions to download and process. Accessions will be selected from a pre-randomised list that was fetched during Preliminary_assembly.py run and stored in in the logs.json file.\
    Default set to 500.")
    parser.add_argument("-kr", "--k_range", type=str, metavar="", default="2:20",
    help = "Specifies the range of k (number of clusters) to iterate through during k-means clustering. Lower and upper limit seperated by colon(:). \
    Set to 2:20 by default. Increasing the upper limit for high-resolution clustering is advised when downloading and processing > 500 accessions or if accessions are expected to have various experimental permutations (genotypes, pertubations, tissue-types) ." )    
    parser.add_argument("-kr", "--k_range", type=str, metavar="", default="2:20",
    
    parser.add_argument("-con", "--conti", action="store_true",
    help = "Resume incomplete run based on output directory. Only requires -o to run.")
    
    #parse args
    args=parser.parse_args()

    #assign output_dir and conti to variables
    outputdir= args.output_dir
    conti=args.conti
    #check if outputdir and its subdir from previous step exists. If not, exit.
    fastqdir=os.path.join(outputdir,"fastq")
    ssadir=os.path.join(outputdir, "ssa")
    if not os.path.exists(fastqdir) or not os.path.exists(ssadir):
        sys.exit("Error in output directory. Either directory or its required sub-directory are not found. Exiting...")
    
    #create logfile object, exit if step 1 not completed.
    logfile=misc.logfile(os.path.join(outputdir,"logs.json"))
    if  logfile["prelim"]["status"] != "complete":
        sys.exit("Prior Preliminary_assembly.py run (step 1) was either incomplete or not started at all. \nNOTE: This is step 2 of 3 in the LSTrAP-denovo pipeline.\nExiting...")
    
    #assigning arguments to variables, writing to log OR fetching variables from log
    if conti==False:
        pseudoalignment_threshold = args.pseudoalignment_threshold
        filesizelimit= args.filesizelimit
        threadpool= args.threads
        workers=args.workers
        download_method= args.download_method
        accessions_limit= args.accessions_limit
        k_range= args.k_range
        
        #check if threads/worker arguments are logical. Correct if necessary.
        if threadpool % workers != 0:
            print(f"Specified thread pool of {threadpool} is not divisible by number of workers.")
            threadpool= threadpool - (threadpool % workers)
            print(f"Using thread pool of {threadpool} instead.\n")
        threads=int(threadpool/workers)