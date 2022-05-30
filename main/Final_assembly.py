#setting sys.path for importing modules
import os
import sys
if __name__ == "__main__":
        abspath= os.getcwd()
        parent_module= os.path.join(abspath.split("LSTrAP-denovo")[0], "LSTrAP-denovo")
        sys.path.insert(0, parent_module)

import argparse
import concurrent.futures
import numpy as np
from time import sleep
from datetime import datetime
from tqdm import tqdm

from assembly import soapdenovo, misc , postprocess, report
from setup import constants
from preprocess import trim

#retry limit determines the number of retries before aborting
retrylimit=2


#definitions
def parallel_fastp(threadpool, accessions, sourcedir, targetdir, cluster):
    '''wrapper to praallelize fastp trimming of accessions. All input should be strings'''
    with concurrent.futures.ProcessPoolExecutor(max_workers=threadpool) as executor:
        progress_bar= tqdm(total=len(accessions), desc= f"Cluster {cluster}: Trimming accessions", unit="Acsn", leave=True)
        #results = [executor.submit(misc.run_with_retries, retrylimit, trim.launch_fastp,[os.path.join(sourcedir, accession+".fastq.gz"), os.path.join(targetdir, accession + ".fastq"), 1], f"{accession}: Fastp failed. retrying", "") for accession in accessions]#
        results = [executor.submit(trim.launch_fastp, os.path.join(sourcedir, accession+".fastq.gz"), os.path.join(targetdir, accession + ".fastq"), 1) for accession in accessions]
        for f in concurrent.futures.as_completed(results):
            print(f.result())
            progress_bar.update(1)


if __name__ == "__main__":      
    #arguments
    parser= argparse.ArgumentParser(description="LSTrAP-denovo.Final_assembly: Modular multi-kmer de novo transcriptome assembly of accession clusters for CDS mining.\n\
    NOTE: This is step 3 of 3 in the LSTrAP-denovo pipeline. Requires prior run of step 1: Preliminary_assembly.py and step 2: Cluster_accessions.py. \n\
    Refer to https://github.com/pengkenlim/LSTrAP-denovo for more information on pipeline usage and implmentation")
    parser.add_argument("-o", "--output_dir", type=str, metavar= "", required=True,
    help= "Directory for data output. Directory needs to be same as for step 1 (Preliminary_assembly.py) and step 2 (Cluster_accessions.py).")
    parser.add_argument("-t", "--threads", type=int, metavar="", default=16, 
    help = "Total threads for running CLI tools. Also serve as thread pool for parallel Fastp trimming.")
    parser.add_argument("-g", "--gene_code", type=int, metavar="", default=1, choices=range(1, 31), 
    help= "Genetic code (codon table) passed to ORFfinder during ORF extraction. Set to 1 (universal) by default. Refer to https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi for more information.")
    parser.add_argument("-sc", "--start_codon", type=int, metavar="", default=0, choices=range(0, 2+1),
    help= "ORF start codon passed to ORFfinder during ORF extraction. Set to 0 (ATG only) by default. Refer to ORFfinder usage https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/ORFfinder/USAGE.txt for more information")
    parser.add_argument("-ml", "--min_len", type=int, metavar="", default=300, choices=range(30, 500),
    help= "Minimal ORF length (nt) passed to ORFfinder during ORF extraction. Set to 300 by default.")
    parser.add_argument("-kmr", "--kmer_range", type=str, metavar="", default="25:49:6",
    help= "K-mer lengths for assembly. Lower limit , upper limit(inclusive) and step size seperated by colon(:). \
    Set to 25:49:6 by default (k-mer lengths of 25 , 31, 37, 43 and 49). NOTE: upper and lower limit must be odd integers while step size be an even integer. Range window (upper - lower) must be divisible by step size. Users are not advised to change this parameter.")
    parser.add_argument("-con", "--conti", action="store_true",
    help = "Resume incomplete run based on output directory. Only requires -o to run.")
    parser.add_argument("-f", "--force", action="store_true",
    help = "Delete data from previous Cluster_accessions.py run.")
    parser.add_argument("-rf", "--remap_filter" , default=1, choices=range(0, 1+1),
    help= "Remapping filter to get rid of exogenous CDSs. Enabled (1) by default. Set 0 to disable.")
    
    #parse args
    args=parser.parse_args()
    outputdir= args.output_dir
    conti=args.conti
    force=args.force
    
    #check if outputdir and its subdir from previous step exists. If not, exit.
    C_fastqdir=os.path.join(outputdir,"cluster","fastq")
    if not os.path.exists(C_fastqdir):
        sys.exit("Error in output directory. Either directory or its required sub-directory are not found. Exiting...")
        
    #check if logfile is present in output directory. exit if not there.
    if not os.path.exists(os.path.join(outputdir,"logs.json")):
        sys.exit("Error in output directory. logs.json file not found. Exiting...")
    
    #create logfile object, exit if step 1 not completed.
    logfile=misc.logfile(os.path.join(outputdir,"logs.json"))
    if  logfile.contents["cluster"]["status"] != "completed":
        sys.exit("Step 2 (Cluster_accessions.py) is either incomplete or not started at all. \nNOTE: This is step 3 of 3 in the LSTrAP-denovo pipeline.\nExiting...")
        
    #assigning arguments to variables, writing to log OR fetching variables from log
    if conti==False:
        if logfile.contents["final"]["run_info"].get("init_time") != None:
            if force != True:
                temp_string = logfile.contents["final"]["run_info"].get("init_time")
                sys.exit(f"Previous incomplete run started on {temp_string} has been detected in logfile.\nPlease use --conti to continue incomplete run or --force to clear previous run data. Exiting...")
            
            #Trigger run reset if force argument is used
            if force == True:
                print("\n--force argument has been specified by user. Deleting data from previous run and clearing logs....\n")
                logfile.clear("final")
                if os.path.exists(os.path.join(outputdir, "final")):
                    temp_string= os.path.join(outputdir, "final")
                    os.system(f"rm -r {temp_string}")
    
        #assign argument vars
        kmer_range = args.kmer_range
        threadpool= args.threads
        orfminlen=args.min_len
        startcodon=args.start_codon
        geneticcode=args.gene_code
        remap_filter= args.remap_filter
        #logic to check if kmer_range is valid
        kmer_min , kmer_max , kmer_step = kmer_range.split(":")
        kmer_min , kmer_max , kmer_step = int(kmer_min) , int(kmer_max) , int(kmer_step)
        if kmer_max - kmer_min <= kmer_step or  kmer_min%2 == 0 or kmer_max%2 ==0 or kmer_step%2 !=0 or kmer_min <21 or (kmer_max - kmer_min)%kmer_step != 0 : 
            sys.exit(f"Specified k-mer range {kmer_range} is invalid. Use -h for help information. Exiting...")
            
        #write run info (for final report) to logfile
        logfile.contents["final"]["run_info"]["init_time"]=  datetime.now().strftime("%d/%m/%Y %H:%M:%S")
        logfile.contents["final"]["run_info"]["command_issued"]= " ".join(sys.argv)
        
        #write run variables to logfile
        logfile.contents["final"]["run_var"]= {
        "kmer_range":kmer_range,
        "threadpool": threadpool,
        "orfminlen": orfminlen,
        "startcodon": startcodon,
        "geneticcode" : geneticcode,
        "remap_filter": remap_filter}
        logfile.update()
        
    elif conti==True:
        if logfile.contents["final"]["run_info"]["init_time"]==None:
            sys.exit(f"\nNo previous run initiation detected in {outputdir}. Exiting...")
            
        if logfile.contents["final"]["status"]== "completed":
            sys.exit(f"\nPrevious run initiated in {outputdir} has fully completed. There is nothing to run. Use --force to delete all previous run data in order to restart run.")
        
        print("\n--conti specified by user. Resuming previous incomplete run...\n")
        #inherit run variables when continuing run
        kmer_range, threadpool, orfminlen, startcodon, geneticcode, remap_filter = logfile.contents["final"]["run_var"].values()
        kmer_min , kmer_max , kmer_step = kmer_range.split(":")
        kmer_min , kmer_max , kmer_step = int(kmer_min) , int(kmer_max) , int(kmer_step)
        
    #get relavant info from previous steps
    #get accessions and sci name from step 1
    #cluster information from step 2
    _, scientific_name, _, _, _ = logfile.contents["prelim"]["run_info"].values()
    taxid, _, _, _, _, _, _ , _, _, _ ,_, _ = logfile.contents["prelim"]["run_var"].values()
    cluster_assignment_dict = logfile.contents["cluster"]["kmeans"].get("cluster_assignment_dict")
    total = logfile.contents["cluster"]["qc"].get("total")
    passed = logfile.contents["cluster"]["qc"].get("passed")
    clusters = list(cluster_assignment_dict.keys())
    clusters.sort()
    
    print(f"Organism name: {scientific_name} (NCBI TaxID: {taxid})\n")
    print(f"{len(total)} accessions downloaded and pseudoaligned in step 2. \n{len(passed)} accessions that passed QC has been assigned into {len(clusters)} clusters\n")
    for cluster in clusters:
        #fetch list of accessions
        accessions= cluster_assignment_dict.get(cluster)
        
        #create cluster dir and fastq dir
        cluster_dir = os.path.join(outputdir, "final" ,f"cluster_{cluster}")
        if not os.path.exists(cluster_dir):
            os.makedirs(cluster_dir)
        
        cluster_fastqdir = os.path.join(cluster_dir, "fastq")
        if not os.path.exists(cluster_fastqdir):
            os.makedirs(cluster_fastqdir)
        
        #make subdirs to hold assemblt intermediate files
        cluster_assemblydir = os.path.join(cluster_dir, "assembly")
        if not os.path.exists(os.path.join(cluster_assemblydir, "raw")):
            os.makedirs(os.path.join(cluster_assemblydir, "raw"))
        if not os.path.exists(os.path.join(cluster_assemblydir, "combined")):
            os.makedirs(os.path.join(cluster_assemblydir, "combined"))
        if not os.path.exists(os.path.join(cluster_assemblydir, "CPC2")):
            os.makedirs(os.path.join(cluster_assemblydir, "CPC2"))
        if not os.path.exists(os.path.join(cluster_assemblydir,"remap")):
            os.makedirs(os.path.join(cluster_assemblydir, "remap"))
        
        #tempdir to hold truncated but un-trimmed reads
        cluster_tempdir= os.path.join(cluster_dir, "temp")
        if not os.path.exists(cluster_tempdir):
            os.makedirs(cluster_tempdir)
        
        
        ##check size. Truncate libary if >10gb
        if cluster not in logfile.contents["final"]["progress"]["size_check"].keys():
            print(f"Cluster {cluster}: Checking size of total libarary...\n")
            sizes_to_truncate , cap_size = misc.get_truncate_sizes(accessions, C_fastqdir, 10737418240 ) #10gb total lib limit per cluster
            if sizes_to_truncate == "NA":  
                logfile.contents["final"]["progress"]["size_check"][cluster] = "pass"
            else:
                print(f"Cluster {cluster}: Truncating accession read files to fit library limit...\n")
                for acc , size in zip(accessions, sizes_to_truncate):
                    sourcepath= os.path.join(C_fastqdir,f"{acc}.fastq.gz")
                    targetpath= os.path.join(cluster_tempdir,f"{acc}.fastq.gz")
                    os.system(f"cp {sourcepath} {targetpath} && truncate -s {size} {targetpath}")
                logfile.contents["final"]["progress"]["size_check"][cluster] = "truncated"
            logfile.update()
        print(f"Cluster {cluster}: Library size check completed.\n")
            
        #parallel adapter trimming by fastp
        if cluster not in logfile.contents["final"]["progress"]["fastp"].keys():
            if logfile.contents["final"]["progress"]["size_check"].get(cluster) == "truncated":
                parallel_fastp(threadpool, accessions, cluster_tempdir, cluster_fastqdir, cluster)
                os.system(f"rm -r {cluster_tempdir}")
                print(f"Cluster {cluster}: concatenating trimmed reads into single library...\n")
                os.system(f"cat {cluster_fastqdir}/*.fastq > " + os.path.join(cluster_fastqdir, "concat.fq"))
                logfile.contents["final"]["progress"]["fastp"][cluster] = "done"
            
            if logfile.contents["final"]["progress"]["size_check"].get(cluster) == "pass":
                parallel_fastp(threadpool, accessions, C_fastqdir , cluster_fastqdir, cluster)
                os.system(f"rm -r {cluster_tempdir}")
                print(f"Cluster {cluster}: concatenating trimmed reads into single library...\n")
                os.system(f"cat {cluster_fastqdir}/*.fastq > " + os.path.join(cluster_fastqdir, "concat.fq"))
                logfile.contents["final"]["progress"]["fastp"][cluster] = "done"
            logfile.update()    
        print(f"\nCluster {cluster}: Fastp trimming completed.\n")
        if cluster not in logfile.contents["final"]["progress"]["ORNA"].keys():
            print(f"Cluster {cluster}: Normalising reads using ORNA....\nNOTE: This will take a while so this step will be verbose\n\n\n")
            sleep(5)
            if threadpool == 16:
                #for some reason, ORNA cannot multi-thread at 16 threads. bump it down to 15 instead
                trim.launch_ORNA(os.path.join(cluster_fastqdir, "concat.fq"), os.path.join(cluster_fastqdir, "concat_norm"), 15)
            else:
                trim.launch_ORNA(os.path.join(cluster_fastqdir, "concat.fq"), os.path.join(cluster_fastqdir, "concat_norm"), threadpool)
            logfile.contents["final"]["progress"]["ORNA"][cluster] = "done"
            logfile.update()
        print(f"\nCluster {cluster}: Read normalization completed.\n")
        
        if cluster not in logfile.contents["final"]["progress"]["assembly"].keys():
            soapdenovo.make_config(os.path.join(cluster_fastqdir, "concat_norm.fq"), os.path.join(cluster_assemblydir,"raw" ,"configfile.config"))
            logfile.contents["final"]["progress"]["assembly"][cluster] = {}
            logfile.update()  
        
        #loop for each kemer in kmer range    
        for kmer in range(kmer_min, kmer_max +1 , kmer_step):
            if str(kmer) not in logfile.contents["final"]["progress"]["assembly"][cluster].keys():
                #dnt assembly
                soapdenovo.launch_soap_verbose(os.path.join(cluster_assemblydir, "raw","configfile.config"), 
                kmer, 
                os.path.join(cluster_assemblydir, "raw", f"c{cluster}k{str(kmer)}"),
                threadpool)
                #change name of dnt output
                os.system(f"mv " + os.path.join(cluster_assemblydir, "raw" ,f"c{cluster}k{str(kmer)}_temp.scafSeq") + " " 
                + os.path.join(cluster_assemblydir,"raw" ,f"c{cluster}k{str(kmer)}.fasta"))
                #delete intermediate files
                os.system(f"rm "+os.path.join(cluster_assemblydir, "raw","*temp*"))
                #extract ORF
                soapdenovo.extract_orf(os.path.join(cluster_assemblydir, "raw" ,f"c{cluster}k{str(kmer)}.fasta"),
                os.path.join(cluster_assemblydir, "raw" ,f"c{cluster}k{str(kmer)}_cds.fasta"),
                orfminlen,
                startcodon,
                geneticcode)
                
                logfile.contents["final"]["progress"]["assembly"][cluster][str(kmer)]= "done"
                logfile.update()
                print(f"\nCluster {cluster}: CDS assembly at {kmer} k-mer length completed.\n")
        print(f"Cluster {cluster}: Assembly completed at all specified k-mer lengths.\n")
        
        #combine using CD hit retain rep of CDS clusters with hits from more than 1 kmer length
        if "combine" not in logfile.contents["final"]["progress"]["assembly"][cluster].keys():
            print(f"Cluster {cluster}: Combining all k-mer CDSs into a single multi-kmer CDS assembly...\n")
            #concat and rename 
            postprocess.concat_rename_assemblies(os.path.join(cluster_assemblydir, "raw"), os.path.join(cluster_assemblydir, "combined", f"c{cluster}_mk_red_cds.fasta"))
            #launch CDhit
            postprocess.launch_cdhit(os.path.join(cluster_assemblydir, "combined", f"c{cluster}_mk_red_cds.fasta"), 
            0.98,
            os.path.join(cluster_assemblydir, "combined", f"c{cluster}_mk_nr_cds.fasta"),
            threadpool)
            #extract Clstr file
            seqtoretain = postprocess.cluster_seq_extractor(1,os.path.join(cluster_assemblydir, "combined", f"c{cluster}_mk_nr_cds.fasta.clstr"))
            #subset nr fasta from CDHIT output
            postprocess.fasta_subset(os.path.join(cluster_assemblydir, "combined", f"c{cluster}_mk_nr_cds.fasta"),
            os.path.join(cluster_assemblydir, "combined", f"c{cluster}_mk_cds.fasta"),
            seqtoretain)
            n_cds, avg_cds_len, GC = misc.get_assembly_stats(os.path.join(cluster_assemblydir, "combined", f"c{cluster}_mk_cds.fasta"))
            logfile.contents["final"]["progress"]["assembly"][cluster]["combine"] = [n_cds, avg_cds_len, GC]
            logfile.update()
        print(f"Cluster {cluster}: Multi-kmer CDS assembly generated.\n")
        
        if cluster not in logfile.contents["final"]["progress"]["CPC2"].keys():
            print(f"Cluster {cluster}: Extracting CDS with coding potential using CPC2...\n")
            
            postprocess.launch_CPC2(os.path.join(cluster_assemblydir, "combined", f"c{cluster}_mk_cds.fasta"),
            os.path.join(cluster_assemblydir, "CPC2", "CPC2_output"))
            
            cds_list , ncds_list = postprocess.parse_CPC2_output(os.path.join(cluster_assemblydir, "CPC2", "CPC2_output.txt"))
            
            postprocess.fasta_subset(os.path.join(cluster_assemblydir, "combined", f"c{cluster}_mk_cds.fasta"),
            os.path.join(cluster_assemblydir, "CPC2", f"c{cluster}_CPC2_cds.fasta"),
            cds_list)
            
            n_cds, avg_cds_len, GC = misc.get_assembly_stats(os.path.join(cluster_assemblydir, "CPC2", f"c{cluster}_CPC2_cds.fasta"))
            
            logfile.contents["final"]["progress"]["CPC2"][cluster] = [n_cds, avg_cds_len, GC]
            logfile.update()
        print(f"Cluster {cluster}: CDS with coding potential extracted.\n")
        print(f"Cluster {cluster}: CDS-mining partially complete.\n")
        #if cluster not in logfile.contents["final"]["progress"]["CPC2"].keys():
            #print(f"Cluster {cluster}: Remapping reads from ")
        
        
        
        #print(f"Cluster {cluster}: CDS-mining complete.\n")
