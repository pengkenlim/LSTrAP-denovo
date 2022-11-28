#Replace with 
import os
import sys
import numpy as np
import concurrent.futures
import subprocess
from datetime import datetime
import argparse


def extract_ORFs(filepath,dirname):
    return_code= os.system(f"cd {working_dir}; " + os.path.join(transdecoder_bin_dir, "TransDecoder.LongOrfs") + " " + 
    f"-t {filepath} --output_dir {dirname} -G {genetic_code} -m {str(min_prot_len)} > {working_dir}/{dirname}.logs")
    return return_code
    
def predict_ORFs(filepath,dirname, domtbloutpath):
    return_code= os.system(f"cd {working_dir}; " + os.path.join(transdecoder_bin_dir, "TransDecoder.Predict") + " " + 
    f"-t {filepath} --output_dir {dirname} -G {genetic_code} --retain_pfam_hits {domtbloutpath} >> {working_dir}/{dirname}.logs")
    return return_code


def Pfam_hmmsearch(outdir, domtbloutpath):
    input_pep = os.path.join(outdir, "longest_orfs.pep") #path/to/splitfile_partx/longest_orfs.pep
    logpath = os.path.join(outdir, "PfamHMM.log") #path/to/splitfile_partx/PfamHMM.log
    return_code= os.system(f"{hmmsearch_bin} --cpu {worker_cpu} --domtblout {domtbloutpath} {pathtoPfamHMM} {input_pep} >{logpath}")
    return return_code
    

def run_job(file_name):
    print(f"{file_name}: Running Transdecoder.LongOrfs...")
    filepath = os.path.join(working_dir ,file_name) #path/to/splitfile_partx.fasta
    outdir = os.path.join(working_dir ,file_name.split(".fasta")[0]) #path/to/splitfile_partx
    
    #extract ORF
    return_code = extract_ORFs(filepath, file_name.split(".fasta")[0])
    if return_code ==0:
        print(f"{file_name}: Transdecoder.LongOrfs done. Running hmmsearch for Pfam domains...")
    else:
        return f"{file_name}: ERROR. Transdecoder.LongOrfs Failed."
    
    domtbloutpath = os.path.join(outdir, "longest_orfs.domtblout")
    #Pfam hmmsearch
    return_code= Pfam_hmmsearch(outdir, domtbloutpath)
    if return_code ==0:
        print(f"{file_name}: hmmsearch done. Running Transdecoder.Predict...")
    else:
        return f"{file_name}: ERROR. hmmsearch Failed."
    
    #flipping targets and queries in domtblouts output by hmmsearch so that it is consistent with hmmscan
    flipped_domtbloutpath = os.path.join(outdir, "longest_orfs_flipped.domtblout")
    os.system("awk \'BEGIN{OFS=FS=\" \"} NR<=3{print}; NR>3{tmp=$1; $1=$4; $4=tmp; tmp=$2; $2=$5; $5=tmp; print}\'" + f" {domtbloutpath} > {flipped_domtbloutpath}")
    
    #predcit ORFs hmmsearch
    return_code= predict_ORFs(filepath,dirname, flipped_domtbloutpath)
    endtime = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
    if return_code ==0:
        return f"{file_name} completed at {endtime}"
    else:
        return f"{file_name}: ERROR. Transdecoder.Predict Failed."
    
    
if __name__ == "__main__":
    #specify arguments
    parser= argparse.ArgumentParser(description="AnnotatePredictORFs.py: Helper script that uses Transdecoder's ORF prediction to extract Coding sequences (CDS) and annotate protein function based on Pfam domains.\n \
    NOTE: This is to be run following completion of the HSS-Trans pipeline and Trinity transcriptome assembly (Run_Trinity).\n\
    Please make sure binaries for TransDecoder, hmmsearch, hmmpress and its associated dependencies are installed and added into $PATH.\n\
    \n\
    Refer to https://github.com/pengkenlim/HSS-Trans for more information on pipeline usage and implmentation.\n\
    Refer to ... for more information on Transdecoder.\n\
    Refer to ... for more information on hmmsearch.\n\
    Refer to ... for more information on hmmpress.")
    
    parser.add_argument("-o", "--output_dir", type=str, metavar= "", required=True,
    help= "Directory for data output. Directory needs to be same as for step 1 and step 2 of HSS-Trans pipeline (i.e. MakeDraftCDS.py, SelectAccessions.py).")
    
   
    parser.add_argument("-t", "--threads", type=int, metavar="", default=4, 
    help = "Total thread pool for workers. Needs to be divisible by number of workers.")
    
    parser.add_argument("-w", "--workers", type=int, metavar="", default=2, 
    help= "Specifies the maximum workers for running Transdecoder and hmmsearch in parallel. Reccomended to set to (thread pool)/3. Set to 2 by default.")
    
    parser.add_argument("-G", "--genetic_code", type=str, metavar="", default="Universal", 
    help= "Genetic code passed to TransDecoder.LongOrfs and Transdecoder.Predict. Set to \"Universal\" by default. Refer to ... for more information.")
    
    parser.add_argument("-m", "--min_prot_len", type=int, metavar="", default=100,
    help= "Minimal protein length (aa) passed to TransDecoder.LongOrfs. Set to 100 by default.")
    
    parser.add_argument("-pdir","--pfam_dir", type=str, metavar="", required=True,
    help="Path to directory containing Pfam hmm model (Pfam-A.hmm). If directory specified do not contain Pfam-A.hmm, it will be downloaded auotmatically." )
    
    parser.add_argument("-hb", "--hmmsearch_bin", type=str, metavar= "", default="hmmsearch",
    help= "Path to hmmsearch binary. Not required if hmmsearch binary directory has been added to $PATH.")
    
    parser.add_argument("-hpb", "--hmmpress_bin", type=str, metavar= "", default="hmmpress",
    help= "Path to hmmpress binary. Not required if hmmpress binary directory has been added to $PATH.")
    
    parser.add_argument("-tbdir", "--transdecoder_bin_dir", type=str, metavar= "", default="",
    help= "Path to directory containing transdecoder binaries (Transdecoder.LongORFs and Transdecoder.Predict). Not required if hmmsearch binary directory has been added to $PATH.")
    
    #parse args
    args=parser.parse_args()

    output_dir= args.output_dir
    threadpool =args.threads
    workers =args.workers
    genetic_code = args.genetic_code
    min_prot_len= args.min_prot_len
    pfam_dir= args.pfam_dir
    hmmsearch_bin = args.hmmsearch_bin
    hmmpress_bin = args.hmmpress_bin
    transdecoder_bin_dir = args.transdecoder_bin_dir
    

    #check for trinity_dir, findout name of assembly fasta. set path for assembly fasta.
    Trinity_dir= os.path.join(output_dir, "Trinity_output")
    if not os.path.exists(Trinity_dir):
         sys.exit("Error in output directory. Trinity directory not found. Exiting...")
    if "Trinity_over-assembly_nr.fasta" not in os.listdir(Trinity_dir) and "Trinity_all_samples.fasta" not in os.listdir(Trinity_dir):
        sys.exit("Error in output directory. Trinity transcriptome assembly fasta files not found. Exiting...")
    if "Trinity_over-assembly_nr.fasta" in os.listdir(Trinity_dir):
        fastapath = os.path.join(Trinity_dir, "Trinity_over-assembly_nr.fasta")
    else:
        fastapath = os.path.join(Trinity_dir, "Trinity_all_samples.fasta")
    print(f"Assembly file detected in {fastapath}.")
    
    #check and set threads. Check if threadpool can be divided by number of workers. Else set to appropriate threads
    if threadpool%workers != 0:
        print(f"Threadpool of {threadpool} specified by user is not divisible by {workers} workers. Threadpool corrected to {int((threadpool - (threadpool%workers)))}.")
        threadpool = int((threadpool - (threadpool%workers)))
    
    threads=int(threadpool/workers)
    
    #check if path to hmmsearch, hmmpress and transdecoder binaries are valid.
    if hmmsearch_bin!= "hmmsearch":
        if not os.path.exists(hmmsearch_bin):
            sys.exit(f"Error: hmmsearch not found at {hmmsearch_bin}. Exiting...")
    
    if hmmpress_bin!= "hmmpress":
        if not os.path.exists(hmmpress_bin):
            sys.exit(f"Error: hmmpress_bin not found at {hmmpress_bin}. Exiting...")
    
    if transdecoder_bin_dir!= "":
        if not os.path.exists(os.path.join(transdecoder_bin_dir, "TransDecoder.LongOrfs")):
            sys.exit(f"Error: TransDecoder.LongOrfs not found in {transdecoder_bin_dir}. Exiting...")
        if not os.path.exists(os.path.join(transdecoder_bin_dir, "TransDecoder.Predict")):
            sys.exit(f"Error: TransDecoder.Predict not found in {transdecoder_bin_dir}. Exiting...")

    #check if Pfam hmm is downloaded in directory

    if not os.path.exists(pfam_dir):
        print(f"User specified pfam directory ({pfam_dir}) is created as it cannot be found in the system.")
        os.makedirs(pfam_dir)
    pathtoPfamHMM =os.path.join(pfam_dir, "Pfam-A.hmm")
    if not os.path.exists(pathtoPfamHMM):
        print(f"Pfam-A.hmm not detected in user specifed pfam directory. Attempting to download from ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam35.0/Pfam-A.hmm.gz...")
        return_code = os.system(f"wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam35.0/Pfam-A.hmm.gz --directory-prefix={pfam_dir}")
        if return_code != 0:
            sys.exit("Error in downloading Pfam-A.hmm. Please download and decompress manaully into Pfam directory . Exiting...")
        print("un-gunzipping Pfam-A.hmm ... ")
        return_code = os.system("gunzip " + os.path.join(pfam_dir,"Pfam-A.hmm.gz"))
        if return_code != 0:
            sys.exit("Error in un-gunzipping Pfam-A.hmm. Please download and decompress manaully into Pfam directory . Exiting...")
    if not os.path.exists(pathtoPfamHMM):
        sys.exit("Error in un-gunzipping Pfam-A.hmm. Please download and decompress manaully into Pfam directory . Exiting...")
    
    print(f"Generating hmm database from Pfam-A.hmm using hmmpress...")
    return_code = os.system(f"{hmmpress_bin} {pathtoPfamHMM}")
    if return_code != 0:
        sys.exit("Error in running hmmpress on Pfam-A.hmm. Exiting...")
    
    #make working directory if not already exists
    working_dir= os.path.join(output_dir, "AnnotatePredictORFs")
    if not os.path.exists(working_dir):
        os.makedirs(working_dir)
    
    #read content of fasta file
    with open(fastapath, "r") as f:
       contents = f.read()
       contents= contents.split(">")
    
    print(f"Total number of transcripts in assembly fasta: {len(contents)}")
    n_seq_per_file = int((len(contents) - (len(contents)%workers))/workers)
    
    print(f"Number of specified workers= {workers} \n\
    {len(contents)}Sequences in assembly fasta file will be split into {workers}splitfiles containing approx. {n_seq_per_file} sequences each.\n")
    
    #create each seq_chunk in working_dir
    file_names=[]
    for i in range(0,workers):
        if i == workers-1:
            towrite= contents[(i*n_seq_per_file):]
        else:    
            towrite= contents[(i*n_seq_per_file):((i+1)*n_seq_per_file)]
        towrite=">".join(towrite)
        with open(os.path.join(working_dir, f"splitfile_part{i+1}.fasta"), "w")as f:
            f.write(towrite)
        file_names += [f"splitfile_part{i+1}.fasta"]
        print(f"splitfile_part{i+1}.fasta created")
    
    with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as executor:
        results= [executor.submit(run_job, file_name) for file_name in file_names]
        for f in concurrent.futures.as_completed(results):
            print(f.result())

    print("script completed", datetime.now().strftime("%d/%m/%Y %H:%M:%S"))

