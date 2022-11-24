#this is a throw away version of the code
import os
import sys
import numpy as np
import concurrent.futures
import subprocess
from datetime import datetime


#concat_fastapath="/shared/scratch/ken/pipeline_output/TAXID_4565/Step_2/selected_accessions/concat_renamed_wo_all_2.fasta"
concat_fastapath="/shared/scratch/ken/pipeline_output/TAXID_4565/Step_2/selected_accessions/AnnotatePredictORFs/concat_test.fasta"
tempdir="/shared/scratch/ken/pipeline_output/TAXID_4565/Step_2/selected_accessions/AnnotatePredictORFs"
pathtotransdecoderdir="~/LSTrAP-denovo/programs/TransDecoder"
pathtohmmsearch="hmmsearch"
pathtoPfamHMM= "~/my_interproscan/interproscan-5.59-91.0/data/pfam/Pfam-A.hmm"
workers= 1
worker_cpu= 4
genetic_code = "Universal"
min_prot_len= 100

def extract_ORFs(filepath,outdir):
    returncode=subprocess.run([os.path.join(pathtotransdecoderdir, "TransDecoder.LongOrfs"),
    "-t", filepath, "--outputdir", outdir, "--genetic_code", genetic_code, "-m", str(min_prot_len)]#,
    #stdout=subprocess.DEVNULL, 
    #stderr=subprocess.STDOUT
    )
    return returncode.returncode

def Pfam_hmmsearch(outdir, domtbloutpath):
    input_pep = os.path.join(outdir, "longest_orfs.pep") #path/to/splitfile_partx/longest_orfs.pep
    logpath = os.path.join(outdir, "PfamHMM.log") #path/to/splitfile_partx/PfamHMM.log
    os.system(f"{pathtohmmsearch} --cpu {worker_cpu} --domtblout {domtbloutpath} {pathtoPfamHMM} >{logpath}")
    

def run_job(file_name):
    print(f"{file_name}: Running Transdecoder.LongOrfs...")
    filepath = os.path.join(tempdir ,file_name) #path/to/splitfile_partx.fasta
    outdir = os.path.join(tempdir ,file_name.split(".fasta")[0]) #path/to/splitfile_partx
    returncode = extract_ORFs(filepath, outdir)
    print(f"{file_name}: Transdecoder.LongOrfs done. Running hmmsearch for Pfam domains...")
    domtbloutpath = os.path.join(outdir, "longest_orfs.domtblout")
    Pfam_hmmsearch(outdir, domtbloutpath)
    endtime = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
    return f"{file_name} completed at {endtime}"
    
if __name__ == "__main__":
    print(datetime.now().strftime("%d/%m/%Y %H:%M:%S"))   
    if not os.path.exists(tempdir):
        os.system(f"mkdir {tempdir}")

    #read content
    with open(concat_fastapath, "r") as f:
       contents = f.read()
       contents= contents.split(">")

    print(f"Total number of transcripts in fasta: {len(contents)}")
    n_seq_per_chunk = int((len(contents) - (len(contents)%workers))/workers) #<--- number of sequences in each split file


    #create each seq_chunk in tempdir
    file_names=[]
    for i in range(0,workers):
        if i == workers-1:
            towrite= contents[(i*n_seq_per_chunk):]
        else:    
            towrite= contents[(i*n_seq_per_chunk):((i+1)*n_seq_per_chunk)]
        towrite=">".join(towrite)
        with open(os.path.join(tempdir, f"splitfile_part{i+1}.fasta"), "w")as f:
            f.write(towrite)
        file_names += [f"splitfile_part{i+1}.fasta"]
        print(f"splitfile_part{i+1}.fasta created")
        

    with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as executor:
        results= [executor.submit(run_job, file_name) for file_name in file_names]
        for f in concurrent.futures.as_completed(results):
            print(f.result())

    print("script completed", datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
