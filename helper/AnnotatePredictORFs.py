#this is a throw away version of the code
import os
import sys
import numpy as np
import concurrent.futures
import subprocess

concat_fastapath="/shared/scratch/ken/pipeline_output/TAXID_4565/Step_2/selected_accessions/concat_renamed_wo_all_2.fasta"
tempdir="/shared/scratch/ken/pipeline_output/TAXID_4565/Step_2/selected_accessions/AnnotatePredictORFs"
pathtotransdecoderdir="~/LSTrAP-denovo/programs/TransDecoder"
pathtohmmsearch="hmmsearch"
pathtoPfamHMM= "~/my_interproscan/interproscan-5.59-91.0/data/pfam/Pfam-A.hmm"
workers= 10
worker_cpu= 4
genetic_code = "Universal"
min_prot_len= 100

def extract_ORFs(filepath,outpathdir):
    returncode=subprocess.run([os.path.join(pathtotransdecoder, "TransDecoder.LongOrfs"),
    "-t", filepath, "--outputdir", outdir, "--genetic_code", genetic_code, "-m", str(min_prot_len)],
    stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
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
    
    
if not os.path.exists(tempdir):
    os.system(f"mkdir {tempdir}")

#read content
with open(fastapath) as f:
   contents = f.read()
   contents= contents.split(">")

n_seq_per_chunk = int((len(contents) - (len(contents)%workers))/workers) #<--- number of sequences in each split file

#create each seq_chunk in tempdir
file_names=[]
for i in range(0,workers):
    towrite= contents[(i*n_seq_per_chunk):((i+1)*n_seq_per_chunk)]
    towrite=">".join(towite)
    with open(os.path.join(tempdir, f"splitfile_part{i}.fasta"))as f:
        f.write(towrite)
    file_names += ["splitfile_part{i}.fasta"]
    

with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as executor:
    pass