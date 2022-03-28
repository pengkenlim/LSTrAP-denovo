#setting sys.path for importing modules
import os
import sys
if __name__ == "__main__":
        abspath= os.getcwd()
        parent_module= os.path.join(abspath.split("LSTrAP-denovo")[0], "LSTrAP-denovo")
        sys.path.insert(0, parent_module)


import subprocess

from setup import constants

def make_config(fastqpath,configoutpath ):
    with open(configoutpath, "w") as f:
        f.write(f"max_rd_len=100\n[LIB]\nrd_len_cutof=100\navg_ins=200\nreverse_seq=0\nasm_flags=3\nmap_len=32\nq={fastqpath}")

def extract_orf(fastainputpath, fastaoutputpath, orfminlen, startcodon, geneticcode):
    returncode= subprocess.run([constants.orffinderpath, "-in", fastainputpath, "-out", fastaoutputpath, "-ml", str(orfminlen), "-s", str(startcodon), "-g", str(geneticcode)],
    stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
    return returncode.returncode
    
    
def launch_soap(configoutpath, kmerlen, outputpath_prefix, threads):
    returncode=subprocess.run([constants.soappath, "all", "-s", configoutpath, "-o", outputpath_prefix+"_temp", "-K", str(kmerlen), "-p", str(threads)],
    stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
    #renaming transcript file to keep and deleting others
    os.system(f"mv {outputpath_prefix}_temp.scafSeq {outputpath_prefix}.fasta")
    os.system(f"rm -r {outputpath_prefix}_temp*")
    return returncode.returncode

