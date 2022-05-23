#setting sys.path for importing modules
import os
import sys
if __name__ == "__main__":
        abspath= os.getcwd()
        parent_module= os.path.join(abspath.split("LSTrAP-denovo")[0], "LSTrAP-denovo")
        sys.path.insert(0, parent_module)
        
import subprocess
from setup import constants
import shutil
import pandas as pd

def launch_kallisto_index(fastapath,indexpath):
    """Contruct a kallisto index using assemnbly fasta"""
    returncode=subprocess.run([constants.kallistopath, "index", "-i", indexpath , fastapath],
    stdout=subprocess.DEVNULL, 
    stderr=subprocess.STDOUT)
    return returncode.returncode

def launch_kallisto_quant(threads,indexpath, outputdir, fastqpath):
    """ launch kallisto quant in single-ended mode """
    returncode=subprocess.run([constants.kallistopath, "quant", "-t", str(threads), "-i", indexpath, "-o", outputdir, "--single", "-l", "200", "-s" , "20" , fastqpath ],
    stdout=subprocess.DEVNULL, 
    stderr=subprocess.STDOUT) 
    return returncode.returncode

def write_quant_info(accession ,outdirpath, tpm_matpath):
    """extract data from quant output and write to tpm matrix. Returns maprate (pseudoalignment %)
    NOTE: it wil auto-remove kallisto output dir when done."""
    target_id = list(pd.read_csv(os.path.join(outdirpath,"abundance.tsv"), sep="\t")["target_id"])
    tpm_data = list(pd.read_csv(os.path.join(outdirpath,"abundance.tsv"), sep="\t")["tpm"])
    map_rate = open(os.path.join(outdirpath,"run_info.json")).read().split("\"p_pseudoaligned\": ")[1].split(",")[0]
    
    #if file exist, create with header before writing info. Else, append info to file
    if not os.path.exists(tpm_matpath):
        with open(tpm_matpath, "w") as f:
            f.write("accession\t"+"\t".join(target_id)+"\n")
            f.write(f"{accession}\t"+"\t".join([str(tpm) for tpm in tpm_data])+ "\n")
    else:
        with open(tpm_matpath, "a") as f:
            f.write(f"{accession}\t"+"\t".join([str(tpm) for tpm in tpm_data])+ "\n")
    
    #shutil.rmtree(outdirpath) ###might include this to remove kallisto output directory
    return map_rate


__all__=["launch_kallisto_index", "launch_kallisto_quant", "write_quant_info"]