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
import numpy as np

def launch_kallisto_index(fastapath,indexpath):
    """Contruct a kallisto index using assemnbly fasta"""
    returncode=subprocess.run([constants.kallistopath, "index", "-i", indexpath , fastapath],
    stdout=subprocess.DEVNULL, 
    stderr=subprocess.STDOUT)
    return returncode.returncode

def launch_kallisto_quant(threads,indexpath, outputdir, fastqpath):
    """ launch kallisto quant in single-ended mode """
    try:
        returncode=subprocess.run([constants.kallistopath, "quant", "-t", str(threads), "-i", indexpath, "-o", outputdir, "--single", "-l", "200", "-s" , "20" , fastqpath ],
        stdout=subprocess.DEVNULL, 
        stderr=subprocess.STDOUT)
    except:
        return 1
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
    return map_rate


def re_mapping(assemblydir, cluster, cluster_list, threads, basedir, rankingthreshold): #not in use
    '''outputdir = assemblydir (i.e. final/cluster_{cluster}/assembly/) . Cluster the assembly to map against. Cluster_list holds the list of clusters
    basdir = .../final/'''
    if not os.path.exists(os.path.join(assemblydir, "remap")):
        os.makedirs(os.path.join(assemblydir, "remap"))
    #build index    
    launch_kallisto_index(os.path.join(assemblydir, "CPC2", f"c{cluster}_CPC2_cds.fasta"), os.path.join(assemblydir, "remap", f"c{cluster}_CPC2_cds.index"))
    #map every single concat read from each accession cluster to the index
    
    for i in cluster_list:
        launch_kallisto_quant(threads, os.path.join(assemblydir, "remap", f"c{cluster}_CPC2_cds.index") , os.path.join(assemblydir, "remap", f"AC{i}"), os.path.join(basedir, f"cluster_{i}" ,"fastq", "concat.fq"))
    tpm_df = pd.DataFrame()
    for i in cluster_list:
        tpm_data= pd.read_csv(os.path.join(assemblydir, "remap", f"AC{i}", "abundance.tsv") , sep="\t").set_index("target_id")["tpm"]
        tpm_df[f"AC{i}"] = tpm_data   
    tpm_df.to_csv(os.path.join(assemblydir, "remap", "tpm_mat.tsv") , sep= "\t")
    print(f"Cluster {cluster}: Removing CDSs with poor mapping...\n ")
    seqtoretain=[]
    total=len(tpm_df)
    tpm_df = tpm_df[tpm_df[f"AC{cluster}"] != 0.0] #remove CDSs with TPM 0 when read lib of respective AC is aligned
    total2 =len(tpm_df)
    for column in list(tpm_df.columns):
        tpm_df.sort_values(by=[column], axis=0, ascending=False, inplace=True) #sort CDS from highest TPM to lowest for each liibrary
        sorted_genes= list(tpm_df.index)
        seqtoretain += sorted_genes[:int(np.round(len(sorted_genes)*((rankingthreshold)/100)))]
    seqtoretain = list(set(seqtoretain))
    total3= len(seqtoretain)
    return seqtoretain , total, total2, total3
    


__all__=["launch_kallisto_index", "launch_kallisto_quant", "write_quant_info"]