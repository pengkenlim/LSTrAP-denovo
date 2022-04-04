#setting sys.path for importing modules
import os
import sys
if __name__ == "__main__":
        abspath= os.getcwd()
        parent_module= os.path.join(abspath.split("LSTrAP-denovo")[0], "LSTrAP-denovo")
        sys.path.insert(0, parent_module)
        

import json
import numpy as np    
def run_with_retries(retry_limit,func,arg_list,retry_message,run_message):
    """wrapper function to re-run shell-command functions that failed (based on subprocess.run()'s returncode)"""
    retries=0
    while True:
        if retries!=0:
            print(retry_message)
        print(run_message)
        returncode=func(*arg_list)
        if returncode==0:
            break
        else:
            retries+=1
        if retries > retry_limit:
            return "failed"        
    return "succeeded"

class logfile:
    ''' logfile object to store, retrieve and update checkpoint information'''
    def __init__(self, path):
        self.path= path
        self.template={"prelim":{"run_info":{"taxid":None, "sci_name":None, "n_total_acc":None, "command_issued": None, "init_time": None}, 
        "run_var":None ,"total_acc": None,"processed_acc":None, "consensus":{"CDS":None,"stats":{"CT":None, "path":None, "n_CDS": None ,"CDS_len":None, "GC": None }}, "status":"incomplete"}, "cluster":None, "modular_assembly":None}
        if not os.path.exists(path):
            with open(path, "w") as f:
                json.dump(self.template,f, indent=2)
                self.contents = self.template
        else:   
            with open(path, "r") as f:
               self.contents= json.load(f)
    
    def update(self):
        path= self.path
        with open(path, "w") as f:
            json.dump(self.contents,f, indent=2)
    
    def load(self):
        path= self.path
        with open(path, "r") as f:
            self.contents= json.load(f)
    
    def clear(self,step):
        path= self.path
        if step == "prelim":
            self.contents["prelim"]=self.template["prelim"]
        elif step == "cluster":
            self.contents["cluster"]= self.template["cluster"]
        elif step == "modular_assembly":
            self.contents["modular_assembly"]= self.template["modular_assembly"]
        with open(path, "w") as f:
            json.dump(self.contents,f, indent=2)


def get_assembly_stats(pathtofasta):
    with open(pathtofasta, "r") as f:
        fasta= f.read()
    n_cds=fasta.count(">")
    fasta_lines= fasta.split("\n")
    seq_concat= "".join([k  for k in fasta_lines if ">" not in k])
    avg_cds_len = int(len(seq_concat)/n_cds)
    GC= np.round((seq_concat.count("G") + seq_concat.count("C"))/len(seq_concat)*100,2)
    
    return n_cds, avg_cds_len, GC
        
    