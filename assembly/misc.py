#setting sys.path for importing modules
import os
import sys
if __name__ == "__main__":
        abspath= os.getcwd()
        parent_module= os.path.join(abspath.split("HSS-Trans")[0], "HSS-Trans")
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
    return ""

class logfile:
    ''' logfile object to store, retrieve and update checkpoint information'''
    def __init__(self, path):
        self.path= path
        self.template={"Step_1":{"run_info":{"taxid":None, "sci_name":None, "n_total_acc":None, "command_issued": None, "init_time": None}, 
        "run_var":None ,
        "total_acc": None,
        "processed_acc":None, 
        "consensus":{"stats":None, "optimal":None}, 
        "status":"incomplete"},
        "Step_2":{"run_info":{"command_issued":None, "init_time":None},
        "run_var":None,
        "processed_acc": None,
        "qc":{"threshold":None, "passed":None, "failed":None, "total":None},
        "kmeans":{"s_coeficient":None, "cluster_assignment_dict":None},
        "selected_accessions": {},
        "status": "incomplete"}}

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
        if step == "step_1":
            self.contents["step_1"]=self.template["step_1"]
        elif step == "step_2":
            self.contents["step_2"]= self.template["step_2"]
        with open(path, "w") as f:
            json.dump(self.contents,f, indent=2)


def get_assembly_stats(pathtofasta):
    with open(pathtofasta, "r") as f:
        fasta= f.read()
    n_cds=fasta.count(">")
    fasta_lines= fasta.split("\n")
    seq_concat= "".join([k  for k in fasta_lines if ">" not in k])
    if n_cds !=0:
        avg_cds_len = int(len(seq_concat)/n_cds)
        GC= np.round((seq_concat.count("G") + seq_concat.count("C"))/len(seq_concat)*100,2)
    else:
        avg_cds_len= "NA"
        GC= "NA"
    return n_cds, avg_cds_len, GC
        

def get_truncate_sizes(accessions, fastqdir ,lib_limit):
    """get sizes to truncate files based on target lib size
    Returns array of target sizes to truncate"""
    sizes = [os.path.getsize(os.path.join(fastqdir, f"{acc}.fastq.gz")) for acc in accessions]
    cap = max(sizes)
    if sum(sizes) < lib_limit:
        return "NA" , "NA"
    while sum(sizes) > lib_limit:
        cap+= -10485760 #10mb
        sizes= [min([size, cap]) for size in sizes]
        
    return sizes, cap
    
class logfile_expmat:
    ''' logfile object to store, retrieve and update checkpoint information'''
    def __init__(self, path):
        self.path= path
        self.template={"run_info":{"taxid":None, "sci_name":None, "n_total_acc":None, "command_issued": None, "init_time": None}, 
        "run_var":None ,
        "total_acc": None,
        "processed_acc":None,
        "status":"incomplete"}

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
    
    def clear(self): # reset log file
        self.contents = self.template
        path= self.path
        with open(path, "w") as f:
            json.dump(self.contents,f, indent=2)
