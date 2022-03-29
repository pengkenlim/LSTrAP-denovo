#setting sys.path for importing modules
import os
import sys
if __name__ == "__main__":
        abspath= os.getcwd()
        parent_module= os.path.join(abspath.split("LSTrAP-denovo")[0], "LSTrAP-denovo")
        sys.path.insert(0, parent_module)
        

import json

#basiclly a wrapper function to re-run shell-command functions that failed (based on subprocess.run()'s returncode)       
def run_with_retries(retry_limit,func,arg_list,retry_message,run_message):
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

#log file to read and write progress
class logfile:
    def __init__(self, path):
        self.path= path
        if not os.path.exists(path):
            with open(path, "w") as f:
                empty_dict={"prelim":{"cmd_args":{},"processed":[],"status":"incomplete"}}
                json.dump(empty_dict,f, indent=1)
                self.contents = empty_dict
        else:
            
            with open(path, "r") as f:
               self.contents= json.load(f)
    def update(self):
        path= self.path
        with open(path, "w") as f:
            json.dump(self.contents,f, indent=1)

