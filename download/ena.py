#setting sys.path for importing modules
import os
import sys
if __name__ == "__main__":
        abspath= os.getcwd()
        parent_module= os.path.join(abspath.split("LSTrAP-denovo")[0], "LSTrAP-denovo")
        sys.path.insert(0, parent_module)

import requests
import json

def get_runs(taxid):
    req= requests.get(f"https://www.ebi.ac.uk/ena/portal/api/search?result=read_run&query=tax_tree({taxid})%20AND%20library_source=%22TRANSCRIPTOMIC%22&format=json")
    run_list = [k["run_accession"]for k in req.json()]
    return run_list

def get_sciname(taxid):
    req= requests.get(f"https://www.ebi.ac.uk/ena/portal/api/search?result=taxon&query=tax_eq({taxid})&fields=scientific_name&format=json")
    try:
        return req.json()
    except:
        return None

__all__ = ['get_runs', 'get_sciname']
