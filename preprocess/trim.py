#setting sys.path for importing modules
import os
import sys
if __name__ == "__main__":
        abspath= os.getcwd()
        parent_module= os.path.join(abspath.split("HSS-Trans")[0], "HSS-Trans")
        sys.path.insert(0, parent_module)
        

import subprocess

from setup import constants

def launch_fastp(inputpath,outputpath,threads):
    returncode=subprocess.run([constants.fastppath, "-w", str(threads), "--in1",  inputpath, "--out1", outputpath],
    stdout=subprocess.DEVNULL, 
    stderr=subprocess.STDOUT)
    return returncode.returncode

def launch_ORNA(inputpath, outputpathprefix, threads):
    returncode=subprocess.run([constants.ORNApath, "-type","fastq", "-input", inputpath, "-nb-cores", str(threads), "-kmer", "29", "-output",outputpathprefix])
    return returncode.returncode
    

    


__all__=["launch_fastp"]