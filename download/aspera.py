#setting sys.path for importing modules
import os
import sys
if __name__ == "__main__":
        abspath= os.getcwd()
        parent_module= os.path.join(abspath.split("LSTrAP-denovo")[0], "LSTrAP-denovo")
        sys.path.insert(0, parent_module)
 

import subprocess


from setup import constants

#returns paths of ftp and ascp directories where reads for a particular SRR accession can be found
def get_download_path(accession):
    if len(accession) ==9:
        path =f"/vol1/fastq/{accession[:6]}/{accession}/"
    elif len(accession) ==10:
        path =f"/vol1/fastq/{accession[:6]}/00{accession[-1]}/{accession}/"
    elif len(accession) == 11:
        path= f"/vol1/fastq/{accession[:6]}/0{accession[-2:]}/{accession}/"
    elif len(accession) == 12:
        path= f"/vol1/fastq/{accession[:6]}/{accession[-3:]}/{accession}/"
    http_dirpath=  "http://ftp.sra.ebi.ac.uk" + path
    http_fullpath_1 = os.path.join(http_dirpath,accession+".fastq.gz")
    http_fullpath_2 = os.path.join(http_dirpath,accession+"_1.fastq.gz")
    filesize_1 = check_filesize(http_fullpath_1)
    filesize_2 = check_filesize(http_fullpath_2)
    if  filesize_1 > filesize_2:
        filesize= filesize_1
        ascp_fullpath=  os.path.join("era-fasp@fasp.sra.ebi.ac.uk:" +path, accession+".fastq.gz")
        ftp_fullpath= os.path.join("ftp://ftp.sra.ebi.ac.uk" + path, accession+".fastq.gz")
    elif filesize_1 < filesize_2:
        filesize= filesize_2
        ascp_fullpath=  os.path.join("era-fasp@fasp.sra.ebi.ac.uk:" +path, accession+"_1.fastq.gz")
        ftp_fullpath= os.path.join("ftp://ftp.sra.ebi.ac.uk" + path, accession+"_1.fastq.gz")
    else:
        filesize= filesize_1
        ascp_fullpath=  os.path.join("era-fasp@fasp.sra.ebi.ac.uk:" +path, accession+".fastq.gz")
        ftp_fullpath= os.path.join("ftp://ftp.sra.ebi.ac.uk" + path, accession+".fastq.gz")
    return ascp_fullpath, ftp_fullpath, filesize

def check_filesize(http_fullpath):
    header= subprocess.check_output(f"curl -sI {http_fullpath}", shell=True).decode("utf-8")
    if "404 Not Found" in header:
        return 0
    elif "200 OK" in header:
        filesize= int(header.split("Content-Length: ")[1].split("\r")[0])
        return filesize
        
def launch_ascp(ascp_fullpath, outputdir, filesizelimit=1500000000):
    subprocess.run([constants.asperapath, "-QT", "-l", "300m", "-P33001", "-@",f"0:{filesizelimit}" ,"-i", constants.aspera_ssh_key , ascp_fullpath, outputdir], 
    stdout=subprocess.DEVNULL, 
    stderr=subprocess.STDOUT)
    
def launch_curl(ftp_fullpath, outputdir, filesizelimit=1500000000):
    subprocess.run(["curl", "-r", f"0-{str(filesizelimit)}","-o", outputdir, ftp_fullpath ])
    #stdout=subprocess.DEVNULL,
    #stderr=subprocess.STDOUT)
__all__=["get_download_path", "check_filesize", "launch_ascp","launch_curl" ]