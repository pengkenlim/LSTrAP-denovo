#setting sys.path for importing modules
import os
import sys
if __name__ == "__main__":
        abspath= os.getcwd()
        parent_module= os.path.join(abspath.split("HSS-Trans")[0], "HSS-Trans")
        sys.path.insert(0, parent_module)
 

import subprocess
import json

from setup import constants

#returns paths of ftp and ascp directories where reads for a particular SRR accession can be found
def get_download_path(accession):
    """deprecated, replaced by get_download_path_ffq()"""
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
    if filesize_1 is None:
        filesize_1=0
    filesize_2 = check_filesize(http_fullpath_2)
    if filesize_1 is None:
        filesize_2 =0
    if  filesize_1 > filesize_2:
        filesize= filesize_1
        ascp_fullpath=  os.path.join("era-fasp@fasp.sra.ebi.ac.uk:" +path, accession+".fastq.gz")
        ftp_fullpath= os.path.join("ftp://ftp.sra.ebi.ac.uk" + path, accession+".fastq.gz")
        #ftp_fullpath= os.path.join("http://ftp.sra.ebi.ac.uk" + path, accession+".fastq.gz")#remove 
    elif filesize_1 < filesize_2:
        filesize= filesize_2
        ascp_fullpath=  os.path.join("era-fasp@fasp.sra.ebi.ac.uk:" +path, accession+"_1.fastq.gz")
        ftp_fullpath= os.path.join("ftp://ftp.sra.ebi.ac.uk" + path, accession+"_1.fastq.gz")
        #ftp_fullpath= os.path.join("http://ftp.sra.ebi.ac.uk" + path, accession+"_1.fastq.gz")#remove
    else:
        filesize= filesize_1
        ascp_fullpath=  os.path.join("era-fasp@fasp.sra.ebi.ac.uk:" +path, accession+".fastq.gz")
        ftp_fullpath= os.path.join("ftp://ftp.sra.ebi.ac.uk" + path, accession+".fastq.gz")
        #ftp_fullpath= os.path.join("http://ftp.sra.ebi.ac.uk" + path, accession+".fastq.gz")#remove
    return ascp_fullpath, ftp_fullpath, filesize

def check_filesize(http_fullpath):
    """deprecated, replaced by launch_ffq_ftp()"""
    header= subprocess.check_output(f"curl -sI {http_fullpath}", shell=True).decode("utf-8")
    if "404 Not Found" in header:
        return 0
    elif "200 OK" in header:
        filesize= int(header.split("Content-Length: ")[1].split("\r")[0])
        return filesize
        
def launch_ascp(ascp_fullpath, outputdir, filesizelimit=1500000000):
    returncode= subprocess.run([constants.asperapath, "-QT", "-l", "300m", "-P33001", "-@",f"0:{filesizelimit}" ,"-i", constants.aspera_ssh_key , ascp_fullpath, outputdir], 
    stdout=subprocess.DEVNULL, 
    stderr=subprocess.STDOUT)
    return returncode.returncode

def launch_curl(ftp_fullpath, outputdir, filesizelimit=1500000000):
    returncode= subprocess.run(["curl", "-r", f"0-{str(filesizelimit)}","-o", outputdir, ftp_fullpath ],
    stdout=subprocess.DEVNULL,
    stderr=subprocess.STDOUT)
    return returncode.returncode
    
    

def launch_ffq_ftp(accession):
    """CLI command to fetch accession metadata. This function captures the stdout stream, decode into json and load as list"""
    completedprocess = subprocess.run(["ffq" ,"--ftp", accession], capture_output=True)
    returncode= completedprocess.returncode
    stdout= completedprocess.stdout.decode("utf-8")
    return json.loads(stdout)

def get_download_path_ffq(accession):
    """wrapper for launch_ffq_ftp(). Parse json output to return download links of largest read file"""
    ftp_metadata = launch_ffq_ftp(accession)
    #hard-coded to retry fetching metadata ## this is a quick fix to occational failure of ffq
    while type(ftp_metadata) != list:
        ftp_metadata = launch_ffq_ftp(accession)
    if ftp_metadata == []:
        return "NOT_FOUND", "NOT_FOUND",0   
    filesizes = [datadict.get("filesize") for datadict in ftp_metadata if datadict.get("filetype") == "fastq"]
    if ftp_metadata == []:
        return "NOT_FOUND", "NOT_FOUND",0
    ftp_fullpath =  ftp_metadata[filesizes.index(max(filesizes))].get("url")
    ascp_fullpath = ftp_fullpath.replace("ftp://ftp.sra.ebi.ac.uk/", "era-fasp@fasp.sra.ebi.ac.uk:")
    return ascp_fullpath , ftp_fullpath, max(filesizes)
    
__all__=["get_download_path", "check_filesize", "launch_ascp","launch_curl" , "get_download_path_ffq"]