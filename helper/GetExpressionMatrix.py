#this is a throw away version of the code
import os
import sys
import numpy as np
import concurrent.futures
import subprocess
from datetime import datetime
from time import sleep
import argparse


#setting sys.path for importing modules
if __name__ == "__main__":
        abspath= os.getcwd()
        parent_module= os.path.join(abspath.split("HSS-Trans")[0], "HSS-Trans")
        sys.path.insert(0, parent_module)

#quick and easy filepaths
cds_fastapath="/shared/scratch/ken/pipeline_output/TAXID_4565/Step_2/selected_accessions/concat_renamed_wo_all_2_nr.fasta.transdecoder.cds" #In reality will be {outputdir}/Annotations/cds_from_transcripts.fasta
tempdir="/shared/scratch/ken/pipeline_output/TAXID_4565/GetExpressonMatrix"
fastqdir = "/shared/scratch/ken/pipeline_output/TAXID_4565/Step_2/fastq"
pathtologfile= "/shared/scratch/ken/pipeline_output/TAXID_4565/logs.json"
kaldir=os.path.join("/shared/scratch/ken/pipeline_output/TAXID_4565/GetExpressonMatrix/kallisto")
tpm_matpath= os.path.join(kaldir,"Draft_CDS_exp_mat.tsv")
workers= 30
threadpool= 30
threads = 1
retrylimit = 1

from setup import constants
from preprocess import read_map, classify
from assembly import misc

def PS_job(accession,index):
    try:
        #to unsync workers
        if index < workers:
            sleep((index%workers)*2)
            
        #check if processed
        if type(processed_acc_dict.get(accession)) == float:
            return accession, index , "Already processed"
            
        fastqpath= os.path.join(fastqdir, accession + ".fastq.gz")
        #check if fastq file exists
        if not os.path.exists(fastqpath):
            return accession , index , "failed because fastq file not found"
        #Run kallisto pseudoalignment
        kaloutdir= os.path.join(kaldir, accession)
        result= misc.run_with_retries(retrylimit,
                read_map.launch_kallisto_quant,
                [threads, indexpath , kaloutdir , fastqpath],
                f"{accession}: Kallisto pseudoalignment failed. Retrying...",
                f"{accession}: Kallisto pseudoalignment of accession reads against cds_from_transcripts.fasta\n")
        if result == "failed" or not os.path.exists(kaloutdir):
            return accession , index , "PS failed"
        
        #write quantification information to tpm matrix
        map_rate = read_map.write_quant_info(accession, kaloutdir, tpm_matpath)
        print(f"{accession}: Pseudoalignment completed.")
        return accession , index , float(map_rate)
    except:
        accession, index , "ERROR"
    
def parallel_job():
    '''Wrapper to parallelize download and pseudolaignment jobs for each accession'''
    with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as executor:
        progress_bar= tqdm(total=len(passed_accessions), desc="Accessions processed", unit="Acsn", leave=True)
        results= [executor.submit(PS_job, accession, index) for index, accession in enumerate(passed_accessions)]
        for f in concurrent.futures.as_completed(results):
            accession , index , map_rate = f.result()
            if map_rate== "Already processed":
                msg = f"{accession} already processed."
            else:
                with open(pathtoprocessed, "a") as f:
                    f.write(f"{accession}\t{map_rate}\n")
                if map_rate == "failed because fastq file not found":
                    msg= f"{accession}: Aborted. File not found."
                elif map_rate == "PS failed":
                    msg= f"{accession}: Aborted after {retrylimit} retries."
                elif type(map_rate) == float:
                    msg= f"{accession}: processed."
                else:
                    msg= f"{accession}: ERROR. Unknown exception occurred."
            progress_bar.update(1)
            progress_bar.set_postfix_str(s=msg)
            
            print("\n")
        progress_bar.close()



logfile=misc.logfile(os.path.join(outputdir,"logs.json"))
passed_accessions = logfile.contents["Step_2"].get("qc").get("passed")

if not os.path.exists(kaldir):
    os.makedirs(kaldir)
indexpath= os.path.join(kaldir, fastapath.split("/")[-1] + ".index")

if not os.path.exists(indexpath):
    result= misc.run_with_retries(retrylimit,
    read_map.launch_kallisto_index,
    [cds_fastapath, indexpath],
    f"Kallisto index failed. Retrying...",
    f"\nContructing Kallisto index for cds_from_transcripts.fasta...\n")
    if result == "failed":
        sys.exit(f"Error: Kallisto indexing failed after {retrylimit} retries. Exiting...")

pathtoprocessed= os.path.join(tempdir, "processed.tsv")    
    if not os.path.exists(pathtoprocessed):
         with open(pathtoprocessed, "w") as f:
                f.write("Accession\tMap_rate\n")
    #load processed files into dict
    with open(pathtoprocessed, "r") as f:
        processed_acc_dict = {chunk.split("\t")[0]: chunk.split("\t")[1] for chunk in f.read().split("\n") if chunk != "Accession\tMap_rate" and chunk != ""}
        processed_acc_dict = {key : value if "failed" in value or "exception" in value else float(value) for key, value in processed_acc_dict}

parallel_job()



    





kaloutdir= os.path.join(kaldir, accession)
            result= misc.run_with_retries(retrylimit,
            read_map.launch_kallisto_quant,
            [threads, indexpath , kaloutdir , fastqpath],
            f"{accession}: Kallisto pseudoalignment failed. Retrying...",
            f"{accession}: Kallisto pseudoalignment of accession reads against draft CDSs...\n")
            if result == "failed" or not os.path.exists(kaloutdir):
                #os.system(f"rm {fastqpath}") add in final build
                return accession , index , "PS failed"
            map_rate = read_map.write_quant_info(accession, kaloutdir, tpm_matpath)
            #os.system(f"rm {fastqpath}") add in in final build
            print(f"{accession}: Pseudoalignment completed.")
            return accession , index , float(map_rate)
        return accession , index , "Unknown exception"
    except:
        return accession , index , "Unknown exception"


####################################################


def download_PS_job(accession, index):
    '''Job to validate download path -> download via FTP/ascp -> Peudoalignment(PS) by Kallisto'''
    try:
    #to slow down jobs
        sleep(0.2)
        logfile.load()
        #to unsync workers
        if index < workers:
            sleep((index%workers)*5)
        if type(logfile.contents["Step_2"]["processed_acc"].get(accession)) == float:
            return accession, index , "Already processed"
        if logfile.contents["Step_2"]["processed_acc"].get(accession) == "PS failed":
            return accession, index , "PS failed"
        if accession not in logfile.contents["Step_2"]["processed_acc"].keys() or logfile.contents["Step_2"]["processed_acc"].get(accession) == "Download failed" or logfile.contents["Step_2"]["processed_acc"].get(accession) == "Download failed because link not found" or logfile.contents["Step_2"]["processed_acc"].get(accession) == "Unknown exception" : #check if accession has been downloaded/processed. Proceed with download if not.
            ascp_fullpath, ftp_fullpath, filesize = aspera.get_download_path_ffq(accession)
            fastqpath =os.path.join(C_fastqdir,accession+".fastq.gz")
            #remove fastq if exists
            if os.path.exists(fastqpath):
                os.system(f"rm {fastqpath}")
            if ascp_fullpath == "NOT_FOUND":
                return accession , index , "Download failed because link not found"
            #download
            if download_method == "ascp":
                result= misc.run_with_retries(retrylimit,
                aspera.launch_ascp,
                [ascp_fullpath,fastqpath,min([filesizelimit,filesize])],
                f"{accession}: Download failed. Retrying...",
                f"{accession}: Downloading file via ascp...\n")
            elif download_method == "ftp":
                result= misc.run_with_retries(retrylimit,
                aspera.launch_curl,
                [ftp_fullpath,fastqpath,min([filesizelimit,filesize])],
                f"{accession}: Download failed. Retrying...",
                f"{accession}: Downloading file via ftp...\n")
            if result == "failed":
                return accession , index , "Download failed"
            logfile.contents["Step_2"]["processed_acc"][accession]= "Downloaded"

        if  logfile.contents["Step_2"]["processed_acc"].get(accession)== "Downloaded" and os.path.exists(fastqpath):
            if os.path.getsize("fastqpath")==0:
                os.system(f"rm {fastqpath}")
                return accession , index , "Download failed"
            fastqpath =os.path.join(C_fastqdir,accession+".fastq.gz")
            kaloutdir= os.path.join(kaldir, accession)
            result= misc.run_with_retries(retrylimit,
            read_map.launch_kallisto_quant,
            [threads, indexpath , kaloutdir , fastqpath],
            f"{accession}: Kallisto pseudoalignment failed. Retrying...",
            f"{accession}: Kallisto pseudoalignment of accession reads against draft CDSs...\n")
            if result == "failed" or not os.path.exists(kaloutdir):
                #os.system(f"rm {fastqpath}") add in final build
                return accession , index , "PS failed"
            map_rate = read_map.write_quant_info(accession, kaloutdir, tpm_matpath)
            #os.system(f"rm {fastqpath}") add in in final build
            print(f"{accession}: Pseudoalignment completed.")
            return accession , index , float(map_rate)
        return accession , index , "Unknown exception"
    except:
        return accession , index , "Unknown exception"
        
def parallel_job(workers):
    '''Wrapper to parallelize download and pseudolaignment jobs for each accession'''
    #logfile.load()
    with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as executor:
        progress_bar= tqdm(total=len(accessions), desc="Accessions processed", unit="Acsn", leave=True)
        #results= [executor.submit(download_PS_job, accession, index) for index, accession in enumerate(accessions)]
        results= [executor.submit(runjob,download_PS_job, accession, index, 1200) for index, accession in enumerate(accessions)]
        for f in concurrent.futures.as_completed(results): # 20 minutes
            accession , index , map_rate = f.result()
            if map_rate== "Already processed":
                msg = f"{accession} already processed."
            else:
                with open(pathtoprocessed, "a") as f:
                    f.write(f"{accession}\t{map_rate}\n")
                if map_rate == "Download failed because link not found":
                    msg= f"{accession}: Aborted. Download link not found."
                elif map_rate == "Download failed":
                    msg= f"{accession}: Aborted. Download failed."
                elif map_rate == "PS failed":
                    msg= f"{accession}: Aborted after {retrylimit} retries."
                elif type(map_rate) == float:
                    msg= f"{accession}: processed."
                else:
                    msg= f"{accession}: ERROR. Unknown exception occurred."
            progress_bar.update(1)
            progress_bar.set_postfix_str(s=msg)
            
            print("\n")
        progress_bar.close()
        #logfile.load()