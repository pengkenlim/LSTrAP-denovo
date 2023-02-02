import os
import sys
if __name__ == "__main__":
        abspath= os.getcwd()
        parent_module= os.path.join(abspath.split("LSTrAP-denovo")[0], "LSTrAP-denovo")
        sys.path.insert(0, parent_module)

import numpy as np
import concurrent.futures
import subprocess
from datetime import datetime
import argparse
import pandas as pd
from setup import constants


        
from assembly import misc

def extract_ORFs(filepath,dirname):
    return_code= os.system(f"cd {working_dir}; " + os.path.join(transdecoder_bin_dir, "TransDecoder.LongOrfs") + " " + 
    f"-t {filepath} --output_dir {dirname} -G {genetic_code} -m {str(min_prot_len)} 2>{working_dir}/{dirname}.logs")
    return return_code
    
def predict_ORFs(filepath,dirname, domtbloutpath):
    return_code= os.system(f"cd {working_dir}; " + os.path.join(transdecoder_bin_dir, "TransDecoder.Predict") + " " + 
    f"-t {filepath} --output_dir {dirname} -G {genetic_code} --retain_pfam_hits {domtbloutpath} 2>>{working_dir}/{dirname}.logs")
    return return_code


def Pfam_hmmsearch(outdir, domtbloutpath):
    input_pep = os.path.join(outdir, "longest_orfs.pep") #path/to/splitfile_partx/longest_orfs.pep
    logpath = os.path.join(outdir, "PfamHMM.log") #path/to/splitfile_partx/PfamHMM.log
    return_code= os.system(f"{hmmsearch_bin} --cpu {threads} --domtblout {domtbloutpath} {pathtoPfamHMM} {input_pep} > {logpath}")
    return return_code
    
def swap_target_query(input_domtblout, output_domtblout):
    with open(input_domtblout, "r") as f:
        contents = f.readlines()
    with open(output_domtblout, "w")as f:
        for line in contents:
            if "#" in line:
                f.write(line)
            else:
                f.write(line[38:76]+line[:38]+line[76:])


def run_job(file_name):
    print(f"{file_name}: Running Transdecoder.LongOrfs...")
    filepath = os.path.join(working_dir ,file_name) #path/to/splitfile_partx.fasta
    outdir = os.path.join(working_dir ,file_name.split(".fasta")[0]) #path/to/splitfile_partx
    
    #extract ORF
    return_code = extract_ORFs(filepath, file_name.split(".fasta")[0])
    if return_code ==0:
        print(f"{file_name}: Transdecoder.LongOrfs done. Running hmmsearch for Pfam domains...")
    else:
        return f"{file_name}: ERROR. Transdecoder.LongOrfs Failed."
    
    domtbloutpath = os.path.join(outdir, "longest_orfs.domtblout")
    
    #Pfam hmmsearch
    if not os.path.exists(outdir+".__checkpoints"): #check if not done
        return_code= Pfam_hmmsearch(outdir, domtbloutpath)
        if return_code ==0:
            print(f"{file_name}: hmmsearch done. Running Transdecoder.Predict...")
        else:
            return f"{file_name}: ERROR. hmmsearch Failed."
    else:
        print(f"{file_name}: hmmsearch done. Running Transdecoder.Predict...")
    
    #flipping targets and queries in domtblouts output by hmmsearch so that it is consistent with hmmscan
    flipped_domtbloutpath = os.path.join(outdir, "longest_orfs_flipped.domtblout")
    #os.system("awk \'BEGIN{OFS=FS=\" \"} NR<=3{print}; NR>3{tmp=$1; $1=$4; $4=tmp; tmp=$2; $2=$5; $5=tmp; print}\'" + f" {domtbloutpath} > {flipped_domtbloutpath}")
    swap_target_query(domtbloutpath, flipped_domtbloutpath)
    
    #predcit ORFs hmmsearch
    return_code= predict_ORFs(filepath,file_name.split(".fasta")[0], flipped_domtbloutpath)
    endtime = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
    if return_code ==0:
        return f"{file_name} completed at {endtime}"
    else:
        return f"{file_name}: ERROR. Transdecoder.Predict Failed."
    
def combine():
    #combine cds extracted from each splitfile
    os.system("cat "+ os.path.join(working_dir,"*.transdecoder.cds")+">" + os.path.join(annot_dir,"cds_from_transcripts.fasta"))
    print(os.path.join(annot_dir,"cds_from_transcripts.fasta") + " generated")

    #combine pep
    os.system("cat "+ os.path.join(working_dir,"*.transdecoder.pep")+">" + os.path.join(annot_dir,"translated_cds.fasta"))
    print(os.path.join(annot_dir,"translated_cds.fasta") + " generated")
    
    #combine gff3
    os.system("cat "+ os.path.join(working_dir,"*.transdecoder.gff3")+">" + os.path.join(annot_dir,"transcripts.gff3"))
    print(os.path.join(annot_dir,"transcripts.gff3") + " generated")
    
    #copy trinity output assembly inot annot dir
    os.system("cp "+ fastapath + " " + os.path.join(annot_dir,"transcripts.fasta"))
    print(os.path.join(annot_dir,"transcripts.fasta") + " generated")
    
    #combine domtblout from hmmsearches. Establish headers.
    os.system("head -n 3 " + os.path.join(working_dir,"splitfile_part1", "longest_orfs_flipped.domtblout") + ">" + os.path.join(annot_dir,"translated_cds.domtblout"))
    os.system("awk \'!/#/\' " + os.path.join(working_dir,"splitfile_part*", "longest_orfs_flipped.domtblout") + ">>" + os.path.join(annot_dir,"translated_cds.domtblout"))
    print(os.path.join(annot_dir,"translated_cds.domtblout") + " generated")
    
    #cd-hit to cluster isoforms
    pathtocds = os.path.join(annot_dir,"cds_from_transcripts.fasta")
    outputpath = os.path.join(annot_dir, "Transcript_isoforms")
    os.system(f"{constants.cdhitpath} -p 1 -b 3 -t 10 -T {threadpool} -M 0 -c 0.98 -i {inputpath} -o {outputpath}")
    os.system(f"rm {outputpath}")
    
    #Select primary transcripts to make cds_from_primary_transcripts.fasta
    ##read cds_from_transcripts.fasta and extract coding score
    with open(pathtocds, "r") as f:
        seqinfo = f.read().split(">")[1:]
    seq_info_dict = {}
    for chunk in seqinfo:
        seqid = chunk.split("\n")[0].split(" ")[0]
        score = float(chunk.split("\n")[0].split("score=")[1].split(" ")[0]. split(",")[0])
        seq_info_dict[seqid]={"fasta_content" : chunk ,"Coding_score" : score}
    
    #read cluster information 
    pathtoclstr = outputpath+ ".clstr"
    with open(pathtoclstr, "r") as f:
        clstr_contents = f.read().split(">Cluster ")[1:]
    cluster_dict = {}
    for chunk in clstr_contents:
        clusterid= chunk.split("\n")[0]
        seqid = chunk.split(">")[1:]
        seqid= [ line.split("...")[0] for line in seqid]
        cluster_dict[clusterid] = seqid
    
    #Select isoform with highest coding score as cds_from_primary_transcripts.fasta    
    primary_cds_path = os.path.join(annot_dir,"cds_from_primary_transcripts.fasta")
    with open(path_to_fasta_subset, "w") as f:
        for clusterid , seqid_list in cluster_dict.items():
            scoretobeat = -10000
            seqtokeep= ""
            for seqid in seqid_list:
                score = seq_info_dict.get(seqid).get("Coding_score")
                if score > scoretobeat:
                    scoretobeat=score
                    seqtokeep=seqid
            f.write(">"+ seq_info_dict.get(seqtokeep).get("fasta_content"))
        
    


def parse_domtblout(dombloutpath):
    with open(dombloutpath, "r") as f:
        contents= f.read()
    contents= [line for line in contents.split("\n") if "#" not in line and line != "" ]
    cds_annot_dict = {}
    for line in contents:
        target = line[21:31].strip()
        query = line[38:59].strip()
        i_E_val = float(line[117:127].strip())
        if query not in cds_annot_dict.keys():
            cds_annot_dict[query] = {target: i_E_val}
        else:
            if target not in cds_annot_dict[query].keys():
                cds_annot_dict[query][target] = i_E_val
            else:
                if cds_annot_dict[query][target] > i_E_val:
                    cds_annot_dict[query][target] = i_E_val
    return cds_annot_dict
    
def parse_interpro2go(interpro2gopath):
    with open(interpro2gopath, "r") as f:
        interpro2go= f.read().split("\n")
        interpro2go_dict = {}
    for line in interpro2go:
        if "InterPro:" in line:
            acc = line.split("InterPro:")[1].split(" ")[0]
            GO_term= line.split("; ")[1]
            if acc not in interpro2go_dict.keys():
                interpro2go_dict[acc]= [GO_term]
            else:
                interpro2go_dict[acc] += [GO_term]
    return interpro2go_dict
    

def create_annotation_desc():
    cds_annot_dict = parse_domtblout(os.path.join(annot_dir,"translated_cds.domtblout"))
    Pfam_descriptions = pd.read_csv(os.path.join(data_dir, "Pfam_descriptions.tsv"), sep = "\t")
    Pfam2desc_dict = Pfam_descriptions[["Accession", "Name"]].set_index("Accession").to_dict()["Name"]
    Pfam2interpro_dict = Pfam_descriptions[["Accession", "Integrated Into"]].set_index("Accession").to_dict()["Integrated Into"]
    interpro2go_dict = parse_interpro2go(os.path.join(data_dir, "interpro2go"))
    annotation_df = pd.DataFrame()
    cds_col = cds_annot_dict.keys()
    pfam_col= []
    des_col= []
    interpro_col = []
    GO_col=[]
    for cds in cds_col:
        Eval_pfam_tuple = [(value, key) for key, value in cds_annot_dict[cds].items()]
        Eval_pfam_tuple.sort
        pfam_list = [pfam for Eval, pfam in Eval_pfam_tuple]
        interpro_list =[]
        des_list= []
        GO_list=[]
        for pfam in pfam_list:
            des_list += [Pfam2desc_dict[pfam.split(".")[0]]]
            interpro_list += [Pfam2interpro_dict[pfam.split(".")[0]]]
            GO_list += interpro2go_dict.get(Pfam2interpro_dict[pfam.split(".")[0]], [])
        interpro_list = [interpro for interpro in interpro_list if type(interpro)!=float]
        pfam_col+= [pfam_list]
        interpro_col += [list(set(interpro_list))]
        des_col += [list(set(des_list))]
        GO_col += [list(set(GO_list))]
    #make dataframe
    annotation_df["Protein"] = cds_col
    annotation_df["Descriptions"] = des_col
    annotation_df["Pfam Domains"] = pfam_col
    annotation_df["Interpro Entries"] = interpro_col
    annotation_df["Go Terms"] = GO_col
    annotation_df = annotation_df.set_index("Protein")
    return annotation_df
        
        
    
    

    
if __name__ == "__main__":
    #specify arguments
    parser= argparse.ArgumentParser(description="LSTrAP-denovo.AnnotatePredictORFs.py: Helper script that uses Transdecoder's ORF prediction to extract Coding sequences (CDS) and annotate protein function based on Pfam domains.\n \
    NOTE: This is to be run following completion of the LSTrAP-denovo pipeline and Trinity transcriptome assembly (Run_Trinity).\n\
    Please make sure binaries for TransDecoder, hmmsearch, hmmpress are installed and added into $PATH.\n\
    \n\
    Refer to https://github.com/pengkenlim/LSTrAP-denovo for more information on pipeline usage and implementation .\n\
    Refer to https://github.com/TransDecoder/TransDecoder/wiki for more information on Transdecoder.\n\
    Refer to http://hmmer.org/documentation.html for more information on hmmsearch and hmmpress.\n\
    Refer to https://www.ebi.ac.uk/interpro/ for more information pfam domains.")
    
    parser.add_argument("-o", "--output_dir", type=str, metavar= "", required=True,
    help= "Directory for data output. Directory needs to be same as for step 1 and step 2 of LSTrAP-denovo pipeline (i.e. MakeDraftCDS.py, SelectAccessions.py).")
    
   
    parser.add_argument("-t", "--threads", type=int, metavar="", default=4, 
    help = "Total thread pool for workers. Needs to be divisible by number of workers.")
    
    parser.add_argument("-w", "--workers", type=int, metavar="", default=2, 
    help= "Specifies the maximum workers for running Transdecoder and hmmsearch in parallel. Reccomended to set to (thread pool)/3. Set to 2 by default.")
    
    parser.add_argument("-G", "--genetic_code", type=str, metavar="", default="Universal", 
    help= "Genetic code passed to TransDecoder.LongOrfs and Transdecoder.Predict. Set to \"Universal\" by default. Refer to ... for more information.")
    
    parser.add_argument("-m", "--min_prot_len", type=int, metavar="", default=100,
    help= "Minimal protein length (aa) passed to TransDecoder.LongOrfs. Set to 100 by default.")
    
    parser.add_argument("-pdir","--pfam_dir", type=str, metavar="", required=True,
    help="Path to directory containing Pfam hmm model (Pfam-A.hmm). If directory specified do not contain Pfam-A.hmm, it will be downloaded auotmatically." )
    
    parser.add_argument("-hbdir", "--hmm_bin_dir", type=str, metavar= "", default="",
    help= "Path to directory containing hmmsearch and hmmpress binaries. Not required if directory has been added to $PATH.")
    
    parser.add_argument("-tbdir", "--transdecoder_bin_dir", type=str, metavar= "", default="",
    help= "Path to directory containing transdecoder binaries (Transdecoder.LongORFs and Transdecoder.Predict). Not required if hmmsearch binary directory has been added to $PATH.")
    
    #banner
    misc.print_logo("AnnotatePredictORFs.py")
    
    
    #parse args
    args=parser.parse_args()

    output_dir= args.output_dir
    threadpool =args.threads
    workers =args.workers
    genetic_code = args.genetic_code
    min_prot_len= args.min_prot_len
    pfam_dir= args.pfam_dir
    hmm_bin_dir = args.hmm_bin_dir
    transdecoder_bin_dir = args.transdecoder_bin_dir
    hmmsearch_bin= os.path.join(hmm_bin_dir, "hmmsearch")
    hmmpress_bin= os.path.join(hmm_bin_dir, "hmmpress")

    print("\nAnnotatePredictORFs.py started running on ", datetime.now().strftime("%d/%m/%Y %H:%M:%S")+ "\n")
    
    #check for trinity_dir, findout name of assembly fasta. set path for assembly fasta.
    Trinity_dir= os.path.join(output_dir, "Trinity_output")
    fastapath=os.path.join(Trinity_dir, "transcripts.fasta")
    if not os.path.exists(Trinity_dir):
         sys.exit("Error in output directory. Trinity directory not found. Exiting...")
    if not os.path.exists(fastapath):
        sys.exit(f"Error in output directory. Trinity transcriptome assembly fasta not found in {fastapath}. Exiting...")

    print(f"\nAssembly file detected in {fastapath}.\n")
    
    #check and set threads. Check if threadpool can be divided by number of workers. Else set to appropriate threads
    if threadpool%workers != 0:
        print(f"Threadpool of {threadpool} specified by user is not divisible by {workers} workers. Threadpool corrected to {int((threadpool - (threadpool%workers)))}.\n")
        threadpool = int((threadpool - (threadpool%workers)))
    
    threads=int(threadpool/workers)
    
    #check if path to hmmsearch, hmmpress and transdecoder binaries are valid.
    if hmm_bin_dir!= "":
        if not os.path.exists(hmmsearch_bin):
            sys.exit(f"Error: hmmsearch not found at {hmm_bin_dir}. Exiting...")
        if not os.path.exists(hmmpress_bin):
            sys.exit(f"Error: hmmpress not found at {hmm_bin_dir}. Exiting...")
    
    
    if transdecoder_bin_dir!= "":
        if not os.path.exists(os.path.join(transdecoder_bin_dir, "TransDecoder.LongOrfs")):
            sys.exit(f"Error: TransDecoder.LongOrfs not found in {transdecoder_bin_dir}. Exiting...")
        if not os.path.exists(os.path.join(transdecoder_bin_dir, "TransDecoder.Predict")):
            sys.exit(f"Error: TransDecoder.Predict not found in {transdecoder_bin_dir}. Exiting...")

    #check if LSTrAP-denovo/data folder exists
    data_dir = os.path.join(abspath, "data")
    if not os.path.exists(data_dir):
        sys.exit("Error: Cannot find data directory, please make sure script is run in LSTrAP-denovo directory (i.e. /path/to/LSTrAP-denovo/). Exiting...")
    
    #check if Pfam hmm is downloaded in directory

    if not os.path.exists(pfam_dir):
        print(f"User specified pfam directory ({pfam_dir}) is created as it cannot be found in the system.")
        os.makedirs(pfam_dir)
    pathtoPfamHMM =os.path.join(pfam_dir, "Pfam-A.hmm")
    if not os.path.exists(pathtoPfamHMM):
        print(f"Pfam-A.hmm not detected in user specifed pfam directory. Attempting to download from ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam35.0/Pfam-A.hmm.gz...")
        return_code = os.system(f"wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam35.0/Pfam-A.hmm.gz --directory-prefix={pfam_dir}")
        if return_code != 0:
            sys.exit("Error in downloading Pfam-A.hmm. Please download and decompress manaully into Pfam directory . Exiting...")
        print("un-gunzipping Pfam-A.hmm ... ")
        return_code = os.system("gunzip " + os.path.join(pfam_dir,"Pfam-A.hmm.gz"))
        if return_code != 0:
            sys.exit("Error in un-gunzipping Pfam-A.hmm. Please download and decompress manaully into Pfam directory . Exiting...")
    if not os.path.exists(pathtoPfamHMM):
        sys.exit("Error in un-gunzipping Pfam-A.hmm. Please download and decompress manaully into Pfam directory . Exiting...")
    
    if not os.path.exists(os.path.join(pfam_dir,"Pfam-A.hmm.h3i")):
        print(f"Generating hmm database from Pfam-A.hmm using hmmpress...")
        return_code = os.system(f"{hmmpress_bin} {pathtoPfamHMM}")
        if return_code != 0:
            sys.exit("Error in running hmmpress on Pfam-A.hmm. Exiting...")
    
    #make working directory if not already exists
    working_dir= os.path.join(output_dir, "AnnotatePredictORFs")
    if not os.path.exists(working_dir):
        os.makedirs(working_dir)
    
    #read content of fasta file
    with open(fastapath, "r") as f:
       contents = f.read()
       contents= contents.split(">")
    
    print(f"\nTotal number of transcripts in assembly fasta: {len(contents)}")
    n_seq_per_file = int((len(contents) - (len(contents)%workers))/workers)
    
    print(f"Number of specified workers= {workers} \n{len(contents)} Sequences in assembly fasta file will be split into {workers} splitfiles each containing approx. {n_seq_per_file} sequences.\n")
    
    #create each seq_chunk in working_dir
    file_names=[]
    for i in range(0,workers):
        if i == workers-1:
            towrite= contents[(i*n_seq_per_file):]
        else:    
            towrite= contents[(i*n_seq_per_file):((i+1)*n_seq_per_file)]
        towrite=">".join(towrite)
        with open(os.path.join(working_dir, f"splitfile_part{i+1}.fasta"), "w")as f:
            f.write(towrite)
        file_names += [f"splitfile_part{i+1}.fasta"]
        print(f"splitfile_part{i+1}.fasta created")
    
    print("\nFinished splitting files. Annotating splitfiles in parallel...\n")
    with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as executor:
        results= [executor.submit(run_job, file_name) for file_name in file_names]
        for f in concurrent.futures.as_completed(results):
            print(f.result())
    
    print("\nParallel anotation of splitfiles completed.\n")
    #make dir to hold annotations
    annot_dir = os.path.join(output_dir, "Annotations") 
    if not os.path.exists(annot_dir):
        os.makedirs(annot_dir)
    
    #combine transdecoder output to form annotation files
    combine()

    #write cds, descriptions, pfam , interopro entries and Go terms
    annotation_df = create_annotation_desc()
    annotation_df.to_csv(os.path.join(annot_dir,"cds_annotations.tsv"), sep="\t")
    print("\nDescriptions, pfam domains, interopro entries and Go terms associated with each Coding sequence written to" + os.path.join(annot_dir,"cds_annotations.tsv"))

    print("\nAnnotatePredictORFS.py finished running on ", datetime.now().strftime("%d/%m/%Y %H:%M:%S"))

