#setting sys.path for importing modules
import os
import sys
if __name__ == "__main__":
        abspath= os.getcwd()
        parent_module= os.path.join(abspath.split("LSTrAP-denovo")[0], "LSTrAP-denovo")
        sys.path.insert(0, parent_module)

import argparse
from time import sleep

from setup import constants
from assembly import misc

if __name__ == "__main__":
    #arguments
    parser= argparse.ArgumentParser(description="Run_Trinity.py: Helper script for de novo transcriptome assembly using Trinity.\n \
    NOTE: This is to be run following completion of the HSS-Trans pipeline. Please make sure that binaries for Trinity, its associated dependencies (e.g. bowtie) are installed and added into $PATH.\
    Refer to https://github.com/pengkenlim/LSTrAP-denovofor more information on pipeline usage and implmentation\n\
    Refer to https://github.com/trinityrnaseq/trinityrnaseq/wiki for more information on Trinity installation, usage and implmentation")
    
    parser.add_argument("-o", "--output_dir", type=str, metavar= "", required=True,
    help= "Directory for data output. Directory needs to be same as for step 1 and step 2 of LSTrAP-denovo pipeline (i.e. MakeDraftCDS.py, SelectAccessions.py).")
    
    parser.add_argument("-trin_args", "--trinity_arguments", type=str, metavar="", default="--full_cleanup --max_memory 10G", 
    help = "Arguments in quotes to be passed verbatim into Trinity (e.g. \"--CPU 32 --max_memory 32G\" ). \
    Do not include arguments for output directory and sample paths.\
    If not specified, Trinity will run on memory and CPU as per its default settings with full_cleanup enabled.\n")
    
    parser.add_argument("-no", "--no_overassembly", action="store_true",
    help = "Do not assemble multiple cluster-specific transcript assemblies and concatenate them to generate an over-assembly. Instead, just feed all samples into Trinity (vanilla).")
    
    
    #banner
    #banner
    misc.print_logo("RunTrinity.py")
    
    args=parser.parse_args()
    outputdir= args.output_dir
    trinity_arguments = args.trinity_arguments
    no_overassembly = args.no_overassembly
    pathtotrinity= outputdir
    s_fastqdir = os.path.join(outputdir, "Step_2", "selected_accessions") #dir holding selected accessions
    TSV_file = os.path.join(outputdir,"Samples_for_trinity.tsv")
    #checking outputdir...
    if not os.path.exists(TSV_file):
        sys.exit("Error in output directory. Samples_for_trinity.tsv file not found. Exiting...")
    if not os.path.exists(s_fastqdir):
         sys.exit("Error in output directory. Downloaded accessions from HSS-Trans pipeline not found. Exiting...")
    
    #make dir to hold Trinity ouput
    Trinity_dir= os.path.join(outputdir, "Trinity_output")
    if not os.path.exists(Trinity_dir):
        os.makedirs(Trinity_dir)
        
    print("\nChecking if paths in Samples_for_trinity.tsv are valid...")
    with open(TSV_file, "r") as f:
        TSV_contents = f.read().split("\n")[:-1]
    
    Download = True
    for line in TSV_contents:
        acc_file_path = line.split("\t")[2]
        if not os.path.exists(acc_file_path):
            print(f"WARNING:{acc_file_path} does not exists.")
            Download = False
        acc_file_path = line.split("\t")[3]    
        if not os.path.exists(acc_file_path):
            print(f"WARNING:{acc_file_path} does not exists.")
            Download = False
    if Download == False:
        sys.exit("Some files are missing. Exiting...")
    else:
        print("All paths are valid.\n\n")
    
    if trinity_arguments is None:
        trinity_arguments=""
        
    #vanilla-Trinity
    if no_overassembly:
        trinity_outpath = os.path.join(Trinity_dir, "Trinity_all_samples")
        leftstr= ",".join([line.split("\t")[2] for line in TSV_contents])
        rightstr=",".join([line.split("\t")[3] for line in TSV_contents])
        cmdstr= f"Trinity {trinity_arguments} --left {leftstr} --right {rightstr} --output {trinity_outpath} --seqType fq"
        pathtofinal=os.path.join(Trinity_dir, "transcripts.fasta")
        
        #check if assembly has been prviously completed
        if not os.path.exists(trinity_outpath+ ".fasta"):
            print(f"Running Trinity on all downloaded files using command: \n{cmdstr}\n\n")
            os.system(cmdstr)
            #rename trinity output fasta
            os.system(f"cp {trinity_outpath}.fasta {pathtofinal}")
        else:
            sys.exit("Assembly has been previously completed. Skipping...")
        print(f"\nRun_Trinity.py completed. Trinity assembly is available at:\n{pathtofinal}")
            
    #Cluster-specific assembly
    else:
        print("\nAssembling multiple cluster-specific transcript assemblies for over-assembly...\n")
        #get total number of clusters
        clusters = list({line.split("\t")[0] for line in TSV_contents})
        print(f"Number of clusters-specific assemblies to assemble:{len(clusters)}")
        sleep(2)
        for cluster in clusters:
            
            leftstr= ",".join([line.split("\t")[2] for line in TSV_contents if line.split("\t")[0] == cluster])
            rightstr=",".join([line.split("\t")[3] for line in TSV_contents if line.split("\t")[0] == cluster])
            trinity_outpath = os.path.join(Trinity_dir, "Trinity_" + cluster)
            cmdstr = f"Trinity {trinity_arguments} --left {leftstr} --right {rightstr} --output {trinity_outpath} --seqType fq"
            #check if already cluster is already assembled
            if not os.path.exists(trinity_outpath+ ".Trinity.fasta"):
                print(f"Running Trinity on downloaded files of cluster {cluster} using command: \n{cmdstr}\n\n")
                os.system(cmdstr)
            else:
                print(f"Cluster {cluster} has been previously assembled. Moving on to next...")
        
        print("\nAssembly done. Renaming transcripts and concatenating assemblies...")
        fasta_files = [file for file in os.listdir(Trinity_dir) if "Trinity.fasta" in file and "gene_trans_map" not in file]
        pathtoconcat = os.path.join(Trinity_dir, "concat.fasta")
        concat_content=[]
        for file in files:
            count=1
            assembly_id= file.split("_")[0]
            with open(os.path.join(workdir,file), "r") as f:
                contents = f.read().split("\n")
            for line in contents:
                if ">" in line:
                    concat_content += [f">{assembly_id}_{count}"]
                    count += 1
                elif line == "":
                    pass
                else:
                    concat_content += [line]
            print(file,"parsed. Number of transcripts:", count)
        with open(pathtoconcat, "w") as f:
            f.write("\n".join(concat_content))
        
        print("\nLaunching CD-HIT-EST to remove redundant sequences...\n")
        pathtofinal=os.path.join(Trinity_dir, "transcripts.fasta")
        os.system(f"{constants.cdhitpath} -p 1 -b 3 -t 10 -T 0 -M 0 -c 0.98 -i {pathtoconcat} -o {pathtofinal}")
        print(f"\nRun_Trinity.py completed. Trinity over-assembly is available at:\n{pathtofinal}")
        
    
            
            
        
    
    
    