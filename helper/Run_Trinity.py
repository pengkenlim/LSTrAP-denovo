#setting sys.path for importing modules
import os
import sys
if __name__ == "__main__":
        abspath= os.getcwd()
        parent_module= os.path.join(abspath.split("HSS-Trans")[0], "HSS-Trans")
        sys.path.insert(0, parent_module)

import argparse
from time import sleep

from setup import constants

if __name__ == "__main__":
    #arguments
    parser= argparse.ArgumentParser(description="Run_Trinity.py: Helper script for de novo transcriptome assembly using Trinity.\n \
    NOTE: This is to be run following completion of the HSS-Trans pipeline. Please make sure that Trinity and its associated dependencies are installed and added into $PATH.\
    Refer to https://github.com/pengkenlim/HSS-Trans for more information on pipeline usage and implmentation\n\
    Refer to https://github.com/trinityrnaseq/trinityrnaseq/wiki for more information on Trinity installation, usage and implmentation")
    parser.add_argument("-o", "--output_dir", type=str, metavar= "", required=True,
    help= "Directory for data output. Directory needs to be same as for step 1 and step 2 of HSS-Trans pipeline (i.e. MakeDraftCDS.py, SelectAccessions.py).")
    parser.add_argument("-ta", "--trinity_arguments", type=str, metavar="", default="--full_cleanup --max_memory 10G", 
    help = "Arguments in quotes to be passed verbatim into Trinity (e.g. \"--CPU 32 --max_memory 32G\" ). \
    Do not include arguments for output directory and sample paths.\
    If not specified, Trinity will run on memory and CPU as per its default settings with full_cleanup enabled.\n")
    parser.add_argument("-no", "--no_overassembly", action="store_true",
    help = "Do not assemble multiple cluster-specific transcript assemblies and then concatenate them to generate an over-assembly. Instead just use feed all samples into Trinity (vanilla).")
    
    #banner
    print("\n \
    _  _ ___ ___    _____\n\
    | || / __/ __|__|_   _| _ __ _ _ _  ___\n\
    | __ \__ \__ \___|| || '_/ _` | ' \(_-<\n\
    |_||_|___/___/    |_||_| \__,_|_||_/__/ Run_Trinity.py\n")
    
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
    Trinity_dir= os.path.join(outputdir, "trinity_output")
    if not os.path.exists(Trinity_dir):
        os.makedirs(Trinity_dir)
        
    print("Checking if paths in Samples_for_trinity.tsv are valid...")
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
        trinity_path = os.path.join(Trinity_dir, "Trinity_all_samples")
        leftstr= ",".join([line.split("\t")[2] for line in TSV_contents])
        rightstr=",".join([line.split("\t")[3] for line in TSV_contents])
        cmdstr= f"Trinity {trinity_arguments} --left {leftstr} --right {rightstr} --output {trinity_path} --seqType fq"
            
        
        print(f"Running Trinity on all downloaded files using command: \n{cmdstr}\n\n")
        os.system(cmdstr)
    else:
        print("Assembling multiple cluster-specific transcript assemblies for over-assembly...")
        #get total number of clusters
        clusters = list({line.split("\t")[0] for line in TSV_contents})
        print(f"Number of clusters-specific assemblies to assemble:{len(clusters)}")
        sleep(2)
        for cluster in clusters:
            #generate a new TSV file
            #cs_TSV_file = os.path.join(Trinity_dir,cluster+"_samples.tsv")
            #cs_TSV_content = "\n".join([line for line in TSV_contents if line.split("\t")[0] == cluster])+"\n"
            #with open(cs_TSV_file, "w") as f:
            #    f.write(cs_TSV_content)
            #trinity_path = os.path.join(Trinity_dir, "Trinity_" + cluster)
            #cmdstr = f"Trinity {trinity_arguments} --samples_file {cs_TSV_file} --output {trinity_path}"
            #print(f"Running Trinity on all downloaded files using command: \n{cmdstr}\n\n")
            #os.system(cmdstr)
            
            leftstr= ",".join([line.split("\t")[2] for line in TSV_contents if line.split("\t")[0] == cluster])
            rightstr=",".join([line.split("\t")[3] for line in TSV_contents if line.split("\t")[0] == cluster])
            trinity_path = os.path.join(Trinity_dir, "Trinity_" + cluster)
            cmdstr = f"Trinity {trinity_arguments} --left {leftstr} --right {rightstr} --output {trinity_path} --seqType fq"
            print(f"Running Trinity on all downloaded files using command: \n{cmdstr}\n\n")
            os.system(cmdstr)
        
        print("Assembly done. Renaming transcripts and concatenating assemblies...")
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
        
        print("Launching CD-HIT-EST to remove redundant sequences...")
        pathtofinal=os.path.join(Trinity_dir, "Trinity_over-assembly_nr.fasta")
        os.system(f"{constants.cdhitpath} -p 1 -b 3 -t 10 -T 0 -M 0 -c 0.98 -i {pathtoconcat} -o {pathtofinal}")
        
        print(f"Run_Trinity.py completed. Trinity over-assembly is available at:\n{pathtofinal}")
        
            
            
        
    
    
    