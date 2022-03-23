#setting sys.path for importing modules
import os
import sys
if __name__ == "__main__":
        abspath= os.getcwd()
        parent_module= os.path.join(abspath.split("LSTrAP-denovo")[0], "LSTrAP-denovo")
        sys.path.insert(0, parent_module)

import argparse
from download import ena

if __name__ == "__main__":
	#parse arguments
	parser= argparse.ArgumentParser(description="LSTrAP-denovo.preliminary_assembly: Assemble a reduced but high-confidence assembly from public RNA-seq data")
	parser.add_argument("-o", "--output_dir", type=str , metavar= "", required=True, help= "Directory for data output.")
	parser.add_argument("-k", "--kmer_len", type=int, metavar="", default=35, help = "Specify K-mer length for assembly using Soapdenovo-Trans. K-mer length will be set to 35 by default.")
	parser.add_argument("-ct", "--consensus_threshold", type=int ,metavar="", default=0 , help = "Specify consensus threshold during filtering. Threshold will be determined automatically by default.")
	parser.add_argument("-s","--file_size" , type=int, default=1500, help="Specify the size limit(mb) of accession read files to partially download. Limit set to 1500 by default.")
	parser.add_argument("-t", "--threads", type=int, required=True, help = "Total thread pool for workers. Needs to be divisible by number of workers.")
	parser.add_argument("-w", "--workers", type=int, default=2, help= "Specify the maximum workers for running multiple download-assembly jobs in parellel. Set to 2 by default.")

	ME_group_1 = parser.add_mutually_exclusive_group(required=True)
	ME_group_1.add_argument("-i", "--id", type=int, metavar="", help= "NCBI TaxID of organism for fetching SRA run accessions.")
	ME_group_1.add_argument("-a", "--accessions", type=str, metavar="", help= "User-defined list of SRA run accessions to fetch. Requires at least 10 accessions. E.g.: SRR123456,SRR654321,ERR246810,...")
	args=parser.parse_args()
	#assigning arguments to variables
	taxid= args.id
	accessions= args.accessions
	outputdir= args.output_dir
	consensus_threshold= args.consensus_threshold
	file_size= args.file_size
	threads= args.threads
	workers=args.workers
	#check if threads pool is divisible by number of workers
	if threads % workers != 0:
		print(f"Specified thread pool of {args.threads} is not divisible by number of workers.")
		threads= threads - (threads % workers)
		print(f"Using thread pool of {threads} instead.")
	#create outputdir if not found
	#if not os.path.exists(outputdir):
		#os.makedirs(outputdir)

	#check if accessions are given
	if accessions is not None:
		accessions = accessions.split(",")
		if len(accessions) < 10:
			sys.exit("Not enough accessions provided. Refer to --help for more information.")
	#check if taxid is given
	elif taxid is not None:
		scientific_name= ena.get_sciname(taxid)
		if type(scientific_name) is  not list:
			sys.exit("TaxID {taxid} is invalid/not found.")
		elif len(scientific_name) > 1:
			sys.exit("More than one organism found for TaxID {taxid}.")
		else:
			scientific_name= scientific_name[0]["scientific_name"]
			print(f"Fetching RNA-seq accessions of {scientific_name}, NCBI TaxID {taxid} form ENA..")
			accessions = ena.get_runs(taxid)
			print(f"Total accessions: {len(accessions)}")
