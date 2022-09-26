# HSS-Trans
 **High-throughput Sample Selection for Transcriptome assembly**  
 An automated pipeline to quality control and select public RNA-seq accessions for the transcriptome assembly of species without reference genomes using unsupervised machine learning
 
 **Features**  
 - Simplest implementation of HSS-Trans only requires taxonomic id for species of interest
 - Construction of high confidence draft coding sequences (CDSs) using an efficient *de novo* assembler (SOAPdenovo-Trans) coupled with a novel consensus-based approach to retain common CDSs assembled independently in multiple RNA-seq accessions
 - Rapid identification of low quality RNA-seq accessions (quality control) without the need of reference genome annotations
 - High-speed download of RNA-seq data using IBM Aspera file transfer framework. Option to download via FTP using cURLis also available.
 - Parallel pseudoalignment of accessions against draft CDS to obtain gene expression data using Kallisto
 - selection of representative accessions gene expression dta using unsupervised machine learning (k-means clustering)
 - Automatic download of selected accessions at the end of pipeline
 - Option to seamlessly run popular *de novo* transcriptome assembler, Trinity at the end of the pipeline
 - Informative HTML report for run details  

For more information refer to: *Place-holder for publication DOI*

# Setup
**clone repository to local machine**
```
git clone https://github.com/pengkenlim/HSS-Trans.git
```
**Create environment and install relavent packages (grab the latest version of ffq)**
```
cd HSS-Trans
virtualenv -p python3 <MY_ENV>
source ./<MY_ENV>/bin/activate
pip install --upgrade pip
pip install -r ./setup/requirements.txt
pip install ffq
```
**Install dependencies into programs sub-directory **
```
python3 ./setup/install.py
```

# Step 1. Assembling Draft CDSs  
**Simplest implementation**
```
python3 ./main/MakeDraftCDS.py --output_dir <working directory> -i <NCBI TaxID> -g <Genetic code>
```

Full options:
```
     _  _ ___ ___    _____
    | || / __/ __|__|_   _| _ __ _ _ _  ___
    | __ \__ \__ \___|| || '_/ _` | ' \(_-<
    |_||_|___/___/    |_||_| \__,_|_||_/__/ MakeDraftCDS.py

usage: MakeDraftCDS.py [-h] -o  [-k] [-s] [-t] [-w] [-g] [-sc] [-ml] [-dm] [-na] [-a] (-i  | -con)

HSS-Trans.MakeDraftCDS.py: Assemble a reduced but high-confidence assembly from public RNA-seq data

optional arguments:
  -h, --help            show this help message and exit
  -o , --output_dir     Directory for data output.
  -k , --kmer_len       Specifies K-mer length (odd integer only) for assembly using Soapdenovo-Trans. K-mer length will be set to 35
                        by default.
  -s , --filesizelimit
                        Specifies the parital download limit/ file size requirement(mb) of accession read files. Limit set to 1500
                        (mb) by default.
  -t , --threads        Total thread pool for workers. Needs to be divisible by number of workers.
  -w , --workers        Specifies the maximum workers for running multiple download-assembly jobs in parallel. Set to 2 by default.
  -g , --gene_code      Genetic code (codon table) passed to ORFfinder during ORF extraction. Set to 1 (universal) by default. Refer
                        to https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi for more information.
  -sc , --start_codon   ORF start codon passed to ORFfinder during ORF extraction. Set to 0 (ATG only) by default. Refer to ORFfinder
                        usage https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/ORFfinder/USAGE.txt for more information
  -ml , --min_len       Minimal ORF length (nt) passed to ORFfinder during ORF extraction. Set to 300 by default.
  -dm , --download_method
                        Method to download accession runs. ftp/ascp. Set to aspera download (ascp) by default.
  -na , --n_accessions
                        Number of single-accession-assemblies to combine in order to generate the preliminary assembly.
  -a , --accessions     User-defined list of SRA run accessions to fetch for preliminary assembly. If insufficient accessions
                        provided, run will be supplemented with other public accessions. E.g.: SRR123456,SRR654321,ERR246810,...
  -i , --id             NCBI TaxID of organism for fetching SRA run accessions.
  -con, --conti         Resume incomplete run based on output directory. Only requires -o to run.
```
**Alternative download methods**\
For high-speed download using IBM Aspera file transfer framework:
```
python3 ./main/MakeDraftCDS.py --output_dir <output directory> -i <NCBI TaxID> -g <Genetic code> -dm ascp
```
For download via FTP using cURL:
```
python3 ./main/MakeDraftCDS.py --output_dir <output directory> -i <NCBI TaxID> -g <Genetic code> -dm ftp
```
**Continuing an interrupted run / Overwriting a previous run in the same directory**

Running the main pipeline (steps 1 and 2) might take quite long especially if there are many accessions to download for the organism of interest or if the user has limited internet bandwidth.
Therefore, a logging system (in a logs.json file) has been incorporated into the pipeline for users to continue an interrupted run using the --conti option.
Simply state the output directory and the pipeline will parse the log file, inherit arguments from the interrupted run and pickup from where it left off.
```
python3 ./main/MakeDraftCDS.py --output_dir <output directory> --conti

python3 ./main/SelectAccessions.py --output_dir <output directory> --conti 
```
To prevent corruption of data files, the pipeline will not allow users to re-run pipeline steps that has previously been initialized in the output directory implicitly.
To re-run Step 1, users must first delete all contents within the output directory recursively:
```
rm <output directory>/* -r
python3 ./main/MakeDraftCDS.py --output_dir <output directory> -i <NCBI TaxID>
```
To re-run Step 2, users can delete data files specific to Step 2 and roll-back logs using the --force option:
```
python3 ./main/MakeDraftCDS.py --output_dir <output directory> --force
``` 

 