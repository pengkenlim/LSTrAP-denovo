# LSTrAP-denovo
 Large Scale Transcriptome Assembly Pipeline
 Note: Pipeline is still under active development.

# Setup
**clone repository to local machine**
```
git clone https://github.com/pengkenlim/LSTrAP-denovo.git
```
**Create environment**
```
cd LSTrAP-denovo
virtualenv -p python3 <MY_ENV>
source ./<MY_ENV>/bin/activate
pip install --upgrade pip
pip install -r ./setup/requirements.txt
```
**Installing programs used in the pipeline**
```
python3 ./setup/install.py
```

# Running the pipeline

**Step 1. Generating a preliminary assembly (reduced but high-confidence assembly)**

Note: Not in active development. All functionalities are fully fleshed-out.

Simplest implementation is as follows
```
python main/Preliminary_assembly.py -o <OUTPUT_DIR> -t <THREAD_POOL> -i <NCBI_TAXID> -dm ftp
```

Full options:
```
usage: Preliminary_assembly.py [-h] -o  [-k] [-ct] [-s] [-t] [-w] [-g] [-sc] [-ml] [-dm] [-na] [-a] (-i  | -con)

LSTrAP-denovo.Preliminary_assembly: Assemble a reduced but high-confidence assembly from public RNA-seq data

optional arguments:
  -h, --help            show this help message and exit
  -o , --output_dir     Directory for data output.
  -k , --kmer_len       Specifies K-mer length (odd integer only) for assembly using Soapdenovo-Trans. K-mer length will be set to 35 by default.
  -ct , --consensus_threshold
                        Specifies consensus threshold during filtering. Threshold will be determined automatically by default.
  -s , --filesizelimit
                        Specifies the parital download limit/ file size requirement(mb) of accession read files. Limit set to 1500 (mb) by default.
  -t , --threads        Total thread pool for workers. Needs to be divisible by number of workers.
  -w , --workers        Specifies the maximum workers for running multiple download-assembly jobs in parellel. Set to 2 by default.
  -g , --gene_code      Genetic code (codon table) passed to ORFfinder during ORF extraction. Set to 1 (universal) by default. Refer to https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi for more information.
  -sc , --start_codon   ORF start codon passed to ORFfinder during ORF extraction. Set to 0 (ATG only) by default. Refer to ORFfinder usage https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/ORFfinder/USAGE.txt for more information
  -ml , --min_len       Minimal ORF length (nt) passed to ORFfinder during ORF extraction. Set to 300 by default.
  -dm , --download_method
                        Method to download accession runs. ftp/ascp.
  -na , --n_accessions
                        Number of single-accession-assemblies to combine in order to generate the preliminary assembly.
  -a , --accessions     User-defined list of SRA run accessions to fetch for preliminary assembly. If insufficient accessions provided, run will be supplemented with other public accessions. E.g.: SRR123456,SRR654321,ERR246810,...
  -i , --id             NCBI TaxID of organism for fetching SRA run accessions.
  -con, --conti         Resume incomplete run based on output directory. Only requires -o to run.
```

**Step 2. Large-scale download , Quality Control and transcriptome-profile-based clustering of accessions** 

Note: Currently in active development

Simplest implementation is as follows
```
<INSERT CLI IMPLEMENTATION HERE>
```


**Step 3. Moduluar concatenation of multi-kmer and transcriptome-profile-specific assemblies.** 

Note: Not yet in development

(organ/condition-specific assembly, remap filtering to get bona fide genes and concatenation)
