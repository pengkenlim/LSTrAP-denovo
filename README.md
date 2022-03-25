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
virtualenv -p python <MY_ENV>
source ./<MY_ENV>/bin/activate
pip install --upgrade pip
pip install -r ./setup/requirements.txt
```
**Installing programs used in the pipeline**
```
python3 ./setup/install.py
```

# Running the pipeline

**Generating a preliminary assembly**

Note: this functionality is still under development

Simplest implementation is as follows
```
python main/preliminary_assembly.py -o <OUTPUT_DIR> -t <THREAD_POOL> -i <NCBI_TAXID> -dm ftp
```

Full options:
```
usage: preliminary_assembly.py [-h] -o  [-k] [-ct] [-s] -t  [-w] [-g] [-sc] [-ml] -dm  (-i  | -a )

LSTrAP-denovo.preliminary_assembly: Assemble a reduced but high-confidence assembly from public RNA-seq data

optional arguments:
  -h, --help            show this help message and exit
  -o , --output_dir     Directory for data output.
  -k , --kmer_len       Specify K-mer length (odd integer only) for assembly using Soapdenovo-Trans. K-mer length will be set to 35 by
                        default.
  -ct , --consensus_threshold
                        Specify consensus threshold during filtering. Threshold will be determined automatically by default.
  -s , --filesizelimit
                        Specify the size limit(mb) of accession read files to partially download. Limit set to 1500 by default.
  -t , --threads        Total thread pool for workers. Needs to be divisible by number of workers.
  -w , --workers        Specify the maximum workers for running multiple download-assembly jobs in parellel. Set to 2 by default.
  -g , --gene_code      Genetic code (codon table) passed to ORFfinder during ORF extraction. Set to 1 (universal) by default. Refer to
                        https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi for more information.
  -sc , --start_codon   ORF start codon passed to ORFfinder during ORF extraction. Set to 1 (ATG only) by default. Refer to ORFfinder usage
                        https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/ORFfinder/USAGE.txt for more information
  -ml , --min_len       Minimal ORF length (nt) passed to ORFfinder during ORF extraction. Set to 300 by default.
  -dm , --download_method
                        Method to download accession runs. ftp/ascp.
  -i , --id             NCBI TaxID of organism for fetching SRA run accessions.
  -a , --accessions     User-defined list of SRA run accessions to fetch. Requires at least 10 accessions. E.g.:
                        SRR123456,SRR654321,ERR246810,...
```

**Place holder for next step** 

**(sample download, QC and k-means classification)**

Place holder description

**Place holder for next step** 

**(organ/condition-specific assembly, remap filtering to get bona fide genes and concatenation)**

Place holder description