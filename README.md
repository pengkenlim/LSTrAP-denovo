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

# Setting up
**clone repository to local machine**
```
$ git clone https://github.com/pengkenlim/HSS-Trans.git
```
**Create environment, install packages and grab the latest version of FFQ**

```
$ cd HSS-Trans
$ virtualenv -p python3 <MY_ENV>
$ source ./<MY_ENV>/bin/activate
$ pip install --upgrade pip
$ pip install -r ./setup/requirements.txt
$ pip install ffq
```
**Install dependencies into programs sub-directory**
```
$ python3 ./setup/install.py
```

# Step 1. Assembling Draft CDSs  
**Simplest implementation**
```
$ python3 ./main/MakeDraftCDS.py --output_dir <output directory> -i <NCBI TaxID> -g <Genetic code>
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

# Step 2. Selecting representative acccessions for transcriptome assembly
**Simplest implementation**
```
python3 ./main/SelectAccessions.py --output_dir <output directory> --accessions_limit 1000
```
Full options:
```
     _  _ ___ ___    _____
    | || / __/ __|__|_   _| _ __ _ _ _  ___
    | __ \__ \__ \___|| || '_/ _` | ' \(_-<
    |_||_|___/___/    |_||_| \__,_|_||_/__/ SelectAccessions.py

usage: SelectAccessions.py [-h] -o  [-ps] [-s] [-t] [-w] [-dm] [-al] [-kr] [-ct] [-clib] [-con] [-f]

HSS-Trans.SelectAccessions.py: Selection of representative accessions for transcriptome assembly. NOTE: This is step 2 of 2 in the
HSS-Trans pipeline. Requires prior run of step 1: MakeDraftCDS.py. Refer to https://github.com/pengkenlim/HSS-Trans for more
information on pipeline usage and implmentation

optional arguments:
  -h, --help            show this help message and exit
  -o , --output_dir     Directory for data output. Directory needs to be same as for step 1 (MakeDraftCDS.py).
  -ps , --pseudoalignment_threshold
                        Specifies reads pseudoaligned (%PS) threshold for quality control. Expression data of accessions that do not
                        meet this threshold will be excluded as features for clustering. Set to 0 by default where %PS threshold will
                        be set to be the lower bound of the %PS distribution (Q1 - 1.5 *IQR) or 20%, whichever is higher.
  -s , --filesizelimit
                        Specifies the size limit(mb) for the partial download of accession read files (gzip compressed) for
                        pseudoalignment. Limit set to 500 (mb) by default. Size maybe decreased to improve the overall runtime
                        (Download and Psudoalignment) and storage used in this step of the pipeline. However, doing so might
                        compromise accurate gene expression quantification as a result of limited sequencing depth.
  -t , --threads        Total thread pool for workers. Needs to be divisible by number of workers.
  -w , --workers        Specify the maximum workers for running multiple download-pseudoalignment jobs in parallel. Set to 4 by
                        default.
  -dm , --download_method
                        Method to download accession runs. ftp/ascp. Set to aspera download (ascp) by default.
  -al , --accessions_limit
                        Specifies the upper limit for number of accessions to download and process. Accessions will be selected from a
                        pre-randomised list that was fetched during MakeDraftCDS.py run and stored in in the logs.json file. Default
                        set to 500.
  -kr , --k_range       Specifies the range of k (number of clusters) to iterate through during clustering. Lower and upper limit
                        seperated by colon(:). Set to auto(5:20) by default. Optimal k will be chosen from this range based on
                        silhouette coefficient, a clustering performance metric. As such, please set range within expectations based
                        on heterogenity of expression data for that organism.
  -ct , --consensus_threshold
                        Specifies consensus threshold of preliminary assembly. Default set to 0 where optimal threshold determined
                        automatically in step 1 will be used.
  -clib , --cluster_lib_size
                        Specifies the minimum library size (mb) for each cluster to guide the selection of representative accessions
                        to download. Total file sizes of read files (gzipped) to download for each cluster is used as an approximation
                        of libary size instead of number of reads.
  -con, --conti         Resume incomplete run based on output directory. Only requires -o to run.
  -f, --force           Delete data from previous SelectAccessions.py run and start a fresh SelectAccessions.py in output directory.
```

# Continuing an interrupted run in the same directory
**Continuing an interrupted run**

Running the main pipeline (steps 1 and 2) might take quite long especially if there are many accessions to download for the organism of interest or if the user has limited internet bandwidth.
Therefore, a logging system (in a logs.json file) has been incorporated into the pipeline for users to continue an interrupted run using the `--conti` option.
Simply state the output directory and the pipeline will parse the log file, inherit arguments from the interrupted run and pickup from where it left off.
```
$ python3 ./main/MakeDraftCDS.py --output_dir <output directory> --conti

$ python3 ./main/SelectAccessions.py --output_dir <output directory> --conti 
```
**Overwriting a previous run in the same directory**

To prevent corruption of data files, the pipeline will not allow users to re-run pipeline steps that has previously been initialized in the output directory implicitly.
To re-run Step 1, users must first delete all contents within the output directory recursively:
```
$ rm <output directory>/* -r
$ python3 ./main/MakeDraftCDS.py --output_dir <output directory> -i <NCBI TaxID>
```
To re-run Step 2, users can delete data files specific to Step 2 and roll-back logs using the `--force` option:
```
$ python3 ./main/MakeDraftCDS.py --output_dir <output directory> --force
``` 

# Data downloading
HSS-Trans downloads RNA-seq data directly from the servers of the European Nucleotide Archive (ENA). This is favoured over other popular methods of fetching fastq files through fastq-dump / fasterq-dump (SRA-tools) as it facillitates the partial (truncated) download of deep-sequencing accessions in compressed (gunzipped) FASTQ files.

**Alternative download methods**

HSS-Trans can be download RNA-seq FASTQ files in one of two methods. In the event that one method doesn't work due to reasons server-side, users are encouraged to use the other alternative. Download method can be specified using the `--download_method` option.

For high-speed download using IBM Aspera file transfer framework:
```
$ python3 ./main/MakeDraftCDS.py --output_dir <output directory> -i <NCBI TaxID> -g <Genetic code> --download_method ascp

$ python3 ./main/SelectAccessions.py --output_dir <output directory> --download_method ascp
```
For download via FTP using cURL:
```
$ python3 ./main/MakeDraftCDS.py --output_dir <output directory> -i <NCBI TaxID> -g <Genetic code> --download_method ftp

$ python3 ./main/SelectAccessions.py --output_dir <output directory> --download_method ftp
```
 
**Download Parallelization**

In HSS-Trans, downloading and processing of accessions are coupled into individual jobs and spawned as parallel worker processes. This coupling allows for downloaded data to be immediately fed into compute-heavy tasks (assembly in MakeDraftCDS.py / pseudoalignment in SelectAccessions.py) without waiting for all downloads to finish. In the event where the user has high download speed, this parallelization can allow user to overcome connection bandwidth limitations through multiple concurrent download connections.
In both steps of the pipeline, the number of workers for parallelization can be specified with the `--workers` option.

**Limit the number of accessions to download and process in step 2**

RNA-seq accessions of the specified taxa are fetched from ENA. However the number of RNA-seq accessioons for some taxa (especially model species; >80k for *A. thaliana*) might be too large to be practically processed by the pipeline. As such, step 2 of the pipeline is limited to process up to 500 accessions (chosen randomly) by default. Users can adjust this upper limit using the `--accessions_limit` option

# File size options 

Note: integers arguments after `--filesizelimit` and `--cluster_lib_size` correspond to the size in megabytes (mb) of gunzipped FASTQ files. 

**Step 1. (MakeDraftCDS.py): Specifying required library size for Single-sample assemblies (SSAs)**

To make sure that Single-sample assemblies are assembled from RNA-seq libraries of similar sequencing depth, RNA-seq gunzipped FASTQs are partially downloaded to the same file size before assembly. Using the `--filesizelimit` option, user can specify the file size to partially download. In addition, this limit also serves the minimal size requirement for accessions to be considered for SSA. File size to be used depends on the transcriptome complexity of the taxa. Refer to DOI for more details.


**Step 2. (SelectAccessions.py): Specifying library size limit for gene expression estimation**

Due to the large size of some RNA-seq libraries, RNA-seq gunzipped FASTQs are partially downloaded to the same file size before pseudoalignment using Kallisto to improve runtime. Using the `--filesizelimit` option, user can specify the file size to partially download. Accessions that do not exceed the specifed limit will be downloaded in full.

**Step 2. (SelectAccessions.py): Specifying total library size of representative accessions to download for each cluster**

# Quality control and clustering options for Step 2 (SelectAccessions.py)
**Pseudoalignment rate for quality control of RNA-seq accessions**

`--pseudoalignment_threshold`

**Number of clusters for k-medoids clustering**

Mention k-range and silhouette coefficient

``--k_range``



# After running HSS-Trans (Step 1 and Step 2) 

**interpreting report**

<insert link to example HTML file>

**Assembling *de novo* transcriptome using downloaded accessions from the pipeline**

Helper script for running trinity

Install instructions

two options: vanilla assembly and over-assembly (explain)


