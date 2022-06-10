# HSS-Trans
 **High-throughput Sample Selection for Transcriptome assembly**  
 An automated pipeline to quality control and select public RNA-seq accessions for transcriptome assembly of species without reference genomes using unsupervised machine learning
 
 **Features**  
 - Simplest implementation of HSS-Trans only requires taxonomic id for species of interest
 - Construction of high confidence draft coding sequences (CDSs) using an efficient *de novo* assembler (SOAPdenovo-Trans) coupled with a novel consensus-based approach to retain common CDSs assembled independently in multiple RNA-seq accessions
 - Rapid identification of low quality RNA-seq accessions (quality control) without the need of reference genome annotations
 - High-speed download of RNA-seq data using IBM Aspera file transfer framework. Option to download using FTP is also available.
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
**Create environment**
```
cd HSS-Trans
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
The pipeline is split into two steps.  
**Step 1. Assembling Draft CDSs (reduced but high-confidence assembly)**



Simplest implementation is as follows
```

```

Full options:
```

```

