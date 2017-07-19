# scRNA-seq
  scRNA-seq Repo on Chan Zuckerberg Science Initiative
## About This Repo
  This repo serves two purposes.

 1. **AWS STAR**: The first purpose is to provide tools for bio researchers to run the STAR alignment pipeline in a parallelized fashion on AWS cloud. The input is the fastq files and the output is the genome-mapped BAM files and gene expression counts (from running htseq-count).

 2. **NCBI SRA Pipeline**: The second purpose is to download publicly available scRNA-seq data from the [NCBI SRA archive](https://www.ncbi.nlm.nih.gov/sra), put them on Amazon S3, use tools from 1 to run the assembly pipeline and generate gene
 cell table for all these public datasets. The rawdata, metadata and outputs for all these datasets are made accessible for everyone through Amazon S3.

## AWS STAR
### Data Requirements
  In order to use our tools,  you need to have AWS S3 access and have all your data under a primary S3 path (i.e. s3://czi-hca/data). Under this path, gzipped fastq files for a sample should be placed under ``` <S3_BUCKET>/<EXPERIMENT_ID>/rawdata/<SAMPLE_ID>/```

```
,where
    <S3_BUCKET>: the primary S3_PATH
<EXPERIMENT_ID>: experiment id of your choosing
    <SAMPLE_ID>: sample id under the experiment
```

### Tool Installation
 1. Install [aegea](https://github.com/kislyuk/aegea). This is the tool we use to launch the STAR jobs on the AWS cloud.
 2. Clone this [repo](https://github.com/chanzuckerberg/cloud-rnaseq).
 3. Generate dockerimage for our STAR jobs.

 It looks like the following in code:

```
git clone https://github.com/kislyuk/aegea.git
cd aegea; make; make install; cd ..;
# Read the documentation if the above failed.
git clone https://github.com/chanzuckerberg/hca.git
cd hca/batch; make dockerimage;  cd ../..;
```


### Running the pipeline
After installation, here is how you run your pipeline. Right now we support two reference genomes hg38-plus (for homo sapiens) and mm10-plus (for mus musculus). The -plus part is the spiked ERCC RNA references so if you spike your data, you should be able to map those reads. You also have to decide how many jobs you want to run at the same time. If you have 400 samples and 10 jobs, each job would run through 40 samples.

Now, run the following command

```cd hca/batch; ./aws_star.sh ```

The above will give you the usage of this command:

```
Usage: ./aws_star.sh <S3_BUCKET> <EXP_IDS> <NUM_PARTITIONS> <TAXON>
   <S3_BUCKET>: S3 Path. Example: s3://czi-hca/data
   <EXP_IDS>: experiment directories under <S3_BUCKET> seperated by ','. Example: 200067835,200057872
   <NUM_PARTITIONS>: Number of partitions you want to run. Example: 10 means you will run in 10 different jobs
   <TAXON>: mus or homo. mus will use mm10 as reference genome, homo will use hg38 as reference genome
```

Now, if you are ready you can run a real command such as

``` ./aws_star.sh s3://czi-hca/data 200067835,200057832 3 homo ```

This would generate a shell script like the following:

```
aegea batch submit --execute sra_download/run_star_and_htseq.py --storage /mnt=500 --ecr-image sra_download --environment S3_BUCKET=s3://czi-hca/data EXP_IDS=200067835,200057832 NUM_PARTITIONS=3 TAXON=homo PARTITION_ID=0 --memory 64000
sleep 100
aegea batch submit --execute sra_download/run_star_and_htseq.py --storage /mnt=500 --ecr-image sra_download --environment S3_BUCKET=s3://czi-hca/data EXP_IDS=200067835,200057832 NUM_PARTITIONS=3 TAXON=homo PARTITION_ID=1 --memory 64000
sleep 100
aegea batch submit --execute sra_download/run_star_and_htseq.py --storage /mnt=500 --ecr-image sra_download --environment S3_BUCKET=s3://czi-hca/data EXP_IDS=200067835,200057832 NUM_PARTITIONS=3 TAXON=homo PARTITION_ID=2 --memory 64000
sleep 100
```


Save the script and then run it under the ```cloud-rnaseq/batch``` directory. This would launch a bunch of AWS batch jobs after a few minutes.

You can then drink some coffee, watch and wait the jobs to finish. If you are really anxious, you can see the list of your jobs by running.  

```aegea batch ls -w99```

You can also watch individual jobs by typing:

```aegea batch watch <job_id>```

See the [aegea documentation](https://github.com/kislyuk/aegea) for more info.

The results will be accessible under ```<S3_BUCKET>/<EXPERIMENT_ID>/results/``` once the pipeline is successfuly run. Each sample would generate a <SAMPLE_ID>.tgz under this directory. The .tgz file includes the htseq-count output, the unmapped reads and the sorted BAM file.


### Generating the gene-cell table
Once all the jobs are done, you can generate the gene-cell table as follows:

```
cd hca/geo_meta/; ./generate_gene_cell_table.py -s <S3_BUCKET> -d <EXPERIMENT_ID> -f <OUT_HTSEQ_CSV> -l <OUT_LOG_CSV>

```

The output htseq-count file would be in  ```<OUT_HTSEQ_CSV>``` and the STAR logs would be in ```<OUT_LOG_CSV>```.

### Notes
  This pipeline doesn't work with 10X genomics data because SRA format doesn't store 10X data correctly. Please go to https://support.10xgenomics.com/single-cell-gene-expression/datasets to download 10X datasets.

## NCBI SRA Pipeline

### How to obtain the gene cell table from the NCBI dataset?

We downloaded about 185 single-cell RNA sequencing SRR datasets from the NCBI database. We converted the data into fastq files and ran STAR pipeline on the samples downloaded.
We only download datasets that contain Homo sapiens and Mus musculus samples. For datasets that contain *Homo sapiens* samples, we ran STAR againt `hg38` for all samples in the datasets. Similarly for datasets that contain *Mus musculus* samples, we ran STAR against `mm10` for all samples in the datasets. We then run htseq-count to get the gene count for each single-cell sample. Finally, we ran a harvest script to aggregate the results into a few files for your ease of use.


The data is currently located under s3://czi-hca/summary-data/20170426/ AWS. This is a publicly readable AWS directory.


The actual (summary) data is directly downloadable at
- https://s3.amazonaws.com/czi-hca/summary-data/20170426/20170426.homo.tgz
- https://s3.amazonaws.com/czi-hca/summary-data/20170426/20170426.mus.tgz

Some interesting statistics about this dataset:

  * Raw data size: 80+TB
  * Summary data size: 1.5GB
  * Sample Size:
    * Homo sapiens: 19416 processed samples out of  51 datasets
    * Mus musculus: 43655 processed samples out of 133 datasets

### How to use this dataset?
5 relevant files are generated after you uncompress the .tar.gz file.
  1. `*.htseq-count.csv`: gene cell table in csv format. The first column is the gene names. The first few rows include some mapping information that would be of your interest. The first row contains SRR id from NCBI. The second row contains the experiment ID from NCBI. The third row indicates taxon. The fourth row contains the GSM ID. The Experiment ID and the GSM ID can be used to mapped to the metadata described below.
  2. `*.log.csv.gz`: STAR log files for each sample in csv format. The first column is the description of the STAR output fields. It includes # reads, #/% of mapped reads, etc by sample. Again, the first few rows include mapping information to the metadata
  3. `*.metadata.series.csv`: Important metadata fields by dataset in csv format. EXP_ID row in htseq-count references doc_id column in this file.  
  4. `*.metadata.samples.csv`: Important metadata fields by cell in csv format. GSM row in htseq-count references geo_accession column in this file.  
  5. `*.metadata.tgz`: compressed full metadata dataset, once you decompressed it ("tar xvfz <filename>"), you would find a directory structure of <DIR>/<EXPERIMENT_ID>.metadata.json. You can extract full metadata for a sample by doing the following:
      # Python Code below. (exp_id, gsm_id): Experiment ID and GSM ID as described in 1.
      
      ```python
      import json
      meta_file = open("%s.metadata.json" % exp_id.rstrip(), 'rb')
      metadata_hash = json.loads(meta_file.read())
      metadata = metadata_hash[gsm_id]
      ```
      
    `metadata` is structured as follows:
      * `metadata['series_data']`: dataset(experiment) specific data. same for samples in the same dataset
      * `metadata['platform_data']`: sequencing platform specific data. same for samples in the same dataset
      * `metadata['sample_data']`:  sample specific metaadata

### Appendix
 * STAR pipeline script: script used from fastq files to htseq-count
   https://github.com/chanzuckerberg/cloud-rnaseq/blob/master/batch/sra_download/run_star_and_htseq.py
 * Harvest script: script to aggregate all the samples to a few big files
   https://github.com/chanzuckerberg/cloud-rnaseq/blob/master/util/big_gc_table.py




