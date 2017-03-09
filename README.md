# HCA
  Human Cell Atlas Repo on Chan Zuckerberg Science Initiative
## About This Repo
  This repo serves two purposes.

 1. **AWS STAR**: The first purpose is to provide tools for bio researchers to run the STAR assembly pipeline in a parallee
lized fashion on AWS cloud. The input is the fastq files and the output is the mapped BAM files and gene expression countss
 (from running htseq-count).

 2. **NCBI SRA Pipeline**: The second purpose is to download publicly available scRNA-seq data from the [NCBI SRA archive](https://www.ncbi.nlm.nih.gov/sra), put them on Amazon S3, use tools from 1 to run the assembly pipeline and generate genee
 cell table for all these public datasets. The rawdata, metadata and outputs for all these datasets are made accessible foo
r everyone through Amazon S3.

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
 2. Clone this [repo](https://github.com/chanzuckerberg/hca).
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


Save the script and then run it under the ```hca/batch``` directory. This would launch a bunch of AWS batch jobs after a few minutes. 

You can then drink some coffee, watch and wait the jobs to finish. If you are really anxious, you can see the list of your jobs by running.  

```aegea batch ls -w99```

You can also watch individual jobs by typing: 

```aegea batch watch <job_id>```

See the [aegea documentation](https://github.com/kislyuk/aegea) for more info. 
### Generating the gene-cell table
Once all the jobs are done, you can generate the gene-cell table as follows: 

```
cd hca/geo_metadata/; ./generate_gene_cell_table.py -s <S3_BUCKET> -d <EXPERIMENT_ID> -f <OUT_HTSEQ_CSV> -l <OUT_LOG_CSV> 

```

The output htseq-count file would be in  ```<OUT_HTSEQ_CSV>``` and the STAR logs would be in ```<OUT_LOG_CSV>```. 


## NCBI SRA Pipeline

Coming soon~~~

