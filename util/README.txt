What is the data in the directory? 
We downloaded about 180 single-cell RNA sequencing SRR datasets from the NCBI database. We converted the data into fastq files and run STAR pipeline on the samples downloaded. 
We only download datasets that contain Homo sapiens and Mus musculus samples. For datasets that contain homo samples, we ran STAR againt HG38 for all samples in the datasets. Similaryly for datasets that contain mus samples, we ran STAR against MM10 for all samples in the datasets. We then run htseq-count to get the gene count for each single-cell sample. Finally, we ran a harvest script to aggregate the results into a few files for your ease of use. 

Some interesting statistics about this dataset: 
  * Raw data size: 60+ TB
  * Summary data size: 1.5GB
  * Sample Size:
    ** Homo sapiens: 21968 processed samples out of  50 datasets
    ** Mus musculus: 47102 processed samples out of 129 datasets

How to use this data ? 
3 relevant files are generated for Homo sapiens and Mus musculus respectively as follows:
  1. [mus|homo].mus.htseq-count.csv.gz: gzipped gene cell table in csv format. The first column is the gene names. The first few rows include some mapping information that would be of your interest. The first row contains SRR id from NCBI. The second row contains the experiment ID from NCBI. The third row indicates taxon. The fourth row contains the GSM ID. The Experiment ID and the GSM ID can be used to mapped to the metadata described below. 
  2. [mus|homo]_meta.tgz: compressed metadata file, once you decompressed it ("tar xvfz <filename>"), you would find a directory structure of [mus|homo]_meta/<EXPERIMENT_ID>.metadata.json. You can extract metadata for a sample by doing the following:
      # Python Code below. (exp_id, gsm_id): Experiment ID and GSM ID as described in 1. 
      import json
      meta_file = open("%s.metadata.json" % exp_id.rstrip(), 'rb')
      metadata_hash = json.loads(meta_file.read())
      metadata = metadata_hash[gsm_id]
    metadata is structured as follows: 
      * metadata['series_data']: dataset(experiment) specific data. same for samples in the same dataset 
      * metadata['platform_data']: sequencing platform specific data. same for samples in the same dataset
      * metadata['sample_data']:  sample specific metaadata
   3. [mus|homo].log.csv.gz: gzipped STAR log files for each sample in csv format. The first column is the description of the STAR output fields. It includes # reads, #/% of mapped reads, etc by sample. Again, the first few rows include mapping information to the metadata

Appendix:
 * STAR pipeline script: script used from fastq files to htseq-count
   https://github.com/chanzuckerberg/cloud-rnaseq/blob/master/batch/sra_download/run_star_and_htseq.py
 * Harvest script: script to aggregate all the samples to a few big fiels
   https://github.com/chanzuckerberg/cloud-rnaseq/blob/master/util/big_gc_table.py 

