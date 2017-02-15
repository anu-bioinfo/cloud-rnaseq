dockerimage:
	aegea-build-image-for-mission --mission-dir . starline starline
launch_job_with_image:
	aegea batch submit --execute 1pass_and_htseq.sh --storage /mnt=500 --ecr-image starline --environment READS1=s3://science-sample-datasets/ncbi-sra-rnaseq/SRR1813952_1.fastq.gz READS2=s3://science-sample-datasets/ncbi-sra-rnaseq/SRR1813952_2.fastq.gz --watch --memory 24000
