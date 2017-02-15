dockerimage:
	aegea-build-image-for-mission --mission-dir . sra_download sra_download
launch_job_with_image:
	GDS_ID=200067835 NUM_PARTITIONS=10 PARTITION_ID=0 RDB_S3_PATH=s3://czi-hca/config/dump.rdb.gz ./sra_download.py
