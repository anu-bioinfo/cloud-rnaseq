#!/bin/bash -ex
# replace each of the following with correct paths

REF_PACK=s3://chapmanb/genomes_automated_test.tar.gz
aws s3 cp $REF_PACK - | tar xvz -C /mnt

STAR=STAR
HTSEQ=htseq-count
SAMTOOLS=samtools
GENOME_DIR=/mnt/genomes/hg19/star
GENOME_FASTA=/mnt/genomes/hg19/snpeff/hg19/sequences.fa
SJDB_GTF=/mnt/genomes/hg19/snpeff/hg19/genes.gtf
COMMON_PARS="--runThreadN 12 --sjdbGTFfile ${SJDB_GTF} --outFilterType BySJout \
--outFilterMultimapNmax 20 \
--alignSJoverhangMin 8 \
--alignSJDBoverhangMin 1 \
--outFilterMismatchNmax 999 \
--outFilterMismatchNoverLmax 0.04 \
--alignIntronMin 20 \
--alignIntronMax 1000000 \
--alignMatesGapMax 1000000 \
--outSAMstrandField intronMotif \
--outSAMtype BAM Unsorted \
--outSAMattributes NH HI NM MD \
--readFilesCommand zcat"
#COMMON_PARS="--runThreadN 12 --sjdbGTFfile ${SJDB_GTF} --outSAMtype SAM"
#if [ $? -eq 0 ]; then touch z.txt; fi

OUTPUT_DIR=/mnt
#READS1=$1
#READS2=$2
aws s3 cp $READS1 /mnt/
aws s3 cp $READS2 /mnt/
READS="/mnt/$(basename $READS1) /mnt/$(basename $READS2)"

#if [ $# -lt 2 ]; then 
#  echo "Usage: $0 <outputdir> <read_file_1> [<read_file_2>]"
#  exit 1
#fi

echo "Output dir: ${OUTPUT_DIR}"
echo "Input files: ${READS}"

# cd to output dir
mkdir -p ${OUTPUT_DIR}
cd ${OUTPUT_DIR}

# run 1st pass
echo "Runnning STAR"
mkdir -p Pass1
cd Pass1
$STAR $COMMON_PARS --genomeDir $GENOME_DIR --readFilesIn $READS
cd ..

# sam tools sort
$SAMTOOLS sort -n -m 6000000000 -o Pass1/Aligned.out.sorted.bam Pass1/Aligned.out.bam

# htseq-count 
echo "Running htseq"
$HTSEQ -s no -f bam -m intersection-nonempty  ./Pass1/Aligned.out.sorted.bam $SJDB_GTF > htseq-count.txt

ls /mnt
