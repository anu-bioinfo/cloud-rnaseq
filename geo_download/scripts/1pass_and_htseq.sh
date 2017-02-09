#!/bin/bash
# replace each of the following with correct paths

STAR=/usr/local/bin/STAR
HTSEQ=/usr/local/bin/htseq-count
SAMTOOLS=/usr/local/bin/samtools
GENOME_DIR=/data/jwang/genome/STAR/HG38/
GENOME_FASTA=/data/jwang/genome/hg38/hg38.fa
SJDB_GTF=/data/jwang/genome/hg38/hg38.gtf
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
--outSAMattributes NH HI NM MD"
#COMMON_PARS="--runThreadN 12 --sjdbGTFfile ${SJDB_GTF} --outSAMtype SAM"
OUTPUT_DIR=$1
READS="${2} ${3}"
#Reads="/home/dobin/STARruns/STARtests/2pass/LID16627_FC61U30AAXX_3_1.txt /home/dobin/STARruns/STARtests/2pass/LID16627_FC61U30AAXX_3_2.txt"

if [ $# -lt 2 ]; then 
  echo "Usage: $0 <outputdir> <read_file_1> [<read_file_2>]"
  exit 1
fi

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
