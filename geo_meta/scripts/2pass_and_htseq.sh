#!/bin/bash
# replace each of the following with correct paths

STAR=/usr/local/bin/STAR
HTSEQ=/usr/local/bin/htseq-count
GENOME_DIR=/data/jwang/genome/STAR/HG38/
GENOME_FASTA=/data/jwang/genome/hg38/hg38.fa
SJDB_GTF=/data/jwang/genome/hg38/hg38.gtf
COMMON_PARS="--runThreadN 12 --outSAMattributes All --genomeLoad LoadAndKeep"
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
echo "Runnning first pass"
mkdir -p Pass1
cd Pass1
$STAR $COMMON_PARS --genomeDir $GENOME_DIR --readFilesIn $READS
cd ..

# make splice junctions database file out of SJ.out.tab, filter out non-canonical junctions
echo "Generating Genome based on first results"
mkdir -p GenomeForPass2
cd GenomeForPass2
awk 'BEGIN {OFS="\t"; strChar[0]="."; strChar[1]="+"; strChar[2]="-";} {if($5>0){print $1,$2,$3,strChar[$4]}}' ../Pass1/SJ.out.tab > SJ.out.tab.Pass1.sjdb

# generate genome with junctions from the 1st pass
$STAR --genomeDir ./ --runMode genomeGenerate --genomeFastaFiles $GENOME_FASTA --sjdbFileChrStartEnd SJ.out.tab.Pass1.sjdb --sjdbOverhang 100 --runThreadN 12
cd ..

# run 2nd pass with the new genome
echo "Running second pass"
mkdir -p Pass2
cd Pass2
$STAR $COMMON_PARS --genomeDir ../GenomeForPass2 --readFilesIn $READS
cd ..

# htseq-count
echo "Running htseq"
$HTSEQ -s no -r pos -f sam -m intersection-nonempty  ./Pass2/Aligned.sam $SJDB_GTF > htseq-count.txt
