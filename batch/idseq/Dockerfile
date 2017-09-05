# Set the base image
FROM ubuntu:16.04

# File Author / Maintainer
MAINTAINER Yun-Fang Juan yunfang@chanzuckerberg.com

# Add packages, update image, and clear cache
RUN apt-get update && apt-get install -y build-essential python-scipy python-redis gdebi-core

RUN pip install htseq==0.6.1p1
RUN pip install redis biopython pysam

WORKDIR /tmp
# install STAR
RUN curl -L https://github.com/alexdobin/STAR/archive/2.5.2b.tar.gz | tar xz
RUN mv STAR-2.5.2b/bin/Linux_x86_64_static/* /usr/local/bin
RUN rm -rf STAR-2.5.2b"

# install Samtools
RUN curl -L https://sourceforge.net/projects/samtools/files/samtools/1.3/samtools-1.3.tar.bz2/download > samtools-1.3.tar.bz2
RUN tar xjf samtools-1.3.tar.bz2; cd samtools-1.3; make; make install; samtools --version

# Compile and install bowtie
RUN wget https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.3.2/bowtie2-2.3.2-source.zip
RUN unzip bowtie2-2.3.2-source.zip
WORKDIR /tmp/bowtie2-2.3.2
RUN make NO_TBB=1 # we may want tbb (for "superior thread scaling") later, so will need to find a way to compile without "NO_TBB=1"
RUN mv -t /usr/local/bin/ bowtie2 bowtie2-align-l bowtie2-align-s bowtie2-build bowtie2-build-l bowtie2-build-s bowtie2-inspect bowtie2-inspect-l bowtie2-inspect-s

WORKDIR /tmp

# Cleanup
RUN rm -rf /tmp/*

WORKDIR /
