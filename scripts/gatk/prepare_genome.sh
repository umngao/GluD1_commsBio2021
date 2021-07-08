#!/bin/bash

# index reference genome for bwa, create fasta indexes (fai and dict)

TMP=/bulk/lianggao/tmp
GATKDIR=/homes/lianggao/software/gatk-4.2.0.0
export PATH=$GATKDIR:$PATH

ref='/bulk/jpoland/genome/wheat-tauschii/jagger_AL8/index/jagger.v1.1.rm.D.tauschii.AL8.chr1to7.fa'
# Genome summary files needed and by GATK tools
gatk CreateSequenceDictionary -R $ref -O $ref.dict
#samtools faidx genome.fa
#
## index for BWA alignment
#bwa index genome.fa
#
## index image file needed by some Spark-based tools (if used)
#gatk --java-options "-Djava.io.tmpdir=$TMP" BwaMemIndexImageCreator \
#     -I genome.fa \
#     -O genome.fa.img
