#!/bin/bash

REFFASTA=/bulk/jpoland/genome/Ae_tauschii/AL8_78/index/Aet_v4.fa
#REFFASTA=/bulk/jpoland/genome/wheat-tauschii/jagger_AL8/index/jagger.v1.1.rm.D.tauschii.AL8.chr1to7.fa
GATKDIR=/homes/lianggao/software/gatk-4.2.0.0
TMP=/bulk/lianggao/tmp

export PATH=$GATKDIR:$PATH

#export JAVA_HOME=/usr/local/jdk1.8.0_121
#export PATH=$JAVA_HOME/bin:$PATH

# We will run the genotyping on one chromosome only.
# Other chromosomes clould be handlen in separate runs, 
# possibly in parallel..

REGION="CM008368.1:419304488-419368495"

ACC=$1

# multi-threading does not work well here - they recommend using Queue to
# parallelize over regions, or just manually run a few in parallel...

echo Genotyping for $ACC started
date

gatk --java-options "-Xmx4g"  HaplotypeCaller \
      --tmp-dir $TMP \
     -R $REFFASTA \
     -I $ACC  \
     -L $REGION \
     -ERC GVCF \
     --native-pair-hmm-threads 4 \
     --minimum-mapping-quality 10 \
     -O $ACC.g.vcf

echo Run ended
date

