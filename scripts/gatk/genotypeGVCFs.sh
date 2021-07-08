#!/bin/bash

REFFASTA=/bulk/jpoland/genome/Ae_tauschii/AL8_78/index/Aet_v4.fa
#REFFASTA=/bulk/jpoland/genome/wheat-tauschii/jagger_AL8/index/jagger.v1.1.rm.D.tauschii.AL8.chr1to7.fa
GATKDIR=/homes/lianggao/software/gatk-4.2.0.0
TMP=/bulk/lianggao/tmp

export PATH=$GATKDIR:$PATH


# the gVCFG files obyained before need to be combined to be used with GenotypeGVCFs tool

REGION=CM008368.1


echo Joint genotyping started
date

gatk  GenotypeGVCFs \
   -R $REFFASTA \
   --tmp-dir $TMP \
   -V all.317.vcf \
   -O all.vcf \
   -stand-call-conf 5

echo Run ended
date

