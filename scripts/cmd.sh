## link genome.fa and genes.gff to AET_V4 folder under software/snpEff/data
## modify .config file to include 
cd ~/HMW_glutenin/data/
java -jar ~/software/snpEff/snpEff.jar build  -gff3 -v AET_V4 > aet.v4.db.build.out 2> aet.v4.db.build.err &
java -jar ~/software/snpEff/snpEff.jar AET_V4  hmw.glutenin.tauschii.vcf > hmw.glutenin.tauschii.ann.vcf
grep "CHROM" hmw.glutenin.tauschii.ann.vcf | tr '\t' '\n' | tail -n +10 > ids.txt
cat ids.txt| parallel --will-cite -j 5 "samtools faidx /bulk/jpoland/genome/Ae_tauschii/AL8_78/index/GCA_002575655.1_Aet_v4.0_genomic.fna CM008368.1:419306988-419309556 | bcftools consensus hmw.glutenin.tauschii.vcf.gz -s {} | sed 's/>.*/>{}/' > seq/ind/{}.subunit.x.fa"
minimap2 -t 5 -2 -I 20G -K5G -x asm5 arina.chr1D.first100mb.fa aet.v4.chr1D.first100mb.format.id.fa | gzip > 210518.aet.on.arina.1d.paf.gz
cat hmw.glutenin.tauschii.ann.vcf |grep -Fv './.' | ~/owwc/scripts/vcf2long.v1.sh > hmw.glutenin.tauschii.ann.vcflong
cat hmw.glutenin.tauschii.ann.vcf | ~/owwc/scripts/vcf2long.v1.sh | grep -Fv './.'> hmw.glutenin.tauschii.ann.vcflong.correct.grep

##
cd /bulk/lianggao/HMW_glutenin/align-to-insilico/Lancer/
samtools merge -h rg1.txt -R CM008368.1:419304488-419368495 lancer.merged.2.5k.bam  *.bam > merge.Region.out 2> merge.Region.err &
##

cd /bulk/lianggao/HMW_glutenin/var_call/
prefix="210528.insilico.jagger.aet.10genomes.glud1"
crams="/bulk/lianggao/HMW_glutenin/align-to-insilico/bamlist.10genomes.extra.sorted.txt"
ref="/bulk/jpoland/genome/wheat-tauschii/jagger_AL8/index/jagger.v1.1.rm.D.tauschii.AL8.chr1to7.fa"
/homes/lianggao/scripts/Variant_call/call_bcftools_leagacy.zsh $ref $crams $prefix CM008368.1:419304488-419368495

prefix='210528_v2'
/homes/lianggao/HMW_glutenin/scripts/call_bcftools_all.sites.legacy.zsh $ref $crams $prefix CM008368.1:419304488-419368495 ### needs  output gvcf ...?



cd /bulk/lianggao/HMW_glutenin/var_call/210528.insilico.jagger.aet.10genomes.glud1/
bcftools merge -r CM008368.1:419304488-419368495 owwc.raw.legacy.vcf.gz 210528.insilico.jagger.aet.10genomes.glud1_CM008368.1:419304488-419368495.vcf.gz | bgzip > 210528.merged.breadwheat.owwc.vcf.gz
zcat 210528.merged.breadwheat.owwc.vcf.gz | /homes/lianggao/scripts/vcf_filtering/filter_vcf_rm_path.awk -v dphom=2 -v dphet=4 -v minqual=40 -v mindp=100 -v minhomn=1 -v minhomp=0.9 -v tol=0.2 -v minmaf=0.01 -v minpresent=0.1 | bgzip > 210528.merged.breadwheat.owwc.filtered.vcf.gz
tabix -Cp vcf 210528.merged.breadwheat.owwc.filtered.vcf.gz

cd /bulk/lianggao/HMW_glutenin/var_call/210528_v2
bcftools merge -r CM008368.1:419304488-419368495 owwc.raw.legacy.vcf.gz 210528_v2_CM008368.1:419304488-419368495.vcf.gz | bgzip > 210528.merged.breadwheat.owwc.vcf.gz
zcat 210528.merged.breadwheat.owwc.vcf.gz | /homes/lianggao/scripts/vcf_filtering/filter_vcf_rm_path.awk -v dphom=2 -v dphet=4 -v minqual=40 -v mindp=100 -v minhomn=1 -v minhomp=0.9 -v tol=0.2 -v minmaf=0.01 -v minpresent=0.1 | bgzip > 210528.merged.breadwheat.owwc.filtered.vcf.gz
tabix -Cp vcf 210528.merged.breadwheat.owwc.filtered.vcf.gz


cat ids.210530.txt | parallel --will-cite -j 5 "samtools faidx /bulk/jpoland/genome/Ae_tauschii/AL8_78/index/GCA_002575655.1_Aet_v4.0_genomic.fna CM008368.1:419306988-419309556 | bcftools consensus -M N -te -n iqtree -c bioconda  iqtreeH R 210528.merged.breadwheat.owwc.filtered.vcf.gz  -s {} | sed 's/>.*/>{}/'" > 210530.recall.glud1x.M.HR.fa

cat ids.210530.txt | parallel --will-cite -j 5 "samtools faidx /bulk/jpoland/genome/Ae_tauschii/AL8_78/index/GCA_002575655.1_Aet_v4.0_genomic.fna CM008368.1:419364015-419365995 | bcftools consensus -M N -H R 210528.merged.breadwheat.owwc.filtered.vcf.gz  -s {} | sed 's/>.*/>{}.glud1y/'" > 210530.recall.glud1y.M.HR.true.fa

conda create -n iqtree -c bioconda  iqtree
conda activate iqtree
mkdir test.18.dna
iqtree -s 210530.sel18.glud1x.fa -m TEST -bb 1000 -alrt 1000 > test.out 2> test.err &
iqtree -s 210530.recall.glud1x.M.HR.aa.fa -st AA -m TEST -bb 1000 -alrt 1000 > glud1x.runiqtree.out 2> glud1x.runiqtree.err &
iqtree -s 210530.sel18.glud1x.aa.fa -st AA -m TEST -bb 1000 -alrt 1000 > glud1x.runiqtree.out 2> glud1x.runiqtree.err &
iqtree -s 210530.recall.glud1y.M.HR.true.aa.fa -st AA -m TEST -bb 1000 -alrt 1000 > glud1y.runiqtree.out 2> glud1y.runiqtree.err &

## bed.N.C.terminal.R
seqtk subseq 210530.recall.glud1x.M.HR.aa.fa glud1x.N.terminal.aa.bed > 210530.recall.glud1x.M.HR.aa.Nterminal.fa
seqtk subseq 210530.recall.glud1x.M.HR.aa.fa glud1x.CRR.aa.bed > 210530.recall.glud1x.M.HR.aa.CRR.fa
seqtk subseq 210530.recall.glud1x.M.HR.aa.fa glud1x.Cterminal.aa.bed > 210530.recall.glud1x.M.HR.aa.Cterminal.fa
seqtk subseq  210530.recall.glud1y.M.HR.true.aa.fa  glud1y.Cterminal.aa.bed > 210530.recall.glud1y.M.HR.aa.Cterminal.fa
seqtk subseq  210530.recall.glud1y.M.HR.true.aa.fa  glud1y.CRR.aa.bed > 210530.recall.glud1y.M.HR.aa.CRR.fa
seqtk subseq  210530.recall.glud1y.M.HR.true.aa.fa  glud1y.N.terminal.aa.bed > 210530.recall.glud1y.M.HR.aa.N.terminal.fa


## Cys grep
cat  hmw.glutenin.tauschii.ann.vcf | cut -f 8 | tr ";" "\n" | grep ANN= | tr "," "\n" | sed 's/ANN=//' | grep misse | tr "|" "\t" | cut -f1-3,10-15 | column -t| grep -i Cys

perl query_fa_extract_based_on_column.pattern.batch.cmb.pl --fasta_file /homes/lianggao/HMW_glutenin/data/210530.recall.glud1x.M.HR.fa --fasta_file2 /homes/lianggao/HMW_glutenin/data/phylogenetic_trees/glud1y/210530.recall.glud1y.M.HR.true.fa --ids /homes/lianggao/HMW_glutenin/archive/210603.ids.link.vcf.iTOL.txt  --column 0 --chop1 no > ~/HMW_glutenin/data/210604.cmb.glud1x.glud1y.fa &

bcftools filter -r CM008368.1:419306988-419307120,CM008368.1:9154-9556,CM008368.1:419364015-419364147,CM008368.1:419365551-419365995  210528.merged.breadwheat.owwc.filtered.vcf.gz|bgzip > 210528.merged.breadwheat.owwc.filtered.glud1x.glud1y..n.and.c.terminals.vcf.gz &

prefix='210608';mkdir /homes/lianggao/HMW_glutenin/data/$prefix; 
crams='/homes/lianggao/HMW_glutenin/data/210608.bamlist.owwc.bread.wheat.txt'
ref="/bulk/jpoland/genome/wheat-tauschii/jagger_AL8/index/jagger.v1.1.rm.D.tauschii.AL8.chr1to7.fa"
/homes/lianggao/HMW_glutenin/scripts/call_bcftools_all.sites.legacy.zsh $ref $crams $prefix CM008368.1:419304488-419368495

cd /bulk/lianggao/HMW_glutenin/sub_bams/
cat 210608.bamlist.owwc.bread.wheat.txt |  parallel --will-cite -j 10 --plus  'samtools view -hb {} CM008368.1:419304488-419368495 > {/}'
ls *.bam |  parallel --will-cite -j 10 "picard AddOrReplaceReadGroups I={} O={.}_RG.bam RGID={} RGLB=LIB1 RGPL=ILLUMINA RGPU=unit1 RGSM={} > logs/{.}.out 2> logs/{.}.err" &
ls *_RG.bam | parallel --will-cite -j 8 'samtools index -c {}' &
ls *_RG.bam | parallel --will-cite -j 8 "/homes/lianggao/HMW_glutenin/scripts/gatk/hc.sh {} > {/.}.hc.out 2>{/.}.hc.err" &
ls *.g.vcf | parallel --will-cite -j 8 'grep -v NWV {} > gvcfs.cleanedup.chrun/{}'
cd gvcfs.cleanedup.chrun
ls *.g.vcf | grep -v all | xargs -I {} echo '--variant {} \' > variants.txt
sh combineGVCFs.sh > log/combine.out 2>&1 &
sh genotypeGVCFs.sh ##

## some random commands to merge bams and concat same id header fasta
samtools merge -h rg1.txt -R CM008368.1:419304488-419368495 mace.merged.2.5k.bam  *.bam > merge.Region.out 2> merge.Region.err &
seqkit concat phylogenetic_trees/glud1x/210530.recall.glud1x.M.HR.aa.simplify.name.fa phylogenetic_trees/glud1y/210530.recall.glud1y.M.HR.true.aa.simplify.name.fa > 210605.all.glud1x.glud1y.aa.cmb.fa



cat /bulk/lianggao/HMW_glutenin/sub_bams/gvcfs.cleanedup.chrun/all.vcf |  ~/scripts/vcf_filtering/filter_vcf_DP.AD.awk  -v dphom=2 -v dphet=4 -v minqual=40 -v mindp=100 -v minhomn=1 -v minhomp=0.9 -v tol=0.2 -v minmaf=0.01 -v minpresent=0.01 | bgzip > 210610.all.filtered.ad.min01.vcf.gz

