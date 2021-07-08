library(data.table)
# install.packages('gggenes')
library(ggplot2)
library(gggenes)

# ### I. first version of vcf and vcf annotation file
# fread('grep -v "#" /homes/lianggao/HMW_glutenin/data/hmw.glutenin.tauschii.ann.vcf | cut -f 1-9') -> vcf.cut
# cbind(vcf.cut[,1:2],vcf.cut[,tstrsplit(V8,'\\|')][,2]) -> vcf.cut.split.eff
# setnames(vcf.cut.split.eff, c('chr','pos', 'eff'))
# vcf.cut.split.eff[,.(chr=1, start=pos-1, end=pos, eff)] -> vcf.cut.split.eff.v2
# data.table(molecule='genome1',gene='Glu-D1x', start=419306988, end=419309556, strand='reverse', orientation=-1) -> test.gene
# data.table(molecule='genome1',gene='Glu-D1y', start=419364015, end=419365995, strand='reverse', orientation=-1) -> test.gene2
# vcf.cut.split.eff.v2[start %between% c(test.gene$start, test.gene$end)] -> vcse1
# vcf.cut.split.eff.v2[start %between% c(test.gene2$start, test.gene2$end)] -> vcse2
# data.table(test.gene, from=vcse1$start, to=vcse1$end, eff=vcse1$eff) -> test.gene.sub; test.gene.sub[,subgene:=paste0('segX',1:.N)]
# data.table(test.gene2, from=vcse2$start, to=vcse2$end, eff=vcse2$eff) -> test.gene.sub2; test.gene.sub2[,subgene:=paste0('segY',1:.N)]
# rbind(test.gene, test.gene2) -> glu.d1.gene
# glu.d1.gene[grepl('D1y', gene), molecule:='genome2']
# rbind(test.gene.sub, test.gene.sub2) -> glu.d1.snps
# glu.d1.snps[grepl('D1y', gene), molecule:='genome2']
# saveRDS(list(glu.d1.gene,glu.d1.snps), file='glu.d1.gene.snps.RDS')
# 
# ggplot(glu.d1.gene, aes(xmin = start, xmax = end, y = molecule, fill = gene, forward = orientation)) +
#     geom_gene_arrow(arrow_body_height = grid::unit(10,"mm"),
#                     arrowhead_height = grid::unit(15,'mm'),
#                     arrowhead_width = grid::unit(15,'mm')) +
#     facet_wrap(~ molecule, scales = "free", ncol = 2)+
#     scale_fill_brewer(palette = "Set3")+
#     geom_subgene_arrow(data = glu.d1.snps,
#                        arrow_body_height = grid::unit(10,"mm"),
#                        arrowhead_height = grid::unit(15,'mm'),
#                        arrowhead_width = grid::unit(15,'mm'),
#                        aes(xmin = start, xmax = end, y = molecule, fill = gene,
#                            xsubmin = from, xsubmax = to), color='black') +
#     geom_subgene_arrow(data = glu.d1.snps[grepl('missen', eff)],
#                        arrow_body_height = grid::unit(10,"mm"),
#                        arrowhead_height = grid::unit(15,'mm'),
#                        arrowhead_width = grid::unit(15,'mm'),
#                        aes(xmin = start, xmax = end, y = molecule, fill = gene,
#                            xsubmin = from, xsubmax = to), color='#cc4c02') +
#     theme_genes() + theme(axis.text=element_text(size=16, face = 'bold', angle=90), axis.text.y=element_blank())
# dev.off()
# 
# fread('/homes/lianggao/HMW_glutenin/data/hmw.glutenin.tauschii.ann.vcflong.correct.grep') -> g
# setnames(g, c("chr", "pos", "ref", "alt", "qual", "sample.id", "gt", "dp", "dv"))
# g[gt == '0/0', gt := paste(0)]
# g[gt == '0/1', gt := paste(1)]
# g[gt == '1/1', gt := paste(2)]
# g[sample.id %in% c('Stanley_5+10'), .(chr, pos, sample1=sample.id, gt1=gt)] -> a
# g[sample.id %in% c('TA10103', 'TA10104','TA2582', 'TA2580', 'TA10928'), .(chr, pos, sample2=sample.id, gt2=gt)] -> b
# a[b, on=c("chr", "pos"), allow.cartesian=T][, d := abs(as.integer(gt1)-as.integer(gt2))][]->ab
# copy(ab) -> ab.stanley.against.l3
# dcast(ab.stanley.against.l3[, .N, key=.(sample1, sample2, d)], sample1 + sample2 ~ d, value.var="N", fill=0) -> y.stanley
# vcf.cut.split.eff[ pos %in% ab.stanley.against.l3[d>0, pos]]
# vcf.cut.split.eff[ pos %in% ab.stanley.against.l3[sample2 =='TA2580' & d>1, pos]]
# 
# 
# g[, .(chr, pos, sample1=sample.id, gt1=gt)] -> a
# g[, .(chr, pos, sample2=sample.id, gt2=gt)] -> b
# a[b, on=c("chr", "pos"), allow.cartesian=T][, d := abs(as.integer(gt1)-as.integer(gt2))][]->ab
# ab[sample1 != sample2, table(d)/sum(table(d))]
# gene1.start=419306988; gene1.end=419309556; gene2.start=419364015; gene2.end=419365995
# ab[pos %between% c(gene1.start, gene1.end),gene:='glu.d1x'][gene=='glu.d1x', gene.start:=gene1.start][gene=='glu.d1x',gene.end:=gene1.end]
# ab[pos %between% c(gene2.start, gene2.end),gene:='glu.d1y'][gene=='glu.d1y', gene.start:=gene2.start][gene=='glu.d1y',gene.end:=gene2.end]
# ab[!is.na(gene),aa.pos:=(gene.end-pos)/3]
# ab[sample1 == 'Stanley_5+10' & sample2 == 'TA2580' & d>0]
# 
# dcast(ab[, .N, key=.(sample1, sample2, d)], sample1 + sample2 ~ d, value.var="N", fill=0) -> y
# # y[sample1 == 'Stanley_5+10'][order(`2`)][1:20]
# # y[sample1 == 'Stanley_5+10', boxplot(`2`)]
# 
# setnames(y, c("0","1", "2"), c("ibs2", "ibs1", "ibs0"))[, n := ibs2 + ibs0]
# y[, pdiff := ibs0/n]
# y[sample1 != sample2 & n>100, summary(n)]
# y[sample1 != sample2 & n>100, summary(pdiff)]
# dcast(y, sample1~sample2, value.var='pdiff') -> m
# as.matrix(m[,-c(1)]) -> m2
# list(ibs=y, matrix=m) -> d
# # saveRDS(d, file='~/HMW_glutenin/210525.ibs.matrix.rds')
# 
# library('gplots')
# pdf('/homes/lianggao/HMW_glutenin/data/210525.heatmap.based.on.glu.locus.pdf', width=100, height=100)
# heatmap.2(sqrt(m2), trace='none', keysize=0.5, key.title=NA, key.xlab=NA)
# dev.off()
# 

###############. ####
#### II. new genotyping data based on Liang's genotype recall gvcf, merge with owwc... etc.
### read in annotation
fread('grep -v "#" /homes/lianggao/HMW_glutenin/data/210528.merged.breadwheat.owwc.filtered.ann.vcf | cut -f 1-9') -> vcf.cut
fread('grep -v "#" /homes/lianggao/HMW_glutenin/data/210609.all.filtered.ad.ann.vcf | cut -f 1-9') -> vcf.cut
cbind(vcf.cut[,.(V1,V2,V4,V5)],vcf.cut[,tstrsplit(V8,'\\|')][,2]) -> vcf.cut.split.eff
setnames(vcf.cut.split.eff, c('chr','pos', 'ref','alt', 'eff'))
gene1.start=419306988; gene1.end=419309556; gene2.start=419364015; gene2.end=419365995
vcf.cut.split.eff[pos %between% c(gene1.start, gene1.end),gene:='glu.d1x'][gene=='glu.d1x', gene.start:=gene1.start][gene=='glu.d1x',gene.end:=gene1.end]
vcf.cut.split.eff[pos %between% c(gene2.start, gene2.end),gene:='glu.d1y'][gene=='glu.d1y', gene.start:=gene2.start][gene=='glu.d1y',gene.end:=gene2.end]
c.term.glud1x = c(419306988, 419307120); vcf.cut.split.eff[pos %between% c.term.glud1x,motif:='glu.d1x.cterm']
crr.glud1x =c(419307120, 419309154); vcf.cut.split.eff[pos %between% crr.glud1x, motif:='glud1x.crr']
n.term.glud1x = c(419309154, 419309556); vcf.cut.split.eff[pos %between% n.term.glud1x, motif:='glud1x.nterm']
c.term.glud1y = c(419364015, 419364147); vcf.cut.split.eff[pos %between% c.term.glud1y,motif:='glu.d1y.cterm']
crr.glud1y =c(419364147, 419365551); vcf.cut.split.eff[pos %between% crr.glud1y, motif:='glud1y.crr']
n.term.glud1y = c(419365551, 419365995); vcf.cut.split.eff[pos %between% n.term.glud1y, motif:='glud1y.nterm']
vcf.cut.split.eff[!is.na(gene),aa.pos:=(gene.end-pos)/3]
vcf.cut.split.eff[!is.na(gene),table(eff, motif)]

# write.csv(vcf.cut.split.eff, file='~/HMW_glutenin/archive/210609.vcf.effects.gene.pos.csv', row.names=F)
vcf.cut.split.eff[pos %in% vcf.cut[grepl('Cys', V8), V2],]


# 66666
# # fread('/bulk/lianggao/HMW_glutenin/var_call/210528_v2/210528.merged.breadwheat.owwc.filtered.vcf2long') -> g
# # setnames(g, c("chr", "pos", "ref", "alt", "qual", "sample.id", "gt", "dp", "dv"))
# # g[gt == '0/0', gt := paste(0)]
# # g[gt == '0/1', gt := paste(1)]
# # g[gt == '1/1', gt := paste(2)]
# # g[sample.id %in% c('stanley.merged.Region.2.5k.bam'), .(chr, pos, sample1=sample.id, gt1=gt)] -> a
# # g[sample.id %in% c('TA10103', 'TA10104','TA2582', 'TA2580', 'TA10928'), .(chr, pos, sample2=sample.id, gt2=gt)] -> b
# # a[b, on=c("chr", "pos"), allow.cartesian=T][, d := abs(as.integer(gt1)-as.integer(gt2))][]->ab
# # copy(ab) -> ab.stanley.against.l3
# # dcast(ab.stanley.against.l3[, .N, key=.(sample1, sample2, d)], sample1 + sample2 ~ d, value.var="N", fill=0) -> y.stanley
# # vcf.cut.split.eff[ pos %in% ab.stanley.against.l3[d>1, pos]]
# # 
# # 
# # g[, .(chr, pos, sample1=sample.id, gt1=gt)] -> a
# # g[, .(chr, pos, sample2=sample.id, gt2=gt)] -> b
# # a[b, on=c("chr", "pos"), allow.cartesian=T][, d := abs(as.integer(gt1)-as.integer(gt2))][]->ab
# # ab[sample1 != sample2, table(d)/sum(table(d))]
# # gene1.start=419306988; gene1.end=419309556; gene2.start=419364015; gene2.end=419365995
# # ab[pos %between% c(gene1.start, gene1.end),gene:='glu.d1x'][gene=='glu.d1x', gene.start:=gene1.start][gene=='glu.d1x',gene.end:=gene1.end]
# # ab[pos %between% c(gene2.start, gene2.end),gene:='glu.d1y'][gene=='glu.d1y', gene.start:=gene2.start][gene=='glu.d1y',gene.end:=gene2.end]
# # ab[!is.na(gene),aa.pos:=(gene.end-pos)/3]
# # ab[sample1 == 'stanley.merged.Region.2.5k.bam' & sample2 == 'TA2580' & d>0]
# # 
# # dcast(ab[, .N, key=.(sample1, sample2, d)], sample1 + sample2 ~ d, value.var="N", fill=0) -> y
# # setnames(y, c("0","1", "2"), c("ibs2", "ibs1", "ibs0"))[, n := ibs2 + ibs0]
# # y[, pdiff := ibs0/n]
# # y[sample1 != sample2 & n>100, summary(n)]
# # y[sample1 != sample2 & n>100, summary(pdiff)]
# # y[sample1=='stanley.merged.Region.2.5k.bam', ][order(pdiff)]
# # list(vcf.eff=vcf.cut.split.eff, vcf.long.ab=ab, sample2sample.ident=y) -> results.210528
# # saveRDS(results.210528, file='~/HMW_glutenin/data/210528.results.vcf.eff_long.ab_sample2sample.ident.y.RDS')

dcast(y, sample1~sample2, value.var='pdiff') -> m
as.matrix(m[,-c(1)]) -> m2
library('gplots')
pdf('/homes/lianggao/HMW_glutenin/data/210528.heatmap.based.on.glu.locus.pdf', width=100, height=100)
heatmap.2(sqrt(m2), trace='none', keysize=0.5, key.title=NA, key.xlab=NA)
dev.off()

### checking vcf2long call differences between gvcf merged version and -mv merged version
# fread('/bulk/lianggao/HMW_glutenin/var_call/210528.insilico.jagger.aet.10genomes.glud1/210528.merged.breadwheat.owwc.filtered.vcf2long') -> g
# setnames(g, c("chr", "pos", "ref", "alt", "qual", "sample.id", "gt", "dp", "dv"))
# g[,N:=.N, by=sample.id]
# g[!duplicated(sample.id), .(sample.id, N)] 
# g[!duplicated(sample.id), hist(N, main='normal vcf', col='#41ab5d', xlim=c(500,2000))]
# fread('/bulk/lianggao/HMW_glutenin/var_call/210528_v2/210528.merged.breadwheat.owwc.filtered.vcf2long') -> g2
fread('/homes/lianggao/HMW_glutenin/data/210609.all.filtered.ad.ann.vcf2long') -> g2
setnames(g2, c("chr", "pos", "ref", "alt", "qual", "sample.id", "gt", "dp", "dv"))
# g2[,N:=.N, by=sample.id]
# g2[!duplicated(sample.id), .(sample.id, N)]
# g2[!duplicated(sample.id), hist(N, main='gvcf', col='#4292c6', xlim=c(500,2000))]
# par(mfrow=c(2,1))
# g2[!duplicated(sample.id), .(sample.id, N)]
# g[g2, on=c('chr','pos','ref','sample.id')] -> gg
# gg[sample.id=='cs.merged.Region.2.5k.bam' & alt == i.alt]


vcf.cut.split.eff[,.(chr=1, start=pos-1, end=pos, eff)] -> vcf.cut.split.eff.v2
data.table(molecule='genome1',gene='Glu-D1x', start=419306988, end=419309556, strand='reverse', orientation=-1) -> test.gene
data.table(molecule='genome1',gene='Glu-D1y', start=419364015, end=419365995, strand='reverse', orientation=-1) -> test.gene2
vcf.cut.split.eff.v2[start %between% c(test.gene$start, test.gene$end)] -> vcse1 ## glud1x eff
vcf.cut.split.eff.v2[start %between% c(test.gene2$start, test.gene2$end)] -> vcse2 ## glud1y eff
data.table(test.gene, from=vcse1$start, to=vcse1$end, eff=vcse1$eff) -> test.gene.sub; test.gene.sub[,subgene:=paste0('segX',1:.N)]
data.table(test.gene2, from=vcse2$start, to=vcse2$end, eff=vcse2$eff) -> test.gene.sub2; test.gene.sub2[,subgene:=paste0('segY',1:.N)]
rbind(test.gene, test.gene2) -> glu.d1.gene
glu.d1.gene[grepl('D1y', gene), molecule:='genome2']
rbind(test.gene.sub, test.gene.sub2) -> glu.d1.snps
glu.d1.snps[grepl('D1y', gene), molecule:='genome2']g2[pos==419309204, table(gt)]

# saveRDS(list(glu.d1.gene,glu.d1.snps), file='~/HMW_glutenin/data/210601_glu.d1.gene.snps.RDS')

ggplot(glu.d1.gene, aes(xmin = start, xmax = end, y = molecule, fill = gene, forward = orientation)) +
    geom_gene_arrow(arrow_body_height = grid::unit(10,"mm"),
                    arrowhead_height = grid::unit(15,'mm'),
                    arrowhead_width = grid::unit(15,'mm')) +
    facet_wrap(~ molecule, scales = "free", ncol = 2)+
    scale_fill_brewer(palette = "Set3")+
    geom_subgene_arrow(data = glu.d1.snps,
                       arrow_body_height = grid::unit(10,"mm"),
                       arrowhead_height = grid::unit(15,'mm'),
                       arrowhead_width = grid::unit(15,'mm'),
                       aes(xmin = start, xmax = end, y = molecule, fill = gene,
                           xsubmin = from, xsubmax = to), color='black') +
    geom_subgene_arrow(data = glu.d1.snps[grepl('missen', eff)],
                       arrow_body_height = grid::unit(10,"mm"),
                       arrowhead_height = grid::unit(15,'mm'),
                       arrowhead_width = grid::unit(15,'mm'),
                       aes(xmin = start, xmax = end, y = molecule, fill = gene,
                           xsubmin = from, xsubmax = to), color='#cc4c02') +
    theme_genes() + theme(axis.text=element_text(size=16, face = 'bold', angle=90), axis.text.y=element_blank())
dev.off()
g2[,sample.id:=sub('.bam','', sample.id)]
ids[taxa %in% g2[,sample.id]]
View(g2[pos=='419311695' & gt=='1/1' ])
(g2[pos==419311695 & gt=='1/1' & grepl('TA', sample.id) ])
ids[g2[pos==419311695], on=c(taxa='sample.id')][,table(lineage, gt)] 

vcf.cut.split.eff[eff=='stop_gained']
vcf.cut.split.eff[eff=='stop_gained', pos]
g2[pos %in% vcf.cut.split.eff[eff=='stop_gained', pos] & gt=='1/1']

vcf.cut.split.eff[!is.na(gene),pos] -> pos
dist.vec=c()
for (i in 1:length(pos)){
    dist.vec=c(dist.vec,pos[i+1]-pos[i])
}
table(dist.vec)

fread('~/HMW_glutenin/data/find.cys.out.txt') ->find.cys
unique(find.cys[,2:11])
