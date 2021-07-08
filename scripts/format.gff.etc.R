
library(data.table)
chr.format=fread('/bulk/jpoland/genome/Ae_tauschii/AL8_78/assembly/GCA_002575655.1_Aet_v4.0.chr.num.agp.format.txt')
chr.format[,.(chr,chr_gff=paste0('Chr',1:7))] -> c
fread('/bulk/lianggao/genome/Aet_v4.0/annotation_Aet_v4/AET_High_confidence_gene.gff3') -> gff
gff[c, on=c(V1='chr_gff')] -> gffc

gffc[,V1:=chr]
gffc[,1:9] -> gff.new
write.table(gff.new, file='/bulk/lianggao/genome/Aet_v4.0/annotation_Aet_v4/AET_High_confidence_gene.new.chrid.gff3', sep='\t', quote=F, row.names=F, col.names=F)





