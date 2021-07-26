library(data.table)
library(openxlsx)
library(ggplot2)
library(scales)


data.table(read.xlsx('/homes/lianggao/HMW_glutenin/data/10Genomes_HMW-GS-Loci_positions_Mar24_ED2.xlsx')) -> d
d[, unique(Genome)][c(1,3,5,6,7,8,9,12,13,14,16)] -> genome10
d[Genome %in% genome10 & allele %in% c('2+12', '5+10', '2.2+12'),] -> d.sub10
d.sub10[allele=='2.2+12', allele:='2+12']
d.sub10[,]


## read in kmer  assignments matrix
admix.tabs=list()
for (i in 1:11) {
        admix.tabs[[i]] = data.table(read.xlsx('/homes/lianggao/HMW_glutenin/data/lineage_admixture_raw_counts.xlsx', sheet=i))
        colnames(admix.tabs[[i]]) = c('chr','start', 'end', 'l1.kmer','l2.kmer','l3.kmer','all')
        admix.tabs[[i]][, l3.pct:=l3.kmer/all][,l2.pct:=l2.kmer/all][,l1.pct:=l1.kmer/all]
}

## assign lineage based on simple stats
genome.names.11 = c('ArinaLrFor', 'Chinese_Spring', 'Jagger', 'Julius', 'LongReach_Lancer', 'CDC_Landmark', 'Mace', 'SY_Mattis', 'Norin_61', 'PI190962_Spelt','CDC_Stanley')
for (i in 1:11){
        admix.tabs[[i]][l3.kmer-l2.kmer>10 & l3.kmer-l1.kmer > 10 & l3.pct>0.2, lineage.assign:='L3'][,genome:=genome.names.11[i]]
        admix.tabs[[i]][l2.kmer-l1.kmer>10 & l2.kmer-l3.kmer > 10 & l2.pct>0.2, lineage.assign:='L2'][,genome:=genome.names.11[i]]
        admix.tabs[[i]][l1.kmer-l2.kmer>10 & l1.kmer-l3.kmer > 10 & l1.pct>0.2, lineage.assign:='L1'][,genome:=genome.names.11[i]]
}
rbindlist(admix.tabs) -> admix.tabs
admix.tabs[is.na(lineage.assign), lineage.assign:='undecided']
admix.tabs[,table(lineage.assign)]

chr1d.400 = data.table(molecule=genome.names.11, gene='genome', start=350e6, end=450e6, strand='forward', orientation=1)
admix.tabs[chr=='chr1D', max(end), by=genome][,.(genome, chr=1:11,start=0, length=V1)] -> chr1d.all
admix.tabs[chr=='chr1D', ] -> admix.tabs.1d
admix.tabs.1d[chr1d.all, on ='genome'] -> admix.tabs.1d.formats



genomes.key=data.table(x=unique(d.sub10$Genome), y=genome.names.11[c(1,6,3,4,5,7,9,10,11,8,2)])
d.sub10[genomes.key, on = c(Genome='x')] -> d.sub10.formats
chr1d.all[d.sub10.formats,on=c(genome='y')] -> d.sub10.formats.more

pdf('~/HMW_glutenin/data/publications/210609.kmer.regions.10genomes.pdf', width = 10, height=6)
ggplot()+geom_bar(data=chr1d.all, aes(chr,length), width=0.5, stat='identity', fill='grey70')+
        geom_segment(data=admix.tabs.1d.formats[lineage.assign=='L2'], aes(x=as.integer(i.chr)-.25, xend=as.integer(i.chr)+.25, y=start, yend=start), color='#0073B3')+
        geom_segment(data=admix.tabs.1d.formats[lineage.assign=='L3'], aes(x=as.integer(i.chr)-.25, xend=as.integer(i.chr)+.25, y=start, yend=start), color='#D45F14')+
        geom_segment(data=admix.tabs.1d.formats[lineage.assign=='L1'], aes(x=as.integer(i.chr)-.25, xend=as.integer(i.chr)+.25, y=start, yend=start), color='#F0E442')+
        geom_point(data=d.sub10.formats.more[allele=='2+12'], aes(x=as.integer(chr)+.5,y=end), fill='#0073B3', shape=25, size=3)+
        geom_point(data=d.sub10.formats.more[allele=='5+10'], aes(x=as.integer(chr)+.5,y=end), fill='#D45F14', shape=25, size=3)+
        # scale_y_continuous(name='Chr1D physical position (Mb)', labels=c(0,100,200,300,400,500))+
        geom_hline(yintercept=c(390e6,430e6), size=1)+
        scale_y_continuous(name='', labels=c(0,100,200,300,400,500), oob=rescale_none)+
        ggtitle('Lineage specific kmers distributions on 10+ genomes') +
        theme(axis.title=element_text(face='bold', size=14), axis.text=element_text(face='bold', size=12), plot.title=element_text(face='bold',size=16, hjust=0.5), plot.margin=unit(c(2,2,2,2), 'cm'))+
        scale_x_discrete(name='', limits=genome.names.11) +coord_flip() 
dev.off()

pdf('~/HMW_glutenin/data/publications/210609.kmer.regions.10genomes.zoom.in.pdf', width = 10, height=6)
ggplot()+geom_bar(data=chr1d.all, aes(chr,length), width=0.5, stat='identity', fill='grey70')+
        geom_segment(data=admix.tabs.1d.formats.zoom[lineage.assign=='L2'], aes(x=as.integer(i.chr)-.25, xend=as.integer(i.chr)+.25, y=start, yend=end), color='#0073B3')+
        geom_segment(data=admix.tabs.1d.formats[lineage.assign=='L3'], aes(x=as.integer(i.chr)-.25, xend=as.integer(i.chr)+.25, y=start, yend=start), color='#D45F14')+
        geom_segment(data=admix.tabs.1d.formats[lineage.assign=='L1'], aes(x=as.integer(i.chr)-.25, xend=as.integer(i.chr)+.25, y=start, yend=start), color='#F0E442')+
        geom_point(data=d.sub10.formats.more[allele=='2+12'], aes(x=as.integer(chr)+.5,y=end), fill='#0073B3', shape=25, size=3)+
        geom_point(data=d.sub10.formats.more[allele=='5+10'], aes(x=as.integer(chr)+.5,y=end), fill='#D45F14', shape=25, size=3)+
        geom_segment(data=d.sub10.formats.more[allele=='2+12'], aes(x=as.integer(chr)-.25, xend=as.integer(chr)+.25, y=i.start, yend=end), color='black', shape=25, size=1)+
        geom_segment(data=d.sub10.formats.more[allele=='5+10'],  aes(x=as.integer(chr)-.25, xend=as.integer(chr)+.25, y=i.start, yend=end), color='black', shape=25, size=1)+
        scale_y_continuous(name='Chr1D physical position (Mb)', limits=c(390e6,430e6), labels=c(390,400,410,420, 430), oob = rescale_none)+
        # ggtitle('Lineage specific kmers distributions on 10+ genomes') +
        theme(axis.title=element_text(face='bold', size=14), axis.text=element_text(face='bold', size=12), plot.title=element_text(face='bold',size=16, hjust=0.5), plot.margin=unit(c(2,2,2,2), 'cm'))+
        scale_x_discrete(name='', limits=genome.names.11) +coord_flip() 
dev.off()     
         


# ggplot()+geom_bar(data=chr1d.all, aes(chr,length), width=0.5, stat='identity', fill='grey')+
#         geom_segment(data=admix.tabs.1d.formats[lineage.assign=='L2'], aes(x=as.integer(i.chr)-.25, xend=as.integer(i.chr)+.25, y=start, yend=end), color='#0073B3')+
#         geom_segment(data=admix.tabs.1d.formats[lineage.assign=='L3'], aes(x=as.integer(i.chr)-.25, xend=as.integer(i.chr)+.25, y=start, yend=end), color='#D45F14')+
#         geom_segment(data=admix.tabs.1d.formats[lineage.assign=='L1'], aes(x=as.integer(i.chr)-.25, xend=as.integer(i.chr)+.25, y=start, yend=end), color='#F0E442')+
#         geom_point(data=d.sub10.formats.more[allele=='2+12'], aes(x=as.integer(chr)+.5,y=end), fill='#0073B3', shape=25, size=3)+
#         geom_point(data=d.sub10.formats.more[allele=='5+10'], aes(x=as.integer(chr)+.5,y=end), fill='#D45F14', shape=25, size=3)+
#         geom_segment(data=d.sub10.formats.more[allele=='2+12'], aes(x=as.integer(chr)-.25, xend=as.integer(chr)+.25, y=i.start, yend=end), color='black',  size=3)+
#         geom_segment(data=d.sub10.formats.more[allele=='5+10'],  aes(x=as.integer(chr)-.25, xend=as.integer(chr)+.25, y=i.start, yend=end), color='black', size=3)+
#         scale_y_continuous(name='Chr1D physical position (Mb)', limits=c(400e6,420e6), labels=c(400, 405,410,415,420))+
#         ggtitle('Lineage specific kmers distributions on 10+ genomes') +
#         theme(axis.title=element_text(face='bold', size=14), axis.text=element_text(face='bold', size=12), plot.title=element_text(face='bold',size=16, hjust=0.5), plot.margin=unit(c(2,2,2,2), 'cm'))+
#         scale_x_discrete(name='cultivar', limits=genome.names.11) +coord_flip() 
# 
