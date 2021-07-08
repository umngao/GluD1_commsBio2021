library(ggtree)
library(treeio)
library(data.table)

fread('~/HMW_glutenin/data/phylogenetic_trees/201215.owwc.306.lineages.based.on.PC.txt') -> ids.lineage
# ids.lineage[,plot(PC1, PC2)]
# ids.lineage[!grepl('L', lineage),points(PC1, PC2, col='red')]
ids.lineage[,.(taxa,lineage)][data.table(labels=tree$tip.label), on=c(taxa='labels')] -> ids
ids[is.na(lineage), lineage:='wheat']
ids[!grepl('L', lineage)]
ids[lineage=='wheat' & grepl('TA', taxa), lineage:='hybrid']
ids[, new.id:=taxa]
ids[grepl('Stanley', taxa, ignore.case=T), new.id:='Stanley5+10']
ids[grepl('cs', taxa), new.id:='CS2+12']
ids[grepl('BW971', taxa, ignore.case=T), new.id:='Landmark']
ids[grepl('landmark', taxa, ignore.case=T),new.id:='Landmark2+12']
ids[grepl('lancer', taxa, ignore.case=T),new.id:='Lancer2+12']
ids[grepl('CS.PE470', taxa), new.id:='cs.PE0601']
ids[, idx:=1:.N, by=new.id]
ids[!grepl('L', lineage), new.id:=paste(new.id,idx, sep='_')]
# write.table(ids, file='~/HMW_glutenin/data/phylogenetic_trees/210531.ids.w.lineages.txt', sep='\t')

## read in trees
tree = read.tree('~/HMW_glutenin/data/phylogenetic_trees/glud1x/glud1x.all.aa/210530.recall.glud1x.M.HR.aa.fa.treefile_copy.txt')
groupInfo=list(L1=ids[lineage=='L1', taxa], 
               L2=ids[lineage=='L2', taxa],
               L3=ids[lineage=='L3', taxa],
               wheat=ids[lineage=='wheat', taxa],
               hybrid=ids[lineage=='hybrid', taxa]
               )
tree=groupOTU(tree, groupInfo)
tree$tip.label = ids$new.id
library(ggrepel)
ggtree(tree, layout='equal_angle')+geom_tiplab(aes(color=group))
as.treedata(tree) -> tt
as.phylo(tt) -> tt
write.tree(tt, file='~/HMW_glutenin/data/phylogenetic_trees/210601.tree.nwk')
rbind(ids[grepl('L1', lineage),.(new.id, 'label', '#F0E442', 'bold', 3)],
ids[grepl('L2', lineage),.(new.id, 'label', '#0073B3', 'bold', 3)],
ids[grepl('L3', lineage),.(new.id, 'label', '#D45F14', 'bold', 3)],
ids[lineage=='wheat', .(new.id, 'label', '#009E73', 'bold', 3)],
ids[lineage=='hybrid', .(new.id, 'label', 'grey', 'bold', 3)])  -> dd
# write.table(dd, file='~/HMW_glutenin/data/phylogenetic_trees/label.color2.txt', sep='\t', quote=F, row.names=F, col.names=F)



tree = read.tree('~/HMW_glutenin/data/phylogenetic_trees/glud1y/glud1y.all.aa/210530.recall.glud1y.M.HR.true.aa.fa.treefile')
data.table(labels=tree$tip.label)[,labels.y:=sub('.glud1y_5', '', labels)][,labels.y:=sub('.bam', '', labels.y)][] -> labs.y
ids[labs.y, on=c(taxa='labels.y')] -> ids2
groupInfo=list(L1=ids2[lineage=='L1', labels], 
               L2=ids2[lineage=='L2', labels],
               L3=ids2[lineage=='L3', labels],
               wheat=ids2[lineage=='wheat', labels],
               hybrid=ids2[lineage=='hybrid', labels]
)
tree=groupOTU(tree, groupInfo)
tree$tip.label = ids2$new.id
ggtree(tree, layout='equal_angle')+geom_tiplab(aes(color=group))
as.treedata(tree) -> tt
as.phylo(tt) -> tt
# write.tree(tt, file='~/HMW_glutenin/data/phylogenetic_trees/210601.tree.glud1y.aa.nwk')
# write.table(dd, file='~/HMW_glutenin/data/phylogenetic_trees/label.color2.txt', sep='\t', quote=F, row.names=F, col.names=F)



