library(reshape2)
library(ggplot2)
#read in te pav data which has been filtered for population structure as well as DMR data:
tes.input=read.delim('te.pol.fixed2.txt',head=T)
met.input=read.delim('mc_dmrs.tsv',head=T)
#subset methylation accessions to remove duplicates and fix names for bud tissue samples
met.filt=met.input[,intersect(grep('Ag_0_bud',names(met.input),invert=T),intersect(grep('Bor_4_bud',names(met.input),invert=T),intersect(grep('Ema_1_bud',names(met.input),invert=T),intersect(grep('Er_0_bud',names(met.input),invert=T),intersect(grep('Fr_2_bud',names(met.input),invert=T),intersect(grep('Is_0_bud',names(met.input),invert=T),intersect(grep('Kin_0_bud',names(met.input),invert=T),intersect(grep('Litva_bud',names(met.input),invert=T),intersect(grep('Pu2_23_bud',names(met.input),invert=T),intersect(grep('Ragl_1_bud',names(met.input),invert=T),intersect(grep('Uk_1_bud',names(met.input),invert=T),intersect(grep('Zdr_1_bud',names(met.input),invert=T),intersect(grep('vi_0',names(met.input),invert=T),intersect(grep('HR_5',names(met.input),invert=T),grep('Ler_1',names(met.input),invert=T)))))))))))))))]
names(met.filt)=gsub("_bud","",names(met.filt))
#perform a first pass to identify samples with matching names between the input TE and DMR files
havename=intersect(names(tes.input),names(met.filt))
#subset the TE data to those with methylation data
tes.with.met.data=tes.input[,colnames(tes.input)%in%havename]
met.with.te.data=met.filt[,colnames(met.filt)%in%havename]
#order the accessions the same way in both TE and DMR sets
tes.with.met.data=tes.with.met.data[,order(names(tes.with.met.data))]
met.with.te.data=met.with.te.data[,order(names(met.with.te.data))]
#remove the conserved 'stop' column!
tes.with.met.data=tes.with.met.data[,c(1:102,104:122)]
met.with.te.data=met.with.te.data[,c(1:102,104:122)]
#There are three accessions that are present in both sets, but have typos in their names. We will add these three samples manualy to the end of the file (same order in both)
tes.with.met.data=cbind.data.frame(tes.input[,1:7],tes.with.met.data,tes.input$Da1_12,tes.input$Li_2.1,tes.input$Rubeznhoe_1)
met.with.te.data=cbind.data.frame(met.filt[,1:3],met.with.te.data,met.filt$Da_1_12,met.filt$Li_2_1,met.filt$Rubezhnoe_1)
names(tes.with.met.data)[(ncol(tes.with.met.data)-2):ncol(tes.with.met.data)]=c('Da_1_12','Li_2_1','Rubeznhoe_1')
names(met.with.te.data)[(ncol(met.with.te.data)-2):ncol(met.with.te.data)]=c('Da_1_12','Li_2_1','Rubeznhoe_1')
#We can now confirm that columns have the same names. We omit the first three columns of dmr data (chr,start,stop) and the 7 metadata columns in the TE dataset
namecheck=cbind(names(met.with.te.data[,4:ncol(met.with.te.data)]),names(tes.with.met.data[,8:ncol(tes.with.met.data)]))
namecheck[,2]=gsub("\\.","_",namecheck[,2])
table(ifelse(namecheck[,1]==namecheck[,2],T,F))
#There should be no non-uniform columns:
namecheck[ifelse(namecheck[,1]==namecheck[,2],F,T),]
#rename dmr set to make things easier on the eyes
dmrs=met.with.te.data
#there are a couple of DMRs which have no variation across samples. Can remove them with this
dmrs=subset(dmrs,rowSums(dmrs[,4:ncol(dmrs)])!=0)
#same for tes:
tes=tes.with.met.data
#determine the number of present TE calls across all samples. Do not include TE-sample values that are not present (NAs)
tes$count_present=rowSums(tes[,8:ncol(tes)],na.rm=T)
#determine the number of samples that have any call at all
tes$count_notna=ncol(tes[,8:ncol(tes)]) - apply(tes[,8:(ncol(tes)-1)],MARGIN=1,FUN=function(x) length(x[is.na(x)]))
############
#set a hard limit to only include TE PAV in which we end up with at least 100 samples with data
tes=subset(tes,tes$count_notna > 100)
#now make a column filtering TE PAV requiring the MAF to be > 3% from the set of samples that have data (1: PASS / 0:FAIL)
tes$filter=ifelse((tes$count_present >= tes$count_notna*0.03) & (tes$count_present <= (tes$count_notna*0.97)),1,0)
#filter down to tes passing our MAF cutoff
tes=subset(tes,tes$filter==1)
#how many tes are we going to work with
dim(tes)
#File from Tim providing coordinates for 5 centromeres in A. thaliana
centromere=read.delim('centromeres.txt',head=F)
#grab the chr,start,stop information for the TEs and DMRs:
te.info=tes[,1:3]
dmr.info=dmrs[,1:3]
#create matrix to fill with plot coordinates for geom_rect() on our plot. Values are NOT actual chr,start,stop for plotting.
te.call=matrix(NA,nrow(te.info),ncol=5)
dmr.call=matrix(NA,nrow(dmr.info),ncol=5)
#for each element, determine if it falls within the percentromeric region (+/- 300kb around centromere boundaries)
for(i in 1:nrow(te.info)){for(q in 1:5){te.call[i,q]=ifelse(te.info[i,1]==centromere[q,1] & te.info[i,2]>=centromere[q,2]-300000 & te.info[i,3]<=centromere[q,3]+300000,1,0)}}
for(i in 1:nrow(dmr.info)){for(q in 1:5){dmr.call[i,q]=ifelse(dmr.info[i,1]==centromere[q,1] & dmr.info[i,2]>=centromere[q,2]-300000 & dmr.info[i,3]<=centromere[q,3]+300000,1,0)}}
#collapse into a single vector across the set
te.call=rowSums(te.call)
dmr.call=rowSums(dmr.call)
#identify the changepoints going in and out of perecentromeric regions
te.perecentromeric.breakpoints=c(1,1+which(diff(te.call)!=0))
dmr.perecentromeric.breakpoints=c(1,1+which(diff(dmr.call)!=0))
#create a dataframe in use for ggplot geom_rect() that will provide plot positions for each chromosomes perecentromeric boundaries
te.x=cbind.data.frame(rep(1:5,each=2),te.perecentromeric.breakpoints[-1],rep(c('xmin','xmax')))
te.y=cbind.data.frame(rep(1:5,each=2),dmr.perecentromeric.breakpoints[-1],rep(c('ymin','ymax')))
colnames(te.x)=c('chrom','value','class')
colnames(te.y)=c('chrom','value','class')
centromere.reg=rbind(te.x,te.y)
centromere.reg=dcast(centromere.reg,chrom ~ class,value.var='value')
#looks like this:
centromere.reg
#Also for plotting, identify where chromosome breaks are within the file for visualisation:
te.chr.breaks=as.vector(cumsum(table(te.info[,1])))[1:4]
dmr.chr.breaks=as.vector(cumsum(table(dmr.info[,1])))[1:4]
#push out dmrs and tes for finding closest DMR (cis band analysis)
write.table(dmrs,'dmrs.bed',sep='\t',row.names=F,col.names=F,quote=F)
write.table(tes,'tes.bed',sep='\t',row.names=F,col.names=F,quote=F)
#### RUN OUTSIDE OF R
#bedtools closest -D "ref" -a tes.bed -b dmrs.bed -t first > local_data.out
####
#read in cis data
local.out=read.delim('local_data.out',head=F)
#setup output matrix
out.local.r2=matrix(NA,ncol=2,nrow=nrow(local.out))
#test real data
for(i in 1:nrow(out.local.r2)){out.local.r2[i,1]=cor(t(local.out[i,8:131]),t(local.out[i,138:261]),use='pairwise.complete.obs')}
#test permuted data
set.seed(42)
for(i in 1:nrow(out.local.r2)){out.local.r2[i,2]=cor(t(local.out[i,8:131]),t(local.out[i,sample(138:261)]),use='pairwise.complete.obs')}
#make r2 values
out.local.r2=(out.local.r2)^2
#format for plotting
colnames(out.local.r2)=c('real','permuted')
out.melt=melt(out.local.r2)
########## plot with 1% permuted line highlighted
pdf('FIGURE6B.pdf')
ggplot(out.melt,aes(value)) + geom_density(aes(fill=Var2),alpha=0.5) + coord_cartesian(ylim=c(0,10),xlim=c(0,0.5)) + scale_fill_manual(values=c('#ff9933','#999999')) +geom_vline(xintercept=quantile(out.local.r2[,2],0.99),col='red',alpha=0.8) + theme_classic() + xlab('r2 value')
dev.off()
##########
#How enriched is real cis r2 compared to permuted?
sum(out.local.r2[,1] >= quantile(out.local.r2[,2],0.99)) / sum(out.local.r2[,2]>=quantile(out.local.r2[,2],0.99))
#########
#back to trans data
#########
#create a matrix of only the actual data
te.matrix=tes[,8:(ncol(tes)-3)]
dmr.matrix=dmrs[,4:ncol(dmrs)]
#the correlation is as simple as transposing our TE PAV and DMR datasets. We are using 'pairwise.complete.obs' to indicate a direct pairwise comparison for complete observations only
test=cor(t(te.matrix),t(dmr.matrix),use='pairwise.complete.obs')
#remove names at this point just to prevent factor changes in ggplot
colnames(test)=NULL
rownames(test)=NULL
#develop r2 values from the correlation values
r2.values=test^2
#we also want a permuted dataset. We will create this by randomizing our DMR samples. Seed is set to always return same permuted set.
set.seed(42)
dmrs.permute=dmr.matrix[,sample(ncol(dmr.matrix))]
#perform the same correlation and steps as above
test.permute=cor(t(te.matrix),t(dmrs.permute),use='pairwise.complete.obs')
colnames(test.permute)=NULL
rownames(test.permute)=NULL
r2.values.permute=test.permute^2
#as we are plotting in ggplot2, it wants long-form dataframes. Lets make some
r2.melt=melt(r2.values)
r2.melt.permuted=melt(r2.values.permute)
#add and identifyer column 'class' for both
r2.melt$class=rep('real')
r2.melt.permuted$class=rep('permuted')
#stick them together and make one massive long-form dataframe for ggplotting
r2.melt.all=rbind(r2.melt,r2.melt.permuted)
#determine the top 1% quantile cutoff for the permuted data. This will be our threshold to see how many r2 values are above the top 1% of r2 values in the permuted dataset
r2.quantile.cutoff=quantile(r2.values.permute,0.99,na.rm=T)
#so what r2 value is deemed good?
r2.quantile.cutoff
##########
pdf('FIGURE6C.pdf')
ggplot(r2.melt.all,aes(value)) + geom_density(aes(fill=class),alpha=0.5) + coord_cartesian(ylim=c(0,1),xlim=c(0,0.5)) + scale_fill_manual(values=c('#999999','#ff9933')) + geom_vline(xintercept=r2.quantile.cutoff,col='red',alpha=0.8) + theme_classic() + xlab('r2 value')
dev.off()
##########
#using our 0.1% threshold of permuted data as a cutoff, we can make our matrix a binary state
r2.values.binary=r2.values
r2.values.binary[r2.values.binary < r2.quantile.cutoff] <- 0
r2.values.binary[r2.values.binary >= r2.quantile.cutoff] <- 1
r2.values.permute.binary=r2.values.permute
r2.values.permute.binary[r2.values.permute.binary < r2.quantile.cutoff] <- 0
r2.values.permute.binary[r2.values.permute.binary >= r2.quantile.cutoff] <- 1
binary.melt=melt(r2.values.binary)
binary.permute.melt=melt(r2.values.permute.binary)
#add and identifyer column 'class' for both
binary.melt$class=rep('real')
binary.permute.melt$class=rep('permuted')
#stick them together and make one massive long-form dataframe for ggplotting
binary.melt.all=rbind(binary.melt,binary.permute.melt)
#determine the top 1% quantile cutoff for the permuted data. This will be our threshold to see how many COUNT values are above the top 1% of r2 values in the permuted dataset
binary.quantile.cutoff=quantile(r2.values.permute.binary,0.99,na.rm=T)
#what is our binary cutoff value?
binary.quantile.cutoff
#For each TE PAV, perform a binary sum across the DMRs (i.e. the number of DMRs that pass the threshold)
real.te.count.binary=rowSums(r2.values.binary,na.rm=T)
perm.te.count.binary=rowSums(r2.values.permute.binary,na.rm=T)
te.rowcount.binary=cbind(real.te.count.binary,perm.te.count.binary)
sum.binary.quantile.cutoff=quantile(perm.te.count.binary,0.99,na.rm=T)
##########
pdf('FIGURE6D.pdf')
ggplot(melt(te.rowcount.binary),aes(value)) + geom_density(aes(fill=Var2),alpha=0.5) +  theme_classic() + xlab('Sum of DMRs above 1% threshold') + ylab('Density') + scale_fill_manual(values=c('#ff9933','#999999')) + geom_vline(xintercept=sum.binary.quantile.cutoff,col='red',alpha=.5)
dev.off()
##########
#heatmap plot of all pairwise comparisons
##########
png('FIGURE6A.png',width=9000,height=9000)
ggplot(data = binary.melt, aes(x=Var1, y=Var2, fill=value)) + geom_tile() + scale_fill_gradient2(mid='white',high='black') + geom_vline(xintercept=te.chr.breaks,col='black') + geom_hline(yintercept=dmr.chr.breaks,col='black') + geom_rect(data=centromere.reg,inherit.aes=F,aes(xmax=xmax,ymax=ymax,xmin=xmin,ymin=ymin),fill='cyan',alpha=0.1) + theme_classic() + xlab('Transposable Element Polymorphisms') + ylab('Differentially Methylated Regions')
dev.off()
##########
#will also create same plot with permuted data to show how cis-band goes away
##########
png('FIGURE6_SUPPLEMENT_1.png',width=9000,height=9000)
ggplot(data = binary.permute.melt, aes(x=Var1, y=Var2, fill=value)) + geom_tile() + scale_fill_gradient2(mid='white',high='black') + geom_vline(xintercept=te.chr.breaks,col='black') + geom_hline(yintercept=dmr.chr.breaks,col='black') + geom_rect(data=centromere.reg,inherit.aes=F,aes(xmax=xmax,ymax=ymax,xmin=xmin,ymin=ymin),fill='cyan',alpha=0.1) + theme_classic() + xlab('Transposable Element Polymorphisms') + ylab('Differentially Methylated Regions')
dev.off()
##########
#For each TE PAV, perform a r2 sum across the DMRs
real.te.r2.count=rowSums(r2.values,na.rm=T)
perm.te.r2.count=rowSums(r2.values.permute,na.rm=T)
te.rowcount=cbind(real.te.r2.count,perm.te.r2.count)
#determine top 1% cutoff value of permuted data r2 sums
sum.r2.quantile.cutoff=quantile(perm.te.r2.count,0.99,na.rm=T)

#create summary dataframe containg our trans-band sum r2 and binary counts for real and permuted data
te.sums=cbind.data.frame(
  seq(1:length(rowSums(r2.values))),
  real.te.r2.count,
  perm.te.r2.count,
  real.te.count.binary,
  perm.te.count.binary)
#tack on te metadata
te.sums=cbind.data.frame(tes[,1:7],te.sums)
#name our columns
colnames(te.sums)=c('chr','start','stop','teid','te_names_at_site','num_pres','num_abs','pos','real_r2sum','permuted_r2sum','real_binarysum','permuted_binarysum')
#add column indicating a punative trans-effect TE PAV if binary sum is greater than permuted binary sum threshold
##########
perm.te.count.binary.threshold=quantile(perm.te.count.binary,0.99)
##########
te.sums$punative_trans_te=ifelse(te.sums$real_binarysum >= perm.te.count.binary.threshold, 1, 0)
#quick check how many pass our threshold:
table(te.sums$punative_trans_te)
#grab our LD calls for TE PAV from prior analysis
age.call=read.delim('dec3_outlist_snpassociation.txt',head=T)
#merge with our summary datafraome
te.sums.merge=merge(te.sums,age.call,by.x=c('chr','start','teid'),by.y=c('chrom','pos','TEID'),all.x=T)
#as not all TE PAV in this analysis were tested in the LD analysis, we will note this in the LD flag column as 'not_tested'
te.sums.merge$flag=ifelse(is.na(te.sums.merge$flag)==T,'not_tested',as.character(te.sums.merge$flag))
#make a varient version of the flag column in which we only provide LD information if the trans-band is above our permute-data-based threshold (perm.te.count.binary.threshold)
te.sums.merge$flag_pass=ifelse(te.sums.merge$punative_trans_te==0,'under_threshold',as.character(te.sums.merge$flag))

##########
# Figure X-X: Barplots of TE PAV trans-band sums
##########
pdf('FIGURE6E.pdf',width=20,height=5)
ggplot(te.sums.merge,aes(pos,real_binarysum)) + geom_bar(stat='identity',aes(fill=flag_pass),width=1) + theme_classic() + expand_limits(y=c(0,1200)) + scale_y_reverse() + ggtitle('count above threshold real data') + geom_hline(yintercept=perm.te.count.binary.threshold,color='red',alpha=0.5) + scale_fill_manual(values=c('#00BA38','darkgrey','#619CFF','black','#F8766D'))
dev.off()
##########


##########
# FILE OUTPUT: Summary of Trans-band analysis and punative sites
##########
te.sums.merge=te.sums.merge[,c(1,2,4,3,5:ncol(te.sums.merge))]
write.table(te.sums.merge,'te_trans.txt',sep='\t',row.names=F,quote=F)
##########

#outside_of_R
################
#sed 's/^c/C/g' te_trans.txt > te_trans.bed
#tail -n +2 te_trans.bed > te_trans.nohead.bed
#bedtools closest -D "ref" -a te_trans.nohead.bed -b TAIR10_gene.sorted.bed > te_trans_nearest_gene.bed
#wget https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_functional_descriptions
################

#read in annotation descriptions
info=read.delim('TAIR10_functional_descriptions',head=T)
#read back in summary table
data=read.delim('te_trans_nearest_gene.bed',head=F)
#fix column names
colnames(data)=c(names(te.sums.merge),'gene_chr','gene_start','gene_end','geneID','type','strand','dot','class','distance')
#flag TE PAV that intersect with genes
data$intersect=ifelse(data$distance==0,1,0)
#get geneID from info set parsed
info$gene=matrix(unlist(strsplit(as.character(info[,1]),'\\.')),ncol=2,byrow=T)[,1]
#merge te pav summary with annotation descriptions, keep all te pavs
data.merge=merge(data,info,by.x=c('geneID'),by.y=c('gene'),all.x=T)
#the merge duplicates anotations in which the gene has multiple isoforms. Try to keep things to the primary transcript as best as possible
data.merge.unique=unique(data.merge[,c(1:36,38,39:40)])
data.merge.unique.transonly=subset(data.merge.unique,data.merge.unique$punative_trans_te==1)
#write out summary table
write.table(data.merge.unique.transonly,'FIGURE6_SUPPLEMENT_2.txt',sep='\t',row.names=F,quote=F)
