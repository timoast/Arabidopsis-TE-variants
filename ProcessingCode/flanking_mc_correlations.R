# From Steve's RMarkdown script
# Won't run, needs some input data files and paths changed, but shows what was done to get the correlation values

ins.met=read.delim('insertion_mc.tsv',head=T)
ins.met.filt=ins.met[,intersect(grep('Ag_0_bud',names(ins.met),invert=T),intersect(grep('Bor_4_bud',names(ins.met),invert=T),intersect(grep('Ema_1_bud',names(ins.met),invert=T),intersect(grep('Er_0_bud',names(ins.met),invert=T),intersect(grep('Fr_2_bud',names(ins.met),invert=T),intersect(grep('Is_0_bud',names(ins.met),invert=T),intersect(grep('Kin_0_bud',names(ins.met),invert=T),intersect(grep('Litva_bud',names(ins.met),invert=T),intersect(grep('Pu2_23_bud',names(ins.met),invert=T),intersect(grep('Ragl_1_bud',names(ins.met),invert=T),intersect(grep('Uk_1_bud',names(ins.met),invert=T),intersect(grep('Zdr_1_bud',names(ins.met),invert=T),intersect(grep('vi_0',names(ins.met),invert=T),intersect(grep('HR_5',names(ins.met),invert=T),grep('Ler_1',names(ins.met),invert=T)))))))))))))))]
names(ins.met.filt)=gsub("_bud","",names(ins.met.filt))

insertions.matrix=read.delim('../1_file_formatting/TE_insertion_matrix_genotyped.txt',head=T)
insertions.matrix.subset=insertions.matrix[,c(1:7,8,9,10,11,13,14,15,17,18,19,20,21,22,25,28,29,30,31,33,34,35,36,37,39,40,42,43,46,48,49,50,51,52,54,55,63,64,65,68,69,70,71,72,74,75,76,78,79,80,82,83,84,85,86,87,89,95,97,98,99,101,103,106,108,109,110,112,115,116,117,118,119,120,121,123,126,128,129,131,132,135,137,139,141,144,145,146,148,149,151,152,155,156,157,158,159,160,161,162,179,165,166,168,171,173,174,176,180,182,183,184,186,187,188,189,191,192,193,194,195,196,200,203,204,206,207,208,210,211,213,214,217,219,220,221,222,224)]

namecheck=cbind(names(ins.met.filt)[4:140],names(insertions.matrix.subset[,8:144]))
namecheck[,2]=gsub("\\.","_",namecheck[,2])
table(ifelse(namecheck[,1]==namecheck[,2],T,F))

namecheck[ifelse(namecheck[,1]==namecheck[,2],F,T),]
namecheck[,2]=gsub('Da_1__12','Da_1_12',namecheck[,2])
namecheck[,2]=gsub('RMX_A02','Rmx_A02',namecheck[,2])
namecheck[ifelse(namecheck[,1]==namecheck[,2],F,T),]

colnames(insertions.matrix.subset)[8:144]=namecheck[,2]

insertions.matrix.subset$filt_poly_count=rowSums(insertions.matrix.subset[,8:144],na.rm=T)
insertions.matrix.subset$filt_nopoly_count=rowSums(insertions.matrix.subset[,8:144]==0,na.rm=T)
insertions.matrix.subset$filt_na_count=rowSums(is.na(insertions.matrix.subset[,8:144])==T,na.rm=T)
insertions.matrix.subset$filt_min_freq=apply(insertions.matrix.subset[,145:146],1,min)
insertions.matrix.subset$maf=insertions.matrix.subset$filt_min_freq/(insertions.matrix.subset[,145]+insertions.matrix.subset[,146])

insertions.matrix.subset.maf=subset(insertions.matrix.subset,insertions.matrix.subset$maf>=0.03)
ins.met.filt.maf=subset(ins.met.filt,insertions.matrix.subset$maf>=0.03)

correlation.values=rep(NA,nrow(insertions.matrix.subset.maf))
for(i in 1:nrow(insertions.matrix.subset.maf)){
  correlation.values[i]=cor(unlist(ins.met.filt.maf[i,4:140]), unlist(insertions.matrix.subset.maf[i,8:144]),use="pairwise.complete.obs")
}

out=cbind(insertions.matrix.subset.maf,correlation.values)
write.table("../ProcessedData/correlation_insertions_regenerated.tsv")

#############################

deletions.matrix=read.delim('../1_file_formatting/TE_deletion_matrix_genotyped.txt',head=T)
deletions.matrix=deletions.matrix[,-4]
deletions.matrix.subset=deletions.matrix[,c(1:7,8,9,10,11,13,14,15,17,18,19,20,21,22,25,28,29,30,31,33,34,35,36,37,39,40,42,43,46,48,49,50,51,52,54,55,63,64,65,68,69,70,71,72,74,75,76,78,79,80,82,83,84,85,86,87,89,95,97,98,99,101,103,106,108,109,110,112,115,116,117,118,119,120,121,123,126,128,129,131,132,135,137,139,141,144,145,146,148,149,151,152,155,156,157,158,159,160,161,162,179,165,166,168,171,173,174,176,180,182,183,184,186,187,188,189,191,192,193,194,195,196,200,203,204,206,207,208,210,211,213,214,217,219,220,221,222,224)]

del.met=read.delim('deletions_mc.tsv',head=T)
del.met.filt=del.met[,intersect(grep('Ag_0_bud',names(del.met),invert=T),intersect(grep('Bor_4_bud',names(del.met),invert=T),intersect(grep('Ema_1_bud',names(del.met),invert=T),intersect(grep('Er_0_bud',names(del.met),invert=T),intersect(grep('Fr_2_bud',names(del.met),invert=T),intersect(grep('Is_0_bud',names(del.met),invert=T),intersect(grep('Kin_0_bud',names(del.met),invert=T),intersect(grep('Litva_bud',names(del.met),invert=T),intersect(grep('Pu2_23_bud',names(del.met),invert=T),intersect(grep('Ragl_1_bud',names(del.met),invert=T),intersect(grep('Uk_1_bud',names(del.met),invert=T),intersect(grep('Zdr_1_bud',names(del.met),invert=T),intersect(grep('vi_0',names(del.met),invert=T),intersect(grep('HR_5',names(del.met),invert=T),grep('Ler_1',names(del.met),invert=T)))))))))))))))]
names(del.met.filt)=gsub("_bud","",names(del.met.filt))

namecheck=cbind(names(del.met.filt)[4:140],names(deletions.matrix.subset[,8:144]))
namecheck[,2]=gsub("\\.","_",namecheck[,2])
table(ifelse(namecheck[,1]==namecheck[,2],T,F))
namecheck[ifelse(namecheck[,1]==namecheck[,2],F,T),]

namecheck[,2]=gsub('Da_1__12','Da_1_12',namecheck[,2])
namecheck[,2]=gsub('RMX_A02','Rmx_A02',namecheck[,2])
namecheck[ifelse(namecheck[,1]==namecheck[,2],F,T),]
colnames(deletions.matrix.subset)[8:144]=namecheck[,2]

deletions.matrix.subset$filt_poly_count=rowSums(deletions.matrix.subset[,8:144],na.rm=T)
deletions.matrix.subset$filt_nopoly_count=rowSums(deletions.matrix.subset[,8:144]==0,na.rm=T)
deletions.matrix.subset$filt_na_count=rowSums(is.na(deletions.matrix.subset[,8:144])==T,na.rm=T)
deletions.matrix.subset$filt_min_freq=apply(deletions.matrix.subset[,145:146],1,min)
deletions.matrix.subset$maf=deletions.matrix.subset$filt_min_freq/(deletions.matrix.subset[,145]+deletions.matrix.subset[,146])

deletions.matrix.subset.maf=subset(deletions.matrix.subset,deletions.matrix.subset$maf>=0.03)
del.met.filt.maf=subset(del.met.filt,deletions.matrix.subset$maf>=0.03)
correlation.values=rep(NA,nrow(deletions.matrix.subset.maf))
for(i in 1:nrow(deletions.matrix.subset.maf)){
  correlation.values[i]=cor(unlist(del.met.filt.maf[i,4:140]), unlist(deletions.matrix.subset.maf[i,8:144]),use="pairwise.complete.obs")
}

out=cbind(deletions.matrix.subset.maf,correlation.values)
write.table("../ProcessedData/correlation_deletions_regenerated.tsv")
