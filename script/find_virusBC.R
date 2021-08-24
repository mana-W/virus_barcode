# find Virus barcode 
# input reads file & barcode fa & cell type file 
# reads file header is Read.ID Read.Seq Cell.BC UMI
# cell type header : Cell.BC Cluster

# need to in put virus.fa (random seq : NNNNNNNNNN...)

source("function.R")

args=commandArgs(T)
command = args[1]
reads=args[2]
BCseq=args[3]

mat  =  nucleotideSubstitutionMatrix(match = 1, mismatch = -3)
BCseq=ReadFasta(BCseq)
BClength=str_count(as.character(BCseq$VBC),"N")
reads<-read.table(reads,header=T,stringsAsFactors=F)
alig_score<-(length(BCseq)-floor(length(BCseq)/10)-(3*floor(length(BCseq)/10))-3*BClength)
#environment() <- .GlobalEnv

#find virus barcode and create output file
print("Finding barcode ...")
library(parallel)
cl <- makeCluster(16)
clusterEvalQ(cl,library(Biostrings))
clusterEvalQ(cl,library(stringr))
clusterExport(cl,c('BCseq','BClength','alig_score','mat','find_barcode','reads'))
barcode<-parLapply(cl,reads$Read.Seq,find_barcode)
stopCluster(cl)
reads$Virus.BC<-unlist(barcode)
reads<-na.omit(reads)
reads<-reads[reads$Virus.BC!="unknown",]
reads$VBC_len<-nchar(reads$Virus.BC)
pdf("reads_VBC_length.pdf",width=3,height=3)
plot(density(reads$VBC_len),main="Length density")
dev.off()
reads<-reads[reads$VBC_len <= (BClength+5) & reads$VBC_len >= (BClength-5),]
write.table(reads,"CB_UMI_barcode.tsv",row.names=F,quote=F)
print("Finished")

#filter noise & confirm true barcode
print("Filter reads ...")
print("raw")
print(c("reads :",dim(reads)[1]))
print(c("cell barcodes :",length(unique(as.character(reads$Cell.BC)))))

passreads<-system(
		"sed '1d' CB_UMI_barcode.tsv | awk '{print$3$4$5}' | sort | uniq -c | awk '$1>=5 {print $2}'",
		intern = TRUE
		)
reads<-reads[paste0(reads$Cell.BC,reads$UMI,reads$Virus.BC) %in% passreads,]
print("more than 5 reads")
print(c("reads :",dim(reads)[1]))
print(c("cell barcodes :",length(unique(as.character(reads$Cell.BC)))))
write.table(reads,"CB_UMI_barcode_reads.tsv",row.names=F,quote=F)

passcells<-system(
		"awk 'NR>1 {print $3,$4}' CB_UMI_barcode_reads.tsv | sort -u | awk '{print$1}' | uniq -c | awk '$1>1 {print$2}'",
		intern = TRUE
		)
reads<-reads[reads$Cell.BC %in% passcells,]
print("cells with more than 2 UMIs")
print(c("reads :",dim(reads)[1]))
print(c("cell barcodes :",length(unique(as.character(reads$Cell.BC)))))
write.table(reads,"CB_UMI_barcode_UMI.tsv",row.names=F,quote=F)
system("awk 'NR>1 {print$3$5}' CB_UMI_barcode_UMI.tsv > CB_VB_for_cluster.tsv")
reads$CV<-paste0(reads$Cell.BC,reads$Virus.BC)

#barcode cluster (based on similarity)
system("starcode -d 5 --print-clusters -s CB_VB_for_cluster.tsv > CB_VB_cluster.tsv")
cv_clu<-read.table("CB_VB_cluster.tsv")
cbvb<-apply(cv_clu,1,split_clu)
cbvb<-data.frame(matrix(unlist(cbvb),ncol=2,byrow =T))
j_cbvb<-data.frame(matrix(unlist(apply(cbvb,1,estimate_clu)),ncol=9,byrow=T))
j_cbvb$rcbvb<-apply(j_cbvb,1,r_cbvb_stat)
j_cbvb<-j_cbvb[,c(1,10)]
names(j_cbvb)<-c("CV","rCV")
reads<-merge(reads,j_cbvb,by="CV",all.x=T)
newdata<-reads[,c(5,8)]
newdata$Cell.BC<-substring(newdata$rCV,1,16)
newdata$Virus.BC<-substring(newdata$rCV,17,)
print("cluster")
print(c("reads :",dim(newdata)[1]))
print(c("cell barcodes :",length(unique(as.character(newdata$Cell.BC)))))
write.table(newdata,"real_CV.tsv",row.names=F,sep="\t",quote=F)

#reform clone and match celltye
celltype=args[4]
celltype<-read.table(celltype,header=T,stringsAsFactors=F)
newdata<-newdata[,c(3,4)]
newdata<-newdata[!duplicated(newdata),]
clone_raw<-data.frame(Cell.BC=unique(as.character(newdata$Cell.BC)),Virus.BC=unlist(lapply(unique(as.character(newdata$Cell.BC)),clone_find)))
clone_raw$num<-unlist(lapply(as.character(clone_raw$Virus.BC),function(x){length(strsplit(x,split=";")[[1]])}))
write.table(clone_raw,"clone_raw.tsv",row.names=F,quote=F)
clone<-merge(clone_raw,celltype,by="Cell.BC")
print(c("final cells :",length(clone$Cell.BC)))
write.table(clone,"clone_final.tsv",row.names=F,quote=F)

#calculate similarity (with in one group)
if(length(args)<5){
	thresh=1
}else{
	thresh=as.numeric(args[5])
}

clone_tab<-data.frame(table(clone$Virus.BC))
clone_tab$Var1<-as.character(clone_tab$Var1)
clone_tab$num<-clone$num[match(clone_tab$Var1,clone$Virus.BC)]
clone_tab<-clone_tab[order(clone_tab$Freq,clone_tab$num),]
#clone$time<-unlist(lapply(clone$Cluster,function(x){strsplit(x,"_")[[1]][1]}))

barcode_split<-lapply(as.character(clone_tab$Var1),barcode_split_stat)

all_length<-length(barcode_split)
library(parallel)
cl <- makeCluster(4)
clusterEvalQ(cl,library(jaccard))
clusterExport(cl,c('clone_tab','barcode_split','all_length','jaccard_stat'))
jac_matirx<-parLapply(cl,clone_tab$Var1,jaccard_stat)
stopCluster(cl)

names(jac_matirx)<-clone_tab$Var1
clone_names<-names(jac_matirx)

for(n in c(1:(length(jac_matirx)-1))){
  i=data.frame(jac_matirx[[n]])
  if(max(i$dist)<thresh){
    next;
  }else{
    ind<-which(i$dist==max(i$dist))
    test<-i[ind,]
    ind<-max(as.numeric(row.names(test)))
    clone_names[ind+n]<-paste0(clone_names[ind+n],"_",clone_names[n])
    clone_names[n]<-"None"
  }
}
rm(n,i,ind,test)

clone_names<-lapply(c(1:length(clone_names)),change_form_stat)
clone_names<-do.call("rbind",clone_names)
clone<-merge(clone,clone_names,by="Virus.BC")
write.csv(clone,"clone.csv",row.names = F,quote=F)
clone<-clone[clone$clone %in%  as.character(data.frame(table(clone$clone))$Var1[data.frame(table(clone$clone))$Freq>1]),]
write.csv(clone,"clone_2.csv",row.names = F,quote=F)

pdf("clone_size_log2.pdf",width = 4,height = 3)
plot(density(log2(table(clone$clone))))
dev.off()

jac_sample<-list()
for(t in c(1:500)){
  clone$clone_sample<-sample(clone$clone,length(clone$clone))
  clone_tab_sample_sub<-data.frame(acast(clone,clone_sample~Cluster))
  clu<-names(clone_tab_sample_sub)
  jac_sample_sub<-data.frame(matrix(unlist(lapply(clu,all_jac_stat1,clone_stat_data=clone_tab_sample_sub)),ncol=length(clu)))
  names(jac_sample_sub)<-clu
  row.names(jac_sample_sub)<-clu
  jac_sample<-c(jac_sample,list(jac_sample_sub))
}
saveRDS(jac_sample,"jac_sample.rds")

clone_tab<-data.frame(acast(clone,clone~Cluster))
clu<-names(clone_tab)
all_jac<-data.frame(matrix(unlist(lapply(clu,all_jac_stat1,clone_stat_data=clone_tab)),ncol=length(clu)))
names(all_jac)<-clu
row.names(all_jac)<-clu
write.csv(all_jac,"all_jac.csv",quote=F)

all_p<-lapply(clu,barplot_stat)
pdf("similarity_box.pdf",width = 3.5,height = 2)
for(p in all_p){
  print(p)
}
dev.off()

pdf("clust.pdf",width = 4,height = 4)
plot(hclust(dist(all_jac)))
dev.off()

all_ob_p<-lapply(clu,all_op_stat1,jac=all_jac,jac_sample=jac_sample,clu=clu)
all_ob<-do.call("rbind",lapply(all_ob_p,function(x){unlist(x)[c(1:length(unlist(x)))[c(1:length(unlist(x)))%%2==1]]}))
all_ob<-data.frame(all_ob)
all_p<-do.call("rbind",lapply(all_ob_p,function(x){unlist(x)[c(1:length(unlist(x)))[c(1:length(unlist(x)))%%2==0]]}))
all_p<-data.frame(all_p)
names(all_ob)<-clu
row.names(all_ob)<-clu
names(all_p)<-clu
row.names(all_p)<-clu

write.csv(all_ob,"all_obpre.csv",quote=F)
write.csv(all_p,"all_pvalue.csv",quote=F)

pheatmap(all_ob,cluster_rows = F,cluster_cols = F,color = rev(c(colorRampPalette(colors = c("#D73027","#FDAE61"))(200),colorRampPalette(colors =c("#FDAE61","white"))(30),colorRampPalette(colors = c("white","#4575B4"))(20))),file="obspre.pdf",width = 3.3,height = 3)
dev.new()
pheatmap(all_p,cluster_rows = F,cluster_cols = F,color = c(colorRampPalette(colors = c("#D73027","#FDAE61"))(20),colorRampPalette(colors =c("#FDAE61","white"))(30),colorRampPalette(colors = c("white","#4575B4"))(200)),file="pvalue.pdf",width = 3.3,height = 3)
dev.new()
pheatmap(all_jac,cluster_rows = F,cluster_cols = F,color = rev(c(colorRampPalette(colors = c("#D73027","#FDAE61"))(200),colorRampPalette(colors =c("#FDAE61","white"))(30),colorRampPalette(colors = c("white","#4575B4"))(20))),file="jaccard.pdf",width = 3.3,height = 3)
dev.new()

all_jac$Cluster<-row.names(all_jac)
all_jac<-gather(all_jac,Cluster2,jaccard,-Cluster)
all_jac<-unique(all_jac)

all_p$Cluster<-row.names(all_p)
all_p<-gather(all_p,Cluster2,pvalue,-Cluster)
all_p<-unique(all_p)
all_p$pvalue[is.na(all_p$pvalue)]<-0

dat<-merge(all_jac,all_p,by=c("Cluster","Cluster2"))

pdf("pop_all.pdf",width = 5.5,height = 4)
ggplot(dat,aes(Cluster,Cluster2))+
  geom_point(aes(size=jaccard,color=-log2(pvalue+0.0001)))+
  scale_colour_gradient(low="white",high="brown4")+
  labs(color=expression(-log2(pvalue+0.0001)),size="jaccard similarity",x="",y="",title="")+
  scale_size(limits = c(0.0001,max(dat$jaccard)))+
  theme_bw()+theme(axis.text.x = element_text(angle = 90,vjust = 0.5))
dev.off()

dat$jaccard[dat$pvalue>0.05]<-0
dat$pvalue[dat$pvalue>0.05]<-1

pdf("pop_sig.pdf",width = 5.5,height = 4)
ggplot(dat,aes(Cluster,Cluster2))+
  geom_point(aes(size=jaccard,color=-log2(pvalue+0.0001)))+
  scale_colour_gradient(low="white",high="brown4")+
  labs(color=expression(-log2(pvalue+0.0001)),size="jaccard similarity",x="",y="",title="")+
  scale_size(limits = c(0.0001,max(dat$jaccard)))+
  theme_bw()+theme(axis.text.x = element_text(angle = 90,vjust = 0.5))
dev.off()

































