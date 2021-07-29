# input 1 is clonal threshold (0.5 / 0.7 / 1 ...)
# input 2 is celltype1 (F:/you/D30/clone_final) # title need
# input 3 is whether need filtering at first

library(jaccard)
library(reshape2)
library(tidyr)
library(pheatmap)
library(ggplot2)

barcode_split_stat<-function(x){
  return(unlist(strsplit(x,";")))
}

sample_p_stat2<-function(y,x,num,clone_sample){
  clone_real<-num[x,y]
  clone_pred<-as.numeric(unlist(lapply(clone_sample,function(z){z[x,y]})))
  zscore<-(clone_real-mean(clone_pred))/sd(clone_pred)
  p<-pnorm(zscore,lower.tail = F)
  return(c(p))
}
sample_p_stat1<-function(x,num,clone_sample,clu){
  #x is clone_name
  return(lapply(clu,sample_p_stat2,x=x,num=num,clone_sample=clone_sample))
}

jaccard_stat<-function(x){
  ind<-which(clone_tab$Var1==x)
  jac_dist<-unlist(lapply(barcode_split[(ind+1):all_length],function(i){length(intersect(barcode_split[[ind]],i))/length(union(barcode_split[[ind]],i))}))
  jac_clonenum<-clone_tab$num[(ind+1):all_length]
  jac_cellnum<-clone_tab$Freq[(ind+1):all_length]
  return(data.frame(cell_num=jac_cellnum,tag_num=jac_clonenum,dist=jac_dist))
}

change_form_stat<-function(x){
  VBC=unlist(strsplit(clone_names[x],"_"))
  return(data.frame(Virus.BC=VBC,clone=rep(x,length(VBC))))
}

all_jac_stat2<-function(l,VBC1,clone_stat_data){
  VBC2<-clone_stat_data[,l]
  VBC2[VBC2>0]<-1
  return(jaccard(VBC1,VBC2))
}
all_jac_stat1<-function(i,clone_stat_data){
  VBC1<-clone_stat_data[,i]
  VBC1[VBC1>0]<-1
  return(unlist(lapply(clu,all_jac_stat2,VBC1=VBC1,clone_stat_data=clone_stat_data)))
}

all_op_stat2<-function(y,x,jac,jac_sample){
  jac_real<-jac[x,y]
  jac_pred<-as.numeric(unlist(lapply(jac_sample,function(z){z[x,y]})))
  ob<-jac_real/mean(jac_pred)
  zscore<-(jac_real-mean(jac_pred))/sd(jac_pred)
  p<-pnorm(zscore,lower.tail = F)
  return(c(ob,p))
}
all_op_stat1<-function(x,jac,jac_sample,clu){
  return(lapply(clu,all_op_stat2,x=x,jac=jac,jac_sample=jac_sample))
}

barplot_stat<-function(x){
  idx<-which(clu==x)
  dat<-do.call("rbind",jac_sample)[seq(idx,(500*length(clu)),by=length(clu)),-idx]
  dat<-dat %>% gather(time,similarity)
  dat_true<-all_jac[idx,-idx]
  dat_true<-dat_true %>% gather(time,similarity)
  p<-ggplot()+
    geom_boxplot(data=dat,aes(x=time,y=similarity),outlier.shape = NA,fill="slategray3",color=NA,fatten=NULL)+
    geom_point(data=dat_true,aes(time,y=similarity),color="hotpink3")+
    theme_bw()+xlab(x)+
    theme(axis.text.x = element_text(angle = 90,vjust = 0.5))
  return(p)
}


args=commandArgs(T)
clone<-read.table(args[2],header=T,stringsAsFactors = F)
thresh=as.numeric(args[1])
Filt=args[3]

clone_tab<-data.frame(table(clone$Virus.BC))
clone_tab$Var1<-as.character(clone_tab$Var1)
clone_tab$num<-clone$num[match(clone_tab$Var1,clone$Virus.BC)]
clone_tab<-clone_tab[order(clone_tab$Freq,clone_tab$num),]
clone$time<-unlist(lapply(clone$Cluster,function(x){strsplit(x,"_")[[1]][1]}))

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
#clone_names<-clone_names[grep("_",clone_names,ignore.case = T)]
clone_names<-lapply(c(1:length(clone_names)),change_form_stat)
clone_names<-do.call("rbind",clone_names)
clone<-merge(clone,clone_names,by="Virus.BC")
write.csv(clone,"clone.csv",row.names = F,quote=F)
clone<-clone[clone$clone %in%  as.character(data.frame(table(clone$clone))$Var1[data.frame(table(clone$clone))$Freq>1]),]


if(Filt=="T"){
  #clone<-clone[clone$clone %in%  as.character(data.frame(table(clone$clone))$Var1[data.frame(table(clone$clone))$Freq>1]),]
  write.csv(clone,"clone_2_raw.csv",row.names = F,quote=F)
  clone_tab<-data.frame(acast(clone,clone~Cluster))
  clone_sample<-list()
  for(t in c(1:500)){
    for(i in names(table(clone$time))){
      clone$clone_sample[clone$time==i]<-sample(clone$clone[clone$time==i],length(clone$clone[clone$time==i]))
    }
    clone_tab_sample_sub<-data.frame(acast(clone,clone_sample~Cluster))
    clone_sample<-c(clone_sample,list(clone_tab_sample_sub))
  }
  clone_name<-row.names(clone_tab)
  clone_p<-lapply(clone_name,sample_p_stat1,num=clone_tab,clone_sample=clone_sample,clu=names(clone_tab))
  all_p<-data.frame(do.call("rbind",lapply(clone_p,function(x){unlist(x)[c(1:length(unlist(x)))]})))
  row.names(all_p)<-clone_name
  colnames(all_p)<-names(clone_tab)
  all_p$clone<-row.names(all_p)
  all_p<-gather(all_p,Cluster,pvalue,-clone)
  clone<-merge(clone,all_p,by=c("Cluster","clone"))
  clone<-clone[clone$pvalue<=0.05,c("Virus.BC","Cell.BC","num","Cluster","time","clone")]
  clone<-clone[clone$clone %in%  as.character(data.frame(table(clone$clone))$Var1[data.frame(table(clone$clone))$Freq>1]),]
}

write.csv(clone,"clone_2.csv",row.names = F,quote=F)

pdf("clone_size_log2.pdf",width = 4,height = 3)
plot(density(log2(table(clone$clone))))
dev.off()

jac_sample<-list()
for(t in c(1:500)){
  for(i in names(table(clone$time))){
    clone$clone_sample[clone$time==i]<-sample(clone$clone[clone$time==i],length(clone$clone[clone$time==i]))
  }
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

pdf("clust_p.pdf",width = 4,height = 4)
plot(hclust(dist(all_p)))
dev.off()

write.csv(all_ob,"all_obpre.csv",quote=F)
write.csv(all_p,"all_pvalue.csv",quote=F)

#annotation_row = data.frame(
#  Time = factor(unlist(lapply(names(all_jac),function(x){strsplit(x,"_")[[1]][1]})))
#)
#rownames(annotation_row) = names(all_jac)


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


