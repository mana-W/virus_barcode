# input 1 is clonal threshold (0.5 / 0.7 / 1 ...)
# input 2 is celltype1 (F:/you/D30/clone_final) # title need
# input 3 is celltype2 (F:/you/D44/clone_final)

library(jaccard)
library(reshape2)
library(tidyr)
library(pheatmap)
library(ggplot2)

barcode_split_stat<-function(x){
  return(unlist(strsplit(x,";")))
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
  VBC2<-row.names(clone_stat_data)[clone_stat_data[,l]>0]
  return(length(intersect(VBC1,VBC2))/length(VBC1))
}
all_jac_stat1<-function(i,clone_stat_data,clust){
  VBC1<-row.names(clone_stat_data)[clone_stat_data[,i]>0]
  return(unlist(lapply(clust,all_jac_stat2,VBC1=VBC1,clone_stat_data=clone_stat_data)))
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

barplot_stat<-function(x,times){
  idx<-which(clu==x)
  dat<-t(do.call("cbind",jac_sample))[seq(idx,(500*length(clu)),by=length(clu)),times]
  dat<-data.frame(dat) %>% gather(time,similarity)
  dat_true<-data.frame(t(all_jac))[idx,times]
  dat_true<-data.frame(dat_true) %>% gather(time,similarity)
  p<-ggplot()+
    geom_boxplot(data=dat,aes(x=time,y=similarity),outlier.shape = NA,fill="slategray3",color=NA,fatten=NULL)+
    geom_point(data=dat_true,aes(time,y=similarity),color="hotpink3")+
    theme_bw()+xlab(x)+
    theme(axis.text.x = element_text(angle = 90,vjust = 0.5))
  return(p)
}

args=commandArgs(T)
clone_1<-read.table(args[2],header=T,stringsAsFactors = F)
clone_2<-read.table(args[3],header=T,stringsAsFactors = F)
clone<-rbind(clone_1,clone_2)
thresh=as.numeric(args[1])
prothresh = as.numeric(args[4])

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

clonetime1 <- unique(clone$time[clone$Cluster %in% clone_1$Cluster])
clonetime2 <- unique(clone$time[clone$Cluster %in% clone_2$Cluster])
clone1 <- clone[clone$time == clonetime1,]
tab1 <- acast(clone1,clone~Cluster)
tab1 <- t(apply(tab1,1,function(x){x/sum(x)}))
tab1 <- data.frame(tab1)
tab1$clone <- rownames(tab1)
tab1 <- melt(tab1)
names(tab1)[2]<-"Cluster"
tab1 <- tab1[tab1$value>prothresh,]

clone2 <- clone[clone$time == clonetime2,]
tab2 <- acast(clone2,clone~Cluster)
tab2 <- t(apply(tab2,1,function(x){x/sum(x)}))
tab2 <- data.frame(tab2)
tab2$clone <- rownames(tab2)
tab2 <- melt(tab2)
names(tab2)[2]<-"Cluster"
tab2 <- tab2[tab2$value>prothresh,]
taball <- rbind(tab1,tab2)
taball <- taball[,c(1,2)]
clone <- merge(clone,taball,by=c("clone","Cluster"))
clone <- clone[,c(3,4,5,2,6,1)]

timeclu<-unlist(lapply(unique(paste(clone$time,clone$clone,sep = "_")),function(x){strsplit(x,"_")[[1]][2]}))
clone<-clone[clone$clone %in%  as.character(data.frame(table(timeclu))$timeclu[data.frame(table(timeclu))$Freq>1]),]
write.csv(clone,"clone_2.csv",row.names = F,quote=F)



pdf("clone_size_log2.pdf",width = 4,height = 3)
plot(density(log2(table(clone$clone))))
dev.off()
time1clu<-names(table(clone$Cluster))[grep(names(table(clone$time))[1],names(table(clone$Cluster)))]
time2clu<-names(table(clone$Cluster))[grep(names(table(clone$time))[2],names(table(clone$Cluster)))]
clu<-names(clone_tab)

jac_sample<-list()
for(t in c(1:500)){
  for(i in names(table(clone$time))){
    clone$clone_sample[clone$time==i]<-sample(clone$clone[clone$time==i],length(clone$clone[clone$time==i]))
  }
  clone_tab_sample_sub<-data.frame(acast(clone,clone_sample~Cluster))
  clu<-names(clone_tab_sample_sub)
  jac_sample_sub<-data.frame(matrix(unlist(lapply(clu,all_jac_stat1,clone_stat_data=clone_tab_sample_sub,clust=clu)),ncol=length(clu)))
  names(jac_sample_sub)<-clu
  row.names(jac_sample_sub)<-clu
  jac_sample<-c(jac_sample,list(jac_sample_sub))
}
saveRDS(jac_sample,"jac_sample.rds")

#row names is denominator
clone_tab<-data.frame(acast(clone,clone~Cluster))
clu<-names(clone_tab)
all_jac<-data.frame(matrix(unlist(lapply(clu,all_jac_stat1,clone_stat_data=clone_tab,clust=clu)),ncol=length(clu)))
names(all_jac)<-clu
row.names(all_jac)<-clu
write.csv(all_jac,"all_jac.csv",quote=F)

all_p_1<-lapply(time1clu,barplot_stat,times=time2clu)
pdf("similarity_box_2_1.pdf",width = 3.5,height = 2)
for(p in all_p_1){
  print(p)
}
dev.off()
pdf("clust_2_1.pdf",width = 4,height = 4)
plot(hclust(dist(all_jac[time2clu,time1clu])))
dev.off()
all_p_2<-lapply(time2clu,barplot_stat,times=time1clu)
pdf("similarity_box_1_2.pdf",width = 3.5,height = 2)
for(p in all_p_2){
  print(p)
}
dev.off()
pdf("clust_1_2.pdf",width = 4,height = 4)
plot(hclust(dist(all_jac[time1clu,time2clu])))
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

annotation_row = data.frame(
  Time = factor(unlist(lapply(names(all_jac),function(x){strsplit(x,"_")[[1]][1]})))
)
rownames(annotation_row) = names(all_jac)

pheatmap(all_ob,cluster_rows = F,cluster_cols = F,annotation_col =annotation_row,annotation_row = annotation_row,color = rev(c(colorRampPalette(colors = c("#D73027","#FDAE61"))(200),colorRampPalette(colors =c("#FDAE61","white"))(30),colorRampPalette(colors = c("white","#4575B4"))(20))),file="all_obspre.pdf",width = 5,height = 4)
dev.new()
pheatmap(all_p,cluster_rows = F,cluster_cols = F,color = c(colorRampPalette(colors = c("#D73027","#FDAE61"))(20),colorRampPalette(colors =c("#FDAE61","white"))(30),colorRampPalette(colors = c("white","#4575B4"))(200)),annotation_col =annotation_row,annotation_row = annotation_row,file="all_pvalue.pdf",width = 5,height = 4)
dev.new()
pheatmap(all_jac,cluster_rows = F,cluster_cols = F,annotation_col =annotation_row,annotation_row = annotation_row,color = rev(c(colorRampPalette(colors = c("#D73027","#FDAE61"))(200),colorRampPalette(colors =c("#FDAE61","white"))(30),colorRampPalette(colors = c("white","#4575B4"))(20))),file="all_jaccard.pdf",width = 5,height = 4)
dev.new()


pheatmap(all_ob[time1clu,time2clu],cluster_rows = F,cluster_cols = F,color = rev(c(colorRampPalette(colors = c("#D73027","#FDAE61"))(200),colorRampPalette(colors =c("#FDAE61","white"))(30),colorRampPalette(colors = c("white","#4575B4"))(20))),file="data1_2_obspre.pdf",width = 5,height = 4)
dev.new()
pheatmap(all_p[time1clu,time2clu],cluster_rows = F,cluster_cols = F,color = c(colorRampPalette(colors = c("#D73027","#FDAE61"))(20),colorRampPalette(colors =c("#FDAE61","white"))(30),colorRampPalette(colors = c("white","#4575B4"))(200)),file="data1_2_pvalue.pdf",width = 5,height = 4)
dev.new()
pheatmap(all_jac[time1clu,time2clu],cluster_rows = F,cluster_cols = F,color = rev(c(colorRampPalette(colors = c("#D73027","#FDAE61"))(200),colorRampPalette(colors =c("#FDAE61","white"))(30),colorRampPalette(colors = c("white","#4575B4"))(20))),file="data1_2_jaccard.pdf",width = 5,height = 4)
dev.new()

pheatmap(all_ob[time2clu,time1clu],cluster_rows = F,cluster_cols = F,color = rev(c(colorRampPalette(colors = c("#D73027","#FDAE61"))(200),colorRampPalette(colors =c("#FDAE61","white"))(30),colorRampPalette(colors = c("white","#4575B4"))(20))),file="data2_1_obspre.pdf",width = 5,height = 4)
dev.new()
pheatmap(all_p[time2clu,time1clu],cluster_rows = F,cluster_cols = F,color = c(colorRampPalette(colors = c("#D73027","#FDAE61"))(20),colorRampPalette(colors =c("#FDAE61","white"))(30),colorRampPalette(colors = c("white","#4575B4"))(200)),file="data2_1_pvalue.pdf",width = 5,height = 4)
dev.new()
pheatmap(all_jac[time2clu,time1clu],cluster_rows = F,cluster_cols = F,color = rev(c(colorRampPalette(colors = c("#D73027","#FDAE61"))(200),colorRampPalette(colors =c("#FDAE61","white"))(30),colorRampPalette(colors = c("white","#4575B4"))(20))),file="data2_1_jaccard.pdf",width = 5,height = 4)
dev.new()

all_jac$Cluster<-row.names(all_jac)
all_jac<-gather(all_jac,Cluster2,jaccard,-Cluster)
all_jac<-unique(all_jac)

all_p$Cluster<-row.names(all_p)
all_p<-gather(all_p,Cluster2,pvalue,-Cluster)
all_p<-unique(all_p)
all_p$pvalue[is.na(all_p$pvalue)]<-0

dat<-merge(all_jac,all_p,by=c("Cluster","Cluster2"))
pdf("pop_all_data2_1.pdf",width = 5.5,height = 4)
ggplot(dat[dat$Cluster %in% time2clu & dat$Cluster2 %in% time1clu,],aes(Cluster,Cluster2))+
  geom_point(aes(size=jaccard,color=-log2(pvalue+0.0001)))+
  scale_colour_gradient(low="white",high="brown4")+
  labs(color=expression(-log2(pvalue+0.0001)),size="jaccard similarity",x="",y="",title="")+
  scale_size(limits = c(0.0001,max(dat$jaccard)))+
  theme_bw()+theme(axis.text.x = element_text(angle = 90,vjust = 0.5))
dev.off()
pdf("pop_all_data1_2.pdf",width = 5.5,height = 4)
ggplot(dat[dat$Cluster %in% time1clu & dat$Cluster2 %in% time2clu,],aes(Cluster,Cluster2))+
  geom_point(aes(size=jaccard,color=-log2(pvalue+0.0001)))+
  scale_colour_gradient(low="white",high="brown4")+
  labs(color=expression(-log2(pvalue+0.0001)),size="jaccard similarity",x="",y="",title="")+
  scale_size(limits = c(0.0001,max(dat$jaccard)))+
  theme_bw()+theme(axis.text.x = element_text(angle = 90,vjust = 0.5))
dev.off()

dat$jaccard[dat$pvalue>0.05]<-0
dat$pvalue[dat$pvalue>0.05]<-1

pdf("pop_sig_data2_1.pdf",width = 5.5,height = 4)
ggplot(dat[dat$Cluster %in% time2clu & dat$Cluster2 %in% time1clu,],aes(Cluster,Cluster2))+
  geom_point(aes(size=jaccard,color=-log2(pvalue+0.0001)))+
  scale_colour_gradient(low="white",high="brown4")+
  labs(color=expression(-log2(pvalue+0.0001)),size="jaccard similarity",x="",y="",title="")+
  scale_size(limits = c(0.0001,max(dat$jaccard)))+
  theme_bw()+theme(axis.text.x = element_text(angle = 90,vjust = 0.5))
dev.off()
pdf("pop_sig_data1_2.pdf",width = 5.5,height = 4)
ggplot(dat[dat$Cluster %in% time1clu & dat$Cluster2 %in% time2clu,],aes(Cluster,Cluster2))+
  geom_point(aes(size=jaccard,color=-log2(pvalue+0.0001)))+
  scale_colour_gradient(low="white",high="brown4")+
  labs(color=expression(-log2(pvalue+0.0001)),size="jaccard similarity",x="",y="",title="")+
  scale_size(limits = c(0.0001,max(dat$jaccard)))+
  theme_bw()+theme(axis.text.x = element_text(angle = 90,vjust = 0.5))
dev.off()


