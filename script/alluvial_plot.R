#input 1-2 celltype 3-4 pvalue & similarity
library(ggplot2)
library(ggalluvial)
library(tidyr)

args=commandArgs(T)
celltype_1<-read.table(args[1],header=T,stringsAsFactors = F)
celltype_2<-read.table(args[2],header=T,stringsAsFactors = F)
pvalue_dat<-read.csv(args[3],header=T,row.names = 1)
simi_dat<-read.csv(args[4],header=T,row.names = 1)

simi_dat$Cluster<-row.names(simi_dat)
simi_dat<-gather(simi_dat,Cluster2,jaccard,-Cluster)
simi_dat<-unique(simi_dat)

pvalue_dat$Cluster<-row.names(pvalue_dat)
pvalue_dat<-gather(pvalue_dat,Cluster2,pvalue,-Cluster)
pvalue_dat<-unique(pvalue_dat)
pvalue_dat$pvalue[is.na(pvalue_dat$pvalue)]<-0

dat<-merge(simi_dat,pvalue_dat,by=c("Cluster","Cluster2"))
dat<-dat[unlist(lapply(dat$Cluster,function(x){strsplit(x,"_")[[1]][1]}))!=unlist(lapply(dat$Cluster2,function(x){strsplit(x,"_")[[1]][1]})),]

data<-dat[1,]
data<-data[-1,]
for(i in which(dat$pvalue<=0.05)){
  clu=dat$Cluster[i]
  clu2=dat$Cluster2[i]
  data<-rbind(data,dat[dat$Cluster==clu & dat$Cluster2==clu2,])
  data<-rbind(data,dat[dat$Cluster==clu2 & dat$Cluster2==clu,])
}
data<-unique(data)
data<-data[order(data$Cluster2),]

celltype_1<-data.frame(table(celltype_1$Cluster)/dim(celltype_1)[1])
celltype_2<-data.frame(table(celltype_2$Cluster)/dim(celltype_2)[1])
celltype_1$Freq<-300*(celltype_1$Freq)
celltype_2$Freq<-300*(celltype_2$Freq)
if else
celltype_1$Freq<-30
celltype_2$Freq<-30

data$group<-unlist(lapply(data$Cluster2,function(x){strsplit(x,"_")[[1]][1]}))
which(duplicated(data[,c(1,2)]))
data$connect<-NA
data$connect[c(1:(dim(data)[1]/2))]<-c(1:(dim(data)[1]/2))

for(i in which(!is.na(data$connect))){
  clu1<-data$Cluster[i]
  clu2<-data$Cluster2[i]
  data$connect[data$Cluster==clu2 & data$Cluster2==clu1]<-data$connect[i]
}

for(clus in unique(data$Cluster2)){
  data$jaccard[data$Cluster2==clus]<-data$jaccard[data$Cluster2==clus]/sum(data$jaccard[data$Cluster2==clus])
}

data<-data[,c(2,3,5,6)]
names(data)[1]<-"name"

celltype<-rbind(celltype_1,celltype_2)
names(celltype)[1]<-"name"

data<-rbind(data,data.frame(name=setdiff(celltype$name,data$name),
                            jaccard=1,
                            group=unlist(lapply(setdiff(celltype$name,data$name),function(x){strsplit(x,"_")[[1]][1]})),
                            connect=c(1:length(setdiff(celltype$name,data$name)))+50
                            ))

data$value<-NA
for(clus in celltype$name){
  data$value[data$name==clus]=celltype$Freq[celltype$name==clus]*data$jaccard[data$name==clus]
}

pdf("sanky.pdf",width = 6,height = 5)
ggplot(data = data,aes(x = group, stratum = name, alluvium = connect,y = value,label=name))+
  geom_stratum(fill="lightcyan3")+
  geom_flow()+
  geom_text(stat = "stratum",size=2.5)+
  ylab("")+
  theme_bw()+
  NoLegend()+
  scale_y_discrete(breaks=NULL)
dev.off()
