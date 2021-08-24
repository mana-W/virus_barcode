#grapth connectivity

#3 method

library(Seurat)
args=commandArgs(T)

rds<-readRDS(args[1])
if(length(grep("Low",rds@active.ident))>0){
  rds<-subset(rds,cells=row.names(rds@meta.data)[-grep("Low",rds@active.ident)])
}

#aver exp
av.exp <- AverageExpression(rds)$RNA
cor.exp <- as.data.frame(cor(av.exp))
hc = hclust(dist(cor.exp))
pdf("exp_cor_cluster.pdf",width = 6.5,height = 5.5)
plot(hc, hang = -1,col = "#487AA1", col.main = "#45ADA8", col.lab = "#7C8071", 
     col.axis = "#487AA1")
dev.off()
cor.exp$x <- rownames(cor.exp)
cor.df <- tidyr::gather(data = cor.exp, y, correlation,-x)
pdf("exp_cor.pdf",width = 6.5,height = 5.5)
ggplot(cor.df, aes(x, y, fill = correlation)) +
  scale_fill_gradient(low="#e2eff1",high="#233b6e")+
  geom_tile()+
  xlab("")+ylab("")+
  theme_bw()+theme(axis.text.x = element_text(angle = 90,vjust = 0.5))
dev.off()

#aver dist
umap_df<-data.frame(rds@reductions$umap@cell.embeddings)
umap_df$Type<-rds@active.ident
clus<-data.frame(clus=unique(umap_df$Type))
clus$x<-unlist(lapply(clus$clus,
                      function(x){mean(umap_df$UMAP_1[umap_df$Type==x])}
                      ))
clus$y<-unlist(lapply(clus$clus,
                      function(x){mean(umap_df$UMAP_2[umap_df$Type==x])}
))
clus_dist<-data.frame(as.matrix(dist(clus[,c(2,3)],method="euclidean",upper =T)))
row.names(clus_dist)<-clus$clus
colnames(clus_dist)<-clus$clus
hc = hclust(dist(clus_dist))
pdf("dist_cluster.pdf",width = 6.5,height = 5.5)
plot(hc, hang = -1,col = "#487AA1", col.main = "#45ADA8", col.lab = "#7C8071", 
     col.axis = "#487AA1")
dev.off()
clus_dist$x <- rownames(clus_dist)
cor.df <- tidyr::gather(data = data.frame(clus_dist), y, dist,-x)
pdf("dist_cor.pdf",width = 6.5,height = 5.5)
ggplot(cor.df, aes(x, y, fill = dist)) +
  scale_fill_gradient(high="#e2eff1",low="#233b6e")+
  geom_tile()+
  xlab("")+ylab("")+
  theme_bw()+theme(axis.text.x = element_text(angle = 90,vjust = 0.5))
dev.off()

#knn
library(kknn)
train<-umap_df[c(1,2),]
for(i in as.character(unique(umap_df$Type))){
  dat<-umap_df[sample(which(umap_df$Type==i),300,replace = T),]
  train<-rbind(train,dat)
}
test<-umap_df[c(1,2),]
for(i in as.character(unique(umap_df$Type))){
  dat<-umap_df[sample(which(umap_df$Type==i),300,replace = T),]
  test<-rbind(test,dat)
}
rm(dat)
pre1<- kknn(Type~., train, test,k=200, kernel = "rectangular")
fit <- fitted(pre1) 
table(fit,test$Type)
df<-data.frame(apply(table(fit,test$Type),2,function(x){x/sum(x)}))
hc = hclust(dist(df))
pdf("knn_cluster.pdf",width = 6.5,height = 5.5)
plot(hc, hang = -1,col = "#487AA1", col.main = "#45ADA8", col.lab = "#7C8071", 
     col.axis = "#487AA1")
dev.off()
df$x <- rownames(df)
cor.df <- tidyr::gather(data = data.frame(df), y, dist,-x)
pdf("knn_cor.pdf",width = 6.5,height = 5.5)
ggplot(cor.df, aes(x, y, fill = dist)) +
  scale_fill_gradient(high="#e2eff1",low="#233b6e")+
  geom_tile()+
  xlab("")+ylab("")+
  theme_bw()+theme(axis.text.x = element_text(angle = 90,vjust = 0.5))
dev.off()












