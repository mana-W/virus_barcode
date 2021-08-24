#input1 is rds 2 is "RNA" or "integrated"
#cluster must be int like : 0,1,2,3 ...

#SubMap input prepare

library(Seurat)
args<-commandArgs(T)

rds<-readRDS(args[1])

if(args[2]=="RNA"){
	dat<-rds@assays$RNA@scale.data
	}else{
	dat<-rds@assays$integrated@scale.data
	}
	
write.table(dat,"gene.tsv",quote=F,row.names = TRUE)
write.table(data.frame(rds@active.ident),"cluster.tsv",row.names = F,col.names = F,quote=F)

