library(Seurat)
args<-commandArgs(T)
rds<-readRDS(args[1])
genes<-VariableFeatures(rds)
Idents(rds)<-"type"
diff<-FindMarkers(rds,ident.1="cell_neu_bias",ident.2="cell_prog_bias")
diff<-diff[which(row.names(diff)%in% genes),]
diff$p.adj.fdr<-p.adjust(diff$p_val,method ="fdr",n=dim(diff)[1])
diff$p.adj.bonferroni<-p.adjust(diff$p_val,method ="bonferroni",n=dim(diff)[1])
diff$type<-"cell_neu_bias"
diff$type[diff$avg_log2FC<0]<-"cell_prog_bias"
diff$genes<-row.names(diff)
write.csv(diff,"diff.csv",row.names = F,quote=F)
