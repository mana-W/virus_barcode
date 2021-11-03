#functions
library(Biostrings)
library(stringr)
library(stringdist)
library(ggplot2)
library(jaccard)
library(reshape2)
library(tidyr)
library(pheatmap)
library(ggplot2)
library(parallel)

ReadFasta = function(filename){
  sv = read.table(filename)
  VBC = DNAString(sv[2,1])
  return(c("VBC" = VBC))
}

testFunction <- function (data_in) {
  return(tryCatch(data_in, error=function(e) "unknown"))
}

find_barcode<-function(x){
	tryCatch({
	usereads<-DNAString(as.character(x))
	alig<-pairwiseAlignment(BCseq,usereads,substitutionMatrix = mat,type="global-local",gapOpening = 5, gapExtension = 2)
	if(score(alig)<alig_score){
		return(NA)
	}else{
		BC<-strsplit(as.character(subject(alig)),split="")[[1]][which(strsplit(as.character(pattern(alig)),split="")[[1]]=="N")]
		BC<-paste0(BC,collapse="")
		BC<-gsub("-","",BC)
		return(BC)
	}
	},error=function(e) "unknown")	
}
		 
find_barcode_parallel<-function(dat,n){
	cl <- makeCluster(n)
	clusterEvalQ(cl,library(Biostrings))
	clusterEvalQ(cl,library(stringr))
	clusterExport(cl,c('BCseq','BClength','alig_score','mat','find_barcode','reads'))
	barcode<-parLapply(cl,dat$Read.Seq,find_barcode)
	stopCluster(cl)
	dat$Virus.BC<-unlist(barcode)
	dat<-na.omit(dat)
	dat<-dat[dat$Virus.BC!="unknown",]
	dat$VBC_len<-nchar(dat$Virus.BC)
	return(dat)
}
		 
passnreads<-function(file,reads_threshold){
	comm=paste0("sed '1d' ",file," | awk '{print$3$4$5}' | sort | uniq -c | awk '$1>=",reads_threshold," {print $2}'")
	dat=system(comm,intern = TRUE)
	return(dat)
}

passncells<-function(file,UMI_threshold){
	comm=paste0("awk 'NR>1 {print $3,$4}' ", file," | sort -u | awk '{print$1}' | uniq -c | awk '$1>",UMI_threshold," {print$2}'")
	dat=system(comm,intern = TRUE)
	return(dat)
}
		 
starcode_cluster<-function(file,starcode_path){
	dirpath=paste0(dirname(file),"/")
	system(paste0("awk 'NR>1 {print$3$5}' ",file," > ",dirpath,"CB_VB_for_cluster.tsv"))
	system(paste0(starcode_path," -d 5 --print-clusters -s ",dirpath,"CB_VB_for_cluster.tsv > ",dirpath,"CB_VB_cluster.tsv"))
	dat<-read.table(paste0(dirpath,"CB_VB_cluster.tsv"))
	return(dat)
}

split_clu<-function(x){
	cbvb<-unlist(strsplit(as.character(x[3]),split=","))
	n<-length(cbvb)
	r_cbvb=rep(as.character(x[1]),n)
	as.character(t(matrix(c(cbvb,r_cbvb),ncol=2)))
}

estimate_clu<-function(x){
	a<-substring(x[1],1,16)
	b<-substring(x[1],17,)
	c<-substring(x[2],1,16)
	d<-substring(x[2],17,)
	e<-as.numeric(stringdist(a,c))
	f<-as.numeric(stringdist(b,d))
	if(e<=4 & f<=4){
		j<-"clu"
	}else{
		j<-"non"
	}
	return(c(x[1],x[2],a,b,c,d,e,f,j))
}
r_cbvb_stat<-function(x){
	if(x[9]=="clu"){
		return(x[2])
	}else{
		return(x[1])
	}
}

final_barcode<-function(cv_clu,reads_data){
	cbvb<-apply(cv_clu,1,split_clu)
	cbvb<-data.frame(matrix(unlist(cbvb),ncol=2,byrow =T))
	j_cbvb<-data.frame(matrix(unlist(apply(cbvb,1,estimate_clu)),ncol=9,byrow=T))
	j_cbvb$rcbvb<-apply(j_cbvb,1,r_cbvb_stat)
	j_cbvb<-j_cbvb[,c(1,10)]
	names(j_cbvb)<-c("CV","rCV")
	reads_data<-merge(reads_data,j_cbvb,by="CV",all.x=T)
	newdata<-reads_data[,c(5,8)]
	newdata$Cell.BC<-substring(newdata$rCV,1,16)
	newdata$Virus.BC<-substring(newdata$rCV,17,)
	return(newdata)
}

clone_find<-function(x){
	use<-unique(as.character(newdata$Virus.BC[which(as.character(newdata$Cell.BC)==x)]))
	use<-sort(use)
	paste0(use,collapse=";")
}
		 
clone_call<-function(x){
	dat=data.frame(Cell.BC=unique(as.character(x$Cell.BC)),Virus.BC=unlist(lapply(unique(as.character(x$Cell.BC)),clone_find)))
	dat$num<-unlist(lapply(as.character(dat$Virus.BC),function(i){length(strsplit(i,split=";")[[1]])}))
	return(dat)
}

barcode_split_stat<-function(x){
  return(unlist(strsplit(x,";")))
}

jaccard_stat<-function(x,clone_tab,barcode_split){	 
  ind<-which(clone_tab$Var1==x)
  jac_dist<-unlist(lapply(barcode_split[(ind+1):length(barcode_split)],function(i){length(intersect(barcode_split[[ind]],i))/length(union(barcode_split[[ind]],i))}))
  jac_clonenum<-clone_tab$num[(ind+1):length(barcode_split)]
  jac_cellnum<-clone_tab$Freq[(ind+1):length(barcode_split)]
  return(data.frame(cell_num=jac_cellnum,tag_num=jac_clonenum,dist=jac_dist))
}

change_form_stat<-function(x,clone_names){
  VBC=unlist(strsplit(clone_names[x],"_"))
  return(data.frame(Virus.BC=VBC,clone=rep(x,length(VBC))))
}

all_jac_stat2<-function(l,VBC1,clone_stat_data){
  VBC2<-clone_stat_data[,l]
  VBC2[VBC2>0]<-1
  return(jaccard(VBC1,VBC2))
}
all_jac_stat1<-function(i,clone_stat_data,clu){
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

relation_call<-function(clone,thresh){
	clone_tab<-data.frame(table(clone$Virus.BC))
	clone_tab$Var1<-as.character(clone_tab$Var1)
	clone_tab$num<-clone$num[match(clone_tab$Var1,clone$Virus.BC)]
	clone_tab<-clone_tab[order(clone_tab$Freq,clone_tab$num),]
	#clone$time<-unlist(lapply(clone$Cluster,function(x){strsplit(x,"_")[[1]][1]}))

	barcode_split<-lapply(as.character(clone_tab$Var1),barcode_split_stat)

	all_length<-length(barcode_split)
	#library(parallel)
	#cl <- makeCluster(4)
	#clusterEvalQ(cl,library(jaccard))
	#clusterExport(cl,c('clone_tab','barcode_split','all_length','jaccard_stat'))
	jac_matirx<-lapply(clone_tab$Var1,jaccard_stat,clone_tab=clone_tab,barcode_split=barcode_split)

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

	clone_names<-lapply(c(1:length(clone_names)),change_form_stat,clone_names=clone_names)
	clone_names<-do.call("rbind",clone_names)
	clone<-merge(clone,clone_names,by="Virus.BC")
	write.csv(clone,"res/clone.csv",row.names = F,quote=F)
	clone<-clone[clone$clone %in%  as.character(data.frame(table(clone$clone))$Var1[data.frame(table(clone$clone))$Freq>1]),]
	write.csv(clone,"res/clone_2.csv",row.names = F,quote=F)

	pdf("res/clone_size_log2.pdf",width = 4,height = 3)
	plot(density(log2(table(clone$clone))))
	dev.off()

	jac_sample<-list()
	for(t in c(1:500)){
  		clone$clone_sample<-sample(clone$clone,length(clone$clone))
  		clone_tab_sample_sub<-data.frame(acast(clone,clone_sample~Cluster,fun.aggregate = sum))
  		clu<-names(clone_tab_sample_sub)
  		jac_sample_sub<-data.frame(matrix(unlist(lapply(clu,all_jac_stat1,clone_stat_data=clone_tab_sample_sub,clu=clu)),ncol=length(clu)))
  		names(jac_sample_sub)<-clu
  		row.names(jac_sample_sub)<-clu
  		jac_sample<-c(jac_sample,list(jac_sample_sub))
	}
	saveRDS(jac_sample,"res/jac_sample.rds")

	clone_tab<-data.frame(acast(clone,clone~Cluster,fun.aggregate = sum))
	clu<-names(clone_tab)
	all_jac<-data.frame(matrix(unlist(lapply(clu,all_jac_stat1,clone_stat_data=clone_tab,clu=clu)),ncol=length(clu)))
	names(all_jac)<-clu
	row.names(all_jac)<-clu
	write.csv(all_jac,"res/all_jac.csv",quote=F)

	all_p<-lapply(clu,barplot_stat,clu=clu,jac_sample=jac_sample,all_jac=all_jac)
	pdf("res/similarity_box.pdf",width = 3.5,height = 2)
		for(p in all_p){
  			print(p)
		}
	dev.off()

	pdf("res/clust.pdf",width = 4,height = 4)
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

	write.csv(all_ob,"res/all_obpre.csv",quote=F)
	write.csv(all_p,"res/all_pvalue.csv",quote=F)

	pheatmap(all_ob,cluster_rows = F,cluster_cols = F,color = rev(c(colorRampPalette(colors = c("#D73027","#FDAE61"))(200),colorRampPalette(colors =c("#FDAE61","white"))(30),colorRampPalette(colors = c("white","#4575B4"))(20))),file="res/obspre.pdf",width = 3.3,height = 3)
	dev.new()
	pheatmap(all_p,cluster_rows = F,cluster_cols = F,color = c(colorRampPalette(colors = c("#D73027","#FDAE61"))(20),colorRampPalette(colors =c("#FDAE61","white"))(30),colorRampPalette(colors = c("white","#4575B4"))(200)),file="res/pvalue.pdf",width = 3.3,height = 3)
	dev.new()
	pheatmap(all_jac,cluster_rows = F,cluster_cols = F,color = rev(c(colorRampPalette(colors = c("#D73027","#FDAE61"))(200),colorRampPalette(colors =c("#FDAE61","white"))(30),colorRampPalette(colors = c("white","#4575B4"))(20))),file="res/jaccard.pdf",width = 3.3,height = 3)
	dev.new()

	all_jac$Cluster<-row.names(all_jac)
	all_jac<-gather(all_jac,Cluster2,jaccard,-Cluster)
	all_jac<-unique(all_jac)

	all_p$Cluster<-row.names(all_p)
	all_p<-gather(all_p,Cluster2,pvalue,-Cluster)
	all_p<-unique(all_p)
	all_p$pvalue[is.na(all_p$pvalue)]<-0

	dat<-merge(all_jac,all_p,by=c("Cluster","Cluster2"))

	pdf("res/pop_all.pdf",width = 5.5,height = 4)
		ggplot(dat,aes(Cluster,Cluster2))+
  		geom_point(aes(size=jaccard,color=-log2(pvalue+0.0001)))+
  		scale_colour_gradient(low="white",high="brown4")+
  		labs(color=expression(-log2(pvalue+0.0001)),size="jaccard similarity",x="",y="",title="")+
  		scale_size(limits = c(0.0001,max(dat$jaccard)))+
  		theme_bw()+theme(axis.text.x = element_text(angle = 90,vjust = 0.5))
	dev.off()

	dat$jaccard[dat$pvalue>0.05]<-0
	dat$pvalue[dat$pvalue>0.05]<-1

	pdf("res/pop_sig.pdf",width = 5.5,height = 4)
	ggplot(dat,aes(Cluster,Cluster2))+
  		geom_point(aes(size=jaccard,color=-log2(pvalue+0.0001)))+
  		scale_colour_gradient(low="white",high="brown4")+
  		labs(color=expression(-log2(pvalue+0.0001)),size="jaccard similarity",x="",y="",title="")+
  		scale_size(limits = c(0.0001,max(dat$jaccard)))+
  		theme_bw()+theme(axis.text.x = element_text(angle = 90,vjust = 0.5))
	dev.off()
}
		 
barplot_stat<-function(x,clu,jac_sample,all_jac){
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

score_cal<-function(BCseq,BClength){
	return((length(BCseq)-floor(length(BCseq)/10)-(3*floor(length(BCseq)/10))-3*BClength))
}





