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

clone_find<-function(x){
	use<-unique(as.character(newdata$Virus.BC[which(as.character(newdata$Cell.BC)==x)]))
	use<-sort(use)
	paste0(use,collapse=";")
}

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








