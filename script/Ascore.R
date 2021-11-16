ascore_fun<-function(clone,data,type1,type2){
  alltype<-data$Type[data$clone==clone]
  ascore<-length(which(alltype==type1))/(length(which(alltype==type2))+0.5)
  return(ascore)
}

sample_ascore_fun<-function(type1,type2,data){
  dat<-data[sample(which(data$Type==type1),round(length(which(data$Type==type1))/10)),]
  ascore<-unlist(lapply(unique(dat$clone),ascore_fun,data=data,type1=type1,type2=type2))
  return(ascore)
}

library(ggplot2)
args<-commandArgs(T)
data<-read.csv(args[1],header=T)
Prog_1000<-unlist(lapply(rep("Prog",1000),sample_ascore_fun,type2="Neu",data=data))
Neu_1000<-unlist(lapply(rep("Neu",1000),sample_ascore_fun,type2="Prog",data=data))
Prog_1000<-data.frame(Ascore=Prog_1000,Type="Prog")
Neu_1000<-data.frame(Ascore=Neu_1000,Type="Neu")
d_1000<-rbind(Prog_1000,Neu_1000)

pdf("Ascore.pdf",width = 4.5,height = 4)
ggplot(d_1000,aes(color=Type))+
  geom_density(aes(x=Ascore),adjust=2)+
  theme_bw()
dev.off()

#quantile(Prog_1000$Ascore,0.6)
#quantile(Neu_1000$Ascore,0.6)

clone<-data.frame(clone=unique(data$clone))
clone$Ascore_Prog<-unlist(lapply(clone$clone,ascore_fun,data=data,type1="Prog",type2="Neu"))
clone$Ascore_Neu<-unlist(lapply(clone$clone,ascore_fun,data=data,type1="Neu",type2="Prog"))
pdf("Ascore_compare.pdf",width = 4,height = 4)
ggplot(clone,aes(x=Ascore_Prog,y=Ascore_Neu))+geom_point()+theme_bw()
dev.off()
clone$Type_Prog<-"F"
clone$Type_Neu<-"F"
clone$Type_Prog[clone$Ascore_Prog>=quantile(Prog_1000$Ascore,0.65)]<-"T"
clone$Type_Neu[clone$Ascore_Neu>=quantile(Neu_1000$Ascore,0.65)]<-"T"
if(dim(clone[clone$Type_Prog=="T" & clone$Type_Neu=="T",])[1]==0){
  clone$Type<-"None"
  clone$Type[clone$Type_Prog=="T"]<-"Prog"
  clone$Type[clone$Type_Neu=="T"]<-"Neu"
  pdf("Ascore_type.pdf",width = 4.8,height = 4)
  ggplot(clone,aes(x=Ascore_Prog,y=Ascore_Neu,color=Type))+
    geom_point(alpha=0.4,size=3)+
    theme_bw()+
    scale_color_manual(values=c("#ff9de2","gray77","#8c82fc"))
  dev.off()
  write.csv(clone,"clone_type.csv",row.names = F,quote=F)
}else{
  write.csv(clone,"clone_type.csv",row.names = F,quote=F)
  print("Warnning: Clone type undefined！")
}

