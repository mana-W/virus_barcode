# find Virus barcode 
# input reads file & barcode fa & cell type file 
# reads file header is Read.ID Read.Seq Cell.BC UMI
# cell type header : Cell.BC Cluster

# need to in put virus.fa (random seq : NNNNNNNNNN...)

source("function.R")

args=commandArgs(T)

reads=args[1]
BCseq=args[2]
celltype=args[3]

mat  =  nucleotideSubstitutionMatrix(match = 1, mismatch = -3)
BCseq=ReadFasta(BCseq)
BClength=str_count(as.character(BCseq$VBC),"N")
reads<-read.table(reads,header=T,stringsAsFactors=F)
alig_score<-score_cal(BCseq,BClength)
#environment() <- .GlobalEnv

dir.create("res")

#find virus barcode and create output file
print("Finding barcode ...")
print(date())
reads<-find_barcode_parallel(dat=reads,n=16)
print(date())
pdf("res/reads_VBC_length.pdf",width=3,height=3)
plot(density(reads$VBC_len),main="Length density")
dev.off()
reads<-reads[reads$VBC_len <= (BClength+5) & reads$VBC_len >= (BClength-5),]
write.table(reads,"res/CB_UMI_barcode.tsv",row.names=F,quote=F)
print("Finished")

#filter noise & confirm true barcode
print("Filter reads ...")
print("raw")
print(c("reads :",dim(reads)[1]))
print(c("cell barcodes :",length(unique(as.character(reads$Cell.BC)))))

passreads<-passnreads(file="CB_UMI_barcode.tsv",reads_threshold = 5)
reads<-reads[paste0(reads$Cell.BC,reads$UMI,reads$Virus.BC) %in% passreads,]
print("more than 5 reads")
print(c("reads :",dim(reads)[1]))
print(c("cell barcodes :",length(unique(as.character(reads$Cell.BC)))))
write.table(reads,"res/CB_UMI_barcode_reads.tsv",row.names=F,quote=F)

passcells<-passncells(file = "CB_UMI_barcode_reads.tsv",UMI_threshold = 1)
reads<-reads[reads$Cell.BC %in% passcells,]
print("cells with more than 2 UMIs")
print(c("reads :",dim(reads)[1]))
print(c("cell barcodes :",length(unique(as.character(reads$Cell.BC)))))
write.table(reads,"res/CB_UMI_barcode_UMI.tsv",row.names=F,quote=F)
reads$CV<-paste0(reads$Cell.BC,reads$Virus.BC)

#barcode cluster (based on similarity)

cv_clu<-starcode_cluster("res/CB_UMI_barcode_UMI.tsv",system("which starcode",inter=T))
newdata<-final_barcode(cv_clu=cv_clu,reads_data = reads)
print("cluster")
print(c("reads :",dim(newdata)[1]))
print(c("cell barcodes :",length(unique(as.character(newdata$Cell.BC)))))
write.table(newdata,"res/real_CV.tsv",row.names=F,sep="\t",quote=F)

#reform clone and match celltye
celltype<-read.table(celltype,header=T,stringsAsFactors=F)
newdata<-newdata[!duplicated(newdata[,c(3,4)]),c(3,4)]
clone_raw<-clone_call(newdata)
write.table(clone_raw,"res/clone_raw.tsv",row.names=F,quote=F)
clone<-merge(clone_raw,celltype,by="Cell.BC")
print(c("final cells :",length(clone$Cell.BC)))
write.table(clone,"res/clone_final.tsv",row.names=F,quote=F)

#calculate similarity (with in one group)
if(length(args)<4){
	thresh=1
}else{
	thresh=as.numeric(args[5])
}

relation_call(clone,thresh)
































