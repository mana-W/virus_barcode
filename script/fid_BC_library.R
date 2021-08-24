# find Virus barcode 
# input reads file and barcode seq

args=commandArgs(T)
reads=args[1]
BCseq=args[3]

source("/picb/sysgenomics2/projects/wangluyue/script/virus/function.R")

mat  =  nucleotideSubstitutionMatrix(match = 1, mismatch = -3)
BCseq=ReadFasta(BCseq)
BClength=str_count(as.character(BCseq$VBC),"N")
reads<-read.table(reads,header=T,stringsAsFactors=F)
alig_score<-(length(BCseq)-floor(length(BCseq)/10)-(3*floor(length(BCseq)/10))-3*BClength)
#environment() <- .GlobalEnv

#find virus barcode and create output file
print("Finding barcode ...")
library(parallel)
cl <- makeCluster(16)
clusterEvalQ(cl,library(Biostrings))
clusterEvalQ(cl,library(stringr))
clusterExport(cl,c('BCseq','BClength','alig_score','mat','find_barcode','reads'))
barcode<-parLapply(cl,reads$Read.Seq,find_barcode)
stopCluster(cl)
reads$Virus.BC<-unlist(barcode)
