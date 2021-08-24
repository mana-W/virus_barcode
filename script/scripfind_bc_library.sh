#input is R1 R2

#find barcode for library

#first step assamble
#falsh 
mkdir fq
flash $1 $2 -d fq

mkdir reads_BC
#after flash extract reads
cat fq/out.extendedFrags.fastq | sed -n '2~4p' > reads_BC/Reads
awk -F 'GGCTGGATCC||ACCGGTTTGC' '{print$2}' reads_BC/Reads > reads_BC/barcode
paste reads_BC/Reads reads_BC/barcode > reads_barcode.tsv
