#umitools extract cell barcode & UMI

R1=`echo $1`
path=`echo ${R1%/*}`
R2=`echo $2`

umi_tools whitelist --stdin ${R1} --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNNNN --set-cell-number=10000 --log2stderr > ${path}/whitelist.txt

umi_tools extract --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNNNN \
                  --stdin ${R1} \
                  --stdout ${path}/R1_extracted.fastq.gz \
                  --read2-in ${R2} \
                  --read2-out=${path}/R2_extracted.fastq.gz \
                  --filter-cell-barcode \
                  --whitelist=${path}/whitelist.txt
