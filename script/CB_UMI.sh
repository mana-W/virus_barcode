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
                  
mkdir UMI_CB_umitools
zcat ${path}/R2_extracted.fastq.gz | sed -n '1~4p' | awk -F "_" '{print $1,$2,$3}' | awk '{print$1,$2,$3}' > ${path}/UMI_CB_umitools/ID_CB_UMI.tsv
zcat ${path}/R2_extracted.fastq.gz | sed -n '2~4p' > ${path}/UMI_CB_umitools/Reads.tsv
paste ${path}/UMI_CB_umitools/ID_CB_UMI.tsv ${path}/UMI_CB_umitools/Reads.tsv | awk -v OFS="\t" '{print$1,$4,$2,$3}' > ${path}/UMI_CB_umitools/cb.umi.tsv
awk '$3!="" && $4!=""' ${path}/UMI_CB_umitools/cb.umi.tsv > ${path}/UMI_CB_umitools/CB_UMI.tsv
sed -i '1i\Read.ID\tRead.Seq\tCell.BC\tUMI' ${path}/UMI_CB_umitools/CB_UMI.tsv
