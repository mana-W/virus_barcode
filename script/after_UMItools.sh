#input is R2

R2=`echo $1`
path=`echo ${R2%/*}`

mkdir UMI_CB_umitools
zcat ${R2} | sed -n '1~4p' | awk -F "_" '{print $1,$2,$3}' | awk '{print$1,$2,$3}' > ${path}/UMI_CB_umitools/ID_CB_UMI
zcat ${R2} | sed -n '2~4p' > ${path}/UMI_CB_umitools/Reads
paste ${path}/UMI_CB_umitools/ID_CB_UMI ${path}/UMI_CB_umitools/Reads | awk -v OFS="\t" '{print$1,$4,$2,$3}' > ${path}/UMI_CB_umitools/cb.umi.tsv
awk '$3!="" && $4!=""' ${path}/UMI_CB_umitools/cb.umi.tsv > ${path}/UMI_CB_umitools/CB_UMI
sed -i '1i\Read.ID\tRead.Seq\tCell.BC\tUMI' ${path}/UMI_CB_umitools/CB_UMI
awk -F 'GCAAACCGGT||GGATCCAGCC' '{print$2}' ${path}/UMI_CB_umitools/CB_UMI > ${path}/UMI_CB_umitools/barcode
paste -d '\t' ${path}/UMI_CB_umitools/CB_UMI ${path}/UMI_CB_umitools/barcode > ${path}/UMI_CB_umitools/CB_UMI_barcode.txt
awk 'NF==5' ${path}/UMI_CB_umitools/CB_UMI_barcode.txt | sed '1d' > ${path}/UMI_CB_umitools/CB_UMI_barcode_full.txt
Rscript /picb/sysgenomics2/projects/wangluyue/You/you_en1_lmx1a/pipeline/real_CV.R ${path}/UMI_CB_umitools/CB_UMI_barcode_full.txt
/picb/bigdata/project/wly/tools/starcode/starcode -d 5 --print-clusters -s ${path}/UMI_CB_umitools/CB_VB_for_cluster.txt > ${path}/UMI_CB_umitools/CB_VB_cluster.txt
Rscript /picb/sysgenomics2/projects/wangluyue/You/you_en1_lmx1a/pipeline/real_CV2.R ${path}/UMI_CB_umitools
#cells
awk '{print$3}' ${path}/UMI_CB_umitools/real_CV.txt | sort -u > ${path}/UMI_CB_umitools/cells 
Rscript /picb/sysgenomics2/projects/wangluyue/You/you_en1_lmx1a/pipeline/real_CV3.R ${path}/UMI_CB_umitools/real_CV.txt ${path}/UMI_CB_umitools/cells





#Rscript /picb/sysgenomics2/projects/wangluyue/Xie/Dual_2nd/1.find_scar_Ver3.R ${path}/UMI_CB_umitools
#Rscript /picb/sysgenomics2/projects/wangluyue/Xie/Dual/pipeline/1.find_scar_V2.R ${outpath}
#Rscript /picb/sysgenomics2/projects/wangluyue/Xie/Dual_2nd/2.consensus_plot_v1_ver3.R ${path}/UMI_CB_umitools
#Rscript /picb/sysgenomics2/projects/wangluyue/Xie/Dual_2nd/2.consensus_plot_v2_ver3.R ${path}/UMI_CB_umitools

##########s
#Rscript /picb/sysgenomics2/projects/wangluyue/Xie/2.consensus_plot_v1_ver2.R ${path}/v1/UMI_reads_scar_full.txt "GATACGATACGCGCACGCTATGGTTTAGACACGACTCGCGCATACGATGGTTAAGATAGTATGCGTATACGCTATGGTCAAGATATGCATAGCGCATGCTATGGTTATGAGTCGAGACGCTGACGATATGGTCTTGATATGAGACTCGCATGTGATGGTCTAGCGACTGTACGCACACGCGATGGTCATGATACGTAGCACGCAGACTATGGTCAAGACACAGTACTCTCACTCTA" ${outpath}/v1
#Rscript /picb/sysgenomics2/projects/wangluyue/Xie/Dual/pipeline/2.consensus_plot_v2_ver2.R ${path}/v2/UMI_reads_scar_full.txt "GAGTCGAGACGCTGACGATATGGTTTAGATATGAGACTCGCATGTGATGGTTAAGATACGATACGCGCACGCTATGGTCAAGACACGACTCGCGCATACGATGGTTATGATACGTAGCACGCAGACTATGGTCTTGACACAGTACTCTCACTCTATGGTCTAGCGACTGTACGCACACGCGATGGTCATGATAGTATGCGTATACGCTA" ${outpath}/v2


