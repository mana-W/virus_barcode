
**Barlin**
===========

If your data include UMI、cell barcode and exogenous virus barcode <br />
you can use Barlin to extract all of this tags and calculate the intergroup similarity.



## Installation

Barlin can only work on typical Linux systems.

```
git clone https://github.com/mana-W/virus_barcode.git
```

## Dependencies


starcode : https://github.com/gui11aume/starcode <br />
umitools : https://github.com/CGATOxford/UMI-tools

**R**
```
pcks <- c("stringr","stringdist","ggplot2","jaccard","reshape2","tidyr","pheatmap","parallel","ggalluvial")
install.packages(pcks)
```

## Usage


### Data prepare
1. Fastq file (with virus barcode):<br />
R1.fastq.gz<br />
R2.fastq.gz<br />
2. Fasta file of barcode template sequence:<br />
virus.fa<br />
3. Cells annotation (tsv):<br />
celltype.tsv<br />
Contents in column 'Cluster' have to be like: group_annotation
<img src="https://github.com/mana-W/chenlab_you/blob/main/images/celltype.PNG" width = "200" height = "230" align="center" />

<br />

### Running<br />
**Step1**: Extract cell barcode and UMIs, prepare input file for next step.
```
sh CB_UMI.sh R1.fastq.gz R2.fastq.gz
```

**Step2**: Recover virus barcodes of cells and relationship between each pair of clusters.
```
Rscript find_virusBC.R UMI_CB_umitools/CB_UMI.tsv virus.fa celltype.tsv 0.5
```
Output of this step in directory *res*.<br />
The most important file is *res/clone_final.tsv*, include imformation of barcodes.<br />


**Step3 (optional)**: Calculate cells' relationship span multiple groups and results visualization.
```
Rscript similarity.R 0.5 group1/res/clone_final.tsv group2/res/clone_final.tsv
```
Output of this step in directory *spanres*.<br />
<br />
After this step you can also create a sanky plot by:<br />
```
Rscript alluvial_plot.R group1_celltype.tsv group2_celltype.tsv spanres/all_jac.csv spanres/all_pvalue.csv
```
<img src="https://github.com/mana-W/chenlab_you/blob/main/images/sanky.PNG" width = "400" height = "400" align="center" />



