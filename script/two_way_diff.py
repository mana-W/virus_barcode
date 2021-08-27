import os
import numpy as np
from collections import Counter
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib_venn import venn2

os.chdir("/Users/wly/data/chenlab/library_lenti/") 
data1 = pd.read_csv('reads_CB.csv')
data2=pd.read_csv("reads_barcode.tsv",sep="\t",header=None)
data2.columns=["Read.Seq","Virus.BC"]
VB_1=list(data1['Virus.BC'])
VB_2=list(data2['Virus.BC'])
data={'VB_1':VB_1,'VB_2':VB_2}
data=pd.DataFrame(data)


#R script result
d=Counter(data['VB_1']) 
#data[ column_1 ].value_counts()
d=pd.DataFrame(sorted(d.items(),key=lambda x:x[1],reverse=True))
d.columns=["Virus.BC","num"]
d=d[d['num']>1].dropna()
d['logN']=np.log2(d['num'])
 #num distribution
d.logN.plot.density(color='green')
plt.title("Density plot for Virus barcode (log2)")
plt.savefig('density_num.pdf')
plt.show()
 #length distribution
d1=pd.DataFrame.from_dict({'Length':data['VB_1'].str.len()})
d2=pd.DataFrame.from_dict({'Length':data['VB_2'].str.len()})
d1.Length.plot.density(color='green')
plt.title("Barcode length")
plt.savefig('density_length_R.pdf')
plt.show()
d2.Length.plot.density(color='green')
plt.title("Barcode length")
plt.savefig('density_length.pdf')
plt.show()

#two data state
#pie
i=data[data["VB_1"]==data["VB_2"]].shape[0]
a=data.shape[0]
b=data.shape[0]
venn2(subsets=(a-i,b-i,i),set_labels=("align","grep"))
plt.title("Barcode consistency")
plt.savefig('consistency.pdf')
plt.show()

#reshape to motif form
#base & proportion at each site
#select barcode in 32nt (VB_1)
VB_1=data1['Virus.BC'][data1['Virus.BC'].str.len()==32]
VB_1.index=range(VB_1.shape[0])
df={}
for i in range(len(VB_1[1])):
    df["P"+str(i)]={"A":0,"T":0,"C":0,"G":0}


for i in range(VB_1.shape[0]):
    bc=list(VB_1[i])
    for n in range(len(list(VB_1[i]))):
        rawnum=df["P"+str(n)][list(VB_1[i])[n]]
        df["P"+str(n)][list(VB_1[i])[n]]=rawnum+1


df=pd.DataFrame.from_dict(df)
df=df.div(df.sum(axis=0),axis=1)
df.to_csv("motif_data.csv")






