import pandas as pd;
import numpy as np;

fptr1 = open("gene_names.txt", "r");
genes = list(fptr1.readlines());
fptr1.close()

gene = [];
hap_num = [];
above_25 = [];
i = 1;

for x in genes:
	gene_name = x[0:13];
	infile = str(gene_name+".pkl");
	pckl = pd.read_pickle(infile);
	data = pckl[0];
	gene.append(gene_name);
	hap_num.append(len(data[data.Total > 0]));
	above_25.append(len(data[data.Total > 25]));
	print(i,".",gene_name," - done");
	i=i+1;
	
finaldata = pd.DataFrame({'Gene':pd.Series(gene), 'Number_of_Haplotypes':pd.Series(hap_num), 'Number_of_Haplotypes_Above_25_Samples':pd.Series(above_25)});
finaldata.to_csv("pf_haploatlas.csv",index=False);
