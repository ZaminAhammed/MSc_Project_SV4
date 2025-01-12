import pandas as pd;
import numpy as np;

fptr = open("miss_found.txt","r");
miss_found = fptr.readlines();
fptr.close();

gene = [];
hap_num = [];
above_25 = [];
i=1;

for x in miss_found:
	gene_name = x[0:13];
	infile = str("pf-haploatlas-"+gene_name+"_population_summary.csv");
	data = pd.read_csv(infile);
	gene.append(gene_name);
	hap_num.append(len(data[data.Total > 0]));
	above_25.append(len(data[data.Total > 25]));
	print(i,".",gene_name," - done");
	i=i+1;

finaldata = pd.DataFrame({'Gene':pd.Series(gene), 'Number_of_Haplotypes':pd.Series(hap_num), 'Number_of_Haplotypes_Above_25_Samples':pd.Series(above_25)});
finaldata.to_csv("miss_pf_haploatlas.csv",index=False);
