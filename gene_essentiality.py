import pandas as pd;
import numpy as np;

gene_data = pd.read_excel("NIHMS1004827-supplement-Table_S5.xlsx", header=1);
gene_data = gene_data[['Chr','Gene_ID','MIS','MFS']];
gene_data.rename(columns = {'Gene_ID':'Gene'}, inplace=True);
gene_data.to_csv("gene_essentiality.csv",index=False);
