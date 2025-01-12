import pandas as pd;
import numpy as np;

pf = pd.read_csv("pf_haploatlas.csv");
ge = pd.read_csv("gene_essentiality.csv");

final = pd.merge(ge,pf,on='Gene',how="outer");
final.to_csv("integrated.csv",index=False);
