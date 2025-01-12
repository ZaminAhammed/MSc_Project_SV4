import pandas as pd;
import numpy as np;

fptr = open("gene_names.txt","r");
gene_name = fptr.readlines();
fptr.close();

presentgenes = [];
pf = "PF3D7_";
a1 = [];
a2 = [];
start = [];
end = [];
allgenes = [];
missinggenes = [];

#Get the names of genes that are available in GitHub
for x in gene_name:
	presentgenes.append(x[0:13]);

for i in range(1,15):
	if (i < 10):
		a1.append(pf+"0"+str(i));
	else:
		a1.append(pf+str(i));

#Sort the gene that are available to different lists based on chromosome
for x in a1:
	bckup = [];
	for y in presentgenes:
		if (x in y):
			bckup.append(y);
	a2.append(bckup)

#Get the first and last gene name in each chromosome
for i in range(0,14):
	a2[i].sort();
	start.append(a2[i][0]);
	end.append(a2[i][len(a2[i])-1]);

#Find all possible gene names between the first and last gene name in each chromosome
for i in range(0,14):
	if(i<9):
		a = int(start[i][7:11]);
		b = int(end[i][7:11])+1;
		for x in range(a,b):
			allgenes.append(pf+"0"+str(x)+"00");
	else:
		a = int(start[i][6:11]);
		b = int(end[i][6:11])+1;
		for x in range(a,b):
			allgenes.append(pf+str(x)+"00");

#Find the possible genes that could be missing
for x in allgenes:
	if x not in presentgenes:
		missinggenes.append(x);

fptr = open("missing_genes.txt","w");
for x in missinggenes:
	fptr.write(x+"\n");
fptr.close();

print(len(allgenes));
print(len(presentgenes));
print(len(missinggenes));
