# Time-ordered Gene Coexpression Network (TO-GCN) for single time-series transcriptomic data
Pipeline of time-ordered gene coexpression network (TO-GCN) construction from single time-series transcriptomic data

The pipeline contains three steps: (1) Determine the PCC cutoff for TO-GCN, (2) Generate a list of TF genes in TO-GCN level with initial TF seeds, and (2) Generate a list of genes in each TO-GCN level.

## Prepare the gene expression data

Before going to the pipeline, we need to prepare two matrix of genes and expression levels (TF genes and non-TF genes) at different time points. In addition to the data files, you also need to prepare a list of initial TF genes (at lease one).

In the example folder (example_data), there are two data files from the study of "Decoding the differentiation of mesenchymal stem cells into mesangial cells at the transcriptomic level (https://doi.org/10.1186/s12864-020-06868-5)". The data file should be a Tab-separated values (.tsv) format that contains m rows and n columns, where m is the number of genes (TF genes or non-TF genes) and n represents the number of time points.The gene IDs are listed in the first column. For each gene, the expression levels of each time point are listed from the second to (n+1)-th columns. In the example data of TF_gene_matrix.tsv, there are 1086 rows for 1086 TF genes and 5 columns for 5 time points.

## Run the programs of pipeline

As mentioned above, there are three steps for the pipeline. Therefore, we provided a program for each step: (1) Cutoff, (2) TO-GCN, and (3) GeneLevel. You can directly run the program by downloading the corresponding precompiled files for different system platforms, Linux, MacOSX, or Windows. You can also download the C++ source code (.cpp) and compile to an executable one by yourself. For compiling source codes by yourself, you can use the following commands:
```sh
g++ Cutoff_STAR.cpp -o Cutoff
g++ TO-GCN_STAR.cpp -o TO-GCN
g++ GeneLevel_STAR.cpp -o GeneLevel
```
### (1) Cutoff: Determine the PCC cutoff for TO-GCN

First of all, you need postive and negative cutoff values of Pearson’s Correlation Coefficients (PCCs) for constructing the GCN. Our method is to calculate all the PCC values for each TF-gene pair. With all the PCC values, we generate distributions of probability density function (PDF) and cumulative density function (CDF). According to the CDF, we can suggest you the cutoff values with p < 0.05. To run the Cutoff program, you have to give 2 parameters: number of time points and matrix of TF genes. Here is the example of our study:
```sh
Cutoff 5 example_data/TF_gene_matrix.tsv
```
In addition to the suggested cutoff values, the program will also generate a file of PCC value distribution in a file named PCC_histogram.tsv. You can use the file to generate a histogram bar chart by Microsoft Excel or R program. 

### (2) TO-GCN: Generate a list of TF genes in TO-GCN level with initial TF seeds

The second step is to determine the time-order (level) of TF genes in the GCN. The time-order is assigned by the breadth-first search (BFS) algorithm, starting with a set of seed nodes you chose (listed in file of initial_seed.txt). In most case, we will select some genes as seeds that highly expressed in the first time point and lowly expressed in the following time points. In our study, we select a gene with ID, ENSG00000175592, and run the TO-GCN program to assign the time-order (level) in GCN. Besides, we also need the cutoff recommended in the 1st step. Two output files Node_relation.csv and Node_level.tsv can be both imported into the [Cytoscape](http://www.cytoscape.org) to visualize the network. 
```sh
TO-GCN 5 example_data/TF_gene_matrix.tsv example_data/initial_seed.txt 0.9
```
### (2) GeneLevel: Generate a list of genes in each TO-GCN level

In the final step, we want to list a set of genes (TF and non-TF genes) for each level in TO-GCN. So we need 5 parameters, number of time points, matrix of TF genes, matrix of non-TF genes, TF genes with assigned levels (output from 2nd step), and the recommended cutoff. This program output the gene list for each level in two formats. You can run GO or pathway enrichment test for these gene sets to get the dynamic biological processes during different time period. Here is an example:
```sh
GeneLevel 5 example_data/TF_gene_matrix.tsv example_data/Non-TF_gene_matrix.tsv Node_level.tsv 0.9
```
