This repository contains code associated with the following study:

**Temporal multi-omics gene expression data across human embryonic stem cell-derived polyhormonal cell differentiation**

Abdurrahman Keskin<sup>1</sup>, Hani J. Shayya<sup>2</sup>, Achchhe Patel<sup>3</sup>, Dario Sirabella<sup>3</sup>, Barbara Corneo<sup>3</sup> & Marko Jovanovic<sup>1</sup>

<sup>1</sup>Department of Biological Sciences, Columbia University, New York, NY 10027, USA 

<sup>2</sup>Mortimer B. Zuckerman Mind, Brain and Behavior Institute, Columbia University, New York, NY 10027, USA 

<sup>3</sup>Columbia Stem Cell Initiative, Stem Cell Core, Columbia University Irving Medical Center, New York, NY 10032, USA

**Data Availability:**

RNA-seq: NCBI GEO: GSE305933

Ribo-seq: NCBI GEO: GSE305934

**Contents:**

***Create_hg38_gtf.py***: Creates custom GTF/STAR indices used in the study  

***Transcript_to_gene_mapping.R***: Maps transcript IDs to gene names using Ensembl database

***Read_Processing_and_Alignment.sh***: Aligns RNA-seq and Ribo-seq reads to the genome using STAR

***custom_featurecounts.py***: Performs read counting of the RNA-seq and Ribo-seq data

**Note: All scripts were written by Hani J. Shayya**
