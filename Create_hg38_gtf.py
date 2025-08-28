from plastid import *
import numpy as np
from twobitreader import TwoBitFile
import dill
import re

# Import necessary modules
tx = list(BED_Reader("hg38_gencode_v32.bed",return_type=Transcript))
hg38 = TwoBitFile("/seq/hg38/hg38.2bit")

# Filtering criteria: keep coding transcripts localized to chromosomes
coding_tx = [i for i in tx if i.cds_start is not None and i.chrom in hg38.keys()]  # Has CDS, localized
# 108031 transcripts

#Need to remove Tx's that are on unlocalized contigs for STAR.
coding_tx = [i for i in coding_tx if re.search('^chr[0-9XY]+$',i.chrom) is not None] #100171 tx now
len(set([i.attr['ID'] for i in coding_tx])) #still 73 duplicate tx's. Usually X/Y

#Add _1... to duplicate tx_ids to make unique. 
ct_dict = {k:0 for k in set([i.attr['ID'] for i in coding_tx])}
duplicated_ids = []
for tx in coding_tx:
    tx_id = tx.attr['ID']
    if ct_dict[tx_id] != 0:
        duplicated_ids.append(tx_id)
        tx.attr['ID'] = tx_id + "_" + str(ct_dict[tx_id])
    ct_dict[tx_id] += 1

len(set([i.attr['ID'] for i in coding_tx])) #100171 Tx's. 

#Get our mappings of tx_id -> Gene_ids. 
#The _'s we appended shouldn't matter. They modify tx version, which we we strip to match up w/ gene name anyway! 
tx2gene = {}
with open("tx_to_gene_mappings.tsv","r") as fin:
    next(fin)
    for line in fin:
        (tx_id, gene_id) = line.rstrip().split('\t')
        tx_id = tx_id.split('.')[0] #Remove the version number, see RNA_Ectoderm_1_2_pairedend_working.R 
        try: 
            tx2gene[tx_id]
            print 'duplicate key: ' + tx_id
            break
        except KeyError:
            tx2gene[tx_id] = gene_id

#Add the Gene_IDs & Gene_Names to the Transcript Objects
for i in coding_tx:
    i.attr['transcript_id'] = i.attr['ID']
    i.attr['gene_id'] = tx2gene[i.attr['ID'].split('.')[0]]
    i.attr['gene_name'] = tx2gene[i.attr['ID'].split('.')[0]]

#Need all the tx's for a given gene to be on the same chromosome & Strand
gene_tx_dict = {i:list() for i in set([v for k,v in tx2gene.iteritems()])}
for i in coding_tx:
    gene_tx_dict[i.attr['gene_id']].append(i)

multi_chrom_genes = [k for k,v in gene_tx_dict.iteritems() if len(set([i.chrom for i in v])) > 1] #18
multi_strand_genes = [k for k,v in gene_tx_dict.iteritems() if len(set([i.strand for i in v])) > 1] #1
[set([i.chrom for i in gene_tx_dict[g]]) for g in multi_chrom_genes] #all cases on diff chromos with annotations on X and Y chroms
list(set(multi_strand_genes) & set(multi_chrom_genes)) #no overlap between two

#Only keep the tx on X chromosome in multi_chrom cases. 
#Multi-Strans is only one gene (TMSB15B) so just keep +ve strand tx's. -ves are actally for TMSB15A
final_coding_tx = []
tossed_tx = []
for i in coding_tx:
    if i.attr['gene_id'] in multi_chrom_genes:
        if i.chrom == 'chrX':
            final_coding_tx.append(i)
        else:
            tossed_tx.append(i)
    elif i.attr['gene_id'] in multi_strand_genes:
        print 'Kept only +ve Txs for:' + i.attr['gene_id']
        if i.strand == '+':
            final_coding_tx.append(i)
    else:
        final_coding_tx.append(i)

len(final_coding_tx) #100096
#tossed tx turn out to be all the duplicate tx_ids. Makes sense.

#Write the gtf file
#with open("/media/storageA/hani/alignment/hg38_annotations/hg38_gencode_v32_codingonly.gtf","w") as fout:
#    for i in final_coding_tx:
#        fout.write(i.as_gtf())

#Write the Coding .bed file
#with open("/media/storageA/hani/alignment/hg38_annotations/hg38_gencode_v32_codingonly.bed","w") as fout:
#    for i in final_coding_tx:
#        fout.write(i.as_bed())

#Write the .Obj File
#with open("/media/storageA/hani/alignment/hg38_annotations/hg38_gencode_v32_codingonly.obj", "wb") as fout:
#    dill.dump(final_coding_tx, fout)

#Show the GTF2 can be read in correctly
tx_codingonly = [i for i in GTF2_TranscriptAssembler(open("/media/storageA/hani/alignment/hg38_annotations/hg38_gencode_v32_codingonly.gtf"))]
#No Tx rejected now.