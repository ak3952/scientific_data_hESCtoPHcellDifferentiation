#Set things up
library(tidyverse)
library(reshape2)
library(patchwork)

###############################################################
###               Mapping TxIDs -> Gene Names               ###
###############################################################
#Query Ensembl
lookup_table <- read_tsv('/media/storageA/hani/Apo_LUTI_Project/Endoderm_1_initial/SalmonOut_2/quant.sf') %>%
  dplyr::select(Name) %>%
  dplyr::rename('tx_id' = Name) %>%
  mutate('tx_id_noversion' = str_extract(tx_id, '.*(?=\\.)'))

#This changes between ensembl versions, force v101 which was used for all annotations
ensembl <- biomaRt::useEnsembl(biomart = "ensembl",
                               dataset = "hsapiens_gene_ensembl",
                               version = "101")

ens_downloaded <- biomaRt::getBM(attributes = c('ensembl_transcript_id','external_gene_name'),
                                 filters = 'ensembl_transcript_id',
                                 values = lookup_table$tx_id_noversion,
                                 mart = ensembl)

lookup_table <- lookup_table %>% 
  left_join(ens_downloaded, by = c('tx_id_noversion' = 'ensembl_transcript_id')) %>%
  mutate('external_gene_name' = ifelse(is.na(external_gene_name), tx_id, external_gene_name)) %>% #There are a very few transcripts (~500) that don't match even with this pretty degenerate search. Treat those as their own genes, not so useful for DTU anyways
  dplyr::select(tx_id, external_gene_name)

#write_tsv(lookup_table, path = '/media/storageA/hani/Apo_LUTI_Project/hg38_annotations/tx_to_gene_mappings.tsv')

