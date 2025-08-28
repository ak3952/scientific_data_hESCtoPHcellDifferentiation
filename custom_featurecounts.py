#!/usr/bin/env python

#Import Relevent Packages
from __future__ import division
from plastid import GTF2_TranscriptAssembler, BAMGenomeArray, VariableFivePrimeMapFactory, FivePrimeMapFactory, CenterMapFactory
import twobitreader as twobit
import numpy as np
import re 
from collections import Counter, defaultdict
import itertools
from plastid.util.io.filters import CommentReader

#Import Additional Relevent Packages
import dill
import pandas as pd
from pathos.multiprocessing import ProcessingPool as Pool 
import copy
import os 
import argparse

################################
######  Define Functions  ######
################################

def meta_roi(segmentchain_list):
    """An old function that I wrote and have adopted here. Takes a list of segment chains, breaks them apart into constituent genomic segments, and throws them back together into a segment chain representing a meta-window. 
    
    Function was updated on 12/14/18 to deepcopy starter segment for construcitng the meta window. This ensures that input segmentchain list is not altered by the function.
    
    --Input--
    segmentchain_list: a list of plastid segmentchains or plastid transcripts. 
    
    --Output--
    A single plastid segmentchain representing the union of all input segmentchains. Attributes are propagated from the first segmentchain in the input list. 
    """
    if type(segmentchain_list) is str:
        return(segmentchain_list) 
    
    elif (len(segmentchain_list) == 1):
        uORF_meta = segmentchain_list[0]
        return uORF_meta
    
    else:
        # Figure out how to chain the genomic segments together
        segments_list = []
        for i in segmentchain_list[1:len(segmentchain_list)]:
            segments_list.append(i.segments)
        
        segments_flattened = list(itertools.chain.from_iterable(segments_list))
        uORF_meta = copy.deepcopy(segmentchain_list[0]) #implemented deepcopy approach here which is necessary if you want to call this again.
        uORF_meta.add_segments(*segments_flattened) 
        return uORF_meta

def get_counts_and_lengths_masked(input_segment_chain, mapped_read_array, masked_logical, keep_true):
    """A basic method for getting counts and lengths from a plastid segmentchain and a plastid BAMGenome Array. Advantage over get_counts is that it provides methods for handling masked regions and functionality for keeping either masked OR unmasked region read/length counts.
    
    --Input--
    input_segment_chain: plastid segment_chain_object to calculate parameters over
    mapped_read_array: plastid BAMGenomeArray with mapping set, containing reads that need to be counted
    masked_logical: (True|False) do you want to look at masks that exist on the input segment_chain? Masks must have been previously added with .add_masks method
    keep_true: (string: 'yes'|'no') applicable only when masked_logical=True. Do you want to count reads/length for masked region (yes) or the NON-masked region (no)
    
    --Output--
    a tuple: (counts, length) for input_segment_chain given mapped_read_array
    """
    if masked_logical:
        masked_counts = input_segment_chain.get_masked_counts(mapped_read_array)
        if keep_true == 'yes':
            roi_masked_counts = input_segment_chain.get_counts(mapped_read_array)[masked_counts.mask] #keep only reads in masked region
            my_counts = np.nansum(roi_masked_counts)
            my_length = len(roi_masked_counts)
        #
        else:
            roi_masked_counts = input_segment_chain.get_counts(mapped_read_array)[np.invert(masked_counts.mask)] #Trick b/c get_masked_counts doesnt count the masked regions
            my_counts = np.nansum(roi_masked_counts)
            my_length = len(roi_masked_counts)
        #
    #
    else:
        my_counts = np.nansum(input_segment_chain.get_counts(mapped_read_array))
        my_length = input_segment_chain.get_length()
    #
    return my_counts, my_length

def genewise_counter_riboprof(run_mode, gene_id, full_transcript_list, uORF_list, path_to_ribprof_bam, path_to_psite, path_to_rnaseq, rnamapfactory):
    """A genewise counter function for ribosome profiling and RNA-seq data. When run in rnaseq mode, returns cds and full transcript counts as well as lengths. When run in riboprof mode, returns ribosome profiling counts and lengths for various regions of the meta-transcript (UTR/CDS/uORFs) under various sets of assumptions as detailed below. When run in 'both' mode, returns both ribosome profiling and RNAseq read counts, along with a translational efficiency that is calculated internally. 
    
    --Input--
    run_mode: str('riboprof' | 'rnaseq' | 'both') run mode to use for counting. See output below. 
    gene_ID: string representing gene ID to be quantified
    full_transcript_list: list of complete plastid transcript objects to be concatenated for the given gene_id. Ensure these are complete objects, the UTR/CDS parsing will happen internally. 
    uORF_list: list of all uORFs to be included in for the given gene_id. These should be plastid segmentchains. The meta roi will be constructed internally. If there are no uORFs, you must pass an empty list to this paramater (ie. uORF_list = [])
    path_to_ribprof_bam: string giving path to ribosome profiling CHX reads. Required if mode is "riboprof" or "both". Otherwise pass 'None'.
    path_to_psite: string giving path to file containing table of read lengths:psite offsets (generated by plastid psite script). Required if run_mode is "riboprof" or "both". Otherwise pass 'None'.
    path_to_rnaseq: string giving path to mRNA-seq reads. Required if mode is "rnaseq" or "both". Otherwise pass 'None'.
    rnamapfactory: str('fiveprime'|'center') detailing which map factory to use for counting rnaseq data. 'fiveprime' will use plastid's FivePrimeMapFactory which will assign a value of 1 to the 5'-most base of the read. 'center' will use plastid's CenterMapFactory which assigns a value of 1/L to every position on read of length L. Use 'fiveprime' for tools like DESeq2. 
    
    --Output--
    'rnaseq' mode: a list of length 6- [gene_id, gene_name, reads in meta-cds, length of meta-cds, reads in meta-transcript, length of meta-transcript]. Gene name is taken from attributes of input transcript list, meta-transcript and meta-cds regions are constructed via internal calls to meta-roi using the full transcript list and the annotated CDS from that full transcript list respectively. 
    'riboprof' mode: a list of length 12- [gene_id, gene_name, utr5_counts, utr5_len, cds_counts, cds_len, utr3_counts, utr3_len, meta_uORF_counts, meta_uORF_len, meta_uORF_counts_in5utronly_noCDS, meta_uORF_len_in5utronly_noCDS]. 
        utrs are defined as the meta-region constructed from 5utr and 3utr annotations of input transcript list. These regions are masked by the meta-CDS, to ensure that reads are not counted twice. 
        'meta_uORF' parameters refer to meta-regions constructed via a call to meta-uORF with input uORF_list. No masks are applied.
        'meta_uORF_in5utronly_noCDS' parameters are the result of masking meta-uORF region with the 5'utr, which is itself masked by the CDS. These parameters refer to the region of the meta-uORF which is in the 5'UTR and does not overlap CDS on any transcript for that gene. 
    'both' mode: a list of length 15- [gene_id, gene_name, ribosome profiling utr5 counts, utr5_length, ribosome profiling cds_counts, cds_len, ribosome profiling utr3 counts, utr3_len, ribosome profiling meta_uORF_counts, meta_uORF_len, ribosome profiling meta_uORF_counts_in5utronly_noCDS, meta_uORF_len_in5utronly_noCDS, rnaseq counts cds, rnaseq counts full meta-transcript, length full meta-transcript]
        ribosome profiling parameters are precisely as documented above for 'riboprof' mode
        rnaseq parameters are as described above for 'rnaseq' mode
        full transcript length gives the length of the full meta-transcript constructed via a call to meta-roi using the input full_transcript_list. This is useful for calculating rpkm for the full transcript. 
    """
    
    print 'Working on ' + str(gene_id)
    #
    #Test run mode. Now implemented at top level in wrapper function. 
    #if run_mode not in ['riboprof', 'rnaseq', 'both']:
    #    raise Exception("Did not understand the value for run_mode")
    #
    #Set up meta-window annotations for given gene using full transcript annotations
    CDS_list = [i.get_cds() for i in full_transcript_list]
    utr5_list = [i.get_utr5() for i in full_transcript_list]
    utr3_list = [i.get_utr3() for i in full_transcript_list]
    meta_CDS = meta_roi(CDS_list)
    meta_transcript = meta_roi(full_transcript_list)
    
    #Import the reads for rna
    if run_mode == 'rnaseq' or run_mode == 'both':
        rnaseq_reads = BAMGenomeArray(path_to_rnaseq)
        if rnamapfactory == 'fiveprime':
            rnaseq_reads.set_mapping(FivePrimeMapFactory())
        elif rnamapfactory == 'center':
            rnaseq_reads.set_mapping(CenterMapFactory())
        else: 
            raise Exception("Did not understand the value for rnamapfactory.")
        #
    #
    
    #Import the reads for ribosome profiling
    if run_mode == 'riboprof' or run_mode == 'both':
        CHX_reads = BAMGenomeArray(path_to_ribprof_bam)
        CHX_reads.set_mapping(VariableFivePrimeMapFactory.from_file(open(path_to_psite)))
    #
    
    #RNA-seq read counting
    if run_mode == 'rnaseq' or run_mode == 'both':
        cds_rnaseq_counts = np.nansum(meta_CDS.get_counts(rnaseq_reads))
        full_transcript_rnaseq_counts = np.nansum(meta_transcript.get_counts(rnaseq_reads))
        full_transcript_meta_length = meta_transcript.get_length()
        if run_mode == 'rnaseq':
            output_line = [gene_id, full_transcript_list[0].attr['gene_name'], cds_rnaseq_counts, meta_CDS.get_length(), full_transcript_rnaseq_counts, full_transcript_meta_length] 
            print '...done!'
            return output_line
        #
    #
    
    #Ribosome profiling analysis, used for run_mode = 'both' or 'riboprof'
    meta_utr5 = meta_roi(utr5_list)
    meta_utr3 = meta_roi(utr3_list)
    
    #Mask CDS overlap of 5'UTR and 3'UTR windows. For calculating 5'UTR/CDS/3'UTR ratios most efficiently
    meta_utr5.add_masks(*meta_CDS.segments)
    meta_utr3.add_masks(*meta_CDS.segments)
    
    #Import Reads and Set mapping
    CHX_reads = BAMGenomeArray(path_to_ribprof_bam)
    CHX_reads.set_mapping(VariableFivePrimeMapFactory.from_file(open(path_to_psite))) #doesn't read string like it should, so my workaround is to give it the open file handle instead
    
    #Calculate Everything Except for uORF Parameters
    utr5_counts_maskedbyCDS, utr5_len_maskedbyCDS = get_counts_and_lengths_masked(input_segment_chain=meta_utr5, mapped_read_array=CHX_reads, masked_logical=True, keep_true='no')
    cds_counts, cds_len = get_counts_and_lengths_masked(input_segment_chain=meta_CDS, mapped_read_array=CHX_reads, masked_logical=False, keep_true='no')
    utr3_counts_maskedbyCDS, utr3_len_maskedbyCDS = get_counts_and_lengths_masked(input_segment_chain=meta_utr3, mapped_read_array=CHX_reads, masked_logical=True, keep_true='no')
    
    #Build in if statements b/c may not have uORFs 
    if len(uORF_list)!= 0:
        meta_uORF = meta_roi(uORF_list)
        meta_uORF_counts, meta_uORF_len = get_counts_and_lengths_masked(input_segment_chain=meta_uORF, mapped_read_array=CHX_reads, masked_logical=False, keep_true='no') #whole meta_uORF used for uORF/CDS ratio count. 
        meta_uORF.add_masks(*meta_utr5.segments) 
        meta_uORF_in5utr = meta_uORF.get_masks_as_segmentchain() #takes the intersection of uORF and 5'UTR
        meta_uORF_in5utr.add_masks(*meta_CDS.segments) #masks the CDs
        meta_uORF_counts_in5utronly_noCDS, meta_uORF_len_in5utronly_noCDS = get_counts_and_lengths_masked(input_segment_chain=meta_uORF_in5utr, mapped_read_array=CHX_reads, masked_logical=True, keep_true='no') #used for uORF/5'UTR 
    else:
        meta_uORF_counts, meta_uORF_len = (np.nan, np.nan)
        meta_uORF_counts_in5utronly_noCDS, meta_uORF_len_in5utronly_noCDS = (np.nan, np.nan)
        #
    #
    
    #Output for 'riboprof mode'
    if run_mode == 'riboprof':
        output_line = [gene_id, full_transcript_list[0].attr['gene_name'], utr5_counts_maskedbyCDS, utr5_len_maskedbyCDS, cds_counts, cds_len, utr3_counts_maskedbyCDS, utr3_len_maskedbyCDS, meta_uORF_counts, meta_uORF_len, meta_uORF_counts_in5utronly_noCDS, meta_uORF_len_in5utronly_noCDS]
        print '...done!'
        return output_line
    #
    if run_mode == 'both':
        output_line = [gene_id, full_transcript_list[0].attr['gene_name'], utr5_counts_maskedbyCDS, utr5_len_maskedbyCDS, cds_counts, cds_len, utr3_counts_maskedbyCDS, utr3_len_maskedbyCDS, meta_uORF_counts, meta_uORF_len, meta_uORF_counts_in5utronly_noCDS, meta_uORF_len_in5utronly_noCDS, cds_rnaseq_counts, full_transcript_rnaseq_counts, full_transcript_meta_length]
        print '...done!'
        return output_line
    #

def feature_counter(run_mode, list_of_transcripts, list_of_uORFs, path_to_ribprof_bam, path_to_psite, path_to_rnaseq, rnamapfactory, nthreads):
    """Wrapper function for genewise_counter_riboprof that parses input transcripts/uORF into groups by gene, generates genewise lists that are fed into genewise_counter_riboprof in a multithreaded way, and creates a clean pandas dataframe with count paramaters for all genes analyzed. Calculates relevent parameters like rpkm for rnaseq and translational efficiency, which require all data in aggreate. 
    
    --Input--
    run_mode: str('riboprof' | 'rnaseq' | 'both') run mode to use for counting. Passed to genewise_counter_riboprof internally. See output below. 
    list_of_transcripts: list of all annotated transcripts in genome to be quantified (eg. all coding transcripts from mm10). Transcripts should be represented as plastid transcript objects, with a 'gene_id' attribute that will be used to group transcripts by gene for meta-region construction and analysis. Transcripts should have a 'gene_name' attribute which will be included in the output table (assuming that gene_id is a UCSC or Ensembl ID that is not particuarily helpful). 
    list_of_uORFs: list of all uORFs to be included in the analysis. uORFs should be plastid segmentchains, with a gene_id that can be cross-referenced with the list_of_transcripts. 
    path_to_ribprof_bam: string giving path to ribosome profiling CHX reads. Required if mode is "riboprof" or "both". Otherwise pass 'None'.
    path_to_psite: string giving path to file containing table of read lengths:psite offsets (generated by plastid psite script). Required if run_mode is "riboprof" or "both". Otherwise pass 'None'. 
    path_to_rnaseq: string giving path to mRNA-seq reads. Required if mode is "rnaseq" or "both". Otherwise pass 'None'. 
    rnamapfactory: str('fiveprime'|'center') detailing which map factory to use for counting rnaseq data. 'fiveprime' will use plastid's FivePrimeMapFactory which will assign a value of 1 to the 5'-most base of the read. 'center' will use plastid's CenterMapFactory which assigns a value of 1/L to every position on read of length L. Use 'fiveprime' for tools like DESeq2. 
    nthreads: number of threads to be used for multiprocessing. 
    
    --Output--
    'rnaseq' mode: a pandas dataframe of shape (number_genes_input, 7). Columns: [gene_id, gene_name, reads in meta-cds, length of meta-cds, reads in meta-transcript, length of meta-transcript, rpkm]. 
        Gene name is taken from attributes of input transcript list
        Meta-transcript and meta-cds regions are constructed via internal calls to meta-roi using the full transcript list and the annotated CDS from that full transcript list respectively. 
        RPKM is definied for the total meta-transcript: (reads in meta-transcript) / (length of meta-transcript in kb * total number of reads in any annotated meta-transcript) * 1,000,000. 
    'riboprof' mode: a pandas dataframe of shape (number_genes_input, 12). Columns: [gene_id, gene_name, utr5_counts, utr5_len, cds_counts, cds_len, utr3_counts, utr3_len, meta_uORF_counts, meta_uORF_len, meta_uORF_counts_in5utronly_noCDS, meta_uORF_len_in5utronly_noCDS]. 
        utrs are defined as the meta-region constructed from 5utr and 3utr annotations of input transcript list. These regions are masked by the meta-CDS, to ensure that reads are not counted twice. 
        'meta_uORF' parameters refer to meta-regions constructed via a call to meta-uORF with input uORF_list. No masks are applied.
        'meta_uORF_in5utronly_noCDS' parameters are the result of masking meta-uORF region with the 5'utr, which is itself masked by the CDS. These parameters refer to the region of the meta-uORF which is in the 5'UTR and does not overlap CDS on any transcript for that gene. 
    'both' mode: a pandas dataframe of shape (number_genes_input, 16). Columns: [gene_id, gene_name, ribosome profiling utr5 counts, utr5_length, ribosome profiling cds_counts, cds_len, ribosome profiling utr3 counts, utr3_len, ribosome profiling meta_uORF_counts, meta_uORF_len, ribosome profiling meta_uORF_counts_in5utronly_noCDS, meta_uORF_len_in5utronly_noCDS, rnaseq counts cds, rnaseq counts full meta-transcript, length full meta-transcript, translation_efficiency]
        ribosome profiling parameters are precisely as documented above for 'riboprof' mode
        rnaseq parameters are as described above for 'rnaseq' mode
        full transcript length gives the length of the full meta-transcript constructed via a call to meta-roi using the input full_transcript_list. 
        translational efficiency follows the definition of Ingolia et al, 2011: (riboprof_reads_meta_cds / [length meta_CDS in kb * total number of ribosome reads in any region of any annotated transcript]) / (rnaseq_reads_meta_cds / [length meta_CDS in kb * total number of rnaseq reads in any region of any annotated transcript])
            Note that the Ingolia definition of TE normalizes to reads in any region of any annotated transcript (UTR or CDS), rather than the total number of reads that STAR could align. This avoids confounding by rRNA or small ncRNA, which STAR can align, but may differ between experiments.  
    """
    #Test run mode.  
    if run_mode not in ['riboprof', 'rnaseq', 'both']:
        raise Exception("Did not understand the value for run_mode")
    #
    print "Creating Dictionaries..."
    #Create full transcript dictionary
    full_transcript_dict = defaultdict(list)
    for transcript in list_of_transcripts:
        full_transcript_dict[transcript.get_gene()].append(transcript)
    #
    #Create meta-uORF dictionary
    uORF_dict = defaultdict(list)
    for uORF in list_of_uORFs:
        uORF_dict[uORF.attr['gene_id']].append(uORF)
    #
    #Merged dictionary has structure gene_id: list of 2 [full transcripts list, uORFs list]
    gene_id_list = []
    full_transcript_list_of_lists = []
    uORF_list_of_lists = []
    for gene_id in full_transcript_dict.iterkeys():
        gene_id_list.append(gene_id)
        full_transcript_list_of_lists.append(full_transcript_dict[gene_id])
        try:
            uORF_list_of_lists.append(uORF_dict[gene_id])
        except KeyError:
            uORF_list_of_lists.append([])
        #
    
    #Run the loop
    print "Counting Features..."
    output_list = Pool(nthreads).map(genewise_counter_riboprof, itertools.repeat(run_mode), gene_id_list, full_transcript_list_of_lists, uORF_list_of_lists, itertools.repeat(path_to_ribprof_bam), itertools.repeat(path_to_psite), itertools.repeat(path_to_rnaseq), itertools.repeat(rnamapfactory))
    
    print "Concatenating Dataframes..."
    if run_mode == 'rnaseq':
        output_df = pd.DataFrame(output_list, columns = ["gene_id", "gene_name", "rnaseq_cds_counts", "cds_len", "rnaseq_total_counts", "total_len"])
        output_df['rnaseq_total_rpkm'] = output_df.rnaseq_total_counts / ((output_df.total_len/1000) * np.nansum(output_df.rnaseq_total_counts)) * 1000000
        print "...done!"
        return output_df
    elif run_mode == 'riboprof':
        output_df = pd.DataFrame(output_list, columns = ["gene_id", "gene_name", "ribo_utr5_counts", "utr5_len", "ribo_cds_counts", "cds_len", "ribo_utr3_counts", "utr3_len", "ribo_uORF_total_counts", "uORF_total_len", "ribo_uORF_utr5_counts", "uORF_utr5_len"])
        print "...done!"
        return output_df
    else:
        output_df = pd.DataFrame(output_list, columns = ["gene_id", "gene_name", "ribo_utr5_counts", "utr5_len", "ribo_cds_counts", "cds_len", "ribo_utr3_counts", "utr3_len", "ribo_uORF_total_counts", "uORF_total_len", "ribo_uORF_utr5_counts", "uORF_utr5_len", "rnaseq_cds_counts", "rnaseq_total_counts", "total_transcript_len"])
        total_ribo_reads = np.nansum(output_df.ribo_utr5_counts) + np.nansum(output_df.ribo_cds_counts) + np.nansum(output_df.ribo_utr3_counts)
        output_df['translation_efficiency'] = (output_df.ribo_cds_counts / total_ribo_reads) / (output_df.rnaseq_cds_counts / np.nansum(output_df.rnaseq_total_counts)) #cds_length divides out
        print '...done!'
        return output_df
    #

def main():
    """
    Description
    ------------------------------------------------------------------------
    A read counting program for RNA-seq or ribosome profiling applications. 
    
    --'rnaseq' mode --
    Returns a tab-delimited table of shape (number_genes_input, 7). Columns: [gene_id, gene_name, reads in meta-cds, length of meta-cds, reads in meta-transcript, length of meta-transcript, rpkm]. 
        Meta-transcripts are compiled for a given gene_id (attribute of plastid obj or field 9 of GTF) by grouping all full transcripts (--transcript_annotations) with that gene_id and collapsing into a single maximal spanning window.
        RPKM is definied for the total meta-transcript: (reads in meta-transcript) / (length of meta-transcript in kb * total number of reads in any annotated meta-transcript) * 1,000,000. 
    
    --'riboprof' mode--
    Returns a tab-delimited table of shape (number_genes_input, 12). Columns: [gene_id, gene_name, utr5_counts, utr5_len, cds_counts, cds_len, utr3_counts, utr3_len, meta_uORF_counts, meta_uORF_len, meta_uORF_counts_in5utronly_noCDS, meta_uORF_len_in5utronly_noCDS]. 
        utrs are defined as the meta-region constructed from 5utr and 3utr annotations of input transcript list. These regions are masked by the meta-CDS, to ensure that reads are not counted twice. 
        'meta_uORF' parameters refer to meta-regions constructed via a call to meta-uORF with input uORF_list. No masks are applied.
        'meta_uORF_in5utronly_noCDS' parameters are the result of masking meta-uORF region with the 5'utr, which is itself masked by the CDS. These parameters refer to the region of the meta-uORF which is in the 5'UTR and does not overlap CDS on any transcript for that gene. 
    
    --'both' mode--
    Returns a tab-delimited table of shape (number_genes_input, 16). Columns: [gene_id, gene_name, ribosome profiling utr5 counts, utr5_length, ribosome profiling cds_counts, cds_len, ribosome profiling utr3 counts, utr3_len, ribosome profiling meta_uORF_counts, meta_uORF_len, ribosome profiling meta_uORF_counts_in5utronly_noCDS, meta_uORF_len_in5utronly_noCDS, rnaseq counts cds, rnaseq counts full meta-transcript, length full meta-transcript, translation_efficiency]
        ribosome profiling parameters are precisely as documented above for 'riboprof' mode
        rnaseq parameters are as described above for 'rnaseq' mode
        full transcript length gives the length of the full meta-transcript constructed via a call to meta-roi using the input full_transcript_list. 
        translational efficiency follows the definition of Ingolia et al, 2011: (riboprof_reads_meta_cds / [length meta_CDS in kb * total number of ribosome reads in any region of any annotated transcript]) / (rnaseq_reads_meta_cds / [length meta_CDS in kb * total number of rnaseq reads in any region of any annotated transcript])
            Note that the Ingolia definition of TE normalizes to reads in any region of any annotated transcript (UTR or CDS), rather than the total number of reads that STAR could align. This avoids confounding by rRNA or small ncRNA, which STAR can align, but may differ between experiments.  
    ----------------------------------------------------------------------
    """
    parser = argparse.ArgumentParser(description=main.__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--run_mode', required=True, choices=['riboprof', 'rnaseq', 'both'], help='Required. str(riboprof|rnaseq|both). Run mode to use for counting.')
    parser.add_argument('--transcript_annotations', help='Required. Path to file containing transcript annotations in GTF2 format or as a list of plastid transcript objects that have been pickled using dill. The gene_id attribute (column 9 of GTF) will be used to group transcripts by gene- so please ensure this is correct. Ensembl annotations are preferred for this reason.', required=True)
    parser.add_argument('--transcript_format', choices=['GTF2','Obj'], help='Required. str(GTF2|Obj) specifying if transcript annotations are in GTF2 format or a serialized list of plastid transcripts. Usual security concerns apply when using pickled objects', required=True)
    parser.add_argument('--uORF_annotations', default=[], help='Path to file containing uORF annotations in GTF2 format or as a list of plastid segmentchain objects. Special keyword "none" (default) can be used if no uORF annotations are desired. If using a GTF file, be especially careful that the gene_id attribute from column 9 of the GTF is correct and agrees with the gene ID of the appropriate transcripts in your transcript_annotations, as this will be used to group the features and pair uORFs/full transcript annotations')
    parser.add_argument('--uORF_format', choices = ['GTF2','Obj'], help='Required only if uORF_annotations != "none", str(GTF2|Obj) specifying if uORF annotations are in GTF2 format or are a serilarized list of plastid segmentchain objects. Usual security concerns apply with de-pickling so be careful.')
    parser.add_argument('--path_to_ribprof_bam', default=None, help='String specifying path to BAM file containing ribosome profiling read alignments. Required if mode is "riboprof" or "both".')
    parser.add_argument('--path_to_psite', default=None, help='String specifying path to string giving path to file containing table of read lengths:psite offsets (generated by plastid psite script). Required if run_mode is "riboprof" or "both".')
    parser.add_argument('--path_to_rnaseq', default=None, help='String specifying path to BAM file containing RNA-seq read alignments. Required if mode is "rnaseq" or "both".')
    parser.add_argument('--rnamapfactory', default='fiveprime', choices=['fiveprime','center'], help="str('fiveprime'|'center') detailing which map factory to use for counting rnaseq data. 'fiveprime' (default) will use plastid's FivePrimeMapFactory, assigning a count of 1 to the 5'-most base of the read. 'center' will use plastid's CenterMapFactory which assigns a value of 1/L to every position on read of length L. You must use 'fiveprime' for tools that require integer values for counts, like DESeq2.")
    parser.add_argument('--nthreads', default= 4, type=int, help="Number of threads to be used for multiprocessing. Default is 4. For optimal performance set as max processors - 1. ")
    parser.add_argument('--outfile_prefix', required=True, help='Required. Prefix for disk location where output data table should be saved. The final .tsv will be appended by the program.')
    
    #Parse Args into Namespace object
    args = parser.parse_args()
    
    #Check that specified file output can safely be written.
    if os.path.exists(args.outfile_prefix + '.tsv'):
        raise ValueError('Output file already exists and will not be overwritten. Please specify a new filename.')
    
    #Write arguements to output file as commented lines. 
    fout = open(args.outfile_prefix + '.tsv', 'w')
    for arg, val in vars(args).iteritems():
        fout.write('#' + str(arg) + ': ' + str(val) + '\n')
    
    #Build Transcript annotations list from input GTF or load from serialized object
    if args.transcript_format == 'GTF2':
        print 'Assembling Transcripts from GTF2 File: ' + args.transcript_annotations
        annotation_list = list(GTF2_TranscriptAssembler(args.transcript_annotations))
        print '...done!'
    else:
        annotation_list = dill.load(open(args.transcript_annotations))
        print 'Successfullly assembled transcripts from serialized object: ' + args.transcript_annotations
    
    #Build uORF annotations list. uORF_annotations are not required, default will be list of zero. If string (ie. filepath) specified, parse based on format. 
    if isinstance(args.uORF_annotations,(str,)):
        if args.uORF_format == 'GTF2':
            print 'Assembling uORFs from GTF2 File: ' + args.uORF_annotations
            list_of_uORFs = list(GTF2_TranscriptAssembler(args.uORF_annotations))
            print '...done!'
        else:
            list_of_uORFs = dill.load(open(args.uORF_annotations))
            print 'Successfullly assembled uORFs from serilarized object: ' + args.uORF_annotations
    else:
        list_of_uORFs = args.uORF_annotations
    
    #Run the program
    output_frame = feature_counter(run_mode=args.run_mode, list_of_transcripts=annotation_list, list_of_uORFs=list_of_uORFs, path_to_ribprof_bam=args.path_to_ribprof_bam, path_to_psite=args.path_to_psite, path_to_rnaseq=args.path_to_rnaseq, rnamapfactory=args.rnamapfactory, nthreads=args.nthreads)
    
    #Save file
    output_frame.to_csv(fout, sep='\t', header=True, index=False)
    fout.close()
    
    print 'Program Complete.'

if __name__ == '__main__':
    main()


#When we get this working
#parser.add_argument('--check_overlap', default=False, action='store_true', help="When this flag is specified, two additional columns are added to output specifying overlap between a given gene's maximum spanning window and any other genes in annotation set. See output section above for description. Default is to skip this functionality" )

def cluster_segmentchain_families(meta_gene_list):
    """A function which finds the minimal number of non-overlapping clusters to describe a list of input plastid transcript objects. This is done iteratively: two overlapping objects are identified, grouped into a famlily, and family boundaries are expanded recursively to maximize family size. Optimized to find overlapping genes (represented as meta-transcripts) in read counting applications.  
    
    --Input--
    meta_gene_list: list of plastid transcript objects representing meta-genes (result of calling meta_roi on list of transcripts for a given gene). 
    
    --Output--
    A dictionary {key:value}. key is a string, 'family_chr_start_end_strand' based on maximum spanning genomic segment for the family. Value is a list of meta-genes within that family (using gene_id attribute). 
    """
    annotation_dictionary = {i.get_gene():i for i in meta_gene_list}
    genome_hash = GenomeHash(meta_gene_list) #block off spanning segments in genome for each known transcript
    families = defaultdict(list)
    assigned_to_family = []
    for my_metagene in meta_gene_list:
        print 'Working on gene ' + my_metagene.get_gene()
        
        #If gene has already been included in a family we can safely ignore it b/c family guarenteed to be bounded
        if my_metagene.get_name() in assigned_to_family:
            continue
        
        else:
            #check if spanning segment for my_metagene overlaps anything in the hash
            maximum_segment_family = my_metagene.spanning_segment
            continue_to_search = True
            while continue_to_search:
                #Find overlapping features in genome
                my_family = genome_hash.get_overlapping_features(maximum_segment_family, stranded=True)
                
                #We are searching for any features bigger than self to know if we should expand the family. 
                #THIS IS TOTALLY BROKEN!! CHECK STRAND, THEN COMPARE BOUNDS OF SPANNING SEGMENT TO FIND MAX. APPEND SPANNING SEGMENT END/START AS NECESSARY
                larger_features = [i for i in my_family if len(i.spanning_segment)>len(my_family)]
                
                #If no features larger than maximum_segment_family our search is done. Break the while loop.
                if len(larger_features) > 0:
                    continue_to_search = False
                
                #Otherwise adjust maximum_segment_family to reflect the largest spanning segment and continue search
                else: 
                    maximum_segment_family = max(larger_features, key=lambda item: len(item.spanning_segment))
                #
            #
            #At this point my_family contains all the transcripts in the family! Output to dictionary named on overlap 
            family_name = 'family_' + maximum_segment_family.chrom + '_' + str(maximum_segment_family.start) + '_' + str(maximum_segment_family.end) + '_' + maximum_segment_family.strand
            families[family_name] = my_family
            #
            #Reserve names of family to assigned_to_family list to narrow search space
            assigned_to_family.extend([i.get_name() for i in my_family])
        #
    #
    return families
