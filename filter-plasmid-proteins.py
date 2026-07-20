#!/usr/bin/env python

'''
filter-plasmid-proteins.py by Rohan Maddamsetti.

This script filters for plasmid proteins.

Usage: python filter-plasmid-proteins.py

'''

import os
from Bio import SeqIO
from tqdm import tqdm
import argparse
import gzip

outf = "../results/filtered-plasmid-proteins.csv"
print("outfile =>", outf)
    
## use the data in chromosome-plasmid-table.csv to look up replicon type,
## based on Annotation_Accession and then NCBI Nucleotide ID.

replicon_type_lookup_table = {}
chromosome_plasmid_tbl = "../data/Maddamsetti2025-results/chromosome-plasmid-table.csv"
with open(chromosome_plasmid_tbl, 'r') as chromosome_plasmid_fh:
    for i, line in enumerate(chromosome_plasmid_fh):
        if i == 0: continue ## skip the header
        line = line.strip()
        fields = line.split(',')
        my_annot_accession = fields[0]
        rep_type = fields[-2] 
        rep_id = fields[-3]
        if my_annot_accession in replicon_type_lookup_table:
            replicon_type_lookup_table[my_annot_accession][rep_id] = rep_type
        else:
            replicon_type_lookup_table[my_annot_accession] = {rep_id : rep_type}
            

with open(outf, 'w') as outfh:
    header = "SeqIndex,AnnotationAccession,SeqID,SeqType,product,sequence\n"
    outfh.write(header)
    seq_index = 1 ## index for sequences.
    for gbk in tqdm(os.listdir("../data/Maddamsetti2025-results/gbk-annotation")):
        if not gbk.endswith(".gbff.gz"): continue
        annotation_accession = gbk.split("_genomic.gbff.gz")[0]
        ## The next line is a consistency check to make sure that
        ## chromosome-plasmid-table.csv and the data in ../results/gbk-annotation are consisten
        if annotation_accession not in replicon_type_lookup_table: continue
        infile = os.path.join("../data/Maddamsetti2025-results/gbk-annotation/", gbk)
        with gzip.open(infile, "rt") as genome_fh:
            for replicon in SeqIO.parse(genome_fh, "gb"):
                """ hack to make the replicon id from biopython parsing. """
                replicon_id = ".".join([replicon.annotations['accessions'][0],str(replicon.annotations['sequence_version'])])
                replicon_type = "NA"                
                if replicon_id in replicon_type_lookup_table[annotation_accession]:
                    replicon_type = replicon_type_lookup_table[annotation_accession][replicon_id]
                else: ## replicon is not annotated as a plasmid or chromosome
                    ## in the NCBI Genome report, prokaryotes.txt.
                    ## assume that this is an unassembled contig or scaffold.
                    replicon_type = "contig"
                ## now skip replicons that are not plasmids.
                if replicon_type != "plasmid":
                    continue
                for feat in replicon.features:
                    if feat.type != "CDS": continue
                    try:
                        prot_seq = feat.qualifiers['translation'][0]
                    except:
                        continue
                    try: ## replace commas with semicolons! otherwise csv format is messed up.
                        prot_product = feat.qualifiers['product'][0].replace(',',';')
                    except:
                        prot_product = "NA"
                    ## now print to file.
                    row = ','.join([str(seq_index), annotation_accession, replicon_id, replicon_type, prot_product, prot_seq])
                    row = row + '\n'
                    outfh.write(row)
                    seq_index += 1 ## increment the seq_index after printing.

                        
