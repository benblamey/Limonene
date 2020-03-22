#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 15 14:06:03 2020

@author: Thomsn
"""
def progressinator(loop_count, process_atom, num_taxa):
    progress_criterion = -( - num_taxa // process_atom)
    if (loop_count % progress_criterion) == 1:
        percentage_factor = process_atom/100
        return f' - {int(loop_count // (progress_criterion*percentage_factor)) - 1} % -'
    else:
        return

def entreztive(my_mail,
            my_API_key,
            search_limit,
            query,
            e_type,
            ncbi_database,
            acc_id = None,
            sort_criterion = None):
    
    from Bio import Entrez
    Entrez.email = my_mail
    Entrez.api_key = my_API_key
    
    if e_type == 'esearch':   
        with Entrez.esearch(db=ncbi_database, term=query, retmax=search_limit, sort=sort_criterion) as handle:
            record = Entrez.read(handle)
              
    elif e_type == 'esummary':
        with Entrez.esummary(db=ncbi_database, id=query, retmax=search_limit) as handle:
            record = Entrez.read(handle)[0]
            
    else:
        with Entrez.efetch(db=ncbi_database, id=query, retmode='text', rettype="gb") as handle:
            record = SeqIO.read(handle, "genbank")
            
    return record

    
# Database search for gene and taxon:
def query_gen(my_mail,
              my_API_key,
              coding_sequence, 
              gene_name, 
              txid, 
              exclude_query, 
              ncbi_database,
              search_limit, 
              sort_criterion,
              ref_seq,
              exless_ref_seq):
    if coding_sequence:
        query = f'{gene_name}[Gene Name]+txid{str(txid)}'
    else:
        query = f'{gene_name}+txid{str(txid)}'
    if exclude_query:
        exclude_string = f'+NOT+{exclude_query}'
        query = f'{query}{exclude_string}'
    if ref_seq:
        rsquery = query
        query = f'{query}+srcdb​_refseq[prop]'
        rsq_perf = True
    gene_record = entreztive(my_mail, my_API_key, search_limit, query, 'esearch', ncbi_database, sort_criterion)
    
    if (ref_seq and gene_record['Count'] == 0 and not exless_ref_seq):
        gene_record = entreztive(my_mail, my_API_key, search_limit, rsquery, 'esearch', ncbi_database, sort_criterion)
        rsq_perf = False
    return gene_record, rsq_perf

# Collect data and safe as fasta:
def gene_hoover(features,
             record,
             full_taxon_list,
             data_list,
             genus_ID,
             ncbi_database,
             gene_info,
             newpath,
             taxon_query,
             gene_name,
             taxon_level,
             entries_per_tax,
             entries_per_genus):
    from Bio import Entrez
    Entrez.email = my_mail
    Entrez.api_key = my_API_key
    sample_locality = features[0].qualifiers.get('country')
    try:
        referenza = record.annotations['references'][0]
    except TypeError:
        referenza = ['']
    selection = full_taxon_list[full_taxon_list['taxID'] == genus_ID]
    data_list.loc[len(data_list)] = [genus_ID,
                        record.id,
                        len(record),
                        ncbi_database,
                        record.annotations['date'],
                        record.annotations['organism'],
                        referenza,
                        gene_info,
                        sample_locality,
                        selection['outgroup'].iat[0], # how safe is that with 0 (could be only datatrafo)??? Maybe I just do not undestand it (15/03/20)
                        selection['genus'].iat[0],
                        selection['epithet'].iat[0],
                        selection['refseq'].iat[0]]
    with open(f'{newpath}/{taxon_query}_{gene_name}.fasta', 'a') as finalfasta:
        rawseq = str(record.seq)
        if taxon_level == 'species':
            fasta_head = str(record.annotations['organism'])
            fasta_head = fasta_head.replace(' ', '_').replace('.', '').replace('-', '_')
        else:
            fasta_head = genus_ID
        if entries_per_tax == 1 and entries_per_genus == 1:
            finalfasta.write(f'>{fasta_head}\n')
        else:
            finalfasta.write(f'>{fasta_head}_{record.id}\n')
        finalfasta.write(f'{rawseq}\n\n')
    return data_list

# This function talks to you:
def printowaz(essay_num,
          entry_name = None,
          i = None,
          newpath = None,
          tree_resolution = None,
          taxonlist_path = None,
          taxon_query = None,
          gene_name = None,
          org = None,
          process_time = None):
    import pandas as pd
    param_name, restart, proti, prosi = [None]*4
    if entry_frame:
        param_name = entry_frame.index[i]
        restart = entry_frame[i]
    if process_time:
        proti = process_time.seconds // 60
        prosi = process_time.seconds % 60
    
    ESSAY = pd.Series([f'The parameter "{param_name}" was unfortunately set to a value \
which is out of the lumberjacks`range of possibilities. Please change it to one of the following and restart: \
{restart}.',
'The protein database is too tired today to sort by sequence length as you asked for. Please feel happy with default order.',
f'Step 0:\n >> Directory {newpath} was created.',
f'Mining data with tree resolution on the taxonomic level \'{tree_resolution}\' \
is not possible with this script. Try \'genus\' (default) or \'species\'.',
f'It was not possible to find input file {taxonlist_path}, please check \
the path and restart \'ncbi_miner\'.',
f'The csv-file {taxonlist_path} was loaded. Step 1 will be \
skipped. Follow up steps will use those taxa to search for sequences.',
f'Step 2:\n >> Esearch for gene entries of all species in the taxon {taxon_query}. \
This may take time, so keep the internet connection, chill down, drink a tea. \n The progress is ...',
'New Error, but RUN Forrest RUN!',
f' >> Done - entry database successfully established. Summary was saved as \
{taxon_query}_{gene_name}_stat.csv',
f'The input file {taxonlist_path} does not fit with the conditions (header).\
Please change and restart \'ncbi_miner\'.',
'Step 1:\n >> Done - taxon database successfully established.',
f'Step 2:\n >> Esearch for gene entries of all species in the taxon {taxon_query}. \
This may take time, so keep the internet connection, chill down, drink a tea. \n The progress is ...',
'New Error. Nemas Problemas!',
'New Error, but bro, stay seated, we skip it and keep searching...',
'New Error! But nothing big to worry about.',
'There are no entries of the gene for outgroup {taxon}.',
'Warning: You have no outgroups!',
f' >> Done - entry database successfully established. Summary (without outgroups) was saved as \
{taxon_query}_{gene_name}_stat_ingroup.csv',
f'Step 3:\n >> Efetch for each species of the given genera with the most entries \
for the given gene. This may take some time as well, I would answer some mails in the meantime ;-) \n The progress is ...',
'New Error, but chill down, everything is soft.',
'Updating the outgroup data...',
f'Error, did not find {org}.',
'Errore furore, no problemo spaghetto!',
'New Error! But nothing big to worry about.',
'Ulala, there was an error, check the outgroups in the stat/data-csv.files afterwards.',
f' >> Done - most recent fasta sequences were collected and successfully concatenated. \
It was saved in the file {taxon_query}_{gene_name}.fasta and is proove for nexus conversion. \
Summary of used sequences was saved as {taxon_query}_{gene_name}_data.csv. In addition the \
outgroup was added to the stat-file and saved as: {taxon_query}_{gene_name}_stat.csv.',
f' >> NCBI mining finished. It took {proti} min, \
{prosi} sec. The files are stored in the directory {newpath}.'])
    return ESSAY[essay_num]
    


query = f'{gene_name}[Gene Name]+txid{str(txid)}'

from Bio import Entrez
Entrez.email = my_mail
Entrez.api_key = my_API_key
    with Entrez.efetch(db='nucleotide', id='XM_001461340', retmode='text', rettype="gb") as handle:
        record = SeqIO.read(handle, "genbank")

record.features

