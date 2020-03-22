#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 27 12:18:17 2020

@author: Thomsn
"""

__author__ = 'thomas m. huber'
__email__ = ['thomas.huber@evobio.eu']

def ncbi_miner(taxon_query,
               gene_name,
               gene_length,
               my_mail,
               outgroups = [],
               coding_sequence = True,
               exclude_query = None,
               my_API_key = None,
               tree_resolution = 'genus',
               ncbi_database = 'nuccore',
               sortby = 'SLEN',
               ref_seq = True,
               exless_ref_seq = False,
               length_tolerance_percent = 20,
               upper_tolerance = 2,
               search_limit = 10000,
               entries_per_genus = 1,
               entries_per_tax = 1,
               taxonlist_path = None,
               #random_mining = False, I should implement this into sortby!
               strict_search = True):
    import os
    import pandas as pd
    from Bio import SeqIO
    from datetime import datetime
    from urllib.error import HTTPError
    if random_mining:
        from random import shuffle
    
    from funzioni import progressinator
    from funzioni import entreztive
    from funzioni import query_gen
    from funzioni import gene_hoover
    from funzioni import printowaz
        
##### FUNCTIONS #####

# Constants used
    TAX_COLUMNS = ['taxID', 'genus', 'epithet', 'entry_UIDs', 'count', 'outgroup', 'refseq']

    ENTRY_FRAME = pd.Series([['genus', 'species', 'GENUS', 'SPECIES'],
               ['nuccore', 'protein', 'NUCCORE', 'PROTEIN'],
               [None, 'slen', 'pdat', 'random', 'SLEN', 'PDAT', 'RANDOM']])
    ENTRY_FRAME.index = ['tree_resolution', 'ncbi_database', 'sortby']
    

# Catching the time/date
    start_time = datetime.now()
    
    
# Avoid unpleasant choices    
    for i, entry in enumerate([tree_resolution, ncbi_database, sortby]):
        if entry not in ENTRY_FRAME[i]:
            printowaz(0, entry_frame = ENTRY_FRAME, i = i)
            return
    if ncbi_database in ['protein', 'PROTEIN']:
        if sortby in ['slen', 'SLEN']:
            sortby = None
            printowaz(1)

# Assign a new folder for the data, with the format date_taxon_gene. Non overwriting.
    newpath = f'{datetime.now().strftime("%Y_%m_%d")}_{taxon_query}_{gene_name}'
    if not os.path.exists(newpath):
        os.makedirs(newpath)
        printowaz(2, newpath = newpath)

# Filters for gene_length are set. 
# Lower border: gene_length - (gene_length * length_tolerance_percent/100)
# Upper border: gene_length + (gene_length * length_tolerance_percent/100) * upper_tolerance
    gl_lowend = gene_length * (100 - length_tolerance_percent)/100
    gl_highend = gene_length * (100 + upper_tolerance * length_tolerance_percent)/100
    
# Will the tree be build up by species or (species reprecenting) genera?
    tree_resolution = str(tree_resolution).lower()


    if tree_resolution == 'genus':
        taxon_level = 'Genus'
    elif tree_resolution == 'species':
        taxon_level = 'species'
    else:
        print(3, tree_resolution = tree_resolution)

# Taxonlist was already previously established, so some steps are going to be skipped.
    if taxonlist_path:
        try:
            old_taxon_list = pd.read_csv(taxonlist_path)
        except FileNotFoundError:
            printowaz(4, taxonlist_path = taxonlist_path)
            return
        else:
            if all([(any(old_taxon_list.keys() == i )) for i in TAX_COLUMNS]):
                printowaz(5, taxonlist_path = taxonlist_path)

                Entrez.email = my_mail
                Entrez.api_key = my_API_key
                all_taxaIDs = list(old_taxon_list['taxID'])
                        
                printowaz(6, taxon_query = taxon_query)
                taxon_list = pd.DataFrame(columns = TAX_COLUMNS)
                PROCESS_ATOM = 20
                for i, taxon in enumerate(all_taxaIDs):
                    progressinator(i, PROCESS_ATOM, len(all_taxaIDs))
                    try: 
                        gene_record, rsq_perf = query_gen(my_mail, my_API_key, coding_sequence, gene_name, taxon, exclude_query, ncbi_database, search_limit, sortby, ref_seq)
                    except HTTPError:
                        printowaz(7)
                    else:
                        count = gene_record['Count']
                        if count != '0':
                            taxon_list.loc[taxon] = [taxon,
                                                 old_taxon_list['genus'].iloc[i],
                                                 old_taxon_list['epithet'].iloc[i],
                                                 gene_record['IdList'],
                                                 int(gene_record['Count']),
                                                 old_taxon_list['outgroup'].iloc[i],
                                                 rsq_perf]
                taxon_list.to_csv(f'{newpath}/{taxon_query}_{gene_name}_stat.csv')
                og_list = taxon_list.iloc[:,:][taxon_list.outgroup]
                taxon_list = taxon_list.drop(og_list.index)
                printowaz(8, taxon_query = taxonquery, gene_name = gene_name)

            else:
                printowaz(9, taxonlist_path = taxonlist_path)
                return

    else:
        #### > STEP 1 < ####

        que = f'{taxon_query}[orgn]'
        try:
            record = entreztive(my_mail, my_API_key, search_limit, que, 'esearch', 'taxonomy')
        except HTTPError:
            return str('Database error, try later...')
        else:
            all_taxaIDs = record['IdList']

            printowaz(10)
        
    
        #### > STEP 2 < ####
        
        PROCESS_ATOM = 20

        printowaz(6, taxon_query = taxon_query)
        taxon_list = pd.DataFrame(columns = TAX_COLUMNS)
        for i, taxon in enumerate(all_taxaIDs, 1):
            progressinator(i, PROCESS_ATOM, len(all_taxaIDs))
            
            try:
                record = entreztive(my_mail, my_API_key, search_limit, taxon, 'esummary', 'taxonomy')
            except IndexError:
                printowaz(11)
                break
            except HTTPError:
                printowaz(12)
            else:
                if record['Rank'] == 'species' and record['Genus'] != '':
                    try: 
                        gene_record, rsq_perf = query_gen(my_mail, my_API_key, coding_sequence, gene_name, taxon, exclude_query, ncbi_database, search_limit, sortby, ref_seq)
                    except HTTPError:
                        printowaz(13)
                    else:
                        count = gene_record['Count']
                        if count != '0':
                            taxon_list.loc[taxon] = [taxon,
                                                 record['Genus'],
                                                 record['Species'],
                                                 gene_record['IdList'],
                                                 int(gene_record['Count']),
                                                 'False',
                                                 rsq_perf]
                else:
                    pass
        og_list = pd.DataFrame(columns = TAX_COLUMNS)
        if len(outgroups) > 0:
            for taxon in outgroups:
                try: 
                    gene_record, rsq_perf = query_gen(my_mail, my_API_key, coding_sequence, gene_name, taxon, exclude_query, ncbi_database, search_limit, sortby, ref_seq)
                except HTTPError:
                    printowaz(13)
                else:
                    count = gene_record['Count']
                    if count == '0':
                        printowaz(14, taxon)
                    else:
                        og_list.loc[len(og_list)] = [taxon,
                                            f'{taxon}_gen',
                                            f'{taxon}_sp',
                                            gene_record['IdList'],
                                            int(gene_record['Count']),
                                            'True',
                                            rsq_perf]

        else:
            printowaz(15)
        #for taxon in enumerate(all_taxaIDs, 1):
            taxon_list.to_csv(f'{newpath}/{taxon_query}_{gene_name}_stat_ingroup.csv')
        printowaz(16, taxon_query = taxon_query, gene_name = gene_name)
        
    
    #### > STEP 3 < ####

    printowaz(17)
    data_list = pd.DataFrame(columns = ['taxID',
                                        'accession',
                                        'length',
                                        'database'
                                        'date',
                                        'organism',
                                        'reference',
                                        'gene_information',
                                        'sampling_locality',
                                        'outgroup',
                                        'genus',
                                        'epithet',
                                        'refseq'])

    full_taxon_list = pd.concat([taxon_list, og_list])
    if taxon_level == 'Genus':
        genera = full_taxon_list['genus'].unique()
    else:
        genera = full_taxon_list['taxID'].unique()

    PROCESS_ATOM = 10
    for j, genus_ID in enumerate(genera):
        progressinator(i, PROCESS_ATOM, len(all_taxaIDs))
        if taxon_level == 'Genus':
          all_species = full_taxon_list[full_taxon_list['genus'] == genus_ID]
        else:
          all_species = full_taxon_list[full_taxon_list['taxID'] == genus_ID]
          
        all_species = all_species.sort_values(by=['count'], ascending=False)
        if len(list(all_species['entry_UIDs'])) < entries_per_genus:
            epg = len(list(all_species['entry_UIDs']))
        else:
            epg = entries_per_genus
            
        for entry in list(range(epg)):
            entry_list = list(all_species['entry_UIDs'])[entry]
            if sortby in ['random', 'RANDOM']:
                shuffle([1,2,3])
                shuffle(entry_list)
                    
            entry_counter = 0
            for i, acc_id in enumerate(entry_list):
                try:
                    record = entreztive(my_mail, my_API_key, search_limit, acc_id, 'efetch', ncbi_database)
                except HTTPError:
                    printowaz(18)
                else:
                    if gl_lowend < len(record) < gl_highend:
                        features = record.features
                        filterframe = [x.type == 'gene' for x in features]
                        if any(filterframe) == True:
                            keyword = 'gene'
                        else:
                            filterframe = [x.type != 'source' for x in features]
                            keyword = 'product'
                        gene_features = [x for i, x in enumerate(features) if filterframe[i]==True]
                        try:
                            gene_info = [(x.qualifiers.get(keyword))[0] for x in gene_features]
                        except TypeError:
                            gene_info = ['']
                        if strict_search:
                            if (len(gene_info) == 1 and gene_name.lower() in gene_info[0].lower()):
                                gene_hoover(features,
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
                                             entries_per_genus)
                                entry_counter += 1 #attention! What happens with multiple entries per species????
                                if entry_counter == entries_per_tax:
                                    break
                                
                        else:
                            gene_hoover(features,
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
                                         entries_per_genus)
                            entry_counter += 1 #attention! What happens with multiple entries per species????
                            if entry_counter == entries_per_tax:
                                break
                    else:
                        pass
    printowaz(19)
    out_selection = data_list[data_list['outgroup'] == 'True']
    for i, acc in enumerate(out_selection['accession']):
        org = out_selection.loc[out_selection['accession'] == acc, ['organism']]
        org = ''.join(str(org.values[0])).replace(' ','_').replace("'",'')
        que = f'{org}'
        try:
            record = entreztive(my_mail, my_API_key, search_limit, que, 'esearch', 'taxonomy')
        except HTTPError:
            printowaz(20, org = org)
        else:
            txid = record['IdList'][0]
            try:
                spef_record = entreztive(my_mail, my_API_key, search_limit, txid, 'esummary', 'taxonomy')
            except IndexError:
                printowaz(11)
                break
            except HTTPError:
                printowaz(21)
            else:
                try:
                    gene_record, rsq_perf = query_gen(my_mail, my_API_key, coding_sequence, gene_name, txid, exclude_query, ncbi_database, search_limit, sortby, ref_seq)
                except HTTPError:
                    printowaz(22)
                else:
                    count = gene_record['Count']
                    if count == '0':
                        printowaz(23)
                        pass
                    else:
                        og_gen = spef_record['Genus']
                        og_sp = spef_record['Species']
                        data_list.at[data_list['accession'] == acc, ['taxID']] = txid
                        data_list.at[data_list['accession'] == acc, ['genus']] = og_gen
                        data_list.at[data_list['accession'] == acc, ['epithet']] = og_sp
                        data_list.at[data_list['accession'] == acc, ['refseq']] = og_rs
                        taxon_list.loc[txid] = [txid,
                                              og_gen,
                                              og_sp,
                                              gene_record['IdList'],
                                              int(gene_record['Count']),
                                              'True',
                                              'og_rs']

    data_list.to_csv(f'{newpath}/{taxon_query}_{gene_name}_data.csv')
    taxon_list.to_csv(f'{newpath}/{taxon_query}_{gene_name}_stat.csv')
    printowaz(24, taxon_query = taxon_query, gene_name = gene_name)
    stop_time = datetime.now()
    process_time = stop_time - start_time
    printowaz(25, newpath = newpath, process_time = process_time)

