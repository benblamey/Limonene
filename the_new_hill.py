#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 16 10:04:58 2020

@author: Thomsn
"""

import os
import python_cipres.client as cip
import pandas as pd
from Bio import AlignIO
from Bio import SeqIO
import numpy as np

home = os.path.expanduser('~')

os.chdir(home)

url = 'https://cipresrest.sdsc.edu/cipresrest/v1'
appname = 'treepinkler'
appid = 'treepinkler-EE06FD418FD147B99EE41EA2B3B40E4B'
username = 'biotomme'
password = 'vacc1Niumm'

date_name = '2019_12_22'
gene_names = ['rpb1', 'rpb2', '28s']
taxon_name = 'Boletales'

pycipres_access = f'''URL={url}
APPNAME={appname}
APPID={appid}
USERNAME={username}
PASSWORD={password}
VERBOSE='''


#### Activate your CIPRES API ####
def wake_up_cipres(pycipres_access,
                   appname,
                   appid,
                   username,
                   password,
                   url):
    with open('pycipres.conf', 'w') as file:
        file.write(pycipres_access)
    global verbose
    properties = cip.Application().getProperties()
    client = cip.Client(appname, appid, username, password, url)
    if properties.VERBOSE:
        verbose = True
    return client


#### Finally navigate to your main directory with ncbi_lumberjack output folders ####


#### Do you have additional sequences ####

add_seq_choice = 'True'


#### Add additional sequences - place sequences in folder 'additional' main directory and name '{gene_name}.fasta' ####
def sequence_additor(add_seq_choice,
                     gene_name,
                     working_dir,
                     date_name,
                     taxon_name,
                     username,
                     my_mail):
    if add_seq_choice:
        gene_dir = f'{working_dir}/{date_name}_{taxon_name}_{gene_name}'
        records = list(SeqIO.parse(f'{gene_dir}/{taxon_name}_{gene_name}.fasta', 'fasta'))
        taxon_list = pd.read_csv(f'{gene_dir}/{taxon_name}_{gene_name}_data.csv')
        try:
            additional_records = list(SeqIO.parse(f'{working_dir}/additional/{gene_name}.fasta', 'fasta'))
        except FileNotFoundError:
            additional_records = []
            print(f'> {gene_name}: WARNING: There were no additional files for the gene {gene_name}.')
        all_records = records + additional_records
        with open(f'{gene_dir}/{taxon_name}_{gene_name}_wadd.fasta', 'w') as output_handle:
            SeqIO.write(all_records, output_handle, 'fasta')
        for add_rec in additional_records:
            taxon_list.loc[add_rec.id] = ['NA',
                                          'NA',
                                          len(add_rec.id),
                                          ncbi_database,
                                          'NA',
                                          add_rec.id,
                                          'Current user: {username}, {my_mail}',
                                          'This is a manually added sequence. You could have more information by yourself.',
                                          'NA',
                                          False,  # Only for ingroups!!!!
                                          'NA',
                                          'NA',
                                          False]
        taxon_list.to_csv(f'{gene_dir}/{taxon_name}_{gene_name}_data_wadd.csv')
        print(f'> {gene_name}: {len(additional_records)} sequences (names: {[ent.id for ent in additional_records]}) \
were added to the fasta-file and saved with postfix \'_wadd\'. Furthermore, \
some information of the additional sequences was added to the data.csv file and stored in \
{gene_dir}/{taxon_name}_{gene_name}_data_wadd.csv. Chillig!')


#### MAFFT every sequence ####
def the_big_maffti(client,
                   gene_names,
                   working_dir,
                   add_seq_choice,
                   date_name,
                   taxon_name):
    mafft_progress_cipres = pd.Series()
    for gene_name in gene_names:  # starts the CIPRES run
        gene_dir = f'{working_dir}/{date_name}_{taxon_name}_{gene_name}'
        if add_seq_choice:
            fill_in_file = f'{gene_dir}/{taxon_name}_{gene_name}_wadd.fasta'
        else:
            fill_in_file = f'{gene_dir}/{taxon_name}_{gene_name}.fasta'

        mafft_progress_cipres[gene_name] = client.submitJob({'toolId': 'MAFFT_TG',
                                                             'runtime_': '.1',
                                                             'adjust_direction_': '1'},
                                                            {'infile_': fill_in_file},
                                                            {'statusEmail': 'true'},
                                                            validateOnly=False)
        cipres_lyrics(mafft_progress_cipres[gene_name],
                      'mafft',
                      gene_name,
                      'NA',
                      '',
                      fill_in_file)
        print(f'> Job for MAFFT alignment of {gene_name} was successfully \
submitted to CIPRES.')
    return mafft_progress_cipres


## Is the alignment process finished? Do you want to queue in the waiting process? ##
def mafft_checker(gene_names,
                  mafft_progress_cipres,
                  wait_until_finished=False):
    for gene_name in gene_names:
        print(f'------------------------- {gene_name} -------------------------')
        mafft_progress_cipres.at[gene_name].update()
        print(mafft_progress_cipres.at[gene_name].show(messages=True))
        if wait_until_finished:
            print(f'Please wait, we are queued until the MAFFT alignment for {gene_name} is finished. \
It is only about minutes...')
            mafft_progress_cipres.at[gene_name].waitForCompletion()
            print(f'> CIPRES successfully aligned your {gene_name} sequences - \
you can follow the next steps to download it.')


## Download alignments from CIPRES ##
def dload_maffto(gene_names,
                 working_dir,
                 date_name,
                 taxon_name,
                 mafft_progress_cipres):
    for gene_name in gene_names:
        print(f'------------------------- {gene_name} -------------------------')
        gene_dir = f'{working_dir}/{date_name}_{taxon_name}_{gene_name}'
        newpath = f'{gene_dir}/{taxon_name}_{gene_name}_mafft'
        if not os.path.exists(newpath):
            os.makedirs(newpath)
            print(f'>> Directory {newpath} was created.')
        if not mafft_progress_cipres.at[gene_name].terminalStage:
            return f'The MAFFT process for {gene_name} is not yet finished. Sing a song \
about todays weather and come back afterwards to find your alignment...'

        mafft_progress_cipres.at[gene_name].downloadResults(directory=newpath)
        print(f'CIPRES files were downloaded. The alignment is stored in \
directory {newpath} and can be found named \'output.fasta\'.')


## How large is the alignment?
def ntax_nchar(alignment):
    ntax = int(len(alignment))
    nchar = int(len(alignment[0].seq))
    return ntax, nchar


## Load the alignment form directory for usage
def spit_aligned(gene_names,
                 working_dir,
                 date_name,
                 taxon_name,
                 ncbi_database):
    from Bio.Alphabet import IUPAC

    alignments = pd.Series()
    for gene_name in gene_names:
        gene_dir = f'{working_dir}/{date_name}_{taxon_name}_{gene_name}'
        proc_mafft_file = f'{gene_dir}/{taxon_name}_{gene_name}_mafft_no_outliers.fasta'
        mafft_file = f'{gene_dir}/{taxon_name}_{gene_name}_mafft/output.mafft'
        try:
            input_handle = open(proc_mafft_file, "r")
        except FileNotFoundError:
            input_handle = open(mafft_file, "r")
            print(f'> {gene_name}: No processed fasta file was found, therefore the raw file \
{mafft_file} used for the alignment.')
        else:
            print(f'> {gene_name}: The processed fasta {proc_mafft_file} file was found and used for the alignment.')
        if ncbi_database[0].lower() == 'nuccore':
            alph = IUPAC.extended_dna
        elif ncbi_database[0].lower() == 'protein':
            alph = IUPAC.extended_protein
        else:
            return f'The database / sequence type (\'ncbi_database\' = {ncbi_database}) you specified was not found. \
Try \'nuccore\' for DNA or \'protein\' for protein sequences'
        alignments.loc[gene_name] = AlignIO.read(input_handle, "fasta", alphabet=alph)
        print(f'> Alignment for {gene_name} was loaded into memory.')
    return alignments


#### Countlister
def countlister(alignment):
    count_list = []
    _, nchar = ntax_nchar(alignment)
    for i in range(int(nchar)):
        base_count = 0
        for record in alignment:
            if record.seq[i] != '-':
                base_count += 1
        count_list.append(base_count)
    return count_list


#### Plot alignment ####
def alignoplot(gene_names,
               alignments):
    import matplotlib.pyplot as plt

    nchar = pd.Series(index=gene_names)
    ntax = pd.Series(index=gene_names)
    base_proportion = pd.Series()
    for gene_name in gene_names:
        ntax.loc[gene_name], nchar.loc[gene_name] = ntax_nchar(alignments.loc[gene_name])
        count_list = countlister(alignments.loc[gene_name])
        base_proportion[gene_name] = [i / ntax.loc[gene_name] for i in count_list]
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        ax.set_xlabel(f'Base site index (total length = {nchar.loc[gene_name]})', fontsize=8)
        ax.set_ylabel(f'Proportion of sequences with base (n = {ntax.loc[gene_name]})', fontsize=8)
        ax.set_title(f'Mafft alignment of {gene_name} in {taxon_name}', fontsize=10)
        print('WARNING: Please check the alignment also manually, the follow up automated analysis is not bulletproof!')
        plt.plot(base_proportion[gene_name])


## visualize the alignment
def visualign(alignments,
              gene_names,
              ncbi_database,
              working_dir,
              date_name,
              taxon_name):
    import plotly.graph_objects as go

    if ncbi_database[0].lower() != 'nuccore':
        return 'This function only supports visualization of DNA alignments - sorry.'
    else:
        print('You initiated a visualisation of the alignment. Please be patient, this may take some minutes to plot...')
        olignments = alignments.copy()
        for gene_name in gene_names:
            print(f'------------------------- {gene_name} -------------------------')
            alignment = alignments.loc[gene_name]
            ntax, nchar = ntax_nchar(alignment)
            id_list = [seq.id for seq in alignment]

            encolored = pd.DataFrame()
            colors = {'A': [255, 0, 0],
                      'T': [0, 255, 0],
                      'G': [0, 0, 255],
                      'C': [255, 0, 255],
                      '-': [255, 255, 255]}
            texto = []
            x_axis = []
            y_axis = []
            for j, sequ in enumerate(alignment):
                texto += [i.upper() for x, i in enumerate(sequ)]
                x_axis += [i for i, x in enumerate(sequ)]
                y_axis += [j for i in sequ]
                encolored[sequ.id] = [colors[i.upper()] if i.upper() in ['A', 'T', 'G', 'C', '-'] else [255, 255, 255] for i in sequ]

            encolored = encolored.transpose()
            fig = go.Figure(data=go.Image(z=encolored))
            fig.add_trace(go.Scatter(x=x_axis,
                                     y=y_axis,
                                     opacity=0,
                                     text=texto,
                                     hovertemplate='Sequence %{y} at pos %{x}: %{text}',
                                     name='Base'))
            fig.show()
            olignments.loc[gene_name] = eliminoutliers(gene_name, alignments, working_dir, date_name, taxon_name, outlier_sequences=id_list, delete_history=False)
    return olignments


#### Cut alignments etc. ####
def alignmentrimmer(gene_names,
                    alignments,
                    cut_alignment=True,
                    treshold_percentage=20):
    nchar = pd.Series(index=gene_names)
    ntax = pd.Series(index=gene_names)
    olignments = alignments.copy()
    for gene_name in gene_names:
        ntax.loc[gene_name], nchar.loc[gene_name] = ntax_nchar(alignments.loc[gene_name])
        if cut_alignment:
            count_list = countlister(alignments.loc[gene_name])
            base_proportion = [i / ntax.loc[gene_name] for i in count_list]
            precut = int(next(i for i, base in enumerate(base_proportion) if (base * 100) >= treshold_percentage))
            postcut = int(nchar.loc[gene_name] - next(i for i, base in enumerate(reversed(base_proportion)) if (base * 100) >= treshold_percentage))
            olignments.loc[gene_name] = alignments.loc[gene_name][:, precut:postcut]
            print(f'The alignment was trimmed in the beginning until base {precut} \
and in the end from base {postcut}. This means at least {treshold_percentage} % of all sequences \
have a present base there.')
    return olignments


## Where is the BAD FIT? ###
def weird_comparison(gene_names,
                     alignments,
                     working_dir,
                     date_name,
                     taxon_name,
                     count_list,
                     percentage=2,
                     quantile_tune=.005):
    import matplotlib.pyplot as plt
    from adjustText import adjust_text

    nchar = pd.Series(index=gene_names)
    ntax = pd.Series(index=gene_names)
    percentage = percentage / 100
    outlier_sequences = pd.Series()

    for gene_name in gene_names:
        print(f'------------------------- {gene_name} -------------------------')
        count_list = countlister(alignments.loc[gene_name])
        ntax.loc[gene_name], nchar.loc[gene_name] = ntax_nchar(alignments.loc[gene_name])
        base_proportion = [i / ntax.loc[gene_name] for i in count_list]
        base_count_list = pd.DataFrame(columns=['Base_count', 'Antibase_count'])
        antibase_count_list = []

        for record in alignments.loc[gene_name]:
            antibase_count = []
            base_count = []
            for i in range(int(nchar.loc[gene_name])):
                if record.seq[i] == '-':
                    antibase_count.append(base_proportion[i])
                else:
                    base_count.append(base_proportion[i])
            base_count_list.loc[record.id] = [np.quantile(base_count, .5),
                                              # Median occurence proportion of sequences per base site where there is none for the present record - expected low for low quality fit
                                              np.quantile(antibase_count,
                                                          .5)]  # Median occurence proportion of sequences per base site where there is none for the present record - expected high for low quality fit

        fig = plt.figure(figsize=(8, 8))
        ax = fig.add_subplot(1, 1, 1)
        ax.set_xlabel(f'Base median', fontsize=15)
        ax.set_ylabel(f'Antibase median', fontsize=15)
        ax.set_title(f'Quality of {gene_name} alignment fit', fontsize=20)
        texts = []
        ax.scatter(base_count_list['Base_count'],
                   base_count_list['Antibase_count'],
                   s=5)
        outliers = np.quantile(base_count_list['Base_count'], quantile_tune)  # np.quantile(base_count_list['Base_count'], .01)
        outlier_sequences[gene_name] = pd.Series()
        outcount = 0
        if sum(base_count_list.Base_count <= outliers) < (ntax.loc[gene_name] * percentage):
            for i, spe in enumerate(base_count_list.index):
                if base_count_list['Base_count'].iloc[i] <= outliers:  # 'Base'
                    outlier_sequences[gene_name].loc[outcount] = spe
                    texts.append(ax.annotate(f'{outcount}_{spe}',
                                             ((base_count_list['Base_count'].iloc[i]),
                                              (base_count_list['Antibase_count'].iloc[i]))))
                    outcount += 1
            adjust_text(texts,
                        arrowprops=dict(arrowstyle="-", color='black', lw=1))
        outliers = np.quantile(base_count_list['Antibase_count'], 1 - quantile_tune)  # np.quantile(base_count_list['Base_count'], .01)
        if sum(base_count_list.Antibase_count >= outliers) < (ntax.loc[gene_name] * percentage):
            for i, spe in enumerate(base_count_list.index):
                if (base_count_list['Antibase_count'].iloc[i] >= outliers and len(outlier_sequences[gene_name].values) == 0):
                    outlier_sequences[gene_name].loc[outcount] = spe
                    texts.append(ax.annotate(f'{outcount}_{spe}',
                                             ((base_count_list['Base_count'].iloc[i]),
                                              (base_count_list['Antibase_count'].iloc[i]))))
                    outcount += 1
                elif (base_count_list['Antibase_count'].iloc[i] >= outliers and all(outlier_sequences[gene_name].values != spe)):
                    outlier_sequences[gene_name].loc[outcount] = spe
                    texts.append(ax.annotate(f'{outcount}_{spe}',
                                             ((base_count_list['Base_count'].iloc[i]),
                                              (base_count_list['Antibase_count'].iloc[i]))))
                    outcount += 1
            adjust_text(texts,
                        arrowprops=dict(arrowstyle="-", color='black', lw=1))
            print(f'''The outliers for {gene_name} are: \n{outlier_sequences[gene_name]}''')

            alignment = eliminoutliers(gene_name, alignments, working_dir, date_name, taxon_name, outlier_sequences.at[gene_name], delete_history=False)
            print('WARNING: Please check the alignment also manually, the follow up automated analysis is not bulletproof!')
    return alignment


# Get rid of outliers
def eliminoutliers(gene_name,
                   alignments,
                   working_dir,
                   date_name,
                   taxon_name,
                   outlier_sequences,
                   delete_history=False):
    from Bio.Align import MultipleSeqAlignment
    input_list = []
    print('WARNING: Please check the alignment also manually, this automated analysis is not bulletproof!')
    ip = input('>> Do you want to delete any of the outlier sequences from your alignment?\n \
If yes, tell me the number. Otherwise press enter to skip and proceed\n')
    if ip:
        try:
            input_list += [int(ip)]
        except ValueError:
            print('Error, you did not enter an integer (= number without fractional component). Retry...')
        except TypeError:
            print('Error, you did not enter an integer (= number without fractional component). Retry...')
        else:
            print(f'-- Input {ip} accepted --')
    else:
        print('No entry, the alignment was not changed.')
    while ip:
        ip = input('-->> Do you want to delete another number? Please enter it or press enter to skip and proceed.\n')
        if ip:
            try:
                input_list += [int(ip)]
            except ValueError:
                print('Error, you did not enter an integer (= number without fractional component). Retry...')
            except TypeError:
                print('Error, you did not enter an integer (= number without fractional component). \
Please be so kind to try it again')
            else:
                print(f'-- Input {ip} accepted --')
    else:
        print(f'Aquisation of sequences to delete stopped. These numbers were entered: {input_list}.')

    if input_list:
        ids = [seq.id for seq in alignments.loc[gene_name]]
        outis = [outlier_sequences[i] for i in input_list]
        delete_it = [[i for i, x in enumerate(ids) if x == outi][0] for outi in outis]
        gene_dir = f'{working_dir}/{date_name}_{taxon_name}_{gene_name}'
        print('Do you really want to delete the following sequences?')
        for i, record in enumerate(alignments.loc[gene_name]):
            if i in delete_it:
                print(record)
        deal = input('Type in \'Y\' or \'y\' to finally take these entries out of your alignment. Input anything else to keep them.\n')
        if deal:
            if deal in 'Yy':
                new_alignment = []
                for i, record in enumerate(alignments.loc[gene_name]):
                    if not i in delete_it:
                        new_alignment.append(alignments.loc[gene_name][i])
                seq_data = pd.read_csv(f'{gene_dir}/{taxon_name}_{gene_name}_data.csv')
                selection = [x.replace('-', '_').replace(' ', '_').replace('.', '') for x in seq_data.organism]
                delete_data = [i for i, x in enumerate(selection) if x in outis]
                seq_data.iloc[delete_data, 1:].to_csv(f'{gene_dir}/{taxon_name}_{gene_name}_data_deleted.csv')
                seq_data.drop(delete_data).iloc[:, 1:].to_csv(f'{gene_dir}/{taxon_name}_{gene_name}_data_selected.csv')
                alignments.loc[gene_name] = MultipleSeqAlignment(new_alignment)
                SeqIO.write(alignments.loc[gene_name], f'{gene_dir}/{taxon_name}_{gene_name}_mafft_no_outliers.fasta', 'fasta')
                print(f'The sequences were successfully deleted. Further analysis will be conducted without those.')
    return alignments.loc[gene_name]


## Concatenate all genes
def concaterpillar(alignments,
                   gene_names,
                   working_dir,
                   date_name,
                   taxon_name,
                   ncbi_database):
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq
    from Bio.Alphabet import IUPAC
    from Bio.Align import MultipleSeqAlignment

    print(f'------------------------- {" - ".join(gene_names)} -------------------------')
    ntax = pd.Series(index=gene_names)
    nchar = pd.Series(index=gene_names)
    al_ind = []
    for gene_name in gene_names:
        ntax.loc[gene_name], nchar.loc[gene_name] = ntax_nchar(alignments.loc[gene_name])
        al_ind += [rec.id for rec in alignments.loc[gene_name]]

    al_ind_cat = list(dict.fromkeys(al_ind))

    alignment_cat = None
    for taxon in al_ind_cat:
        seq_cat = ''
        total_gene_name = ''
        for gene_name in gene_names:
            total_gene_name += f'{gene_name}_'
            if len([rec for rec in alignments.loc[gene_name] if rec.id == taxon]) >= 1:
                curseq = [rec for rec in alignments.loc[gene_name] if rec.id == taxon][0]
            else:
                if ncbi_database[0].lower() == 'nuccore':
                    alph = IUPAC.extended_dna
                elif ncbi_database[0].lower() == 'protein':
                    alph = IUPAC.extended_protein
                curseq = SeqRecord(seq=Seq('-' * int(nchar.loc[gene_name]), alph),
                                   id=taxon,
                                   name=taxon,
                                   description=taxon,
                                   dbxrefs=[])
            seq_cat += curseq
        if alignment_cat:
            alignment_cat.append(seq_cat)
        else:
            alignment_cat = MultipleSeqAlignment([seq_cat])

    total_align_file = f'{working_dir}/{taxon_name}_{total_gene_name}mafft'
    SeqIO.write(alignment_cat, f'{total_align_file}.fasta', 'fasta')
    print(f'> The total alignment was produced by concatenating the alignments \
of the genes {", ".join(gene_names)}. It was saved as {taxon_name}_{total_gene_name}mafft.fasta.')
    return alignment_cat


#### Fasta to NEXUS ####
def to_nexus(alignment,
             working_dir,
             taxon_name,
             gene_names,
             concat=True):
    if concat:
        total_gene_name = '_'.join(gene_names) + '_'
        total_align_file = f'{working_dir}/{taxon_name}_{total_gene_name}mafft'
        nexus_file = to_nexus_help(alignment,
                                   working_dir,
                                   taxon_name,
                                   gene_names,
                                   total_align_file)
    else:
        for gene_name in gene_names:
            total_align_file = f'{working_dir}/{date_name}_{taxon_name}_{gene_name}/{taxon_name}_{gene_name}mafft'
            nexus_file = to_nexus_help(alignment.loc[gene_name],
                                       working_dir,
                                       taxon_name,
                                       gene_names,
                                       total_align_file)
    return


def to_nexus_help(alignment,
                  working_dir,
                  taxon_name,
                  gene_names,
                  total_align_file):
    nexus_file = f'{total_align_file}.nex'
    ntax_tot, nchar_tot = ntax_nchar(alignment)

    longest_name = max([len(x.id.replace('.', '').replace('-', '_')) for x in alignment])
    pre_nexus = f'''#NEXUS
    
BEGIN DATA;
    	DIMENSIONS  NTAX={ntax_tot} NCHAR={nchar_tot};
    	FORMAT DATATYPE = DNA GAP = - MISSING = ?;
    	MATRIX
'''
    post_nexus = f'''
;
    
END;'''
    with open(nexus_file, 'w') as file:
        file.write(pre_nexus)
        for record in alignment:
            header = '\t' + record.id.replace('.', '').replace('-', '_')
            num_of_blanks = longest_name - len(header) + 2
            blank_string = ' ' * num_of_blanks
            file.write(f'''{header}{blank_string}{str(record.seq).upper()}
''')
        file.write(post_nexus)
    return nexus_file


## MRBAYES postfix ##
def to_mrbayes_nexus_help(alignments,
                          gene_names,
                          working_dir,
                          taxon_name,
                          date_name,
                          total_align_file,
                          mrb_generations,
                          save_each,
                          burnin,
                          jmod):
    from shutil import copyfile
    partition_str = ''
    start = 1
    end = 0
    old = 0
    part_gen = str(gene_names).replace("\'", "").replace("]", "").replace("[", "")
    part_count = str([x + 1 for x in range(len(gene_names))]).replace("\'", "").replace("]", "").replace("[", "")
    partition_genes = f'\tpartition genes = {len(gene_names)}: {part_gen};\n\n'
    partition_genes += f'\tset partition = genes;\n\n'
    partition_genes += f'\tset applyto=(all);'
    priors = ''
    if len(gene_names) > 1:
        for i, gene_name in enumerate(gene_names):
            _, nchar = ntax_nchar(alignments.loc[gene_name])
            start += old
            end += nchar
            old = nchar
            partition_str += f'\tcharset {gene_name} = {start} - {end};\n'
            if jmod:
                jmodel_stdout = f'{working_dir}/{date_name}_{taxon_name}_{gene_name}/jmodeltest2_out/stdout.txt'
                prset = jmt_to_mrbayes_priors(jmodel_stdout)
                priors += f'\tprset applyto=({i+1})\n{prset}\n'
    else:
        _, nchar = ntax_nchar(alignments)
        start += old
        end += nchar
        old = nchar
        partition_str += f'\tcharset {gene_name} = {start} - {end};\n'
        if jmod:
            jmodel_stdout = f'{working_dir}/{date_name}_{taxon_name}_{gene_name}/jmodeltest2_out/stdout.txt'
            prset = jmt_to_mrbayes_priors(jmodel_stdout)
            priors = f'\tprset applyto=({i+1})\n{prset}\n'
    partition_str += f'\n{partition_genes}'
    nex_file = f'{total_align_file}.nex'
    mrb_nex_file = f'{total_align_file}_mrb.nex'

    mrb_credits = f'''
begin mrbayes;
    
\tlog start filename=Combined_run_{taxon_name}.log;
    
{partition_str}

{priors}
\tlset nst=6 rates=invgamma;
\tunlink statefreq=(all) revmat=(all) shape=(all) pinvar=(all);
    
\tprset applyto=(all) ratepr=variable;
    	
\tmcmc ngen={mrb_generations} nruns=2 printfreq={save_each} savebrlens=yes;
    
\tsump burnin={burnin};
\tsumt burnin={burnin};
    
quit;
end;'''

    copyfile(nex_file, mrb_nex_file)

    with open(mrb_nex_file, 'a') as file:
        file.write(mrb_credits)
    return mrb_nex_file


def to_mrbayes_nexus(alignments,
                     gene_names,
                     working_dir,
                     taxon_name,
                     mrb_generations=50000000,
                     save_each=1000,
                     burnin=25000,
                     jmod=False,
                     concat=False):
    total_gene_name = '_'.join(gene_names) + '_'
    if concat:
        total_align_file = f'{working_dir}/{taxon_name}_{total_gene_name}mafft'
        mrb_nex_file = to_mrbayes_nexus_help(alignments,
                                             gene_names,
                                             working_dir,
                                             taxon_name,
                                             date_name,
                                             total_align_file,
                                             mrb_generations,
                                             save_each,
                                             burnin,
                                             jmod)
    else:
        mrb_nex_file = pd.Series()
        for gene_name in gene_names:
            total_align_file = f'{working_dir}/{date_name}_{taxon_name}_{gene_name}/{taxon_name}_{gene_name}_mafft'
            mrb_nex_file.loc[gene_name] = to_mrbayes_nexus_help(alignments,
                                                                gene_name,
                                                                working_dir,
                                                                taxon_name,
                                                                date_name,
                                                                total_align_file,
                                                                mrb_generations,
                                                                save_each,
                                                                burnin,
                                                                jmod)
    return mrb_nex_file


#### MrBayes on CIPRES ####
def cipro_bayes_help(client,
                     gene_names,
                     working_dir,
                     taxon_name,
                     total_align_file,
                     runtime,
                     concat):
    mrb_nex_file = f'{total_align_file}_mrb.nex'
    global verbose
    properties = cip.Application().getProperties()
    if properties.VERBOSE:
        verbose = True

    mrbob = client.submitJob({"toolId": "MRBAYES_XSEDE", "runtime_": runtime, "mrbayesblockquery_": 1},
                             {'infile_': mrb_nex_file},
                             {'statusEmail': 'true'}, validateOnly=False)
    cipres_lyrics(mrbob,
                  'mrbayes',
                  gene_names,
                  concat,
                  '',
                  mrb_nex_file)
    mrbob.show(messages=True)
    mrbob.update()
    return mrbob


def cipro_bayes(client,
                gene_names,
                working_dir,
                taxon_name,
                runtime=168,
                concat=True):
    if concat:
        total_gene_name = '_'.join(gene_names) + '_'
        total_align_file = f'{working_dir}/{taxon_name}_{total_gene_name}mafft'
        mrbob = cipro_bayes_help(client,
                                 gene_names,
                                 working_dir,
                                 taxon_name,
                                 total_align_file,
                                 runtime=runtime,
                                 concat)
    else:
        mrbob = pd.DataFrame()
        for gene_name in gene_names:
            total_align_file = f'{working_dir}/{date_name}_{taxon_name}_{gene_name}/{taxon_name}_{gene_name}_mafft'
            mrbob[gene_name] = cipro_bayes_help(client,
                                                gene_names,
                                                working_dir,
                                                taxon_name,
                                                total_align_file,
                                                runtime=runtime,
                                                concat)
    return mrbob


# mrbob.delete()
# mrbob.show(messages=True)
# mrbob.waitForCompletion()

## If you lost your job (retrieves last job)
def retreive_last_cipro(client):
    return client.listJobs()[-1]


## Delete a job:
def delete_cipro(client,
                 job):
    job.delete()


## Delete a job:
def queue_cipro(client,
                job):
    job.update()
    job.show(messages=True)
    job.waitForCompletion()


## Download MrBayes Results ##
def dl_cipro_bayes(client,
                   gene_names,
                   working_dir,
                   taxon_name):
    total_gene_name = '_'.join(gene_names) + '_'
    treepath = f'{working_dir}/{taxon_name}_{total_gene_name}_mrbayes'
    if not os.path.exists(treepath):
        os.makedirs(treepath)
        print(f'>> Directory {treepath} was created.')
    mrbob.downloadResults(directory=treepath)
    return f'> MrBayes results were downloaded to {treepath}.'


## To Phylip
def to_phylip_help(alignment,
                   working_dir,
                   taxon_name,
                   gene_names,
                   date_name,
                   ncbi_database,
                   total_align_file):
    partofile = f'{total_align_file}_partitions.txt'
    phylip_file = f'{total_align_file}.phy'
    start = 1
    end = 0
    old = 0
    partition_str = ''
    if ncbi_database[0].lower() == 'nuccore':
        dcode = 'DNA'
    else:
        dcode = 'AA'

    if len(gene_names) > 1:
        for i, gene_name in enumerate(gene_names):
            _, nchar = ntax_nchar(alignment.loc[gene_name])
            start += old
            end += nchar
            old = nchar
            partition_str += f'{dcode}, {gene_name}={start}-{end}\n'
    else:
        _, nchar = ntax_nchar(alignment)
        start += old
        end += nchar
        old = nchar
        partition_str += f'{dcode}, {gene_name}={start}-{end}\n'

    with open(partofile, 'w') as file:
        file.write(partition_str)
    # Phylip file
    totalignment = concaterpillar(alignment,
                                  gene_names,
                                  working_dir,
                                  date_name,
                                  taxon_name,
                                  ncbi_database)
    longest_name = max([len(x.id.replace('.', '').replace('-', '_')) for x in totalignment])
    ntax, nchar = ntax_nchar(alignment.loc[gene_name])
    pre_phylip = f''' {ntax} {nchar}
'''
    with open(phylip_file, 'w') as file:
        file.write(pre_phylip)
        for record in totalignment:
            header = record.id.replace('.', '').replace('-', '_')
            num_of_blanks = longest_name - len(header) + 1
            blank_string = ' ' * num_of_blanks
            file.write(f'''{header}{blank_string}{str(record.seq).upper()}
''')
    return phylip_file


def to_phylip(alignment,
              working_dir,
              taxon_name,
              gene_names,
              date_name,
              ncbi_database,
              concat=True):
    if concat:
        total_gene_name = '_'.join(gene_names) + '_'
        total_align_file = f'{working_dir}/{taxon_name}_{total_gene_name}mafft'
        phylip_file = to_phylip_help(alignment,
                                     working_dir,
                                     taxon_name,
                                     gene_names,
                                     date_name,
                                     ncbi_database,
                                     total_align_file)
    else:
        phylip_file = pd.Dataframe
        for gene_name in gene_names:
            total_align_file = f'{working_dir}/{date_name}_{taxon_name}_{gene_name}/{taxon_name}_{gene_name}mafft'
            phylip_file.loc[gene_name] = to_phylip_help(alignment,
                                                        working_dir,
                                                        taxon_name,
                                                        gene_names,
                                                        date_name,
                                                        ncbi_database,
                                                        total_align_file)
    return phylip_file


## RaxML
def cipro_rax_help(client,
                   gene_names,
                   working_dir,
                   taxon_name,
                   ncbi_database,
                   total_align_file,
                   runtime,
                   concat):
    total_gene_name = '_'.join(gene_names) + '_'
    phylip_file = f'{total_align_file}.phy'
    partofile = f'{total_align_file}_partitions.txt'
    global verbose
    properties = cip.Application().getProperties()
    if properties.VERBOSE:
        verbose = True
    if ncbi_database[0].lower() == 'nuccore':
        dcode = 'dna'
    else:
        dcode = 'protein'
    # Submit the job
    raxbob = client.submitJob({"toolId": "RAXMLHPC8_REST_XSEDE",
                               "runtime_": runtime,
                               "datatype_": dcode},
                              {'input.infile_': phylip_file,
                               'input.partition_': partofile},
                              {'statusEmail': 'true'}, validateOnly=False)
    cipres_lyrics(raxbob,
                  'raxml',
                  gene_names,
                  concat,
                  dcode,
                  f'infile: {phylip_file}; partfile: {partofile}')
    raxbob.show(messages=True)
    raxbob.update()
    return raxbob


def cipro_rax(client,
              gene_names,
              working_dir,
              taxon_name,
              date_name,
              ncbi_database,
              runtime=48,
              concat=True):
    if concat:
        total_gene_name = '_'.join(gene_names) + '_'
        total_align_file = f'{working_dir}/{taxon_name}_{total_gene_name}mafft'
        raxjob_file = cipro_rax_help(client,
                                     gene_names,
                                     working_dir,
                                     taxon_name,
                                     ncbi_database,
                                     total_align_file,
                                     runtime,
                                     concat)
    else:
        raxjob_file = pd.Dataframe
        for gene_name in gene_names:
            total_align_file = f'{working_dir}/{date_name}_{taxon_name}_{total_gene_name}/{taxon_name}_{gene_name}mafft'
            raxjob_file.loc[gene_name] = cipro_rax_help(client,
                                                        gene_names,
                                                        working_dir,
                                                        taxon_name,
                                                        ncbi_database,
                                                        total_align_file,
                                                        runtime,
                                                        concat)
    return raxjob_file


def cipres_lyrics(job,
                  ciprosoft,
                  gene_name,
                  concat,
                  additional_information,
                  file_name):
    HEADER = ['date_submitted',
              'software',
              'job_handle',
              'job_object',
              'gene_name',
              'concatenated_set',
              'additional_information',
              'intput_file']
    job.update()
    cipres_memory = pd.DataFrame(data=[str(job.dateSubmitted),
                                       ciprosoft,
                                       str(job.jobHandle),
                                       [job],
                                       str(gene_name),
                                       concat,
                                       additional_information,
                                       file_name],
                                 index=HEADER)
    csv_filepath = f'{working_dir}/CIPRES_memory.csv'
    try:
        cipres_csv = pd.read_csv(csv_filepath,
                                 sep=',',
                                 header='infer')
        cipres_memory = pd.concat([cipres_csv, cipres_memory])
    except FileNotFoundError:
        pass

    cipres_memory.to_csv(csv_filepath)
    print(f'CIPRES run was documented in {csv_filepath}.')
    return


#### JModelTest 2 on CIPRES ####
def cipro_jmod(client,
               gene_names,
               working_dir,
               date_name,
               taxon_name,
               concat_input=False,
               information_criterion='AIC',
               runtime=5):
    print(f'------------------------- {" - ".join(gene_names)} -------------------------')
    COLS = ['date_released', 'job_name', 'handle_for_use', 'gene', 'information_criterion']
    cipres_memory = pd.DataFrame(columns=COLS)
    if concat_input:
        print('One jModeltest will be run for the concatenated alignment. \
Not recommended -> indepentent evolution of different loci.')
        total_gene_name = '_'.join(gene_names) + '_'
        nex_file = f'{working_dir}/{taxon_name}_{total_gene_name}mafft.nex'
        i = len(cipres_memory)
        job = client.submitJob(
            [
                ("toolId", "JMODELTEST2_XSEDE"),
                ("runtime_", runtime),
                ("criteria_1_", f'-{information_criterion}')
            ],
            {"infile_": nex_file},
            {"statusEmail": "true"}, validateOnly=False)
        job.update()
        cipres_memory.loc[1] = [str(job.dateSubmitted),
                                str(job.jobHandle),
                                [job],  # str(job)
                                str(gene_name),
                                information_criterion]
    else:
        print('Subsequent jModeltests will be run for each of the single alignments. \
Recommended -> indepentent evolution of different loci.')
        for gene_name in gene_names:
            gene_dir = f'{working_dir}/{date_name}_{taxon_name}_{gene_name}'
            nex_file = f'{gene_dir}/{taxon_name}_{gene_name}_mafft.nex'
            i = len(cipres_memory)
            job = client.submitJob(
                [
                    ("toolId", "JMODELTEST2_XSEDE"),
                    ("runtime_", "5"),
                    ("criteria_1_", f'-{information_criterion}')
                ],
                {"infile_": nex_file},
                {"statusEmail": "true"}, validateOnly=False)
            cipres_lyrics(job,
                          'jmodeltest',
                          gene_name,
                          concat_input,
                          information_criterion,
                          nex_file)
    print(f'The jModeltest analysis for the genes {gene_names} was started. A summary \
can be found in {working_dir}/CIPRES_memory.csv. You can refer to that \
list if you lose the overview.')
    return cipres_memory


## Download the jModeltest results
def dlcipro_jmod(client,
                 gene_names,
                 working_dir,
                 date_name,
                 taxon_name):
    for gene_name in gene_names:
        jmod_dir = f'{working_dir}/{date_name}_{taxon_name}_{gene_name}/jmodeltest2_out'
        cipres_memory = pd.read_csv(f'{working_dir}/CIPRES_stat_jmodeltest2.csv',
                                    sep=',',
                                    header='infer')
        handle = cipres_memory.loc[cipres_memory.loc[:, 'gene'] == gene_name, 'job_name'].iloc[0]
        job = client.getJobStatus(handle)
        if job.terminalStage:
            if not os.path.exists(jmod_dir):
                os.makedirs(jmod_dir)
                print(f'>> Directory {jmod_dir} was created.')
            job.downloadResults(directory=jmod_dir)
            return f'The jModeltest results for the genes {gene_names} were downloaded \
to the specific directory.'
        else:
            return f'jModeltest analysis is not yet finished.'


## Use the jModeltest 2 results for mrbayes block
def jmt_to_mrbayes_priors(jmodel_stdout):
    file = open(jmodel_stdout, mode='r')
    lines = file.readlines()
    i, k = [False] * 2
    statefreqpr = '\tstatefreqpr=fixed('
    for j, line in enumerate(lines):
        if 'INFORMATION CRITERION (' in line:
            i = True
            y = 0
            p = 0
        if i:
            if y > 8:
                numbers = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, '.']
                num = ''.join([x for x in line if x in '0123456789.'])
                statefreqpr += num
                p += 1
                if p == 4:
                    i = False
                    k = True
                    statefreqpr += ')\n\t'
                    statefreqpr += 'prset revmatpr=fixed('
                else:
                    statefreqpr += ','
            y += 1
        if k:
            if y > 13:
                num = ''.join([x for x in line if x in '0123456789.'])
                p += 1
                if p == 10:
                    statefreqpr += num
                elif p == 13:
                    statefreqpr += ');\n'
                    k = False
                elif p > 10:
                    if 'p-inv' in line:
                        statefreqpr += ')\n\tpinvarpr=fixed(' + num
                    if 'gamma shape =' in line:
                        statefreqpr += ')\n\tshapepr=exponential(' + num
                else:
                    statefreqpr += num + ','
            y += 1
    return statefreqpr


#### @ future me: try this sometimes!!!! Function of function
def concat_or_not(alignment,
                  working_dir,
                  taxon_name,
                  date_name,
                  gene_names,
                  function,
                  concat=True):
    if concat:
        total_gene_name = '_'.join(gene_names) + '_'
        total_align_file = f'{working_dir}/{taxon_name}_{total_gene_name}mafft'
        output_file = function
    else:
        output_file = pd.Dataframe
        for gene_name in gene_names:
            total_align_file = f'{working_dir}/{date_name}_{taxon_name}_{total_gene_name}/{taxon_name}_{gene_name}mafft'
            output_file.loc[gene_name] = function
    return output_file


def foo(some_func, gene_names, x, y, z, concat):
    if concat:
        result = some_func(gene_names, x, y, z)
    else:
        result = pd.DataFrame()
        for gene_name in gene_names:
            result[gene_name] = some_func(x, y, z)
    return result



# client.getJobStatus('key')
def unfinished_jobs(client):
    COLS = ['date_released', 'job_name', 'handle_for_use']
    cipres_unfinished = pd.DataFrame(columns=COLS)
    for i, job in enumerate(client.listJobs()):
        if not job.terminalStage:
            cipres_unfinished.loc[i] = [str(job.dateSubmitted),
                                        str(job.jobHandle),
                                        str(job)]
    return cipres_unfinished

#### If you lost yout job (retrieves last job): ####
# mrbob = client.listJobs()[-1]
# mrbob = client.getJobStatus('NGBW-JOB-MRBAYES_XSEDE-370B383830DA4A8D98233EF1585B5FB9')
