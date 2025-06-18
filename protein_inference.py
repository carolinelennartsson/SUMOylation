import os
import pandas as pd
import numpy as np
import re

import biotite.database.entrez as entrez
import biotite.sequence as seq
from biotite.sequence import ProteinSequence
import biotite.sequence.io.fasta as fasta

def find_uni(df, annos): 
    for index, row in df.iterrows(): 
        test_proteins = row['Proteins_list']
        filtered_df = annos[annos['Entry'].isin(test_proteins)]

        sorted_df = filtered_df.sort_values(by=['Reviewed', 'Total Count', 'Annotation', 'Entry', 'Entry name len'], ascending=[True, False, False, True, False])
        sorted_df = sorted_df.reset_index(drop=True)
        if not sorted_df.empty: 
            df.at[index, 'Uniprot_cannonised'] = sorted_df['Entry'][0] 
            df.at[index, 'Gene_name'] = sorted_df['Gene Names (primary)'][0] 
            df.at[index, 'Gene Ontology (biological process)'] = sorted_df['Gene Ontology (biological process)'][0] 
            df.at[index, 'Gene Ontology (cellular component)'] = sorted_df['Gene Ontology (cellular component)'][0] 
            df.at[index, 'Gene Ontology (molecular function)'] = sorted_df['Gene Ontology (molecular function)'][0] 
    return df


# extract isoforms 
def extract_following_symbols(string, list_of_strings):
    if not isinstance(string, str): 
        return []
    
    if not isinstance(list_of_strings, list):  
        return []
    
    result = []
    for item in list_of_strings:
        if isinstance(item, str) and string in item:
            index = item.find(string)
            # Extract the two characters following the match
            following_symbols = item[index + len(string):index + len(string) + 2]
            if following_symbols.startswith('-'):
                result.append(following_symbols)
    return result


def combine_string_and_list(string, list_of_strings):

    if not isinstance(string, str):  # Ensure the string is valid
        string = ""
    if not isinstance(list_of_strings, list):  # Ensure the list is valid
        return []
    
    combined = [string + item for item in list_of_strings]
    combined.insert(0, string)  
    return combined


# used to check the reult 
def match_any_col2(list1, list2):
    match_results = [any(item1 in item2 for item1 in list1) for item2 in list2]
    matched_peptides = [item2 for item2, match in zip(list2, match_results) if match]

    return match_results, matched_peptides


    

def count_go_and_keywords(row):
    go_terms_count = len(row['Gene Ontology IDs'].split(";")) if row['Gene Ontology IDs'] else 0
    keywords_count = len(row['Keywords'].split(";")) if row['Keywords'] else 0
    return go_terms_count + keywords_count
    

def process_data_msms_MQ(data, species = "human"):

    def find_matching_keys(peptide):
        return [key for key, seq in string_sequences.items() if peptide in seq]

    if species == 'human': 
        f_path = "../HUMAN.fasta"
        annos = "../uniprotkb_proteome_UP000005640_2025_03_18.tsv.gz"

    if species == 'mouse': 
        annos = "../uniprotkb_proteome_UP000000589_2025_05_07.tsv.gz"
        f_path = "../MOUSE.fasta"
        
    with open(f_path, "r") as file:
        fasta_file = fasta.FastaFile.read(file)
    
    sequence_dict = fasta.get_sequences(fasta_file)
    sequence_dict = {header.split('|')[1]: seq for header, seq in sequence_dict.items()}
    string_sequences = {key: str(seq) for key, seq in sequence_dict.items()}

    uni_annos = pd.read_csv(annos, delimiter="\t")
    uni_annos['Gene Ontology IDs'] = uni_annos['Gene Ontology IDs'].fillna('')
    
    # Apply the function and add a new column
    uni_annos['Total Count'] = uni_annos.apply(count_go_and_keywords, axis=1)
    uni_annos['Entry name len'] = uni_annos['Entry'].str.len()

    data["matching_proteins"] = data["Sequence"].apply(find_matching_keys)
    
    data['Proteins_list'] = data['Proteins'] #.apply(lambda x: x.split(';') if pd.notnull(x) else [])
    data['Proteins_list'] = data['Proteins_list'].apply(lambda sublist: [item.split('|')[1] if '|' in item else item for item in sublist])
    data['Proteins_list'] = data['Proteins_list'].apply(lambda sublist: [item.split('-')[0] if '-' in item else item for item in sublist])

    data['Uniprot_cannonised'] = "" 
    data['Gene_name'] = ""
    data['Gene Ontology (biological process)'] = ''
    data['Gene Ontology (cellular component)'] = ''
    data['Gene Ontology (molecular function)'] = ''

    data_cannonised = find_uni(data, uni_annos)

    data_cannonised['isoforms'] = data_cannonised.apply(
        lambda row: extract_following_symbols(row['Uniprot_cannonised'], row['Proteins']), axis=1
    )
    
    data_cannonised['uniprot_w_isoforms'] = data_cannonised.apply(
        lambda row: combine_string_and_list(row['Uniprot_cannonised'], row['isoforms']), axis=1
    )

    data_cannonised['isoform_seqs'] = data_cannonised['uniprot_w_isoforms'].apply(
        lambda key_list: [str(sequence_dict.get(key, None)) for key in key_list]
    )
    
    return data_cannonised
    

    


def protein_inference_pg(data, species = 'human'):

    if species == 'human': 
        f_path = "../HUMAN.fasta"
        annos = "../uniprotkb_proteome_UP000005640_2025_03_18.tsv.gz"

    if species == 'mouse': 
        annos = "../uniprotkb_proteome_UP000000589_2025_05_07.tsv.gz"
        f_path = "../MOUSE.fasta"
        
    with open(f_path, "r") as file:
        fasta_file = fasta.FastaFile.read(file)
    
    sequence_dict = fasta.get_sequences(fasta_file)
    sequence_dict = {header.split('|')[1]: seq for header, seq in sequence_dict.items()}
    string_sequences = {key: str(seq) for key, seq in sequence_dict.items()}

    uni_annos = pd.read_csv(annos, delimiter="\t")
    uni_annos['Gene Ontology IDs'] = uni_annos['Gene Ontology IDs'].fillna('')
    uni_annos['Total Count'] = uni_annos.apply(count_go_and_keywords, axis=1)
    uni_annos['Entry name len'] = uni_annos['Entry'].str.len()
    
    data['Proteins_list'] = data['Protein IDs'].apply(lambda x: x.split(';') if pd.notnull(x) else [])
    data['Proteins_list'] = data['Proteins_list'].apply(lambda sublist: [item.split('|')[1] if '|' in item else item for item in sublist])
    data['Uniprot_cannonised'] = "" 
    data['Gene_name'] = "" 
    data['Gene Ontology (biological process)'] = ''
    data['Gene Ontology (cellular component)'] = ''
    data['Gene Ontology (molecular function)'] = ''

    data_cannonised = find_uni(data, uni_annos)
    
    return data_cannonised



