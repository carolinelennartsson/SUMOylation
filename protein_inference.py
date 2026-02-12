import os
import pandas as pd
import numpy as np
import re

import biotite.database.entrez as entrez
import biotite.sequence as seq
from biotite.sequence import ProteinSequence
import biotite.sequence.io.fasta as fasta


def find_uni(df, annos):
    """
    Canonicalize UniProt entries for each row in a dataframe based on UniProt annotation priority.

    For each row in `df`, this function:
    - Extracts a list of candidate UniProt protein ids from the
      `Proteins_list` column.
    - Filters the annotation dataframe (`annos`) to only those entries.
    - Sorts the filtered annotations by a defined priority:
        1. Reviewed status (reviewed entries preferred)
        2. Total Count (higher preferred)
        3. Annotation score (higher preferred)
        4. UniProt Entry ID (alphabetical)
        5. Entry name length (longer preferred)
    - Selects the top-ranked UniProt entry as the canonical protein.
    - Populates UniProt, gene, Gene Ontology, and keyword fields in `df`
      using the selected annotation.

    Parameters
    ----------
    df : pandas.DataFrame
        Input dataframe containing a `Proteins_list` column with UniProt
        accession IDs per row. The dataframe is modified in place.

    annos : pandas.DataFrame
        Annotation dataframe containing UniProt metadata. Must include
        the following columns:
        - Entry
        - Reviewed
        - Total Count
        - Annotation
        - Entry name len
        - Gene Names (primary)
        - Gene Ontology (biological process)
        - Gene Ontology (cellular component)
        - Gene Ontology (molecular function)
        - Keywords

    Returns
    -------
    pandas.DataFrame
        The input dataframe with additional columns populated:
        - Uniprot_cannonised
        - Gene_name
        - Gene Ontology (biological process)
        - Gene Ontology (cellular component)
        - Gene Ontology (molecular function)
        - Keywords

    Notes
    -----
    - Rows with no matching annotations in `annos` are left unchanged.
    - If multiple candidate proteins exist, only the highest-ranked
      annotation is retained.
    """
    
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
            df.at[index, 'Keywords'] = sorted_df['Keywords'][0] 
    return df


def extract_following_symbols(string, list_of_strings):
    """
    Extract symbol pairs following a target substring within a list of strings.

    This function searches each element in `list_of_strings` for occurrences
    of `string`. When found, it extracts the two characters immediately
    following the matched substring. Only symbol pairs that begin with a
    hyphen ('-') are retained.

    Parameters
    ----------
    string : str
        Substring to search for within each item of `list_of_strings`.

    list_of_strings : list of str
        List of strings to be searched.

    Returns
    -------
    list of str
        A list of two-character substrings that directly follow `string`
        and start with '-'.

    Notes
    -----
    - Non-string inputs return an empty list.
    - If `string` is not found in an item, that item is ignored.
    - Only the first occurrence per item is considered.
    """

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
    """
    Combine a base string with each element of a list of strings.

    This function prepends `string` to every element in `list_of_strings`
    and also includes the original `string` as the first element of the
    returned list.

    Parameters
    ----------
    string : str
        Base string to prepend. If not a string, it is replaced with
        an empty string.

    list_of_strings : list of str
        List of strings to be combined with `string`.

    Returns
    -------
    list of str
        A list where the first element is `string`, followed by
        `string + item` for each item in `list_of_strings`.


    """

    if not isinstance(string, str):  # Ensure the string is valid
        string = ""
    if not isinstance(list_of_strings, list):  # Ensure the list is valid
        return []
    
    combined = [string + item for item in list_of_strings]
    combined.insert(0, string)  
    return combined


# used to controlcheck the result  
def match_any_col2(list1, list2):
    """
    Identify elements in a list that contain any substring from another list.

    This function checks whether any element in `list1` is a substring
    of each element in `list2`. It returns both a boolean match mask and
    the list of matching elements from `list2`.

    """

    match_results = [any(item1 in item2 for item1 in list1) for item2 in list2]
    matched_peptides = [item2 for item2, match in zip(list2, match_results) if match]

    return match_results, matched_peptides


    

def count_go_and_keywords(row):
    """
    Count Gene Ontology terms and keywords for a dataframe row.

    This function calculates the total number of Gene Ontology (GO) terms
    and UniProt keywords associated with a single row by counting
    semicolon-separated entries in the relevant columns.

    Parameters
    ----------
    row : pandas.Series
        A row from a pandas DataFrame containing the following fields:
        - 'Gene Ontology IDs' : str or None
        - 'Keywords' : str or None

        Values are expected to be semicolon-separated strings.

    Returns
    -------
    int
        Total count of Gene Ontology terms and keywords combined.

    Notes
    -----
    - Empty or missing values result in a count of 0 for that field.
    - The function assumes ';' as the delimiter for both columns.
    """

    go_terms_count = len(row['Gene Ontology IDs'].split(";")) if row['Gene Ontology IDs'] else 0
    keywords_count = len(row['Keywords'].split(";")) if row['Keywords'] else 0
    return go_terms_count + keywords_count
    

def process_data_msms_MQ(data, species="human"):
    """
    Process MaxQuant msms.txt peptide-level data and annotate canonical UniProt proteins.

    This function is designed specifically to operate on the MaxQuant
    `msms.txt` output file. It maps identified peptide sequences to
    UniProt protein entries using FASTA sequence matching, enriches the
    data with UniProt annotations, selects a canonical protein per peptide,
    and resolves associated isoforms and protein sequences.

    The workflow includes:
    - Loading species-specific UniProt FASTA and annotation files
    - Mapping peptide sequences from msms.txt to UniProt protein sequences
    - Normalizing MaxQuant protein identifiers
    - Selecting a canonical UniProt entry based on annotation priority
    - Extracting isoform identifiers and corresponding protein sequences

    Parameters
    ----------
    data : pandas.DataFrame
        DataFrame created from the MaxQuant `msms.txt` file. Must include
        at least the following columns:
        - 'Sequence' : str
            Peptide sequence identified by MS/MS
        - 'Proteins' : list of str
            UniProt protein accessions assigned by MaxQuant

    species : str, optional
        Species identifier used to select UniProt reference data.
        Supported values:
        - 'human' (default)
        - 'mouse'

    Returns
    -------
    pandas.DataFrame
        The input msms.txt DataFrame enriched with UniProt-based
        annotations and isoform information, including:
        - 'matching_proteins'
        - 'Proteins_list'
        - 'Uniprot_cannonised'
        - 'Gene_name'
        - Gene Ontology annotations
        - 'Keywords'
        - 'isoforms'
        - 'uniprot_w_isoforms'
        - 'isoform_seqs'

    Notes
    -----
    - This function assumes MaxQuant-style formatting and column names
      from `msms.txt`.
    - UniProt FASTA and annotation paths are currently hard-coded.
    - Canonical protein selection relies on `find_uni()`.
    - Peptide-to-protein mapping is based on exact substring matching
      within FASTA sequences.
    - The input DataFrame is modified during processing.

    See Also
    --------
    find_uni
    extract_following_symbols
    combine_string_and_list
    count_go_and_keywords
    """


    def find_matching_keys(peptide):
        return [key for key, seq in string_sequences.items() if peptide in seq]


    if species == 'human': 
        f_path = "../uniprot/HUMAN.fasta"
        annos = "../uniprot/uniprotkb_proteome_UP000005640_2025_03_18.tsv.gz"

    if species == 'mouse': 
        annos = "../uniprot/uniprotkb_proteome_UP000000589_2025_05_07.tsv.gz"
        f_path = "../uniprot/MOUSE.fasta"
  
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
    data['Keywords'] = ''

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
    



def protein_inference_pg(data, species='human'):
    """
    Perform protein inference and canonical UniProt annotation for MaxQuant protein groups.

    This function is designed to process MaxQuant `proteinGroups.txt` output.
    It normalizes protein group identifiers, annotates UniProt proteins using
    reference FASTA and UniProt metadata, and selects a canonical UniProt
    accession per protein group based on annotation priority.

    The workflow includes:
    - Loading species-specific UniProt FASTA and annotation files
    - Normalizing protein identifiers from `Majority protein IDs`
    - Computing Gene Ontology and keyword annotation counts
    - Selecting a canonical UniProt entry per protein group
    - Annotating gene names, Gene Ontology terms, and keywords

    Parameters
    ----------
    data : pandas.DataFrame
        DataFrame created from the MaxQuant `proteinGroups.txt` file.
        Must include at least the following column:
        - 'Majority protein IDs' : str
            Semicolon-separated UniProt protein identifiers assigned
            to the protein group.

    species : str, optional
        Species identifier used to select UniProt reference data.
        Supported values:
        - 'human' (default)
        - 'mouse'

    Returns
    -------
    pandas.DataFrame
        The input proteinGroups DataFrame enriched with UniProt-based
        protein inference results, including:
        - 'Proteins_list'
        - 'Uniprot_cannonised'
        - 'Gene_name'
        - 'Gene Ontology (biological process)'
        - 'Gene Ontology (cellular component)'
        - 'Gene Ontology (molecular function)'
        - 'Keywords'

    Notes
    -----
    - This function assumes MaxQuant-style formatting and column names
      from `proteinGroups.txt`.
    - UniProt FASTA and annotation file paths are currently hard-coded.
    - Canonical protein selection is performed by `find_uni()`.
    - The input DataFrame is modified during processing.

    See Also
    --------
    find_uni
    count_go_and_keywords
    """


    if species == 'human': 
        f_path = "../uniprot/HUMAN.fasta"
        annos = "../uniprot/uniprotkb_proteome_UP000005640_2025_03_18.tsv.gz"

    if species == 'mouse': 
        annos = "../uniprot/uniprotkb_proteome_UP000000589_2025_05_07.tsv.gz"
        f_path = "../CL/uniprot/MOUSE.fasta"
        
    with open(f_path, "r") as file:
        fasta_file = fasta.FastaFile.read(file)
    
    sequence_dict = fasta.get_sequences(fasta_file)
    sequence_dict = {header.split('|')[1]: seq for header, seq in sequence_dict.items()}
    string_sequences = {key: str(seq) for key, seq in sequence_dict.items()}

    uni_annos = pd.read_csv(annos, delimiter="\t")
    uni_annos['Gene Ontology IDs'] = uni_annos['Gene Ontology IDs'].fillna('')
    uni_annos['Total Count'] = uni_annos.apply(count_go_and_keywords, axis=1)
    uni_annos['Entry name len'] = uni_annos['Entry'].str.len()
    
    data['Proteins_list'] = data['Majority protein IDs'].apply(lambda x: x.split(';') if pd.notnull(x) else [])
    data['Proteins_list'] = data['Proteins_list'].apply(lambda sublist: [item.split('|')[1] if '|' in item else item for item in sublist])
    data['Uniprot_cannonised'] = "" 
    data['Gene_name'] = "" 
    data['Gene Ontology (biological process)'] = ''
    data['Gene Ontology (cellular component)'] = ''
    data['Gene Ontology (molecular function)'] = ''
    data['Keywords'] = ''

    data_cannonised = find_uni(data, uni_annos)
    
    return data_cannonised




