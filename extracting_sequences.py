# function that takes a dataframe
# check if any of the peptides in 'Peptide sequences' matches w the string in 'seq'
# if K is at the end of the 'peprtide sequences', make sure it's followed by a E or D in seq
# takes the positions of the K in the string match in seq 
# takes the ['SUMO2/3-DVFQQQTGG site positions'] # column with list of numbers to see if the numbers match , return the overlapping numbers 
# return the matched seuqences +-7 aroung the site positions, it's in the beginning or end add -- to account for the missing letters

# function that takes a dataframe
# check if any of the peptides in 'Peptide sequences' matches w the string in 'seq'
# if K is at the end of the 'peprtide sequences', make sure it's followed by a E or D in seq
# takes the positions of the K in the string match in seq 
# takes the ['SUMO2/3-DVFQQQTGG site positions'] # column with list of numbers to see if the numbers match , return the overlapping numbers 
# return the matched seuqenceßs +-7 aroung the site positions, it's in the beginning or end add -- to account for the missing letters

# todo, add if k pos 
# add preprosessing steps 

import re
import pandas as pd

# not 0 indexed
def find_k_positions(seq):
    """
    Identify lysine (K) positions in a protein or peptide sequence.

    This function scans a sequence and returns the positions of all
    lysine ('K') residues using 1-based indexing (biological convention).

    Parameters
    ----------
    seq : str
        Amino acid sequence to be scanned.

    Returns
    -------
    list of int
        List of positions (1-based) where the residue is 'K'.

    Notes
    -----
    - Indexing follows biological convention (1-based).
    - The input sequence is converted to a string before processing.
    """
    return [i - 1 for i, letter in enumerate(str(seq)) if letter == 'K']


def extract_sublist(lst, position, number): # , placeholder=None
    """
    Extract a fixed-width window around a position in a list.

    This function returns a sublist centered at a given position,
    extending `number` elements to the left and right. If the window
    extends beyond the list boundaries, a placeholder value is used.

    Parameters
    ----------
    lst : list
        Input list from which elements are extracted.

    position : int
        Central position around which to extract the sublist.
        Expected to follow 0-based indexing.

    number : int
        Number of elements to include on each side of `position`.

    placeholder : any, optional
        Value used when the window exceeds list boundaries.
        Defaults to '-' if not provided.

    Returns
    -------
    list
        Sublist of length `(2 * number + 1)` centered at `position`.

    Notes
    -----
    - If indices fall outside the list bounds, '-' is inserted.
    - The returned list always has a fixed length.
    """

    start = position - number
    end = position + number + 1

    result = []
    for i in range(start, end):
        if 0 <= i < len(lst): 
            result.append(lst[i])
        else: 
            result.append("-")
    return result

    

def extract_matching_sequences(df, seqs_to_use, results_column, pos_column = 'Positions within proteins', res_pdllt_column = 'extracted_pddlt',  AF = False):
    import pandas as pd

    def extract_surrounding_sequence(seq, position, seq_len_tp_extract):
        """Helper function to extract the sequence with padding around the position."""
        start = max(0, position - seq_len_tp_extract)
        end = min(len(seq), position + seq_len_tp_extract + 1)
        substring = seq[start:end]
    
        # Pad with '-' if we're at the beginning or end
        pad_before = max(0, seq_len_tp_extract - position)
        pad_after = max(0, seq_len_tp_extract - (len(seq) - position - 1))
        padded_substring = '-' * pad_before + substring + '-' * pad_after

        return padded_substring

    def process_row(row, seq_len_tp_extract, seqs_to_use, AF):
        """Processes a single row to extract sequences and match peptides."""
        peptide = row['Sequence']
        #print(peptide)
        #print(row[pos_column])
        
        site_positions = row[pos_column]
        seqs = row[seqs_to_use]
        if 'sse' in row and 'plddt' in row:
            row_sse = row['sse']
            row_pdllt = row['plddt']
        else:
            row_sse = None
            row_pdllt = None
            
        pos_peptide = row['Position in peptide']

        matched_sequences = []
        sses = []
        pdllts = []

        for site_pos in site_positions:
            for seq in seqs:
                match_found = False  
                surrounding_sequence = extract_surrounding_sequence(seq, site_pos, seq_len_tp_extract)
                middle_index = seq_len_tp_extract - 1
                
                if surrounding_sequence[middle_index] == 'K':
                    #print(surrounding_sequence[middle_index])
                    #print(surrounding_sequence[middle_index + 1])
                    
                    if pos_column == 'Positions within proteins': 
                        if peptide.endswith('K') and pos_peptide == len(peptide):
                            if surrounding_sequence[middle_index + 1] == 'D' or surrounding_sequence[middle_index + 1] == 'E':
                                if peptide in surrounding_sequence:
                                    #print('------- MATCH ---------')
                                    matched_sequences.append(surrounding_sequence)
                                    match_found = True
                                    #break  
                        else:   
                            if peptide in surrounding_sequence:
                                #print('------- MATCH ---------')
                                matched_sequences.append(surrounding_sequence)
                                match_found = True
                                #break 
                                
                        if AF and match_found:
                            sse = extract_sublist(row_sse, site_pos, seq_len_tp_extract)
                            sses.append(sse)
                            pdllt = extract_sublist(row_pdllt, site_pos, seq_len_tp_extract)
                            #print(pdllt)
                            pdllts.append(pdllt)
                            
                    elif pos_column == 'k_positions': 
                        matched_sequences.append(surrounding_sequence)
                        pdllt = extract_sublist(row_pdllt, site_pos, seq_len_tp_extract)
                        #print(pdllt)
                        pdllts.append(pdllt)
                        match_found = False
                        
                if match_found: 
                    break
        
        #print(matched_sequences)
        return matched_sequences, sses, pdllts 

    seq_len_tp_extract = int(df['Length'].max() + 1)

    df[[results_column, 'extracted_sse', res_pdllt_column]] = df.apply(
    lambda row: pd.Series(process_row(row, seq_len_tp_extract, seqs_to_use, AF)), axis=1
)
    return df

    
#big_df['k_positions'] = big_df['AF_seq'].apply(find_k_positions)

    
# Call the funciton like so to get different output 
#testtt = extract_matching_sequences(df_sumo_unique.head(100), max_length, 'isoform_seqs', 'matched_sequences') 
#result_ttt = extract_matching_sequences(df_sumo_unique.head(1000), max_length, 'AF_seq', 'all_matched_sequences_AF', 'k_positions', 'all_surr_pddlt', AF = True) 
#result_ttt = extract_matching_sequences(df_sumo_unique.head(100), max_length, 'AF_seq', 'matched_sequences_AF', AF = True) 



# extract af info from file 


# af stuff 

#exposure = "/Users/hmt128/Library/CloudStorage/OneDrive-UniversityofCopenhagen/Science/Projects/Ubi_like_modifiers/identification/data/exposure_sites/sites_exposure_w_FDR.txt"
# exposure comes from extra script 
#exposure_info = pd.read_csv(exposure, delimiter=";") # old file
#exposure_info['seq'] = exposure_info['sequence'].str.extract(r'"(.*?)"')

def convert_to_number_list(str_value):
    str_value = str_value.strip('[]').replace('\n', '') #.replace(' ', ',') #.replace('.', '') 
    return [float(x) for x in str_value.split(',') if x.strip() != '']  # Convert only non-empty parts

def convert_to_string_list(str_value):
    str_value = str_value.strip('[]').replace('\n', '') #.replace(' ', ',')
    return [x.strip().strip("'").strip('"') for x in str_value.split(',') if x.strip() != ''] 

#exposure_info['plddt'] = exposure_info['plddt'].apply(convert_to_number_list)
#exposure_info['sse'] = exposure_info['sse'].apply(convert_to_string_list)

#big_df = pd.merge(big_df, exposure_info, left_on='Uniprot_cannonised', right_on='outer_key', how='left') 
#big_df['AF_seq'] = big_df['seq'].fillna('').astype(str)
#big_df['AF_seq'] = big_df['AF_seq'].apply(lambda x: [x]) 





# getting the sequences into a flat list for plotting 

def extract_trimmed_sequences(df, seq_column='matched_sequences', seq_len_tp_extract=7):
    """
    Extract and trim matched sequences from specified datasets in a DataFrame.

    Parameters:
    - df: pandas DataFrame containing the data.
    - dataset_names: list of dataset names to extract sequences from.
    - seq_column: name of the column containing sequences.
    - seq_len_tp_extract: number of amino acids to extract on each side of the center.

    Returns:
    - dict of {dataset_name: list of trimmed sequences}.
    """
    df['matched_sequence'] = df['matched_sequences'].apply(lambda x: x[0] if x else None)
    dataset_names = df['Dataset'].unique()
    
    def trim_sequence(sequence):
        """Trim the sequence to take the part +-7 amino acids around the center."""
        
        middle_index = len(sequence) // 2  # Center of the sequence
        start = max(0, middle_index - seq_len_tp_extract)
        end = min(len(sequence), middle_index + seq_len_tp_extract + 1)
        return sequence[start:end]

    result = {}
    for dataset in dataset_names:
        seqs = df[df['Dataset'] == dataset]['matched_sequence'].dropna().astype(str)
        trimmed_seqs = [trim_sequence(seq) for seq in seqs if isinstance(seq, str)]
        result[dataset] = [s for s in trimmed_seqs if isinstance(s, str) and len(s) > 0]
        
    return result




def calculate_motif_adherence(sequence_list, anywhere=False):
    """
    Calculate the percentage of sequences adhering to a lysine-centered motif.

    This function evaluates a list of peptide or protein sequences for the
    presence of a defined lysine (K) motif and reports the percentage of
    sequences that match the motif.

    Two motif modes are supported:
    - Fixed-position motif (default):
        Six residues upstream of K followed by any residue and E
        (pattern: .{6}K\\wE)
    - Anywhere motif:
        K[su]XE motif anywhere in the sequence
        (pattern: K\\[su\\][A-Z]E)

    Parameters
    ----------
    sequence_list : list of str
        List of sequences to be evaluated. Duplicate sequences are removed
        prior to analysis.

    anywhere : bool, optional
        If False (default), the fixed-position motif is used.
        If True, the motif may occur anywhere in the sequence.

    Returns
    -------
    float
        Percentage of unique sequences that match the selected motif.

    Notes
    -----
    - Sequences are deduplicated using a set before analysis.
    - Motif matching is performed using regular expressions.
    - The function prints the total number of sequences and the number
      of matching sequences as a side effect.
    - The returned value is expressed as a percentage.
    """

    pattern = r'.{6}K\wE'  
    anywhere_pattern = r'K\[su\][A-Z]E'

    sequence_list = set(sequence_list)
    
    if anywhere == False:
        total_sequences = len(sequence_list)
        print(total_sequences)
        matched_sequences = sum(1 for seq in sequence_list if re.search(pattern, seq))
        print(matched_sequences)
    if anywhere == True: 
        total_sequences = len(sequence_list)
        print(total_sequences)
        matched_sequences = sum(1 for seq in sequence_list if re.search(anywhere_pattern, seq))
        print(matched_sequences)
        
    return (matched_sequences / total_sequences) * 100



def find_motif(df, column_name):
    """
    Identify sequences containing a fixed-position lysine-centered motif.

    This function scans sequences in a specified DataFrame column for the
    presence of a lysine (K) motif defined as six residues upstream of K,
    followed by any residue and glutamate (E).

    Motif pattern:
        .{6}K\\wE

    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame containing sequence data.

    column_name : str
        Name of the column in `df` that contains sequences to be analyzed.

    Returns
    -------
    pandas.Series
        Boolean Series indicating whether each sequence contains the motif.
        Entries that are not strings return the string 'None'.

    Notes
    -----
    - Motif detection is performed using regular expressions.
    - Non-string entries are explicitly labeled as 'None' rather than False.
    - The returned Series aligns with the index of the input DataFrame.

    -To run:  df['pattern_match'] = find_motif(df, 'matched_sequence')
    """

    pattern = r'.{6}K\wE'
    return df[column_name].apply(
        lambda seq: bool(re.search(pattern, seq)) if isinstance(seq, str) else 'None'
    )



    
