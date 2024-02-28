#!/usr/bin/env python
# coding: utf-8

# **Obtaining Files for Analysis**

# In[194]:


file_names = [
    "aro_categories.tsv",
    "aro_categories_index.tsv",
    "aro_index.tsv",
    "CARD-Download-README.txt",
    "card.json",
    "nucleotide_fasta_protein_homolog_model.fasta",
    "nucleotide_fasta_protein_knockout_model.fasta",
    "nucleotide_fasta_protein_overexpression_model.fasta",
    "nucleotide_fasta_protein_variant_model.fasta",
    "nucleotide_fasta_rRNA_gene_variant_model.fasta",
    "protein_fasta_protein_homolog_model.fasta",
    "protein_fasta_protein_knockout_model.fasta",
    "protein_fasta_protein_overexpression_model.fasta",
    "protein_fasta_protein_variant_model.fasta",
    "shortname_antibiotics.tsv",
    "shortname_pathogens.tsv",
    "snps.txt",
]


# In[195]:


def extract_fasta_files(file_names):
    # Create a dictionary to store FASTA files based on their suffix
    fasta_files_dict = {}
    
    for file in file_names:
        # Extract the suffix of the file name
        suffix = file.split('.')[-1]
        
        # Check if the suffix is 'fasta'
        if suffix == 'fasta':
            # Check if the suffix is already a key in the dictionary
            if suffix in fasta_files_dict:
                fasta_files_dict[suffix].append(file)
            else:
                # If not, create a new key and store the file name in a list
                fasta_files_dict[suffix] = [file]
    
    return fasta_files_dict

result = extract_fasta_files(file_names)
print(result)


# **Analysis of CARD using Biopython**

# In[196]:


pip install biopython


# As described in Biopython: The combination of homolog fasta files below uses the Bio.SeqIO.parse(...) function, which is used to loop over all the records in a file as SeqRecord objects. TheS eqRecord (SequenceRecord) class is defined in the Bio.SeqRecordmodule. This class allows higher level features such as identifiers and features to be associated with a sequence (see Chapter3), and is the basic data type for the Bio.SeqIO sequence input/output interface (seeChapter5). https://biopython.org/DIST/docs/tutorial/Tutorial.pdf

# In[224]:


# Obtaining tools
from Bio import SeqIO
import pandas as pd

# Defining a funciton to convert fasta files into a dataframe
def fasta_to_dataframe(fasta_file):
    records = SeqIO.parse(fasta_file, "fasta")
    data = []
    for record in records:
        parts = record.id.split('|')
        aro_id = parts[-1].split(':')[0]
        other_identifiers = '|'.join(parts[:-1])
        data.append((aro_id, other_identifiers, str(record.seq)))
    df = pd.DataFrame(data, columns=['ARO', 'Other_identifiers', 'Sequence'])
    return df

# Defining a function to merge converted fasta dataframe
def merge_fasta_files(file1, file2):
    df1 = fasta_to_dataframe(file1)
    df2 = fasta_to_dataframe(file2)
    merged_df = pd.merge(df1, df2, on='ARO', how='outer')
    
    # Drop NaN values not needed since dataframe contains homolog (complete) data --> useful for other forms of Fasta files. 
   # merged_df.dropna(inplace=True)
    return merged_df

# File paths for homolog Fasta files
file1 = "nucleotide_fasta_protein_homolog_model.fasta"
file2 = "protein_fasta_protein_homolog_model.fasta"

# Merge the files
merged_df = merge_fasta_files(file1, file2)


# In[225]:


# Organize dataframe and rename the columns
merged_df.rename(columns={'Other_identifiers_x': 'Other_identifiers_nucleotide',
                          'Other_identifiers_y': 'Other_identifiers_protein',
                          'Sequence_x': 'Sequence_nucleotide',
                          'Sequence_y': 'Sequence_protein'}, inplace=True)

# Renaming the ARO term to match CARD literature, see the 
dict = {'ARO':'CARD_Short_Name'}
merged_df.rename(columns=dict,
                 inplace=True)


# In[226]:


merged_df.head(5)


# In[228]:


# Using Python documentation to develop pattern used to extract specific identifiers for proteins and nucleotides: https://docs.python.org/3/library/re.html and https://www.regular-expressions.info/

# Extraction pattern for nucleotide identifiers
pattern_nucleotide = r'(\w+)\|(\w+)\.(\d+)\|([+-])\|(\d+-\d+)'

# Extraction pattern for protein identifiers
pattern_protein = r'(\w+)\|(\w+)\.(\d+)\|ARO:(\d+)'

# Extracting information from the "Other_identifiers_nucleotide" column
merged_df[['Database_nucleotide', 'Accession_Number_Nucleotide', 'Sequence_Version_Nucleotide', 'Strand_Orientation_Nucleotide', 'Sequence_Position_Nucleotide']] = merged_df['Other_identifiers_nucleotide'].str.extract(pattern_nucleotide)

# Extracting information from the "Other_identifiers_protein" column
merged_df[['Database_protein', 'Accession_Number_protein', 'Sequence_Version_protein', 'ARO']] = merged_df['Other_identifiers_protein'].str.extract(pattern_protein)
merged_df.head()


# In[229]:


# Dropping the Other_identifiers column, the data within these columns is already present in the dataframe
merged_df.drop(['Other_identifiers_nucleotide','Other_identifiers_protein'], axis = 1, inplace = True)
merged_df.head()


# **Reviewing Content** using CblA-1 (row 0) as an example
# 
# gb: This indicates the type of identifier. In this case, it refers to GenBank, a widely used nucleotide sequence database.
# 
# GQ343019.1: This is the accession number, a unique identifier assigned to a sequence record in the GenBank database. It allows users to quickly locate and retrieve specific sequences.
# 
# +: This symbol indicates the strand orientation of the sequence. In this case, the plus sign (+) typically denotes the forward or sense strand.
# 
# 132-1023: This represents the position of the sequence within the genome or the length of the sequence. In this example, it likely indicates that the sequence starts at position 132 and ends at position 1023.
# 
# ARO:3002999: This part of the identifier is specific to the Antibiotic Resistance Ontology (ARO). It provides information about antibiotic resistance genes or mechanisms associated with the sequence. 
# _______________________________________________________________________________________________________________________________________________________________________________________________________________________
# Investigation of the second term in the identifier "gb|GQ343019.1|+|132-1023|ARO:3002999" is "GQ343019.1". This part refers to the accession number of the sequence. Here's what each component of this term typically signifies:
# 
# GQ: This part of the accession number usually represents the GenBank division where the sequence is stored. In this case, "GQ" could indicate a specific category or source of the sequence within the GenBank database.
# 
# 343019: This is a unique numerical identifier assigned to the sequence within its division or category. It's used to distinguish this sequence from others within the same division.
# 
# .1: The decimal portion of the accession number often indicates the version of the sequence. When a sequence undergoes updates or revisions, a new version is assigned to it. The version number allows users to track changes to the sequence over time.
# _______________________________________________________________________________________________________________________________________________________________________________________________________________________
# Information regarding pattern used for extraction of identifiers from columns:
# 
# (\w+): Matches one or more word characters (alphanumeric or underscore), capturing the GenBank division.
# 
# \|: Matches the pipe character.
# 
# (\w+)\.(\d+): Matches the alphanumeric characters followed by a dot and then more digits, capturing the numerical identifier and version.
# 
# \|: Matches the pipe character.
# 
# ([+-]): Matches either a plus or minus sign, capturing the strand information.
# 
# \|: Matches the pipe character.
# 
# (\d+-\d+): Matches digits followed by a hyphen and then more digits, capturing the position range.
# 
# \|ARO:(\d+): Matches "|ARO:" followed by digits, capturing the ARO identifier.

# **Investigating the Dataframe**

# In[230]:


#determine the type
merged_df.dtypes


# In[231]:


merged_df.describe


# In[232]:


merged_df["CARD_Short_Name"].value_counts() 


# The value counts demonstrates how there is one row deticated to each CARD_Short_Name term

# In[233]:


from Bio.Data import CodonTable

# Access the standard/universal codon table
standard_table = CodonTable.unambiguous_dna_by_name["Standard"]

print(standard_table)


# This table shows the three letter combinations of nucleotides that make up protein sequences. The percentage of amino acids and their orientation within a protein denote some to the proteins characteristics.
# **Can be used as a figure to demonstrate conversion of nucleotides to amino acid sequence**

# In[240]:


def calculate_lengths(row):
    lengths = {}
    lengths['Sequence_nucleotide'] = len(row['Sequence_nucleotide'])
    lengths['Sequence_protein'] = len(row['Sequence_protein'])
    return lengths

# Iterating through each row
for index, row in merged_df.iterrows():
    lengths = calculate_lengths(row)
    print(f"CARD_Short_Name: {row['CARD_Short_Name']}")
    for key, value in lengths.items():
        print(f"{key}: {value}")
    print()


# The output shows the protein and amino acid sequence length for each row

# In[ ]:




