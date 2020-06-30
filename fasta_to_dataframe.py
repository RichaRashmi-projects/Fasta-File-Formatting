# -*- coding: utf-8 -*-
"""
Created on Tue Jun  9 16:36:31 2020

@author: Richa Rashmi
"""
 ## Please instal the packages biopython and pandas beforehand
 ## command is "pip install biopython" & "pip install pandas"
 
from Bio import SeqIO
import pandas as pd

fasta_dict = dict()

## Enter path to the fasta file here
fasta_path = "F:/scripts/PF05485_full_length_sequences.fasta" 

## Enter path to the csv file here
table_path = "F:/scripts/PF05485_full_length_sequences.csv" 

## Covert fasta file into dictionary
for record in SeqIO.parse(fasta_path, 'fasta'):
    #print('>{}\t{}'.format(record.description, record.seq))
    fasta_dict[record.description] = record.seq

## this command creates table from fasta dictionary
fasta_df = pd.DataFrame(fasta_dict.items())
fasta_df.columns = ["header", "sequence"]

## Formating header
## Protein Id, organism name, protein name and protein sequence format
fasta_df['header'] = fasta_df.header.str.replace(r"\(.*\)","")


fasta_df[['Protein_ID', 'Organism']] = fasta_df.header.str.split('_',expand=True)
## for removing duplicated based on Pritein ID column
df = fasta_df[['Protein_ID', 'Organism', 'sequence']].drop_duplicates(subset=['Protein_ID'])

# Create a column with sequence length
df['Sequence_Length'] = df['sequence'].str.len()


## save the tabular file
df.to_csv(table_path)

##############################################################################

from Bio import SeqIO
import pandas as pd

fasta_dict = dict()

## Enter path to the fasta file here
fasta_path = "F:/scripts/sequence.fasta" 

## Enter path to the csv file here
table_path = "F:/scripts/sequence.csv" 

## Covert fasta file into dictionary
for record in SeqIO.parse(fasta_path, 'fasta'):
    #print('>{}\t{}'.format(record.description, record.seq))
    fasta_dict[record.description] = record.seq

## ths command creates table from fasta dictionary
fasta_df = pd.DataFrame(fasta_dict.items())
fasta_df.columns = ["header", "sequence"]

## excluding isoforms, uncharacterized, partial, predicted and hypothetical proteins from the dataframe
exclude = ['hypothetical',
           'Hypothetical',
           'predicted',
           'uncharacterized',
           'isoform',
           'PREDICTED',
           'Uncharacterized',
           'partial']

fasta_df = fasta_df[~fasta_df['header'].str.contains('|'.join(exclude))]

## Formating header
## Protein Id, organism name, protein name and protein sequence format

fasta_df[['temp1','temp2','temp3','temp4','temp5' ]] = fasta_df.header.str.split('|',expand=True)

fasta_df[['info','organism']] = fasta_df.temp5.str.split('\[|\]', expand=True).iloc[:,[0,1]]

## Selecting and renaming the columns needed
df = (fasta_df[['temp4','organism', 'info', 'sequence']]
      .rename(columns={"temp4": "Protein_ID", "info": "Protein_Name"}))

## creating column with sequence length
df['Sequence_Length'] = df['sequence'].str.len()


## save the tabular file
df.to_csv(table_path)
