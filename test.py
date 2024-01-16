# -*- coding: utf-8 -*-
"""
testing the get_gff_file function and make_annotations_table
"""
import pandas as pd
# import csv



def make_annotations_table(file_path):
    attr_dict_list = []
    columns = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
    annotation_df = pd.read_csv(file_path, sep='\t', comment='#', header=None, names=columns)
    for entry in annotation_df['attributes'].str.split(';'):
        attr_dict = {}
        for attr in entry:
            if(attr != '' and '.gtf' in file_path):
                key, value = attr.strip().split(' ')
                attr_dict[key] = value.strip('"')
            elif ('.gff3' in file_path):
                key, value = attr.strip().split('=')
                attr_dict[key] = value.strip()
        attr_dict_list.append(attr_dict)
    attributes_df = pd.DataFrame(attr_dict_list)
    final_dataset = pd.concat([annotation_df,attributes_df],axis=1)
    final_dataset.drop('attributes', axis=1, inplace=True) 
    # print("Number of unique transcripts = ",final_dataset['transcript_id'].nunique()) 
    return final_dataset



def get_gff_file(gff_table, transcript_list = None):
    if(transcript_list != None):
        gff_table = gff_table[gff_table['ID'].str.contains('|'.join(transcript_list), case=False)]
    gff_table = gff_table.reset_index(drop=True)
    final_table = gff_table.loc[:,'seqid':'phase']
    attributes = gff_table.loc[:, 'ID':] 
    attributes = attributes.where(pd.notna(attributes), None)
    attributes = attributes.to_dict()
    attribute_column = []
    for i in range(len(attributes['ID'])):
        column_value_pairs = []
        for column in attributes.keys():
            if (attributes[column][i] is not None):
                column_value = f'{column}={attributes[column][i]}'
                column_value_pairs.append(column_value)
        attribute_column.append(column_value_pairs)
    attribute_column = [';'.join(inner_list) for inner_list in attribute_column]
    final_table['attributes'] = attribute_column
    return final_table

def get_gtf_file (gtf_table):
    gtf_table = gtf_table.reset_index(drop=True)
    final_table = gtf_table.loc[:,'seqid':'phase']    
    attributes = gtf_table.loc[:, 'transcript_id':]
    attributes = attributes.where(pd.notna(attributes), None)
    attributes = attributes.to_dict()
    attribute_column = []
    for i in range(len(attributes['transcript_id'])):
        column_value_pairs = []
        for column in attributes.keys():
            if (attributes[column][i] is not None):
                column_value = f'{column} "{attributes[column][i]}"'
                column_value_pairs.append(column_value)
        attribute_column.append(column_value_pairs)
    attribute_column = [';'.join(inner_list) for inner_list in attribute_column]
    #print(attribute_column[0])
    final_table['attributes'] = attribute_column
    # final_table.to_csv('test.gtf', sep='\t', index=False, header=False, quoting=csv.QUOTE_NONE)    




    
#%%
'''
Testing the feature specific gff extraction
'''

annotation_df = make_annotations_table('final_optimized_dupremoved_sorted.gff3')
gff_table = annotation_df[annotation_df['type'].str.contains('utr',case=False)]
gff_table = gff_table.reset_index(drop=True)
final_table = gff_table.loc[:,'seqid':'phase']
attributes = gff_table.loc[:, 'ID':]
attributes = attributes.where(pd.notna(attributes), None)
attributes = attributes.to_dict()
attribute_column = []
for i in range(len(attributes['ID'])):
    # print(i)
    column_value_pairs = []
    for column in attributes.keys():
        # print(attributes[column][i])
        if (attributes[column][i] is not None):
            column_value = f'{column}={attributes[column][i]}'
            column_value_pairs.append(column_value)
    attribute_column.append(column_value_pairs)
attribute_column = [';'.join(inner_list) for inner_list in attribute_column]
final_table['attributes'] = attribute_column

#%%

'''Testing the feature counts'''

table = make_annotations_table('consensus_noduplicates.gtf')
table = gff_table.drop_duplicates(subset='ID')
pivot_df = table.groupby(['type']).size().reset_index()
pivot_df.rename(columns={0: 'Count'}, inplace=True)
pivot_df.rename(columns={'type': 'Feature'}, inplace=True)
print(pivot_df)

#%%

'''Removing the duplicates from consensus annotation files'''

gtf_file = 'inputs/consensus.gtf'
gff_file = 'inputs/final_optimized_sorted.gff3'

gtf_table = make_annotations_table(gtf_file)
gff_table = make_annotations_table(gff_file)

pattern1 = r'SRAE_1000189300\.([2-9]|10|11|12)'
pattern2 = r'SRAE_1000277600\.[2-7]'
pattern3 = r'SRAE_1000336700\.([2-9]|10|11|12|13|14)'

# import re
# pattern1 = r'SRAE_1000189300\.([2-9]|10)'
# string = 'ALL_batch_14_SRAE_1000189300.10.1'
# print(re.search(pattern1, string))


gtf_filtered = gtf_table[~gtf_table['transcript_id'].str.contains(pattern1, regex=True)]
gtf_filtered = gtf_filtered[~gtf_filtered['transcript_id'].str.contains(pattern2, regex=True)]
gtf_filtered = gtf_filtered[~gtf_filtered['transcript_id'].str.contains(pattern3, regex=True)]

get_gtf_file(gtf_filtered)

#%%


'''
creating a new gff3 file from the exisiting gtf and the cds/utr info
'''

gtf_file = make_annotations_table('consensus_noduplicates.gtf')
tmap_table = pd.read_csv('../Sratti_isoform/tmap_files/gffcmp.consensus.gtf.tmap', sep='\t', comment='#')

# creating a list of genes
genes = gtf_file[gtf_file['type'] == 'transcript']
# dup_genes = genes[genes.duplicated(['start', 'end'], keep=False)]
genes = genes.merge(tmap_table, left_on='transcript_id', right_on='qry_id', how='left')
genes.loc[genes['source'] == 'FLAIR', 'gene_id'] = genes['ref_gene_id']
mask = ~genes['gene_id'].str.contains('Gene:')
genes.loc[mask, 'gene_id'] = 'Gene:' + genes.loc[mask, 'gene_id']

genes = genes.rename(columns = {'gene_id':'ID'})
genes = genes.rename(columns = {'ref_gene_id':'Alias'})
genes = genes.drop(columns =  ['transcript_id','gene_name','xloc','class_code_x','cmp_ref','tss_id','exon_number','contained_in'])
genes = genes.iloc[:,:10]
genes.loc[(genes['Alias'] == '-'),'Alias'] = genes['ID']
genes['Parent'] = None
genes['type'] = 'gene'
genes = genes.drop_duplicates(subset='ID')

# creating list of transcripts
transcripts = gtf_file[gtf_file['type'] == 'transcript']
transcripts = transcripts.merge(tmap_table, left_on='transcript_id', right_on='qry_id', how='left')
transcripts.loc[transcripts['source'] == 'FLAIR', 'gene_id'] = transcripts['ref_gene_id']
mask = ~transcripts['gene_id'].str.contains('Gene')
transcripts.loc[mask,'gene_id'] = 'Gene:' + transcripts.loc[mask,'gene_id']

transcripts = transcripts.rename(columns = {'gene_id':'Parent'})
transcripts = transcripts.rename(columns = {'cmp_ref':'Alias'})
mask = ~transcripts['transcript_id'].str.contains('Transcript:')
transcripts.loc[mask, 'transcript_id'] = 'Transcript:' + transcripts.loc[mask, 'transcript_id']
transcripts.loc[transcripts['Alias'].isna(), 'Alias'] = transcripts['transcript_id']
transcripts = transcripts.drop(columns = ['gene_name','xloc','class_code_x','tss_id','exon_number','contained_in'])
transcripts = transcripts.iloc[:,:11]
transcripts = transcripts.rename(columns = {'transcript_id':'ID'})
transcripts = transcripts.drop_duplicates(subset='ID')

# changing the exon feature to be unique
for index,row in gtf_file.iterrows():
    if (row['type'] == 'transcript'):
        ID = row['transcript_id']
        if 'Transcript:' in ID:
            ID = ID.replace('Transcript:', '')
        c = 0
        #print(row['type'])
        #print(ID)
    elif(row['type'] == 'exon' and ID is not None):
        c += 1
        gtf_file.at[index, 'transcript_id'] = 'Exon:'+ ID + '.' + str(c)

# creating list of exons
exons = gtf_file[gtf_file['type'] == 'exon']
exons['Parent'] = exons['transcript_id'].str.replace(r'^Exon:|(\.\d+)?$', '', regex=True)
mask = ~exons['Parent'].str.contains('Transcript:')
exons.loc[mask,'Parent'] = 'Transcript:' + exons.loc[mask, 'Parent']
exons = exons.drop(columns = ['gene_id','gene_name','xloc','cmp_ref','class_code','tss_id','exon_number','contained_in'])
exons = exons.rename(columns = {'transcript_id':'ID'})
exons['Alias'] = None

# checking the columns are matchin in genes, transcripts and exons

result = genes['ID'].isin(transcripts['Parent'])
result = transcripts['ID'].isin(exons['Parent'])
print(result.value_counts())

# merging the exons, transcripts and genes

final_gff3 = pd.concat([genes, transcripts], ignore_index=True)
final_gff3 = pd.concat([final_gff3,exons], ignore_index=True) 

get_gff_file(final_gff3, 'test_gff3.gff3')

# merging cds and utrs

gff3_table = make_annotations_table('inputs/final_optimized_sorted.gff3') 

utr3_table = gff3_table[gff3_table['type'] == '3UTR_Length']
utr5_table = gff3_table[gff3_table['type'] == '5UTR_Length']
cds_table = gff3_table[gff3_table['type'] == 'CDS']

# removing duplicate entries
cds_table = cds_table[~cds_table.duplicated(['start', 'end'], keep=False)]
utr3_table = utr3_table[~utr3_table.duplicated(['start', 'end'], keep=False)]
utr5_table = utr5_table[~utr5_table.duplicated(['start', 'end'], keep=False)]

cds_table = cds_table.reset_index(drop=True)
cds_table['Parent'] = cds_table['Note'] 
cds_table['Note'] = cds_table['Note'].str.replace('Transcript:','')
cds_table['ID'] = cds_table['ID'] + ':' + cds_table['Note']
mask = ~cds_table['Parent'].str.contains('Transcript:')
cds_table.loc[mask,'Parent'] = 'Transcript:' + cds_table.loc[mask, 'Parent']
cds_table = cds_table.drop(cds_table.columns[9:15], axis = 1)

final_gff3 = pd.concat([final_gff3,cds_table], ignore_index=True)

get_gff_file(final_gff3, 'test_gff3.gff3')

utr3_table['Parent'] = utr3_table['Note']
utr3_table['Note'] = utr3_table['Note'].str.replace('Transcript:','')
utr3_table['ID'] = utr3_table['ID'] + ':' + utr3_table['Note']
mask = ~utr3_table['Parent'].str.contains('Transcript:')
utr3_table.loc[mask,'Parent'] = 'Transcript:' + utr3_table.loc[mask, 'Parent']
utr3_table = utr3_table.drop(utr3_table.columns[9:15], axis = 1) 

utr5_table['Parent'] = utr5_table['Note']
utr5_table['Note'] = utr5_table['Note'].str.replace('Transcript:','')
utr5_table['ID'] = utr5_table['ID'] + ':' + utr5_table['Note']
mask = ~utr5_table['Parent'].str.contains('Transcript:')
utr5_table.loc[mask,'Parent'] = 'Transcript:' + utr5_table.loc[mask, 'Parent']
utr5_table = utr5_table.drop(utr5_table.columns[9:15], axis = 1) 


final_gff3 = pd.concat([final_gff3,utr3_table], ignore_index=True)
final_gff3 = pd.concat([final_gff3,utr5_table], ignore_index=True)

final_gff3 = final_gff3.drop_duplicates(subset = 'ID')

get_gff_file(final_gff3)



result = utr5_table['Parent'].isin(transcripts['ID'])
print(result.value_counts())






#%%
from subprocess import call
from pathlib import Path
call('echo "AJ"', shell=True)

file_path = "/path/to/your/file/filename.txt"
file_name = Path(file_path).name

print(file_name)






























