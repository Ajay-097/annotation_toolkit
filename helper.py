import pandas as pd
import csv
from pathlib import Path
from Bio import SeqIO
import seaborn as sns
import matplotlib.pyplot as plt


class Annotation :
    
    def __init__(self, file_path:Path, output_path:Path = None):
        self.file_path = file_path
        self.output_path = output_path
    
    def validate_file (self, file_path) :
        return True
    

    def make_anotations_table (self, file_path):
        """
        The following finction makes the annotation file (gtf/gff) into a large dataframe
        including spearate columns for each of the attrributes
        """
        
        if not self.validate_file(file_path):
            print("File error")
        
        attr_dict_list = []
        columns = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
        annotation_df = pd.read_csv(self.file_path, sep='\t', comment='#', header=None, names=columns)
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
        final_dataset.to_csv(self.output_path+".csv", index=False)
    
    def make_transcript_list(transcript_file):
        """
        Function returns a list of transcript_ids from the transcript_ids text file
        """
        transcript_list = []
        with open(transcript_file, 'r') as file:
            for line in file:
                transcript_list.append(line.strip())
        return transcript_list
    
    def get_tmap_stats(tmap_file):
        tmap_table = pd.read_csv(tmap_file, sep='\t', comment='#')
        print('Query transcript: ', tmap_table['qry_id'].nunique())
        print('Reference transcript: ',tmap_table['ref_id'].nunique())
        print('Query gene: ', tmap_table['qry_gene_id'].nunique())
        print('Reference gene: ', tmap_table['ref_gene_id'].nunique())
        
    def get_lenghtDistribution_from_fasta (fasta_file):
        sequence_lengths = []
        with open(fasta_file, "r") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                length = len(record.seq)
                sequence_lengths.append(length)
        sns.boxplot(sequence_lengths)#, showfliers = False)
        plt.ylabel('Sequence lenghts (nt)')
        plt.title('Transcript lenght distributions')
        plt.show()
        


