import pandas as pd
import subprocess
import sys
import csv
from Bio import SeqIO
import seaborn as sns
import matplotlib.pyplot as plt


class Annotation :
    
    def __init__(self, file_path, output_path):
        self.file_path = file_path
        self.output_path = output_path
        

    def make_annotations_table (self, file_path):
        """
        The following finction makes the annotation file (gtf/gff) into a large dataframe
        including spearate columns for each of the attrributes
        """
        # if not self.validate_file(file_path):
        #     print("File error")
        #     return pd.DataFrame(),0
        # print (file_path)
        attr_dict_list = []
        columns = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
        annotation_df = pd.read_csv(self.file_path, sep='\t', comment='#', header=None, names=columns)
        for entry in annotation_df['attributes'].str.split(';'):
            attr_dict = {}
            for attr in entry:
                if(attr != '' and file_path.endswith('.gtf')):
                    key, value = attr.strip().split(' ')
                    attr_dict[key] = value.strip('"')
                elif file_path.endswith('.gff3') or file_path.endswith('.gff'):
                    key, value = attr.strip().split('=')
                    attr_dict[key] = value.strip()
            attr_dict_list.append(attr_dict)
        attributes_df = pd.DataFrame(attr_dict_list)
        final_dataset = pd.concat([annotation_df,attributes_df],axis=1)
        final_dataset.drop('attributes', axis=1, inplace=True) 
        # print("Number of unique transcripts = ",final_dataset['transcript_id'].nunique()) 
        return final_dataset
    
    def extract_transcripts(self, file_path, transcript_file):
        transcript_list = []
        with open(transcript_file, 'r') as file:
            for line in file:
                transcript_list.append(line.strip())
        if '.gtf' in file_path:
            gtf_table = self.make_annotations_table(file_path)
            gtf_file = self.get_gtf_file(gtf_table, transcript_list)
            return gtf_file
        elif '.gff3' in file_path:
            gff_table = self.make_annotations_table(file_path)
            gff_file = self.get_gff_file(gff_table, transcript_list)
            return gff_file
        else:
            return pd.DataFrame()
    
    def get_gff_file (self, gff_table, transcript_list=None):
        """
        When a set of transcript Ids are provided the function returns a gff file with
        just the annotations for the provided list of transcripts.
        """
        if(transcript_list != None):
            search_str = '|'.join(transcript_list)
            search_str = search_str.replace('.','\.')
            gff_table = gff_table[gff_table['ID'].str.contains(search_str, case=False)]
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
        
    def get_gtf_file (self, gtf_table, transcript_list=None):
        """
        When a set of transcript Ids are provided the function returns a gtf file with
        just the annotations for the provided list of transcripts.
        """
        if(transcript_list!=None):
            search_str = '|'.join(transcript_list)
            search_str = search_str.replace('.','\.')
            gtf_table = gtf_table[gtf_table['transcript_id'].str.contains(search_str, case=False)]
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
        return final_table
    
    def get_feature_counts (self, file_path):
        """
        Function returns the feature counts of all the features present in a annotations file.
        """
        annotation_table = self.make_annotations_table(file_path)
        pivot_df = annotation_table.groupby('type').size().reset_index()
        pivot_df.rename(columns={0: 'Count'}, inplace=True)
        pivot_df.rename(columns={'type': 'Feature'}, inplace=True)
        print(pivot_df)
    
    def compare_annotation_files(self, reference_file, input_file, tmap_file_path):
        '''
        Function runs gff compare on the two annotation files and generates a tmap file
        The tmap file is used to calculate the matching statistics
        '''
        try:
            print('Running gff comapre...\n')
            result = subprocess.run(f'gffcompare -r {reference_file} -R {input_file}',shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            if result.returncode !=0 :
                print(f' Error: {result.stderr}')
                sys.exit(1)
            # subprocess.call(f'mkdir {self.output_path}\gff_comapre', shell=True)
        except Exception as e:
            print(f'An error occured: {e}')
            sys.exit(1)
        
        tmap_table = pd.read_csv(tmap_file_path, sep='\t', comment='#')
        print('\n----------------------------------------------------------------')
        print('No. of transcripts in input file: ', tmap_table['qry_id'].nunique() )
        print('No. of matching transcripts (including partial matches): ',tmap_table['ref_id'].nunique())
        print('No. of exact matching transcripts (classcode =): ', tmap_table[tmap_table['class_code'] == '=']['qry_id'].nunique())
        print('No. of novel transcripts: ', tmap_table[tmap_table['class_code'] == 'u']['qry_id'].nunique())
        
        print('No. of genes in input file: ', tmap_table['qry_gene_id'].nunique())
        print('No. of matching genes: ', tmap_table['ref_gene_id'].nunique())
        print('No. of novel genes: ', tmap_table[tmap_table['class_code']=='u']['qry_gene_id'].nunique())
        print('----------------------------------------------------------------\n')
        # subprocess.call(f'mv {os.getcwd()}\gffcmp* {self.output_path}\gff_compare', shell=True)
    
    def get_fasta_file(self, reference, feature, input_file):
        
        annotation_df = self.make_annotations_table(input_file)
        if feature == 'transcript' :
            try:
                result = subprocess.run(f'gffread -w transcript.fa -g {reference} {input_file}', shell=True, 
                                        stdout=subprocess.PIPE, 
                                        stderr=subprocess.PIPE, text=True
                                        )
                if(result.returncode != 0):
                    print(result.stderr)
                    sys.exit(1)
            except Exception as e:
                print(f'An error occured: {e}')
                sys.exit(1)
            self.get_lenghtDistribution_from_fasta('transcript.fa')
        else:
            feature_df = annotation_df[annotation_df['type'].str.contains(feature, case=False)]
            if(feature_df.empty):
                print(' The given feature is not present in the annotation file')
                sys.exit(1)
            feature_gff = self.get_gff_file(feature_df)
            feature_gff.to_csv(f'{feature}.gff3' , sep='\t', index=False, header=False, quoting=csv.QUOTE_NONE)
            try:
                result = subprocess.run(f'gffread -w {feature}.fa -g {reference} {feature}.gff3', shell=True, 
                                        stdout=subprocess.PIPE, 
                                        stderr=subprocess.PIPE, text=True
                                        )
                if(result.returncode != 0):
                    print(result.stderr)
                    sys.exit(1)
            except Exception as e:
                print(f'An error occured: {e}')
                sys.exit(1)
            self.get_lenghtDistribution_from_fasta(f'{feature}.fa')
    
    def get_lenghtDistribution_from_fasta (self, fasta_file):
        sequence_lengths = []
        with open(fasta_file, "r") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                length = len(record.seq)
                sequence_lengths.append(length)
        sns.set_palette("viridis")
        sns.boxplot(sequence_lengths, color='skyblue', showfliers = False, showmeans=True)
        plt.ylabel('Sequence lenghts (nt)')
        plt.title('Lenght distribution box plot')
        # output_file_path = os.path.join(self.output_path, f'{fasta_file}.LDplot.png')
        plt.savefig(f'{fasta_file}.LDplot.png')
        #  plt.show()
        


