import os
import shutil
import argparse
import csv
from pathlib import Path
from helper import Annotation

def main():
    parser = argparse.ArgumentParser(description = "Hunt Lab Annotation toolkit helps in manipulation of gff and gtf annotation files. It can also generate annotation statistics and extract of features of interest from provided annotation files.")
    # Positional arguments for CLI
    parser.add_argument("path", action = "store", type=str, help = "Path to input annotation (gff3/gtf) file")
    
    # Optional arguments
    parser.add_argument("-r", "--reference" , type=str, nargs="?", help="Path to reference genome")
    parser.add_argument("-m","--make_csv", action="store_true" , help="create the annotation file into a table in csv format")
    parser.add_argument("-e" , "--extract", type=str, help="extract the set of transcript_ids provided in a txt file and create a new gtf file for provided transcripts")
    parser.add_argument("-c" , "--count", action="store_true", help="counts all the features in the annotation file")
    parser.add_argument("-i", "--compare", type=str , help="compare the input annotation file with provided file and identify novel transcripts and genes")
    # parser.add_argument("-s", "--sort", action="store_true", help="sorts the provided input annotation file")
    parser.add_argument("-f", "--get_fasta", choices=["transcript", "utr", "cds", "exon"] , help="Generate fasta file for the provided feature, calculate the sequence length distributions and create a length distribution plot")
    parser.add_argument("-p", "--output_prefix", type=str, default="output" )
    
    args = parser.parse_args()
    
    print('********RUNNING HL-ANNOTATION TOOLKIT********')
    
    output_directory = os.path.join(os.getcwd(), 'Annotation_toolkit_outputs')
    if os.path.exists(output_directory):
        shutil.rmtree(output_directory)
    os.makedirs(output_directory)
    
    annotation = Annotation(file_path=args.path, 
                            output_path=output_directory)
    
    file_path = Path(annotation.file_path)
    file_name = file_path.name
    extension = file_path.suffix
    
    if extension not in ['.gtf','.gff3', '.gff']:
        parser.exit(status=0, message="Input error: .gtf or .gff3 input required, but received other")
    
    if(args.make_csv):
        # print(annotation.output_path)
        print('Generating csv file..')
        annotations_table = annotation.make_annotations_table(annotation.file_path)
        output_file_path = os.path.join(annotation.output_path, f"{args.output_prefix}_table.csv")
        annotations_table.to_csv(output_file_path, index=False)
        print('Done')
    
    if(args.extract):
       print('Extracting the provided transcript features...') 
       transcript_file = args.extract
       table = annotation.extract_transcripts(annotation.file_path, transcript_file)
       if table.empty:
           parser.exit(status=0, message="Input error: .gtf or .gff3 input required, but received other")
       else:
           output_file_path = os.path.join(annotation.output_path, f"{args.output_prefix}_extracted{extension}")
           table.to_csv(output_file_path , sep='\t', index=False, header=False, quoting=csv.QUOTE_NONE)
           print('Done')
    
    if(args.count):
        print('\n----------------------------------------------------------------')
        annotation.get_feature_counts(annotation.file_path)
        print('-----------------------------------------------------------------\n')
    
    if(args.compare):
        print('Comparing input annotation files...')
        comparison_file = Path(args.compare)
        if comparison_file.suffix not in ['.gtf','.gff3','.gff']:
            parser.exit(status=0, message="Comparison file error: .gtf or .gff3 input required, but received other")
        tmap_file = f'gffcmp.{file_name}.tmap'
        annotation.compare_annotation_files(args.compare, annotation.file_path, tmap_file)
        
        gffcmp_directory = os.path.join(annotation.output_path, 'gffcmp_files')
        os.makedirs(gffcmp_directory)
        
        for filename in os.listdir(file_path.parent):
            if filename.startswith('gffcmp'):
                source_path = os.path.join(file_path.parent, filename)
                destination_path = os.path.join(gffcmp_directory, filename)
                shutil.move(source_path, destination_path)
        print('Done')
    
    if(args.get_fasta):
        if (args.reference):
            reference_path = Path(args.reference)
            if reference_path.suffix not in ['.fa','.fasta','.fai']:
                parser.exit(status=0, message="Fasta file input required for refrence genome, but received other")
            print('Creating fasta file...')
            feature = args.get_fasta
            annotation.get_fasta_file(args.reference, feature, annotation.file_path)
            for filename in os.listdir(os.getcwd()):
                if filename.startswith(feature):
                    source_path = os.path.join(os.getcwd(), filename)
                    destination_path = os.path.join(output_directory,filename)
                    shutil.move(source_path, output_directory)
            print('Done')
        else:
            parser.exit(status=0, message="Error: Need additional input. Reference genome input required through -r tag")
               
                    
if __name__ == '__main__':
    main()