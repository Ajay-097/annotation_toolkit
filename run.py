import os
import argparse
from helper import Annotation

def main():
    parser = argparse.ArgumentParser()
    # Positional arguments for CLI
    parser.add_argument("path", action = "store", type=argparse.FileType('r'), help = "Path to input annotation (gff/gtf) file")
    
    # Optional arguments
    parser.add_argument("-r", "--reference" , type=str, nargs="?", help="Path to reference genome")
    parser.add_argument("-m","--make_csv", action="store_true" , help="create the annotation file into a table in csv format")
    parser.add_argument("-e" , "--transcript_list", type=str, help="extract the set of transcript_ids provided in a txt file and create a new gtf file for provided transcripts")
    parser.add_argument("-c" , "--count", action="store_true", help="counts all the features in the annotation file")
    parser.add_argument("-i", "--identify_novel", type=str , help="compare the input annotation file with provided file and identify novel transcripts and genes")
    parser.add_argument("-s", "--sort", action="store_true", help="sorts the provided input annotation file")
    parser.add_argument("-p", "--output_prefix", type=str, default="output" )
    parser.add_argument("-o" , "--output_directory", type=str, help="output directory path", default = os.getcwd())
    
    args = parser.parse_args()

    annotation = Annotation(file_path=args.path, output_path=args.output_directory+f"\{args.output_prefix}")
    
    if(args.make_csv):
        print(annotation.output_path)
        annotation.make_anotations_table(annotation.file_path)
        


if __name__ == '__main__':
    main()