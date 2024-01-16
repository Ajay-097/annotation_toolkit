# HL_annotation_toolkit
The Hunt Lab Annotation toolkit helps in manipulation of gff and gtf annotation files. It can also generate annotation statistics and extract of features of interest from provided annotation files. 

# Installation
The scripts can be accessed by simply cloning the git repository to your local machine.

```bash
git clone <copied url>
```

# Pre-requisites
Make sure gffcompare and Python is already isntalled in your machine

# Running the toolkit
The toolkit can be run by executing the following command
```bash
python run.py <optional_arguments> <path/to/annotation_file>

# Optional arguments are as follows
-r --reference <path/to/reference_genome>
-m --make_csv #create the annotation file into a table in csv format
-e --extracts <path/to/text_file_with_transcript_ids> #extract the set of transcript_ids provided in a txt file and create a new gtf file for provided transcripts
-c --count #displays total counts of all different features in the annotation file
-i --compare <path/to/comparison_annotation_file> #compare the input annotation file with provided comparison file and displays counts of novel transcripts and genes
-f --get_fasta {transcript,utr,cds,exon} #Generate fasta file for the provided feature calculate the sequence length distributions and createa length distribution plot. (reference genome has to be provided) 
-p --output_prefix <prefix_string> #default:output
```
Note that the toolkit only accepts .gtf, .gff3 or .gff inputs that are not gzipped.<br> 
All the outputs generated by the toolkit are stored in the 'Annotation_toolkit_outputs' directory. The directory gets over written everytime the toolkit is run. <br>

# Usage
The HL annotation toolkit currently has the following functionalities

1. Make the annotation file into a csv file including separate columns for each of the attributes
```bash
python run.py -m <path/to/annotation_file>
```  

2. Creating a new annotation file (gtf) with only certain transcripts of interest. The list of transcript IDs should be provided as a txt file.<br>
Extracted transcripts are in gtf format in the output dirtectory with name <output_prefix_table.csv>
```bash
python run.py -e  <path/to/text_file_with_transcript_ids> <path/to/annotation_file>
```

3. Total counts of each feature present in the annotation file are displayed with the --count tag.
```bash
python run.py -c <path/to/text_file_with_transcript_ids> <path/to/annotation_file>
```

4. The total transcripts and genes present in an annotation file can be compared with a comparison file to identify the novel and matching features.<br>
Novel denotes - Transcripts that are present in the input file but not present in the comparison file <br>
Matching denotes - Both exact matches and partial matches are considered matching transcripts <br>
```bash
python run.py -i <path/to/comparison_annotation_file> <path/to/annotation_file>
```

5. The --get_fasta option can be used to fetch the sequences of required features in fasta format.Currently the tollkit can fetch utrs, exons, trasnscripts and cds sequences. 
Reference genome input is mandatory for this option.
```bash
python run.py -r <path/to/reference_genome> -f <transcript/utr/cds/exon> <path/to/annotation_file>
``` 
