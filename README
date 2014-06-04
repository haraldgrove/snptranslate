
The project 'snptranslate' is for converting between the various genotype formats.
The script will read one and one file from the input file and translate it
as necessary before writing it to the output file.

Contact: harald.grove@nmbu.no

*HELP TEXT*

Running the script with the -h option gives the following help text
(the list of formats will automatically update when new formats are added):

usage: convSNP.py [-h] [-i GENOSFILE] [-n INFORM] [-p PEDIGREEFILE]
                  [-m MARKERFILE] [-o OUTPUTFILE] [-u OUTFORM] [-v VERBOSE]

Converts between the following SNP-formats:
 Plink DMU Geno LD Alphaimpute Simplegeno

optional arguments:
  -h, --help            show this help message and exit
  -i GENOSFILE, --genotypes GENOSFILE
                        Input file
  -n INFORM, --informat INFORM
                        Input file format
  -p PEDIGREEFILE, --pedigree PEDIGREEFILE
                        Pedigree file
  -m MARKERFILE, --markers MARKERFILE
                        Marker file
  -o OUTPUTFILE, --output OUTPUTFILE
                        Output file
  -u OUTFORM, --outformat OUTFORM
                        Output file format
  -v VERBOSE, --verbose VERBOSE
                        Prints runtime info

*EXAMPLE*

Translating a file called 'input.txt' from DMU to Plink format and store it in a new file 'output.txt',
with marker information located in a file 'marker.txt':

convSNP.py -i input.txt -n DMU -o output.txt -u Plink -m marker.txt

