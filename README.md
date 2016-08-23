# annotate_vcf
A simple VCF annotator.
The program will parse an input vcf as well as fetch annotation information from ExAC server.

The input is a VCF file, and the output is an annotated tab-delimited file with VCF header removed.

You can use Excel (or text editor) to open the output file.

To run the program, you need to have the following environment setup:

(1) python 2.6 or 2.7;

(2) python packages: argparse, urllib2, json;

(3) Put annotate_vcf.py and variant.py into a same directory.

Type "python annotate_vcf.py -h" for help.

