import os
import sys
import argparse
from variant import VARIANT
# A simple VCF annotator
# Author: Zhengdeng Lei (zlei2@uic.edu), last update 8/23/2016

class ANNOTATE_VCF:
	"""
	VCF annotation class, which is used to parse a VCF file and create an annotated tab-delimited file with VCF header removed.
	The user can use Excel (or text editor) to open the output file.
	"""
	def __init__(self, in_vcf, out_vcf):
		"""
		ANNOTATE_VCF constructor, initialize the input/output file handles
		"""

		self.in_vcf	= in_vcf
		self.out_vcf	= out_vcf
		self.in_vcf_fh	= None
		self.out_vcf_fh	= None

		# Assign output VCF file name if the name is not provided in command line
		if self.out_vcf is None:
			self.out_vcf = self.in_vcf + ".annotated.txt"

		# Make sure the input VCF file is readable and the output VCF file is writable
		try:
			self.in_vcf_fh	= open(self.in_vcf)
		except IOError as ex:
			self.__io_error(ex, self.in_vcf)
		try:
			self.out_vcf_fh	= open(self.out_vcf, "wt")
		except IOError as ex:
			self.__io_error(str(ex), self.out_vcf)
		self.parse_vcf()		


	def __io_error(self, ex_message, filename):
		"""
		Error in opening input/output file
		"""
		print 'Error: open %s failed.\n%s' % (filename, ex_message)
		if self.in_vcf_fh:
			self.in_vcf_fh.close()
		if self.out_vcf_fh:
			self.out_vcf_fh.close()
		sys.exit(1)

	
	def __del__(self):
		"""
		Close input/output file handles
		"""
		if self.in_vcf_fh:
			self.in_vcf_fh.close()
		if self.out_vcf_fh:
			self.out_vcf_fh.close()


	def parse_vcf(self):
		"""
		Parse the input VCF file, for each variant, it will be processed in VARIANT class
		Output will be an annotated tab-delimited file with VCF header removed.
		"""		
		header1 = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "sample"]
		header2 = ["Variant type", "Total coverage (TC)", "Total reads supporting variant (TR)", "TR/TC", "Allele freq", "SYMBOL", "SIFT", "PolyPhen", "major_consequence", "Existing_variation"]
		header = header1 + header2
		self.out_vcf_fh.write("\t".join(header) + os.linesep)

		# parse input vcf file, ignore the comment lines
		variant_count = 1
		for line in self.in_vcf_fh:
			if not line.startswith("#"):
				print "processing variant %d" % variant_count
				line = line.rstrip(os.linesep)
				v = VARIANT(line)
				out_str = "\t".join([line, v.variant_type, v.tr, v.tc, v.pct_tr, v.allele_freq, v.gene_symbol, v.sift, v.polyphen, v.major_consequence, v.existing_variation, os.linesep])
				self.out_vcf_fh.write(out_str)
				variant_count = variant_count + 1
				

		
def get_arguments():
	main_description = """\
	A simple VCF annotator
	Author: Zhengdeng Lei (zlei2@uic.edu)
	"""
	help_help = """\
	show this help message and exit\
	"""
	version_help = """\
	show the version of this program\
	"""
	input_vcf_help = """\
	input vcf file
	"""

	output_vcf_help = """\
	output vcf file with VCF header removed
	default is <input_vcf>.annotated.txt
	"""

	arg_parser = argparse.ArgumentParser(description=main_description, formatter_class=argparse.RawTextHelpFormatter, add_help=False)
	# register boolean type of input argument
	# arg_parser.register('type','bool', lambda s: str(s).lower() in ['true', '1', 't', 'y', 'yes']) # add type keyword to registries

	###############################
	#    1. required arguments    #
	###############################
	required_group = arg_parser.add_argument_group("required arguments")
	required_group.add_argument("-i", dest="input_vcf",  action="store", required=True, default=None, help=input_vcf_help)

	###############################
	#    2. optional arguments    #
	###############################
	optional_group = arg_parser.add_argument_group("optional arguments")
	optional_group.add_argument("-o", dest="output_vcf", default=None, help=output_vcf_help)
	optional_group.add_argument("-h", "--help", action="help", help=help_help)
	optional_group.add_argument("-v", "--version", action="version", version="%(prog)s: version 0.1", help=version_help)

	args = arg_parser.parse_args()
	args.input_vcf = os.path.abspath(args.input_vcf)

	if not os.path.exists(args.input_vcf):
		arg_parser.print_help()
		print "\n\nError: input vcf does not exist!\n"
		sys.exit(1)

	return args



def main():
	args = get_arguments()
	anno_vcf = ANNOTATE_VCF(args.input_vcf, args.output_vcf)
	
		
if __name__ == "__main__":
	main()
