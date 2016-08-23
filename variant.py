import os
import sys
import time
import urllib2
import json

class VARIANT:
	"""
	VARIANT class, which is used to parse one variant line in the VCF file
	A VARIANT instance will be create in ANNOTATE_VCF class for each variant line in the input VCF file
	"""
	def __init__(self, variant_line):
		"""
		Initialize the variable members for VARIANT instance
		"""
		self.variant_type	= "" # Variant type: deletion, insertion, or substitution
		self.alt_order		= "" # In ALT field, there could be multiple alts separated by ",", we want to record the number for selected variant_type
		self.tc			= "" # Total coverage at this locus
		self.tr			= "" # Total number of reads containing this variant; refer to self.alt_order to get selected number
		self.pct_tr		= "" # pct_tr = tr/tc

		# The following information can be fetched from ExAC server
		self.allele_freq	= "" # Allele frequency
		self.gene_symbol	= "" # Gene symbol
		self.sift		= "" # Prediction of the effects of coding non-synonymous variants on protein function
		self.polyphen		= "" # Possible impact of an amino acid substitution
		self.major_consequence	= "" # Variant type: synonymous, missense, intron, inframe_deletion, 3_prime_UTR etc
		self.existing_variation	= "" # variant ID, e.g. dbSNP ID
		
		# Here we assume there is only one sample. Please see below:
		# CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample
		self.chrom		= None 
		self.pos		= None
		self.myid		= None
		self.ref		= None
		self.alt		= None
		self.qual		= None
		self.myfilter		= None
		self.info		= None
		self.format		= None
		self.sample		= None
		self.chrom, self.pos, self.myid, self.ref, self.alt, self.qual, self.myfilter, self.info, self.format, self.sample = variant_line.split("\t")
		
		self.variant_type, self.alt_order	= self.get_variant_type()
		self.tc, self.tr			= self.get_tc_tr()
		try:
			self.pct_tr	= str(float(self.tr)/float(self.tc))
		except Exception as ex:
			print "Error in calculation of pct_tr, tc=%s, tr=%s: %s" % (self.tc, self.tr, str(ex))

		valid_chroms = ["X", "Y"]
		for i in xrange(1,23):
			valid_chroms.append(str(i))
		if self.chrom in valid_chroms:
			self.allele_freq, self.gene_symbol, self.sift, self.polyphen, self.major_consequence, self.existing_variation = self.get_anno_from_exac()
			
	def get_anno_from_exac(self):
		"""
		Fetch annotation information from ExAC server
		The parameter is VARIANT_STR with its format: CHROM-POS-REF-ALT
		"""
		URL_template = "http://exac.hms.harvard.edu/rest/variant/variant/VARIANT_STR"
		# example: http://exac.hms.harvard.edu/rest/variant/variant/22-46615880-T-C

		variant_str = "-".join([self.chrom, self.pos, self.ref, self.alt])
		url = URL_template.replace("VARIANT_STR", variant_str)
		the_page = self.get_page(url)
		variant_obj = json.loads(the_page)

		# anno_list is a return list for annotation fields found from ExAC
		# If no annotation found from ExAC, the field with be assigned with empty string ""
		anno_list = []
		if "allele_freq" in variant_obj:
			anno_list.append(str(variant_obj["allele_freq"]))
		else:
			anno_list.append("")
			
		# We only use the first transcript for simplicity, i.e. variant_obj["vep_annotations"][0]
		vep_anno_list = ["SYMBOL", "SIFT", "PolyPhen", "major_consequence", "Existing_variation"]
		for vep_field in vep_anno_list:
			vep_field_found = False
			if "vep_annotations" in variant_obj:
				if len(variant_obj["vep_annotations"]) > 0:
					if vep_field in variant_obj["vep_annotations"][0]:
						anno_list.append(variant_obj["vep_annotations"][0][vep_field])
						vep_field_found = True
			if not vep_field_found:
				anno_list.append("")

		return anno_list



	def get_page(self, url):
		"""
		Do a HTTP GET request and get the page
		We will try 20 times, because occasionally there could be some network problem.
		"""
		request = urllib2.Request(url)
		response = None
		TRY_TIMES = 20
		for try_count in xrange(1, TRY_TIMES):
			try: 
				#print try_times
				response = urllib2.urlopen(request, timeout=180)
				if response:
					break				
			except Exception as ex:
				print "Exception found for %s time(s)!\n%s" % (try_count, str(ex))
				time.sleep(10)
				continue
		if response:
			return response.read()
		else:
			print "After trying %d times, fail to get the page for %s\nExit now!\n" % (TRY_TIMES, url)
			sys.exit(1)



	
	def get_variant_type(self):
		"""
		Parse to obtain variant type by comparing REF and ALT.
		Variant type can be one of these: Deletion, Insertion, or Substitution
		"""
		# By comparing REF and alts (seperated by "," in ALT field), we can assign self.variant_type with one of following types:
		# Deletion, Insertion, or Substitution
		# If multiple variants type found in ALT column, then we only select one and the priority order is:
		# Deletion > Insertion > Substitution, which is exactly the same as the alphabetical order.
		variant_type_set = set()
		# mapping from variant_type to alt order, so that we know which alt is selected
		variant_type_to_alt_order = {} 
		ref	= self.ref # just use a short name in the following codes
		alts	= self.alt.split(",")
		for i in xrange(len(alts)):
			alt = alts[i]
			if (alt in ref) and (ref not in alt):
				variant_type_set.add("Deletion")
				variant_type_to_alt_order["Deletion"] = i
			elif (ref in alt) and (alt not in ref):
				variant_type_set.add("Insertion")
				variant_type_to_alt_order["Insertion"] = i
			else:
				variant_type_set.add("Substitution")
				variant_type_to_alt_order["Substitution"] = i
		variant_type = sorted(list(variant_type_set))[0]
		return [variant_type, variant_type_to_alt_order[variant_type]]


		
	def get_tc_tr(self):
		"""
		Parse to get variant type by comparing REF and ALT.
		Variant type can be one of these: Deletion, Insertion, or Substitution
		"""
		tc = ""
		tr = ""
		for info_field in self.info.split(";"):
			name, value = info_field.split("=")
			if name == "TC":
				tc = value
			elif name == "TR":
				tr = value.split(",")[self.alt_order]
		return [tc, tr]
				
				
			

			
		
		
		

