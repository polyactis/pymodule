#!/usr/bin/env python
"""
2011-8-28
	class to store functions related to next-gen sequencing
"""
import os, sys, csv, re
import utils

def getPEInputFiles(input_dir, isPE=True):
	"""
	2011-8-28
		copied from MpiBWA.py
	2011-8-5
		add argument isPE, which flags whether input_dir contains PE or single-end reads
		become a classmethod
	2011-2-7
		for paired-end files, sequence_628BWAAXX_1_1.fastq.gz and sequence_628BWAAXX_1_2.fastq.gz
			are regarded as one pair of two files.
	"""
	sys.stderr.write("Pair input files from %s ..."%input_dir)
	pairedEndPrefix2FileLs = {}
	files = os.listdir(input_dir)
	no_of_fastq_files = 0
	for fname in files:
		fname_prefix, fname_suffix = utils.getRealPrefixSuffixOfFilenameWithVariableSuffix(fname)
		if fname_suffix!='.fastq':		#skip non-fastq files
			continue
		no_of_fastq_files += 1
		if isPE==True:
			pairedEndPrefix = fname_prefix[:-2]
			pairedEndOrder = fname_prefix[-2:]
			
			if pairedEndPrefix not in pairedEndPrefix2FileLs:
				pairedEndPrefix2FileLs[pairedEndPrefix] = ['', '']
			
			if pairedEndOrder=='_1':	#the first file
				pairedEndPrefix2FileLs[pairedEndPrefix][0] = fname
			else:
				pairedEndPrefix2FileLs[pairedEndPrefix][1] = fname
		else:
			pairedEndPrefix2FileLs[fname_prefix] = [fname]	#single End
	no_of_files = len(files)
	no_of_pairedEndPrefix = len(pairedEndPrefix2FileLs)
	if no_of_pairedEndPrefix>0:
		avg_no_of_files_per_prefix = no_of_fastq_files/float(no_of_pairedEndPrefix)
	else:
		avg_no_of_files_per_prefix = 0.0
	sys.stderr.write("%.2f files per one pairedEnd prefix. %s fastq files. %s total files. Done.\n"%\
					(avg_no_of_files_per_prefix, no_of_fastq_files, no_of_files))
	return pairedEndPrefix2FileLs

def isFileNameVCF(inputFname, includeIndelVCF=False):
	"""
	2011-11-11
	"""
	isVCF = False
	if (inputFname[-3:]=='vcf' or inputFname[-6:]=='vcf.gz'):
		isVCF=True
		if not includeIndelVCF and inputFname.find('indel')!=-1:	#exclude indel vcf
			isVCF =False
	return isVCF

def isVCFFileEmpty(inputFname, checkContent=False):
	"""
	2011-11-11
		function to test if the input VCF has any locus.
		empty VCF file could still have headers.
	"""
	import os
	if not os.path.isfile(inputFname):
		return True
	fileSize = os.path.getsize(inputFname)
	if fileSize==0:
		return True
	if checkContent:
		if inputFname[-2:]=='gz':
			import gzip
			inf = gzip.open(inputFname)
		else:
			inf = open(inputFname)
		
		fileIsEmpty = True
		for line in inf:
			if line[0]!='#':
				fileIsEmpty=False
				break
		del inf
		return fileIsEmpty
	else:
		return False

if __name__ == '__main__':
	pass