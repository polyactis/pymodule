#!/usr/bin/env python
"""
Examples:
	%s -i input.vcf -o -o selectedStKitts.vcf --country_id_ls 144 --tax_id_ls 60711
	
	%s -i input.vcf -o -o selectedNevis.vcf --country_id_ls 148 --tax_id_ls 60711

	%s  -a ~/NetworkData/vervet/db/individual_sequence/524_superContigsMinSize2000.fasta
		-i foldermap/Contig1_54079_VCF_52496_VCF_49977_VCF_Contig1_splitVCF_u1_popNevis.vcf
		-o /tmp/Contig1_54079_VCF_52496_VCF_49977_VCF_Contig1_splitVCF_u1_popNevis_flankSeq.tsv
	

Description:
	2012.10.5 program that extracts samples from a VCF and form a new VCF.
		need to re-calculate the AC/AF values of each variant.
"""

import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0], sys.argv[0])


sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import csv, re
from pymodule import ProcessOptions, getListOutOfStr, PassingData, getColName2IndexFromHeader, figureOutDelimiter, SNPData
from pymodule.utils import getColName2IndexFromHeader, getListOutOfStr, figureOutDelimiter
from pymodule import yh_matplotlib, GenomeDB, utils
from pymodule.MatrixFile import MatrixFile
from pymodule import SNP

#used in getattr(individual_site_id_set, '__len__', returnZeroFunc)()
from pymodule.utils import returnZeroFunc
from pymodule.VCFFile import VCFFile
from pymodule.Genome.FastaFile import FastaFile
from AbstractVCFMapper import AbstractVCFMapper


class ExtractFlankingSequenceForVCFLoci(AbstractVCFMapper):
	__doc__ = __doc__
	option_default_dict = AbstractVCFMapper.option_default_dict
	option_default_dict.update({
						('refFastaFname', 1, ): [None, 'a', 1, 'path to the reference sequence fasta file', ],\
						('flankingLength', 1, int): [24, '', 1, 'number of flanking bases on either side of the locus.\n\
	length of flanking = 2*flankingLength+locusLength', ],\
						('outputFormatType', 1, int): [1, '', 1, 'output format type. 1: fasta, 2: fastq (standard quality=H (72-33) format)'],\
					})
	def __init__(self, inputFnameLs=None, **keywords):
		"""
		"""
		AbstractVCFMapper.__init__(self, inputFnameLs=inputFnameLs, **keywords)
	
	#2012.10.8 a regular expression pattern to parse the title of fasta sequences that are outputted by this program
	sequenceTitlePattern = 	re.compile(r'^([a-zA-Z0-9]+)_(\d+)_(\d+)_([a-zA-Z\-]*)_([a-zA-Z\-,]*)_positionInFlank(\d+)$')

	def extractFlankingSequence(self, inputFname=None, refFastaFname=None, outputFname=None, flankingLength=24,\
							outputFormatType=1):
		"""
		2012.10.10
			added argument outputFormatType. 1: fasta, 2: fastq
		2012.10.8
		"""
		sys.stderr.write("Extracting flanking sequences of loci from %s, based on ref-sequence of %s ...\n"%\
						(inputFname, refFastaFname))
		vcfFile = VCFFile(inputFname=inputFname)
		outf = open(outputFname, 'w')
		refFastaFile = FastaFile(inputFname=refFastaFname)
		
		counter = 0
		for vcfRecord in vcfFile:
			refBase = vcfRecord.refBase
			stopPos = vcfRecord.pos + len(refBase) -1
			
			SNP_ID = '%s_%s_%s_%s_%s'%(vcfRecord.chr, vcfRecord.pos, stopPos, vcfRecord.refBase, vcfRecord.altBase)
			fastaTitle = '%s_positionInFlank%s'%(SNP_ID, flankingLength+1)	#positionInFlank is 1-based.
			flankSeqStart = max(1, vcfRecord.pos-flankingLength)
			flankSeqStop = stopPos + flankingLength
			flankingSequence = refFastaFile.getSequence(vcfRecord.chr, start=flankSeqStart, stop=flankSeqStop)
			if flankingSequence:
				if outputFormatType==1:
					outf.write(">%s\n"%(fastaTitle))
					outf.write('%s\n'%(flankingSequence))
				else:
					outf.write("@%s\n"%(fastaTitle))
					outf.write('%s\n'%(flankingSequence))
					outf.write("+\n")
					outf.write("%s\n"%('H'*len(flankingSequence)))
					
			
			counter += 1
		
		del outf
		vcfFile.close()
		refFastaFile.close()
		sys.stderr.write("%s loci.\n"%(counter))
	
	def run(self):
		"""
		2012.10.5
		"""
		
		if self.debug:
			import pdb
			pdb.set_trace()
		
		self.extractFlankingSequence(inputFname=self.inputFname, refFastaFname=self.refFastaFname, \
									outputFname=self.outputFname, flankingLength=self.flankingLength,\
									outputFormatType=self.outputFormatType)
		

if __name__ == '__main__':
	main_class = ExtractFlankingSequenceForVCFLoci
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(po.arguments, **po.long_option2value)
	instance.run()