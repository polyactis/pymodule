#!/usr/bin/env python3
"""
Examples:
	%s 
	
	%s -i ~namtran/panasas/Experiment/RNA-seq/Freimer/Developmental/ASE/Variant/MultipleSampleCalling/genome.algn.split.part17/samtools.var.filt.vcf.gz 
		-o /tmp/ -m 2

Description:
	2012.5.10
		
"""

import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import csv
from palos import ProcessOptions, getListOutOfStr, PassingData, utils
from palos import VCFFile
from palos.pegasus.mapper.AbstractVCFMapper import AbstractVCFMapper
from vervet.src import VervetDB

class ConvertAlignmentReadGroup2UCLAIDInVCF(AbstractVCFMapper):
	__doc__ = __doc__
	option_default_dict = AbstractVCFMapper.option_default_dict.copy()
	option_default_dict.pop(("chromosome", 0, ))
	option_default_dict.pop(("chrLength", 1, int))
	option_default_dict.pop(('outputFnamePrefix', 0, ))
	option_default_dict.update({
						('includeIndels', 0, ): [0, 'n', 0, 'By default, indels are excluded. INFO tag contains INDEL.', ],\
						})
	def __init__(self,  **keywords):
		"""
		"""
		AbstractVCFMapper.__init__(self, **keywords)
		import re
		self.chr_pattern = re.compile(r'(\w+\d+).*')
		self.contig_number_pattern = re.compile(r'Contig(?P<contigNumber>\d+)')

	def convertAlignmentReadGroup2UCLAIDInVCF(self, inputFname, outputFname, minDepth=1, includeIndels=False,\
											maxContigNumber=None):
		"""
		2012.5.10
		"""
		sys.stderr.write("Converting %s from VCF to EigenStrat ...\n"%(inputFname))
		
		vcfFile = VCFFile(inputFname=inputFname, minDepth=minDepth)
		#replace Variant/PooledTissues/2002053/genome.algn.split.part17/5tissues.pooled.rmdup.bam with just monkey ID
		
		newSampleIDHeader = []
		for sampleID in vcfFile.sampleIDHeader:
			readGroupData = VervetDB.VervetDB.parseAlignmentReadGroupWithoutDB(sampleID)
			UCLAID = readGroupData.individual_code
			newSampleIDHeader.append(UCLAID)
		#new header for every output contig
		newHeader = vcfFile.header[:vcfFile.sampleStartingColumn] + newSampleIDHeader
		
		counter = 0
		real_counter = 0
		outVCFFile = VCFFile(outputFname=outputFname)
		outVCFFile.metaInfoLs = vcfFile.metaInfoLs
		outVCFFile.header = newHeader
		outVCFFile.writeMetaAndHeader()
		for vcfRecord in vcfFile.parseIter():
			counter += 1
			if not includeIndels and (len(vcfRecord.refBase)!=1 or len(vcfRecord.altBase)!=1):
				#it's an indel if refBase or altBase is not just one base
				continue
			
			chr = vcfRecord.chr
			if maxContigNumber:
				contigNumber = int(self.contig_number_pattern.search(chr).group('contigNumber'))
				if contigNumber>maxContigNumber:
					continue
			real_counter += 1
			# set genotype whose depth is below minDepth to ./. (=missing)
			for i in xrange(1, len(vcfRecord.data_row)):	#[0] is the ref base
				callData = vcfRecord.data_row[i]
				if callData is None or callData.get('DP',0)<minDepth:
					sampleColumnIndex = i+vcfFile.sampleStartingColumn-1
					vcfRecord.row[sampleColumnIndex] = './.'
			outVCFFile.writeVCFRecord(vcfRecord)
		
		vcfFile.close()
		#close all output files
		outVCFFile.close()
		
		sys.stderr.write("%s (out of %s) loci.\n"%(real_counter, counter))
	
	def run(self):
		"""
		"""
		
		if self.debug:
			import pdb
			pdb.set_trace()
		
		self.convertAlignmentReadGroup2UCLAIDInVCF(self.inputFname, self.outputFname, minDepth=self.minDepth,\
												includeIndels=self.includeIndels)

if __name__ == '__main__':
	main_class = ConvertAlignmentReadGroup2UCLAIDInVCF
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()