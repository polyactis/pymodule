#!/usr/bin/env python3
"""
2012.8.13
	this program outputs chr, pos, gap (distance to the next SNP).

Examples:
	%s -i gatk/Contig799.vcf.gz -l 1000000 -c Contig799 -o /tmp/Contig799_siteGap.tsv

"""

import sys, os, math
__doc__ = __doc__%(sys.argv[0])
import csv
from palos import ProcessOptions, getListOutOfStr, PassingData, utils
from palos.ngs.io.VCFFile import VCFFile
from palos.mapper.AbstractVCFMapper import AbstractVCFMapper

class OutputVCFSiteGap(AbstractVCFMapper):
	__doc__ = __doc__
	option_default_dict = AbstractVCFMapper.option_default_dict.copy()

	def __init__(self,  **keywords):
		"""
		"""
		AbstractVCFMapper.__init__(self, **keywords)
	
	
	def calculateSiteGap(self, inputFname, outputFname, chromosome=None, chrLength=None, minDepth=1):
		"""
		2011-11-2
			given a VCF file, count the number of homo-ref, homo-alt, het calls
			
		"""
		sys.stderr.write("Calculate the distances between sites of %s .\n"%(inputFname))
		writer = csv.writer(open(outputFname, 'w'), delimiter='\t')
		writer.writerow(['chromosome', 'position', 'length', "distanceToNextSite"])
		vcfFile = VCFFile(inputFname=inputFname, minDepth=minDepth)
		
		no_of_total = 0.
		minStart = None
		previousPosition = None
		for vcfRecord in vcfFile.parseIter():
			chr = vcfRecord.chr
			pos = vcfRecord.pos
			pos = int(pos)
			if previousPosition is not None:
				distanceToNextSite = pos-previousPosition
				data_row = [chr, previousPosition, chrLength, distanceToNextSite]
				writer.writerow(data_row)
			previousPosition = pos
		del writer
		sys.stderr.write("Done.\n")
	
	def run(self):
		"""
		"""
		
		if self.debug:
			import pdb
			pdb.set_trace()
		
		
		self.calculateSiteGap(self.inputFname, self.outputFname, chromosome=self.chromosome, \
												chrLength=self.chrLength, minDepth=self.minDepth)

if __name__ == '__main__':
	main_class = OutputVCFSiteGap
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()