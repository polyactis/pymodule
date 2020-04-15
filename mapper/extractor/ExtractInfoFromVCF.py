#!/usr/bin/env python3
"""
2012.9.17 program that extracts certain key in the INFO column of VCF (could be multiple input files) and outputs them in tsv format.

Examples:
	%s -i ~/NetworkData/vervet/db/genotype_file/method_27/*Contig0.vcf.gz
		-o /tmp/Contig0_HaplotypeScore.tsv.gz 
		-k HaplotypeScore
		~/NetworkData/vervet/db/genotype_file/method_27/*.vcf.gz
	
"""
import sys, os, math
__doc__ = __doc__%(sys.argv[0])

from palos import ProcessOptions
from palos.io.AbstractMatrixFileWalker import AbstractMatrixFileWalker

ParentClass=AbstractMatrixFileWalker

class ExtractInfoFromVCF(ParentClass):
	__doc__ = __doc__
	option_default_dict = ParentClass.option_default_dict.copy()
	option_default_dict.update({
						('infoKey', 1, ): ['HaplotypeScore', 'k', 1, 'the key of the INFO to be extracted'],\
						}
						)
	option_default_dict[('inputFileFormat', 0, int)][0] = 4
	
	def __init__(self, inputFnameLs=None, **keywords):
		"""
		2011-7-12
		"""
		ParentClass.__init__(self, inputFnameLs=inputFnameLs, **keywords)
	
	
	def processRow(self, row=None, pdata=None):
		"""
		2013.09.05
		"""
		vcfRecord = row
		real_counter = 0
		value = vcfRecord.info_tag2value.get(self.infoKey)
		if value is not None:
			data_row = [vcfRecord.chr, vcfRecord.pos, value]
			self.invariantPData.writer.writerow(data_row)
			real_counter += 1
		return real_counter
	
	def processHeader(self, header=None, pdata=None, rowDefinition=None):
		"""
		2013.09.05
		"""
		header = ['CHROM', 'POS', self.infoKey]
		self.invariantPData.writer.writerow(header)


if __name__ == '__main__':
	main_class = ExtractInfoFromVCF
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(po.arguments, **po.long_option2value)
	instance.run()