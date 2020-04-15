#!/usr/bin/env python
"""
2012.3.14
	count the number of reads/bases of fastq or fasta input file.
"""

import sys, os
import csv
from palos import ProcessOptions, PassingData, utils
from palos.mapper.AbstractMapper import AbstractMapper
from palos import ngs

class CountFastqReadBaseCount(AbstractMapper):
	__doc__ = __doc__
	option_default_dict = AbstractMapper.option_default_dict.copy()
	option_default_dict.update({
							('isq_id', 0, int): [0, 'q', 1, 'IndividualSequence.id'],\
							('isqf_id', 0, int): [0, 'f', 1, 'IndividualSequenceFile.id if applicable'],\
							})
	def __init__(self, inputFnameLs, **keywords):
		"""
		"""
		AbstractMapper.__init__(self, **keywords)
		self.inputFnameLs = inputFnameLs
		#2012.3.19 not used
	
	def run(self):
		"""
		"""
		
		if self.debug:
			import pdb
			pdb.set_trace()
		
		baseCountData = ngs.getReadBaseCount(self.inputFname)
		
		writer = csv.writer(open(self.outputFname, 'w'), delimiter='\t')
		header = ['isq_id', 'isqf_id', 'read_count', 'base_count']
		writer.writerow(header)
		data_row = [self.isq_id, self.isqf_id, baseCountData.read_count, baseCountData.base_count]
		writer.writerow(data_row)
		del writer
	
if __name__ == '__main__':
	main_class = CountFastqReadBaseCount
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(po.arguments, **po.long_option2value)
	instance.run()