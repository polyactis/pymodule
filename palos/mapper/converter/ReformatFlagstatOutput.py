#!/usr/bin/env python
"""
Examples:
	%s 
	
	%s -i /tmp/outputStat.txt.gz -a 1043 -o /tmp/outputStat.tsv.gz

Description:
	2012.4.3
		reformat output of samtools flagstat so that it could be inserted into db.
		both input and output files could be either gzipped or not.

	
"""

import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.expanduser('~/script'))

import csv, copy, re
from pymodule import ProcessOptions, getListOutOfStr, PassingData, utils
from pymodule.mapper.AbstractMapper import AbstractMapper as ParentClass

class ReformatFlagstatOutput(ParentClass):
	__doc__ = __doc__
	option_default_dict = copy.deepcopy(ParentClass.option_default_dict)
	option_default_dict.pop(('outputFnamePrefix', 0, ))
	option_default_dict.update({
						('alignmentID', 1, ): [None, 'a', 1, 'ID of this alignment from which all the stats are extracted.'],\
						
							})
	def __init__(self,  **keywords):
		"""
		"""
		ParentClass.__init__(self, **keywords)
	
	def getNumberOutOfFlagStatLine(self, line=None, grabPattern=None):
		"""
		20170602 added argument grabPattern
		2012.4.3
		"""
		searchResult = grabPattern.search(line)
		if searchResult:
			n1 = int(searchResult.group(1))
			n2 = int(searchResult.group(2))
			return n1+n2
		else:
			sys.stderr.write("Error: could not parse numbers out of this line (%s).\n"%(line))
			sys.exit(2)
			return None

	def run(self):
		"""
		"""
		
		if self.debug:
			import pdb
			pdb.set_trace()
		"""
		2012.4.3
		the output of samtools flagstat looks like:

20170602 new flagstat output

470131994 + 0 in total (QC-passed reads + QC-failed reads)
63918054 + 0 secondary
0 + 0 supplementary
3001858 + 0 duplicates
460732266 + 0 mapped (98.00% : N/A)
406213940 + 0 paired in sequencing
203106970 + 0 read1
203106970 + 0 read2
391157952 + 0 properly paired (96.29% : N/A)
394571382 + 0 with itself and mate mapped
2242830 + 0 singletons (0.55% : N/A)
2443798 + 0 with mate mapped to a different chr
1751451 + 0 with mate mapped to a different chr (mapQ>=5)

		"""
		
		inf = utils.openGzipFile(self.inputFname, openMode='r')
		writer = csv.writer(utils.openGzipFile(self.outputFname, openMode='w'), delimiter='\t')
		header = ['alignmentID', 'total_no_of_reads', 'perc_secondary', 'perc_supplementary', \
				'perc_reads_mapped', 'perc_duplicates', 'perc_paired', 'perc_properly_paired', \
				'perc_both_mates_mapped', 'perc_singletons',\
				'perc_mapped_to_diff_chrs', 'perc_mapq5_mapped_to_diff_chrs']
		writer.writerow(header)
		
		#float total_no_of_reads now so that no "float" upon division
		total_no_of_reads = float(self.getNumberOutOfFlagStatLine(line=inf.next(), \
									grabPattern=re.compile(r'^(\d+) \+ (\d+) in total')))
		no_of_secondary = self.getNumberOutOfFlagStatLine(line=inf.next(),\
									grabPattern=re.compile(r'^(\d+) \+ (\d+) secondary'))
		no_of_supplementary = self.getNumberOutOfFlagStatLine(line=inf.next(),\
									grabPattern=re.compile(r'^(\d+) \+ (\d+) supplementary'))
		no_of_duplicates = self.getNumberOutOfFlagStatLine(line=inf.next(),\
									grabPattern=re.compile(r'^(\d+) \+ (\d+) duplicates'))
		no_of_mapped = self.getNumberOutOfFlagStatLine(line=inf.next(),\
									grabPattern=re.compile(r'^(\d+) \+ (\d+) mapped'))
		no_of_paired = self.getNumberOutOfFlagStatLine(line=inf.next(),\
									grabPattern=re.compile(r'^(\d+) \+ (\d+) paired in sequencing'))
		no_of_read1 = self.getNumberOutOfFlagStatLine(line=inf.next(),\
									grabPattern=re.compile(r'^(\d+) \+ (\d+) read1'))
		no_of_read2 = self.getNumberOutOfFlagStatLine(line=inf.next(),\
									grabPattern=re.compile(r'^(\d+) \+ (\d+) read2'))
		no_of_properly_paired = self.getNumberOutOfFlagStatLine(line=inf.next(),\
									grabPattern=re.compile(r'^(\d+) \+ (\d+) properly paired'))
		no_of_both_mates_mapped = self.getNumberOutOfFlagStatLine(line=inf.next(),\
									grabPattern=re.compile(r'^(\d+) \+ (\d+) with itself and mate mapped'))
		no_of_singletons = self.getNumberOutOfFlagStatLine(line=inf.next(),\
									grabPattern=re.compile(r'^(\d+) \+ (\d+) singletons'))
		no_of_mates_mapped_to_diff_chrs = self.getNumberOutOfFlagStatLine(line=inf.next(),\
									grabPattern=re.compile(r'^(\d+) \+ (\d+) with mate mapped to a different chr\n'))
		no_of_mates_mapped_to_diff_chrs_mapQAbove5 = self.getNumberOutOfFlagStatLine(line=inf.next(),\
									grabPattern=re.compile(r'^(\d+) \+ (\d+) with mate mapped to a different chr \(mapQ>=5\)'))
		#
		del inf
		
		data_row = [self.alignmentID, total_no_of_reads, no_of_secondary/total_no_of_reads*100, \
				no_of_supplementary/total_no_of_reads*100, \
				no_of_mapped/total_no_of_reads*100, no_of_duplicates/total_no_of_reads*100,\
				no_of_paired/total_no_of_reads*100, no_of_properly_paired/total_no_of_reads*100,\
				no_of_both_mates_mapped/total_no_of_reads*100, no_of_singletons/total_no_of_reads*100,\
				no_of_mates_mapped_to_diff_chrs/total_no_of_reads*100,
				no_of_mates_mapped_to_diff_chrs_mapQAbove5/total_no_of_reads*100]
		writer.writerow(data_row)
		del writer
		

if __name__ == '__main__':
	main_class = ReformatFlagstatOutput
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()