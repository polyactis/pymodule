#!/usr/bin/env python3
"""
2013.2.10 This program replaces the ms command path in the inputFname (output of history2ms.pl) with correct msPath
	so that the inputFname could be run as a shell script.   

Examples:
	#2013.2.11 add --replaceTheHengLiOutputFlagAsWell to get rid of "-l" argument so that msHOT or ms could run as well
	%s -i 1534_788_2009098_GA_vs_524.ms_command.sh -o 1534_788_2009098_GA_vs_524.msHOT-lite.output.traditional_output.sh
		--msPath ~/script/lh3_foreign/msHOT-lite/msHOT-lite
		#--replaceTheHengLiOutputFlagAsWell

"""

import sys, os, math
__doc__ = __doc__%(sys.argv[0])
import re
from palos import ProcessOptions, utils
from palos.mapper.AbstractMapper import AbstractMapper

class ReplaceMSPathInMSCommandFile(AbstractMapper):
	__doc__ = __doc__
	option_default_dict = AbstractMapper.option_default_dict.copy()
	#option_default_dict.pop(('inputFname', 1, ))
	option_default_dict.update({
		('oldMSPath', 0, ): ['msHOT-lite', '', 1, 'path of the ms program in inputFname '],\
		('msPath', 0, ): [None, '', 1, 'path to the ms or msHOT, msHOT-lite program, '],\
		('replaceTheHengLiOutputFlagAsWell', 0, int): [0, '', 0, "Heng Li's msHOT-lite has a '-l' flag to output in a succinct format."
			" Toggle this to get rid of '-l'. "],\
		})
	def __init__(self, inputFnameLs=None, **keywords):
		"""
		"""
		AbstractMapper.__init__(self, inputFnameLs=inputFnameLs, **keywords)	#self.connectDB() called within its __init__()
			
	def run(self):
		"""
		"""
		if self.debug:
			import pdb
			pdb.set_trace()
		
		if not os.path.isfile(self.inputFname):
			sys.stderr.write("Error: file, %s,  is not a file.\n"%(self.inputFname))
			sys.exit(3)
			
		inf = utils.openGzipFile(self.inputFname, 'r')
		outf = open(self.outputFname, 'w')
		for line in inf:
			newLine = re.sub(r'%s'%(self.oldMSPath), r'%s'%(self.msPath), line)
			if self.replaceTheHengLiOutputFlagAsWell:
				newLine = newLine.replace(" -l", "")	#it's global and exhaustive, any " -l " will be replaced.
			outf.write(newLine)
		inf.close()
		outf.close()

	
if __name__ == '__main__':
	main_class = ReplaceMSPathInMSCommandFile
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()