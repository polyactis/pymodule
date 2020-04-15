#!/usr/bin/env python3
"""
2012.3.13
	abstract mapper for Sunset mappers/reducers.

"""
from . AbstractDBJob import AbstractDBJob
import sys, os
sys.path.insert(0, os.path.expanduser('~/script'))
sys.path.insert(0, os.path.expanduser('~/src'))
import copy
from palos import ProcessOptions
from palos.db import SunsetDB

class AbstractSunsetMapper(AbstractDBJob):
	__doc__ = __doc__
	option_default_dict = copy.deepcopy(AbstractDBJob.option_default_dict)
	#option_default_dict.pop(('inputFname', 0, ))
	
	def __init__(self, inputFnameLs=None, **keywords):
		"""
		"""
		AbstractDBJob.__init__(self, inputFnameLs=inputFnameLs, **keywords)	#self.connectDB() called within its __init__()
		
	
	def connectDB(self):
		"""
		"""
		self.db_main = SunsetDB.SunsetDB(drivername=self.drivername, db_user=self.db_user, db_passwd=self.db_passwd, \
							hostname=self.hostname, dbname=self.dbname, schema=self.schema, port=self.port)
		self.db_main.setup(create_tables=False)
	
	def checkIfAlignmentListMatchMethodDBEntry(self, individualAlignmentLs=[], methodDBEntry=None, session=None):
		"""
		2013.08.23 moved from AddVCFFile2DB.py
		2012.7.18
		"""
		#make sure methodDBEntry.individual_alignment_ls is identical to individualAlignmentLs
		alignmentIDSetInFile = set([alignment.id for alignment in individualAlignmentLs])
		alignmentIDSetInGenotypeMethod = set([alignment.id for alignment in methodDBEntry.individual_alignment_ls])
		if alignmentIDSetInFile!=alignmentIDSetInGenotypeMethod:
			sys.stderr.write("ERROR: alignmentIDSetInFile (%s) doesn't match alignmentIDSetInFile (%s).\n"%\
							(repr(alignmentIDSetInFile), repr(alignmentIDSetInGenotypeMethod)))
			if session:
				session.rollback()
			#delete all target files if there is any
			self.cleanUpAndExitOnFailure(exitCode=2)
			
if __name__ == '__main__':
    main_class = AbstractSunsetMapper
    po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
    instance = main_class(**po.long_option2value)
    instance.run()
