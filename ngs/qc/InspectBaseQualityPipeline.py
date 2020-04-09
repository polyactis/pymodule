#!/usr/bin/env python
"""
Examples:
	
	# 2011-8-16 on condorpool
	%s -o inspectBaseQuality.xml -u yh -i 1-8,15-130
	
	# 2011-8-16 use hoffman2 site_handler
	%s -o inspectBaseQuality.xml -u yh -i 1-8,15-130 
		-l hoffman2 -e /u/home/eeskin/polyacti -t /u/home/eeskin/polyacti/NetworkData/vervet/db
	
	
Description:
	2011-8-16
		construct a pegasus workflow to run InspectBaseQuality.py over a list of individual sequences
"""
import sys, os, copy
__doc__ = __doc__%(sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.expanduser('~/script'))
sys.path.insert(0, os.path.expanduser('~/src'))

from pymodule import ProcessOptions, getListOutOfStr, PassingData, utils, yh_pegasus
from pymodule.db import SunsetDB as DBClass
from Sunset.pegasus.AbstractAccuWorkflow import AbstractAccuWorkflow as ParentClass

class InspectBaseQualityPipeline(ParentClass):
	__doc__ = __doc__
	
	commonOptionDict = copy.deepcopy(ParentClass.option_default_dict)
	#commonOptionDict.pop(('inputDir', 0, ))

	option_default_dict = copy.deepcopy(commonOptionDict)
	option_default_dict.update({
				('ind_seq_id_ls', 1, ): ['', 'i', 1, 'a comma/dash-separated list of IndividualSequence.id. non-fastq entries will be discarded.', ],\
						})

	def __init__(self, **keywords):
		"""
		2011-7-11
		"""
		ParentClass.__init__(self, **keywords)
		if self.ind_seq_id_ls:
			self.ind_seq_id_ls = getListOutOfStr(self.ind_seq_id_ls, data_type=int)

	def registerCustomExecutables(self, workflow=None):
		"""
		"""
		self.addOneExecutableFromPathAndAssignProperClusterSize(
			path=os.path.join(self.thisModulePath, 'mapper/InspectBaseQuality.py'), \
			name='InspectBaseQuality', clusterSizeMultipler=1)
		
			
	def run(self):
		"""
		2011-7-11
		"""
		if self.debug:
			import pdb
			pdb.set_trace()
		
		db_main = self.db_main
		session = db_main.session
		session.begin()
		
		if not self.data_dir:
			self.data_dir = db_main.data_dir
		
		if not self.local_data_dir:
			self.local_data_dir = db_main.data_dir
		
		workflow = self.initiateWorkflow()
		
		self.registerJars()
		self.registerExecutables()
		self.registerCustomExecutables(workflow=workflow)
		
		
		#must use db_vervet.data_dir.
		# If self.data_dir differs from db_vervet.data_dir, this program (must be run on submission host) won't find files.
		individualSequenceID2FilePairLs = db_main.getIndividualSequenceID2FilePairLs(self.ind_seq_id_ls, data_dir=self.data_dir)
		
		for ind_seq_id, FilePairLs in individualSequenceID2FilePairLs.iteritems():
			individual_sequence = db_main.queryTable(DBClass.IndividualSequence).get(ind_seq_id)
			if individual_sequence is not None and individual_sequence.format=='fastq':
				#start to collect all files affiliated with this individual_sequence record 
				inputFilepathLs = []
				for filePair in FilePairLs:
					for fileRecord in filePair:
						relativePath, format, sequence_type = fileRecord[:3]
						filepath = os.path.join(self.data_dir, relativePath)
						inputFilepathLs.append(filepath)
				
				#create jobs
				for filepath in inputFilepathLs:
					prefix, suffix = utils.getRealPrefixSuffixOfFilenameWithVariableSuffix(filepath)
					if suffix=='.fastq':	#sometimes other files get in the way
						inspectBaseQuality_job = self.addGenericDBJob(executable=self.InspectBaseQuality, \
												inputFile=None, \
												inputArgumentOption='-i',\
												outputFile=None, outputArgumentOption=None,\
												parentJobLs=None, extraDependentInputLs=None, \
												extraOutputLs=None, extraArguments=None, \
												transferOutput=False, \
												extraArgumentList=['-i', filepath, '--read_sampling_rate', '0.005', \
														'--quality_score_format', individual_sequence.quality_score_format], \
												objectWithDBArguments=self,\
												job_max_memory=20000, \
												walltime=120)
		
		sys.stderr.write("%s jobs.\n"%(self.no_of_jobs))
		# Write the DAX to stdout
		outf = open(self.outputFname, 'w')
		self.writeXML(outf)
		outf.close()
		if self.commit:
			session.commit()
		else:
			session.rollback()
	
if __name__ == '__main__':
	main_class = InspectBaseQualityPipeline
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()
