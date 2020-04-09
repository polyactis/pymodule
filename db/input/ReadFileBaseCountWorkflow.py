#!/usr/bin/env python
"""
Examples:
	# 2012.5.3 run on hoffman2's condorpool, need sshDBTunnel (-H1)
	%s  -i 963-1346 -o dags/ReadCount/read_count_isq_936_1346.xml -u yh --commit -z localhost
		--pegasusFolderName readcount --needSSHDBTunnel
		-l hcondor -j hcondor 
		-t /u/home/eeskin/polyacti/NetworkData/vervet/db -D /u/home/eeskin/polyacti/NetworkData/vervet/db/
	
	# 2012.3.14 
	%s -i 1-864 -o dags/ReadCount/read_count_isq_1_864.xml -u yh -l condorpool -j condorpool -z uclaOffice
		--pegasusFolderName readCount --commit
	
Description:
	2012.3.14
"""
import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.expanduser('~/script'))
sys.path.insert(0, os.path.expanduser('~/src'))

import subprocess, cStringIO, copy
from pymodule import ProcessOptions, getListOutOfStr, PassingData, utils
from pegaflow.DAX3 import File, Link, PFN, Job
from pymodule.pegasus import yh_pegasus
from pymodule.pegasus.AbstractNGSWorkflow import AbstractNGSWorkflow

from pymodule.db import SunsetDB as DBClass
from Sunset.pegasus.AbstractAccuWorkflow import AbstractAccuWorkflow as ParentClass


class ReadFileBaseCountWorkflow(ParentClass):
	__doc__ = __doc__
	option_default_dict = copy.deepcopy(ParentClass.option_default_dict)
	option_default_dict.update({
						('ind_seq_id_ls', 1, ): ['', 'i', 1, 'a comma/dash-separated list of IndividualSequence.id. \
									non-fastq entries will be discarded.', ],\
						})

	def __init__(self,  **keywords):
		"""
		2011-7-11
		"""
		ParentClass.__init__(self, **keywords)
		
		if self.ind_seq_id_ls:
			self.ind_seq_id_ls = getListOutOfStr(self.ind_seq_id_ls, data_type=int)
	
	def registerCustomExecutables(self, workflow=None):
		"""
		2012.3.14
		"""
		if workflow is None:
			workflow = self
		ParentClass.registerCustomExecutables(self, workflow=workflow)
		
		self.addOneExecutableFromPathAndAssignProperClusterSize(
			path=os.path.join(self.pymodulePath, 'mapper/computer/CountFastqReadBaseCount.py'), \
										name='CountFastqReadBaseCount', clusterSizeMultipler=1)
		
		self.addOneExecutableFromPathAndAssignProperClusterSize(
			path=os.path.join(self.thisModulePath, 'db/input/PutReadBaseCountIntoDB.py'), \
										name='PutReadBaseCountIntoDB', clusterSizeMultipler=0.2)
		
	
	def registerISQFiles(self, workflow=None, db_main=None, ind_seq_id_ls=[], local_data_dir='', pegasusFolderName='', \
						input_site_handler='local'):
		"""
		2012.3.14
		"""
		sys.stderr.write("Finding all ISQ-affiliated files of %s ind seq entries ..."%(len(ind_seq_id_ls)))
		returnData = PassingData(jobDataLs=[])
		counter = 0
		Table = DBClass.IndividualSequence
		query = db_main.queryTable(Table).filter(Table.id.in_(ind_seq_id_ls))
		individual_sequence_id_set = set()
		missed_individual_sequence_id_set = set()
		for individual_sequence in query:
			if individual_sequence.individual_sequence_file_ls:	#not empty
				for individual_sequence_file in individual_sequence.individual_sequence_file_ls:
					absPath = os.path.join(local_data_dir, individual_sequence_file.path)
					if os.path.isfile(absPath):
						inputF = File(os.path.join(pegasusFolderName, individual_sequence_file.path))
						inputF.addPFN(PFN("file://" + absPath, input_site_handler))
						inputF.absPath = absPath
						workflow.addFile(inputF)
						returnData.jobDataLs.append(PassingData(output=inputF, jobLs=[], isq_id=individual_sequence.id,\
															isqf_id=individual_sequence_file.id))
						individual_sequence_id_set.add(individual_sequence.id)
					else:
						missed_individual_sequence_id_set.add(individual_sequence.id)
						sys.stderr.write("Warning: IndividualSequenceFile.id=%s (isq-id=%s) doesn't have any affiliated IndividualSequenceFile entries while its path %s is not a file.\n"%\
									(individual_sequence_file.id, individual_sequence.id, absPath))
			elif individual_sequence.path:
				absPath = os.path.join(local_data_dir, individual_sequence.path)
				if os.path.isfile(absPath):
					inputF = File(os.path.join(pegasusFolderName, individual_sequence.path))
					inputF.addPFN(PFN("file://" + absPath, input_site_handler))
					inputF.absPath = absPath
					workflow.addFile(inputF)
					returnData.jobDataLs.append(PassingData(output=inputF, jobLs=[], isq_id=individual_sequence.id,\
														isqf_id=None))
					individual_sequence_id_set.add(individual_sequence.id)
				else:
					sys.stderr.write("Warning: IndividualSequence.id=%s doesn't have any affiliated IndividualSequenceFile entries while its path %s is not a file.\n"%\
									(individual_sequence.id, absPath))
					missed_individual_sequence_id_set.add(individual_sequence.id)
		
		sys.stderr.write(" %s files registered for %s individual_sequence entries. missed %s individual-sequence entries.\n"%\
						(len(returnData.jobDataLs), len(individual_sequence_id_set), len(missed_individual_sequence_id_set)))
		return returnData
	
	def addPutReadBaseCountIntoDBJob(self, workflow=None, executable=None, inputFileLs=[], \
					logFile=None, commit=False, parentJobLs=[], extraDependentInputLs=[], \
					transferOutput=True, extraArguments=None, \
					job_max_memory=10, sshDBTunnel=1, **keywords):
		"""
		20170502 use addGenericFile2DBJob()
		2012.5.3
			add argument sshDBTunnel
		2012.3.14
		"""
		job = self.addGenericFile2DBJob(executable=executable, \
					inputFile=None, inputArgumentOption="-i", \
					outputFile=logFile, outputArgumentOption="--logFilename", inputFileList=inputFileLs, \
					data_dir=None, commit=commit,\
					parentJobLs=parentJobLs, extraDependentInputLs=extraDependentInputLs, \
					extraOutputLs=None, transferOutput=transferOutput, \
					extraArguments=extraArguments, extraArgumentList=None, \
					job_max_memory=job_max_memory,  sshDBTunnel=sshDBTunnel,\
					key2ObjectForJob=None, objectWithDBArguments=self, **keywords)
		return job
	
	
	def addCountFastqReadBaseCountJob(self, workflow=None, executable=None, inputFile=None, \
								outputFile=None, isq_id=None, isqf_id=None, \
					parentJobLs=None, extraDependentInputLs=None, transferOutput=True, extraArguments=None, \
					job_max_memory=100, **keywords):
		"""
		20170503 use addGenericJob()
		2012.3.14
		"""
		job = Job(namespace=workflow.namespace, name=executable.name, version=workflow.version)
		job.addArguments("--inputFname", inputFile, "--outputFname", outputFile)
		if not extraArguments:
			extraArguments = ""
		if isq_id:
			extraArguments += " --isq_id %s "%(isq_id)
		if isqf_id:
			extraArguments += " --isqf_id %s "%(isqf_id)
		
		job = self.addGenericJob(executable=executable, \
					inputFile=inputFile, \
					inputArgumentOption="-i", \
					outputFile=outputFile, outputArgumentOption="-o", \
					parentJobLs=parentJobLs, \
					extraDependentInputLs=extraDependentInputLs, \
					transferOutput=transferOutput, \
					extraArguments=extraArguments, \
					job_max_memory=job_max_memory, \
					**keywords)
		return job
	
	def addJobs(self, workflow=None, inputData=None, pegasusFolderName="", needSSHDBTunnel=0):
		"""
		2012.3.14
		"""
		
		sys.stderr.write("Adding read counting jobs on %s input ..."%(len(inputData.jobDataLs)))
		returnJobData = PassingData()
		
		no_of_jobs = 0
		
		topOutputDir = pegasusFolderName
		if topOutputDir:
			topOutputDirJob = self.addMkDirJob(outputDir=topOutputDir)
			no_of_jobs += 1
		else:
			topOutputDirJob = None
		
		finalReduceFile = File(os.path.join(topOutputDir, 'read_base_count.tsv'))
		
		readBaseCountMergeJob = self.addStatMergeJob(workflow, statMergeProgram=workflow.mergeSameHeaderTablesIntoOne, \
						outputF=finalReduceFile, transferOutput=True, extraArguments=None, parentJobLs=[topOutputDirJob])
		
		logFile = File(os.path.join(topOutputDir, 'PutReadBaseCountIntoDB.log'))
		putCountIntoDBJob = self.addPutReadBaseCountIntoDBJob(workflow, executable=workflow.PutReadBaseCountIntoDB, \
					inputFileLs=[finalReduceFile], \
					logFile=logFile, commit=self.commit, parentJobLs=[readBaseCountMergeJob], \
					extraDependentInputLs=[], transferOutput=True, \
					extraArguments=None, \
					job_max_memory=10, sshDBTunnel=needSSHDBTunnel)
		no_of_jobs += 2
		for jobData in inputData.jobDataLs:
			#add the read count job
			outputFile = File(os.path.join(topOutputDir, 'read_count_isq_%s_isqf_%s.tsv'%(jobData.isq_id, jobData.isqf_id)))
			readCountJob = self.addCountFastqReadBaseCountJob(workflow, executable=workflow.CountFastqReadBaseCount, \
								inputFile=jobData.output, outputFile=outputFile, isq_id=jobData.isq_id, \
								isqf_id=jobData.isqf_id, \
								parentJobLs=jobData.jobLs + [topOutputDirJob], extraDependentInputLs=None, \
								transferOutput=False, extraArguments=None, \
								job_max_memory=10, no_of_cpus=4)
			
			no_of_jobs += 1
			self.addInputToStatMergeJob(workflow, statMergeJob=readBaseCountMergeJob, \
								inputF=readCountJob.output, parentJobLs=[readCountJob])
			
		sys.stderr.write("%s jobs.\n"%(no_of_jobs))
		return putCountIntoDBJob
	
	def setup_run(self):
		"""
		2013.04.07 wrap all standard pre-run() related functions into this function.
			setting up for run(), called by run()
		"""
		pdata = AbstractNGSWorkflow.setup_run(self)
		workflow = pdata.workflow
		
		db_main = self.db_main
		session = db_main.session
		session.begin(subtransactions=True)
		"""
		Traceback (most recent call last):
		  File "/u/home/eeskin/polyacti/script/vervet/src/db/ReadFileBaseCountWorkflow.py", line 249, in <module>
		    instance.run()
		  File "/u/home/eeskin/polyacti/script/vervet/src/db/ReadFileBaseCountWorkflow.py", line 232, in run
		    pdata = self.setup_run()
		  File "/u/home/eeskin/polyacti/script/vervet/src/db/ReadFileBaseCountWorkflow.py", line 217, in setup_run
		    session.begin()
		  File "/u/home/eeskin/polyacti/lib/python/sqlalchemy/orm/scoping.py", line 139, in do
		    return getattr(self.registry(), name)(*args, **kwargs)
		  File "/u/home/eeskin/polyacti/lib/python/sqlalchemy/orm/session.py", line 550, in begin
		    "A transaction is already begun.  Use subtransactions=True "
		sqlalchemy.exc.InvalidRequestError: A transaction is already begun.  Use subtransactions=True to allow subtransactions.

		"""
		
		inputData = self.registerISQFiles(workflow=workflow, db_main=db_main, ind_seq_id_ls=self.ind_seq_id_ls, \
										local_data_dir=self.local_data_dir, pegasusFolderName=self.pegasusFolderName,\
										input_site_handler=self.input_site_handler)
		
		registerReferenceData = self.getReferenceSequence()
		
		
		return PassingData(workflow=workflow, inputData=inputData,\
						registerReferenceData=registerReferenceData)
	def run(self):
		"""
		2011-7-11
		"""
		pdata = self.setup_run()
		workflow = pdata.workflow
		
		inputData=pdata.inputData
		
		self.addJobs(workflow, inputData=inputData, pegasusFolderName=self.pegasusFolderName,
					needSSHDBTunnel=self.needSSHDBTunnel)
		
		self.end_run()

if __name__ == '__main__':
	main_class = ReadFileBaseCountWorkflow
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()