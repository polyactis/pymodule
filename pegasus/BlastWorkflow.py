#!/usr/bin/env python
"""
Examples:

	%s -d /Network/Data/vervet/db/individual_sequence/524_superContigsMinSize2000.fasta 
		-i ~/script/vervet/data/OphoffMethylation/DMR330K_ProbeSeq.fasta  -a 2 -l condorpool -j condorpool 
		-C 1 -o workflow/BlastDMR330K_ProbeSeqAgainst524.xml

2011-11-22
	a common class for pegasus workflows that work on NGS (next-gen sequencing) data
"""
import sys, os, math

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

from pymodule import ProcessOptions, getListOutOfStr, PassingData, yh_pegasus, utils
from Pegasus.DAX3 import *
from AbstractWorkflow import AbstractWorkflow

class BlastWorkflow(AbstractWorkflow):
	__doc__ = __doc__
	option_default_dict = AbstractWorkflow.option_default_dict.copy()
	option_default_dict.update({
						("inputFname", 1, ): ["", 'i', 1, 'the input fasta file'],\
						("formatdbPath", 1, ): ["/usr/bin/formatdb", '', 1, 'path to formatdb, index fasta database file'],\
						("blastallPath", 1, ): ["/usr/bin/blastall", '', 1, 'path to blastall'],\
						("blockSize", 1, int): [1000, '', 1, 'how many sequences each blast job handles'],\
						('databaseFname', 1, ): ['', 'd', 1, 'filename of the database to blast against, must be indexed', ],\
						('minNoOfIdentities', 0, int): [None, 'm', 1, 'minimum number of identities between a query and target', ],\
						('maxNoOfMismatches', 0, int): [None, 'a', 1, 'minimum number of mismatches between a query and target', ],\
						('minIdentityPercentage', 0, float): [None, 'n', 1, 'minimum percentage of identities between a query and target', ],\
						
						})
						#('bamListFname', 1, ): ['/tmp/bamFileList.txt', 'L', 1, 'The file contains path to each bam file, one file per line.'],\

	def __init__(self,  **keywords):
		"""
		2012.5.23
		"""
		AbstractWorkflow.__init__(self, **keywords)
	
	
	def registerBlastNucleotideDatabaseFile(self, ntDatabaseFname=None, input_site_handler='local'):
		"""
		2012.5.23
		"""
		return yh_pegasus.registerRefFastaFile(self, ntDatabaseFname, registerAffiliateFiles=True, \
									input_site_handler=input_site_handler,\
									checkAffiliateFileExistence=True, addPicardDictFile=False, \
									affiliateFilenameSuffixLs=['nin', 'nhr', 'nsq'])
	
	
	def getNoOfSequencesFromFasta(self, inputFastaFname=None):
		"""
		2012.5.24
		"""
		sys.stderr.write("Getting number of sequences from %s ..."%(inputFastaFname))
		inf = utils.openGzipFile(inputFastaFname)
		no_of_sequences = 0
		for line in inf:
			if line[0]=='>':
				no_of_sequences += 1
		del inf
		sys.stderr.write("%s sequences.\n"%(no_of_sequences))
		return no_of_sequences
	
	def addSplitFastaFileJob(self, executable=None, inputFile=None, outputFnamePrefix=None, \
						noOfSequencesPerSplitFile=1000, filenameSuffix="", noOfTotalSequences=1000000,\
						parentJobLs=[], extraDependentInputLs=[], transferOutput=False, \
						extraArguments=None, job_max_memory=500, **keywords):
		"""
		2012.5.24
		"""
		job = Job(namespace=self.namespace, name=executable.name, version=self.version)
		
		noOfSplitFiles = int(math.ceil(noOfTotalSequences/float(noOfSequencesPerSplitFile)))
		suffixLength = len(repr(noOfSplitFiles))
		
		job.addArguments("-i", inputFile, "-l %s"%(noOfSequencesPerSplitFile), \
						"-O", outputFnamePrefix, '-f %s'%(filenameSuffix), '-a %s'%(suffixLength))
		if extraArguments:
			job.addArguments(extraArguments)
		job.uses(inputFile, transfer=True, register=True, link=Link.INPUT)
		job.outputList = []
		for i in xrange(noOfSplitFiles):	#start from 0
			splitFname = utils.comeUpSplitFilename(outputFnamePrefix=outputFnamePrefix, suffixLength=suffixLength, fileOrder=i,\
											filenameSuffix=filenameSuffix)
			splitFile = File(splitFname)
			
			job.outputList.append(splitFile)
			job.uses(splitFile, transfer=transferOutput, register=True, link=Link.OUTPUT)
			
		yh_pegasus.setJobProperRequirement(job, job_max_memory=job_max_memory)
		self.addJob(job)
		for parentJob in parentJobLs:
			if parentJob:
				self.depends(parent=parentJob, child=job)
		for input in extraDependentInputLs:
			if input:
				job.uses(input, transfer=True, register=True, link=Link.INPUT)
		return job
	
	def addBlastWrapperJob(self, executable=None, inputFile=None, outputFile=None, databaseFile=None,\
						maxNoOfMismatches=None, minNoOfIdentities=None, minIdentityPercentage=None, blastallPath=None, \
						parentJobLs=[], extraDependentInputLs=[], transferOutput=False, \
						extraArguments=None, job_max_memory=2000, **keywords):
		"""
		2012.5.24
		"""
		extraArgumentList = ['-l %s'%blastallPath, '-d', databaseFile]
		if maxNoOfMismatches:
			extraArgumentList.append('-a %s'%maxNoOfMismatches)
		if minNoOfIdentities:
			extraArgumentList.append("-m %s"%(minNoOfIdentities))
		if minIdentityPercentage:
			extraArgumentList.append("-n %s"%(minIdentityPercentage))
		if extraArguments:
			extraArgumentList.append(extraArguments)
		return self.addGenericJob(executable=executable, inputFile=inputFile, outputFile=outputFile, \
						parentJobLs=parentJobLs, extraDependentInputLs=extraDependentInputLs, transferOutput=transferOutput, \
						extraArgumentList=extraArgumentList, job_max_memory=job_max_memory)
		
	def addJobs(self, inputData=None, outputDirPrefix="", ntDatabaseFileList=None, noOfTotalSequences=None, \
			transferOutput=True, makeBlastDBJob=None):
		"""
		2012.5.24
		"""
		
		sys.stderr.write("Adding blast jobs for %s input ... "%(len(inputData.jobDataLs)))
		no_of_jobs= 0
		
		topOutputDir = "%sBlast"%(outputDirPrefix)
		topOutputDirJob = yh_pegasus.addMkDirJob(self, mkdir=self.mkdirWrap, outputDir=topOutputDir)
		no_of_jobs += 1
		
		allBlastResultFile = File(os.path.join(topOutputDir, 'blast.tsv'))
		allBlastMergeJob = self.addStatMergeJob(statMergeProgram=self.mergeSameHeaderTablesIntoOne, \
							outputF=allBlastResultFile, transferOutput=transferOutput, parentJobLs=[topOutputDirJob])
		no_of_jobs += 1
		
		ntDatabaseFile = ntDatabaseFileList[0]
		returnData = PassingData()
		returnData.jobDataLs = []
		
		for jobData in inputData.jobDataLs:
			inputF = jobData.output
			outputFnamePrefix = os.path.join(topOutputDir, os.path.splitext(os.path.basename(inputF.name))[0])
			
			splitFastaJob = self.addSplitFastaFileJob(executable=self.SplitFastaFile, inputFile=inputF, outputFnamePrefix=outputFnamePrefix, \
						noOfSequencesPerSplitFile=self.blockSize, filenameSuffix=".fasta", noOfTotalSequences=noOfTotalSequences,\
						parentJobLs=jobData.jobLs + [topOutputDirJob], extraDependentInputLs=[], transferOutput=False, \
						extraArguments=None, job_max_memory=500)
			no_of_jobs += 1
			for splitFastaOutput in splitFastaJob.outputList:
				outputFile = File('%s.tsv'%(splitFastaOutput.name))
				blastJob = self.addBlastWrapperJob(executable=self.BlastWrapper, inputFile=splitFastaOutput, outputFile=outputFile, \
								databaseFile=ntDatabaseFile,\
								maxNoOfMismatches=self.maxNoOfMismatches, minNoOfIdentities=self.minNoOfIdentities, \
								minIdentityPercentage=self.minIdentityPercentage, blastallPath=self.blastallPath, \
								parentJobLs=[splitFastaJob, makeBlastDBJob], extraDependentInputLs=ntDatabaseFileList, transferOutput=False, \
								extraArguments=None, job_max_memory=1000)
				
				#add output to some reduce job
				self.addInputToStatMergeJob(statMergeJob=allBlastMergeJob, \
								inputF=blastJob.output, parentJobLs=[blastJob])
				no_of_jobs += 1
		sys.stderr.write("%s jobs. Done.\n"%(no_of_jobs))
		#include the tfam (outputList[1]) into the fileList
		returnData.jobDataLs.append(PassingData(jobLs=[allBlastMergeJob], file=allBlastResultFile, \
											fileList=[allBlastResultFile]))
		return returnData
	
				
		
	def registerCustomExecutables(self):
		"""
		2012.5.23
		"""
		namespace = self.namespace
		version = self.version
		operatingSystem = self.operatingSystem
		architecture = self.architecture
		clusters_size = self.clusters_size
		site_handler = self.site_handler
		
		executableList = []
		blastall = Executable(namespace=namespace, name="blastall", version=version, os=operatingSystem, \
							arch=architecture, installed=True)
		blastall.addPFN(PFN("file://" + os.path.join(self.blastallPath), site_handler))
		executableList.append(blastall)
		
		formatdb = Executable(namespace=namespace, name="formatdb", version=version, os=operatingSystem, \
							arch=architecture, installed=True)
		formatdb.addPFN(PFN("file://" + os.path.join(self.formatdbPath), site_handler))
		executableList.append(formatdb)
		
		BlastWrapper = Executable(namespace=namespace, name="BlastWrapper", version=version, os=operatingSystem, \
							arch=architecture, installed=True)
		BlastWrapper.addPFN(PFN("file://" + os.path.join(self.pymodulePath, 'pegasus/mapper/BlastWrapper.py'), site_handler))
		executableList.append(BlastWrapper)
		
		SplitFastaFile = Executable(namespace=namespace, name="SplitFastaFile", version=version, os=operatingSystem, \
							arch=architecture, installed=True)
		SplitFastaFile.addPFN(PFN("file://" + os.path.join(self.pymodulePath, 'pegasus/mapper/SplitFastaFile.py'), site_handler))
		executableList.append(SplitFastaFile)
		
		for executable in executableList:
			executable.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%self.clusters_size))
			self.addExecutable(executable)
			setattr(self, executable.name, executable)
	
	def addMakeBlastDBJob(self, executable=None, inputFile=None, \
						parentJobLs=[], extraDependentInputLs=[], transferOutput=False, \
						extraArguments=None, job_max_memory=500, **keywords):
		"""
		2012.5.24
			untested
		"""
		job = Job(namespace=self.namespace, name=executable.name, version=self.version)
		job.addArguments("-i", inputFile, "-p F")
		if extraArguments:
			job.addArguments(extraArguments)
		job.uses(inputFile, transfer=True, register=True, link=Link.INPUT)
		job.outputList = []
		for suffix in ['.nin', '.nhr', '.nsq']:	#start from 0
			dbIndexFile = File('%s%s'%(inputFile.name, suffix))
			job.uses(dbIndexFile, transfer=transferOutput, register=True, link=Link.OUTPUT)
			job.outputList.append(dbIndexFile)
			
		yh_pegasus.setJobProperRequirement(job, job_max_memory=job_max_memory)
		self.addJob(job)
		for parentJob in parentJobLs:
			if parentJob:
				self.depends(parent=parentJob, child=job)
		for input in extraDependentInputLs:
			if input:
				job.uses(input, transfer=True, register=True, link=Link.INPUT)
		return job
	
	def run(self):
		"""
		2011-7-11
		"""
		
		if self.debug:
			import pdb
			pdb.set_trace()
		
		self.registerJars()
		self.registerExecutables()
		self.registerCustomExecutables()
		
		
		inputData = PassingData(jobDataLs = [])
		inputFile = self.registerOneInputFile(inputFname=self.inputFname, folderName=self.pegasusFolderName)
		inputData.jobDataLs.append(PassingData(output=inputFile, jobLs=[]))
		noOfTotalSequences= self.getNoOfSequencesFromFasta(inputFastaFname=self.inputFname)
		
		ntDatabaseFileList = self.registerBlastNucleotideDatabaseFile(ntDatabaseFname=self.databaseFname, \
																	input_site_handler=self.input_site_handler)
		ntDatabaseFile = ntDatabaseFileList[0]

		if len(ntDatabaseFileList)<4:	#some nt-database index file is missing
			sys.stderr.write("Adding blast-db-making job...")
			makeBlastDBJob = self.addMakeBlastDBJob(executable=self.formatdb,\
												inputFile=ntDatabaseFile, transferOutput=True)
			#add the index files to the ntDatabaseFileList
			ntDatabaseFileList = [ntDatabaseFile] + makeBlastDBJob.outputList
			sys.stderr.write(".\n")
		else:
			makeBlastDBJob = None
		
		self.addJobs(inputData=inputData, outputDirPrefix=self.pegasusFolderName, ntDatabaseFileList=ntDatabaseFileList, \
					noOfTotalSequences=noOfTotalSequences, \
					transferOutput=True, makeBlastDBJob=makeBlastDBJob)
		# Write the DAX to stdout
		outf = open(self.outputFname, 'w')
		self.writeXML(outf)

if __name__ == '__main__':
	main_class = BlastWorkflow
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()