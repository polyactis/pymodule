#!/usr/bin/env python3
"""
Examples:
	# run the program on crocea and output a local condor workflow
	%s -i ~/NetworkData/vervet/VRC/ -t /u/home/eeskintmp/polyacti/NetworkData/vervet/db/ 
		--bamFname2MonkeyIDMapFname ~/script/vervet/data/VRC_sequencing_20110802.tsv
		-u yh --commit -z dl324b-1.cmb.usc.edu -o /tmp/condorpool.xml

	#2011-8-26	generate a list of all bam file physical paths through find. (doing this because they are not located on crocea)
	find /u/home/eeskintmp/polyacti/NetworkData/vervet/raw_sequence/ -name *.bam  > /u/home/eeskintmp/polyacti/NetworkData/vervet/raw_sequence/bamFileList.txt
	# run the program on the crocea and output a hoffman2 workflow. (with db commit)
	%s -i ~/mnt/hoffman2/u/home/eeskintmp/polyacti/NetworkData/vervet/raw_sequence/bamFileList.txt
		-e /u/home/eeskin/polyacti/
		-t /u/home/eeskin/polyacti/NetworkData/vervet/db/ -u yh
		--bamFname2MonkeyIDMapFname ~/mnt/hoffman2/u/home/eeskintmp/polyacti/NetworkData/vervet/raw_sequence/xfer.genome.wustl.edu/gxfer3/46019922623327/Vervet_12_4X_README.tsv
		-z dl324b-1.cmb.usc.edu -l hoffman2
		-o unpackAndAdd12_2007Monkeys2DB_hoffman2.xml
		--commit 
	
	# 2011-8-28 output a workflow to run on local condorpool, no db commit (because records are already in db)
	%s -i /Network/Data/vervet/raw_sequence/xfer.genome.wustl.edu/gxfer3/
		--bamFname2MonkeyIDMapFname ~/mnt/hoffman2/u/home/eeskintmp/polyacti/NetworkData/vervet/raw_sequence/xfer.genome.wustl.edu/gxfer3/46019922623327/Vervet_12_4X_README.tsv
		 --minNoOfReads 4000000 -u yh --commit
		-z dl324b-1.cmb.usc.edu -j condorpool -l condorpool -o dags/AddReads2DB/DownloadUnpackReads/unpackAndAdd12_2007Monkeys2DB_condor.xml
	
	# 2012.4.30 run on hcondor, to import McGill 1X data (-y2), (-e) is not necessary because it's running on hoffman2 and can recognize home folder.
	#. --needSSHDBTunnel means it needs sshTunnel for db-interacting jobs.
	%s -i ~/NetworkData/vervet/raw_sequence/McGill96_1X/ -z localhost -u yh -j hcondor -l hcondor --commit
		-o dags/AddReads2DB/unpackMcGill96_1X.xml -y2 --needSSHDBTunnel 
		-D /u/home/eeskin/polyacti/NetworkData/vervet/db/ -t /u/home/eeskin/polyacti/NetworkData/vervet/db/
		-e /u/home/eeskin/polyacti
	
	# 2012.6.2 add 18 south-african monkeys RNA read data (-y3), sequenced by Joe DeYoung's core (from Nam's folder),
	# later manually changed its tissue id to distinguish them from DNA data (below) 
	%s -i ~namtran/panasas/Data/HiSeqRaw/Ania/SIVpilot/by.Charles.Demultiplexed/ -z localhost -u yh -j hcondor -l hcondor
		--commit -o dags/AddReads2DB/unpack20SouthAfricaSIVmonkeys.xml -y3 --needSSHDBTunnel
		-D /u/home/eeskin/polyacti/NetworkData/vervet/db/
		-t /u/home/eeskin/polyacti/NetworkData/vervet/db/ -e /u/home/eeskin/polyacti 
		--bamFname2MonkeyIDMapFname ~namtran/panasas/Data/HiSeqRaw/Ania/SIVpilot/by.Charles.Demultiplexed/sampleIds.txt
		--minNoOfReads 4000000
	
	# 2012.6.3 add 24 south-african monkeys DNA read data (-y4), sequenced by Joe DeYoung's core (from Nam's folder)
	# --minNoOfReads 4000000
	%s -i ~namtran/panasas/Data/HiSeqRaw/Ania/SIVpilot/LowpassWGS/Demultiplexed/
		-z localhost -u yh -j hcondor -l hcondor --commit -o dags/AddReads2DB/unpack20SouthAfricaSIVmonkeysDNA.xml
		-y4 --needSSHDBTunnel
		-D /u/home/eeskin/polyacti/NetworkData/vervet/db/ -t /u/home/eeskin/polyacti/NetworkData/vervet/db/
		-e /u/home/eeskin/polyacti/
		--minNoOfReads 4000000
		
	# 2017.04.28 added TCGA sequences (.bam) into db
	%s -i /y/Sunset/tcga/HNSC_TCGA/ --hostname pdc -u huangyu -j ycondor -l ycondor 
		-o dags/unpackTCGAHNSCSamples.xml -y6
		--tissueSourceSiteFname /y/Sunset/tcga/tcga_code_tables/tissueSourceSite.tsv 
		--minNoOfReads 8000000 --dbname pmdb -k xiandao --ref_ind_seq_id 1 --commit
	
Description:
	2011-8-2
		input:
			input (if directory, recursively go and find all bam files; if file, each line is path to bam),
			tsv file (map bam filename to monkey ID),
		
		this program
			1. queries db if individual record for that monkey is already in db, or not. create a new entry in db if not
			2. queries db if individual_sequence record for that monkey is already in db. create one new entry if not.
				update individual_sequence.path if it's none
			3. if commit, the db action would be committed
			4. outputs the whole unpack workflow to an xml output file.
		
		This program has to be run on the pegasus submission host.
		option "--commit" commits the db transaction.
	
	The bamFname2MonkeyIDMapFname contains a map from monkey ID to the bam filename (only base filename is used, relative directory is not used).
	Example ("Library" and "Bam Path" are required):
		FlowCell	Lane	Index Sequence	Library	Common Name	Bam Path	MD5
		64J6AAAXX	1	TGACCA	VCAC-2007002-1-lib1	African Green Monkey	gerald_64J6AAAXX_1.bam	gerald_64J6AAAXX_1.bam.md5
	
	1. Be careful with the db connection setting as it'll be passed to the db-registration job.
		Make sure all computing nodes have access to the db.
	2. The workflow has to be run on nodes where they have direct db and db-affiliated file-storage access.
"""
import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.expanduser('~/script'))
sys.path.insert(0, os.path.expanduser('~/src'))

import copy, re, csv
from pymodule import ProcessOptions, PassingData, MatrixFile, utils
from pymodule.pegasus import yh_pegasus
from pegaflow.DAX3 import File
from . AbstractAccuWorkflow import AbstractAccuWorkflow

class ImportIndividualSequence2DB(AbstractAccuWorkflow):
	__doc__ = __doc__
	option_default_dict = copy.deepcopy(AbstractAccuWorkflow.option_default_dict)
	option_default_dict.update({
						('input', 1, ): ['', 'i', 1, 'if it is a folder, take all .bam/.sam/.fastq files recursively. If it is a file, every line should be a path to the input file.', ],\
						('bamFname2MonkeyIDMapFname', 0, ): ['', '', 1, 'a tsv version of WUSTL xls file detailing what monkey is in which bam file.', ],\
						('minNoOfReads', 1, int): [8000000, '', 1, 'minimum number of reads in each split fastq file. This is the max NoOfReads.\
								 The upper limit in each split file is 2*minNoOfReads.', ],\
						("sequencer_name", 0, ): ["", '', 1, 'sequencing center of TCGA. parsed from TCGA bacode.'],\
						("sequence_type_name", 1, ): ["PairedEnd", '', 1, 'isq.sequence_type_id table column: SequenceType.short_name.'],\
						("sequence_format", 1, ): ["fastq", 'f', 1, 'fasta, fastq, etc.'],\
						("tissueSourceSiteFname", 0, ): ["", '', 1, 'TCGA tissue source site file'],\
						('inputType', 1, int): [1, 'y', 1, 'input type. 1: bam from TCGA 2: import HCC1187', ],\
						})
	#('jobFileDir', 0, ): ['', 'j', 1, 'folder to contain qsub scripts', ],\
	
	def __init__(self,  **keywords):
		"""
		2011-8-3
		"""
		AbstractAccuWorkflow.__init__(self, **keywords)
		self.addJobsDict = {1: self.addJobsToProcessTCGAData,
						2: self.addJobsToProcessHCC1187Data}
	
	def getBamBaseFname2MonkeyID_WUSTLDNAData(self, inputFname, ):
		"""
		2011-8-3
			from WUSTL
			the input looks like this:
			#	FlowCell	Lane	Index Sequence	Library	Common Name	Bam Path	MD5
			1	64J6AAAXX	1	VCAC-2007002-1-lib1	African	Green	Monkey	/gscmnt/sata755/production/csf_111215677/gerald_64J6AAAXX_1.bam	/gscmnt/sata755/production/csf_111215677/gerald_64J6AAAXX_1.bam.md5
			2	64J6AAAXX	2	VCAC-2007006-1-lib1	African	Green	Monkey	/gscmnt/sata751/production/csf_111215675/gerald_64J6AAAXX_2.bam	/gscmnt/sata751/production/csf_111215675/gerald_64J6AAAXX_2.bam.md5
		"""
		sys.stderr.write("Getting bamBaseFname2MonkeyID dictionary ...")
		bamBaseFname2MonkeyID = {}
		reader = csv.reader(open(inputFname), delimiter='\t')
		header = reader.next()
		col_name2index = utils.getColName2IndexFromHeader(header, skipEmptyColumn=True)
		monkeyIDIndex = col_name2index.get("Library")
		if monkeyIDIndex is None:	#2012.6.7
			monkeyIDIndex = col_name2index.get("library")
		bamFnameIndex = col_name2index.get("Bam Path")
		if bamFnameIndex is None:	#2012.2.9
			bamFnameIndex = col_name2index.get("BAM Path")
		if bamFnameIndex is None:	#2012.2.9
			bamFnameIndex = col_name2index.get("BAM")
		if bamFnameIndex is None:	#2012.6.7
			bamFnameIndex = col_name2index.get("bam pathway")
		#monkeyIDPattern = re.compile(r'\w+-(\w+)-\d+-\w+')	# i.e. VCAC-2007002-1-lib1
		monkeyIDPattern = re.compile(r'\w+-(\w+)-\w+-\w+')	# 2012.5.29 i.e. VCAC-VGA00006-AGM0075-lib1 ,
		# VCAC-VZC1014-AGM0055-lib1, VCAC-1996031-VRV0265-lib2a, VCAC-VKD7-361-VKD7-361-lib1 (VKD7 is to be taken),
		for row in reader:
			code = row[monkeyIDIndex]
			pa_search = monkeyIDPattern.search(code)
			if pa_search:
				code = pa_search.group(1)
			else:
				sys.stderr.write("Warning: could not parse monkey ID from %s. Ignore.\n"%(code))
				continue
			bamFname = row[bamFnameIndex]
			bamBaseFname = os.path.split(bamFname)[1]
			bamBaseFname2MonkeyID[bamBaseFname] = code
		sys.stderr.write("%s entries.\n"%(len(bamBaseFname2MonkeyID)))
		return bamBaseFname2MonkeyID
	
	def getAllBamFiles(self, inputDir, bamFnameLs=[]):
		"""
		2011-8-3
			recursively going through the directory to get all bam files
			
		"""
		
		for inputFname in os.listdir(inputDir):
			#get the absolute path
			inputFname = os.path.join(inputDir, inputFname)
			if os.path.isfile(inputFname):
				prefix, suffix = os.path.splitext(inputFname)
				if suffix=='.bam' or suffix=='.sam':
					bamFnameLs.append(inputFname)
			elif os.path.isdir(inputFname):
				self.getAllBamFiles(inputFname, bamFnameLs)
		
	def addIndividualSequence(self, db_main=None, code=None, name=None, isq_comment=None, tax_id=9606,
						tissue_id=None, site_id=None, study_id=None, sequence_batch_id=None,
						sequencer_name='HiSeq', sequence_type_name='PairedEnd', sequence_format='fastq',
						path_to_original_sequence=None, data_dir=None):
		"""
		add individual and then add individual_sequence
		"""
		individual = db_main.getIndividual(code=code, name=name, tax_id=tax_id, site_id=site_id, study_id=study_id)
		individual_sequence = db_main.getIndividualSequence(individual_id=individual.id, sequencer_name=sequencer_name, 
						sequence_type_name=sequence_type_name, sequence_format=sequence_format, 
						path_to_original_sequence=path_to_original_sequence, tissue_name=None, coverage=None,
						subFolder='individual_sequence', data_dir=data_dir, sequence_batch_id=sequence_batch_id, 
						tissue_id=tissue_id,
						filtered=0, version=None, is_contaminated=0, outdated_index=0, comment=isq_comment)
		
		return individual_sequence
	
	def registerCustomExecutables(self, workflow=None):
		"""
		20170503 No pegasus job clustering for these executables 
			because some of these jobs are certain to fail and will affect other good samples's downstream jobs. 
		2012.1.3
		"""
		self.addOneExecutableFromPathAndAssignProperClusterSize(
			path=os.path.join(self.pymodulePath, 'mapper/splitter/SplitReadFileWrapper.sh'), \
			name='SplitReadFileWrapper', clusterSizeMultipler=0)
		self.addOneExecutableFromPathAndAssignProperClusterSize(
			path=os.path.join(self.pymodulePath, 'mapper/computer/VerifyFileMD5Sum.py'), \
			name='VerifyFileMD5Sum', clusterSizeMultipler=0)
		self.addOneExecutableFromPathAndAssignProperClusterSize(
			path=os.path.join(self.thisModulePath, 'db/input/RegisterAndMoveSplitSequenceFiles.py'), \
			name='RegisterAndMoveSplitSequenceFiles', clusterSizeMultipler=0)
		self.addOneExecutableFromPathAndAssignProperClusterSize(
			path=os.path.join(self.thisModulePath, 'db/input/RegisterIndividualSequence2DB.py'), \
			name='RegisterIndividualSequence2DB', clusterSizeMultipler=0)
		self.addOneExecutableFromPathAndAssignProperClusterSize(path=self.javaPath, \
			name='SamToFastqJava', clusterSizeMultipler=0)
		
	def addConvertBamToFastqAndGzipJob(self, workflow=None, executable=None, \
							inputF=None, outputFnamePrefix=None, \
							parentJobLs=None, job_max_memory=8000, walltime = 800, extraDependentInputLs=None, \
							transferOutput=False, **keywords):
		"""
		2013.04.03 use addGenericJavaJob()
		2012.1.3
			walltime is in minutes (max time allowed on hoffman2 is 24 hours).
			The executable should be convertBamToFastqAndGzip.
			
		"""
		if workflow is None:
			workflow = self
		if extraDependentInputLs is None:
			extraDependentInputLs = []
		if executable is None:
			executable = self.SamToFastqJava
		output1 = File("%s_1.fastq.gz"%(outputFnamePrefix))
		output2 = File("%s_2.fastq.gz"%(outputFnamePrefix))
		output_unpaired = File("%s_unpaired.fastq.gz"%(outputFnamePrefix))
		extraOutputLs= [output1, output2, output_unpaired]
		extraArgumentList = ["F=", output1, "F2=", output2, "UNPAIRED_FASTQ=", output_unpaired, \
							"VALIDATION_STRINGENCY=LENIENT", "INCLUDE_NON_PF_READS=true"]
		
		job = self.addGenericJavaJob(executable=self.SamToFastqJava, jarFile=self.PicardJar, \
					inputFile=inputF, inputArgumentOption="INPUT=", \
					outputFile=None, outputArgumentOption="",\
					parentJobLs=parentJobLs, transferOutput=transferOutput, \
					frontArgumentList=["SamToFastq"], \
					extraArgumentList=extraArgumentList, extraOutputLs=extraOutputLs, \
					extraDependentInputLs=extraDependentInputLs, \
					job_max_memory=job_max_memory, walltime=walltime, **keywords)
		
		job.output1 = output1
		job.output2 = output2
		job.output_unpaired = output_unpaired
		return job
	
	def addSplitReadFileJob(self, workflow=None, executable=None, \
							inputF=None, outputFnamePrefix=None, outputFnamePrefixTail="",\
							minNoOfReads=5000000, logFile=None, parentJobLs=None, job_max_memory=4000, walltime = 800, \
							extraDependentInputLs=None, \
							transferOutput=False, **keywords):
		"""
		2013.04.03 use addGenericJob()
		2012.2.9
			argument outputFnamePrefixTail is now useless.
		2012.1.24
			executable is shell/SplitReadFileWrapper.sh
			
			which calls "wc -l" to count the number of reads beforehand to derive a proper minNoOfReads (to avoid files with too few reads).
			
			run SplitReadFile and generate the output directly into the db-affiliated folders
			
			a log file is generated and registered for transfer (so that pegasus won't skip it)
		
			walltime is in minutes (max time allowed on hoffman2 is 24 hours).
			
		"""
		if workflow is None:
			workflow = self
		if executable is None:
			executable = self.SplitReadFileWrapper
		if extraDependentInputLs is None:
			extraDependentInputLs = []
		frontArgumentList = [self.javaPath, repr(job_max_memory), workflow.SplitReadFileJar]
		extraArgumentList = [outputFnamePrefix, repr(minNoOfReads)]
		extraDependentInputLs.append(workflow.SplitReadFileJar)
		if logFile:
			extraArgumentList.append(logFile)
		
		job = self.addGenericJob(executable=executable, inputFile=inputF, inputArgumentOption="", \
					outputFile=None, outputArgumentOption="-o", \
					parentJob=None, parentJobLs=parentJobLs, extraDependentInputLs=extraDependentInputLs, \
					extraOutputLs=[logFile], transferOutput=transferOutput, \
					frontArgumentList=frontArgumentList, \
					extraArgumentList=extraArgumentList, \
					job_max_memory=job_max_memory, \
					walltime=walltime, **keywords)
		return job
	
	def addRegisterIndividualSequence2DBJob(self, workflow=None, executable=None, \
						inputFile=None, individual_id=None, \
						outputFile=None, \
						parentJobLs=None, job_max_memory=100, walltime = 60, commit=0, \
						extraDependentInputLs=None, extraArguments=None, \
						transferOutput=False, sshDBTunnel=1, **keywords):
		"""
		20170428
			walltime is in minutes (max time allowed on hoffman2 is 24 hours).
			
		"""
		if extraDependentInputLs is None:
			extraDependentInputLs = []
		if executable is None:
			executable = self.RegisterIndividualSequence2DB
		extraArgumentList =['--individual_id %s'%individual_id]
		for name, value in keywords.iteritems():
			if value is not None and value != "":
				extraArgumentList.append("--%s %s"%(name, value))
		job = self.addGenericFile2DBJob(executable=executable, \
					inputFile=inputFile, inputArgumentOption="-i", \
					outputFile=outputFile, outputArgumentOption="-o", inputFileList=None, \
					data_dir=None, commit=commit,\
					parentJobLs=parentJobLs, extraDependentInputLs=extraDependentInputLs, \
					extraOutputLs=None, transferOutput=transferOutput, \
					extraArguments=extraArguments, extraArgumentList=extraArgumentList, \
					job_max_memory=job_max_memory,  sshDBTunnel=sshDBTunnel, walltime=walltime,\
					key2ObjectForJob=None, objectWithDBArguments=self)
		return job
	
	def addRegisterAndMoveSplitFileJob(self, workflow=None, inputFile=None, \
							inputDir=None, logFile=None,\
							library=None, mate_id=None, \
							parentJobLs=None, job_max_memory=100, walltime = 60, \
							commit=0, sequence_format='fastq',\
							extraDependentInputLs=None, extraArguments=None, \
							transferOutput=False, sshDBTunnel=1, **keywords):
		"""
		2013.04.03 call addGenericFile2DBJob() instead
		2012.04.30
			add argument extraArguments, sshDBTunnel
		2012.01.24
			walltime is in minutes (max time allowed on hoffman2 is 24 hours).
			
		"""
		if extraDependentInputLs is None:
			extraDependentInputLs = []
		extraArgumentList =['--inputDir', inputDir, \
						'--library', library, '--sequence_format', sequence_format]
		if mate_id:
			extraArgumentList.append('--mate_id %s'%(mate_id))
		
		job = self.addGenericFile2DBJob(executable=self.RegisterAndMoveSplitSequenceFiles, \
					inputFile=inputFile, inputArgumentOption="-i", \
					outputFile=None, outputArgumentOption="-o", inputFileList=None, \
					data_dir=None, logFile=logFile, commit=commit,\
					parentJobLs=parentJobLs, extraDependentInputLs=extraDependentInputLs, \
					extraOutputLs=None, transferOutput=transferOutput, \
					extraArguments=extraArguments, extraArgumentList=extraArgumentList, \
					job_max_memory=job_max_memory,  sshDBTunnel=sshDBTunnel, walltime=walltime,\
					key2ObjectForJob=None, objectWithDBArguments=self)
		return job
	
	def getInputFnameLsFromInput(self, input=None, suffixSet=set(['.fastq']), fakeSuffix='.gz'):
		"""
		2012.4.30
			this function supercedes self.getAllBamFiles() and it's more flexible.
		"""
		inputFnameLs = []
		if os.path.isdir(input):
			self.getFilesWithSuffixFromFolderRecursive(inputFolder=input, suffixSet=suffixSet, fakeSuffix=fakeSuffix, \
												inputFnameLs=inputFnameLs)
		elif os.path.isfile(input):
			inf = open(input)
			for line in inf:
				inputFnameLs.append(line.strip())
			del inf
		else:
			sys.stderr.write("%s is neither a folder nor a file.\n"%(input))
			sys.exit(4)
		return inputFnameLs

	def getMonkeyID2FastqObjectLsForSouthAfricanDNAData(self, fastqFnameLs=None):
		"""
		2012.6.2
			similar to getMonkeyID2FastqObjectLsForMcGillData() (which is for McGill data)
			
			In pairs like this:
				VSAA2015_1.fastq.gz
				VSAA2015_2.fastq.gz
		"""
		sys.stderr.write("Passing monkeyID2FastqObjectLs from %s files ..."%(len(fastqFnameLs)))
		monkeyID2FastqObjectLs = {}
		import random
		filenameSignaturePattern = re.compile(r'(?P<monkeyID>[a-zA-Z0-9]+)_(?P<mateID>\d).fastq')
		counter = 0
		real_counter = 0
		libraryKey2UniqueLibrary = {}	#McGill's library ID , 7_Index-11, is not unique enough.
		for fastqFname in fastqFnameLs:
			counter += 1
			monkeyIDSearchResult = filenameSignaturePattern.search(fastqFname)
			if monkeyIDSearchResult:
				real_counter += 1
				monkeyID = monkeyIDSearchResult.group('monkeyID')
				mateID = monkeyIDSearchResult.group('mateID')
				filenameSignature = (monkeyID)
					
				#concoct a unique library ID
				libraryKey = (monkeyID)	#this combination insures two ends from the same library are grouped together
				if libraryKey not in libraryKey2UniqueLibrary:
					uniqueLibrary = '%s_%s'%(monkeyID, repr(random.random())[2:])
					libraryKey2UniqueLibrary[libraryKey] = uniqueLibrary
				
				uniqueLibrary = libraryKey2UniqueLibrary[libraryKey]
				fastqObject = PassingData(library=uniqueLibrary, monkeyID=monkeyID, mateID=mateID, absPath=fastqFname)
				if monkeyID not in monkeyID2FastqObjectLs:
					monkeyID2FastqObjectLs[monkeyID] = []
				monkeyID2FastqObjectLs[monkeyID].append(fastqObject)
			else:
				sys.stderr.write("Error: can't parse monkeyID, library, mateID out of %s.\n"%fastqFname)
				sys.exit(4)
		sys.stderr.write(" %s monkeys and %s files in the dictionary.\n"%(len(monkeyID2FastqObjectLs), real_counter))
		return monkeyID2FastqObjectLs

	def addJobsToProcessSouthAfricanDNAData(self, workflow=None, db_main=None, bamFname2MonkeyIDMapFname=None, input=None, data_dir=None, \
			minNoOfReads=None, commit=None,\
			sequencer_name=None, sequence_type_name=None, sequence_format=None):
		"""
		2012.6.2
			input fastq files could be gzipped or not. doesn't matter.
			data generated by Joe DeYoung's core, demultiplexed by ICNN
		"""
		fastqFnameLs = self.getInputFnameLsFromInput(input, suffixSet=set(['.fastq']), fakeSuffix='.gz')	#doesn't matter if fastq is not gzipped
		monkeyID2FastqObjectLs = self.getMonkeyID2FastqObjectLsForSouthAfricanDNAData(fastqFnameLs=fastqFnameLs)
		self.addJobsToSplitAndRegisterSequenceFiles(workflow=workflow, db_main=db_main, monkeyID2FastqObjectLs=monkeyID2FastqObjectLs, data_dir=data_dir, \
									minNoOfReads=minNoOfReads, commit=commit,\
									sequencer_name=sequencer_name, sequence_type_name=sequence_type_name, sequence_format=sequence_format)		
	
	def addJobsToProcessUNGCVervetData(self, workflow=None, db_main=None, bamFname2MonkeyIDMapFname=None, input=None, data_dir=None, \
			minNoOfReads=None, commit=None,\
			sequencer_name=None, sequence_type_name=None, sequence_format=None):
		"""
		2013.04.04
			UNGC = UCLA Neuroscience Genomics Core
			
			input fastq files could be gzipped or not. doesn't matter.
			data generated by Joe DeYoung's core, demultiplexed by ICNN
		"""
		sampleID2IndividualData = self.getSampleID2IndividualData_UNGC(bamFname2MonkeyIDMapFname)

		fastqFnameLs = self.getInputFnameLsFromInput(input, suffixSet=set(['.fastq']), fakeSuffix='.gz')	#doesn't matter if fastq is not gzipped
		monkeyID2FastqObjectLs = self.getMonkeyID2FastqObjectLsFromUNGCData(fastqFnameLs=fastqFnameLs, \
																	sampleID2IndividualData=sampleID2IndividualData)
		self.addJobsToSplitAndRegisterSequenceFiles(workflow=workflow, db_main=db_main, monkeyID2FastqObjectLs=monkeyID2FastqObjectLs, \
									data_dir=data_dir, \
									minNoOfReads=minNoOfReads, commit=commit,\
									sequencer_name=sequencer_name, sequence_type_name=sequence_type_name, sequence_format=sequence_format)
	
	def getSampleID2IndividualData_UNGC(self, inputFname=None):
		"""
		2013.04.04
			Format is like this from UNGC  = UCLA Neuroscience Genomics Core:
			FCID	Lane	sample ID	sample code	sample name	Index	Description	SampleProject
			D1HYNACXX	1	Ilmn Human control pool ( 4plex)	IP1			INDEX IS UNKNOWN prepared by Illumina (4 plex pool)	2013-029A
			D1HYNACXX	2	UNGC Human Sample 1	S1	AS001A	ATTACTCG	TruSeq DNA PCR Free beta kit	2013-029A
		"""
		sys.stderr.write("Getting  sampleID2IndividualData from %s ..."%(inputFname))
		sampleID2IndividualData = {}
		
		reader = MatrixFile(inputFname, openMode='r', delimiter=',')
		reader.constructColName2IndexFromHeader()
		sampleIDIndex = reader.getColIndexGivenColHeader("sample ID")
		sampleNameIndex = reader.getColIndexGivenColHeader("sample name")
		libraryIndexIndex = reader.getColIndexGivenColHeader("Index")
		
		for row in reader:
			sampleID = row[sampleIDIndex].replace(' ', '_')	#2013.04.04 stupid quirks
			sampleName = row[sampleNameIndex]
			libraryIndex = row[libraryIndexIndex]
			if sampleID not in sampleID2IndividualData:
				sampleID2IndividualData[sampleID] = PassingData(sampleName=sampleName, libraryIndexList=[])
			if sampleName!=sampleID2IndividualData[sampleID].sampleName:
				sys.stderr.write("Error: sampleID %s is associated with two different sample names (%s, %s).\n"%\
								(sampleID, sampleName, sampleID2IndividualData[sampleID].sampleName))
				raise
			sampleID2IndividualData[sampleID].libraryIndexList.append(libraryIndex)
		
		sys.stderr.write("%s entries.\n"%(len(sampleID2IndividualData)))
		return sampleID2IndividualData
	
	def getMonkeyID2FastqObjectLsFromUNGCData(self, fastqFnameLs=None, sampleID2IndividualData=None):
		"""
		2013.04.05 bugfix. this ensures two ends from same library have same library.
		2013.04.04
			UNGC  = UCLA Neuroscience Genomics Core
			In pairs like this:
				UNGC_Vervet_Sample_11_1.fastq.gz
				UNGC_Vervet_Sample_11_2.fastq.gz
		"""
		sys.stderr.write("Passing monkeyID2FastqObjectLs from %s files ..."%(len(fastqFnameLs)))
		monkeyID2FastqObjectLs = {}
		import re, random
		filenameSignaturePattern = re.compile(r'(?P<sampleID>[\w\d]+)_(?P<mateID>\d).fastq')
		counter = 0
		real_counter = 0
		libraryKey2UniqueLibrary = {}	#McGill's library ID , 7_Index-11, is not unique enough.
		for fastqFname in fastqFnameLs:
			counter += 1
			monkeyIDSearchResult = filenameSignaturePattern.search(fastqFname)
			if monkeyIDSearchResult:
				real_counter += 1
				sampleID = monkeyIDSearchResult.group('sampleID')
				mateID = monkeyIDSearchResult.group('mateID')
				if sampleID in sampleID2IndividualData:
					individualData = sampleID2IndividualData.get(sampleID)
					monkeyID = individualData.sampleName
					library = '_'.join(individualData.libraryIndexList)
					
					libraryKey = (monkeyID, library)
					if libraryKey not in libraryKey2UniqueLibrary:	#2013.04.05 bugfix. this ensures two ends from same library have same library.
						uniqueLibrary = '%s_%s_%s'%(monkeyID, library, repr(random.random())[2:])
						libraryKey2UniqueLibrary[libraryKey] = uniqueLibrary
					else:
						sys.stderr.write("libraryKey %s of %s is already in libraryKey2UniqueLibrary with unique library = %s.\n"%\
										(libraryKey, fastqFname, libraryKey2UniqueLibrary[libraryKey]))
					
					uniqueLibrary = libraryKey2UniqueLibrary[libraryKey]
					fastqObject = PassingData(library=uniqueLibrary, monkeyID=monkeyID, mateID=mateID, absPath=fastqFname)
					if monkeyID not in monkeyID2FastqObjectLs:
						monkeyID2FastqObjectLs[monkeyID] = []
					monkeyID2FastqObjectLs[monkeyID].append(fastqObject)
				else:
					sys.stderr.write("sampleID %s not in sampleID2IndividualData.\n"%(sampleID))
			else:
				sys.stderr.write("Error: can't parse sampleID, mateID out of %s.\n"%fastqFname)
				raise
		sys.stderr.write(" %s monkeys and %s files in the dictionary.\n"%(len(monkeyID2FastqObjectLs), real_counter))
		return monkeyID2FastqObjectLs
	
	
	def getFilenameSignature2MonkeyID_SouthAfricanRNAData(self, inputFname=None):
		"""
		2012.6.2 inputFname is tab-delimited, looks like this
		
			VSAC1012_R      ATCACG  sample2
			VSAF1009_R      ATCACG  sample1
			VSAB2009_RN     TAGCTT  sample2
			VSAB2011_R      TAGCTT  sample1
			VSAB3001_RN     GGCTAC  sample2
			VSAB5004_R      GGCTAC  sample1
			VSAA2015_RN     CTTGTA  sample2
		"""
		sys.stderr.write("Getting filenameSignature2MonkeyID dictionary from %s ..."%(inputFname))
		filenameSignature2MonkeyID = {}
		reader = csv.reader(open(inputFname), delimiter='\t')
		monkeyIDIndex = 0
		folderNameIndex = 1
		subSampleNameIndex = 2
		monkeyIDPattern = re.compile(r'(\w+)_R[\w]*')	# 2012.6.2 VSAA2015_RN or VSAB5004_R
		for row in reader:
			monkeyID = row[monkeyIDIndex]
			pa_search = monkeyIDPattern.search(monkeyID)
			if pa_search:
				monkeyID = pa_search.group(1)
			else:
				sys.stderr.write("Warning: could not parse monkey ID from %s. Ignore.\n"%(monkeyID))
				continue
			folderName = row[folderNameIndex]
			subSampleName = row[subSampleNameIndex]
			filenameSignature = (folderName, subSampleName)
			if filenameSignature in filenameSignature2MonkeyID:
				sys.stderr.write("Error: filenameSignature %s already in filenameSignature2MonkeyID.\n"%(filenameSignature))
				sys.exit(3)
			filenameSignature2MonkeyID[filenameSignature] = monkeyID
		sys.stderr.write("%s entries.\n"%(len(filenameSignature2MonkeyID)))
		return filenameSignature2MonkeyID
	
	def getMonkeyID2FastqObjectLsForSouthAfricanRNAData(self, fastqFnameLs=None, filenameSignature2MonkeyID=None):
		"""
		2012.6.2
			similar to getMonkeyID2FastqObjectLsForMcGillData() (which is for McGill data)
			
			In pairs like this:
			.../GCCAAT/tile_1101_sample1_end1.fastq
			.../GCCAAT/tile_1101_sample1_end2.fastq
		"""
		sys.stderr.write("Passing monkeyID2FastqObjectLs from %s files ..."%(len(fastqFnameLs)))
		monkeyID2FastqObjectLs = {}
		import re, random
		filenameSignaturePattern = re.compile(r'/(?P<folderName>[ACGT]{6})/(?P<library>[\w]+)_(?P<subSampleName>sample[12])_end(?P<mateID>\d).fastq')
		counter = 0
		real_counter = 0
		libraryKey2UniqueLibrary = {}	#McGill's library ID , 7_Index-11, is not unique enough.
		for fastqFname in fastqFnameLs:
			counter += 1
			monkeyIDSearchResult = filenameSignaturePattern.search(fastqFname)
			if monkeyIDSearchResult:
				real_counter += 1
				library = monkeyIDSearchResult.group('library')	#tile_1101
				folderName = monkeyIDSearchResult.group('folderName')
				subSampleName = monkeyIDSearchResult.group('subSampleName')
				mateID = monkeyIDSearchResult.group('mateID')
				filenameSignature = (folderName, subSampleName)
				if filenameSignature in filenameSignature2MonkeyID:
					monkeyID = filenameSignature2MonkeyID.get(filenameSignature)
					
					#concoct a unique library ID
					libraryKey = (folderName, library)	#this combination insures two ends from the same library are grouped together
					if libraryKey not in libraryKey2UniqueLibrary:
						uniqueLibrary = '%s_%s_%s'%(folderName, library, repr(random.random())[2:])
						libraryKey2UniqueLibrary[libraryKey] = uniqueLibrary
					
					uniqueLibrary = libraryKey2UniqueLibrary[libraryKey]
					fastqObject = PassingData(library=uniqueLibrary, monkeyID=monkeyID, mateID=mateID, absPath=fastqFname)
					if monkeyID not in monkeyID2FastqObjectLs:
						monkeyID2FastqObjectLs[monkeyID] = []
					monkeyID2FastqObjectLs[monkeyID].append(fastqObject)
				else:
					sys.stderr.write("%s not in filenameSignature2MonkeyID.\n"%(filenameSignature))
			else:
				sys.stderr.write("Error: can't parse monkeyID, library, mateID out of %s.\n"%fastqFname)
				sys.exit(4)
		sys.stderr.write(" %s monkeys and %s files in the dictionary.\n"%(len(monkeyID2FastqObjectLs), real_counter))
		return monkeyID2FastqObjectLs
	
	def addJobsToProcessSouthAfricanRNAData(self, workflow=None, db_main=None, bamFname2MonkeyIDMapFname=None, input=None, data_dir=None, \
			minNoOfReads=None, commit=None,\
			sequencer_name=None, sequence_type_name=None, sequence_format=None):
		"""
		2012.6.1
			input fastq files could be gzipped or not. doesn't matter.
			data generated by Joe DeYoung's core, demultiplexed by ICNN (Charles in particular)
		"""
		filenameSignature2MonkeyID = self.getFilenameSignature2MonkeyID_SouthAfricanRNAData(bamFname2MonkeyIDMapFname)
		
		fastqFnameLs = self.getInputFnameLsFromInput(input, suffixSet=set(['.fastq']), fakeSuffix='.gz')	#doesn't matter if fastq is not gzipped
		monkeyID2FastqObjectLs = self.getMonkeyID2FastqObjectLsForNamSouthAfricanRNAData(fastqFnameLs=fastqFnameLs, \
																	filenameSignature2MonkeyID=filenameSignature2MonkeyID)
		self.addJobsToSplitAndRegisterSequenceFiles(workflow=workflow, db_main=db_main, monkeyID2FastqObjectLs=monkeyID2FastqObjectLs, data_dir=data_dir, \
									minNoOfReads=minNoOfReads, commit=commit,\
									sequencer_name=sequencer_name, sequence_type_name=sequence_type_name, sequence_format=sequence_format)		
	
	def getMonkeyID2FastqObjectLsForMcGillData(self, fastqFnameLs=None,):
		"""
		2013.1.30 add fixed monkey ID prefixes (VWP, VGA, etc.) into monkeyIDPattern due to new not-all-number monkey IDs.
		
HI.0628.001.D701.VGA00010_R1.fastq.gz  HI.0628.004.D703.VWP00384_R1.fastq.gz  HI.0628.007.D703.VWP10020_R1.fastq.gz
HI.0628.001.D701.VGA00010_R2.fastq.gz  HI.0628.004.D703.VWP00384_R2.fastq.gz  HI.0628.007.D703.VWP10020_R2.fastq.gz
		2012.4.30
			each fastq file looks like 7_Index-11.2006013_R1.fastq.gz, 7_Index-11.2006013_R2.fastq.gz,
				7_Index-10.2005045_replacement_R1.fastq.gz, 7_Index-10.2005045_replacement_R2.fastq.gz
			
			The data from McGill is dated 2012.4.27.
			Each monkey is sequenced at 1X. There are 96 of them. Each library seems to contain 2 monkeys.
				But each monkey has only 1 library.
			
			8_Index_23.2008126_R1.fastq.gz
			8_Index_23.2008126_R2.fastq.gz
			8_Index_23.2009017_R1.fastq.gz
			8_Index_23.2009017_R2.fastq.gz

		"""
		sys.stderr.write("Passing monkeyID2FastqObjectLs from %s files ..."%(len(fastqFnameLs)))
		monkeyID2FastqObjectLs = {}
		import re, random
		#monkeyIDPattern = re.compile(r'(?P<library>[-\w]+)\.(?P<monkeyID>\d+)((_replacement)|(_pool)|())_R(?P<mateID>\d).fastq.gz')
		monkeyIDPattern = re.compile(r'(?P<library>[-\w]+)\.(?P<monkeyID>((VWP)|(VGA)|(VSA)|())\d+)((_replacement)|(_pool)|())_R(?P<mateID>\d).fastq.gz')
		counter = 0
		real_counter = 0
		libraryKey2UniqueLibrary = {}	#McGill's library ID , 7_Index-11, is not unique enough.
		for fastqFname in fastqFnameLs:
			counter += 1
			monkeyIDSearchResult = monkeyIDPattern.search(fastqFname)
			if monkeyIDSearchResult:
				real_counter += 1
				library = monkeyIDSearchResult.group('library')
				monkeyID = monkeyIDSearchResult.group('monkeyID')
				mateID = monkeyIDSearchResult.group('mateID')
				#concoct a unique library ID
				if library not in libraryKey2UniqueLibrary:
					libraryKey2UniqueLibrary[library] = '%s_%s'%(library, repr(random.random())[2:])
				uniqueLibrary = libraryKey2UniqueLibrary[library]
				fastqObject = PassingData(library=uniqueLibrary, monkeyID=monkeyID, mateID=mateID, absPath=fastqFname)
				if monkeyID not in monkeyID2FastqObjectLs:
					monkeyID2FastqObjectLs[monkeyID] = []
				monkeyID2FastqObjectLs[monkeyID].append(fastqObject)
			else:
				sys.stderr.write("Error: can't parse monkeyID, library, mateID out of %s.\n"%fastqFname)
				sys.exit(4)
		sys.stderr.write(" %s monkeys and %s files in the dictionary.\n"%(len(monkeyID2FastqObjectLs), real_counter))
		return monkeyID2FastqObjectLs
	
	def addJobsToSplitAndRegisterSequenceFiles(self, workflow=None, db_main=None, monkeyID2FastqObjectLs=None, data_dir=None, \
			minNoOfReads=None, commit=None,\
			sequencer_name=None, sequence_type_name=None, sequence_format=None):
		"""
		2012.6.2
			split out of addJobsToProcessMcGillData(), used also in addJobsToProcessDeYoungCoreData().
			
		"""
		if workflow is None:
			workflow = self
		sys.stderr.write("Adding split-read & register jobs ...")
		filenameKey2PegasusFile = {}
		for monkeyID, fastqObjectLs in monkeyID2FastqObjectLs.iteritems():
			individual_sequence = self.addIndividualSequence(db_main=db_main, code=monkeyID, \
								sequencer_name=sequencer_name, sequence_type_name=sequence_type_name, \
								sequence_format=sequence_format, data_dir=data_dir)
			
			sequenceOutputDir = os.path.join(data_dir, individual_sequence.path)
			sequenceOutputDirJob = self.addMkDirJob(outputDir=sequenceOutputDir)
			
			splitOutputDir = '%s'%(individual_sequence.id)
			#same directory containing split files from both mates is fine as RegisterAndMoveSplitSequenceFiles could pick up.
			splitOutputDirJob = self.addMkDirJob( outputDir=splitOutputDir)
			for fastqObject in fastqObjectLs:
				library = fastqObject.library
				mateID = fastqObject.mateID
				fastqPath = fastqObject.absPath
				filenameKey = (library, os.path.basename(fastqPath))
				if filenameKey in filenameKey2PegasusFile:
					fastqFile = filenameKey2PegasusFile.get(filenameKey)
					sys.stderr.write("Error: file %s has been registered with monkey %s. Can't happen.\n"%(fastqFile.monkeyID))
					sys.exit(3)
					import pdb
					pdb.set_trace()
					continue
				else:
					fastqFile = self.registerOneInputFile(workflow, fastqPath, folderName=library)
					fastqFile.monkeyID = monkeyID
					fastqFile.fastqObject= fastqObject
					filenameKey2PegasusFile[filenameKey] = fastqFile
			
				splitFastQFnamePrefix = os.path.join(splitOutputDir, '%s_%s_%s'%(individual_sequence.id, library, mateID))
				logFile = File('%s_%s_%s.split.log'%(individual_sequence.id, library, mateID))
				splitReadFileJob1 = self.addSplitReadFileJob(
								inputF=fastqFile, outputFnamePrefix=splitFastQFnamePrefix, \
								outputFnamePrefixTail="", minNoOfReads=minNoOfReads, \
								logFile=logFile, parentJobLs=[splitOutputDirJob], \
								job_max_memory=4000, walltime = 800, \
								extraDependentInputLs=[], transferOutput=True)
				
				logFile = File('%s_%s_%s.register.log'%(individual_sequence.id, library, mateID))
				# 2012.4.30 add '--mate_id_associated_with_bam' to RegisterAndMoveSplitSequenceFiles so that it will be used to distinguish IndividualSequenceFileRaw
				registerJob1 = self.addRegisterAndMoveSplitFileJob(
								inputDir=splitOutputDir, outputDir=sequenceOutputDir, relativeOutputDir=individual_sequence.path, logFile=logFile,\
								individual_sequence_id=individual_sequence.id, bamFile=fastqFile, library=library, mate_id=mateID, \
								parentJobLs=[splitReadFileJob1, sequenceOutputDirJob], job_max_memory=100, walltime = 60, \
								commit=commit, sequence_format=sequence_format, extraDependentInputLs=[], \
								extraArguments='--mate_id_associated_with_bam', \
								transferOutput=True, sshDBTunnel=self.needSSHDBTunnel)
			
		sys.stderr.write("%s jobs.\n"%(self.no_of_jobs))
	
	def addJobsToProcessMcGillData(self, workflow=None, db_main=None, bamFname2MonkeyIDMapFname=None, input=None, data_dir=None, \
			minNoOfReads=None, commit=None,\
			sequencer_name=None, sequence_type_name=None, sequence_format=None):
		"""
		2012.4.30
		"""
		fastqFnameLs = self.getInputFnameLsFromInput(input, suffixSet=set(['.fastq']), fakeSuffix='.gz')
		monkeyID2FastqObjectLs = self.getMonkeyID2FastqObjectLsForMcGillData(fastqFnameLs)
		
		self.addJobsToSplitAndRegisterSequenceFiles(workflow=workflow, db_main=db_main, \
									monkeyID2FastqObjectLs=monkeyID2FastqObjectLs, data_dir=data_dir, \
									minNoOfReads=minNoOfReads, commit=commit,\
									sequencer_name=sequencer_name, sequence_type_name=sequence_type_name, sequence_format=sequence_format)
		
	
	def addJobsToProcessWUSTLData(self, workflow, db_main=None, bamFname2MonkeyIDMapFname=None, input=None, data_dir=None, \
			minNoOfReads=None, commit=None,\
			sequencer_name=None, sequence_type_name=None, sequence_format=None):
		"""
		2012.4.30
		"""
		bamBaseFname2MonkeyID = self.getBamBaseFname2MonkeyID_WUSTLDNAData(bamFname2MonkeyIDMapFname)
		bamFnameLs = self.getInputFnameLsFromInput(input, suffixSet=set(['.bam', '.sam']), fakeSuffix='.gz')
		
		sys.stderr.write("%s total bam files.\n"%(len(bamFnameLs)))
		
		sam2fastqOutputDir = 'sam2fastq'
		sam2fastqOutputDirJob = self.addMkDirJob(outputDir=sam2fastqOutputDir)
		
		for bamFname in bamFnameLs:
			bamBaseFname = os.path.split(bamFname)[1]
			if bamBaseFname not in bamBaseFname2MonkeyID:
				sys.stderr.write("%s doesn't have code affiliated with.\n"%(bamFname))
				continue
			code = bamBaseFname2MonkeyID.get(bamBaseFname)
			individual_sequence = self.addIndividualSequence(db_main=db_main, code=code, \
								sequencer_name=sequencer_name, sequence_type_name=sequence_type_name, \
								sequence_format=sequence_format, data_dir=data_dir)
			#2012.2.10 stop passing path_to_original_sequence=bamFname to self.addMonkeySequence()
			
			"""
			#2012.2.10 temporary, during transition from old records to new ones.
			newISQPath = individual_sequence.constructRelativePathForIndividualSequence()
			newISQPath = '%s_split'%(newISQPath)
			if individual_sequence.path is None or individual_sequence.path !=newISQPath:
				individual_sequence.path = newISQPath
				session.add(individual_sequence)
				session.flush()
			"""
			
			sequenceOutputDir = os.path.join(data_dir, individual_sequence.path)
			sequenceOutputDirJob = self.addMkDirJob(outputDir=sequenceOutputDir)
			
			bamInputF = yh_pegasus.registerFile(workflow, bamFname)
			
			bamBaseFname = os.path.split(bamFname)[1]
			bamBaseFnamePrefix = os.path.splitext(bamBaseFname)[0]
			library = bamBaseFnamePrefix
			
			outputFnamePrefix = os.path.join(sam2fastqOutputDir, '%s_%s'%(individual_sequence.id, library))
			
			convertBamToFastqAndGzip_job = self.addConvertBamToFastqAndGzipJob(executable=workflow.SamToFastqJava, \
							inputF=bamInputF, outputFnamePrefix=outputFnamePrefix, \
							parentJobLs=[sam2fastqOutputDirJob], job_max_memory=6000, walltime = 800, \
							extraDependentInputLs=[], \
							transferOutput=False)
			
			splitOutputDir = '%s_%s'%(individual_sequence.id, library)
			#same directory containing split files from both mates is fine as RegisterAndMoveSplitSequenceFiles could pick up.
			splitOutputDirJob = self.addMkDirJob( outputDir=splitOutputDir)
			
			
			mate_id = 1
			splitFastQFnamePrefix = os.path.join(splitOutputDir, '%s_%s_%s'%(individual_sequence.id, library, mate_id))
			logFile = File('%s_%s_%s.split.log'%(individual_sequence.id, library, mate_id))
			splitReadFileJob1 = self.addSplitReadFileJob(workflow, executable=workflow.SplitReadFileWrapper, \
							inputF=convertBamToFastqAndGzip_job.output1, outputFnamePrefix=splitFastQFnamePrefix, \
							outputFnamePrefixTail="", minNoOfReads=minNoOfReads, \
							logFile=logFile, parentJobLs=[convertBamToFastqAndGzip_job, splitOutputDirJob], \
							job_max_memory=6000, walltime = 800, \
							extraDependentInputLs=[], transferOutput=True)
			
			logFile = File('%s_%s_%s.register.log'%(individual_sequence.id, library, mate_id))
			registerJob1 = self.addRegisterAndMoveSplitFileJob(workflow, executable=workflow.RegisterAndMoveSplitSequenceFiles, \
							inputDir=splitOutputDir, outputDir=sequenceOutputDir, relativeOutputDir=individual_sequence.path, logFile=logFile,\
							individual_sequence_id=individual_sequence.id, bamFile=bamInputF, library=library, mate_id=mate_id, \
							parentJobLs=[splitReadFileJob1, sequenceOutputDirJob], job_max_memory=100, walltime = 60, \
							commit=commit, sequence_format=sequence_format, extraDependentInputLs=[], \
							transferOutput=True, sshDBTunnel=self.needSSHDBTunnel)
			#handle the 2nd end
			mate_id = 2
			splitFastQFnamePrefix = os.path.join(splitOutputDir, '%s_%s_%s'%(individual_sequence.id, library, mate_id))
			logFile = File('%s_%s_%s.split.log'%(individual_sequence.id, library, mate_id))
			splitReadFileJob2 = self.addSplitReadFileJob(workflow, executable=workflow.SplitReadFileWrapper, \
							inputF=convertBamToFastqAndGzip_job.output2, outputFnamePrefix=splitFastQFnamePrefix, \
							outputFnamePrefixTail="", minNoOfReads=minNoOfReads, \
							logFile=logFile, parentJobLs=[convertBamToFastqAndGzip_job, splitOutputDirJob], \
							job_max_memory=6000, walltime = 800, \
							extraDependentInputLs=[], transferOutput=True)
			
			logFile = File('%s_%s_%s.register.log'%(individual_sequence.id, library, mate_id))
			registerJob1 = self.addRegisterAndMoveSplitFileJob(workflow, executable=workflow.RegisterAndMoveSplitSequenceFiles, \
							inputDir=splitOutputDir, outputDir=sequenceOutputDir, relativeOutputDir=individual_sequence.path, logFile=logFile,\
							individual_sequence_id=individual_sequence.id, bamFile=bamInputF, library=library, mate_id=mate_id, \
							parentJobLs=[splitReadFileJob2, sequenceOutputDirJob], job_max_memory=100, walltime = 60, \
							commit=commit, sequence_format=sequence_format, extraDependentInputLs=[], \
							transferOutput=True, sshDBTunnel=self.needSSHDBTunnel)
			
			"""
			
			jobFname = os.path.join(self.jobFileDir, 'job%s.bam2fastq.sh'%(code))
			self.writeQsubJob(jobFname, bamFname, os.path.join(self.data_dir, individual_sequence.path), self.vervet_path)
			commandline = 'qsub %s'%(jobFname)
			if self.commit:	#qsub only when db transaction will be committed.
				return_data = runLocalCommand(commandline, report_stderr=True, report_stdout=True)
			"""
		sys.stderr.write("%s jobs.\n"%(self.no_of_jobs))

	def parseTCGATissueSourceSite(self, db_main=None, inputFname=None):
		"""
		20170412
		"""
		sys.stderr.write("Parsing TCGA tissue source sites from %s... \n"%inputFname)
		reader = MatrixFile(inputFname=inputFname, delimiter="\t")
		header = reader.next()
		tssCode2dbEntry = {}
		for row in reader:
			tss_code, source_site_name, study_name, bcr_code = row[:4]
			if tss_code[0]=='0':	#remove the front 0
				tss_code = tss_code[1:]
			study = db_main.getStudy(short_name=study_name)
			site = db_main.getSite(short_name=tss_code, description=source_site_name, region=bcr_code, study_id=study.id)
			tssCode2dbEntry[tss_code] = site
		sys.stderr.write("%s sites.\n"% len(tssCode2dbEntry))
		
		return tssCode2dbEntry

	def getTCGASamplesFromInputDir(self, db_main=None, inputDir=None, tax_id=None):
		"""
		20170412
		COAD_TCGA/6795b83b-1ca7-4060-a078-6649fa9004bb
			35f008d5-1a15-463e-88e6-955583a33aa6_analysis.xml
			TCGA-AA-3553-10A-01D-1167-02_IlluminaHiSeq-DNASeq_whole.bam
			...
		"""
		sys.stderr.write("Getting TCGA samples from %s... \n"%inputDir)
		counter = 0
		noOfSamplesIntoDB = 0
		noOfSamplesWithoutCenter = 0
		
		tcga_sample_ls = []
		import xml.etree.ElementTree as ET
		tcga_barcode_re = re.compile(r'TCGA-(?P<tss>[0-9A-Z]{2})-(?P<participant>[0-9A-Z]{4})-(?P<tissue>[0-9]{2})[A-Z]-(?P<portion>[0-9]{2}[A-Z])-(?P<plate>[0-9A-Z]{4})(-(?P<center>[0-9]{2})|)')
		for filename in os.listdir(inputDir):
			inputPath = os.path.join(inputDir, filename)
			counter += 1
			if os.path.isdir(inputPath):
				noOfTargets = 0
				tcga_sample = PassingData(uuid=filename)
				for subfilename in os.listdir(inputPath):
					fname_prefix, fname_suffix = os.path.splitext(subfilename)
					if fname_suffix=='.bam':
						bamPath = os.path.join(inputPath, subfilename)
						#tcga_barcode = fname_prefix.split('_')[0], difficult to parse
						# f316e3dff60e757701731d0c0cd94c3d.bam
						# G92908.TCGA-A6-3809-01A-01D-A46W-08.2.bam
						# TCGA-AA-3672-01A-01D-0957-02_IlluminaGAII-DNASeq_whole.bam
						# TCGA-AO-A0JF-01A-11D-A060_130719_SN1120_0270_AC2CVRACXX_s_1_rg.sorted.bam
						tcga_sample.bamPath = bamPath
						tcga_sample.baiPath = "%s.bai"%bamPath
						noOfTargets+=1
					elif fname_suffix=='.xml' and fname_prefix[-9:]=="_analysis":
						#parse to get md5sum
						tree = ET.parse(os.path.join(inputPath, subfilename))
						root = tree.getroot()
						description = root.findall("ANALYSIS_SET/ANALYSIS/DESCRIPTION")[0].text
						#get participant code, tss_id (site.short_name), sample_id (tissue_id), center_id (sequencer.id)
						result = tcga_barcode_re.search(description)
						if result is None:
							sys.stderr.write("\n  ERROR: not enough (!=7) fields in tcga_barcode of %s/%s: %s.\n"%
											(inputPath, subfilename, description))
							sys.exit(-2)
						tcga_barcode = result.group()	##
						split_row = tcga_barcode.split('-')
						
						tss_code = result.group('tss')
						site = db_main.getSite(short_name=tss_code)
						
						participant_code = '-'.join(split_row[:3])
						tissue_id = int(result.group('tissue'))
						
						if result.group('center'):	#some don't have it
							sequencer_id = int(result.group('center'))
						else:
							sequencer_id = None
							noOfSamplesWithoutCenter += 1
						db_entry = db_main.getIndividual(code=participant_code, name=None, sex=None, age=None, \
														site=site, tax_id=tax_id, study=site.study)
						tcga_sample.db_entry = db_entry
						tcga_sample.tissue_id = tissue_id
						tcga_sample.sequencer_id = sequencer_id
						tcga_sample.tcga_barcode = tcga_barcode
						
						#get md5sum
						fileElemLs = root.findall("ANALYSIS_SET/ANALYSIS/DATA_BLOCK/FILES/FILE")
						if len(fileElemLs)!=1:
							sys.stderr.write("\n  ERROR: number of files with MD5SUM in %s/%s is not one. it is %s.\n"%
											(inputPath, subfilename, len(fileElemLs)))
							sys.exit(-2)
						filename = fileElemLs[0].get('filename')
						md5sum = fileElemLs[0].get('checksum')
						if filename is None or md5sum is None:
							import pdb
							pdb.set_trace()
						if filename[-4:]==".bam":
							tcga_sample.md5sum = md5sum
						else:
							sys.stderr.write("\n  ERROR: md5sum is for a non-bam file, %s.\n"%(filename))
							sys.exit(-3)
						noOfTargets+=1
				if noOfTargets==2:
					tcga_sample_ls.append(tcga_sample)
					noOfSamplesIntoDB += 1
		sys.stderr.write("%s samples in this folder %s. %s samples into DB. %s samples without center.\n"%(
			counter, inputDir, len(tcga_sample_ls), noOfSamplesWithoutCenter ))
		return tcga_sample_ls
		

	def addJobsToProcessTCGAData(self, workflow=None, db_main=None, input=None, tissueSourceSiteFname=None, \
					tax_id=9606, data_dir=None, \
					minNoOfReads=None, commit=None,\
					sequencer_name=None, sequence_type_name=None, sequence_format=None, **keywords):
		"""
		2017.4.7
			sequence_type_name is not used in this function. it's determined from tcga barcode.
			input:
				folder of bams
				tissueSourceSite.tsv file
				
		"""
		#parse the tissueSourceSite.tsv file to get all tissue source sites info, save them into table site, and study_id from db
		tssCode2dbEntry = self.parseTCGATissueSourceSite(db_main=db_main, inputFname=tissueSourceSiteFname)
		
		#find all bam in a folder,
		#parse barcode to get participant code, tss_id (site.short_name), sample_id (tissue_id), center_id (sequencer.id),
		#folder name is uuid
		#uuid is for sequence only, add as isq.comment
		#parse analysis.xml to get md5sum 
		#add all individual, tissue, site, study, sequencer, seq_center into db
		tcga_sample_obj_ls = self.getTCGASamplesFromInputDir(db_main=db_main, inputDir=input, tax_id=tax_id)
		sys.stderr.write("%s total TCGA samples.\n"%(len(tcga_sample_obj_ls)))
		
		
		sam2fastqOutputDir = 'sam2fastq'
		sam2fastqOutputDirJob = self.addMkDirJob(outputDir=sam2fastqOutputDir)
		
		for tcga_sample_obj in tcga_sample_obj_ls:
			bamInputF = yh_pegasus.registerFile(workflow, tcga_sample_obj.bamPath)
			bamFileSize = utils.getFileOrFolderSize(tcga_sample_obj.bamPath)
			baiInputF = yh_pegasus.registerFile(workflow, tcga_sample_obj.baiPath)
			bamBaseFname = os.path.split(tcga_sample_obj.bamPath)[1]
			bamBaseFnamePrefix = os.path.splitext(bamBaseFname)[0]
			library = tcga_sample_obj.tcga_barcode
			
			#add VerifyFileMD5Sum, cpu=8 to reduce its IO load
			verifyMD5SumJob = self.addGenericJob(executable=self.VerifyFileMD5Sum, \
					inputFile=bamInputF, \
					inputArgumentOption="-i", \
					outputFile=None, outputArgumentOption="-o", \
					inputFileList=None, argumentForEachFileInInputFileList=None, \
					parentJob=None, parentJobLs=None, extraDependentInputLs=None, \
					extraOutputLs=[bamInputF], transferOutput=False, \
					frontArgumentList=None, \
					extraArguments="--correct_md5sum %s"%tcga_sample_obj.md5sum, extraArgumentList=None, \
					job_max_memory=None, sshDBTunnel=None, \
					key2ObjectForJob=None, objectWithDBArguments=None, no_of_cpus=8, walltime=30)
			
			#no_of_cpus=4 to reduce its IO
			#assume a 10GB file needing a 30GB memory with a cap of 150G
			memoryNeeded = min(max(60000, int(bamFileSize/10000000000.0*30000)), 150000)
			outputFnamePrefix = os.path.join(sam2fastqOutputDir, '%s'%(library))
			convertBamToFastqAndGzip_job = self.addConvertBamToFastqAndGzipJob(executable=workflow.SamToFastqJava, \
					inputF=bamInputF, outputFnamePrefix=outputFnamePrefix, \
					parentJobLs=[sam2fastqOutputDirJob, verifyMD5SumJob], job_max_memory=memoryNeeded, no_of_cpus=4, walltime = 800, \
					extraDependentInputLs=[baiInputF], \
					transferOutput=False)
			
			#job to check if each file is empty or not. If empty, exit non-0. If not, add db isq entry, output isq-id, exit 0.
			
			#uuid is for sequence only, add as isq.comment
			unpaired_seq_db_idFile = File("%s_unpaired_seq_db_id.txt"%outputFnamePrefix)
			registerUnpairedSequence2DBJob =self.addRegisterIndividualSequence2DBJob(
						inputFile=convertBamToFastqAndGzip_job.output_unpaired, \
						outputFile=unpaired_seq_db_idFile, \
						individual_id=tcga_sample_obj.db_entry.id, sequencer_id=tcga_sample_obj.sequencer_id,\
						sequence_type_name="SingleRead", sequence_format=sequence_format, \
						copy_original_file=0, tissue_name=None, tissue_id=tcga_sample_obj.tissue_id, \
						coverage=None, quality_score_format="Standard", filtered=0,\
						parent_individual_sequence_id=None,\
						read_count=None, no_of_chromosomes=None, \
						sequence_batch_id=None, version=None, \
						is_contaminated=0, outdated_index=0, comment=tcga_sample_obj.uuid,\
						data_dir=data_dir,\
						original_sequence_filepath=tcga_sample_obj.bamPath, \
						original_sequence_library=library, \
						original_sequence_mate_id=None, \
						original_sequence_md5sum=tcga_sample_obj.md5sum, \
						parentJobLs=[convertBamToFastqAndGzip_job], \
						job_max_memory=100, walltime = 60, commit=self.commit, \
						extraDependentInputLs=[bamInputF], extraArguments=None, \
						transferOutput=False, sshDBTunnel=self.needSSHDBTunnel)
			
			splitSROutputDir = 'split_%s_singleRead'%(library)
			splitSROutputDirJob = self.addMkDirJob(outputDir=splitSROutputDir)
			
			splitFastQFnamePrefix = os.path.join(splitSROutputDir, '%s_singleRead'%(library))
			logFile = File('%s_split.log'%(splitFastQFnamePrefix))
			splitSRFileJob = self.addSplitReadFileJob(
							inputF=convertBamToFastqAndGzip_job.output_unpaired, \
							outputFnamePrefix=splitFastQFnamePrefix, \
							outputFnamePrefixTail="", minNoOfReads=minNoOfReads, \
							logFile=logFile, \
							parentJobLs=[registerUnpairedSequence2DBJob, splitSROutputDirJob], \
							job_max_memory=4000, walltime = 800, no_of_cpus=4, \
							extraDependentInputLs=[], transferOutput=True)
			
			#mate_id is set to null for singleRead sequences
			logFile = File('%s_register.log'%(splitFastQFnamePrefix))
			registerSRJob = self.addRegisterAndMoveSplitFileJob(
							inputFile=registerUnpairedSequence2DBJob.output,\
							inputDir=splitSROutputDir, logFile=logFile,\
							library=library, mate_id=None, \
							parentJobLs=[splitSRFileJob], job_max_memory=100, walltime = 60, \
							commit=commit, sequence_format=sequence_format, extraDependentInputLs=[], \
							transferOutput=True, sshDBTunnel=self.needSSHDBTunnel)
			
			#job to check if PE sequence file is empty or not. but only check one mate file.
			pairedEnd_seq_db_idFile = File("%s_pairedEnd_seq_db_id.txt"%(outputFnamePrefix))
			registerPairedEndSequence2DBJob =self.addRegisterIndividualSequence2DBJob(
						inputFile=convertBamToFastqAndGzip_job.output1, \
						outputFile=pairedEnd_seq_db_idFile, \
						individual_id=tcga_sample_obj.db_entry.id, sequencer_id=tcga_sample_obj.sequencer_id,\
						sequence_type_name="PairedEnd", sequence_format=sequence_format, \
						copy_original_file=0, tissue_name=None, tissue_id=tcga_sample_obj.tissue_id, \
						coverage=None, quality_score_format="Standard", filtered=0,\
						parent_individual_sequence_id=None,\
						read_count=None, no_of_chromosomes=None, \
						sequence_batch_id=None, version=None, \
						is_contaminated=0, outdated_index=0, comment=tcga_sample_obj.uuid,\
						data_dir=data_dir,\
						original_sequence_filepath=tcga_sample_obj.bamPath, \
						original_sequence_library=library, \
						original_sequence_mate_id=None, \
						original_sequence_md5sum=tcga_sample_obj.md5sum, \
						parentJobLs=[convertBamToFastqAndGzip_job], \
						job_max_memory=100, walltime = 60, commit=self.commit, \
						extraDependentInputLs=[bamInputF], extraArguments=None, \
						transferOutput=False, sshDBTunnel=self.needSSHDBTunnel)
			
			splitOutputDir = 'split_%s_PE'%(library)
			#same directory containing split files from both mates is fine as RegisterAndMoveSplitSequenceFiles could pick up.
			splitOutputDirJob = self.addMkDirJob(outputDir=splitOutputDir)
			
			mate_id = 1
			splitFastQFnamePrefix = os.path.join(splitOutputDir, '%s_%s'%(library, mate_id))
			logFile = File('%s_split.log'%(splitFastQFnamePrefix))
			splitReadFileJob1 = self.addSplitReadFileJob(
							inputF=convertBamToFastqAndGzip_job.output1, outputFnamePrefix=splitFastQFnamePrefix, \
							outputFnamePrefixTail="", minNoOfReads=minNoOfReads, \
							logFile=logFile, \
							parentJobLs=[registerPairedEndSequence2DBJob, splitOutputDirJob], \
							job_max_memory=4000, walltime = 800, no_of_cpus=4, \
							extraDependentInputLs=[], transferOutput=True)
			
			logFile = File('%s_register.log'%(splitFastQFnamePrefix))
			registerJob1 = self.addRegisterAndMoveSplitFileJob(
							inputFile=registerPairedEndSequence2DBJob.output,\
							inputDir=splitOutputDir, logFile=logFile,\
							library=library, mate_id=mate_id, \
							parentJobLs=[splitReadFileJob1], job_max_memory=100, walltime = 60, \
							commit=commit, sequence_format=sequence_format, extraDependentInputLs=[], \
							transferOutput=True, sshDBTunnel=self.needSSHDBTunnel)
			#handle the 2nd end
			mate_id = 2
			splitFastQFnamePrefix = os.path.join(splitOutputDir, '%s_%s'%(library, mate_id))
			logFile = File('%s_split.log'%(splitFastQFnamePrefix))
			splitReadFileJob2 = self.addSplitReadFileJob(
							inputF=convertBamToFastqAndGzip_job.output2, outputFnamePrefix=splitFastQFnamePrefix, \
							outputFnamePrefixTail="", minNoOfReads=minNoOfReads, \
							logFile=logFile, \
							parentJobLs=[registerPairedEndSequence2DBJob, splitOutputDirJob], \
							job_max_memory=4000, walltime = 800, no_of_cpus=4, \
							extraDependentInputLs=[], transferOutput=True)
			
			logFile = File('%s_register.log'%(splitFastQFnamePrefix))
			registerJob2 = self.addRegisterAndMoveSplitFileJob(
							inputFile=registerPairedEndSequence2DBJob.output,\
							inputDir=splitOutputDir, logFile=logFile,\
							library=library, mate_id=mate_id, \
							parentJobLs=[splitReadFileJob2], job_max_memory=100, walltime = 60, \
							commit=commit, sequence_format=sequence_format, extraDependentInputLs=[], \
							transferOutput=True, sshDBTunnel=self.needSSHDBTunnel)
			
		sys.stderr.write("%s jobs.\n"%(self.no_of_jobs))

	def getHCC1187SamplesFromInputDir(self, db_main=None, inputDir=None, tax_id=None):
		"""
		20170607
		inputDir is /y/home/luozhihui/Downloads/
		
			114G /y/home/luozhihui/Downloads/HCC1187BL_S1.bam
			8.7M /y/home/luozhihui/Downloads/HCC1187BL_S1.bam.bai
			237G /y/home/luozhihui/Downloads/HCC1187C_S1.bam
			8.7M /y/home/luozhihui/Downloads/HCC1187C_S1.bam.bai
		"""
		sys.stderr.write("Getting HCC1187 samples from %s... \n"%inputDir)
		counter = 0
		noOfSamplesIntoDB = 0
		
		tcga_sample_ls = []
		for filename in os.listdir(inputDir):
			inputPath = os.path.join(inputDir, filename)
			counter += 1
			if os.path.isfile(inputPath):
				fname_prefix, fname_suffix = os.path.splitext(filename)
				if fname_suffix=='.bam':
					tcga_sample = PassingData()
					code = fname_prefix.split("_")[0]
					bamPath = inputPath
					tcga_sample.bamPath = bamPath
					tcga_sample.baiPath = "%s.bai"%bamPath
					
					site_id = 10000
					
					if code=="HCC1187C":
						tissue_id = 50	#10 (Blood normal) or 50 (cancer cell line)
					elif code=="HCC1187BL":
						tissue_id = 10
					else:
						sys.stderr.write("\n ERROR: unexpected code %s in getting tissue_id.\n"%
											(code))
						import pdb
						pdb.set_trace()
						
					sequencer_id = 10000
					db_entry = db_main.getIndividual(code=code, name=None, sex=None, age=None, \
											site_id=site_id, tax_id=tax_id, study_id=10000)
					tcga_sample.db_entry = db_entry
					tcga_sample.tissue_id = tissue_id
					tcga_sample.sequencer_id = sequencer_id
					
					tcga_sample_ls.append(tcga_sample)
					noOfSamplesIntoDB += 1
		sys.stderr.write("%s samples in this folder %s. %s samples into DB.\n"%(
			counter, inputDir, len(tcga_sample_ls) ))
		return tcga_sample_ls

	def addJobsToProcessHCC1187Data(self, workflow=None, db_main=None, input=None, \
					tax_id=9606, data_dir=None, \
					minNoOfReads=None, commit=None,\
					sequencer_name=None, sequence_type_name=None, sequence_format=None, **keywords):
		"""
		20170607
			input:
				folder of bams
				
		"""
		
		#find all bam in a folder,
		#parse barcode to get participant code, tss_id (site.short_name), sample_id (tissue_id), center_id (sequencer.id),
		#folder name is uuid
		#uuid is for sequence only, add as isq.comment
		#parse analysis.xml to get md5sum 
		#add all individual, tissue, site, study, sequencer, seq_center into db
		tcga_sample_obj_ls = self.getHCC1187SamplesFromInputDir(db_main=db_main, inputDir=input, tax_id=tax_id)
		sys.stderr.write("%s total TCGA samples.\n"%(len(tcga_sample_obj_ls)))
		
		
		sam2fastqOutputDir = 'sam2fastq'
		sam2fastqOutputDirJob = self.addMkDirJob(outputDir=sam2fastqOutputDir)
		
		for tcga_sample_obj in tcga_sample_obj_ls:
			bamInputF = yh_pegasus.registerFile(workflow, tcga_sample_obj.bamPath)
			bamFileSize = utils.getFileOrFolderSize(tcga_sample_obj.bamPath)
			baiInputF = yh_pegasus.registerFile(workflow, tcga_sample_obj.baiPath)
			bamBaseFname = os.path.split(tcga_sample_obj.bamPath)[1]
			bamBaseFnamePrefix = os.path.splitext(bamBaseFname)[0]
			library = tcga_sample_obj.db_entry.code
			
			#no_of_cpus=4 to reduce its IO
			#assume a 10GB file needing a 30GB memory with a cap of 150G
			memoryNeeded = min(max(60000, int(bamFileSize/10000000000.0*30000)), 150000)
			outputFnamePrefix = os.path.join(sam2fastqOutputDir, '%s'%(library))
			convertBamToFastqAndGzip_job = self.addConvertBamToFastqAndGzipJob(executable=workflow.SamToFastqJava, \
					inputF=bamInputF, outputFnamePrefix=outputFnamePrefix, \
					parentJobLs=[sam2fastqOutputDirJob], job_max_memory=memoryNeeded, no_of_cpus=4, walltime = 800, \
					extraDependentInputLs=[baiInputF], \
					transferOutput=False)
			
			#job to check if each file is empty or not. If empty, exit non-0. If not, add db isq entry, output isq-id, exit 0.
			
			#uuid is for sequence only, add as isq.comment
			unpaired_seq_db_idFile = File("%s_unpaired_seq_db_id.txt"%outputFnamePrefix)
			registerUnpairedSequence2DBJob =self.addRegisterIndividualSequence2DBJob(
						inputFile=convertBamToFastqAndGzip_job.output_unpaired, \
						outputFile=unpaired_seq_db_idFile, \
						individual_id=tcga_sample_obj.db_entry.id, sequencer_id=tcga_sample_obj.sequencer_id,\
						sequence_type_name="SingleRead", sequence_format=sequence_format, \
						copy_original_file=0, tissue_name=None, tissue_id=tcga_sample_obj.tissue_id, \
						coverage=None, quality_score_format="Standard", filtered=0,\
						parent_individual_sequence_id=None,\
						read_count=None, no_of_chromosomes=None, \
						sequence_batch_id=None, version=None, \
						is_contaminated=0, outdated_index=0, comment=None,\
						data_dir=data_dir,\
						original_sequence_filepath=tcga_sample_obj.bamPath, \
						original_sequence_library=library, \
						original_sequence_mate_id=None, \
						original_sequence_md5sum=None, \
						parentJobLs=[convertBamToFastqAndGzip_job], \
						job_max_memory=100, walltime = 60, commit=self.commit, \
						extraDependentInputLs=[bamInputF], extraArguments=None, \
						transferOutput=False, sshDBTunnel=self.needSSHDBTunnel)
			
			splitSROutputDir = 'split_%s_singleRead'%(library)
			splitSROutputDirJob = self.addMkDirJob(outputDir=splitSROutputDir)
			
			splitFastQFnamePrefix = os.path.join(splitSROutputDir, '%s_singleRead'%(library))
			logFile = File('%s_split.log'%(splitFastQFnamePrefix))
			splitSRFileJob = self.addSplitReadFileJob(
							inputF=convertBamToFastqAndGzip_job.output_unpaired, \
							outputFnamePrefix=splitFastQFnamePrefix, \
							outputFnamePrefixTail="", minNoOfReads=minNoOfReads, \
							logFile=logFile, \
							parentJobLs=[registerUnpairedSequence2DBJob, splitSROutputDirJob], \
							job_max_memory=4000, walltime = 800, no_of_cpus=4, \
							extraDependentInputLs=[], transferOutput=True)
			
			#mate_id is set to null for singleRead sequences
			logFile = File('%s_register.log'%(splitFastQFnamePrefix))
			registerSRJob = self.addRegisterAndMoveSplitFileJob(
							inputFile=registerUnpairedSequence2DBJob.output,\
							inputDir=splitSROutputDir, logFile=logFile,\
							library=library, mate_id=None, \
							parentJobLs=[splitSRFileJob], job_max_memory=100, walltime = 60, \
							commit=commit, sequence_format=sequence_format, extraDependentInputLs=[], \
							transferOutput=True, sshDBTunnel=self.needSSHDBTunnel)
			
			#job to check if PE sequence file is empty or not. but only check one mate file.
			pairedEnd_seq_db_idFile = File("%s_pairedEnd_seq_db_id.txt"%(outputFnamePrefix))
			registerPairedEndSequence2DBJob =self.addRegisterIndividualSequence2DBJob(
						inputFile=convertBamToFastqAndGzip_job.output1, \
						outputFile=pairedEnd_seq_db_idFile, \
						individual_id=tcga_sample_obj.db_entry.id, sequencer_id=tcga_sample_obj.sequencer_id,\
						sequence_type_name="PairedEnd", sequence_format=sequence_format, \
						copy_original_file=0, tissue_name=None, tissue_id=tcga_sample_obj.tissue_id, \
						coverage=None, quality_score_format="Standard", filtered=0,\
						parent_individual_sequence_id=None,\
						read_count=None, no_of_chromosomes=None, \
						sequence_batch_id=None, version=None, \
						is_contaminated=0, outdated_index=0, comment=None,\
						data_dir=data_dir,\
						original_sequence_filepath=tcga_sample_obj.bamPath, \
						original_sequence_library=library, \
						original_sequence_mate_id=None, \
						original_sequence_md5sum=None, \
						parentJobLs=[convertBamToFastqAndGzip_job], \
						job_max_memory=100, walltime = 60, commit=self.commit, \
						extraDependentInputLs=[bamInputF], extraArguments=None, \
						transferOutput=False, sshDBTunnel=self.needSSHDBTunnel)
			
			splitOutputDir = 'split_%s_PE'%(library)
			#same directory containing split files from both mates is fine as RegisterAndMoveSplitSequenceFiles could pick up.
			splitOutputDirJob = self.addMkDirJob(outputDir=splitOutputDir)
			
			mate_id = 1
			splitFastQFnamePrefix = os.path.join(splitOutputDir, '%s_%s'%(library, mate_id))
			logFile = File('%s_split.log'%(splitFastQFnamePrefix))
			splitReadFileJob1 = self.addSplitReadFileJob(
							inputF=convertBamToFastqAndGzip_job.output1, outputFnamePrefix=splitFastQFnamePrefix, \
							outputFnamePrefixTail="", minNoOfReads=minNoOfReads, \
							logFile=logFile, \
							parentJobLs=[registerPairedEndSequence2DBJob, splitOutputDirJob], \
							job_max_memory=4000, walltime = 800, no_of_cpus=4, \
							extraDependentInputLs=[], transferOutput=True)
			
			logFile = File('%s_register.log'%(splitFastQFnamePrefix))
			registerJob1 = self.addRegisterAndMoveSplitFileJob(
							inputFile=registerPairedEndSequence2DBJob.output,\
							inputDir=splitOutputDir, logFile=logFile,\
							library=library, mate_id=mate_id, \
							parentJobLs=[splitReadFileJob1], job_max_memory=100, walltime = 60, \
							commit=commit, sequence_format=sequence_format, extraDependentInputLs=[], \
							transferOutput=True, sshDBTunnel=self.needSSHDBTunnel)
			#handle the 2nd end
			mate_id = 2
			splitFastQFnamePrefix = os.path.join(splitOutputDir, '%s_%s'%(library, mate_id))
			logFile = File('%s_split.log'%(splitFastQFnamePrefix))
			splitReadFileJob2 = self.addSplitReadFileJob(
							inputF=convertBamToFastqAndGzip_job.output2, outputFnamePrefix=splitFastQFnamePrefix, \
							outputFnamePrefixTail="", minNoOfReads=minNoOfReads, \
							logFile=logFile, \
							parentJobLs=[registerPairedEndSequence2DBJob, splitOutputDirJob], \
							job_max_memory=4000, walltime = 800, no_of_cpus=4, \
							extraDependentInputLs=[], transferOutput=True)
			
			logFile = File('%s_register.log'%(splitFastQFnamePrefix))
			registerJob2 = self.addRegisterAndMoveSplitFileJob(
							inputFile=registerPairedEndSequence2DBJob.output,\
							inputDir=splitOutputDir, logFile=logFile,\
							library=library, mate_id=mate_id, \
							parentJobLs=[splitReadFileJob2], job_max_memory=100, walltime = 60, \
							commit=commit, sequence_format=sequence_format, extraDependentInputLs=[], \
							transferOutput=True, sshDBTunnel=self.needSSHDBTunnel)
			
		sys.stderr.write("%s jobs.\n"%(self.no_of_jobs))

	def run(self):
		"""
		2011-8-3
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
		
		self.addJobsDict[self.inputType](workflow, db_main=db_main, bamFname2MonkeyIDMapFname=self.bamFname2MonkeyIDMapFname, \
					input=self.input, data_dir=self.data_dir, minNoOfReads=self.minNoOfReads, commit=self.commit,\
					sequencer_name=self.sequencer_name, sequence_type_name=self.sequence_type_name, \
					sequence_format=self.sequence_format, tissueSourceSiteFname=self.tissueSourceSiteFname)
		# Write the DAX to stdout
		outf = open(self.outputFname, 'w')
		self.writeXML(outf)
		if self.commit:
			session.commit()
		else:
			session.rollback()
		outf.close()
		
if __name__ == '__main__':
	main_class = ImportIndividualSequence2DB
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()
