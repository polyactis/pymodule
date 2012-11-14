#!/usr/bin/env python
"""
Examples:
	# 2011-9-29
	%s -a 525 -f 9 -I 8GenomeVsTop156Contigs_GATK_all_bases/call/ -N 0.0 -c 1
		-o 8GenomeVsTop156Contigs_GATK_all_bases_maxNA0_minMAF0_het2NA.xml -j condorpool -l condorpool  -u yh -z uclaOffice
	
	%s  -a 525 -f 9 -I 8GenomeVsTop156Contigs_GATK_all_bases/call/ -N 0.8 -M0
		-c 1 -o 8GenomeVsTop156Contigs_GATK_all_bases_maxNA0.8_minMAF0_het2NA.xml -j condorpool -l condorpool  -u yh -z uclaOffice
	
	#2012.5.11 on hoffman condor, 5 jobs per cluster (-C5), always need db connection on hcondor (-H)
	# add -U 0 -Z 3000 if u want to change the partitioning configuration (how many sites in one blast job)
	%s  -I ~/NetworkData/vervet/db/genotype_file/method_42/
		--oldRefFastaFname ~/NetworkData/vervet/db/individual_sequence/524_superContigsMinSize2000.fasta
		--newRefFastaFname ~/NetworkData/vervet/db/individual_sequence/3231_6483ContigsVervetRef1.0.3.fasta
		--formatdbPath ~/bin/blast/bin/formatdb --blastallPath ~/bin/blast/bin/blastall
		--maxNoOfMismatches 2 -H -C 5
		-j hcondor -l hcondor -D ~/NetworkData/vervet/db/ -t ~/NetworkData/vervet/db/
		-u yh -z localhost -o workflow/polymorphism/FindNewRefCoordinates_Method42_vs_3231_maxContigID2.xml
		#-U 0 -Z 3000

Description:
	2012-10-03
	#. preReduce: split the fasta into difference chromosomes (necessary?) 
    #. mapEachVCF:
    		split VCF into several small ones, N-SNP each
      #. mapEachInterval
         #. extract flanking sequences from the input VCF ():
         	ExtractFlankingSequenceForVCFLoci.py -i vcf -a reference fasta -o output.fasta -f flankingLength (=24)
         #. blast them
         #. run FindSNPPositionOnNewRefFromFlankingBlastOutput.py
             #. where hit length match query length, and no of mismatches <=2 => good => infer new coordinates
         #. output a mapping file between old SNP and new SNP coordinates.
         		#. reduce this thing by combining everything
    #. reduce()
         #. merge all output into one big one
         
"""
import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

from pymodule import ProcessOptions, getListOutOfStr, PassingData, yh_pegasus
from Pegasus.DAX3 import *
from pymodule.pegasus.AbstractVCFWorkflow import AbstractVCFWorkflow
from pymodule.pegasus.BlastWorkflow import BlastWorkflow

class FindNewRefCoordinatesGivenVCFFolderWorkflow(AbstractVCFWorkflow, BlastWorkflow):
	__doc__ = __doc__
	option_default_dict = AbstractVCFWorkflow.option_default_dict.copy()
	option_default_dict.update({
						('oldRefFastaFname', 1, ): ['', '', 1, 'path to the old reference sequence file (on which input VCF is based)', ],\
						("formatdbPath", 1, ): ["/usr/bin/formatdb", 'f', 1, 'path to formatdb, index fasta database file'],\
						("blastallPath", 1, ): ["/usr/bin/blastall", 's', 1, 'path to blastall'],\
						('newRefFastaFname', 1, ): ['', '', 1, 'path to the new reference sequence file (blast db)', ],\
						
						('minNoOfIdentities', 0, int): [None, '', 1, 'minimum number of identities between a query and target', ],\
						('maxNoOfMismatches', 0, int): [None, '', 1, 'minimum number of mismatches between a query and target', ],\
						('minIdentityPercentage', 0, float): [None, '', 1, 'minimum percentage of identities between a query and target', ],\
						('flankingLength', 1, int): [24, '', 1, 'number of flanking bases on either side of the locus.\n\
	length of flanking = 2*flankingLength+locusLength', ],\
						('minAlignmentSpan', 1, int): [10, '', 1, 'minimum length of alignment in blast', ],\
						})
	
	#2012.9.25 no overlap and make the interval a lot smaller, (for VCF file)
	option_default_dict[('intervalOverlapSize', 1, int)][0] = 0
	option_default_dict[('intervalSize', 1, int)][0] = 3000
	
	def __init__(self,  **keywords):
		"""
		2011-7-11
		"""
		AbstractVCFWorkflow.__init__(self, **keywords)
		self.oldRefFastaFname = os.path.abspath(self.oldRefFastaFname)
		self.newRefFastaFname = os.path.abspath(self.newRefFastaFname)
	
	def preReduce(self, workflow=None, outputDirPrefix="", passingData=None, transferOutput=True, **keywords):
		"""
		2012.9.17
		"""
		AbstractVCFWorkflow.preReduce(self, workflow=workflow, outputDirPrefix=outputDirPrefix,\
								passingData=passingData, transferOutput=transferOutput, **keywords)
		
		
		returnData = PassingData(no_of_jobs = 0)
		returnData.jobDataLs = []
		self.oldRefFastaFileList = self.registerBlastNucleotideDatabaseFile(self.oldRefFastaFname, \
							input_site_handler=self.input_site_handler)
		self.newRefFastaFileList = self.registerBlastNucleotideDatabaseFile(self.newRefFastaFname, \
							input_site_handler=self.input_site_handler)
		
		if len(self.newRefFastaFileList)<4:	#some nt-database index file is missing
			sys.stderr.write("Adding blast-db-making job ...")
			makeBlastDBJob = self.addMakeBlastDBJob(executable=self.formatdb,\
										inputFile=self.newRefFastaFileList[0], transferOutput=True)
			#add the index files to the ntDatabaseFileList
			self.newRefFastaFileList = [self.newRefFastaFileList[0]] + makeBlastDBJob.outputList
			sys.stderr.write(".\n")
		else:
			makeBlastDBJob = None
		passingData.makeBlastDBJob = makeBlastDBJob
		return returnData
		
	def registerCustomExecutables(self, workflow=None):
		"""
		2011-11-28
		"""
		if workflow==None:
			workflow=self
		
		AbstractVCFWorkflow.registerCustomExecutables(self, workflow)
		BlastWorkflow.registerCustomExecutables(self, workflow)
		
		namespace = workflow.namespace
		version = workflow.version
		operatingSystem = workflow.operatingSystem
		architecture = workflow.architecture
		clusters_size = workflow.clusters_size
		site_handler = workflow.site_handler
		vervetSrcPath = workflow.vervetSrcPath
		
		executableClusterSizeMultiplierList = []	#2012.8.7 each cell is a tuple of (executable, clusterSizeMultipler (0 if u do not need clustering)		
		
		ExtractFlankingSequenceForVCFLoci = Executable(namespace=namespace, name="ExtractFlankingSequenceForVCFLoci", \
											version=version, \
											os=operatingSystem, arch=architecture, installed=True)
		ExtractFlankingSequenceForVCFLoci.addPFN(PFN("file://" + os.path.join(workflow.pymodulePath, \
												"pegasus/mapper/ExtractFlankingSequenceForVCFLoci.py"), \
										site_handler))
		executableClusterSizeMultiplierList.append((ExtractFlankingSequenceForVCFLoci, 2))
		
		FindSNPPositionOnNewRefFromFlankingBlastOutput = Executable(namespace=namespace, name="FindSNPPositionOnNewRefFromFlankingBlastOutput", \
											version=version, \
											os=operatingSystem, arch=architecture, installed=True)
		FindSNPPositionOnNewRefFromFlankingBlastOutput.addPFN(PFN("file://" + os.path.join(workflow.pymodulePath, \
														"polymorphism/FindSNPPositionOnNewRefFromFlankingBlastOutput.py"), \
										site_handler))
		executableClusterSizeMultiplierList.append((FindSNPPositionOnNewRefFromFlankingBlastOutput, 2))
		
		self.addExecutableAndAssignProperClusterSize(executableClusterSizeMultiplierList, defaultClustersSize=self.clusters_size)
	
	def mapEachInterval(self, workflow=None, \
					VCFFile=None, passingData=None, transferOutput=False, **keywords):
		"""
		2012.10.3
			#. extract flanking sequences from the input VCF (ref sequence file => contig ref sequence)
			#. blast them
			#. run FindSNPPositionOnNewRefFromFlankingBlastOutput.py
				#. where hit length match query length, and no of mismatches <=2 => good => infer new coordinates
			#. output a mapping file between old SNP and new SNP coordinates.
				#. reduce this thing by combining everything
			#. make a new VCF file based on the input split VCF file
				(replace contig ID , position with the new one's, remove the header part regarding chromosomes or replace it)

		"""
		
		returnData = PassingData(no_of_jobs = 0)
		returnData.jobDataLs = []

		topOutputDirJob = passingData.topOutputDirJob
		mapDirJob = passingData.mapDirJob
		reduceOutputDirJob = passingData.reduceOutputDirJob
		
		intervalFnamePrefix = passingData.intervalFnamePrefix
		jobData = passingData.jobData
		
		
		splitVCFJob = passingData.mapEachVCFData.splitVCFJob
		chr = passingData.chromosome
		
		# a flanking sequence extraction job
		
		outputFnamePrefix = os.path.join(mapDirJob.output, '%s_flankSequence'%(intervalFnamePrefix)) 
		outputFile = File('%s.fasta'%(outputFnamePrefix))
		oldRefFastaFile = self.oldRefFastaFileList[0]
		extraArgumentList = ['--flankingLength %s'%(self.flankingLength), '--refFastaFname', oldRefFastaFile]
		extractFlankSeqJob = self.addAbstractMapperLikeJob(workflow=workflow, executable=self.ExtractFlankingSequenceForVCFLoci, \
					inputF=VCFFile, outputF=outputFile, \
					parentJobLs=[mapDirJob, splitVCFJob], transferOutput=transferOutput, job_max_memory=2000,\
					extraArguments=None, extraArgumentList=extraArgumentList, extraDependentInputLs=[oldRefFastaFile])
		
		# a blast job
		blastOutputFnamePrefix = '%s_blast'%(outputFnamePrefix)
		outputFile = File('%s.tsv'%(blastOutputFnamePrefix))
		newRefFastaFile = self.newRefFastaFileList[0]
		blastJob = self.addBlastWrapperJob(executable=self.BlastWrapper, inputFile=extractFlankSeqJob.output, \
										outputFile=outputFile, \
						outputFnamePrefix=blastOutputFnamePrefix, databaseFile=newRefFastaFile,\
						maxNoOfMismatches=self.maxNoOfMismatches, minNoOfIdentities=self.minNoOfIdentities, \
						minIdentityPercentage=self.minIdentityPercentage, blastallPath=self.blastallPath, \
						parentJobLs=[mapDirJob, extractFlankSeqJob, passingData.makeBlastDBJob], \
						extraDependentInputLs=self.newRefFastaFileList, \
						transferOutput=transferOutput, \
						extraArguments=None, job_max_memory=5000)
		
		# a FindSNPPositionOnNewRefFromFlankingBlastOutput job
		outputFile = File('%s_2NewRefCoordinateMap.tsv'%(outputFnamePrefix))
		extraArgumentList = ['--minAlignmentSpan %s'%(self.minAlignmentSpan)]
		if self.maxNoOfMismatches:
			extraArgumentList.append('--maxNoOfMismatches %s'%(self.maxNoOfMismatches))
		if self.minNoOfIdentities:
			extraArgumentList.append("--minNoOfIdentities %s"%(self.minNoOfIdentities))
		if self.minIdentityPercentage:
			extraArgumentList.append("--minIdentityPercentage %s"%(self.minIdentityPercentage))
			
		findNewRefCoordinatesJob = self.addAbstractMapperLikeJob(workflow=workflow, \
					executable=self.FindSNPPositionOnNewRefFromFlankingBlastOutput, \
					inputF=blastJob.output, outputF=outputFile, \
					parentJobLs=[blastJob], transferOutput=transferOutput, job_max_memory=500,\
					extraArguments=None, extraArgumentList=extraArgumentList, extraDependentInputLs=None)
		returnData.findNewRefCoordinatesJob = findNewRefCoordinatesJob
		return returnData
	
	def reduce(self, workflow=None, passingData=None, reduceEachChromosomeDataLs=None, transferOutput=True, **keywords):
		"""
		2012.10.3
			#. merge all output of input jobs (passingData.mapEachIntervalDataLsLs) into one big one
		
		"""
		returnData = PassingData(no_of_jobs = 0)
		returnData.jobDataLs = []
		reduceOutputDirJob = passingData.reduceOutputDirJob
		
		fnamePrefix = os.path.join(reduceOutputDirJob.output, 'SNPID2NewCoordinates.tsv')
		outputFile = File('%s.tsv'%(fnamePrefix))
		reduceJob = self.addStatMergeJob(workflow, \
									statMergeProgram=self.mergeSameHeaderTablesIntoOne, \
									outputF=outputFile, \
									parentJobLs=[reduceOutputDirJob],extraOutputLs=[], \
									extraDependentInputLs=[], transferOutput=transferOutput,)
		returnData.jobDataLs.append(PassingData(jobLs=[reduceJob], file=reduceJob.output, \
											fileList=[reduceJob.output]))
		
		for mapEachIntervalDataLs in passingData.mapEachIntervalDataLsLs:
			for mapEachIntervalData in mapEachIntervalDataLs:
				parentJob = mapEachIntervalData.findNewRefCoordinatesJob
				self.addInputToStatMergeJob(workflow, statMergeJob=reduceJob, \
						parentJobLs=[parentJob])
		
		return returnData

if __name__ == '__main__':
	main_class = FindNewRefCoordinatesGivenVCFFolderWorkflow
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()