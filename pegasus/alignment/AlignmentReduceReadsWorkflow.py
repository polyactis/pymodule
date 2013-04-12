#!/usr/bin/env python
"""
Examples:
	# 2012.9.21 run base quality recalibration on VRC alignments (-S 447), individual_sequence_id from 639-642 (--ind_seq_id_ls ...)
	# filtered sequences (-Q 1), alignment method 2 (-G 2)
	# --contigMaxRankBySize 1000 (top 1000 contigs)
	#  --intervalSize 10000000 (10 million bp for each interval) --intervalOverlapSize 30000 (30kb overlap between intervals),
	%s --inputDir ~/NetworkData/vervet/db/genotype_file/method_17/ --ind_seq_id_ls 639-642
		-S 447 -u yh -z localhost --sequence_filtered 1 --alignment_method_id 2
		-a 524 -o dags/BaseQualityRecalibration/BaseQualityRecalibration_VRC447_vsMethod17.xml
		-l hcondor -j hcondor -z localhost -u yh --contigMaxRankBySize 1000 
		--intervalSize 10000000 --intervalOverlapSize 30000
		--indelVCFFolder ...
		-e /u/home/eeskin/polyacti
		--local_data_dir /u/home/eeskin/polyacti/NetworkData/vervet/db/ --data_dir /u/home/eeskin/polyacti/NetworkData/vervet/db/
		--needSSHDBTunnel -J ~/bin/jdk/bin/java --mask_genotype_method_id 17
		 --commit --skipDoneAlignment

	# 2012.9.18
	%s  --inputDir ~/NetworkData/vervet/db/genotype_file/method_41 --ind_seq_id_ls 633,634,635,636,637,638 
		--ref_ind_seq_id 524
		-o dags/BaseQualityRecalibration/BaseQualityRecalibration_ISQ633_638_vsMethod41.xml -l hcondor
		-j hcondor -z localhost -u yh --intervalSize 10000000 --intervalOverlapSize 30000
		-e /u/home/eeskin/polyacti
		--indelVCFFolder ...
		--local_data_dir /u/home/eeskin/polyacti/NetworkData/vervet/db/ --data_dir /u/home/eeskin/polyacti/NetworkData/vervet/db/
		--clusters_size 5 --needSSHDBTunnel -J ~/bin/jdk/bin/java --mask_genotype_method_id 41
		 --commit --skipDoneAlignment
	
	# 2013.3.19 use sequence coverage to filter alignments
	%s  --inputDir ~/NetworkData/vervet/db/genotype_file/method_41
		--sequence_min_coverage 0 --sequence_max_coverage 2  --ind_seq_id_ls 632-3230
		--ref_ind_seq_id 3280 -o dags/BaseQualityRecalibration/BaseQualityRecalibration_ISQ632_3230_coverage0_2_vsMethod41.xml
		-l hcondor -j hcondor -z localhost -u yh --intervalSize 10000000 --intervalOverlapSize 30000
		-e /u/home/eeskin/polyacti --contigMaxRankBySize 250
		--local_data_dir /u/home/eeskin/polyacti/NetworkData/vervet/db/ --data_dir /u/home/eeskin/polyacti/NetworkData/vervet/db/
		--clusters_size 5 --needSSHDBTunnel -J ~/bin/jdk/bin/java --mask_genotype_method_id 41
		--indelVCFFolder ~/NetworkData/vervet/db/genotype_file/method_88 --commit --skipDoneAlignment
		# --ref_genome_version 2 #(optional, as by default, it gets the outdated_index=0 reference chromosomes from GenomeDB)
		# --ref_genome_outdated_index 1 #to get old reference. incompatible here as alignment is based on 3280, new ref.
		# --needFastaDictJob --needFastaIndexJob
	
Description:
	#2013.04.11 workflow that carries out ReduceReads for all alignments in map-reduce fashion.
		preferably on local_realigned=1 alignments. It is set to 1 by default.
"""
import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

from Pegasus.DAX3 import *
from pymodule import ProcessOptions, getListOutOfStr, PassingData, yh_pegasus, NextGenSeq, \
	figureOutDelimiter, getColName2IndexFromHeader, utils
#from pymodule.pegasus.AbstractVCFWorkflow import AbstractVCFWorkflow
from pymodule import VCFFile
from vervet.src import AbstractVervetAlignmentWorkflow

parentClass = AbstractVervetAlignmentWorkflow
class AlignmentReduceReadsWorkflow(parentClass):
	__doc__ = __doc__
	option_default_dict = parentClass.option_default_dict.copy()
	option_default_dict.update(parentClass.commonAlignmentWorkflowOptionDict.copy())
	option_default_dict.update(parentClass.partitionWorkflowOptionDict.copy())
	option_default_dict.update({
							})
	option_default_dict[('intervalSize', 1, int)][0] = 20000000
	option_default_dict[('local_realigned', 0, int)][0] = 1
	
	def __init__(self,  **keywords):
		"""
		"""
		parentClass.__init__(self, **keywords)
		self.chr2IndelVCFJobData = None	#2013.04.04 mark this variable. setup in setup()
		self.candidateCountCovariatesJob = None	#2013.04.09 this BQSR count-variates job encompasses one of the top big intervals.
			# replacing equivalent jobs for small intervals (not accurate if intervals are too small)
		#AlignmentToCallPipeline.__init__(self, **keywords)
		#self.inputDir = os.path.abspath(self.inputDir)
	
	def mapEachInterval(self, workflow=None, alignmentData=None, intervalData=None, chromosome=None, \
							VCFJobData=None, passingData=None, reduceBeforeEachAlignmentData=None,\
							mapEachChromosomeData=None, transferOutput=False, \
							**keywords):
		"""
		2013.03.31 use VCFJobData to decide whether to add BQSR jobs, called in ShortRead2AlignmentWorkflow.py
		2012.9.17
		"""
		if workflow is None:
			workflow = self
		returnData = PassingData(no_of_jobs = 0)
		returnData.jobDataLs = []
		
		topOutputDirJob = passingData.topOutputDirJob
		
		alignment = alignmentData.alignment
		bamF = alignmentData.bamF
		baiF = alignmentData.baiF
		bamFnamePrefix = passingData.bamFnamePrefix
		
		#SNPVCFFile = VCFJobData.file
		#if SNPVCFFile is None or VCFJobData is None:	#2013.04.09	BQSR requires a VCF input regardless of the chromosome
		#	VCFJobData = self.randomSNPVCFJobDataForBQSR
		
		#SNPVCFFile = VCFJobData.file
		#SNPVCFJobLs = VCFJobData.jobLs
		
		if intervalData.file:
			mpileupInterval = intervalData.interval
			bcftoolsInterval = intervalData.file
		else:
			mpileupInterval = intervalData.interval
			bcftoolsInterval = intervalData.interval
		intervalFnameSignature = intervalData.intervalFnameSignature
		overlapInterval = intervalData.overlapInterval
		overlapFilenameSignature = intervalData.overlapIntervalFnameSignature
		span = intervalData.span
		
		if chromosome is None:
			chromosome = getattr(passingData, 'chromosome', None)
		
		
		median_depth = getattr(alignment, 'median_depth', 4)
		readSpace = median_depth * span
		#base is 4X coverage in 20Mb region => 120 minutes
		reduceReadsJobWalltime = self.scaleJobWalltimeOrMemoryBasedOnInput(realInputVolume=readSpace, \
							baseInputVolume=4*20000000, baseJobPropertyValue=120, \
							minJobPropertyValue=60, maxJobPropertyValue=500).value
		#base is 4X, => 5000M
		reduceReadsJobMaxMemory = self.scaleJobWalltimeOrMemoryBasedOnInput(realInputVolume=median_depth, \
							baseInputVolume=4, baseJobPropertyValue=4000, \
							minJobPropertyValue=4000, maxJobPropertyValue=8000).value
							
		reduceReadsBamFile = File(os.path.join(topOutputDirJob.output, '%s_%s.reduceReads.bam'%\
											(bamFnamePrefix, overlapFilenameSignature)))
		#2013.04.09 GATK generates this file. it is not .bam.bai but just .bai. check if this is True 
		reduceReadsBaiFile = File('%s.bai'%(os.path.splitext(reduceReadsBamFile.name)[0]))
		reduceReadsJob = self.addGATKJob(executable=self.IndelRealignerJava, GenomeAnalysisTKJar=self.GenomeAnalysisTK2Jar, \
					GATKAnalysisType='ReduceReads',\
					inputFile=bamF, inputArgumentOption="-I", refFastaFList=passingData.refFastaFList, inputFileList=None,\
					argumentForEachFileInInputFileList=None,\
					interval=overlapInterval, outputFile=reduceReadsBamFile, \
					parentJobLs=alignmentData.jobLs, transferOutput=False, \
					job_max_memory=reduceReadsJobMaxMemory,\
					frontArgumentList=None, extraArguments=None, \
					extraArgumentList=None, \
					extraOutputLs=[reduceReadsBaiFile], \
					extraDependentInputLs=[baiF], no_of_cpus=None, \
					walltime=reduceReadsJobWalltime)
		indexBamJob = self.addBAMIndexJob(BuildBamIndexFilesJava=self.BuildBamIndexFilesJava, \
										BuildBamIndexJar=self.BuildBamIndexJar, \
					inputBamF=reduceReadsJob.output,\
					parentJobLs=[reduceReadsJob], \
					transferOutput=False, job_max_memory=3000, \
					walltime=max(120, int(reduceReadsJobWalltime/3)))
		passingData.AlignmentJobAndOutputLs.append(PassingData(parentJobLs=[reduceReadsJob, indexBamJob], \
															file=reduceReadsJob.output))
		return returnData
	
	def reduceAfterEachAlignment(self, workflow=None, passingData=None, transferOutput=False, data_dir=None, **keywords):
		"""
		"""
		returnData = PassingData(no_of_jobs = 0)
		returnData.jobDataLs = []
		if workflow is None:
			workflow = self
		AlignmentJobAndOutputLs = getattr(passingData, 'AlignmentJobAndOutputLs', [])
		bamFnamePrefix = passingData.bamFnamePrefix
		topOutputDirJob = passingData.topOutputDirJob
		individual_alignment = passingData.individual_alignment
		reduceOutputDirJob = passingData.reduceOutputDirJob
		
		if len(AlignmentJobAndOutputLs)>0:	#2012.3.29	merge alignment output only when there is something to merge!
			#2013.04.09 create a new child alignment local_realigned =1, etc.
			new_individual_alignment = self.db.copyParentIndividualAlignment(parent_individual_alignment_id=individual_alignment.id,\
										mask_genotype_method_id=self.mask_genotype_method_id,\
										data_dir=self.data_dir, local_realigned=individual_alignment.local_realigned,\
										reduce_reads=1)
			
			#2013.04.09 replace read_group with the new one to each alignment job
			NewAlignmentJobAndOutputLs = []
			for AlignmentJobAndOutput in AlignmentJobAndOutputLs:
				#2012.9.19 add a AddReadGroup job
				alignmentJob, indexAlignmentJob = AlignmentJobAndOutput.parentJobLs
				fname_prefix = os.path.splitext(alignmentJob.output.name)[0]
				outputRGBAM = File("%s.isq_RG.bam"%(fname_prefix))
				addRGJob = self.addReadGroupInsertionJob(workflow=workflow, individual_alignment=new_individual_alignment, \
									inputBamFile=alignmentJob.output, \
									outputBamFile=outputRGBAM,\
									addOrReplaceReadGroupsJava=self.addOrReplaceReadGroupsJava, \
									AddOrReplaceReadGroupsJar=self.AddOrReplaceReadGroupsJar,\
									parentJobLs=[alignmentJob, indexAlignmentJob], extraDependentInputLs=None, \
									extraArguments=None, job_max_memory = 2500, transferOutput=False)
				
				NewAlignmentJobAndOutputLs.append(PassingData(parentJobLs=[addRGJob], file=addRGJob.output))
			#
			mergedBamFile = File(os.path.join(reduceOutputDirJob.output, '%s.merged.bam'%(bamFnamePrefix)))
			alignmentMergeJob, bamIndexJob = self.addAlignmentMergeJob(workflow, AlignmentJobAndOutputLs=NewAlignmentJobAndOutputLs, \
					outputBamFile=mergedBamFile, \
					samtools=workflow.samtools, java=workflow.java, \
					MergeSamFilesJava=workflow.MergeSamFilesJava, MergeSamFilesJar=workflow.MergeSamFilesJar, \
					BuildBamIndexFilesJava=workflow.IndexMergedBamIndexJava, BuildBamIndexJar=workflow.BuildBamIndexJar, \
					mv=workflow.mv, parentJobLs=[reduceOutputDirJob], \
					transferOutput=False)
			#2012.9.19 add/copy the alignment file to db-affliated storage
			#add the metric file to AddAlignmentFile2DB.py as well (to be moved into db-affiliated storage)
			logFile = File(os.path.join(reduceOutputDirJob.output, '%s_2db.log'%(bamFnamePrefix)))
			alignment2DBJob = self.addAddAlignmentFile2DBJob(workflow=workflow, executable=self.AddAlignmentFile2DB, \
								inputFile=alignmentMergeJob.output, otherInputFileList=[],\
								individual_alignment_id=new_individual_alignment.id, \
								mask_genotype_method_id=self.mask_genotype_method_id,\
								logFile=logFile, data_dir=data_dir, \
								parentJobLs=[alignmentMergeJob, bamIndexJob], \
								extraDependentInputLs=[bamIndexJob.output], \
								extraArguments=None, transferOutput=transferOutput, \
								job_max_memory=2000, sshDBTunnel=self.needSSHDBTunnel, commit=True)
			self.no_of_jobs += 1
			returnData.jobDataLs.append(PassingData(jobLs=[alignment2DBJob], file=alignment2DBJob.logFile, \
											fileLs=[alignment2DBJob.logFile]))
		return returnData
	
	
	def isThisAlignmentComplete(self, individual_alignment=None, data_dir=None):
		"""
		2013.04.09 this is more complicated as it tests the local_realigned version of individual_alignment is complete or not.
			not individual_alignment itself
		"""
		new_individual_alignment = self.db.copyParentIndividualAlignment(parent_individual_alignment_id=individual_alignment.id,\
										mask_genotype_method_id=self.mask_genotype_method_id,\
										data_dir=self.data_dir, local_realigned=1)
		return self.db.isThisAlignmentComplete(individual_alignment=new_individual_alignment, data_dir=data_dir)
	
	def registerCustomExecutables(self, workflow=None):
		"""
		2011-11-28
		"""
		parentClass.registerCustomExecutables(self, workflow=workflow)
		
		if workflow is None:
			workflow = self
		
		executableClusterSizeMultiplierList = []	#2012.8.7 each cell is a tuple of (executable, clusterSizeMultipler (0 if u do not need clustering)
		self.addExecutableAndAssignProperClusterSize(executableClusterSizeMultiplierList, defaultClustersSize=self.clusters_size)
		self.setOrChangeExecutableClusterSize(executable=workflow.samtools, clusterSizeMultipler=1)
		
		#self.addOneExecutableFromPathAndAssignProperClusterSize(path=self.javaPath, name='IndelRealignerJava', clusterSizeMultipler=0.2)


if __name__ == '__main__':
	main_class = AlignmentReduceReadsWorkflow
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()
