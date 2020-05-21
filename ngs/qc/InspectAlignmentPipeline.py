#!/usr/bin/env python3
"""
Examples:
    #Change tmpDir (--tmpDir) for AddOrReplaceReadGroups,
    #  no job clustering (--cluster_size 1)
    %s -a 524 -j condorpool -l condorpool -u yh -z uclaOffice
        -o dags/InspectAln1_To_661_RefSeq524Alignments.xml
        --ind_aln_id_ls 1-661
        --tmpDir /Network/Data/vervet/vervetPipeline/tmp/ --cluster_size 1

    #Do perContig depth estimation (--needPerContigJob) and
    #  skip alignments with stats in db already (--skipAlignmentWithStats)
    # need ssh tunnel for db (--needSSHDBTunnel)
    # add --individual_sequence_file_raw_id_type 2 (library-specific alignments,
    # 	 different libraries of one individual_sequence)
    # add --individual_sequence_file_raw_id_type 3
    # 	 (both all-library-fused and library-specific alignments)
    # add "--country_id_ls 135,136,144,148,151" to limit individuals 
    # 	from US,Barbados,StKitts,Nevis,Gambia (AND with -S, )
    %s -a 524 -l hcondor -u yh -z localhost --contigMaxRankBySize 7559
        -o dags/InspectAln1_To_1251_RefSeq524Alignments.xml
        --ind_aln_id_ls 1-1251 --cluster_size 1
        --home_path /u/home/eeskin/polyacti/
        -J ~/bin/jdk/bin/java
        --needPerContigJob --skipAlignmentWithStats --needSSHDBTunnel
        #--individual_sequence_file_raw_id_type 2 
        # --country_id_ls 135,136,144,148,151 --tax_id_ls 60711 #sabaeus
        #--ind_seq_id_ls 632-3230 --site_id_ls 447 --sequence_filtered 1
        #--excludeContaminant	#VRC sequences
        #--sequence_filtered 1 --alignment_method_id  2

    #--min_segment_length
    mLength=100;
    %s
        --alignmentDepthIntervalMethodShortName 16CynosurusRef3488MinLength$mLength
        -a 3488 --ref_genome_tax_id 60711
        --ref_genome_sequence_type_id 1 --ref_genome_version 1
        --contigMaxRankBySize 3000
        --ind_seq_id_ls 632-5000 --excludeContaminant
        --alignment_method_id 6 --sequence_filtered 1
        --tax_id_ls 460675
        --country_id_ls 1,129,130,131,132,133,134,136,144,148,151,152
        --skipAlignmentWithStats
        --min_segment_length $mLength
        -l hcondor -z localhost -u yh
        -o dags/InspectAlignment_RefSeq3488MinLength$mLength\_AlnMethod6.xml
        --cluster_size 1 -J ~/bin/jdk/bin/java
        --needSSHDBTunnel

Description:
    A pegasus workflow that inspects no-of-reads-aligned, infers insert size and etc.
    Use samtools flagstat.
"""
import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0], sys.argv[0])

import copy
import getpass
import logging
from pegaflow.DAX3 import File, Link, Job
from palos import ProcessOptions, utils, PassingData
from palos.ngs.AbstractAlignmentWorkflow import AbstractAlignmentWorkflow

ParentClass = AbstractAlignmentWorkflow

class InspectAlignmentPipeline(ParentClass):
    __doc__ = __doc__
    #	("fractionToSample", 0, float): [0.001, '', 1, 
    # 'fraction of loci to walk through for DepthOfCoverage walker.'],\
    def __init__(self, 
        drivername='postgresql', hostname='localhost',
        dbname='', schema='public', port=None,
        db_user=None,
        db_passwd=None,
        data_dir=None, local_data_dir=None,

        ind_aln_id_ls=None,
        ind_seq_id_ls=None,
        alignment_method_id=None,
        excludeContaminant=False,
        sequence_filtered=None,

        min_segment_length=100,
        needPerContigJob=False,
        skipAlignmentWithStats=False,
        alignmentDepthIntervalMethodShortName=None,
        
        samtools_path="bin/samtools",
        picard_dir="script/picard/dist",
        gatk_path="bin/GenomeAnalysisTK1_6_9.jar",
        gatk2_path="bin/GenomeAnalysisTK.jar",
        picard_path="script/picard.broad/build/libs/picard.jar",

        contigMaxRankBySize=None,
        contigMinRankBySize=None,

        chromosome_type_id=None, 
        ref_genome_tax_id=9606,
        ref_genome_sequence_type_id=1,
        ref_genome_version=15,
        ref_genome_outdated_index=0,

        pegasusFolderName='input',
        site_handler='condor',
        input_site_handler='condor',
        cluster_size=30,
        output_path=None,
        tmpDir='/tmp/',
        max_walltime=4320,
        
        home_path=None,
        javaPath=None,
        pymodulePath="src/pymodule",

        jvmVirtualByPhysicalMemoryRatio=1.2,

        needSSHDBTunnel=False,
        commit=False,
        debug=False,
        report=False):
        """
        """
        ParentClass.__init__(self, 
            drivername=drivername, hostname=hostname,
            dbname=dbname, schema=schema, port=port,
            db_user=db_user, db_passwd=db_passwd,
            data_dir=data_dir, local_data_dir=local_data_dir,

            ind_aln_id_ls=ind_aln_id_ls,
            completedAlignment=1,
            ind_seq_id_ls=ind_seq_id_ls,
            alignment_method_id=alignment_method_id,
            excludeContaminant=excludeContaminant,
            sequence_filtered=sequence_filtered,

            samtools_path=samtools_path,
            picard_dir=picard_dir,
            gatk_path=gatk_path,
            gatk2_path=gatk2_path,
            picard_path=picard_path,

            contigMaxRankBySize=contigMaxRankBySize,
            contigMinRankBySize=contigMinRankBySize,

            chromosome_type_id=chromosome_type_id, 
            ref_genome_tax_id=ref_genome_tax_id,
            ref_genome_sequence_type_id=ref_genome_sequence_type_id,
            ref_genome_version=ref_genome_version,
            ref_genome_outdated_index=ref_genome_outdated_index,

            pegasusFolderName=pegasusFolderName,
            site_handler=site_handler,
            input_site_handler=input_site_handler,
            cluster_size=cluster_size,
            output_path=output_path,
            tmpDir=tmpDir,
            max_walltime=max_walltime, 
            
            home_path=home_path,
            javaPath=javaPath,
            pymodulePath=pymodulePath,
            
            jvmVirtualByPhysicalMemoryRatio=jvmVirtualByPhysicalMemoryRatio,
            needSSHDBTunnel=needSSHDBTunnel,
            commit=commit,
            debug=debug, report=report)
        
        self.min_segment_length = min_segment_length
        self.needPerContigJob = needPerContigJob
        self.skipAlignmentWithStats = skipAlignmentWithStats
        self.alignmentDepthIntervalMethodShortName = \
            alignmentDepthIntervalMethodShortName
        
        self.no_of_alns_with_depth_jobs = 0
        self.no_of_alns_with_flagstat_jobs = 0

        self.alignmentDepthJobDataList=[]
        self.needSplitChrIntervalData = False
    
    def addDepthOfCoverageJob(self,
        DOCWalkerJava=None, GenomeAnalysisTKJar=None,
        refFastaFList=None, bamF=None, baiF=None, DOCOutputFnamePrefix=None,
        fractionToSample=None, minMappingQuality=20, minBaseQuality=20,
        extraArguments="",
        parentJobLs=None, transferOutput=False,
        job_max_memory = 1000, walltime=None, **keywords):
        """
        2013.06.12
        Bugfix, instead of --minBaseQuality, it was --maxBaseQuality passed to GATK.
        Set minMappingQuality (was 30) to 20.
        2013.06.09
            .sample_statistics is new GATK DOC output file,
             (replacing the .sample_interval_summary file).
            ignore argument fractionToSample, not available.
        2013.05.17
            re-activate this because "samtools depth" seems to have 
            trouble working with local-realigned and BQSR-ed bam files
            use addGATKJob()
        2012.5.7
            no longer used, superceded by addSAMtoolsDepthJob()
        2012.4.17
            add --omitIntervalStatistics and --omitLocusTable to the walker
        2012.4.12
            add "--read_filter BadCigar" to GATK to avoid stopping because of malformed Cigar
                malformed: Read starting with deletion. Cigar: 1D65M299S
        2012.4.3
            add argument fractionToSample
        2011-11-25
        """
        sample_summary_file = File('%s.sample_summary'%(DOCOutputFnamePrefix))
        sample_statistics_file = File('%s.sample_statistics'%(DOCOutputFnamePrefix))
        extraOutputLs = [sample_summary_file, sample_statistics_file]
        extraArgumentList = ["-o", DOCOutputFnamePrefix, \
            "-pt sample", "--read_filter BadCigar", \
            "--omitDepthOutputAtEachBase", '--omitLocusTable',
            '--omitIntervalStatistics']
        if minMappingQuality is not None:
            extraArgumentList.append("--minMappingQuality %s"%(minMappingQuality))
        if minBaseQuality is not None:
            extraArgumentList.append("--minBaseQuality %s"%(minBaseQuality))

        #if fractionToSample and fractionToSample>0 and fractionToSample<=1:
        #	extraArgumentList.append("--fractionToSample %s"%(fractionToSample))
        extraDependentInputLs = [baiF]
        job = self.addGATKJob(executable=DOCWalkerJava, \
            GenomeAnalysisTKJar=GenomeAnalysisTKJar, \
            GATKAnalysisType="DepthOfCoverage",\
            inputFile=bamF, inputArgumentOption="-I", 
            refFastaFList=refFastaFList, inputFileList=None,\
            interval=None, outputFile=None, \
            parentJobLs=parentJobLs, transferOutput=transferOutput, \
            job_max_memory=job_max_memory,\
            frontArgumentList=None, extraArguments=extraArguments, 
            extraArgumentList=extraArgumentList, \
            extraOutputLs=extraOutputLs, \
            extraDependentInputLs=extraDependentInputLs,
            no_of_cpus=1, walltime=walltime, **keywords)
        job.sample_summary_file = sample_summary_file
        job.sample_statistics_file = sample_statistics_file
        return job

    def addSAMtoolsDepthJob(self, samtoolsDepth=None, samtools_path=None,\
        bamF=None, outputFile=None, baiF=None, \
        parentJobLs=None, extraOutputLs=None, job_max_memory = 500, \
        extraArguments=None, \
        transferOutput=False, minMappingQuality=None, minBaseQuality=None, \
        walltime=120, **keywords):
        """
        2013.08.27 default minMappingQuality and minBaseQuality are set to None
        2013.3.24 use addGenericJob()
        2012.5.7

        """
        extraArgumentList = []
        if minMappingQuality is not None:
            extraArgumentList.append("%s"%minMappingQuality)
        if minBaseQuality is not None:
            extraArgumentList.append("%s"%minBaseQuality)
        job= self.addGenericJob(executable=samtoolsDepth, \
            frontArgumentList=[samtools_path],\
            inputFile=bamF, inputArgumentOption=None,\
            outputFile=outputFile, outputArgumentOption=None,\
            parentJobLs=parentJobLs, extraDependentInputLs=[baiF], \
            extraOutputLs=extraOutputLs, extraArguments=extraArguments, \
            transferOutput=transferOutput, \
            extraArgumentList=extraArgumentList, \
            key2ObjectForJob=None, job_max_memory=job_max_memory, \
            sshDBTunnel=None, walltime=walltime, **keywords)
        return job


    def addReformatFlagstatOutputJob(self, executable=None, inputF=None,
        alignmentID=None, outputF=None, \
        parentJobLs=None, extraDependentInputLs=None,
        transferOutput=False,
        extraArguments=None, job_max_memory=2000, walltime=20, **keywords):
        """
        2013.3.24 use addGenericJob()
        2012.4.3
            input is output of "samtools flagstat"
        """
        job= self.addGenericJob(executable=executable, \
            inputFile=inputF, inputArgumentOption='-i',\
            outputFile=outputF, outputArgumentOption='-o',\
            parentJobLs=parentJobLs,
            extraDependentInputLs=extraDependentInputLs,
            extraOutputLs=None, extraArguments=extraArguments, \
            transferOutput=transferOutput, \
            extraArgumentList=['-a %s'%(alignmentID)], \
            key2ObjectForJob=None, job_max_memory=job_max_memory, \
            sshDBTunnel=None, walltime=walltime, **keywords)
        return job

    def preReduce(self, passingData=None, transferOutput=True, **keywords):
        """
        2012.9.17
        setup additional mkdir folder jobs, 
            before mapEachAlignment, mapEachChromosome, mapReduceOneAlignment
        """
        returnData = ParentClass.preReduce(self, passingData=passingData, \
                            transferOutput=transferOutput)
        reduceOutputDirJob = passingData.reduceOutputDirJob

        self.logOutputDirJob = self.addMkDirJob(outputDir="%sLog"%(
            passingData.outputDirPrefix))

        depthOfCoverageOutputF = File(os.path.join(reduceOutputDirJob.output,
            'DepthOfCoverage.tsv'))
        passingData.depthOfCoverageOutputMergeJob = self.addStatMergeJob(
            statMergeProgram=self.mergeSameHeaderTablesIntoOne, \
            outputF=depthOfCoverageOutputF, parentJobLs=[reduceOutputDirJob],
            transferOutput=True)

        if self.needPerContigJob:	#need for per-contig job
            depthOfCoveragePerChrOutputF = File(\
                os.path.join(reduceOutputDirJob.output, 'DepthOfCoveragePerChr.tsv'))
            passingData.depthOfCoveragePerChrOutputMergeJob = self.addStatMergeJob(\
                statMergeProgram=self.mergeSameHeaderTablesIntoOne, \
                outputF=depthOfCoveragePerChrOutputF, 
                parentJobLs=[reduceOutputDirJob], transferOutput=True)
        else:
            passingData.depthOfCoveragePerChrOutputMergeJob = None

        flagStatOutputF = File(os.path.join(reduceOutputDirJob.output, 'FlagStat.tsv'))
        passingData.flagStatOutputMergeJob = self.addStatMergeJob(
            statMergeProgram=self.mergeSameHeaderTablesIntoOne, \
            outputF=flagStatOutputF,
            parentJobLs=[reduceOutputDirJob], transferOutput=True)

        passingData.alignmentDataLs = self.addAddRG2BamJobsAsNeeded(
            alignmentDataLs=passingData.alignmentDataLs,
            site_handler=self.site_handler, \
            input_site_handler=self.input_site_handler, \
            AddOrReplaceReadGroupsJar=self.AddOrReplaceReadGroupsJar, \
            BuildBamIndexFilesJava=self.BuildBamIndexFilesJava, 
            BuildBamIndexJar=self.BuildBamIndexJar, \
            mv=self.mv, data_dir=self.data_dir, tmpDir=self.tmpDir)

        passingData.flagStatMapFolderJob = self.addMkDirJob(\
            outputDir="%sFlagStatMap"%(passingData.outputDirPrefix))

        return returnData

    def mapEachAlignment(self, alignmentData=None, passingData=None, \
        ransferOutput=True, **keywords):
        """
        2012.9.22
            similar to reduceBeforeEachAlignmentData() but for mapping
             programs that run on one alignment each.

            passingData.alignmentJobAndOutputLs = []
            passingData.bamFnamePrefix = bamFnamePrefix
            passingData.individual_alignment = alignment
        """
        returnData = PassingData(no_of_jobs = 0)
        returnData.jobDataLs = []

        topOutputDirJob = passingData.topOutputDirJob
        flagStatMapFolderJob = passingData.flagStatMapFolderJob

        refFastaF = passingData.refFastaFList[0]

        alignment = alignmentData.alignment
        parentJobLs = alignmentData.jobLs + [passingData.fastaDictJob, \
            passingData.fastaIndexJob]
        bamF = alignmentData.bamF
        baiF = alignmentData.baiF

        bamFileBasenamePrefix = alignment.getReadGroup()

        #4X coverage alignment => 120 minutes
        realInputVolume = alignment.individual_sequence.coverage
        if realInputVolume is None:	#default is 8X
            realInputVolume = 8
        jobWalltime = self.scaleJobWalltimeOrMemoryBasedOnInput(
            realInputVolume=realInputVolume, \
            baseInputVolume=4, baseJobPropertyValue=120, \
            minJobPropertyValue=60, maxJobPropertyValue=1200).value
        #base is 4X, => 5000M
        #2013.07.26 memory usage of DOC walker java job does not depend on 
        # coverage so much, more on ref genome size.
        jobMaxMemory = self.scaleJobWalltimeOrMemoryBasedOnInput(
            realInputVolume=realInputVolume, \
            baseInputVolume=4, baseJobPropertyValue=20000, \
            minJobPropertyValue=20000, maxJobPropertyValue=96000).value

        if self.skipAlignmentWithStats and alignment.median_depth is not None \
            and alignment.mean_depth is not None and alignment.mode_depth is not None:
            pass
        else:
            """
            depthOutputFile = File(os.path.join(topOutputDirJob.output, '%s_DOC.tsv.gz'%(alignment.id)))
            samtoolsDepthJob = self.addSAMtoolsDepthJob(samtoolsDepth=self.samtoolsDepth, \
                samtools_path=self.samtools_path,\
                bamF=bamF, outputFile=depthOutputFile, baiF=baiF, \
                parentJobLs=[topOutputDirJob] + alignmentData.jobLs, \
                job_max_memory = 500, extraArguments="", \
                transferOutput=False)
            self.addRefFastaJobDependency(job=samtoolsDepthJob, refFastaF=passingData.refFastaF, \
                fastaDictJob=passingData.fastaDictJob, refFastaDictF=passingData.refFastaDictF,\
                fastaIndexJob = passingData.fastaIndexJob, refFastaIndexF=passingData.refFastaIndexF)
            meanMedianModeDepthFile = File(os.path.join(topOutputDirJob.output, \
                "%s_meanMedianModeDepth.tsv"%(alignment.id)))
            meanMedianModeDepthJob = self.addCalculateDepthMeanMedianModeJob(\
                executable=self.CalculateMedianMeanOfInputColumn, \
                inputFile=depthOutputFile, outputFile=meanMedianModeDepthFile, 
                alignmentID=alignment.id, fractionToSample=self.fractionToSample, \
                noOfLinesInHeader=0, whichColumn=2, maxNumberOfSamplings=1E6,\
                parentJobLs=[topOutputDirJob, samtoolsDepthJob], job_max_memory = 500, \
                extraArguments=None, \
                transferOutput=False)
            """
            #2013.05.17 samtools depth + CalculateMedianMeanOfInputColumn is
            #  not working well for realigned and BQSRed alignments
            # use GATK DOC walker
            DOCOutputFnamePrefix = os.path.join(topOutputDirJob.output, '%s_DOC'%(alignment.id))
            DOCJob = self.addDepthOfCoverageJob(DOCWalkerJava=self.DOCWalkerJava,
                refFastaFList=passingData.refFastaFList, bamF=bamF, baiF=baiF, \
                DOCOutputFnamePrefix=DOCOutputFnamePrefix,\
                parentJobLs=parentJobLs + [topOutputDirJob], \
                transferOutput=False,\
                job_max_memory = jobMaxMemory, walltime=jobWalltime)
                #1200 minutes is 20 hours
                #fractionToSample=self.fractionToSample, \
            depthOutputFile = DOCJob.sample_statistics_file
            meanMedianModeDepthFile = File(os.path.join(topOutputDirJob.output, \
                "%s_meanMedianModeDepth.tsv"%(alignment.id)))
            meanMedianModeDepthJob = self.addCalculateDepthMeanMedianModeJob(\
                executable=self.CalculateMedianMeanOfInputColumn, \
                inputFile=depthOutputFile, outputFile=meanMedianModeDepthFile, 
                alignmentID=alignment.id, \
                parentJobLs=[topOutputDirJob, DOCJob], job_max_memory = 500, 
                extraArguments="--inputFileFormat=2", \
                transferOutput=False)

            self.addInputToMergeJob(
                passingData.depthOfCoverageOutputMergeJob, 
                inputF=meanMedianModeDepthFile,\
                parentJobLs=[meanMedianModeDepthJob])
            self.no_of_alns_with_depth_jobs += 1

        if self.skipAlignmentWithStats and alignment.perc_reads_mapped is not None:
            pass
        else:
            #2013.05.17 GATK's flagstat, should be identical to samtools flagstat
            """
    java -jar ~/script/gatk2/GenomeAnalysisTK.jar -T FlagStat
    -I ~/3152_640_1985088_GA_vs_3280_by_method2_realigned1_reduced0_p2312_m87.bam
    -o ~/3152_640_1985088_GA_vs_3280_by_method2_realigned1_reduced0_p2312_m87.flagstat.txt
    --reference_sequence ~/NetworkData/vervet/db/individual_sequence/3280_vervet_ref_6.0.3.fasta

        output (<4 hours) looks like:

            1119300506 in total
            0 QC failure
            186159065 duplicates
            1034122354 mapped (92.39%)
            1119300506 paired in sequencing
            559647234 read1
            559653272 read2
            859005395 properly paired (76.74%)
            949042688 with itself and mate mapped
            85079666 singletons (7.60%)
            80245327 with mate mapped to a different chr
            26716310 with mate mapped to a different chr (mapQ>=5)

            """

            oneFlagStatOutputF = File(os.path.join(flagStatMapFolderJob.output, \
                '%s_flagstat.txt.gz'%(alignment.id)))
            #use jobMaxMemory to reduce the number of running jobs and IO load
            samtoolsFlagStatJob = self.addSamtoolsFlagstatJob(
                executable=self.samtoolsFlagStat, \
                samtoolsExecutableFile=self.samtoolsExecutableFile, 
                inputFile=bamF, outputFile=oneFlagStatOutputF, \
                parentJobLs=parentJobLs + [flagStatMapFolderJob], \
                extraDependentInputLs=[baiF], transferOutput=False, \
                extraArguments=None, job_max_memory=jobMaxMemory/2, \
                walltime=jobWalltime/2)
            self.addRefFastaJobDependency(job=samtoolsFlagStatJob, \
                refFastaF=passingData.refFastaF, \
                fastaDictJob=passingData.fastaDictJob, \
                refFastaDictF=passingData.refFastaDictF,\
                fastaIndexJob = passingData.fastaIndexJob, \
                refFastaIndexF=passingData.refFastaIndexF)
            reformatFlagStatOutputF = File(os.path.join(\
                flagStatMapFolderJob.output, '%s_flagstat.tsv'%(alignment.id)))
            reformatFlagStatOutputJob = self.addReformatFlagstatOutputJob(
                executable=self.ReformatFlagstatOutput, \
                inputF=oneFlagStatOutputF, alignmentID=alignment.id, 
                outputF=reformatFlagStatOutputF, \
                parentJobLs=[flagStatMapFolderJob, samtoolsFlagStatJob],
                extraDependentInputLs=[], \
                transferOutput=False, \
                extraArguments=None, job_max_memory=20, walltime=30)
            self.addInputToMergeJob(passingData.flagStatOutputMergeJob, \
                inputF=reformatFlagStatOutputJob.output, \
                parentJobLs=[reformatFlagStatOutputJob])
            self.no_of_alns_with_flagstat_jobs += 1

        if alignment.path_to_depth_file is None or not \
            os.path.isfile(os.path.join(self.data_dir, alignment.path_to_depth_file)):
            depthOutputFile = File(os.path.join(topOutputDirJob.output, \
                '%s_depth.tsv.gz'%(alignment.id)))
            #use jobMaxMemory to reduce the number of running jobs and IO load
            samtoolsDepthJob = self.addSAMtoolsDepthJob(samtoolsDepth=self.samtoolsDepth, \
                samtools_path=self.samtools_path,\
                bamF=bamF, outputFile=depthOutputFile, baiF=baiF, \
                parentJobLs=[topOutputDirJob] + alignmentData.jobLs, 
                job_max_memory = jobMaxMemory/2, \
                extraArguments=None, \
                transferOutput=False)
            self.addRefFastaJobDependency(job=samtoolsDepthJob, \
                refFastaF=passingData.refFastaF, \
                fastaDictJob=passingData.fastaDictJob, 
                refFastaDictF=passingData.refFastaDictF,\
                fastaIndexJob = passingData.fastaIndexJob, \
                refFastaIndexF=passingData.refFastaIndexF)

            logFile = File(os.path.join(passingData.reduceOutputDirJob.output, \
                "%s_depth_file_2DB.log"%(alignment.id)))
            outputFileRelativePath = "%s_depth.tsv.gz"%(os.path.splitext(alignment.path)[0])
            extraArgumentList = ["--db_entry_id %s"%(alignment.id), \
                "--tableClassName IndividualAlignment", \
                "--filePathColumnName path_to_depth_file",\
                "--fileSizeColumnName depth_file_size", \
                "--outputFileRelativePath %s"%(outputFileRelativePath), \
                "--data_dir %s"%(self.data_dir)]
            depthFile2DBJob = self.addPutStuffIntoDBJob(
                executable=self.AffiliateFile2DBEntry, \
                inputFile=samtoolsDepthJob.output, inputArgumentOption='-i',\
                logFile=logFile, commit=True, \
                parentJobLs=[samtoolsDepthJob, passingData.reduceOutputDirJob], \
                extraDependentInputLs=None, transferOutput=True, extraArguments=None, \
                extraArgumentList=extraArgumentList,\
                job_max_memory=10, sshDBTunnel=self.needSSHDBTunnel)
            pdata = self.constructJobDataFromJob(samtoolsDepthJob)
        else:
            alignmentDepthFile = self.registerOneInputFile(
                os.path.join(self.data_dir, alignment.path_to_depth_file), \
                input_site_handler=None, folderName=self.pegasusFolderName, \
                useAbsolutePathAsPegasusFileName=False,\
                pegasusFileName=None, checkFileExistence=True)
            pdata = PassingData(job=None, jobLs=[], file=alignmentDepthFile, \
                fileLs=[alignmentDepthFile])
        pdata.alignment = alignment
        self.alignmentDepthJobDataList.append(pdata)

        if self.needPerContigJob:	#need for per-contig job
            statOutputDir = 'perContigStatOfAlignment%s'%(alignment.id)
            passingData.statOutputDirJob = self.addMkDirJob(outputDir=statOutputDir)
        else:
            passingData.statOutputDirJob = None

        return returnData

    def mapEachChromosome(self, alignmentData=None, chromosome=None,\
        VCFJobData=None, passingData=None, reduceBeforeEachAlignmentData=None, \
        transferOutput=True, **keywords):
        """
        2012.9.17
        """
        returnData = PassingData(no_of_jobs = 0)
        returnData.jobDataLs = []
        if not self.needPerContigJob:
            #no need for per-contig job
            return returnData

        alignment = alignmentData.alignment

        parentJobLs = alignmentData.jobLs
        bamF = alignmentData.bamF
        baiF = alignmentData.baiF

        bamFnamePrefix = alignment.getReadGroup()

        statOutputDirJob = passingData.statOutputDirJob

        depthOutputFile = File(os.path.join(statOutputDirJob.output,
            '%s_%s_DOC.tsv.gz'%(alignment.id, chromosome)))
        samtoolsDepthJob = self.addSAMtoolsDepthJob(
            samtoolsDepth=self.samtoolsDepth,
            samtools_path=self.samtools_path,\
            bamF=bamF, outputFile=depthOutputFile, baiF=baiF,
            parentJobLs=[statOutputDirJob]+alignmentData.jobLs,
            job_max_memory = 500, extraArguments=None,
            transferOutput=False)
        self.addRefFastaJobDependency(job=samtoolsDepthJob,
            refFastaF=passingData.refFastaF,
            fastaDictJob=passingData.fastaDictJob,
            refFastaDictF=passingData.refFastaDictF,
            fastaIndexJob = passingData.fastaIndexJob, 
            refFastaIndexF=passingData.refFastaIndexF)
        meanMedianModeDepthFile = File(os.path.join(statOutputDirJob.output,
            "%s_%s_meanMedianModeDepth.tsv"%(alignment.id, chromosome)))
        meanMedianModeDepthJob = self.addCalculateDepthMeanMedianModeJob(
            executable=self.CalculateMedianMeanOfInputColumn, \
            inputFile=depthOutputFile, outputFile=meanMedianModeDepthFile, 
            alignmentID="%s-%s"%(alignment.id, chromosome), \
            parentJobLs=[samtoolsDepthJob], job_max_memory = 500,
            extraArguments="-r %s"%(chromosome),
            transferOutput=False)

        self.addInputToMergeJob(
            passingData.depthOfCoveragePerChrOutputMergeJob,
            inputF=meanMedianModeDepthFile,\
            parentJobLs=[meanMedianModeDepthJob])
        """
        DOCOutputFnamePrefix = os.path.join(statOutputDir, 
            '%s_%s_DOC'%(alignment.id, chromosome))
        DOCJob = self.addDepthOfCoverageJob(DOCWalkerJava=ContigDOCWalkerJava, \
                GenomeAnalysisTKJar=GenomeAnalysisTKJar,\
                refFastaFList=refFastaFList, bamF=bamF, baiF=baiF, \
                DOCOutputFnamePrefix=DOCOutputFnamePrefix,\
                parentJobLs=[statOutputDirJob]+alignmentData.jobLs, 
                job_max_memory = perContigJobMaxMemory*3,
                extraArguments="-L %s"%(chromosome),\
                transferOutput=False,\
                fractionToSample=self.fractionToSample)

        reduceDepthOfCoverageJob.addArguments(DOCJob.sample_statistics_file)
        reduceDepthOfCoverageJob.uses(DOCJob.sample_statistics_file, transfer=True,
            register=True, link=Link.INPUT)
        self.depends(parent=DOCJob, child=reduceDepthOfCoverageJob)
        """

        return returnData

    def reduceAfterEachAlignment(self, passingData=None,
        mapEachChromosomeDataLs=None,
        reduceAfterEachChromosomeDataLs=None,\
        transferOutput=True, **keywords):
        """
        """
        returnData = PassingData(no_of_jobs = 0)
        returnData.jobDataLs = []
        returnData.mapEachChromosomeDataLs = mapEachChromosomeDataLs
        returnData.reduceAfterEachChromosomeDataLs = reduceAfterEachChromosomeDataLs
        return returnData

    def reduce(self, passingData=None, reduceAfterEachAlignmentDataLs=None,
            transferOutput=True, **keywords):
        """
        2013.08.14 add 2DB jobs only when their input is not empty
        2012.9.17
        """
        returnData = PassingData(no_of_jobs = 0)
        returnData.jobDataLs = []
        returnData.reduceAfterEachAlignmentDataLs = reduceAfterEachAlignmentDataLs

        reduceOutputDirJob = passingData.reduceOutputDirJob

        if passingData.flagStatOutputMergeJob.inputLs:
            #the merge job's input is not empty or None
            flagStat2DBLogFile = File(os.path.join(reduceOutputDirJob.output, \
                "flagStat2DB.log"))
            flagStat2DBJob = self.addPutStuffIntoDBJob(
                executable=self.PutFlagstatOutput2DB, \
                inputFileList=[passingData.flagStatOutputMergeJob.output], \
                logFile=flagStat2DBLogFile, commit=True, \
                parentJobLs=[reduceOutputDirJob, passingData.flagStatOutputMergeJob],
                extraDependentInputLs=[], transferOutput=True,
                extraArguments=None, \
                job_max_memory=10, sshDBTunnel=self.needSSHDBTunnel)
        if passingData.depthOfCoverageOutputMergeJob.inputLs:
            DOC2DBLogFile = File(os.path.join(reduceOutputDirJob.output, "DOC2DB.log"))
            DOC2DBJob = self.addPutStuffIntoDBJob(
                executable=self.PutDOCOutput2DB,
                inputFileList=[passingData.depthOfCoverageOutputMergeJob.output],
                logFile=DOC2DBLogFile, commit=True, \
                parentJobLs=[reduceOutputDirJob, \
                    passingData.depthOfCoverageOutputMergeJob],
                extraDependentInputLs=[], transferOutput=True,
                extraArguments=None,
                job_max_memory=10, sshDBTunnel=self.needSSHDBTunnel)
        if self.alignmentDepthJobDataList and \
            self.alignmentDepthIntervalMethodShortName:
            if not self.min_segment_length:
                logging.error(f"alignmentDepthIntervalMethodShortName="
                    f"{self.alignmentDepthIntervalMethodShortName} is given but"
                    f" min_segment_length ({self.min_segment_length}) is not.")
                sys.exit(4)
            alignmentIDList = [pdata.alignment.id for pdata in self.alignmentDepthJobDataList]
            alignmentIDListInStr = utils.getSuccinctStrOutOfList(alignmentIDList)
            #job to add an AlignmentDepthIntervalMethod
            logFile = File(os.path.join(self.logOutputDirJob.output,
                'AddAlignmentDepthIntervalMethod2DB.log'))
            addMethod2DBJob = self.addData2DBJob(
                executable=self.AddAlignmentDepthIntervalMethod2DB,
                inputFile=None, inputArgumentOption="-i", \
                outputFile=None, outputArgumentOption="-o", \
                data_dir=self.data_dir, logFile=logFile, commit=True,
                parentJobLs=[self.logOutputDirJob], \
                extraDependentInputLs=None, extraOutputLs=None,
                transferOutput=True, extraArguments=None, \
                extraArgumentList=["--methodShortName %s"%(
                    self.alignmentDepthIntervalMethodShortName),
                    "--alignmentIDList %s"%(alignmentIDListInStr),
                    "--min_segment_length %s"%(self.min_segment_length)],
                job_max_memory=2000, walltime=30,
                sshDBTunnel=self.needSSHDBTunnel)

            logFile = File(os.path.join(self.logOutputDirJob.output,
                'updateMethodNoOfIntervals.log'))
            updateMethodNoOfIntervalsJob = self.addData2DBJob(
                executable=self.UpdateAlignmentDepthIntervalMethodNoOfIntervals,
                data_dir=self.data_dir, logFile=logFile, commit=True,\
                parentJobLs=[self.logOutputDirJob],
                extraDependentInputLs=None, extraOutputLs=None,
                transferOutput=True, extraArguments=None, \
                extraArgumentList=["--methodShortName %s"%\
                    (self.alignmentDepthIntervalMethodShortName) ],
                job_max_memory=2000, walltime=30,
                sshDBTunnel=self.needSSHDBTunnel)

            for chromosome, chromosomeSize in self.chr2size.items():
                #add a ReduceSameChromosomeAlignmentDepthFiles job
                outputFile = File(os.path.join(reduceOutputDirJob.output, \
                    '%s_alignments_chr_%s_depth.tsv.gz'%(\
                        len(self.alignmentDepthJobDataList), chromosome)))
                reduceSameChromosomeAlignmentDepthFilesJob = self.addGenericJob(
                    executable=self.ReduceSameChromosomeAlignmentDepthFiles, \
                    inputFile=None, outputFile=outputFile, \
                    parentJobLs=[reduceOutputDirJob],
                    extraDependentInputLs=None, \
                    extraArgumentList=[f"-w 2 --chromosomePositionColumnIndex 1",
                        f" --chromosomeSize {chromosomeSize}"],
                    extraOutputLs=None,
                    transferOutput=False, \
                    key2ObjectForJob=None, job_max_memory=2000, walltime=60)
                for alignmentDepthJobData in self.alignmentDepthJobDataList:
                    #add a chromosome selection job
                    outputFile = File(os.path.join(passingData.topOutputDirJob.output,
                        '%s_chr_%s.tsv.gz'%(utils.getFileBasenamePrefixFromPath(\
                            alignmentDepthJobData.file.name), chromosome)))
                    selectRowsFromMatrixCCJob = self.addGenericJob(
                        executable=self.SelectRowsFromMatrixCC,
                        inputFile=alignmentDepthJobData.file,
                        outputFile=outputFile,
                        parentJobLs=alignmentDepthJobData.jobLs + \
                            [passingData.topOutputDirJob],
                        extraDependentInputLs=None, \
                        extraArgumentList=["--inputFileSortMode 1 -w 0",
                            " --whichColumnValue %s"%(chromosome)],
                        extraOutputLs=None,\
                        transferOutput=False,
                        key2ObjectForJob=None, job_max_memory=1000, walltime=60)
                    self.addInputToMergeJob(
                        reduceSameChromosomeAlignmentDepthFilesJob,
                        inputF=selectRowsFromMatrixCCJob.output, \
                        inputArgumentOption="-i",
                        parentJobLs=[selectRowsFromMatrixCCJob], \
                        extraDependentInputLs=None)
                #add GADA job
                # add segmentation jobs to figure out intervals at similar
                outputFile = File(os.path.join(reduceOutputDirJob.output, \
                    '%s_alignments_%s_depth_GADAOut_minSegLength%s.tsv.gz'%\
                    (len(self.alignmentDepthJobDataList), chromosome, 
                        self.min_segment_length)))
                #adjust memory based on chromosome size, 135Mb => 21.4g memory
                realInputVolume = chromosomeSize
                jobWalltime = self.scaleJobWalltimeOrMemoryBasedOnInput(
                    realInputVolume=realInputVolume, \
                    baseInputVolume=60000000, baseJobPropertyValue=600, \
                    minJobPropertyValue=60, maxJobPropertyValue=2400).value
                #base is 135M, => 21G
                jobMaxMemory = self.scaleJobWalltimeOrMemoryBasedOnInput(
                    realInputVolume=realInputVolume, \
                    baseInputVolume=135000000, baseJobPropertyValue=25000, \
                    minJobPropertyValue=11000, maxJobPropertyValue=29000).value
                GADAJob = self.addGenericJob(executable=self.GADA, \
                    inputFile=reduceSameChromosomeAlignmentDepthFilesJob.output, 
                    outputFile=outputFile, \
                    parentJobLs=[reduceOutputDirJob, \
                        reduceSameChromosomeAlignmentDepthFilesJob], 
                    extraDependentInputLs=None, \
                    extraArgumentList=[f"--MinSegLen {self.min_segment_length}",
                        '--debug -T 10 -a 0.5'], extraOutputLs=None,\
                    transferOutput=False, \
                    key2ObjectForJob=None, job_max_memory=jobMaxMemory,
                        walltime=jobWalltime)
                """
                GADAJob = self.addGenericJob(executable=self.GADA, \
                    inputFile=reduceSameChromosomeAlignmentDepthFilesJob.output, 
                    outputFile=outputFile, \
                    parentJobLs=[reduceOutputDirJob, \
                        reduceSameChromosomeAlignmentDepthFilesJob], \
                    extraDependentInputLs=None, \
                    extraArgumentList=["-M %s"%(self.min_segment_length)],
                    extraOutputLs=None,\
                    transferOutput=False, \
                    key2ObjectForJob=None, job_max_memory=10000, walltime=200)
                """
                #job that adds AlignmentDepthIntervalFile
                logFile = File(os.path.join(self.logOutputDirJob.output, \
                    'AddAlignmentDepthIntervalFile2DB_chr_%s.log'%(chromosome)))
                addFile2DBJob = self.addData2DBJob(
                    executable=self.AddAlignmentDepthIntervalFile2DB, \
                    inputFile=GADAJob.output, \
                    inputArgumentOption="-i", \
                    inputFileList=None, argumentForEachFileInInputFileList=None,\
                    outputFile=None, outputArgumentOption="-o", \
                    data_dir=self.data_dir, logFile=logFile, commit=True,\
                    parentJobLs=[GADAJob, addMethod2DBJob, self.logOutputDirJob],
                    extraDependentInputLs=None, extraOutputLs=None,
                    transferOutput=True,
                    extraArguments=None,
                    extraArgumentList=[
                        "--methodShortName %s"%(self.alignmentDepthIntervalMethodShortName),
                        "--alignmentIDList %s"%(alignmentIDListInStr), 
                        '--chromosome %s'%(chromosome),\
                        "--format tsv"], \
                    job_max_memory=2000, walltime=30, sshDBTunnel=self.needSSHDBTunnel)
                self.depends(parent=addFile2DBJob, child=updateMethodNoOfIntervalsJob)
        print(f" {self.no_of_jobs} jobs, {self.no_of_alns_with_depth_jobs}"
            f" alignments with depth jobs, {self.no_of_alns_with_flagstat_jobs} "
            f"alignments with flagstat jobs.", flush=True)
        return returnData

    def registerExecutables(self):
        """
        """
        ParentClass.registerExecutables(self)
        self.registerOneExecutable(path=self.javaPath, \
            name='ContigDOCWalkerJava', clusterSizeMultiplier=1)
        self.registerOneExecutable(path=self.javaPath, \
            name='ContigVariousReadCountJava', clusterSizeMultiplier=1)
        self.registerOneExecutable(path=os.path.join(self.pymodulePath, 'GADA/GADA'),
            name='GADA', clusterSizeMultiplier=0.1)
        self.registerOneExecutable(path=os.path.join(self.pymodulePath,
            'db/import/AddAlignmentDepthIntervalMethod2DB.py'),
            name='AddAlignmentDepthIntervalMethod2DB', clusterSizeMultiplier=0)
        self.registerOneExecutable(path=os.path.join(self.pymodulePath,
            'db/import/AddAlignmentDepthIntervalFile2DB.py'),
            name='AddAlignmentDepthIntervalFile2DB', clusterSizeMultiplier=0.3)
        self.registerOneExecutable(path=os.path.join(self.pymodulePath,
             'db/import/UpdateAlignmentDepthIntervalMethodNoOfIntervals.py'),
            name='UpdateAlignmentDepthIntervalMethodNoOfIntervals',
            clusterSizeMultiplier=0)
        self.registerOneExecutable(path=os.path.join(self.pymodulePath,
             'db/import/AffiliateFile2DBEntry.py'), \
            name='AffiliateFile2DBEntry', clusterSizeMultiplier=0.1)
        """
        self.registerOneExecutable(path=os.path.join(self.pymodulePath, 
            'reducer/ReduceDepthOfCoverage.py'), \
            name='ReduceDepthOfCoverage', clusterSizeMultiplier=0)

        self.registerOneExecutable(path=os.path.join(self.pymodulePath, 
            'reducer/ReduceVariousReadCount.py'), \
            name='ReduceVariousReadCount', clusterSizeMultiplier=0)
        """
        self.registerOneExecutable(path=os.path.join(self.pymodulePath,
                'reducer/ReduceSameChromosomeAlignmentDepthFiles'),
            name='ReduceSameChromosomeAlignmentDepthFiles', clusterSizeMultiplier=0.5)

        self.registerOneExecutable(path=os.path.join(self.pymodulePath, 
            "mapper/converter/ReformatFlagstatOutput.py"), \
            name='ReformatFlagstatOutput', clusterSizeMultiplier=1)
        self.registerOneExecutable(path=os.path.join(self.pymodulePath, 
            "polymorphism/mapper/samtoolsDepth.sh"), \
            name='samtoolsDepth', clusterSizeMultiplier=0.1)
        self.registerOneExecutable(path=os.path.join(self.pymodulePath, 
            "polymorphism/mapper/CalculateMedianModeFromSAMtoolsDepthOutput.py"), \
            name='CalculateMedianModeFromSAMtoolsDepthOutput', clusterSizeMultiplier=1)
        self.registerOneExecutable(path=os.path.join(self.pymodulePath, 
            'db/import/PutFlagstatOutput2DB.py'), \
            name='PutFlagstatOutput2DB', clusterSizeMultiplier=0)
        self.registerOneExecutable(path=os.path.join(self.pymodulePath, 
            'db/import/PutDOCOutput2DB.py'), \
            name='PutDOCOutput2DB', clusterSizeMultiplier=0)
        self.registerOneExecutable(path=os.path.join(self.pymodulePath, 
            'shell/pipe2File.sh'), \
            name='samtoolsFlagStat', clusterSizeMultiplier=1)
        
        self.registerOneExecutable(path=os.path.join(self.pymodulePath, \
            "mapper/CalculateMedianMeanOfInputColumn"), 
            name='CalculateMedianMeanOfInputColumn', clusterSizeMultiplier=1)
        #2013.08.23 c++ version of SelectRowsFromMatrix.py
        self.registerOneExecutable(path=os.path.join(self.pymodulePath, 
            'mapper/extractor/SelectRowsFromMatrixCC'), \
            name='SelectRowsFromMatrixCC', clusterSizeMultiplier=1)
        self.registerOneExecutable(path=os.path.join(self.pymodulePath, 
            'mapper/extractor/SelectRowsFromMatrix.py'), \
            name='SelectRowsFromMatrix', clusterSizeMultiplier=1)


if __name__ == '__main__':
    from argparse import ArgumentParser
    ap = ArgumentParser()
    ap.add_argument("--drivername", default="postgresql",
        help='Type of database server (default: %(default)s)')
    ap.add_argument('-z', "--hostname", default="pdc",
        help='name/IP of database server (default: %(default)s)')
    ap.add_argument("--port", default=None,
        help="database server port (default: %(default)s)")
    ap.add_argument("--dbname", default='pmdb',
        help="database name (default: %(default)s)")
    ap.add_argument('-k', "--schema", default='sunset', 
        help="database schema (default: %(default)s)")
    ap.add_argument("-u", "--db_user", help="Database user")
    ap.add_argument("-p", "--db_passwd", required=False,
        help="Password of the database user")

    ap.add_argument('-i', "--ind_aln_id_ls", required=True,
        help='a comma/dash-separated list of IndividualAlignment.id.')
    ap.add_argument("--ind_seq_id_ls",
        help='a comma/dash-separated list of IndividualSequence.id.')
    ap.add_argument("--alignment_method_id", type=int, 
        help='AlignmentMethod.id, to filter alignments.')
    ap.add_argument("--excludeContaminant", action='store_true',
        help='Toggle to exclude sequences from contaminated individuals, '
            '(IndividualSequence.is_contaminated=1)')
    ap.add_argument("--sequence_filtered", action='store_true',
        help='Filter alignments/individual_sequences. '
            'None: no filter, 0: unfiltered sequences, '
            '1: filtered sequences: 2: ...')
    
    ap.add_argument("--min_segment_length", type=int, default=100,
        help='minimal length in segmenting the alignment depth data.')
    ap.add_argument("--needPerContigJob", action='store_true',
        help='Toggle to add DepthOfCoverage and VariousReadCount jobs '
        'for each contig')
    ap.add_argument("--skipAlignmentWithStats", action='store_true',
        help='If an alignment has depth stats filled, no DOC job will be run. '
        'Similar for flagstat job.')
    ap.add_argument("--alignmentDepthIntervalMethodShortName",
        action='store_true',
        help='AlignmentDepthIntervalMethod.short_name, '
        'to store segmented depth intervals from all alignments into db. '
        'This will trigger depth-per-segment jobs.')

    ap.add_argument("-F", "--pegasusFolderName", default='input',
        help='The path relative to the workflow running root. '
        'This folder will contain pegasus input & output. '
        'It will be created during the pegasus staging process. '
        'It is useful to separate multiple sub-workflows. '
        'If empty or None, everything is in the pegasus root.')
    ap.add_argument("-l", "--site_handler", type=str, required=True,
        help="The name of the computing site where the jobs run and "
        "executables are stored. "
        "Check your Pegasus configuration in submit.sh.")
    ap.add_argument("-j", "--input_site_handler", type=str,
        help="It is the name of the site that has all the input files."
        "Possible values can be 'local' or same as site_handler."
        "If not given, it is asssumed to be the same as site_handler and "
        "the input files will be symlinked into the running folder."
        "If input_site_handler=local, the input files will be transferred "
        "to the computing site by pegasus-transfer.")
    ap.add_argument("-C", "--cluster_size", type=int, default=30,
        help="Default: %(default)s. "
        "This number decides how many of pegasus jobs should be clustered "
        "into one job. "
        "Good if your workflow contains many quick jobs. "
        "It will reduce Pegasus monitor I/O.")
    ap.add_argument("-o", "--output_path", type=str, required=True,
        help="The path to the output file that will contain the Pegasus DAG.")
    
    ap.add_argument("--tmpDir", type=str, default='/tmp/',
        help='Default: %(default)s. '
        'A local folder for some jobs (MarkDup) to store temp data. '
        '/tmp/ can be too small for high-coverage sequencing.')
    ap.add_argument("--max_walltime", type=int, default=4320,
        help='Default: %(default)s. '
        'Maximum wall time for any job, in minutes. 4320=3 days. '
        'Used in addGenericJob(). Most clusters have upper limit for runtime.')
    
    ap.add_argument("--home_path",
        help="Path to your home folder. Default is ~.")
    ap.add_argument("--javaPath", default='bin/java',
        help="Path to java. Default is %(default)s.")
    ap.add_argument("--pymodulePath", type=str, default="src/pymodule",
        help="Path to the pymodule code folder. "
        "If relative path, home folder is inserted in the front.")
    
    ap.add_argument("--needSSHDBTunnel", action='store_true',
        help="If all DB-interacting jobs need a ssh tunnel to "
        "access a database that is inaccessible to computing nodes.")
    ap.add_argument("-c", "--commit", action='store_true',
        help="Toggle to commit the db transaction (default: %(default)s)")
    ap.add_argument("--debug", action='store_true',
        help='Toggle debug mode.')
    ap.add_argument("--report", action='store_true',
        help="Toggle verbose mode. Default: %(default)s.")
    args = ap.parse_args()
    if not args.db_user:
        args.db_user = getpass.getuser()
    if not args.db_passwd:
        args.db_passwd = getpass.getpass(f"Password for {args.db_user}:")
    
    instance = InspectAlignmentPipeline(
        drivername = args.drivername,
        hostname = args.hostname,
        port = args.port,
        dbname = args.dbname, schema = args.schema,
        db_user = args.db_user, db_passwd = args.db_passwd,

        ind_aln_id_ls = args.ind_aln_id_ls,
        ind_seq_id_ls = args.ind_seq_id_ls,
        alignment_method_id = args.alignment_method_id,
        excludeContaminant = args.excludeContaminant,
        sequence_filtered = args.sequence_filtered,

        min_segment_length = args.min_segment_length,
        needPerContigJob = args.needPerContigJob,
        skipAlignmentWithStats = args.skipAlignmentWithStats,
        alignmentDepthIntervalMethodShortName = \
            args.alignmentDepthIntervalMethodShortName,
        
        pegasusFolderName = args.pegasusFolderName,
        site_handler = args.site_handler, 
        input_site_handler = args.input_site_handler,
        cluster_size = args.cluster_size,
        output_path = args.output_path,
        tmpDir = args.tmpDir,
        max_walltime = args.max_walltime,

        home_path = args.home_path,
        javaPath = args.javaPath,
        pymodulePath = args.pymodulePath,
        needSSHDBTunnel = args.needSSHDBTunnel,
        commit = args.commit,
        debug = args.debug,
        report=args.report)
    instance.run()
