#!/usr/bin/env python3
"""
A read-alignment workflow.

Examples:
    %s
"""
import sys, os, math
__doc__ = __doc__%(sys.argv[0])

import copy
from pegaflow.DAX3 import Executable, File, PFN, Link, Job
from palos import ProcessOptions, getListOutOfStr, PassingData, utils
import pegaflow
from palos.ngs.AbstractNGSWorkflow import AbstractNGSWorkflow
from palos.ngs.AbstractAlignmentAndVCFWorkflow import AbstractAlignmentAndVCFWorkflow
from . AlignmentReadBaseQualityRecalibrationWorkflow import AlignmentReadBaseQualityRecalibrationWorkflow
ParentClass = AbstractNGSWorkflow

class ShortRead2AlignmentWorkflow(ParentClass, AlignmentReadBaseQualityRecalibrationWorkflow):
    __doc__ = __doc__
    option_default_dict = copy.deepcopy(
        AlignmentReadBaseQualityRecalibrationWorkflow.option_default_dict)
    alignment_option_dict = {
        ('noCheckEmptyReadFile', 0, int):[0, '', 0, 
            "toggle to not check whether each read file is empty (if empty, exclude it). "
            "If IndividualSequenceFile.read_count is null, "
            "it'll try to count them on the fly and take a lot of time. "
            "However, only toggle it if you know every input individual_sequence_file is not empty. "
            "Empty read file fails alignment jobs."],\
        ('additionalArguments', 0, ): ["", '', 1, 
            'a string of additional arguments passed to aln, '
            'not bwasw, add double quote if empty space. i.e. "-q 20"'],\
        ("bwa_path", 1, ): ["%s/bin/bwa", '', 1, 'bwa binary'],\
        ("alignment_method_name", 1, ): ["bwaShortRead", '', 1, 
            'alignment_method.short_name from db. '
            'used only when unable to guess based on individual_sequence.sequencer '
            'and individual_sequence.sequence_type'],\
        ("needRefIndexJob", 0, int): [0, '', 0, 
            'need to add a reference index job by bwa?'],\
        ('no_of_aln_threads', 1, int): [1, '', 1, 
            'number of threads during alignment'],\
        ('alignmentJobClusterSizeFraction', 1, float): [0.01, '', 1, 
            'alignment job cluster size relative to self.cluster_size, '
            'for bwa, PEAlignmentByBWA, LongSEAlignmentByBWA, '
            'AddOrReplaceReadGroupsJava/SortSamFilesJava/samtools jobs'],\
        ('notStageOutFinalOutput', 0, int):[0, '', 0, 
            'toggle to not stage out final output (bam + bam.bai)'],\
        ("coreAlignmentJobWallTimeMultiplier", 1, float): [0.2, '', 1, 
            'The default wall time is 23 hours (for bwa aln, stampy, etc.). '
            'This option controls walltime by multiplying the default with this number.'],
        ("alignmentPerLibrary", 0, int): [0, '', 0, 
            'toggle to run alignment for each library of one individual_sequence'],
        }
    option_default_dict.update(alignment_option_dict.copy())
    #by default not to apply known-indel-free local re-alignment
    option_default_dict[('local_realigned', 0, int)][0] = 0


    def __init__(self,  **keywords):
        """
        2012.3.29
            default to stage out final output.
            Argument stageOutFinalOutput morphs into notStageOutFinalOutput.
        2011-7-11
        """
        ParentClass.__init__(self, **keywords)
        #ParentClass.__init__(self, **keywords)
        if self.ind_seq_id_ls:
            self.ind_seq_id_ls = getListOutOfStr(self.ind_seq_id_ls, data_type=int)
        self.bwa_path =  self.insertHomePath(self.bwa_path, self.home_path)
        if self.notStageOutFinalOutput:
            self.stageOutFinalOutput = False
        else:
            self.stageOutFinalOutput = True

        if self.noCheckEmptyReadFile:
            self.ignoreEmptyReadFile = False
        else:
            self.ignoreEmptyReadFile = True
    
    def registerCustomExecutables(self):
        """
        """
        ParentClass.registerCustomExecutables(self)
        AlignmentReadBaseQualityRecalibrationWorkflow.registerCustomExecutables(self)

        #self.BuildBamIndexFilesJava.addProfile(Profile(Namespace.PEGASUS, 
        # key="clusters.size", value="%s"%self.alignmentJobClusterSizeFraction))
        self.setExecutableClusterSize(executable=self.AddOrReplaceReadGroupsJava,
                clusterSizeMultiplier=self.alignmentJobClusterSizeFraction, \
                defaultClusterSize=None)
        self.setExecutableClusterSize(executable=self.samtools,
                clusterSizeMultiplier=self.alignmentJobClusterSizeFraction, \
                defaultClusterSize=None)
        self.registerOneExecutable(path=self.bwa_path, \
            name="bwa", clusterSizeMultiplier=self.alignmentJobClusterSizeFraction)
        self.registerOneExecutable(
            path=os.path.join(self.pymodulePath, 'mapper/alignment/ShortSEAlignmentByBWA.sh'), \
            name="ShortSEAlignmentByBWA",
            clusterSizeMultiplier=self.alignmentJobClusterSizeFraction)
        self.registerOneExecutable(
            path=os.path.join(self.pymodulePath, 'mapper/alignment/PEAlignmentByBWA.sh'), \
            name="PEAlignmentByBWA",
            clusterSizeMultiplier=self.alignmentJobClusterSizeFraction)
        self.registerOneExecutable(
            path=os.path.join(self.pymodulePath, 'mapper/alignment/LongSEAlignmentByBWA.sh'), \
            name="LongSEAlignmentByBWA",
            clusterSizeMultiplier=self.alignmentJobClusterSizeFraction)
        self.registerOneExecutable(path=self.javaPath, \
            name="SAM2BAMJava", clusterSizeMultiplier=self.alignmentJobClusterSizeFraction)
        #self.registerOneExecutable(
        #	path=os.path.join(self.pymodulePath, 'pegasus/mapper/alignment/BWA_Mem.sh'), \
        #	name="BWA_Mem", clusterSizeMultiplier=self.alignmentJobClusterSizeFraction)
        # 2014.04.04 use generic pipe2File.sh instead
        self.registerOneExecutable(
            path=os.path.join(self.pymodulePath, 'shell/pipe2File.sh'), \
            name="BWA_Mem", clusterSizeMultiplier=self.alignmentJobClusterSizeFraction)
        # 2014.04.04 pipe2File need this file dependency
        self.registerOneExecutableAsFile(pythonVariableName="bwaExecutableFile", path=self.bwa_path)
        self.registerOneExecutable(
            path=os.path.join(self.pymodulePath, "db/import/AddAlignmentFile2DB.py"), \
            name="AddAlignmentFile2DB",
            clusterSizeMultiplier=self.alignmentJobClustersSizeFraction)
        #self.registerOneExecutable(path=self.javaPath, name='RealignerTargetCreatorJava', 
        # clusterSizeMultiplier=0.7)
        #self.registerOneExecutable(path=self.javaPath, name='IndelRealignerJava', 
        # clusterSizeMultiplier=0.2)

    def addBWAReferenceIndexJob(self, refFastaFList=None, \
        refSequenceBaseCount=3000000000, bwa=None,\
        bwaIndexFileSuffixLs = None,\
        transferOutput=True, job_max_memory=4000, walltime=200):
        """
        renamed from addRefIndexJob to addBWAReferenceIndexJob()
        bwaIndexFileSuffixLs is controlled by ParentClass.bwaIndexFileSuffixLs.
        and the suffices of index output has changed ('rbwt', 'rpac', 'rsa', are gone):
            524_superContigsMinSize2000.fasta.bwt
            524_superContigsMinSize2000.fasta.pac
            524_superContigsMinSize2000.fasta.ann
            524_superContigsMinSize2000.fasta.amb
            524_superContigsMinSize2000.fasta.sa
            524_superContigsMinSize2000.fasta.nhr
            524_superContigsMinSize2000.fasta.nin
            524_superContigsMinSize2000.fasta.nsq
        """
        if bwaIndexFileSuffixLs is None:
            bwaIndexFileSuffixLs = self.bwaIndexFileSuffixLs

        if refSequenceBaseCount is None:
            #default
            index_algorithm = 'bwtsw'
        elif refSequenceBaseCount<500000000:
            #500 million
            index_algorithm = "is"
        else:
            index_algorithm = "bwtsw"

        refFastaFile = refFastaFList[0]
        extraArgumentList = ["index", "-a", index_algorithm, refFastaFile]
        extraOutputLs = []
        for suffix in bwaIndexFileSuffixLs:
            file = File("%s.%s"%(refFastaFile.name, suffix))
            extraOutputLs.append(file)
        extraDependentInputLs = refFastaFList
        bwa_index_job = self.addGenericJob(executable=bwa, inputFile=None,
            outputFile=None, \
            parentJobLs=None, extraDependentInputLs=extraDependentInputLs, \
            extraOutputLs=extraOutputLs,\
            transferOutput=transferOutput, \
            extraArguments=None, extraArgumentList=extraArgumentList, \
            key2ObjectForJob=None,\
            job_max_memory=job_max_memory, walltime=walltime)
        return bwa_index_job

    def addStampyGenomeIndexHashJob(self, executable=None, refFastaFList=None, \
        parentJobLs=None, job_max_memory=100, walltime = 60, \
        extraDependentInputLs=None, \
        transferOutput=True, **keywords):
        """
        2012.10.10 use addGenericJob()
        2012.2.23
        """
        """
        stampyGenomeIndexJob= Job(namespace=self.namespace, name=executable.name, version=self.version)
        refFastaFile = refFastaFList[0]

        genomeIndexFile = File('%s.stidx'%(refFastaFile.name))
        stampyGenomeIndexJob.addArguments("-G ", refFastaFile, refFastaFile)
        stampyGenomeIndexJob.uses(refFastaFile, transfer=True, register=True, link=Link.INPUT)
        #genomeIndexFile is output of this stampyGenomeIndexJob but ignore it as output 
        # otherwise it'll get auto-cleaned by pegasus
        #stampyGenomeIndexJob.uses(genomeIndexFile, transfer=transferOutput, register=True, link=Link.OUTPUT)
        self.addJob(stampyGenomeIndexJob)
        for input in extraDependentInputLs:
            stampyGenomeIndexJob.uses(input, transfer=True, register=True, link=Link.INPUT)
        for parentJob in parentJobLs:
            self.depends(parent=parentJob, child=stampyGenomeIndexJob)
        pegaflow.setJobResourceRequirement(stampyGenomeIndexJob, job_max_memory=job_max_memory,
            walltime=walltime)

        stampyGenomeHashJob = Job(namespace=self.namespace, name=executable.name, version=self.version)
        stampyGenomeHashJob.addArguments("-g", refFastaFile, "-H", refFastaFile)

        stampyGenomeHashJob.uses(refFastaFile, transfer=True, register=True, link=Link.INPUT)
        genomeHashFile = File('%s.sthash'%(refFastaFile.name))
        #genomeIndexFile is input to this stampyGenomeHashJob 
        # but mark it as output otherwise it'll get auto-cleaned by pegasus
        stampyGenomeHashJob.uses(genomeIndexFile, transfer=transferOutput, register=True, link=Link.OUTPUT)
        stampyGenomeHashJob.uses(genomeHashFile, transfer=transferOutput, register=True, link=Link.OUTPUT)
        stampyGenomeHashJob.outputLs = [genomeIndexFile, genomeHashFile]

        self.addJob(stampyGenomeHashJob)
        self.depends(parent=stampyGenomeIndexJob, child=stampyGenomeHashJob)
        pegaflow.setJobResourceRequirement(stampyGenomeHashJob, job_max_memory=job_max_memory, 
            walltime=walltime)
        """

        refFastaFile = refFastaFList[0]
        genomeIndexFile = File('%s.stidx'%(refFastaFile.name))
        extraOutputLs = [genomeIndexFile]
        extraArgumentList = [refFastaFile]
        stampyGenomeIndexJob = self.addGenericJob(executable=executable,
            inputFile=refFastaFile,
            inputArgumentOption="-G", outputFile=None, \
            parentJobLs=parentJobLs,
            extraDependentInputLs=extraDependentInputLs,
            extraOutputLs=extraOutputLs,\
            transferOutput=transferOutput, \
            extraArguments=None, extraArgumentList=extraArgumentList,
            key2ObjectForJob=None,\
            job_max_memory=job_max_memory, walltime=walltime)

        genomeHashFile = File('%s.sthash'%(refFastaFile.name))
        extraOutputLs = [genomeHashFile, genomeIndexFile]
        extraArgumentList = ["-H", refFastaFile]
        stampyGenomeHashJob = self.addGenericJob(executable=executable,
            inputFile=refFastaFile,
            inputArgumentOption="-g", outputFile=None,
            parentJobLs=[stampyGenomeIndexJob],
            extraDependentInputLs=[genomeIndexFile],
            extraOutputLs=extraOutputLs,\
            transferOutput=transferOutput, \
            extraArguments=None, extraArgumentList=extraArgumentList, \
            key2ObjectForJob=None,\
            job_max_memory=job_max_memory, walltime=walltime)
        return stampyGenomeHashJob

    def registerISQFileObjLsToWorkflow(self, fileObjectLs=None, data_dir=None):
        '''
        2012-2.24
        '''

        newFileObjLs = []
        for fileObject in fileObjectLs:
            relativePath = fileObject.db_entry.path
            fastqF = File(relativePath)
            fastqF.addPFN(PFN("file://" + fileObject.path, self.input_site_handler))
            self.addFile(fastqF)
            fileObject.fastqF = fastqF
            newFileObjLs.append(fileObject)
        return newFileObjLs


    def addRefIndexJobAndItsOutputAsParent(self, refIndexJob=None, childJob=None):
        """
        2012.10.18
            updated by using addJobDependency() and self.addJobUse()
        2012.2.24
        """
        self.addJobDependency(parentJob=refIndexJob, childJob=childJob)
        for output in refIndexJob.outputLs:
            self.addJobUse(childJob, file=output, transfer=True,
                register=True, link=Link.INPUT)
            #childJob.uses(output, transfer=True, register=True, link=Link.INPUT)

    def addStampyAlignmentJob(self, fileObjectLs=None, \
        refFastaFList=None, bwa=None, additionalArguments=None, \
        samtools=None, refIndexJob=None, parentJobLs=None,\
        alignment_method=None, outputDir=None, \
        PEAlignmentByBWA=None, ShortSEAlignmentByBWA=None,
        LongSEAlignmentByBWA=None,
        java=None, SortSamFilesJava=None, SortSamJar=None,\
        AddOrReplaceReadGroupsJava=None, AddOrReplaceReadGroupsJar=None,\
        no_of_aln_threads=3, stampy=None, transferOutput=False, **keywords):
        """
        2014.4.4
        update for stampy 1.0.17
            which added option --gatkcigarworkaround removing 
            adjacent I/D events from CIGAR strings, which trips up GATK
        2012.3.5
        added "--overwrite" to stampy.py to overwrite partial output alignment file
        added two stampy options
            --baq          (SAM format) Compute base-alignment quality (BAQ; BQ tag)
            --alignquals   (SAM format) Compute posterior alignment probabilities (YQ tag)
        2012.2.26
            handle read files with quality_score_format="Illumina" (post 1.3 solexa).

        2012.2.24
            alignment job for stampy
            no_of_aln_threads is only for bwa.
        2011-9-13
            add argument java & SortSamJar
        2011-9-9
            two steps:
                1. aln doesn't use pipe, outputs to sai files.
                2. sampe/samse, convert, sort => connected through pipe
        """
        aln_job_max_memory = 6000	#in MB, bwa needs 3G. stampy needs 3G.
        #bwa: memory 2.5GB is enough for 4.5G gzipped fastq versus 480K contigs
        #  (total size~3G)

        bwasw_job_max_memory = 3800	#in MB, same ref, bwasw needs more memory
        samse_job_max_memory = 4500	#in MB
        sampe_job_max_memory = 6000	#in MB
        addRGJob_max_memory = 2500	#in MB
        #2013.3.1 change the walltime by multiplying it
        baseAlnJobWalltime = int(1380*self.coreAlignmentJobWallTimeMultiplier)
        #1380 is 23 hours, because all reads are stored in chunks of 5-million-read files

        refFastaFile = refFastaFList[0]
        firstFileObject = fileObjectLs[0]
        read_count = firstFileObject.db_entry.read_count
        aln_job_walltime = self.scaleJobWalltimeOrMemoryBasedOnInput(
            realInputVolume=read_count, \
            baseInputVolume=4000000, baseJobPropertyValue=baseAlnJobWalltime, \
            minJobPropertyValue=baseAlnJobWalltime*2/3,
            maxJobPropertyValue=baseAlnJobWalltime*5).value

        fastqF = firstFileObject.fastqF
        relativePath = fastqF.name
        fileBasenamePrefix = utils.getRealPrefixSuffixOfFilenameWithVariableSuffix(
            os.path.basename(relativePath))[0]
        outputSamFile = File('%s.sam'%(os.path.join(outputDir, fileBasenamePrefix)))

        alignmentJob = Job(namespace=self.namespace, name=stampy.name,
            version=self.version)
        # Make sure to use ', rather than ", 
        #  to wrap the bwaoptions. double-quote(") would disappear during xml translation.
        alignmentJob.addArguments(" --bwa=%s "%(pegaflow.getAbsPathOutOfExecutable(bwa)),
            "--bwaoptions='%s -t%s %s' "%(additionalArguments, 
                no_of_aln_threads, refFastaFile.name),
            "-g", refFastaFile, "-h", refFastaFile, "-o", outputSamFile,
            '--overwrite',
            '--baq', '--alignquals', '--gatkcigarworkaround')
        #Added option --gatkcigarworkaround removing adjacent I/D events
        #  from CIGAR strings, which trips up GATK
        if firstFileObject.db_entry.quality_score_format=='Illumina':
            alignmentJob.addArguments("--solexa")
        alignmentJob.uses(outputSamFile, transfer=transferOutput,
            register=True, link=Link.OUTPUT)
        alignmentJob.output = outputSamFile
        alignmentJob.fileBasenamePrefix = fileBasenamePrefix

        for refFastaFile in refFastaFList:
            alignmentJob.uses(refFastaFile, transfer=True, register=True, link=Link.INPUT)
        pegaflow.setJobResourceRequirement(alignmentJob, job_max_memory=aln_job_max_memory, \
                        no_of_cpus=no_of_aln_threads, walltime=aln_job_walltime)
        self.addJob(alignmentJob)
        #add fastq files after "-M"
        # -M FILE[,FILE] Map fastq file(s). 
        # Use FILE.recaldata for recalibration if available
        alignmentJob.addArguments(" -M ")
        for fileObject in fileObjectLs:
            fastqF = fileObject.fastqF
            relativePath = fastqF.name
            alignmentJob.addArguments(fastqF)
            alignmentJob.uses(fastqF, transfer=True, register=True, link=Link.INPUT)
        if refIndexJob:
            self.addRefIndexJobAndItsOutputAsParent(refIndexJob, childJob=alignmentJob)
        if parentJobLs:
            for parentJob in parentJobLs:
                if parentJob:
                    self.depends(parent=parentJob, child=alignmentJob)
        return alignmentJob

    def addBWAAlignmentJob(self, executable=None, bwaCommand='aln',
        fileObject=None, outputFile=None,\
        refFastaFList=None, no_of_aln_threads=3, \
        maxMissingAlignmentFraction=None, maxNoOfGaps=None,
        additionalArguments=None, \
        refIndexJob=None,\
        parentJobLs=None, extraDependentInputLs=None, extraOutputLs=None,
        transferOutput=False, \
        extraArguments=None, extraArgumentList=None, job_max_memory=2000, \
        key2ObjectForJob=None, walltime=None, \
        **keywords):
        """
        2012.10.10
        """
        """
        saiOutput = File(outputFname)
        alignmentJob = Job(namespace=namespace, name=bwa.name, version=version)
        alignmentJob.addArguments(alignment_method.command, 
            additionalArguments,"-t %s"%no_of_aln_threads, \
                "-f", saiOutput)
        if fileObject.db_entry.quality_score_format=='Illumina':	#2012.4.5
            alignmentJob.addArguments("-I")
        alignmentJob.addArguments(refFastaFList[0], fastqF)
        alignmentJob.uses(fastqF, transfer=True, register=True, link=Link.INPUT)
        for refFastaFile in refFastaFList:
            alignmentJob.uses(refFastaFile, transfer=True, register=True, link=Link.INPUT)
        alignmentJob.uses(saiOutput, transfer=False, register=True, link=Link.OUTPUT)
        pegaflow.setJobResourceRequirement(alignmentJob, job_max_memory=aln_job_max_memory, \
                    no_of_cpus=no_of_aln_threads, walltime=aln_job_walltime)

        self.addJob(alignmentJob)
        """
        if extraArgumentList is None:
            extraArgumentList = []
        extraArgumentList.append(bwaCommand)
        if additionalArguments:
            extraArgumentList.append(additionalArguments)
        if fileObject.db_entry.quality_score_format=='Illumina':
            extraArgumentList.append("-I")
        if no_of_aln_threads:
            extraArgumentList.append("-t %s"%no_of_aln_threads)
        if maxMissingAlignmentFraction is not None:
            extraArgumentList.append('-n %s'%(maxMissingAlignmentFraction))
        if maxNoOfGaps is not None:
            extraArgumentList.append("-o %s"%(maxNoOfGaps))
        extraArgumentList.extend(["-f", outputFile])
        extraArgumentList.append(refFastaFList[0])
        extraArgumentList.append(fileObject.fastqF)

        if extraDependentInputLs is None:
            extraDependentInputLs = []
        extraDependentInputLs.extend(refFastaFList[:])
        extraDependentInputLs.append(fileObject.fastqF)

        if extraOutputLs is None:
            extraOutputLs = []
        extraOutputLs = [outputFile]

        alignmentJob = self.addGenericJob(executable=executable,
            inputFile=None, outputFile=None, \
            parentJobLs=parentJobLs, extraDependentInputLs=extraDependentInputLs, \
            extraOutputLs=extraOutputLs,\
            transferOutput=transferOutput, \
            extraArguments=extraArguments, extraArgumentList=extraArgumentList, \
            key2ObjectForJob=None,\
            job_max_memory=job_max_memory, no_of_cpus=no_of_aln_threads, walltime=walltime)
        if refIndexJob:
            self.addRefIndexJobAndItsOutputAsParent(refIndexJob, childJob=alignmentJob)
        return alignmentJob

    def addBWAAlignmentWrapJob(self, fileObjectLs=None, \
        refFastaFList=None, bwa=None, additionalArguments=None, \
        samtools=None, \
        refIndexJob=None, parentJobLs=None,\
        alignment_method=None, outputDir=None, namespace=None, version=None,\
        PEAlignmentByBWA=None, ShortSEAlignmentByBWA=None, LongSEAlignmentByBWA=None, \
        java=None, 
        AddOrReplaceReadGroupsJava=None, AddOrReplaceReadGroupsJar=None,\
        no_of_aln_threads=3, maxMissingAlignmentFraction=None, maxNoOfGaps=None, \
        transferOutput=False, **keywords):
        """
        2012.10.10
        added argument maxMissingAlignmentFraction, maxNoOfGaps
        each fileObject in fileObjectLs should have 2 attributes
            fastqF: registered pegasus file
            db_entry: an individual_sequence_file db_entry or equivalent object. 2 attributes:
                quality_score_format: 'Standard', 'Illumina'
                individual_sequence: 2 attributes:
                    sequencer:  'GA' , '454'
                    sequence_type: 'PE', 'SR', 'genome'	#(not really used as 2012.10.10)
        alignment_method is also an object (supposedly from vervet db), two attributes.
            short_name (stampy, or bwa-short-read, or bwa-long-read)
            command (stampy.py , aln, bwasw)
        2012.4.5
            handle read files with quality_score_format="Illumina" (post 1.3 solexa).
            input fileObjectLs is different now.
        2012.2.23
            split out of addAlignmentJob()
        2011-9-13
            add argument java & SortSamJar
        2011-9-9
            two steps:
                1. aln doesn't use pipe, outputs to sai files.
                2. sampe/samse, convert, sort => connected through pipe
        """
        if namespace is None:
            namespace = self.namespace
        if version is None:
            version = self.version

        #in MB, 2.5GB is enough for 4.5G gzipped fastq versus 480K contigs (total size~3G)
        aln_job_max_memory = 2600
        bwasw_job_max_memory = 3800	#in MB, same ref, bwasw needs more memory
        samse_job_max_memory = 4500	#in MB
        sampe_job_max_memory = 6000	#in MB
        addRGJob_max_memory = 2500	#in MB

        #2013.3.1 change the walltime by multiplying it
        baseAlnJobWalltime = int(1380*self.coreAlignmentJobWallTimeMultiplier)
            #1380 is 23 hours, because all reads are stored in chunks of 5-million-read files

        memRequirementData = self.getJVMMemRequirment(
            job_max_memory=addRGJob_max_memory, minMemory=2000)
        job_max_memory = memRequirementData.memRequirement
        javaMemRequirement = memRequirementData.memRequirementInStr

        if len(fileObjectLs)==1:	#single end
            fileObject = fileObjectLs[0]
            fastqF = fileObject.fastqF
            relativePath = fastqF.name
            sequence_type = fileObject.db_entry.individual_sequence.sequence_type
            sequencer = fileObject.db_entry.individual_sequence.sequencer
            read_count = fileObject.db_entry.read_count
            aln_job_walltime = self.scaleJobWalltimeOrMemoryBasedOnInput(
                realInputVolume=read_count, \
                baseInputVolume=4000000, baseJobPropertyValue=baseAlnJobWalltime, \
                minJobPropertyValue=baseAlnJobWalltime*2/3,
                maxJobPropertyValue=baseAlnJobWalltime*5).value

            if alignment_method.command.find('aln')>=0 and sequencer.short_name!='454' and \
                (sequence_type and sequence_type.read_length_mean is not None \
                and sequence_type.read_length_mean<150):
                #short single-end read
                fileBasenamePrefix = utils.getRealPrefixSuffixOfFilenameWithVariableSuffix(
                    os.path.basename(relativePath))[0]
                outputFname = os.path.join(outputDir, '%s.sai'%fileBasenamePrefix)
                saiOutput = File(outputFname)
                alignmentJob = self.addBWAAlignmentJob(executable=bwa,
                    bwaCommand=alignment_method.command,
                    fileObject=fileObject, outputFile=saiOutput,\
                    refFastaFList=refFastaFList,
                    no_of_aln_threads=no_of_aln_threads,
                    maxMissingAlignmentFraction=maxMissingAlignmentFraction,
                    maxNoOfGaps=maxNoOfGaps,
                    additionalArguments=additionalArguments, \
                    refIndexJob=refIndexJob, parentJobLs=parentJobLs, \
                    extraDependentInputLs=None, extraOutputLs=None,
                    transferOutput=transferOutput, \
                    extraArguments=None, extraArgumentList=None,
                    job_max_memory=aln_job_max_memory,
                    key2ObjectForJob=None, walltime=aln_job_walltime)

                alignmentSamF = File('%s.sam.gz'%(os.path.join(outputDir, fileBasenamePrefix)))
                sai2samJob = Job(namespace=namespace, name=ShortSEAlignmentByBWA.name, version=version)
                sai2samJob.addArguments(refFastaFList[0], saiOutput, fastqF, alignmentSamF)
                for refFastaFile in refFastaFList:
                    sai2samJob.uses(refFastaFile, transfer=True, register=True, link=Link.INPUT)
                sai2samJob.uses(saiOutput, transfer=True, register=True, link=Link.INPUT)
                sai2samJob.uses(fastqF, transfer=True, register=True, link=Link.INPUT)
                pegaflow.setJobResourceRequirement(sai2samJob,
                    job_max_memory=samse_job_max_memory,\
                    walltime=100)
                self.addJob(sai2samJob)
                self.addJobDependency(parentJob=alignmentJob, childJob=sai2samJob)
                self.no_of_jobs += 1

            elif alignment_method.command.find('bwasw')>=0 or sequencer.short_name=='454' or \
                (sequence_type and sequence_type.read_length_mean is not None \
                and sequence_type.read_length_mean>150):
                #long single-end read
                fileBasenamePrefix = utils.getRealPrefixSuffixOfFilenameWithVariableSuffix(
                    os.path.basename(relativePath))[0]
                alignmentSamF = File('%s.sam.gz'%(os.path.join(outputDir, fileBasenamePrefix)))
                alignmentJob = Job(namespace=namespace,
                    name=LongSEAlignmentByBWA.name, version=version)

                alignmentJob.addArguments(refFastaFList[0], fastqF,
                    alignmentSamF, repr(no_of_aln_threads))
                for refFastaFile in refFastaFList:
                    alignmentJob.uses(refFastaFile, transfer=True,
                        register=True, link=Link.INPUT)
                alignmentJob.uses(fastqF, transfer=True, register=True, link=Link.INPUT)
                pegaflow.setJobResourceRequirement(alignmentJob,
                    job_max_memory=bwasw_job_max_memory, \
                    no_of_cpus=no_of_aln_threads, walltime=aln_job_walltime)
                self.addJob(alignmentJob)
                #fake a sai2samJob
                sai2samJob = alignmentJob
                self.no_of_jobs += 1

        elif len(fileObjectLs)==2:	#paired end
            fileObject = fileObjectLs[0]
            fastqF1 = fileObject.fastqF
            relativePath = fastqF1.name
            sequence_type = fileObject.db_entry.individual_sequence.sequence_type
            sequencer = fileObject.db_entry.individual_sequence.sequencer
            read_count = fileObject.db_entry.read_count
            aln_job_walltime = self.scaleJobWalltimeOrMemoryBasedOnInput(
                realInputVolume=read_count, \
                baseInputVolume=4000000, baseJobPropertyValue=baseAlnJobWalltime, \
                minJobPropertyValue=baseAlnJobWalltime*2/3, 
                maxJobPropertyValue=baseAlnJobWalltime*5).value

            #fastqF1, format, sequence_type = fileObjectLs[0][:3]
            #fastqF2, format, sequence_type = fileObjectLs[1][:3]
            fileBasenamePrefix = utils.getRealPrefixSuffixOfFilenameWithVariableSuffix(
                os.path.basename(relativePath))[0]
            # #2013.3.18 stop discarding the last two characters of filename prefix
            #fileBasenamePrefix = fileBasenamePrefix[:-2]
            # full filename for pair-1 fastq is 13135_3628_20_gerald_81LWDABXX_5_ATCACG_1_3.fastq.gz.
            #  the last 3 is split order.
            # the last 1 means pair-1. Doesn't affect filename uniqueness in either way.
            #  since the file id (13135) is unique.
            alignmentSamF = File('%s.sam.gz'%(os.path.join(outputDir, fileBasenamePrefix)))
            sai2samJob = Job(namespace=namespace, name=PEAlignmentByBWA.name,
                version=version)
            sai2samJob.addArguments(refFastaFList[0])
            pegaflow.setJobResourceRequirement(sai2samJob,
                job_max_memory=sampe_job_max_memory)
            for refFastaFile in refFastaFList:
                sai2samJob.uses(refFastaFile, transfer=True, register=True, link=Link.INPUT)
            self.addJob(sai2samJob)
            self.no_of_jobs += 1
            for fileObject in fileObjectLs:
                fastqF = fileObject.fastqF
                relativePath = fastqF.name
                sequence_type = fileObject.db_entry.individual_sequence.sequence_type
                sequencer = fileObject.db_entry.individual_sequence.sequencer

                #fastqF, format, sequence_type = fileObject[:3]

                #relativePath, format, sequence_type = fileObject[:3]
                tmp_fname_prefix = utils.getRealPrefixSuffixOfFilenameWithVariableSuffix(
                    os.path.basename(relativePath))[0]
                outputFname = os.path.join(outputDir, '%s.sai'%tmp_fname_prefix)
                saiOutput = File(outputFname)

                alignmentJob = self.addBWAAlignmentJob(executable=bwa, 
                    bwaCommand=alignment_method.command, \
                    fileObject=fileObject, outputFile=saiOutput,\
                    refFastaFList=refFastaFList, no_of_aln_threads=no_of_aln_threads, \
                    maxMissingAlignmentFraction=maxMissingAlignmentFraction,
                    maxNoOfGaps=maxNoOfGaps, \
                    additionalArguments=additionalArguments, \
                    refIndexJob=refIndexJob, parentJobLs=parentJobLs, \
                    extraDependentInputLs=None, extraOutputLs=None,
                    transferOutput=transferOutput,
                    extraArguments=None, extraArgumentList=None,
                    job_max_memory=aln_job_max_memory,
                    key2ObjectForJob=None, walltime=aln_job_walltime)

                sai2samJob.addArguments(saiOutput)
                sai2samJob.uses(saiOutput, transfer=True, register=True, link=Link.INPUT)
                self.depends(parent=alignmentJob, child=sai2samJob)

            #add a pair of fastq files to sampe in the end
            for fileObject in fileObjectLs:
                fastqF = fileObject.fastqF
                sai2samJob.addArguments(fastqF)
                sai2samJob.uses(fastqF, transfer=True, register=True, link=Link.INPUT)
            sai2samJob.addArguments(alignmentSamF)
        else:
            sys.stderr.write("addBWAAlignmentWrapJob Error: %s (!=1, !=2) file objects (%s).\n"%\
                            (len(fileObjectLs), fileObjectLs))
            raise

        sai2samJob.uses(alignmentSamF, transfer=transferOutput, register=True, link=Link.OUTPUT)
        sai2samJob.output = alignmentSamF
        sai2samJob.outputLs = [alignmentSamF]
        sai2samJob.fileBasenamePrefix = fileBasenamePrefix
        if refIndexJob:
            self.addRefIndexJobAndItsOutputAsParent(refIndexJob, childJob=sai2samJob)
        if parentJobLs:
            for parentJob in parentJobLs:
                if parentJob:
                    self.addJobDependency(parentJob=parentJob, childJob=sai2samJob)
        return sai2samJob

    def addBWAMemAlignmentJob(self, fileObjectLs=None, \
        refFastaFList=None, bwa=None, additionalArguments=None, \
        samtools=None, refIndexJob=None, parentJobLs=None,\
        alignment_method=None, outputDir=None, namespace=None, version=None,
        PEAlignmentByBWA=None, ShortSEAlignmentByBWA=None,
        LongSEAlignmentByBWA=None,
        java=None, AddOrReplaceReadGroupsJava=None,
        AddOrReplaceReadGroupsJar=None,\
        no_of_aln_threads=3, maxMissingAlignmentFraction=None,
        maxNoOfGaps=None,
        extraArgumentList=None, transferOutput=False, **keywords):
        """
        additionalArguments = additionalArguments.strip()
        #2013.06.20 strip all the spaces around it.
        2013.04.04
new alignment method from Heng Li
raw bwa command:
./bwa mem -M -a 3231_VervetRef1.0.3.fasta 12467_sequence_628BWAAXX_2_1_1.fastq.gz >/tmp/aln-pe.2.sam
in pipe2File:
    pipe2File.sh ./bwa /tmp/aln-pe.2.sam.gz mem -t 1 -M -a 3280.fasta 12457_1.fastq.gz 12457_2.fastq.gz
        """
        if namespace is None:
            namespace = self.namespace
        if version is None:
            version = self.version

        if extraArgumentList is None:
            extraArgumentList = []
        extraArgumentList.append(alignment_method.command)	#add mem first
        additionalArguments = additionalArguments.strip()
        #2013.06.20 strip all the spaces around it.
        if additionalArguments:
            extraArgumentList.append(additionalArguments)
        if no_of_aln_threads:
            extraArgumentList.append("-t %s"%no_of_aln_threads)
        #in MB, 2.5GB is enough for 4.5G gzipped fastq versus 480K contigs (total size~3G)
        aln_job_max_memory = 7000

        #2013.3.1 change the walltime by multiplying it
        baseAlnJobWalltime = int(1380*self.coreAlignmentJobWallTimeMultiplier)
            #1380 is 23 hours, because all reads are stored in chunks of 5-million-read files

        refFastaFile = refFastaFList[0]

        if len(fileObjectLs)==1 or len(fileObjectLs)==2:
            #paired end	#single end
            fileObject = fileObjectLs[0]
            fastqF = fileObject.fastqF
            relativePath = fastqF.name
            read_count = fileObject.db_entry.read_count
            aln_job_walltime = self.scaleJobWalltimeOrMemoryBasedOnInput(
                realInputVolume=read_count, \
                baseInputVolume=4000000, baseJobPropertyValue=baseAlnJobWalltime, \
                minJobPropertyValue=baseAlnJobWalltime*2/3, 
                maxJobPropertyValue=baseAlnJobWalltime*5).value

            fileBasenamePrefix = utils.getRealPrefixSuffixOfFilenameWithVariableSuffix(
                os.path.basename(relativePath))[0]
            alignmentSamF = File('%s.sam.gz'%(os.path.join(outputDir, fileBasenamePrefix)))

            fastqFileList = [fileObject.fastqF for fileObject in fileObjectLs]
            extraArgumentList.extend(["-a -M", refFastaFile] + fastqFileList)
            extraDependentInputLs=fastqFileList + refFastaFList

            alignmentJob = self.addPipeCommandOutput2FileJob(executable=self.BWA_Mem, \
                commandFile=self.bwaExecutableFile, \
                outputFile=alignmentSamF, \
                parentJobLs=parentJobLs, extraDependentInputLs=extraDependentInputLs, \
                extraOutputLs=None, transferOutput=transferOutput, \
                extraArguments=None, extraArgumentList=extraArgumentList, \
                sshDBTunnel=None,\
                job_max_memory=aln_job_max_memory, 
                walltime=aln_job_walltime, no_of_cpus=no_of_aln_threads, \
                **keywords)
        else:
            sys.stderr.write("addBWAMemAlignmentJob Error: %s (!=1, !=2) file objects (%s).\n"%\
                            (len(fileObjectLs), fileObjectLs))
            raise


        alignmentJob.fileBasenamePrefix = fileBasenamePrefix
        if refIndexJob:
            self.addRefIndexJobAndItsOutputAsParent(refIndexJob, childJob=alignmentJob)
        return alignmentJob

    def addAlignmentJob(self, fileObjectLs=None, individual_alignment=None, \
        refFastaFList=None, bwa=None, additionalArguments=None, \
        samtools=None, refIndexJob=None, parentJobLs=None, \
        alignment_method=None, outputDir=None, namespace=None, version=None,\
        PEAlignmentByBWA=None, ShortSEAlignmentByBWA=None, LongSEAlignmentByBWA=None, \
        java=None, SortSamFilesJava=None, SortSamJar=None,\
        AddOrReplaceReadGroupsJava=None, AddOrReplaceReadGroupsJar=None,\
        no_of_aln_threads=3, stampy=None, tmpDir=None,\
        maxMissingAlignmentFraction=None, maxNoOfGaps=None, addBamIndexJob=False,\
        transferOutput=False, **keywords):
        """
        Examples:

        #2012.9.19 individual_alignment is passed as None so that
        #  ReadGroup addition job is not added in addAlignmentJob()
        alignmentJob, alignmentOutput = self.addAlignmentJob(fileObjectLs=newFileObjLs, \
            individual_alignment=None, \
            data_dir=data_dir, refFastaFList=refFastaFList, bwa=bwa, \
            additionalArguments=additionalArguments, samtools=samtools, \
            refIndexJob=refIndexJob, parentJobLs=[refIndexJob, mkdirJob], \
            alignment_method=alignment_method, \
            outputDir=tmpOutputDir, namespace=namespace, version=version,\
            PEAlignmentByBWA=PEAlignmentByBWA, ShortSEAlignmentByBWA=ShortSEAlignmentByBWA, \
            LongSEAlignmentByBWA=LongSEAlignmentByBWA,\
            java=java, SortSamFilesJava=SortSamFilesJava, SortSamJar=SortSamJar,\
            AddOrReplaceReadGroupsJava=AddOrReplaceReadGroupsJava, 
            AddOrReplaceReadGroupsJar=AddOrReplaceReadGroupsJar,\
            no_of_aln_threads=no_of_aln_threads, stampy=stampy)

        2013.04.27 used to be 'mem', now is 'bwamem' in db
        2013.04.04 new alignment method (bwa-mem) from Heng Li
        2012.10.18 add argument addBamIndexJob,
        2012.10.10
        not necessary to add refIndexJob into parentJobLs, because it'll be added as parent job.
        added argument maxMissingAlignmentFraction, maxNoOfGaps for bwa
        each fileObject in fileObjectLs should have 2 attributes
            fastqF: a registered pegasus file
            db_entry: an individual_sequence_file db_entry or equivalent object. 2 attributes:
                quality_score_format: 'Standard', 'Illumina'
                individual_sequence: 2 attributes:
                    sequencer:  'GA' , '454'
                    sequence_type: 'PE', 'SR', 'genome'	#(not really used as 2012.10.10)
        alignment_method is also an object (supposedly from vervet db), two attributes.
            short_name (stampy, or bwa-short-read, or bwa-long-read)
            command (stampy.py , aln, bwasw)
        2012.9.19
            modify it substantially
        2012.4.5
            choose which alignment program to use based on alignment_method.short_name
        2012.2.23
            split the BWA alignment part to addBWAAlignmentWrapJob()
        2011-9-13
            add argument java & SortSamJar
        2011-9-9
            two steps:
                1. aln doesn't use pipe, outputs to sai files.
                2. sampe/samse, convert, sort => connected through pipe
        """

        if alignment_method.short_name=='stampy':
            alignmentJob = self.addStampyAlignmentJob(
                fileObjectLs=fileObjectLs,
                refFastaFList=refFastaFList, bwa=bwa, \
                additionalArguments=additionalArguments, samtools=samtools, \
                refIndexJob=refIndexJob, parentJobLs=parentJobLs, \
                alignment_method=alignment_method, \
                outputDir=outputDir, namespace=namespace, version=version,\
                PEAlignmentByBWA=PEAlignmentByBWA,
                ShortSEAlignmentByBWA=ShortSEAlignmentByBWA,
                LongSEAlignmentByBWA=LongSEAlignmentByBWA,\
                java=java,
                AddOrReplaceReadGroupsJava=AddOrReplaceReadGroupsJava,
                AddOrReplaceReadGroupsJar=AddOrReplaceReadGroupsJar,
                no_of_aln_threads=no_of_aln_threads, stampy=stampy,\
                transferOutput=transferOutput)
        elif alignment_method.short_name.find('bwamem')==0:
            #2013.04.27 used to be 'mem', now is 'bwamem'
            alignmentJob = self.addBWAMemAlignmentJob(
                fileObjectLs=fileObjectLs,
                refFastaFList=refFastaFList, bwa=bwa, \
                additionalArguments=additionalArguments, samtools=samtools, \
                refIndexJob=refIndexJob, parentJobLs=parentJobLs, \
                alignment_method=alignment_method, \
                outputDir=outputDir, namespace=namespace, version=version,\
                PEAlignmentByBWA=PEAlignmentByBWA,
                ShortSEAlignmentByBWA=ShortSEAlignmentByBWA,
                LongSEAlignmentByBWA=LongSEAlignmentByBWA,\
                java=java,
                AddOrReplaceReadGroupsJava=AddOrReplaceReadGroupsJava, 
                AddOrReplaceReadGroupsJar=AddOrReplaceReadGroupsJar,\
                no_of_aln_threads=no_of_aln_threads, 
                maxMissingAlignmentFraction=maxMissingAlignmentFraction, 
                maxNoOfGaps=maxNoOfGaps, transferOutput=transferOutput)
        elif alignment_method.short_name.find('bwa')==0:
            alignmentJob = self.addBWAAlignmentWrapJob(
                fileObjectLs=fileObjectLs, \
                refFastaFList=refFastaFList, bwa=bwa, \
                additionalArguments=additionalArguments, samtools=samtools, \
                refIndexJob=refIndexJob, parentJobLs=parentJobLs, \
                alignment_method=alignment_method, \
                outputDir=outputDir, namespace=namespace, version=version,\
                PEAlignmentByBWA=PEAlignmentByBWA,
                ShortSEAlignmentByBWA=ShortSEAlignmentByBWA,
                LongSEAlignmentByBWA=LongSEAlignmentByBWA,\
                java=java,
                AddOrReplaceReadGroupsJava=AddOrReplaceReadGroupsJava, 
                AddOrReplaceReadGroupsJar=AddOrReplaceReadGroupsJar,\
                no_of_aln_threads=no_of_aln_threads,  
                maxMissingAlignmentFraction=maxMissingAlignmentFraction, 
                maxNoOfGaps=maxNoOfGaps, transferOutput=transferOutput)
        else:
            sys.stderr.write("Alignment method %s is not supported.\n"%(alignment_method.short_name))
            sys.exit(3)
        fileBasenamePrefix = alignmentJob.fileBasenamePrefix

        ## convert sam into bam
        bamOutputF = File(os.path.join(outputDir, "%s.bam"%(fileBasenamePrefix)))
        #2013.3.18
        sam_convert_job = self.addSAM2BAMJob(inputFile=alignmentJob.output,
            outputFile=bamOutputF,
            executable=self.SAM2BAMJava, SamFormatConverterJar=self.PicardJar,
            parentJobLs=[alignmentJob], extraDependentInputLs=None, \
            extraArguments=None, job_max_memory = 3000, walltime=80,
            transferOutput=transferOutput)
        """
        #20130318 samtools is buggy, constantly reporting inconsistency 
        # between sequence length and CIGAR length
        sam_convert_job = self.addGenericJob(executable=samtools, inputFile=None,\
            outputFile=bamOutputF, outputArgumentOption="view -bSh -o", \
            parentJobLs=[alignmentJob],
            extraDependentInputLs=[alignmentJob.output], \
            extraOutputLs=[],\
            transferOutput=transferOutput, \
            extraArgumentList=[alignmentJob.output], \
            job_max_memory=2000)
        """
        if individual_alignment:
            #if not given then , no read-group addition job
            #2012.9.19 add a AddReadGroup job
            outputRGBAM = File(os.path.join(outputDir, "%s.RG.bam"%(fileBasenamePrefix)))
            addRGJob = self.addReadGroupInsertionJob(
                individual_alignment=individual_alignment, \
                inputBamFile=sam_convert_job.output, \
                outputBamFile=outputRGBAM,\
                AddOrReplaceReadGroupsJava=AddOrReplaceReadGroupsJava, \
                AddOrReplaceReadGroupsJar=AddOrReplaceReadGroupsJar,\
                parentJobLs=[sam_convert_job], extraDependentInputLs=None, \
                extraArguments=None, job_max_memory = 2500, walltime=80, \
                transferOutput=transferOutput)
            sortAlnParentJob = addRGJob
        else:
            sortAlnParentJob = sam_convert_job

        """
        # 2010-2-4
            sort it so that it could be used for merge
        """
        bam_output_fname_prefix = '%s.sorted'%(os.path.join(outputDir, fileBasenamePrefix))
        sortBamF = File('%s.bam'%(bam_output_fname_prefix))
        sortAlignmentJob = self.addSortAlignmentJob(
            inputBamFile=sortAlnParentJob.output, \
            outputBamFile=sortBamF,\
            SortSamFilesJava=SortSamFilesJava, SortSamJar=SortSamJar, tmpDir=tmpDir,\
            parentJobLs=[sortAlnParentJob], extraDependentInputLs=None, \
            extraArguments=None, job_max_memory = 2500, walltime=80, \
            transferOutput=transferOutput, needBAMIndexJob=addBamIndexJob)

        if addBamIndexJob:
            #bamIndexJob.parentJobLs[0] is sortAlignmentJob.
            returnJob = sortAlignmentJob.bamIndexJob
        else:
            returnJob = sortAlignmentJob
        return returnJob, returnJob.output

    def addSAM2BAMJob(self, inputFile=None, outputFile=None,\
        executable=None, SamFormatConverterJar=None,\
        parentJobLs=None, extraDependentInputLs=None, \
        extraArguments=None, job_max_memory = 2500, transferOutput=False, 
        walltime=120, **keywords):
        """
        2013.3.18
        """
        if executable is None:
            executable = self.SAM2BAMJava
        memRequirementData = self.getJVMMemRequirment(
            job_max_memory=job_max_memory, minMemory=2000)
        job_max_memory = memRequirementData.memRequirement
        javaMemRequirement = memRequirementData.memRequirementInStr

        extraArgumentList = [memRequirementData.memRequirementInStr, '-jar', 
            SamFormatConverterJar, "SamFormatConverter", \
            "I=", inputFile, "O=", outputFile, "VALIDATION_STRINGENCY=LENIENT"]
        if extraArguments:
            extraArgumentList.append(extraArguments)
        if extraDependentInputLs is None:
            extraDependentInputLs=[]
        extraDependentInputLs.extend([inputFile, SamFormatConverterJar])

        job= self.addGenericJob(executable=executable, inputFile=None,\
            parentJobLs=parentJobLs, extraDependentInputLs=extraDependentInputLs, \
            extraOutputLs=[outputFile],\
            transferOutput=transferOutput, \
            extraArgumentList=extraArgumentList, \
            job_max_memory=memRequirementData.memRequirement, walltime=walltime,\
            **keywords)
        return job

    def addAlignmentMergeBySAMtoolsJob(self, AlignmentJobAndOutputLs=[], 
        outputBamFile=None,samtools=None,\
        java=None, MergeSamFilesJava=None, \
        BuildBamIndexFilesJava=None, BuildBamIndexJar=None, \
        mv=None, transferOutput=False, job_max_memory=2500, **keywords):
        """
        2012-3.22
            samtools merge version of addAlignmentMergeJob(),
                which uses picard's MergeSamFiles.jar
            **untested**
        """
        javaMaxMemory=2500
        if len(AlignmentJobAndOutputLs)>1:
            merge_sam_job = Job(namespace=self.namespace, name=samtools.name,
                version=self.version)
            merge_sam_job.addArguments('-f', outputBamFile)	
            #'overwrite the output BAM if exist'
            merge_sam_job.uses(outputBamFile, transfer=transferOutput,
                register=True, link=Link.OUTPUT)
            pegaflow.setJobResourceRequirement(merge_sam_job,
                job_max_memory=job_max_memory)
            self.addJob(merge_sam_job)
            for AlignmentJobAndOutput in AlignmentJobAndOutputLs:
                alignmentJob, alignmentOutput = AlignmentJobAndOutput[:2]
                merge_sam_job.addArguments(alignmentOutput)
                merge_sam_job.uses(alignmentOutput, transfer=True,
                    register=True, link=Link.INPUT)
                self.depends(parent=alignmentJob, child=merge_sam_job)
        else:	#one input file, no samtools merge. use "mv" to rename it instead
            alignmentJob, alignmentOutput = AlignmentJobAndOutputLs[0][:2]
            merge_sam_job = Job(namespace=self.namespace, name=mv.name,
                version=self.version)
            merge_sam_job.addArguments(alignmentOutput, outputBamFile)
            self.depends(parent=alignmentJob, child=merge_sam_job)
            merge_sam_job.uses(alignmentOutput, transfer=True, register=True,
                link=Link.INPUT)
            merge_sam_job.uses(outputBamFile, transfer=transferOutput,
                register=True, link=Link.OUTPUT)
            self.addJob(merge_sam_job)

        # add the index job on the merged bam file
        bamIndexJob = self.addBAMIndexJob(
            BuildBamIndexFilesJava=BuildBamIndexFilesJava,
            BuildBamIndexJar=BuildBamIndexJar, \
            inputBamF=outputBamFile, parentJobLs=[merge_sam_job], \
            transferOutput=transferOutput, javaMaxMemory=javaMaxMemory,
            walltime=180)
        return merge_sam_job, bamIndexJob

    def addMarkDupJob(self, parentJobLs=[], inputBamF=None, inputBaiF=None,
        outputBamFile=None,\
        MarkDuplicatesJava=None, MarkDuplicatesJar=None, tmpDir=None,\
        BuildBamIndexFilesJava=None, BuildBamIndexJar=None, \
        namespace=None, version=None, transferOutput=True, \
        job_max_memory=4000, walltime=600, no_of_cpus=1):
        """
        2012.3.21
            improve it
            set no_of_cpus=1 (was 2) to avoid thread problem in some linux kernels.
        #2011-11-10 duplicate-marking job
        """
        MarkDupJob = Job(namespace=getattr('namespace', namespace),
            name=MarkDuplicatesJava.name, \
            version=getattr('version', version))
        bamFnamePrefix = os.path.splitext(outputBamFile.name)[0]
        MarkDupOutputF = outputBamFile
        MarkDupOutputMetricF = File('%s.metric'%(bamFnamePrefix))

        memRequirementData = self.getJVMMemRequirment(
            job_max_memory=job_max_memory, minMemory=8000, \
            permSizeFraction=0.2)
        job_max_memory = memRequirementData.memRequirement
        javaMemRequirement = memRequirementData.memRequirementInStr

        MarkDupJob.addArguments(javaMemRequirement, '-jar', MarkDuplicatesJar, 
            "MarkDuplicates", "MAX_FILE_HANDLES=500",\
            "VALIDATION_STRINGENCY=LENIENT", "ASSUME_SORTED=true", "INPUT=", inputBamF, \
            'OUTPUT=', MarkDupOutputF, "M=", MarkDupOutputMetricF, "MAX_RECORDS_IN_RAM=500000",\
            "TMP_DIR=%s"%tmpDir)
        self.addJobUse(MarkDupJob, file=MarkDuplicatesJar, transfer=True, register=True, link=Link.INPUT)
        MarkDupJob.uses(inputBaiF, transfer=True, register=True, link=Link.INPUT)
        MarkDupJob.uses(inputBamF, transfer=True, register=True, link=Link.INPUT)
        MarkDupJob.output = MarkDupOutputF
        MarkDupJob.MarkDupOutputMetricF = MarkDupOutputMetricF
        MarkDupJob.uses(MarkDupOutputF, transfer=transferOutput, register=True, link=Link.OUTPUT)
        MarkDupJob.uses(MarkDupOutputMetricF, transfer=transferOutput, register=True, link=Link.OUTPUT)
        #pass	#don't register the files so leave them there
        self.addJob(MarkDupJob)
        pegaflow.setJobResourceRequirement(MarkDupJob, job_max_memory=job_max_memory, 
            no_of_cpus=no_of_cpus, walltime=walltime)
        for parentJob in parentJobLs:
            if parentJob:
                self.depends(parent=parentJob, child=MarkDupJob)


        # add the index job on the bam file
        bamIndexJob = self.addBAMIndexJob(
            BuildBamIndexFilesJava=BuildBamIndexFilesJava,
            BuildBamIndexJar=BuildBamIndexJar,
            inputBamF=MarkDupOutputF,\
            parentJobLs=[MarkDupJob], namespace=namespace, version=version,\
            transferOutput=transferOutput, walltime=max(180, int(walltime/3)))
        return MarkDupJob, bamIndexJob

    def addSAMtoolsCalmdJob(self, samtoolsCalmd=None, inputBamF=None, \
        refFastaFList=None, outputBamF=None, \
        parentJob=None, \
        BuildBamIndexFilesJava=None, BuildBamIndexJar=None, \
        namespace=None, version=None, transferOutput=True,\
        **keywords):
        """
        2011-11-20
            run "samtools calmd" on the input bam and index the output bam
        """
        job = Job(namespace=namespace, name=BuildBamIndexFilesJava.name, version=version)
        job.addArguments(inputBamF, refFastaFList[0], outputBamF)
        pegaflow.setJobResourceRequirement(job, job_max_memory=1000)
        job.uses(inputBamF, transfer=True, register=True, link=Link.INPUT)
        for refFastaFile in refFastaFList:
            job.uses(refFastaFile, transfer=True, register=True, link=Link.INPUT)
        if transferOutput:
            job.uses(outputBamF, transfer=True, register=True, link=Link.OUTPUT)
        else:
            pass	#don't register the files so leave them there
        self.addJob(job)
        if parentJob:
            self.depends(parent=parentJob, child=job)

        # add the index job on the bam
        return self.addBAMIndexJob(BuildBamIndexFilesJava=BuildBamIndexFilesJava, 
            BuildBamIndexJar=BuildBamIndexJar, \
            inputBamF=outputBamF,\
            parentJobLs=[job], namespace=namespace, version=version,\
            transferOutput=transferOutput)

    def addLocalRealignmentSubWorkflow(self, chr2IntervalDataLs=None,
        registerReferenceData=None, alignmentData=None,\
        inputBamF=None, outputBamF=None,
        parentJobLs=None, \
        outputDirPrefix='localRealignment', transferOutput=False,\
        **keywords):
        """
        2013.03.31
            create a sub-workflow to do local re-alignment for one inputBAM
            for each 2million bp interval
                split a small bam out of input bam (add 1kb on both ends of interval)
                index the small bam
                run RealignerTargetCreator
                run IndelRealigner
                (extract the exact interval out of bam if 1kb is added to 
                    both ends of input interval)
            merge all intervals back into one bam
            index the merged bam

        """
        chrIDSet = set(chr2IntervalDataLs.keys())
        chr2VCFJobData = {}

        if self.report:
            sys.stderr.write("Adding local indel-realignment jobs for %s chromosomes/contigs ..."%\
                            (len(chrIDSet)))
        refFastaFList = registerReferenceData.refFastaFList
        topOutputDir = "%sMap"%(outputDirPrefix)
        topOutputDirJob = self.addMkDirJob(outputDir=topOutputDir)
        reduceOutputDir = "%sReduce"%(outputDirPrefix)
        reduceOutputDirJob = self.addMkDirJob(outputDir=reduceOutputDir)
        returnData = PassingData()
        returnData.jobDataLs = []

        #2012.9.22 AlignmentJobAndOutputLs is a relic.
        #	but it's similar to mapEachIntervalDataLs but designed for addAlignmentMergeJob(),
        #	so AlignmentJobAndOutputLs gets re-set for every alignment.
        # 	mapEachAlignmentDataLs is never reset.
        #	mapEachChromosomeDataLs is reset right after a new alignment is chosen.
        #	mapEachIntervalDataLs is reset right after each chromosome is chosen.
        #	all reduce dataLs never gets reset.
        passingData = PassingData(AlignmentJobAndOutputLs=[], \
            alignmentData = alignmentData,\
            bamFnamePrefix=None, \
            outputDirPrefix=outputDirPrefix, \
            topOutputDirJob=topOutputDirJob,\
            reduceOutputDirJob = reduceOutputDirJob,\
            refFastaFList=refFastaFList, \
            registerReferenceData= registerReferenceData,\
            refFastaF=refFastaFList[0],\
            fastaDictJob = None,\
            refFastaDictF = None,\
            fastaIndexJob = None,\
            refFastaIndexF = None,\
            chrIDSet=chrIDSet,\
            chr2IntervalDataLs=chr2IntervalDataLs,\
            mapEachAlignmentData = None,\
            mapEachChromosomeData=None, \
            mapEachIntervalData=None,\
            reduceBeforeEachAlignmentData = None, \
            reduceAfterEachAlignmentData=None,\
            reduceAfterEachChromosomeData=None,\

            mapEachAlignmentDataLs = [],\
            mapEachChromosomeDataLs=[], \
            mapEachIntervalDataLs=[],\
            reduceBeforeEachAlignmentDataLs = [], \
            reduceAfterEachAlignmentDataLs=[],\
            reduceAfterEachChromosomeDataLs=[],\

            gzipReduceAfterEachChromosomeFolderJob=None,\
            gzipReduceBeforeEachAlignmentFolderJob = None,\
            gzipReduceAfterEachAlignmentFolderJob = None,\
            gzipPreReduceFolderJob = None,\
            gzipReduceFolderJob=None,\
            )
        preReduceReturnData = self.preReduce(passingData=passingData,
            transferOutput=False, **keywords)
        passingData.preReduceReturnData = preReduceReturnData
        alignment = alignmentData.alignment
        parentJobLs = alignmentData.jobLs
        bamF = alignmentData.bamF
        baiF = alignmentData.baiF

        bamFnamePrefix = alignment.getReadGroup()

        passingData.AlignmentJobAndOutputLs = []
        passingData.bamFnamePrefix = bamFnamePrefix
        passingData.individual_alignment = alignment
        passingData.alignmentData = alignmentData


        self.mapReduceOneAlignment(alignmentData=alignmentData,
            passingData=passingData, \
            chrIDSet=chrIDSet, chr2IntervalDataLs=chr2IntervalDataLs,
            chr2VCFJobData=chr2VCFJobData, \
            outputDirPrefix=outputDirPrefix, transferOutput=transferOutput,
            skipChromosomeIfVCFMissing=False)

        """
        2013.03.31 not needed
        reduceAfterEachAlignmentData = self.reduceAfterEachAlignment(\
            mapEachAlignmentData=mapEachAlignmentData,\
            mapEachChromosomeDataLs=passingData.mapEachChromosomeDataLs,\
            reduceAfterEachChromosomeDataLs=passingData.reduceAfterEachChromosomeDataLs,\
            passingData=passingData, \
            transferOutput=False, data_dir=data_dir, **keywords)
        passingData.reduceAfterEachAlignmentData = reduceAfterEachAlignmentData
        passingData.reduceAfterEachAlignmentDataLs.append(reduceAfterEachAlignmentData)

        """
        AlignmentJobAndOutputLs = passingData.AlignmentJobAndOutputLs
        mergedBamFile = File(os.path.join(topOutputDirJob.output, \
            '%s_realigned.bam'%(bamFnamePrefix)))
        alignmentMergeJob, bamIndexJob = self.addAlignmentMergeJob(
            AlignmentJobAndOutputLs=AlignmentJobAndOutputLs,
            outputBamFile=mergedBamFile, \
            samtools=self.samtools, java=self.java, \
            MergeSamFilesJava=self.MergeSamFilesJava,
            MergeSamFilesJar=self.MergeSamFilesJar,
            BuildBamIndexFilesJava=self.IndexMergedBamIndexJava,
            BuildBamIndexJar=self.BuildBamIndexJar, \
            mv=self.mv, parentJobLs=[topOutputDirJob], \
            transferOutput=False)

        if self.report:
            sys.stderr.write("%s jobs.\n"%(self.no_of_jobs))
        return alignmentMergeJob, bamIndexJob

    def mapEachInterval(self, alignmentData=None, intervalData=None,\
        VCFJobData=None, passingData=None, reduceBeforeEachAlignmentData=None,\
        mapEachChromosomeData=None, transferOutput=False, **keywords):
        """
        2013.03.31 copied from AlignmentReadBaseQualityRecalibrationWorkflow.py
            called up by mapReduceOneAlignment(),
             which is in turn called up by addLocalRealignmentJobs
        """
        returnData = AlignmentReadBaseQualityRecalibrationWorkflow.mapEachInterval(
            self,
            alignmentData=alignmentData, \
            intervalData=intervalData,\
            VCFJobData=None, passingData=passingData,
            reduceBeforeEachAlignmentData=reduceBeforeEachAlignmentData,\
            mapEachChromosomeData=mapEachChromosomeData,
            transferOutput=transferOutput, \
            **keywords)
        return returnData
    
    def addAllAlignmentJobs(self, db_main=None,
        individualSequenceID2FilePairLs=None, \
        data_dir=None, \
        isqLs=None,\
        refSequence=None, registerReferenceData=None, \
        chr2IntervalDataLs=None,\
        bwa=None, additionalArguments=None, samtools=None,
        mkdirWrap=None, mv=None,\
        java=None, MergeSamFilesJava=None, MergeSamFilesJar=None, \
        MarkDuplicatesJava=None, MarkDuplicatesJar=None, tmpDir='/tmp',\
        BuildBamIndexFilesJava=None, BuildBamIndexJar=None, \
        SortSamFilesJava=None, SortSamJar=None, \
        AddOrReplaceReadGroupsJava=None, AddOrReplaceReadGroupsJar=None,\
        alignment_method_name='bwaShortRead', alignment_format='bam',\
        namespace='workflow', version='1.0', transferOutput=False,\
        PEAlignmentByBWA=None, ShortSEAlignmentByBWA=None,
        LongSEAlignmentByBWA=None, \
        no_of_aln_threads=3, stampy=None, skipDoneAlignment=False, \
        alignmentPerLibrary=False, outputDirPrefix="", **keywords):
        """
        Scale MarkDuplicates and MergeSam jobs memory & walltime to
            the sequence coverage
        2012.9.19
            isq_id2LibrarySplitOrder2FileLs is replaced by isqLs.
            add argument alignmentPerLibrary
        2012.4.20
            bugfix, pass alignment_method.short_name instead of
                alignment_method_name to db_main.getAlignment()
                because alignment_method might be changed according to sequencer
                    regardless of alignment_method_name.
        2012.4.12
            add argument skipDoneAlignment
        2012.4.5
            fetch the alignment_method directly based on the
                alignment_method_name, except for 454 sequences.
                officially merges "bwa-short-read-SR" (single-read) into "bwa-short-read"
        2012.2.24
            pass data_dir to db_main.getAlignment()
            add stampy part
        2011-9-15
            adjust alignment_method_name according to individual_sequence.sequencer
                and individual_sequence.sequence_type
            only when this is not possible, value of argument alignment_method_name is used.
        2011-9-14
            give different names to different java jobs according to jars
        2011-8-28
        """
        sys.stderr.write("Adding alignment jobs for %s individual sequences ..."%(len(isqLs)))

        no_of_alignment_jobs = 0
        no_of_merging_jobs = 0

        alignmentFolder = "%sAlignment"%(outputDirPrefix)
        alignmentFolderJob = self.addMkDirJob(outputDir=alignmentFolder)

        oneLibraryAlignmentFolder = "%sOneLibAlignment"%(outputDirPrefix)
        oneLibraryAlignmentFolderJob = self.addMkDirJob(
            outputDir=oneLibraryAlignmentFolder)

        refFastaFList = registerReferenceData.refFastaFList
        refIndexJob = None
        if self.needRefIndexJob or registerReferenceData.needBWARefIndexJob or \
            registerReferenceData.needStampyRefIndexJob:
            if self.alignment_method_name.find('bwa')>=0 and \
                registerReferenceData.needBWARefIndexJob:
                refIndexJob = self.addBWAReferenceIndexJob(
                    refFastaFList=refFastaFList,
                    refSequenceBaseCount=refSequence.base_count, bwa=self.bwa)
            elif self.alignment_method_name.find('stampy')>=0 and \
                registerReferenceData.needStampyRefIndexJob:
                refIndexJob = self.addStampyGenomeIndexHashJob(
                    executable=self.stampy, refFastaFList=refFastaFList, \
                        parentJobLs=None, job_max_memory=3000, walltime = 200, \
                        extraDependentInputLs=None, \
                        transferOutput=True)

        for individual_sequence  in isqLs:
            if individual_sequence is not None and \
                    individual_sequence.format=='fastq':
                library2Data = individual_sequence.library2Data
                AlignmentJobAndOutputLs = []
                # get or add alignment method
                if individual_sequence.sequencer.short_name=='454' or \
                        (individual_sequence.sequence_type and \
                        individual_sequence.sequence_type.read_length_mean is not None \
                        and individual_sequence.sequence_type.read_length_mean>150):
                    alignment_method = db_main.getAlignmentMethod("bwaLongRead")
                else:
                    alignment_method = db_main.getAlignmentMethod(alignment_method_name)
                """
                #2012.4.5 have all this commented out
                elif individual_sequence.sequencer.short_name=='GA':
                    if individual_sequence.sequence_type=='SR':	#single-end
                        alignment_method = db_main.getAlignmentMethod("bwa-short-read-SR")
                    else:	#default is PE
                        alignment_method = db_main.getAlignmentMethod("bwa-short-read")
                """
                #alignment for the whole individual_sequence
                individual_alignment = db_main.getAlignment(
                    individual_code=individual_sequence.individual.code, \
                    individual_sequence_id=individual_sequence.id,\
                    path_to_original_alignment=None,
                    sequencer_id=individual_sequence.sequencer_id, \
                    sequence_type_id=individual_sequence.sequence_type_id,
                    sequence_format=individual_sequence.format, \
                    ref_individual_sequence_id=refSequence.id, \
                    alignment_method_name=alignment_method.short_name,
                    alignment_format=alignment_format,\
                    individual_sequence_filtered=individual_sequence.filtered,
                    read_group_added=1,
                    data_dir=data_dir, local_realigned=self.local_realigned)
                skipIndividualAlignment = False
                if individual_alignment.path:
                    if skipDoneAlignment and self.isThisAlignmentComplete(
                        individual_alignment=individual_alignment, data_dir=data_dir):
                        skipIndividualAlignment = True
                #2012.3.29 this folder will store the alignment output by the alignment jbos
                tmpOutputDir = os.path.basename(individual_sequence.path)
                # add a mkdir job
                mkdirJob = None
                for library, pdata in library2Data.items():
                    minIsqFileRawID = min(pdata.isqFileRawID2Index.keys())
                    splitOrder2Index = pdata.splitOrder2Index
                    fileObjectPairLs = pdata.fileObjectPairLs

                    oneLibraryAlignmentJobAndOutputLs = []
                    splitOrderLs = splitOrder2Index.keys()
                    splitOrderLs.sort()
                    oneLibraryCumulativeBaseCount = 0
                    if alignmentPerLibrary:
                        #alignment for this library of the individual_sequence
                        oneLibraryAlignmentEntry = db_main.getAlignment(
                            individual_code=individual_sequence.individual.code, \
                            individual_sequence_id=individual_sequence.id,\
                            path_to_original_alignment=None,
                            sequencer_id=individual_sequence.sequencer_id, \
                            sequence_type_id=individual_sequence.sequence_type_id,
                            sequence_format=individual_sequence.format, \
                            ref_individual_sequence_id=refSequence.id, \
                            alignment_method_name=alignment_method.short_name,
                            alignment_format=alignment_format,\
                            individual_sequence_filtered=individual_sequence.filtered,
                            read_group_added=1,
                            data_dir=data_dir,
                            individual_sequence_file_raw_id=minIsqFileRawID,\
                            local_realigned=self.local_realigned)
                        skipLibraryAlignment = False
                        if oneLibraryAlignmentEntry.path:
                            if skipDoneAlignment and self.isThisAlignmentComplete(
                                individual_alignment=individual_alignment, data_dir=data_dir):
                                # 2013.3.22
                                # file_size is updated in the last of AddAlignmentFile2DB.py.
                                # if it fails in the middle of copying, file_size would be None.
                                skipLibraryAlignment = True
                    else:
                        skipLibraryAlignment = True
                    if skipIndividualAlignment and skipLibraryAlignment:
                        #2012.9.19 if both skip flags are true, then yes
                        continue
                    for splitOrder in splitOrderLs:
                        splitOrderIndex = splitOrder2Index[splitOrder]
                        fileObjectLs = fileObjectPairLs[splitOrderIndex]
                        #for fileObject in fileObjectLs:
                        #	print fileObject,"222222222222"
                        #	print "111111111111111111111111111111111111111111111111111111"
                        #	exit(2)
                        oneLibraryCumulativeBaseCount += sum([fileObject.db_entry.base_count \
                            for fileObject in fileObjectLs])
                        if mkdirJob is None:	#now it's time to add the mkdirJob
                            # add a mkdir job
                            mkdirJob = self.addMkDirJob(outputDir=tmpOutputDir)
                        newFileObjLs = self.registerISQFileObjLsToWorkflow(
                            fileObjectLs=fileObjectLs)
                        #2012.9.19 individual_alignment is passed as None so
                        #  that ReadGroup addition job is not added in addAlignmentJob()
                        alignmentJob, alignmentOutput = self.addAlignmentJob(
                            fileObjectLs=newFileObjLs, \
                            individual_alignment=None, \
                            data_dir=data_dir, refFastaFList=refFastaFList,
                            bwa=bwa,
                            additionalArguments=additionalArguments,
                            samtools=samtools,
                            refIndexJob=refIndexJob,
                            parentJobLs=[refIndexJob, mkdirJob],
                            alignment_method=alignment_method, \
                            outputDir=tmpOutputDir, namespace=namespace,
                            version=version,
                            PEAlignmentByBWA=PEAlignmentByBWA,
                            ShortSEAlignmentByBWA=ShortSEAlignmentByBWA, \
                            LongSEAlignmentByBWA=LongSEAlignmentByBWA,\
                            java=java, SortSamFilesJava=SortSamFilesJava,
                            SortSamJar=SortSamJar,
                            AddOrReplaceReadGroupsJava=AddOrReplaceReadGroupsJava,
                            AddOrReplaceReadGroupsJar=AddOrReplaceReadGroupsJar,
                            no_of_aln_threads=no_of_aln_threads, stampy=stampy,
                            tmpDir=tmpDir)
                        no_of_alignment_jobs += 1

                        fileBasenamePrefix = utils.getRealPrefixSuffixOfFilenameWithVariableSuffix(
                            os.path.basename(alignmentOutput.name))[0]
                        if not skipIndividualAlignment:
                            #2012.9.19 add a AddReadGroup job
                            outputRGBAM = File(os.path.join(tmpOutputDir, \
                                "%s.isq_RG.bam"%(fileBasenamePrefix)))
                            addRGJob = self.addReadGroupInsertionJob(
                                individual_alignment=individual_alignment,
                                inputBamFile=alignmentJob.output, \
                                outputBamFile=outputRGBAM,\
                                AddOrReplaceReadGroupsJava=AddOrReplaceReadGroupsJava,
                                AddOrReplaceReadGroupsJar=AddOrReplaceReadGroupsJar,\
                                parentJobLs=[alignmentJob, mkdirJob],
                                extraDependentInputLs=None,
                                extraArguments=None, job_max_memory = 2500,
                                transferOutput=False)
                            AlignmentJobAndOutputLs.append(
                                PassingData(jobLs=[addRGJob], file=addRGJob.output))
                        if not skipLibraryAlignment:
                            #2012.9.19 add a AddReadGroup job for the library bam file
                            outputRGBAM = File(os.path.join(tmpOutputDir, \
                                "%s.isq_library_%s_RG.bam"%(fileBasenamePrefix, library)))
                            addRGJob = self.addReadGroupInsertionJob(
                                individual_alignment=oneLibraryAlignmentEntry, \
                                inputBamFile=alignmentJob.output, \
                                outputBamFile=outputRGBAM,\
                                AddOrReplaceReadGroupsJava=AddOrReplaceReadGroupsJava,
                                AddOrReplaceReadGroupsJar=AddOrReplaceReadGroupsJar,\
                                parentJobLs=[alignmentJob, mkdirJob],
                                extraDependentInputLs=None,
                                extraArguments=None, job_max_memory = 2500,
                                transferOutput=False)
                            oneLibraryAlignmentJobAndOutputLs.append(
                                PassingData(parentJobLs=[addRGJob], file=addRGJob.output))
                    if alignmentPerLibrary and not skipLibraryAlignment and oneLibraryAlignmentJobAndOutputLs:
                        baseCoverage = 4*3000000000	#baseline
                        minMergeAlignmentWalltime = 240	#in minutes, 4 hours, when coverage is defaultCoverage
                        maxMergeAlignmentWalltime = 2980	#in minutes, 2 days
                        minMergeAlignmentMaxMemory = 16000	#in MB, when coverage is defaultCoverage
                        maxMergeAlignmentMaxMemory = 120000	#in MB

                        mergeAlignmentWalltime = self.scaleJobWalltimeOrMemoryBasedOnInput(
                            realInputVolume=oneLibraryCumulativeBaseCount, \
                            baseInputVolume=baseCoverage,
                            baseJobPropertyValue=minMergeAlignmentWalltime,
                            minJobPropertyValue=minMergeAlignmentWalltime,
                            maxJobPropertyValue=maxMergeAlignmentWalltime).value
                        mergeAlignmentMaxMemory = self.scaleJobWalltimeOrMemoryBasedOnInput(
                            realInputVolume=oneLibraryCumulativeBaseCount, \
                            baseInputVolume=baseCoverage,
                            baseJobPropertyValue=minMergeAlignmentMaxMemory,
                            minJobPropertyValue=minMergeAlignmentMaxMemory,
                            maxJobPropertyValue=maxMergeAlignmentMaxMemory).value
                        markDuplicateWalltime= mergeAlignmentWalltime
                        markDuplicateMaxMemory = mergeAlignmentMaxMemory
                        fileBasenamePrefix = utils.getRealPrefixSuffixOfFilenameWithVariableSuffix(\
                            os.path.basename(oneLibraryAlignmentEntry.constructRelativePath()))[0]
                        mergedBamFile = File(os.path.join(oneLibraryAlignmentFolder, \
                            '%s_%s_merged.bam'%(fileBasenamePrefix, library)))
                        alignmentMergeJob, bamIndexJob = self.addAlignmentMergeJob(
                            AlignmentJobAndOutputLs=oneLibraryAlignmentJobAndOutputLs,
                            outputBamFile=mergedBamFile, \
                            samtools=samtools, java=java, \
                            MergeSamFilesJava=MergeSamFilesJava,
                            MergeSamFilesJar=MergeSamFilesJar,
                            BuildBamIndexFilesJava=self.IndexMergedBamIndexJava, \
                            BuildBamIndexJar=BuildBamIndexJar, \
                            mv=mv, namespace=namespace, version=version, \
                            transferOutput=False, \
                            job_max_memory=mergeAlignmentMaxMemory, walltime=mergeAlignmentWalltime, \
                            parentJobLs=[oneLibraryAlignmentFolderJob])

                        finalBamFile = File(os.path.join(oneLibraryAlignmentFolder, \
                            '%s_%s_dupMarked.bam'%(fileBasenamePrefix, library)))
                        markDupJob, markDupBamIndexJob = self.addMarkDupJob(
                            parentJobLs=[alignmentMergeJob, bamIndexJob], \
                            inputBamF=alignmentMergeJob.output, \
                            inputBaiF=bamIndexJob.output, outputBamFile=finalBamFile,\
                            MarkDuplicatesJava=MarkDuplicatesJava,
                            MarkDuplicatesJar=MarkDuplicatesJar, tmpDir=tmpDir,\
                            BuildBamIndexFilesJava=self.IndexMergedBamIndexJava,
                            BuildBamIndexJar=BuildBamIndexJar, \
                            namespace=namespace, version=version,
                            job_max_memory=markDuplicateMaxMemory, \
                            walltime=markDuplicateWalltime, transferOutput=False)
                        no_of_merging_jobs += 1

                        if self.local_realigned:
                            alignmentData = PassingData(
                                jobLs=[markDupJob, markDupBamIndexJob],
                                bamF=markDupJob.output,
                                baiF=markDupBamIndexJob.output, 
                                alignment=oneLibraryAlignmentEntry)
                            preDBAlignmentJob, preDBAlignmentIndexJob = self.addLocalRealignmentSubWorkflow(
                                chr2IntervalDataLs=chr2IntervalDataLs, \
                                registerReferenceData=registerReferenceData, \
                                alignmentData=alignmentData,\
                                inputBamF=markDupJob.output, \
                                outputBamF=None, \
                                parentJobLs=[markDupJob, markDupBamIndexJob], \
                                outputDirPrefix='%s_%s_localRealignment'%(
                                    fileBasenamePrefix, library),
                                transferOutput=False)
                        else:
                            preDBAlignmentJob = markDupJob
                            preDBAlignmentIndexJob = markDupBamIndexJob
                        #2012.9.19 add/copy the alignment file to db-affliated storage
                        #add the metric file to AddAlignmentFile2DB.py as well
                        #  (to be moved into db-affiliated storage)
                        logFile = File(os.path.join(oneLibraryAlignmentFolder,
                            '%s_%s_2db.log'%(fileBasenamePrefix, library)))
                        alignment2DBJob = self.addAddAlignmentFile2DBJob(
                            executable=self.AddAlignmentFile2DB, \
                            inputFile=preDBAlignmentJob.output, \
                            otherInputFileList=[], \
                            individual_alignment_id=oneLibraryAlignmentEntry.id, \
                            individual_sequence_file_raw_id=minIsqFileRawID,\
                            format=None, local_realigned=self.local_realigned,\
                            logFile=logFile, data_dir=data_dir, \
                            parentJobLs=[preDBAlignmentJob, preDBAlignmentIndexJob], \
                            extraDependentInputLs=[preDBAlignmentIndexJob.output,], \
                            extraArguments=None, transferOutput=transferOutput, \
                            job_max_memory=2000, walltime=max(180, markDuplicateWalltime/2), \
                            sshDBTunnel=self.needSSHDBTunnel, commit=True)

                if AlignmentJobAndOutputLs and not skipIndividualAlignment:
                    baseCoverage = 4	#baseline
                    actualCoverage = getattr(individual_sequence, 'coverage', baseCoverage)
                    minMergeAlignmentWalltime = 240
                    #in minutes, 4 hours, when coverage is defaultCoverage
                    maxMergeAlignmentWalltime = 2880	#in minutes, 2 days
                    minMergeAlignmentMaxMemory = 16000
                    #in MB, when coverage is defaultCoverage
                    maxMergeAlignmentMaxMemory = 120000	#in MB

                    mergeAlignmentWalltime = self.scaleJobWalltimeOrMemoryBasedOnInput(
                        realInputVolume=actualCoverage,
                        baseInputVolume=baseCoverage,
                        baseJobPropertyValue=minMergeAlignmentWalltime,
                        minJobPropertyValue=minMergeAlignmentWalltime,
                        maxJobPropertyValue=maxMergeAlignmentWalltime).value
                    mergeAlignmentMaxMemory = self.scaleJobWalltimeOrMemoryBasedOnInput(
                        realInputVolume=actualCoverage,
                        baseInputVolume=baseCoverage,
                        baseJobPropertyValue=minMergeAlignmentMaxMemory, \
                        minJobPropertyValue=minMergeAlignmentMaxMemory,
                        maxJobPropertyValue=maxMergeAlignmentMaxMemory).value
                    markDuplicateWalltime= mergeAlignmentWalltime
                    markDuplicateMaxMemory = mergeAlignmentMaxMemory

                    #2012.3.29	merge alignment output only when there is something to merge!
                    fileBasenamePrefix = utils.getRealPrefixSuffixOfFilenameWithVariableSuffix(\
                        os.path.basename(individual_alignment.constructRelativePath()))[0]
                    mergedBamFile = File(os.path.join(alignmentFolder, \
                        '%s_merged.bam'%(fileBasenamePrefix)))
                    alignmentMergeJob, bamIndexJob = self.addAlignmentMergeJob(
                        AlignmentJobAndOutputLs=AlignmentJobAndOutputLs, \
                        outputBamFile=mergedBamFile, \
                        samtools=samtools, java=java, \
                        MergeSamFilesJava=MergeSamFilesJava,
                        MergeSamFilesJar=MergeSamFilesJar, \
                        BuildBamIndexFilesJava=self.IndexMergedBamIndexJava,
                        BuildBamIndexJar=BuildBamIndexJar, \
                        mv=mv, namespace=namespace, version=version, \
                        parentJobLs=[alignmentFolderJob],\
                        transferOutput=False, job_max_memory=mergeAlignmentMaxMemory, \
                        walltime=mergeAlignmentWalltime)
                    #relative path in the scratch
                    finalBamFile = File(os.path.join(alignmentFolder, '%s_dupMarked.bam'%(fileBasenamePrefix)))

                    markDupJob, markDupBamIndexJob = self.addMarkDupJob(
                        parentJobLs=[alignmentMergeJob, bamIndexJob],
                        inputBamF=alignmentMergeJob.output, \
                        inputBaiF=bamIndexJob.output, outputBamFile=finalBamFile,\
                        MarkDuplicatesJava=MarkDuplicatesJava,
                        MarkDuplicatesJar=MarkDuplicatesJar, tmpDir=tmpDir,\
                        BuildBamIndexFilesJava=self.IndexMergedBamIndexJava,
                        BuildBamIndexJar=BuildBamIndexJar, \
                        namespace=namespace, version=version,
                        job_max_memory=markDuplicateMaxMemory, \
                        walltime=markDuplicateWalltime, \
                        transferOutput=False)
                    no_of_merging_jobs += 1


                    if self.local_realigned:
                        alignmentData = PassingData(jobLs=[markDupJob, markDupBamIndexJob],
                            bamF=markDupJob.output,
                            baiF=markDupBamIndexJob.output, alignment=individual_alignment)
                        preDBAlignmentJob, preDBAlignmentIndexJob = \
                            self.addLocalRealignmentSubWorkflow(\
                            chr2IntervalDataLs=chr2IntervalDataLs,
                            registerReferenceData=registerReferenceData, \
                            alignmentData=alignmentData,\
                            inputBamF=markDupJob.output, \
                            outputBamF=None, \
                            parentJobLs=[markDupJob, markDupBamIndexJob], \
                            outputDirPrefix='%s_localRealignment'%(fileBasenamePrefix),
                            transferOutput=False)
                    else:
                        preDBAlignmentJob = markDupJob
                        preDBAlignmentIndexJob = markDupBamIndexJob
                    #2012.9.19 add/copy the alignment file to db-affliated storage
                    #add the metric file to AddAlignmentFile2DB.py as well
                    #  (to be moved into db-affiliated storage)
                    logFile = File(os.path.join(alignmentFolder, '%s_2db.log'%(fileBasenamePrefix)))
                    alignment2DBJob = self.addAddAlignmentFile2DBJob(
                        executable=self.AddAlignmentFile2DB, \
                        inputFile=preDBAlignmentJob.output, \
                        otherInputFileList=[],\
                        individual_alignment_id=individual_alignment.id, \
                        format=None, local_realigned=self.local_realigned,\
                        logFile=logFile, data_dir=data_dir, \
                        parentJobLs=[preDBAlignmentJob, preDBAlignmentIndexJob], \
                        extraDependentInputLs=[preDBAlignmentIndexJob.output], \
                        extraArguments=None, transferOutput=transferOutput, \
                        job_max_memory=2000,
                        walltime=max(180, markDuplicateWalltime/2),
                        sshDBTunnel=self.needSSHDBTunnel, commit=True)

        sys.stderr.write("%s alignment jobs; %s merge alignment jobs; %s jobs in total.\n"%(\
            no_of_alignment_jobs, no_of_merging_jobs,\
            self.no_of_jobs))

    def run(self):
        """
        2011-7-11
        """

        pdata = self.setup_run()
        workflow = pdata.workflow

        chr2IntervalDataLs = self.getChr2IntervalDataLsBySplitChrSize(
            chr2size=self.chr2size,
            intervalSize=self.intervalSize, \
            intervalOverlapSize=self.intervalOverlapSize)

        #individualSequenceID2FilePairLs = 
        # db_main.getIndividualSequenceID2FilePairLs(self.ind_seq_id_ls,
        #  data_dir=self.local_data_dir)
        isqLs = self.db_main.getISQDBEntryLsForAlignment(self.ind_seq_id_ls,
            data_dir=self.data_dir, filtered=self.sequence_filtered,
            ignoreEmptyReadFile=self.ignoreEmptyReadFile)
        isqLs = self.db_main.filterIndividualSequenceList(
            individual_sequence_list=isqLs,
            min_coverage=self.sequence_min_coverage,\
            max_coverage=self.sequence_max_coverage, \
            individual_site_id_set=set(self.site_id_ls),
            individual_id_set=None,
            sequence_type_id_set=set(self.sequence_type_id_ls),\
            sequencer_id_set=set(self.sequencer_id_ls),
            sequence_filtered=self.sequence_filtered,\
            sequence_batch_id_set=set(self.sequence_batch_id_ls),
            parent_individual_sequence_id_set=None, \
            version_set=set(self.version_ls),\
            country_id_set=set(self.country_id_ls),
            tax_id_set=set(self.tax_id_ls),
            excludeContaminant=self.excludeContaminant,
            report=True)

        refSequence = self.db_main.queryTable(SunsetDB.IndividualSequence).\
            get(self.ref_ind_seq_id)
        #2011-11-16 new way of registering reference fasta file.
        #  but still dont' want to trasnfer 7Gb of data
        refFastaFname = os.path.join(self.data_dir, refSequence.path)
        registerReferenceData = self.registerRefFastaFile(refFastaFname,
            registerAffiliateFiles=True, 
            input_site_handler=self.input_site_handler,\
            checkAffiliateFileExistence=True)

        self.addAllAlignmentJobs(db_main=self.db_main,
            individualSequenceID2FilePairLs=None, \
            isqLs = isqLs,\
            data_dir=self.data_dir,\
            refSequence=refSequence,
            registerReferenceData=registerReferenceData,
            chr2IntervalDataLs=chr2IntervalDataLs,\
            bwa=self.bwa_path,
            additionalArguments=self.additionalArguments,
            samtools=self.samtools, \
            mkdirWrap=self.mkdirWrap, mv=self.cp, \
            java=self.java, \
            MergeSamFilesJava=self.MergeSamFilesJava,
            MergeSamFilesJar=self.PicardJar, \
            MarkDuplicatesJava=self.MarkDuplicatesJava,
            MarkDuplicatesJar=self.PicardJar, tmpDir=self.tmpDir,
            BuildBamIndexFilesJava=self.BuildBamIndexFilesJava,
            BuildBamIndexJar=self.PicardJar, \
            SortSamFilesJava=self.SortSamFilesJava,
            SortSamJar=self.PicardJar,
            AddOrReplaceReadGroupsJava=self.AddOrReplaceReadGroupsJava,
            AddOrReplaceReadGroupsJar=self.PicardJar,\
            alignment_method_name=self.alignment_method_name,
            alignment_format='bam',\
            namespace=self.namespace, version=self.version,
            transferOutput=self.stageOutFinalOutput,\
            PEAlignmentByBWA=self.PEAlignmentByBWA,
            ShortSEAlignmentByBWA=self.ShortSEAlignmentByBWA,
            LongSEAlignmentByBWA=self.LongSEAlignmentByBWA,
            no_of_aln_threads=self.no_of_aln_threads,
            stampy=None, skipDoneAlignment=self.skipDoneAlignment,
            alignmentPerLibrary=self.alignmentPerLibrary, outputDirPrefix="")

        self.end_run()

if __name__ == '__main__':
    main_class = ShortRead2AlignmentWorkflow
    po = ProcessOptions(sys.argv, main_class.option_default_dict,
        error_doc=main_class.__doc__)
    instance = main_class(**po.long_option2value)
    instance.run()
