#!/usr/bin/env python3
"""
abstract class for pegasus workflows that work on alignment files (already aligned).
"""
import sys, os, math
import copy
import logging
from pegaflow.DAX3 import Executable, File, PFN, Link, Job
import pegaflow
from palos import ProcessOptions, getListOutOfStr, PassingData, utils
from . MapReduceGenomeFileWorkflow import MapReduceGenomeFileWorkflow

ParentClass = MapReduceGenomeFileWorkflow
class AbstractAlignmentWorkflow(ParentClass):
    __doc__ = __doc__
    option_default_dict = copy.deepcopy(ParentClass.option_default_dict)
    commonAlignmentWorkflowOptionDict = {
        ('ind_seq_id_ls', 0, ): ['', 'i', 1, 
            'a comma/dash-separated list of IndividualSequence.id.'],\
        ('ind_aln_id_ls', 0, ): ['', '', 1, 
            'a comma/dash-separated list of IndividualAlignment.id. '
            'This overrides ind_seq_id_ls.', ],\
        ('alignment_outdated_index', 0, int): [0, '', 1, 
            'filter based on value of IndividualAlignment.outdated_index.', ],\
        ("alignment_method_id", 0, int): [None, 'G', 1,
            'To filter alignments. None: whatever; integer: AlignmentMethod.id'],\
        ("local_realigned", 0, int): [None, '', 1, 
            'To filter which input alignments to fetch from db (i.e. '
            'AlignmentReadBaseQualityRecalibration.py)'
            'OR to instruct whether local_realigned should be applied '
            '(i.e. ShortRead2Alignment.py)'],\
        ('defaultSampleAlignmentDepth', 1, int): [10, '', 1, 
            "when database doesn't have median_depth info for one alignment, "
            "use this number instead.", ],\
        ('individual_sequence_file_raw_id_type', 1, int): [1, '', 1, 
            "1: only all-library-fused libraries,\n"
            "2: only library-specific alignments,\n"
            "3: both all-library-fused and library-specific alignments", ],\
        }
    option_default_dict.update(commonAlignmentWorkflowOptionDict)
    partitionWorkflowOptionDict= {
        ("selectedRegionFname", 0, ): ["", 'R', 1, 
            'the file is in bed format, tab-delimited, chr start stop. '
            'used to restrict SAMtools/GATK to only make calls at this region. '
            'start and stop are 0-based. i.e. start=0, stop=100 means bases from 0-99. '
            'This overrides the contig/chromosome selection approach defined'
            ' by --contigMaxRankBySize and --contigMinRankBySize. '
            'This file would be split into maxNoOfRegionsPerJob lines.'],\
        ('maxNoOfRegionsPerJob', 1, int): [5000, 'K', 1, 
            'Given selectedRegionFname, this dictates the maximum number of '
            'regions each job would handle. '
            'The actual number could be lower because the regions are first '
            'grouped into chromosomes. '
            'If one chromosome has <maxNoOfRegionsPerJob, then that job handles less.', ],\
        }
    option_default_dict.update(partitionWorkflowOptionDict)

    def __init__(self,
        pegasusFolderName='input',
        output_path=None,

        drivername='postgresql', hostname='localhost',
        dbname='', schema='public', port=None,
        db_user=None,
        db_passwd=None,
        data_dir=None, local_data_dir=None,

        ref_ind_seq_id=None,

        ind_seq_id_ls=None,
        ind_aln_id_ls=None,
        local_realigned=0,
        completedAlignment=None,
        skipDoneAlignment=False,
        excludeContaminant=False,
        sequence_filtered=None,

        alignment_outdated_index=0,
        alignment_method_id=None,
        defaultSampleAlignmentDepth=10,
        individual_sequence_file_raw_id_type=1,

        selectedRegionFname=None,
        maxNoOfRegionsPerJob=5000,
        
        samtools_path="bin/samtools",
        picard_dir="script/picard/dist",
        gatk_path="bin/GenomeAnalysisTK1_6_9.jar",
        gatk2_path="bin/GenomeAnalysisTK.jar",
        picard_path="script/picard.broad/build/libs/picard.jar",
        tabixPath="bin/tabix",
        vcftoolsPath="bin/vcftools/vcftools",
        ligateVcfPerlPath="bin/ligateVcf.pl",
        maxContigID=None,
        minContigID=None,
        contigMaxRankBySize=2500,
        contigMinRankBySize=1,

        chromosome_type_id=None, 
        ref_genome_tax_id=9606,
        ref_genome_sequence_type_id=1,
        ref_genome_version=15,
        ref_genome_outdated_index=0,

        mask_genotype_method_id=None, 
        checkEmptyVCFByReading=False,

        needFastaIndexJob=False,
        needFastaDictJob=False,
        reduce_reads=None,

        site_id_ls="",
        country_id_ls="",
        tax_id_ls="9606",
        sequence_type_id_ls="",
        sequencer_id_ls="",
        sequence_batch_id_ls="",
        version_ls="",
        sequence_min_coverage=None,
        sequence_max_coverage=None,
        
        alignmentDepthIntervalMethodShortName=None,
        minAlignmentDepthIntervalLength=1000,
        alignmentDepthMaxFold=2,
        alignmentDepthMinFold=0.1,

        intervalOverlapSize=500000,
        intervalSize=5000000,
        
        defaultGATKArguments=\
        " --unsafe ALL --validation_strictness SILENT --read_filter BadCigar ",
        
        tmpDir='/tmp/',
        max_walltime=4320,
        
        home_path=None,
        javaPath=None,
        pymodulePath="src/pymodule",
        thisModulePath=None,
        jvmVirtualByPhysicalMemoryRatio=1.2,
        
        site_handler='condor',
        input_site_handler='condor',
        cluster_size=30,

        needSSHDBTunnel=False,
        commit=False,
        debug=False, report=False):
        """
        """
        ParentClass.__init__(self,
            pegasusFolderName=pegasusFolderName,
            output_path=output_path,
            
            drivername=drivername, hostname=hostname,
            dbname=dbname, schema=schema, port=port,
            db_user=db_user, db_passwd=db_passwd,
            data_dir=data_dir, local_data_dir=local_data_dir,

            ref_ind_seq_id=ref_ind_seq_id,

            completedAlignment=completedAlignment,
            skipDoneAlignment=skipDoneAlignment,
            excludeContaminant=excludeContaminant,
            sequence_filtered=sequence_filtered,

            samtools_path=samtools_path,
            picard_dir=picard_dir,
            gatk_path=gatk_path,
            gatk2_path=gatk2_path,
            picard_path=picard_path,
            tabixPath=tabixPath,
            vcftoolsPath=vcftoolsPath,
            ligateVcfPerlPath=ligateVcfPerlPath,

            maxContigID=maxContigID,
            minContigID=minContigID,
            contigMaxRankBySize=contigMaxRankBySize,
            contigMinRankBySize=contigMinRankBySize,

            chromosome_type_id=chromosome_type_id, 
            ref_genome_tax_id=ref_genome_tax_id,
            ref_genome_sequence_type_id=ref_genome_sequence_type_id,
            ref_genome_version=ref_genome_version,
            ref_genome_outdated_index=ref_genome_outdated_index,
            
            mask_genotype_method_id=mask_genotype_method_id, 
            checkEmptyVCFByReading=checkEmptyVCFByReading,

            needFastaIndexJob=needFastaIndexJob,
            needFastaDictJob=needFastaDictJob,
            reduce_reads=reduce_reads,

            site_id_ls=site_id_ls,
            country_id_ls=country_id_ls,
            tax_id_ls=tax_id_ls,
            sequence_type_id_ls=sequence_type_id_ls,
            sequencer_id_ls=sequencer_id_ls,
            sequence_batch_id_ls=sequence_batch_id_ls,
            version_ls=version_ls,

            sequence_min_coverage=sequence_min_coverage,
            sequence_max_coverage=sequence_max_coverage,
            alignmentDepthIntervalMethodShortName=alignmentDepthIntervalMethodShortName,
            minAlignmentDepthIntervalLength=minAlignmentDepthIntervalLength,
            
            alignmentDepthMaxFold=alignmentDepthMaxFold,
            alignmentDepthMinFold=alignmentDepthMinFold,

            intervalOverlapSize=intervalOverlapSize,
            intervalSize=intervalSize,
            
            defaultGATKArguments=defaultGATKArguments,

            tmpDir=tmpDir,
            max_walltime=max_walltime, 
            home_path=home_path,
            javaPath=javaPath,
            pymodulePath=pymodulePath,
            thisModulePath=thisModulePath,
            jvmVirtualByPhysicalMemoryRatio=jvmVirtualByPhysicalMemoryRatio,
            
            site_handler=site_handler,
            input_site_handler=input_site_handler,
            cluster_size=cluster_size,

            needSSHDBTunnel=needSSHDBTunnel,
            commit=commit,
            debug=debug, report=report)
        
        self.ind_aln_id_ls = ind_aln_id_ls
        self.ind_seq_id_ls = ind_seq_id_ls
        self.local_realigned = local_realigned
        self.alignment_outdated_index = alignment_outdated_index
        self.alignment_method_id = alignment_method_id
        self.defaultSampleAlignmentDepth = defaultSampleAlignmentDepth
        self.individual_sequence_file_raw_id_type = individual_sequence_file_raw_id_type

        self.selectedRegionFname = selectedRegionFname
        self.maxNoOfRegionsPerJob = maxNoOfRegionsPerJob

        listArgumentName_data_type_ls = [('ind_seq_id_ls', int), 
            ("ind_aln_id_ls", int)]
        ProcessOptions.processListArguments(listArgumentName_data_type_ls, 
            emptyContent=[], class_to_have_attr=self)

        self.needSplitChrIntervalData = True
        self.mapReduceType = 1
        # type 1: split VCF with fixed number of sites
        # type 2: SelectVariants from VCF with fixed-size windows        
        # child classes can turn it off
        self.needGzipPreReduceReturnData = True
        self.needGzipReduceReturnData = True

    def addAlignmentAsInputToJobLs(self, alignmentDataLs=None, jobLs=[],
        jobInputOption=""):
        """
        Used in addGenotypeCallJobs() to add alignment files as input
            to calling jobs.
        """
        for alignmentData in alignmentDataLs:
            alignment = alignmentData.alignment
            parentJobLs = alignmentData.jobLs
            bamF = alignmentData.bamF
            baiF = alignmentData.baiF
            for job in jobLs:
                if jobInputOption:
                    job.addArguments(jobInputOption)
                job.addArguments(bamF)
                #it's either symlink or stage-in
                job.uses(bamF, transfer=True, register=True, link=Link.INPUT)
                job.uses(baiF, transfer=True, register=True, link=Link.INPUT)
                for parentJob in parentJobLs:
                    if parentJob:
                        self.addJobDependency(parentJob=parentJob, childJob=job)

    def addAlignmentAsInputToPlatypusJobLs(self, alignmentDataLs=None, jobLs=[],
        jobInputOption="--bamFiles"):
        """
        2013.05.21 bugfix: pegasus/condor would truncate long single-argument.
        2013.05.16 different from addAlignmentAsInputToJobLs, that platypus is like this:
                --bamFiles=1.bam,2.bam,3.bam
            used in addGenotypeCallJobs() to add alignment files as input to calling jobs
        """
        for job in jobLs:
            if jobInputOption:
                job.addArguments(jobInputOption)
            fileArgumentLs = []
            alignmentFileFolder = None
            for alignmentData in alignmentDataLs:
                alignment = alignmentData.alignment
                parentJobLs = alignmentData.jobLs
                bamF = alignmentData.bamF
                baiF = alignmentData.baiF

                fileArgumentLs.append(bamF.name)
                if alignmentFileFolder is None:
                    alignmentFileFolder = os.path.split(bamF.name)[0]
                #it's either symlink or stage-in
                job.uses(bamF, transfer=True, register=True, link=Link.INPUT)
                job.uses(baiF, transfer=True, register=True, link=Link.INPUT)
                for parentJob in parentJobLs:
                    if parentJob:
                        self.addJobDependency(parentJob=parentJob, childJob=job)
            #if alignmentFileFolder:
            # #2013.05.21 pegasus/condor would truncate long single-argument.
            #	job.addArguments('%s/*.bam'%(alignmentFileFolder))
            #else:
            job.addArguments(','.join(fileArgumentLs))

    def addAddRG2BamJobsAsNeeded(self, alignmentDataLs=None,
        tmpDir="/tmp"):
        """
        2011-9-15
            add a read group only when the alignment doesn't have it according to db record
            DBVervet.pokeBamReadGroupPresence() from misc.py helps to fill in db records if it's unclear.
        2011-9-14
            The read-group adding jobs will have a "move" part that overwrites
                the original bam&bai if site_handler and input_site_handler is same.
            For those alignment files that don't need to. It doesn't matter.
             pegasus will transfer/symlink them.
        """
        print(f"Adding add-read-group2BAM jobs for {len(alignmentDataLs)} "
            f"alignments if read group is not detected ... ", flush=True)
        job_max_memory = 3500	#in MB
        javaMemRequirement = "-Xms128m -Xmx%sm"%job_max_memory
        indexJobMaxMem = 2500

        addRG2BamDir = None
        addRG2BamDirJob = None

        no_of_rg_jobs = 0
        returnData = []
        for alignmentData in alignmentDataLs:
            alignment = alignmentData.alignment
            parentJobLs = alignmentData.jobLs
            bamF = alignmentData.bamF
            baiF = alignmentData.baiF
            if alignment.read_group_added!=1:
                if addRG2BamDir is None:
                    addRG2BamDir = "addRG2Bam"
                    addRG2BamDirJob = self.addMkDirJob(outputDir=addRG2BamDir)
                # add RG to this bam
                sequencer = alignment.individual_sequence.sequencer
                read_group = alignment.getReadGroup()
                if sequencer=='454':
                    platform_id = 'LS454'
                elif sequencer=='GA':
                    platform_id = 'ILLUMINA'
                else:
                    platform_id = 'ILLUMINA'
                addRGJob = Job(namespace=self.namespace,
                    name=self.AddOrReplaceReadGroupsJava.name,
                    version=self.version)
                outputRGSAM = File(os.path.join(addRG2BamDir,\
                    os.path.basename(alignment.path)))

                addRGJob.addArguments(javaMemRequirement,
                    '-jar', self.AddOrReplaceReadGroupsJar,
                    "INPUT=", bamF,
                    'RGID=%s'%(read_group), 'RGLB=%s'%(platform_id),
                    'RGPL=%s'%(platform_id),
                    'RGPU=%s'%(read_group), 'RGSM=%s'%(read_group),
                    'OUTPUT=', outputRGSAM, 'SORT_ORDER=coordinate',
                    "VALIDATION_STRINGENCY=LENIENT")
                    #(adding the SORT_ORDER doesn't do sorting but it marks the header
                    #  as sorted so that BuildBamIndexJar won't fail.)
                self.addJobUse(addRGJob, file=self.AddOrReplaceReadGroupsJar,
                    transfer=True, register=True, link=Link.INPUT)
                if self.tmpDir:
                    addRGJob.addArguments("TMP_DIR=%s"%self.tmpDir)
                addRGJob.uses(bamF, transfer=True, register=True, link=Link.INPUT)
                addRGJob.uses(baiF, transfer=True, register=True, link=Link.INPUT)
                addRGJob.uses(outputRGSAM, transfer=True, register=True, link=Link.OUTPUT)
                pegaflow.setJobResourceRequirement(addRGJob,
                    job_max_memory=job_max_memory)
                for parentJob in parentJobLs:
                    if parentJob:
                        self.depends(parent=parentJob, child=addRGJob)
                self.addJob(addRGJob)


                index_sam_job = self.addBAMIndexJob(
                    inputBamF=outputRGSAM, parentJobLs=[addRGJob],
                    transferOutput=True, javaMaxMemory=2000)
                newAlignmentData = PassingData(alignment=alignment)
                newAlignmentData.jobLs = [index_sam_job, addRGJob]
                newAlignmentData.bamF = index_sam_job.bamFile
                newAlignmentData.baiF = index_sam_job.baiFile
                no_of_rg_jobs += 1
            else:
                newAlignmentData = alignmentData
            returnData.append(newAlignmentData)
        print(f"{no_of_rg_jobs} alignments need read-group addition.",
            flush=True)
        return returnData

    def preReduce(self, passingData=None, transferOutput=True, **keywords):
        """
        setup additional mkdir folder jobs, before mapEachAlignment,
            mapEachChromosome, mapReduceOneAlignment
        """
        returnData = PassingData(no_of_jobs = 0)
        returnData.jobDataLs = []
        return returnData

    def mapEachChromosome(self, alignmentData=None, chromosome=None,\
        VCFJobData=None, passingData=None,
        reduceBeforeEachAlignmentData=None, transferOutput=True, **keywords):
        """
        2012.9.17
        """
        returnData = PassingData(no_of_jobs = 0)
        returnData.jobDataLs = []
        return returnData

    def map(self, alignmentData=None, intervalData=None,\
        VCFJobData=None, passingData=None,
        mapEachChromosomeData=None, transferOutput=True, **keywords):
        """
        2012.9.17
        """
        returnData = PassingData(no_of_jobs = 0)
        returnData.jobDataLs = []
        return returnData

    def mapEachInterval(self, **keywords):
        """
        2012.9.22 link to map()
        """
        return self.map(**keywords)


    def linkMapToReduce(self, mapEachIntervalData=None,
        preReduceReturnData=None, passingData=None, transferOutput=True, **keywords):
        """
        """
        returnData = PassingData(no_of_jobs = 0)
        returnData.jobDataLs = []
        return returnData

    def mapEachAlignment(self, alignmentData=None,  passingData=None,
        transferOutput=True, **keywords):
        """
        2012.9.22
            similar to reduceBeforeEachAlignmentData() but
             for mapping programs that run on one alignment each.

            passingData.alignmentJobAndOutputLs = []
            passingData.bamFnamePrefix = bamFnamePrefix
            passingData.individual_alignment = alignment
        """
        returnData = PassingData(no_of_jobs = 0)
        returnData.jobDataLs = []

        topOutputDirJob = passingData.topOutputDirJob
        refFastaF = passingData.refFastaFList[0]

        alignment = alignmentData.alignment
        parentJobLs = alignmentData.jobLs
        bamF = alignmentData.bamF
        baiF = alignmentData.baiF

        bamFnamePrefix = alignment.getReadGroup()

        return returnData

    def reduceAfterEachChromosome(self, chromosome=None, passingData=None,
        transferOutput=True,
        mapEachIntervalDataLs=None, **keywords):
        """
        """
        returnData = PassingData(no_of_jobs = 0)
        returnData.jobDataLs = []
        returnData.mapEachIntervalDataLs = mapEachIntervalDataLs
        return returnData

    def reduceBeforeEachAlignment(self, passingData=None,
        transferOutput=True, **keywords):
        """
        2012.9 setup some reduce jobs before loop over all intervals of one alignment begins.
            these reduce jobs will collect stuff from each map() job.
            the link will be established in linkMapToReduce().
        """
        returnData = PassingData(no_of_jobs = 0)
        returnData.jobDataLs = []
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
        2012.9.17
        """
        returnData = PassingData(no_of_jobs = 0)
        returnData.jobDataLs = []
        returnData.reduceAfterEachAlignmentDataLs = reduceAfterEachAlignmentDataLs
        return returnData

    def mapReduceOneAlignment(self, alignmentData=None, passingData=None,
        chrIDSet=None, chrSizeIDList=None, chr2IntervalDataLs=None,
        chr2VCFJobData=None,
        outputDirPrefix=None, transferOutput=False,
        skipChromosomeIfVCFMissing=False, **keywords):
        """
        2013.04.11 moved from AbstractAlignmentAndVCFWorkflow.py
        2013.04.08, added skipChromosomeIfVCFMissing
        2013.1.25
        """
        returnData = PassingData()
        mapEachChromosomeDataLs = passingData.mapEachChromosomeDataLs
        mapEachChromosomeDataLs = []
        reduceBeforeEachAlignmentData = passingData.reduceBeforeEachAlignmentData
        mapEachAlignmentData = passingData.mapEachAlignmentData
        preReduceReturnData = passingData.preReduceReturnData

        for chromosomeSize, chromosome in chrSizeIDList:
            if chr2IntervalDataLs:
                intervalDataLs = chr2IntervalDataLs.get(chromosome, None)
            else:
                intervalDataLs = None
            if chr2VCFJobData:
                VCFJobData = chr2VCFJobData.get(chromosome)
            else:
                VCFJobData = None
            if VCFJobData is None:
                if self.report:
                    logging.warn(f"No VCFJobData for chromosome {chromosome}.")
                if skipChromosomeIfVCFMissing:
                    continue
                VCFJobData = PassingData(job=None, jobLs=[],
                    vcfFile=None, tbi_F=None, file=None, fileLs=[])
                VCFFile = None
            else:
                VCFFile = VCFJobData.file
                if VCFFile is None:
                    if self.report:
                        logging.warn(f"No VCFFile for chromosome {chromosome}.")
                    if skipChromosomeIfVCFMissing:
                        continue
            passingData.chromosome = chromosome
            mapEachChromosomeData = self.mapEachChromosome(
                alignmentData=alignmentData, chromosome=chromosome, \
                VCFJobData=VCFJobData, passingData=passingData,
                reduceBeforeEachAlignmentData=reduceBeforeEachAlignmentData,\
                mapEachAlignmentData=mapEachAlignmentData,\
                transferOutput=False, **keywords)
            passingData.mapEachChromosomeData = mapEachChromosomeData
            mapEachChromosomeDataLs.append(mapEachChromosomeData)

            mapEachIntervalDataLs = passingData.mapEachIntervalDataLs
            mapEachIntervalDataLs = []

            if intervalDataLs:
                for intervalData in intervalDataLs:
                    if intervalData.file:
                        mpileupInterval = intervalData.interval
                        bcftoolsInterval = intervalData.file
                    else:
                        mpileupInterval = intervalData.interval
                        bcftoolsInterval = intervalData.interval
                    intervalFileBasenameSignature = intervalData.intervalFileBasenameSignature
                    overlapInterval = intervalData.overlapInterval
                    overlapFileBasenameSignature = intervalData.overlapIntervalFileBasenameSignature

                    mapEachIntervalData = self.mapEachInterval(
                        alignmentData=alignmentData,
                        intervalData=intervalData,
                        chromosome=chromosome,\
                        VCFJobData=VCFJobData, passingData=passingData,
                        reduceBeforeEachAlignmentData=reduceBeforeEachAlignmentData,
                        mapEachAlignmentData=mapEachAlignmentData,\
                        mapEachChromosomeData=mapEachChromosomeData,
                        transferOutput=False, **keywords)
                    passingData.mapEachIntervalData = mapEachIntervalData
                    mapEachIntervalDataLs.append(mapEachIntervalData)

                    linkMapToReduceData = self.linkMapToReduce(
                        mapEachIntervalData=mapEachIntervalData, \
                        preReduceReturnData=preReduceReturnData, \
                        reduceBeforeEachAlignmentData=reduceBeforeEachAlignmentData,
                        mapEachAlignmentData=mapEachAlignmentData,\
                        passingData=passingData, \
                        **keywords)

            reduceAfterEachChromosomeData = self.reduceAfterEachChromosome(
                chromosome=chromosome,
                passingData=passingData,
                mapEachIntervalDataLs=passingData.mapEachIntervalDataLs,
                transferOutput=False, data_dir=self.data_dir, \
                **keywords)
            passingData.reduceAfterEachChromosomeData = reduceAfterEachChromosomeData
            passingData.reduceAfterEachChromosomeDataLs.append(
                reduceAfterEachChromosomeData)

            gzipReduceAfterEachChromosomeData = self.addGzipSubWorkflow(
                inputData=reduceAfterEachChromosomeData,
                transferOutput=transferOutput,
                outputDirPrefix="%sreduceAfterEachChromosome"%(outputDirPrefix),
                topOutputDirJob=passingData.gzipReduceAfterEachChromosomeFolderJob,
                report=False)
            passingData.gzipReduceAfterEachChromosomeFolderJob = \
                gzipReduceAfterEachChromosomeData.topOutputDirJob
        return returnData

    def setup_chr(self):
        """
        Use self.chr2size to derive chrIDSet.
        Added chrSizeIDList in return
            set chr2VCFJobData to None.
        """
        chrIDSet = set(self.chr2size.keys())
        chrSizeIDList = [(chromosomeSize, chromosome) for chromosome, 
            chromosomeSize in self.chr2size.items()]
        chrSizeIDList.sort()
        chrSizeIDList.reverse()
        #from big to small
        return PassingData(chrIDSet=chrIDSet, chr2VCFJobData=None,
            chrSizeIDList=chrSizeIDList)

    def addAllJobs(self,
        alignmentDataLs=None, chr2IntervalDataLs=None,
        skipDoneAlignment=False,\
        registerReferenceData=None, \
        needFastaIndexJob=False, needFastaDictJob=False, \
        data_dir=None, no_of_gatk_threads = 1, \
        outputDirPrefix="", transferOutput=True, **keywords):
        """
        2012.7.26
        """
        prePreprocessData = self.setup_chr()
        chrIDSet = prePreprocessData.chrIDSet
        chrSizeIDList = prePreprocessData.chrSizeIDList
        chr2VCFJobData = prePreprocessData.chr2VCFJobData

        print(f"Adding jobs that work on {len(alignmentDataLs)} alignments "
            f"(& possibly VCFs) for {len(chrIDSet)} chromosomes/contigs ...",
            flush=True)
        refFastaFList = registerReferenceData.refFastaFList
        refFastaF = refFastaFList[0]

        topOutputDirJob = self.addMkDirJob(outputDir="%sMap"%(outputDirPrefix))
        self.mapDirJob = topOutputDirJob

        plotOutputDirJob = self.addMkDirJob(outputDir="%sPlot"%(outputDirPrefix))
        self.plotOutputDirJob = plotOutputDirJob

        reduceOutputDirJob = self.addMkDirJob(outputDir="%sReduce"%(outputDirPrefix))
        self.reduceOutputDirJob = reduceOutputDirJob

        if needFastaDictJob or registerReferenceData.needPicardFastaDictJob:
            fastaDictJob = self.addRefFastaDictJob(
                refFastaF=refFastaF)
            refFastaDictF = fastaDictJob.refFastaDictF
        else:
            fastaDictJob = None
            refFastaDictF = registerReferenceData.refPicardFastaDictF

        if needFastaIndexJob or registerReferenceData.needSAMtoolsFastaIndexJob:
            fastaIndexJob = self.addRefFastaFaiIndexJob(refFastaF=refFastaF)
            refFastaIndexF = fastaIndexJob.refFastaIndexF
        else:
            fastaIndexJob = None
            refFastaIndexF = registerReferenceData.refSAMtoolsFastaIndexF

        returnData = PassingData()
        returnData.jobDataLs = []

        #2012.9.22 alignmentJobAndOutputLs is a relic.
        #	but it's similar to mapEachIntervalDataLs but
        #   designed for addAlignmentMergeJob(),
        #	so alignmentJobAndOutputLs gets re-set for every alignment.
        # 	mapEachAlignmentDataLs is never reset.
        #	mapEachChromosomeDataLs is reset right after a new alignment is chosen.
        #	mapEachIntervalDataLs is reset right after each chromosome is chosen.
        #	all reduce dataLs never gets reset.
        passingData = PassingData(alignmentJobAndOutputLs=[], \
            alignmentDataLs = alignmentDataLs,\
            bamFnamePrefix=None, \

            outputDirPrefix=outputDirPrefix, \
            topOutputDirJob=topOutputDirJob,\
            plotOutputDirJob=plotOutputDirJob,\
            reduceOutputDirJob = reduceOutputDirJob,\

            refFastaFList=refFastaFList, \
            registerReferenceData= registerReferenceData,\
            refFastaF=refFastaFList[0],\

            fastaDictJob = fastaDictJob,\
            refFastaDictF = refFastaDictF,\
            fastaIndexJob = fastaIndexJob,\
            refFastaIndexF = refFastaIndexF,\

            chromosome=None,\
            chrIDSet=chrIDSet,\
            chrSizeIDList = chrSizeIDList,\
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
        no_of_alignments_worked_on= 0
        for alignmentData in passingData.alignmentDataLs:
            alignment = alignmentData.alignment
            parentJobLs = alignmentData.jobLs + [fastaDictJob, fastaIndexJob]
            bamF = alignmentData.bamF
            baiF = alignmentData.baiF

            bamFnamePrefix = alignment.getReadGroup()

            passingData.alignmentJobAndOutputLs = []
            passingData.bamFnamePrefix = bamFnamePrefix
            passingData.individual_alignment = alignment
            passingData.alignmentData = alignmentData

            if skipDoneAlignment and self.isThisAlignmentComplete(
                individual_alignment=alignment, data_dir=data_dir):
                continue
            no_of_alignments_worked_on += 1
            mapEachAlignmentData = self.mapEachAlignment(
                alignmentData=alignmentData, passingData=passingData, \
                transferOutput=False, \
                preReduceReturnData=preReduceReturnData, **keywords)
            passingData.mapEachAlignmentDataLs.append(mapEachAlignmentData)
            passingData.mapEachAlignmentData = mapEachAlignmentData

            reduceBeforeEachAlignmentData = self.reduceBeforeEachAlignment(
                passingData=passingData,
                preReduceReturnData=preReduceReturnData, transferOutput=False, \
                **keywords)
            passingData.reduceBeforeEachAlignmentData = reduceBeforeEachAlignmentData
            passingData.reduceBeforeEachAlignmentDataLs.append(reduceBeforeEachAlignmentData)


            mapReduceOneAlignmentReturnData = self.mapReduceOneAlignment(
                alignmentData=alignmentData, \
                passingData=passingData, \
                chrIDSet=chrIDSet, chrSizeIDList=chrSizeIDList, \
                chr2IntervalDataLs=chr2IntervalDataLs,
                chr2VCFJobData=chr2VCFJobData,
                outputDirPrefix=outputDirPrefix, transferOutput=transferOutput)

            reduceAfterEachAlignmentData = self.reduceAfterEachAlignment(\
                mapEachAlignmentData=mapEachAlignmentData,\
                mapEachChromosomeDataLs=passingData.mapEachChromosomeDataLs,\
                reduceAfterEachChromosomeDataLs=passingData.reduceAfterEachChromosomeDataLs,\
                passingData=passingData, \
                transferOutput=False, data_dir=data_dir, **keywords)
            passingData.reduceAfterEachAlignmentData = reduceAfterEachAlignmentData
            passingData.reduceAfterEachAlignmentDataLs.append(reduceAfterEachAlignmentData)

            gzipReduceBeforeEachAlignmentData = self.addGzipSubWorkflow(\
                inputData=reduceBeforeEachAlignmentData, transferOutput=transferOutput,\
                outputDirPrefix="%sReduceBeforeEachAlignment"%(outputDirPrefix), \
                topOutputDirJob=passingData.gzipReduceBeforeEachAlignmentFolderJob, report=False)
            passingData.gzipReduceBeforeEachAlignmentFolderJob = \
                gzipReduceBeforeEachAlignmentData.topOutputDirJob

            gzipReduceAfterEachAlignmentData = self.addGzipSubWorkflow(\
                inputData=reduceAfterEachAlignmentData, transferOutput=transferOutput,\
                outputDirPrefix="%sReduceAfterEachAlignment"%(outputDirPrefix), \
                topOutputDirJob=passingData.gzipReduceAfterEachAlignmentFolderJob, \
                report=False)
            passingData.gzipReduceAfterEachAlignmentFolderJob = \
                gzipReduceAfterEachAlignmentData.topOutputDirJob
        reduceReturnData = self.reduce(passingData=passingData, \
            mapEachAlignmentData=passingData.mapEachAlignmentData, \
            reduceAfterEachAlignmentDataLs=passingData.reduceAfterEachAlignmentDataLs,\
            **keywords)
        passingData.reduceReturnData = reduceReturnData


        #2012.9.18 gzip the final output
        newReturnData = self.addGzipSubWorkflow(inputData=preReduceReturnData,
            transferOutput=transferOutput,\
            outputDirPrefix="%sGzipPreReduce"%(outputDirPrefix), \
            topOutputDirJob=passingData.gzipPreReduceFolderJob, \
            report=False)
        passingData.gzipPreReduceFolderJob = newReturnData.topOutputDirJob
        newReturnData = self.addGzipSubWorkflow(inputData=reduceReturnData,
            transferOutput=transferOutput,\
            outputDirPrefix="%sGzipReduce"%(outputDirPrefix), \
            topOutputDirJob=passingData.gzipReduceFolderJob, \
            report=False)
        passingData.gzipReduceFolderJob = newReturnData.topOutputDirJob

        sys.stderr.write("%s alignments to be worked on. %s jobs.\n"%(
            no_of_alignments_worked_on, self.no_of_jobs))
        return returnData


    def registerExecutables(self):
        """
        """
        ParentClass.registerExecutables(self)
        self.registerOneExecutable(path=self.samtools_path,
            name='samtools', clusterSizeMultiplier=0.2)
        self.samtoolsExecutableFile = self.registerOneExecutableAsFile(
            path=self.samtools_path,
            site_handler=self.input_site_handler)

    def setup_run(self):
        """
        Wrap all standard pre-run() related functions into this function.
            setting up for run(), called by run().
        """
        ParentClass.setup_run(self)
        # ParentClass.setup_run() will call getReferenceSequence() to
        #  setup self.registerReferenceData.
        
        if self.needSplitChrIntervalData:
            #2013.06.21 defined in ParentClass.__init__()
            if self.alignmentDepthIntervalMethodShortName and self.db_main and \
                self.db_main.checkAlignmentDepthIntervalMethod(
                    short_name=self.alignmentDepthIntervalMethodShortName):
                #2013.09.01 fetch intervals from db
                #make sure it exists in db first
                chr2IntervalDataLs = self.getChr2IntervalDataLsFromDBAlignmentDepthInterval(
                    db=self.db_main,
                    intervalSize=self.intervalSize,
                    intervalOverlapSize=self.intervalOverlapSize,
                    alignmentDepthIntervalMethodShortName=self.alignmentDepthIntervalMethodShortName,
                    alignmentDepthMinFold=self.alignmentDepthMinFold,
                    alignmentDepthMaxFold=self.alignmentDepthMaxFold,
                    minAlignmentDepthIntervalLength=self.minAlignmentDepthIntervalLength,
                    minContigID=self.minContigID,
                    maxContigID=self.maxContigID)
            else:
                #split evenly using chromosome size
                chr2IntervalDataLs = self.getChr2IntervalDataLsBySplitChrSize(
                    chr2size=self.chr2size,
                    intervalSize=self.intervalSize, \
                    intervalOverlapSize=self.intervalOverlapSize)
            # 2012.8.2 if maxContigID/minContigID is not well defined.
            #  restrictContigDictionry won't do anything.
            chr2IntervalDataLs = self.restrictContigDictionry(
                dc=chr2IntervalDataLs,
                maxContigID=self.maxContigID, minContigID=self.minContigID)
        else:
            logging.warn(f"self.needSplitChrIntervalData="\
                f"{self.needSplitChrIntervalData}, set chr2IntervalDataLs=None.")
            chr2IntervalDataLs = None

        alignmentLs = self.getAlignments()
        alignmentDataLs = self.registerAlignmentAndItsIndexFile(
            alignmentLs=alignmentLs, data_dir=self.data_dir)
        self.alignmentLs = alignmentLs
        self.alignmentDataLs = alignmentDataLs
        self.chr2IntervalDataLs = chr2IntervalDataLs
        return self

    def run(self):
        """
        2013.1.25
        """

        pdata = self.setup_run()
        self.addAllJobs(alignmentDataLs=pdata.alignmentDataLs,
            chr2IntervalDataLs=pdata.chr2IntervalDataLs,
            skipDoneAlignment=self.skipDoneAlignment,
            registerReferenceData=pdata.registerReferenceData,
            needFastaIndexJob=self.needFastaIndexJob,
            needFastaDictJob=self.needFastaDictJob,
            data_dir=self.data_dir, no_of_gatk_threads = 1,
            transferOutput=True,
            outputDirPrefix=self.pegasusFolderName)

        self.end_run()

if __name__ == '__main__':
    from argparse import ArgumentParser
    ap = ArgumentParser()
    ap.add_argument("--completedAlignment", type=int, default=None,
        help='To filter incomplete alignments: '
            '--completedAlignment 0 . '
            '--completedAlignment 1 fetches only the complete alignments. '
            'Default (%(default)s) has no effect.')