#!/usr/bin/env python3
"""
Description:
    2013.11.24 a generic workflow that map-reduces inputs of one or
        multiple genomic files (i.e. multi-chromosome, tabix-indexed )
        parent class is AbstractNGSWorkflow.
    
Examples:
    #2013.11.24
    %s --input_path FindNewRefCoordinates_Method109_vs_3488_BWA/folderReduceLiftOverVCF/
        -H -C 10 -j hcondor -l hcondor
        -D /u/home/p/polyacti/NetworkData/vervet/db
        -t /u/home/p/polyacti/NetworkData/vervet/db/
        -o dags/Method109_vs_3488_BWA_F99.sameSiteConcordance.xml
        --notToUseDBToInferVCFNoOfLoci
        --db_user yh -z localhost
"""
import sys, os, math
__doc__ = __doc__%(sys.argv[0])
import logging
import pegaflow
from pegaflow.DAX3 import Executable, File, PFN
from palos import PassingData, utils
from . AbstractNGSWorkflow import AbstractNGSWorkflow

ParentClass = AbstractNGSWorkflow
class MapReduceGenomeFileWorkflow(ParentClass):
    __doc__ = __doc__
    def __init__(self,
        input_path=None,
        inputSuffixList=None,
        pegasusFolderName='input',
        output_path=None,

        drivername='postgresql', hostname='localhost',
        dbname='', schema='public', port=None,
        db_user=None,
        db_passwd=None,
        data_dir=None, local_data_dir=None,
        
        ref_ind_seq_id=None,

        excludeContaminant=False,
        sequence_filtered=None,
        completedAlignment=None,
        skipDoneAlignment=False,

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
        sequence_max_coverage=None,
        sequence_min_coverage=None,

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
        Default interval is 5Mb.
        Default interval overlap is 500K.
        """
        ParentClass.__init__(self,
            input_path=input_path,
            inputSuffixList=None,
            pegasusFolderName=pegasusFolderName,
            output_path=output_path,

            drivername=drivername, hostname=hostname,
            dbname=dbname, schema=schema, port=port,
            db_user=db_user, db_passwd=db_passwd,
            data_dir=data_dir, local_data_dir=local_data_dir,

            ref_ind_seq_id=ref_ind_seq_id,

            excludeContaminant=excludeContaminant,
            sequence_filtered=sequence_filtered,
            completedAlignment=completedAlignment,
            skipDoneAlignment=skipDoneAlignment,

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

            sequence_max_coverage=sequence_max_coverage,
            sequence_min_coverage=sequence_min_coverage,
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
        
        self.needSplitChrIntervalData = True
        self.mapReduceType = 1
        # type 1: split VCF with fixed number of sites
        # type 2: SelectVariants from VCF with fixed-size windows
        # child classes could change its value in the end of their own __init__()
        
        # child classes can turn it off
        self.needGzipPreReduceReturnData = True
        self.needGzipReduceReturnData = True
    
    def connectDB(self):
        """
        """
        ParentClass.connectDB(self)
        
        self.registerReferenceData = None
        self.refFastaFList= None
    
    def preReduce(self, outputDirPrefix="", passingData=None,
        transferOutput=True, **keywords):
        """
        2013.06.14
            move topOutputDirJob from addAllJobs to here. 
        2012.9.17
        """
        returnData = PassingData(no_of_jobs = 0)
        returnData.jobDataLs = []
        
        self.topOutputDirJob = self.addMkDirJob(
            outputDir="%sRun"%(outputDirPrefix))
        passingData.topOutputDirJob = self.topOutputDirJob
        
        mapDirJob = self.addMkDirJob(outputDir="%sMap"%(outputDirPrefix))
        passingData.mapDirJob = mapDirJob
        returnData.mapDirJob = mapDirJob
        self.mapDirJob = mapDirJob
        
        reduceOutputDirJob = self.addMkDirJob(
            outputDir="%sReduce"%(outputDirPrefix))
        passingData.reduceOutputDirJob = reduceOutputDirJob
        returnData.reduceOutputDirJob = reduceOutputDirJob
        
        self.plotDirJob = self.addMkDirJob(outputDir="%sPlot"%(outputDirPrefix))
        self.statDirJob = self.addMkDirJob(outputDir="%sStat"%(outputDirPrefix))
        self.reduceStatDirJob = self.addMkDirJob(
            outputDir="%sReduceStat"%(outputDirPrefix))
        self.reduceEachInputDirJob = self.addMkDirJob(
            outputDir="%sReduceEachInput"%(outputDirPrefix))
        self.reduceEachChromosomeDirJob = self.addMkDirJob(
            outputDir="%sReduceEachChromosome"%(outputDirPrefix))
        self.reduceOutputDirJob = reduceOutputDirJob
        return returnData
    
    def selectIntervalFromInputFile(self, jobData=None, chromosome=None,\
        intervalData=None, mapEachChromosomeData=None,\
        passingData=None, transferOutput=False,\
        **keywords):
        """
        2013.11.24
        """
        inputSuffix = utils.getRealPrefixSuffix(
            jobData.file.name)[1]
        outputFile = File(os.path.join(self.mapDirJob.output, \
            '%s_%s%s'%(passingData.fileBasenamePrefix, \
            intervalData.overlapInterval, inputSuffix)))
        tabixRetrieveJob = self.addTabixRetrieveJob(
            executable=self.tabixRetrieve, \
            tabixPath=self.tabixPath, \
            inputF=jobData.file, outputF=outputFile, \
            regionOfInterest=intervalData.overlapInterval,
            includeHeader=True,\
            parentJobLs=jobData.jobLs + [self.mapDirJob],
            job_max_memory=100,
            extraDependentInputLs=jobData.fileLs[1:], \
            transferOutput=False)
        return self.constructJobDataFromJob(job=tabixRetrieveJob)
        

    def addAllJobs(self, inputData=None, chr2IntervalDataLs=None,
        data_dir=None,
        intervalSize=3000, intervalOverlapSize=0,
        outputDirPrefix="", passingData=None, \
        transferOutput=True, job_max_memory=2000, **keywords):
        """
        2013.06.14 bugfix regarding noOfUnits,
            which was all inferred from one file.
        2012.7.26
            architect of the whole map-reduce framework
        """
        print(f"Adding jobs for {len(inputData.jobDataLs)} input "
            "genome files ..." , flush=True)
        
        returnData = PassingData()
        returnData.jobDataLs = []
        
        #2012.9.22 
        # 	mapEachAlignmentDataLs is never reset.
        #	mapEachChromosomeDataLs is reset upon new alignment
        #	mapEachIntervalDataLs is reset upon each new chromosome
        #	all reduce lists never get reset.
        #	fileBasenamePrefix is the prefix of input file's basename,
        #    to be used for temporary output files in reduceEachInput()
        #		but not for output files in mapEachInterval()
        passingData = PassingData(\
            fileBasenamePrefix=None, \
            chromosome=None, \
            
            outputDirPrefix=outputDirPrefix, \
            intervalFileBasenamePrefix=None,\
            
            registerReferenceData=None, \
            refFastaFList=None, \
            refFastaF=None,\
            
            fastaDictJob = None,\
            refFastaDictF = None,\
            fastaIndexJob = None,\
            refFastaIndexF = None,\
            
            intervalOverlapSize =intervalOverlapSize,
            intervalSize=intervalSize,
            jobData=None,\
            splitInputFile=None,\
            intervalDataLs=None,\
            preReduceReturnData=None,\
            
            mapEachIntervalData=None,\
            mapEachIntervalDataLs=None,\
            mapEachIntervalDataLsLs=[],\
            mapEachInputData=None,\
            mapEachInputDataLs=None,\
            mapEachInputDataLsLs=[],\
            mapEachChromosomeData=None, \
            mapEachChromosomeDataLs=[], \
            
            chromosome2mapEachIntervalDataLs = {},\
            chromosome2mapEachInputDataLs = {},\
            
            reduceEachInputData=None,\
            reduceEachChromosomeData=None,\
            reduceEachInputDataLs=None,\
            reduceEachInputDataLsLs=[],\
            reduceEachChromosomeDataLs=[],\
            )
        # mapEachIntervalDataLsLs is list of mapEachIntervalDataLs by each Input file.
        # mapEachInputDataLsLs is list of mapEachInputDataLs by each chromosome
        # reduceEachInputDataLsLs is list of reduceEachInputDataLs by each chromosome
        
        preReduceReturnData = self.preReduce(outputDirPrefix=outputDirPrefix,
            passingData=passingData, transferOutput=True,\
            **keywords)
        passingData.preReduceReturnData = preReduceReturnData
        
        #gzip folder jobs (to avoid repeatedly creating the same folder
        gzipReduceEachInputFolderJob = None
        gzipReduceEachChromosomeFolderJob = None
        gzipReduceFolderJob = None
        gzipPreReduceFolderJob = None
        no_of_input_files = 0
        
        firstInterval = True
        
        for chromosome, intervalDataLs in chr2IntervalDataLs.items():
            passingData.chromosome = chromosome
            mapEachChromosomeData = self.mapEachChromosome(
                chromosome=chromosome,
                passingData=passingData,
                transferOutput=False, **keywords)
            passingData.mapEachChromosomeData = mapEachChromosomeData
            passingData.mapEachChromosomeDataLs.append(mapEachChromosomeData)
            
            passingData.mapEachInputDataLsLs.append([])
            #the last one from the double list is the current one
            passingData.mapEachInputDataLs = passingData.mapEachInputDataLsLs[-1]
            passingData.mapEachIntervalDataLs = []
            passingData.chromosome2mapEachIntervalDataLs[chromosome] = []
            
            passingData.reduceEachInputDataLsLs.append([])
            passingData.reduceEachInputDataLs = passingData.reduceEachInputDataLsLs[-1]
            
            for i in range(len(inputData.jobDataLs)):
                jobData = inputData.jobDataLs[i]
                passingData.jobData = jobData
                passingData.inputJobData = jobData
                
                InputFile = jobData.file
                commonFileBasenamePrefix = utils.getFileBasenamePrefixFromPath(
                    InputFile.name)
                passingData.fileBasenamePrefix = commonFileBasenamePrefix
                
                no_of_input_files += 1
                if no_of_input_files%10==0:
                    sys.stderr.write("%s\t%s Inputs."%('\x08'*40, no_of_input_files))
                
                for intervalData in intervalDataLs:
                    selectIntervalJobData = self.selectIntervalFromInputFile(
                        jobData=jobData, chromosome=chromosome,
                        intervalData=intervalData,
                        mapEachChromosomeData=mapEachChromosomeData,
                        passingData=passingData,
                        transferOutput=firstInterval,
                        **keywords)
                    mapEachIntervalData = self.mapEachInterval(
                        inputJobData=jobData,
                        selectIntervalJobData=selectIntervalJobData,
                        chromosome=chromosome,
                        intervalData=intervalData,
                        mapEachChromosomeData=mapEachChromosomeData, \
                        passingData=passingData,
                        transferOutput=firstInterval,
                        **keywords)
                    
                    passingData.mapEachIntervalData = mapEachIntervalData
                    passingData.mapEachIntervalDataLs.append(mapEachIntervalData)
                    passingData.chromosome2mapEachIntervalDataLs[chromosome].append(
                        mapEachIntervalData)
                    
                    linkMapToReduceData = self.linkMapToReduce(
                        mapEachIntervalData=mapEachIntervalData,
                        preReduceReturnData=preReduceReturnData,
                        passingData=passingData,
                        **keywords)
                    if firstInterval==True:
                        firstInterval = False
                reduceEachInputData = self.reduceEachInput(
                    chromosome=chromosome, passingData=passingData,
                    mapEachIntervalDataLs=passingData.mapEachIntervalDataLs,
                    transferOutput=False, data_dir=data_dir, \
                    **keywords)
                passingData.reduceEachInputData = reduceEachInputData
                passingData.reduceEachInputDataLs.append(reduceEachInputData)
                
                gzipReduceEachInputData = self.addGzipSubWorkflow(\
                    inputData=reduceEachInputData,
                    outputDirPrefix="%sReduceEachInput"%(outputDirPrefix),
                    topOutputDirJob=gzipReduceEachInputFolderJob, \
                    transferOutput=transferOutput,
                    report=False)
                gzipReduceEachInputFolderJob = gzipReduceEachInputData.topOutputDirJob
            reduceEachChromosomeData = self.reduceEachChromosome(
                chromosome=chromosome, passingData=passingData, \
                mapEachInputDataLs=passingData.mapEachInputDataLs, \
                chromosome2mapEachIntervalDataLs=passingData.chromosome2mapEachIntervalDataLs,
                reduceEachInputDataLs=passingData.reduceEachInputDataLs,\
                transferOutput=False, data_dir=data_dir, \
                **keywords)
            passingData.reduceEachChromosomeData = reduceEachChromosomeData
            passingData.reduceEachChromosomeDataLs.append(reduceEachChromosomeData)
            
            gzipReduceEachChromosomeData = self.addGzipSubWorkflow(
                inputData=reduceEachChromosomeData,
                outputDirPrefix="%sReduceEachChromosome"%(outputDirPrefix), \
                topOutputDirJob=gzipReduceEachChromosomeFolderJob,
                transferOutput=transferOutput,
                report=False)
            gzipReduceEachChromosomeFolderJob = \
                gzipReduceEachChromosomeData.topOutputDirJob
            
        reduceReturnData = self.reduce(passingData=passingData,
            transferOutput=False,
            mapEachChromosomeDataLs=passingData.mapEachInputDataLs,\
            reduceEachChromosomeDataLs=passingData.reduceEachChromosomeDataLs,
            **keywords)
        passingData.reduceReturnData = reduceReturnData
        
        if self.needGzipPreReduceReturnData:
            gzipPreReduceReturnData = self.addGzipSubWorkflow(
                inputData=preReduceReturnData, transferOutput=transferOutput,
                outputDirPrefix="%sPreReduce"%(outputDirPrefix), \
                topOutputDirJob= gzipPreReduceFolderJob, report=False)
            gzipPreReduceFolderJob = gzipPreReduceReturnData.topOutputDirJob
        
        if self.needGzipReduceReturnData:
            gzipReduceReturnData = self.addGzipSubWorkflow(
                inputData=reduceReturnData, transferOutput=transferOutput,
                outputDirPrefix="%sReduce"%(outputDirPrefix), \
                topOutputDirJob=gzipReduceFolderJob, report=False)
            gzipReduceFolderJob = gzipReduceReturnData.topOutputDirJob
        
        print(f" {no_of_input_files} Input files.", flush=True)
        sys.stderr.write(f"{self.no_of_jobs} jobs.\n")
        return reduceReturnData

    
    def mapEachChromosome(self, chromosome=None,\
        VCFJobData=None, jobData=None, passingData=None,
        transferOutput=True, **keywords):
        """
        2012.9.17
        """
        returnData = PassingData(no_of_jobs = 0)
        returnData.jobDataLs = []
        return returnData
    
    def mapEachInterval(self, inputJobData=None, selectIntervalJobData=None,
        chromosome=None,intervalData=None,\
        mapEachChromosomeData=None, \
        passingData=None, transferOutput=False, **keywords):
        """
        2013.04.08 use inputJobData
        2012.10.3
            #. extract flanking sequences from the input Input
            #   (ref sequence file => contig ref sequence)
            #. blast them
            #. run FindSNPPositionOnNewRefFromFlankingBlastOutput.py
                #. where hit length match query length, and no of
                #  mismatches <=2 => good => infer new coordinates
            #. output a mapping file between old SNP and new SNP coordinates.
                #. reduce this thing by combining everything
            #. make a new Input file based on the input split Input file
                (replace contig ID , position with the new one's,
                 remove the header part regarding chromosomes or replace it)

        """
        returnData = PassingData(no_of_jobs = 0)
        returnData.jobDataLs = []
        #passingData.intervalFileBasenamePrefix
        #passingData.splitInputFile
        #passingData.unitNumber
        """
        ## 2013.06.19 structures available from passingData, specific to the interval
        passingData.splitInputFile = splitInputFile
        passingData.unitNumber = unitNumber
        passingData.intervalFileBasenamePrefix = '%s_%s_splitInput_u%s'%(
            chromosome, commonPrefix, unitNumber)
        passingData.noOfIndividuals = jobData.file.noOfIndividuals
        passingData.span = self.intervalSize + self.intervalOverlapSize*2
        #2013.06.19 for memory/walltime gauging
        """
        return returnData

    def linkMapToReduce(self, mapEachIntervalData=None,
        preReduceReturnData=None, passingData=None, transferOutput=True,
        **keywords):
        """
        """
        returnData = PassingData(no_of_jobs = 0)
        returnData.jobDataLs = []
        return returnData
    
    def reduceEachChromosome(self, chromosome=None, passingData=None,
        mapEachInputDataLs=None, \
        chromosome2mapEachIntervalDataLs=None,\
        reduceEachInputDataLs=None,\
        transferOutput=True, \
        **keywords):
        """
        """
        returnData = PassingData(no_of_jobs = 0)
        returnData.jobDataLs = []
        returnData.mapEachInputDataLs = mapEachInputDataLs
        returnData.reduceEachInputDataLs = reduceEachInputDataLs
        return returnData
    
    def reduceEachInput(self, chromosome=None, passingData=None,
        mapEachIntervalDataLs=None,
        transferOutput=True, **keywords):
        """
        2013.07.10
            #. concatenate all the sub-Inputs into one
        """
        returnData = PassingData(no_of_jobs = 0)
        returnData.jobDataLs = []
        returnData.mapEachIntervalDataLs = mapEachIntervalDataLs
        
        #intervalJobLs = [pdata for pdata in mapEachIntervalDataLs]
        
        """
        realInputVolume = passingData.jobData.file.noOfIndividuals * \
            passingData.jobData.file.noOfLoci
        baseInputVolume = 200*20000
        walltime = self.scaleJobWalltimeOrMemoryBasedOnInput(
            realInputVolume=realInputVolume, \
            baseInputVolume=baseInputVolume, baseJobPropertyValue=60,
            minJobPropertyValue=60, maxJobPropertyValue=500).value
        job_max_memory = self.scaleJobWalltimeOrMemoryBasedOnInput(
            realInputVolume=realInputVolume, \
            baseInputVolume=baseInputVolume, baseJobPropertyValue=5000,
            minJobPropertyValue=5000, maxJobPropertyValue=10000).value
        """
        return returnData
    
    def reduce(self, reduceEachChromosomeDataLs=None, \
        mapEachChromosomeDataLs=None, passingData=None, transferOutput=True,
        **keywords):
        """
        2013.07.18 return each processed-Input job data so that
         followup workflows could carry out map-reduce
        2012.9.17
        """
        returnData = PassingData(no_of_jobs = 0)
        returnData.jobDataLs = []
        returnData.mapEachChromosomeDataLs = mapEachChromosomeDataLs
        returnData.reduceEachChromosomeDataLs = reduceEachChromosomeDataLs
        """
        #2013.07.18 example to return each processed-Input job data
        #  so that followup workflows could carry out map-reduce
        for reduceEachInputDataLs in passingData.reduceEachInputDataLsLs:
            if reduceEachInputDataLs:
                for reduceEachInputData in reduceEachInputDataLs:
                    if reduceEachInputData:
                        returnData.jobDataLs.append(reduceEachInputData.WHATEVERJobData)
        """
        return returnData
    
    def registerExecutables(self):
        """
        """
        ParentClass.registerExecutables(self)
        self.registerOneExecutable(
            path=self.tabixPath, name='tabix', clusterSizeMultiplier=5)
        self.tabixExecutableFile = self.registerOneExecutableAsFile(
            path=self.tabixPath)

    def setup_run(self):
        """
        Wrap all standard pre-run() related functions into this function.
            Setting up for run(), called by run().
        """
        ParentClass.setup_run(self)
        
        #self.chr2size = {}
        #self.chr2size = set(['Contig149'])
        # #temporarily restricting chr2size to Contig149
        if self.needSplitChrIntervalData:
            #self.needSplitChrIntervalData is 
            # defined in __init__().
            chr2IntervalDataLs = self.getChr2IntervalDataLsBySplitChrSize(
                chr2size=self.chr2size,
                intervalSize=self.intervalSize,
                intervalOverlapSize=self.intervalOverlapSize)
        else:
            logging.warn(f"self.needSplitChrIntervalData="\
                f"{self.needSplitChrIntervalData}, set chr2IntervalDataLs=None.")
            chr2IntervalDataLs = None
        inputData = None
        firstInputJobData = None
        if getattr(self, 'input_path', None):
            inputData = self.registerFilesOfInputDir(
                inputDir=self.input_path,
                input_site_handler=self.input_site_handler, \
                pegasusFolderName=self.pegasusFolderName,\
                inputSuffixSet=self.inputSuffixSet,\
                indexFileSuffixSet=set(['.tbi', '.fai']),\
                checkEmptyInputByReading=getattr(
                    self, 'checkEmptyInputByReading', None),\
                maxContigID=self.maxContigID, \
                minContigID=self.minContigID,\
                db_main=getattr(self, 'db_main', None), \
                needToKnowNoOfLoci=getattr(self, 'needToKnowNoOfLoci', True),\
                minNoOfLociInInput=getattr(self, 'minNoOfLociInInput', 10))
            if inputData and inputData.jobDataLs:
                firstInputJobData = inputData.jobDataLs[0]
                firstInputFile = firstInputJobData.file
                print(f"\t Input file {firstInputFile} is chosen as an example"
                    f" input for any job that needs a random Input file.",
                    flush=True)
        
        self.inputData = inputData
        self.chr2IntervalDataLs = chr2IntervalDataLs
        self.firstInputJobData = firstInputJobData
        #self.firstInputFile = firstInputFile
        if self.ref_ind_seq_id:
            self.registerReferenceData = self.getReferenceSequence()
        else:
            self.registerReferenceData = None
        return self
    
    def run(self):
        """
        """
        pdata = self.setup_run()
        inputData = pdata.inputData
        
        if len(inputData.jobDataLs)<=0:
            print(f"No VCF files in this folder, {self.input_path}.", flush=True)
            sys.exit(0)
                
        self.addAllJobs(inputData=inputData,
            chr2IntervalDataLs=self.chr2IntervalDataLs,
            data_dir=self.data_dir,
            intervalSize=self.intervalSize,
            intervalOverlapSize=self.intervalOverlapSize,
            outputDirPrefix=self.pegasusFolderName,
            transferOutput=True,)
        
        self.end_run()

