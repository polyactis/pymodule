#!/usr/bin/env python3
"""
2020/01/29
    an abstract class for pegasus workflows that work on bioinformatic data
"""
import sys, os, math
from palos import ProcessOptions, getListOutOfStr, PassingData, utils
from pegaflow.DAX3 import Executable, File, PFN, Link, Job
from . AbstractWorkflow import AbstractWorkflow

ParentClass = AbstractWorkflow
class AbstractBioinfoWorkflow(ParentClass):
    __doc__ = __doc__
    option_default_dict = ParentClass.option_default_dict.copy()

    def __init__(self,  **keywords):
        """
        20200129
        """
        ParentClass.__init__(self, **keywords)

    def registerPlinkExecutables(self):
        if not workflow:
            workflow = self
        
        self.addExecutableFromPath(path=self.plinkPath, \
            name='plink', clusterSizeMultiplier=1)
        self.addExecutableFromPath(path=self.plinkPath, \
            name='plinkNoClustering', clusterSizeMultiplier=0)

        #2012.8.10 different plinks so that you can differentiate between different types of plink jobs
        self.addExecutableFromPath(path=self.plinkPath, \
            name='plinkMerge', clusterSizeMultiplier=0)
        self.addExecutableFromPath(path=self.plinkPath, \
            name='plinkIBD', clusterSizeMultiplier=0)
        self.addExecutableFromPath(path=self.plinkPath, \
            name='plinkConvert', clusterSizeMultiplier=1)
        self.addExecutableFromPath(path=self.plinkPath, \
            name='plinkLDPrune', clusterSizeMultiplier=1)
        self.addExecutableFromPath(path=self.plinkPath, \
            name='plinkExtract', clusterSizeMultiplier=1)
        #2013.07.24
        self.addExecutableFromPath(path=os.path.join(self.pymodulePath, \
            'mapper/modifier/SplitPlinkLMendelFileSNPIDIntoChrPosition.py'), \
            name='SplitPlinkLMendelFileSNPIDIntoChrPosition', clusterSizeMultiplier=1)
        
        #2013.07.19
        self.addExecutableFromPath(path=os.path.join(self.pymodulePath, \
            'pedigree/CalculateMendelErrorRateGivenPlinkOutput.py'), \
            name='CalculateMendelErrorRateGivenPlinkOutput', clusterSizeMultiplier=1)

    def getTopNumberOfContigs(self, **keywords):
        """
        2013.2.6 placeholder
        """
        pass

    def registerBlastNucleotideDatabaseFile(self, ntDatabaseFname=None,  folderName=""):
        """
        2012.10.8
            moved from BlastWorkflow.py
        2012.5.23
        """
        return self.registerRefFastaFile(refFastaFname=ntDatabaseFname, registerAffiliateFiles=True, \
                                    checkAffiliateFileExistence=True, addPicardDictFile=False, \
                                    affiliateFilenameSuffixLs=['nin', 'nhr', 'nsq'],\
                                    folderName=folderName)

    def registerRefFastaFile(self, refFastaFname=None, registerAffiliateFiles=True, \
                        checkAffiliateFileExistence=True, addPicardDictFile=True,\
                        affiliateFilenameSuffixLs=['fai', 'amb', 'ann', 'bwt', 'pac', 'sa', 'rbwt', 'rpac', 'rsa', \
                        'stidx', 'sthash'], folderName="reference"):
        """
        suffix here doesn't include ".".
        
        2013.08.23 bugfix, check if workflow has a file registered before adding it
        2013.3.26 added refSAMtoolsFastaIndexF, refPicardFastaDictF into returnData
        2013.3.20 deduce needBWARefIndexJob, needSAMtoolsFastaIndexJob, needPicardFastaDictJob, needStampyRefIndexJob from missing suffixes
        2010.10.10 added argument folderName
        2012.5.23
            add an argument "addPicardDictFile" to offer user option to exclude this file (i.e. in registerBlastNucleotideDatabaseFile)
        2012.2.24
            dict is via picard, also required for GATK
            fai is via "samtools faidx" (index reference). also required for GATK
            amb', 'ann', 'bwt', 'pac', 'sa', 'rbwt', 'rpac', 'rsa' are all bwa index.
            stidx is stampy index.
            sthash is stampy hash.
        2012.2.23
            add two suffixes, stidx (stampy index) and sthash (stampy hash)
        2011-11-11
            if needAffiliatedFiles,
                all other files, with suffix in affiliateFilenameSuffixLs, will be registered (symlinked or copied) as well.
        """
        returnData = PassingData(refFastaFList = [], needBWARefIndexJob=False, needSAMtoolsFastaIndexJob=False, \
                                needPicardFastaDictJob=False, needStampyRefIndexJob=False, needBlastMakeDBJob=False,\
                                refPicardFastaDictF=None, refSAMtoolsFastaIndexF=None)
        missingSuffixSet = set()	#2013.3.20
        
        if registerAffiliateFiles:
            refFastaF = File(os.path.join(folderName, os.path.basename(refFastaFname)))	#use relative path, otherwise, it'll go to absolute path
            # Add it into replica only when needed.
            refFastaF.addPFN(PFN("file://" + refFastaFname, self.input_site_handler))
            if not self.hasFile(refFastaF):	#2013.08.12
                self.addFile(refFastaF)
            returnData.refFastaFList.append(refFastaF)
            # If it's not needed, assume the index is done and all relevant files are in absolute path.
            # and no replica transfer
            
            #add extra affiliated files
            suffix2PathToFileLs = {}
            if addPicardDictFile:	#2012.5.23
                picardDictSuffix = 'dict'
                pathToFile = '%s.%s'%(os.path.splitext(refFastaFname)[0], picardDictSuffix)	#remove ".fasta" from refFastaFname
                if checkAffiliateFileExistence and not os.path.isfile(pathToFile):
                    sys.stderr.write("Warning: %s don't exist or not a file on file system. skip registration.\n"%(pathToFile))
                    missingSuffixSet.add(picardDictSuffix)
                    #suffix2PathToFileLs.append(pathToFile)
                else:
                    suffix2PathToFileLs[picardDictSuffix] = pathToFile
            for suffix in affiliateFilenameSuffixLs:
                pathToFile = '%s.%s'%(refFastaFname, suffix)
                if checkAffiliateFileExistence and not os.path.isfile(pathToFile):
                    sys.stderr.write("Warning: %s don't exist or not a file on file system. skip registration.\n"%(pathToFile))
                    missingSuffixSet.add(suffix)
                    continue
                suffix2PathToFileLs[suffix]= pathToFile
            for suffix, pathToFile in suffix2PathToFileLs.items():
                if checkAffiliateFileExistence and not os.path.isfile(pathToFile):
                    sys.stderr.write("Warning: %s don't exist or not a file on file system. skip registration.\n"%(pathToFile))
                    continue
                affiliateF = File(os.path.join(folderName, os.path.basename(pathToFile)))
                #use relative path, otherwise, it'll go to absolute path
                affiliateF.addPFN(PFN("file://" + pathToFile, self.input_site_handler))
                if not self.hasFile(affiliateF):	#2013.08.12
                    self.addFile(affiliateF)
                returnData.refFastaFList.append(affiliateF)
                
                if suffix=='dict':	#2013.3.26
                    returnData.refPicardFastaDictF = affiliateF
                elif suffix=='fai':
                    returnData.refSAMtoolsFastaIndexF = affiliateF
        else:
            refFastaF = File(os.path.join(folderName, os.path.basename(refFastaFname)))
            returnData.refFastaFList.append(refFastaF)
        if 'bwt' in missingSuffixSet or 'pac' in missingSuffixSet:
            returnData.needBWARefIndexJob = True
        if 'fai' in missingSuffixSet:
            returnData.needSAMtoolsFastaIndexJob = True
            returnData.needPicardFastaDictJob = True
        if 'stidx' in missingSuffixSet or 'sthash' in missingSuffixSet:
            returnData.needStampyRefIndexJob = True
        if 'dict' in missingSuffixSet:
            returnData.needPicardFastaDictJob = True
        if 'nin' in missingSuffixSet or 'nhr' in missingSuffixSet or 'nsq' in missingSuffixSet:
            returnData.needBlastMakeDBJob = True
        return returnData

    def addPlinkJob(self, executable=None, inputFileList=None, parentPlinkJob=None,\
                tpedFile=None, tfamFile=None,\
                pedFile=None, famFile=None, mapFile=None, bedFile=None, bimFile=None,\
                inputFnamePrefix=None, inputOption='--file', \
                outputFnamePrefix=None, outputOption='--out',\
                makeBED=False, calculateMendelError=False, checkSex=False, \
                LDPruneWindowSize=100, LDPruneWindowShiftSize=5, LDPruneByPairwiseR2=False, LDPruneMinR2=0.1,\
                LDPruneByRegression=False, LDPruneMinVarianceInflationFactor=2,\
                estimatePairwiseGenomeWideIBD=False, estimatePairwiseGenomeWideIBDFreqFile=None, \
                extractSNPFile=None, recodeOutput=False, recodeTransposeOutput=False, estimateAlleFrequency=False, \
                mergeListFile=None,\
                parentJobLs=None, extraDependentInputLs=None, transferOutput=False, \
                extraArguments=None, extraArgumentList=None, extraOutputLs =None, \
                job_max_memory=2000, **keywords):
        """
        i.e.

            bedFnamePrefix = os.path.join(topOutputDir, '%s_bed'%(commonPrefix))
            convertSingleTPED2BEDJob = self.addPlinkJob(executable=self.plink, inputFileList=[],
                                tpedFile=modifyTPEDJob.output, tfamFile=tfamJob.tfamFile,\
                outputFnamePrefix=bedFnamePrefix, outputOption='--out',\
                makeBED=True, \
                extraDependentInputLs=None, transferOutput=transferOutput, \
                extraArguments=None, job_max_memory=2000,\
                parentJobLs = convertSingleTPED2BEDParentJobLs)


            convertMergedTPED2BEDJob = self.addPlinkJob(executable=self.plink, inputFileList=[tpedFileMergeJob.output, tfamJob.tfamFile], \
                            inputFnamePrefix=mergedPlinkFnamePrefix, inputOption='--tfile', \
                outputFnamePrefix=mergedPlinkBEDFnamePrefix, outputOption='--out',\
                makeBED=True, \
                extraDependentInputLs=None, transferOutput=transferOutput, \
                extraArguments=None, job_max_memory=2000, parentJobLs=[mergedOutputDirJob, tpedFileMergeJob, tfamJob])

            mendelFnamePrefix = os.path.join(setupData.mapDirJob.output, '%s'%(commonPrefix))
            if inputJob.output.name[-4:]=='tped':	#2013.07.25 make sure addPlinkJob could get the right tfamFile
                inputJob.tfamFile = tfamJob.tfamFile
            plinkMendelJob = self.addPlinkJob(executable=self.plink, \
                    parentPlinkJob=inputJob,\
                    outputFnamePrefix=mendelFnamePrefix, outputOption='--out',\
                    calculateMendelError=True, \
                    extraDependentInputLs=None, transferOutput=transferOneContigPlinkOutput, \
                    extraArguments=None, job_max_memory=2000,\
                    parentJobLs =[setupData.mapDirJob, tfamJob]+ jobData.jobLs)

        for plink mendel, LD-prune and other jobs, add extraArguments="--allow-no-sex" to include individuals without sex

        2013.07.25 added parentPlinkJob (returned from this function), and parse input from that job
        2013.07.24 added argument recodeTransposeOutput (--recode --transpose)
        2012.8.28
            add argument
                estimateAlleFrequency, estimate frequency of input file. "--nonfounders" could be added as well.
                estimatePairwiseGenomeWideIBDFreqFile, is the file from which IBD check could draw frequency (rather than estimate from founders)

        2012.8.9
            inputFileList is a list of pegasus Files (.ped, .fam, or .tped, .tfam, etc.) or could be supplied individually.

            inputOption could be, "--file" for .ped .map ; "--tfile" for .tped, .tfam; or '--bfile' for .bed, .fam, .bim

            if extractSNPFile or mergeListFile is given, either recodeOutput or makeBED have to be on. otherwise, no output.
            http://pngu.mgh.harvard.edu/~purcell/plink/index.shtml



        """
        if extraDependentInputLs is None:
            extraDependentInputLs = []
        if inputFileList:
            extraDependentInputLs.extend(inputFileList)

        if extraArgumentList is None:
            extraArgumentList = []
        if extraOutputLs is None:
            extraOutputLs = []
        key2ObjectForJob = {}

        #2013.07.25
        if parentPlinkJob:
            if bedFile is None:
                bedFile = getattr(parentPlinkJob, 'bedFile', None)
            if famFile is None:
                famFile = getattr(parentPlinkJob, 'famFile', None)
            if bimFile is None:
                bimFile = getattr(parentPlinkJob, 'bimFile', None)
            if tpedFile is None:
                tpedFile = getattr(parentPlinkJob, 'tpedFile', None)
            if tfamFile is None:
                tfamFile = getattr(parentPlinkJob, 'tfamFile', None)
            if mapFile is None:
                mapFile = getattr(parentPlinkJob, 'mapFile', None)
            if pedFile is None:
                pedFile = getattr(parentPlinkJob, 'pedFile', None)
            if famFile is None:
                famFile = getattr(parentPlinkJob, 'famFile', None)

        if inputOption and inputFnamePrefix:
            extraArgumentList.extend([inputOption, inputFnamePrefix])
        if tpedFile:
            extraDependentInputLs.append(tpedFile)
            extraArgumentList.extend(["--tped", tpedFile])
        if tfamFile:
            extraDependentInputLs.append(tfamFile)
            extraArgumentList.extend(["--tfam", tfamFile])
        if pedFile:
            extraDependentInputLs.append(pedFile)
            extraArgumentList.extend(["--ped", pedFile])
        if famFile:
            extraDependentInputLs.append(famFile)
            extraArgumentList.extend(["--fam", famFile])
        if mapFile:
            extraDependentInputLs.append(mapFile)
            extraArgumentList.extend(["--map", mapFile])
        if bedFile:
            extraDependentInputLs.append(bedFile)
            extraArgumentList.extend(["--bed", bedFile])
        if bimFile:
            extraDependentInputLs.append(bimFile)
            extraArgumentList.extend(["--bim", bimFile])

        if outputFnamePrefix and outputOption:
            extraArgumentList.extend([outputOption, outputFnamePrefix])
        else:
            outputFnamePrefix = 'plink'


        suffixAndNameTupleList = []	# a list of tuples , in each tuple, 1st element is the suffix. 2nd element is the proper name of the suffix.
            #job.$nameFile will be the way to access the file.
            #if 2nd element (name) is missing, suffix[1:].replace('.', '_') is the name (dot replaced by _)
        if makeBED:
            extraArgumentList.append('--make-bed')
            suffixAndNameTupleList.extend([['.bed',], ('.fam',), ['.bim',]])		#, binary map file, is excluded for now
        if calculateMendelError:
            extraArgumentList.append('--mendel')
            suffixAndNameTupleList.extend([('.mendel',), ('.imendel',), ('.fmendel',), ('.lmendel',)])
            #its output is not tab-delimited. rather it's space (multi) delimited.
        if checkSex:
            extraArgumentList.append('--check-sex')
            suffixAndNameTupleList.extend([('.sexcheck',), ('.hh', )])	#.sexcheck file is accessible as job.sexcheckFile.
                #.hh is heterozygous haplotype genotypes
        if LDPruneByPairwiseR2:
            extraArgumentList.append('--indep-pairwise %s %s %s'%(LDPruneWindowSize, LDPruneWindowShiftSize, LDPruneMinR2))
            suffixAndNameTupleList.extend([('.prune.in',), ('.prune.out',)])	#".prune.in" is accessible as job.prune_inFile
        if LDPruneByRegression:
            extraArgumentList.append('--indep %s %s %s'%(LDPruneWindowSize, LDPruneWindowShiftSize, LDPruneMinVarianceInflationFactor))
            suffixAndNameTupleList.extend([('.prune.in',), ('.prune.out',)])	#".prune.in" is accessible as job.prune_inFile
        if estimatePairwiseGenomeWideIBD:
            extraArgumentList.append('--genome')
            suffixAndNameTupleList.extend([('.genome',)])	#.genome is accessible as job.genomeFile
            if estimatePairwiseGenomeWideIBDFreqFile:	#2012.8.28
                extraArgumentList.extend(['--read-freq', estimatePairwiseGenomeWideIBDFreqFile])
                extraDependentInputLs.append(estimatePairwiseGenomeWideIBDFreqFile)
        if extractSNPFile:
            extraArgumentList.extend(['--extract', extractSNPFile])
            extraDependentInputLs.append(extractSNPFile)
        if recodeOutput:
            extraArgumentList.extend(['--recode',])
            suffixAndNameTupleList.extend([('.ped',), ('.map',)])
        if recodeTransposeOutput:
            extraArgumentList.extend(['--recode', "--transpose"])
            suffixAndNameTupleList.extend([('.tped',), ('.tfam',)])
        if estimateAlleFrequency:	#2012.8.28
            extraArgumentList.append('--freq')
            suffixAndNameTupleList.extend([('.frq',)])

        if mergeListFile:
            extraArgumentList.extend(['--merge-list', mergeListFile])
            extraDependentInputLs.append(mergeListFile)
        if extraArguments:
            extraArgumentList.append(extraArguments)


        self.setupMoreOutputAccordingToSuffixAndNameTupleList(outputFnamePrefix=outputFnamePrefix, \
            suffixAndNameTupleList=suffixAndNameTupleList, \
            extraOutputLs=extraOutputLs, key2ObjectForJob=key2ObjectForJob)
        #2013.07.24 add it in the end
        logFile = File('%s.log'%(outputFnamePrefix))	#2012.8.10 left in the folder dying
        extraOutputLs.append(logFile)

        job= self.addGenericJob(executable=executable, inputFile=None, outputFile=None, \
                parentJobLs=parentJobLs, extraDependentInputLs=extraDependentInputLs, \
                extraOutputLs=extraOutputLs,\
                transferOutput=transferOutput, \
                extraArgumentList=extraArgumentList, key2ObjectForJob=key2ObjectForJob, job_max_memory=job_max_memory, **keywords)
        return job

    def registerExecutables(self):
        """
        """
        if not workflow:
            workflow = self
        ParentClass.registerExecutables(self)

        namespace = self.namespace
        version = self.version
        operatingSystem = self.operatingSystem
        architecture = self.architecture
        cluster_size = self.cluster_size
        site_handler = self.site_handler


        #2013.11.22	#2013.06.25 register tabix
        self.addExecutableFromPath(
            path=self.tabixPath, name='tabix', clusterSizeMultiplier=5)

        #2013.11.22 2011.12.21	for OutputVCFSiteStat.py
        self.addExecutableFromPath(
            path=os.path.join(self.pymodulePath, "mapper/extractor/tabixRetrieve.sh"),
            name='tabixRetrieve', clusterSizeMultiplier=1)

        #2013.11.22 moved from palos/polymorphism/FindNewRefCoordinatesGivenVCFFolderWorkflow.py
        self.addExecutableFromPath(path=os.path.join(self.pymodulePath, \
            "polymorphism/mapper/LiftOverVCFBasedOnCoordinateMap.py"), \
            name='LiftOverVCFBasedOnCoordinateMap', clusterSizeMultiplier=1)

        self.addExecutableFromPath(path=os.path.join(self.pymodulePath, \
            "polymorphism/qc/CalculateLociAndGenomeCoveredAtEachSwitchFrequencyThreshold.py"), \
            name='CalculateLociAndGenomeCoveredAtEachSwitchFrequencyThreshold', clusterSizeMultiplier=0.01)

        self.addExecutableFromPath(path=os.path.join(self.pymodulePath, \
                "mapper/extractor/ExtractFlankingSequenceForVCFLoci.py"), \
            name='ExtractFlankingSequenceForVCFLoci', clusterSizeMultiplier=2)

        self.addExecutableFromPath(path=os.path.join(self.pymodulePath, \
            "polymorphism/mapper/FindSNPPositionOnNewRefFromFlankingBlastOutput.py"), \
            name='FindSNPPositionOnNewRefFromFlankingBlastOutput', clusterSizeMultiplier=2)

        self.addExecutableFromPath(path=os.path.join(self.pymodulePath, \
            "polymorphism/mapper/FindSNPPositionOnNewRefFromFlankingBWAOutput.py"), \
            name='FindSNPPositionOnNewRefFromFlankingBWAOutput', clusterSizeMultiplier=1)


        #2013.08.28
        self.addExecutableFromPath(path=os.path.join(self.pymodulePath, 'Genome/OutputGenomeAnnotation.py'), \
            name='OutputGenomeAnnotation', clusterSizeMultiplier=0.01)
        #2013.07.31
        self.addExecutableFromPath(path=os.path.join(self.pymodulePath, 'statistics/GenomeMovingAverageStatistics.py'), \
            name='GenomeMovingAverageStatistics', clusterSizeMultiplier=0.1)
        #2013.08.23
        self.addExecutableFromPath(path=os.path.join(self.pymodulePath, 'reducer/ReduceSameChromosomeAlignmentDepthFiles'), \
            name='ReduceSameChromosomeAlignmentDepthFiles', clusterSizeMultiplier=0.5)

        #executableList = []
        executableClusterSizeMultiplierList = []	#2012.8.7 each cell is a tuple of (executable, clusterSizeMultiplier (0 if u do not need clustering)
        #noClusteringExecutableSet = set()	#2012.8.2 you don't want to cluster for some jobs.

        PlotLD = Executable(namespace=namespace, name="PlotLD", version=version, os=operatingSystem, arch=architecture, installed=True)
        PlotLD.addPFN(PFN("file://" +  os.path.join(self.pymodulePath, "plot/PlotLD.py"), site_handler))
        executableClusterSizeMultiplierList.append((PlotLD, 0))

        #2012.8.13
        OutputVCFSiteGap = Executable(namespace=namespace, name="OutputVCFSiteGap", \
                            version=version, os=operatingSystem, arch=architecture, installed=True)
        OutputVCFSiteGap.addPFN(PFN("file://" + os.path.join(self.pymodulePath, "mapper/computer/OutputVCFSiteGap.py"), site_handler))
        executableClusterSizeMultiplierList.append((OutputVCFSiteGap, 1))

        
        ConvertBjarniSNPFormat2Yu = Executable(namespace=namespace, name="ConvertBjarniSNPFormat2Yu", \
                                            version=version, os=operatingSystem, arch=architecture, installed=True)
        ConvertBjarniSNPFormat2Yu.addPFN(PFN("file://" +  os.path.join(self.variationSrcPath, "io/ConvertBjarniSNPFormat2Yu.py"), site_handler))
        executableClusterSizeMultiplierList.append((ConvertBjarniSNPFormat2Yu, 1))

        ConvertVCF2BjarniFormat = Executable(namespace=namespace, name="ConvertVCF2BjarniFormat", \
                                            version=version, os=operatingSystem, arch=architecture, installed=True)
        ConvertVCF2BjarniFormat.addPFN(PFN("file://" +  os.path.join(self.pymodulePath, "mapper/converter/ConvertVCF2BjarniFormat.py"), site_handler))
        executableClusterSizeMultiplierList.append((ConvertVCF2BjarniFormat, 1))


        ConvertYuSNPFormat2Bjarni = Executable(namespace=namespace, name="ConvertYuSNPFormat2Bjarni", \
                                            version=version, os=operatingSystem, arch=architecture, installed=True)
        ConvertYuSNPFormat2Bjarni.addPFN(PFN("file://" +  os.path.join(self.variationSrcPath, "io/ConvertYuSNPFormat2Bjarni.py"), site_handler))
        executableClusterSizeMultiplierList.append((ConvertYuSNPFormat2Bjarni, 1))

        ConvertYuSNPFormat2EigenStrat = Executable(namespace=namespace, name="ConvertYuSNPFormat2EigenStrat", \
                                            version=version, os=operatingSystem, arch=architecture, installed=True)
        ConvertYuSNPFormat2EigenStrat.addPFN(PFN("file://" +  os.path.join(self.variationSrcPath, "io/ConvertYuSNPFormat2EigenStrat.py"), site_handler))
        executableClusterSizeMultiplierList.append((ConvertYuSNPFormat2EigenStrat, 1))

        ConvertYuSNPFormat2TPED_TFAM = Executable(namespace=namespace, name="ConvertYuSNPFormat2TPED_TFAM", \
                                            version=version, os=operatingSystem, arch=architecture, installed=True)
        ConvertYuSNPFormat2TPED_TFAM.addPFN(PFN("file://" +  os.path.join(self.variationSrcPath, "io/ConvertYuSNPFormat2TPED_TFAM.py"), site_handler))
        executableClusterSizeMultiplierList.append((ConvertYuSNPFormat2TPED_TFAM, 1))

        CalculatePairwiseDistanceOutOfSNPXStrainMatrix = Executable(namespace=namespace, name="CalculatePairwiseDistanceOutOfSNPXStrainMatrix", \
                            version=version, os=operatingSystem, arch=architecture, installed=True)
        CalculatePairwiseDistanceOutOfSNPXStrainMatrix.addPFN(PFN("file://" + \
                        os.path.join(self.pymodulePath, "mapper/CalculatePairwiseDistanceOutOfSNPXStrainMatrix.py"), \
                        site_handler))
        executableClusterSizeMultiplierList.append((CalculatePairwiseDistanceOutOfSNPXStrainMatrix, 0.5))

        #2013.2.3 use samtools to extract consensus from bam files
        ExtractConsensusSequenceFromAlignment = Executable(namespace=namespace, name="ExtractConsensusSequenceFromAlignment", version=version, \
                        os=operatingSystem, arch=architecture, installed=True)
        ExtractConsensusSequenceFromAlignment.addPFN(PFN("file://" + \
            os.path.join(self.pymodulePath, "mapper/alignment/ExtractConsensusSequenceFromAlignment.sh"), site_handler))
        executableClusterSizeMultiplierList.append((ExtractConsensusSequenceFromAlignment, 1))

        #2013.2.4, wrapper around psmc's splitfa, a program that splits fasta files
        splitfa = Executable(namespace=namespace, name="splitfa", version=version, \
                        os=operatingSystem, arch=architecture, installed=True)
        splitfa.addPFN(PFN("file://" + os.path.join(self.pymodulePath, "mapper/splitter/splitfa.sh"), site_handler))
        executableClusterSizeMultiplierList.append((splitfa, 1))

        self.setExecutablesClusterSize(executableClusterSizeMultiplierList, defaultClusterSize=self.cluster_size)

        self.addExecutableFromPath(path=os.path.join(self.pymodulePath, "plot/PlotVCFtoolsStat.py"), \
                                        name='PlotVCFtoolsStat', clusterSizeMultiplier=0)
        
        #2013.07.19
        self.addExecutableFromPath(path=os.path.join(self.pymodulePath, 'mapper/modifier/AppendExtraPedigreeIndividualsToTPED.py'), \
                                        name='AppendExtraPedigreeIndividualsToTPED', clusterSizeMultiplier=1)

        self.addExecutableFromPath(path=os.path.join(self.pymodulePath, 'mapper/converter/ConvertMSOutput2FASTQ.py'), \
                                        name='ConvertMSOutput2FASTQ', clusterSizeMultiplier=1)

        self.addExecutableFromPath(path=os.path.join(self.pymodulePath, 'mapper/extractor/SelectChromosomeSequences.py'), \
                                        name='SelectChromosomeSequences', clusterSizeMultiplier=0.5)

        #2013.2.11 moved from vervet/src/reduce to pymodule/reducer
        self.addExecutableFromPath(path=os.path.join(self.pymodulePath, 'reducer/MergeGenotypeMatrix.py'), \
                                        name='MergeGenotypeMatrix', clusterSizeMultiplier=0.2)
        self.addExecutableFromPath(path=os.path.join(self.pymodulePath, 'plot/PlotGenomeWideData.py'), \
                                        name='PlotGenomeWideData', clusterSizeMultiplier=1)

        self.bgzipExecutableFile = self.registerOneExecutableAsFile(path=os.path.expanduser("~/bin/bgzip"))	#2013.11.22

if __name__ == '__main__':
    main_class = AbstractBioinfoWorkflow
    po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
    instance = main_class(**po.long_option2value)
    instance.run()