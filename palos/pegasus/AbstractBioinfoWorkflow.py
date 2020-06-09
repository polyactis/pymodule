#!/usr/bin/env python3
"""
2020/01/29
    an abstract class for pegasus workflows that work on bioinformatic data
"""
import sys, os, math
import logging
from palos import ProcessOptions, getListOutOfStr, PassingData, utils
from pegaflow.DAX3 import Executable, File, PFN, Link, Job
from . AbstractWorkflow import AbstractWorkflow

ParentClass = AbstractWorkflow
class AbstractBioinfoWorkflow(ParentClass):
    __doc__ = __doc__
    option_default_dict = ParentClass.option_default_dict.copy()
    option_default_dict.update({
        ("plinkPath", 1, ): ["bin/plink", '', 1, 
            'path to plink, http://pngu.mgh.harvard.edu/~purcell/plink/index.shtml'],
        })
    def __init__(self,
        input_path=None,
        inputSuffixList=None, 
        pegasusFolderName='folder', output_path=None,
        
        tmpDir='/tmp/', max_walltime=4320,
        home_path=None, javaPath=None, 
        pymodulePath="src/pymodule", thisModulePath=None,
        jvmVirtualByPhysicalMemoryRatio=1.2,

        site_handler='condor', input_site_handler='condor', cluster_size=30,
        needSSHDBTunnel=False, commit=False,
        debug=False, report=False):
        """
        20200129
        """
        ParentClass.__init__(self,
            input_path=input_path,
            inputSuffixList=inputSuffixList,
            pegasusFolderName=pegasusFolderName,
            output_path=output_path,
            
            tmpDir=tmpDir, max_walltime=max_walltime, 
            home_path=home_path,
            javaPath=javaPath,
            pymodulePath=pymodulePath, thisModulePath=thisModulePath,
            jvmVirtualByPhysicalMemoryRatio=jvmVirtualByPhysicalMemoryRatio,
            
            site_handler=site_handler,
            input_site_handler=input_site_handler,
            cluster_size=cluster_size,

            needSSHDBTunnel=needSSHDBTunnel, commit=commit,
            debug=debug, report=report)
    
    def registerBlastNucleotideDatabaseFile(self, ntDatabaseFname=None,  folderName=""):
        """
        moved from BlastWorkflow.py
        """
        return self.registerRefFastaFile(refFastaFname=ntDatabaseFname,
            registerAffiliateFiles=True,
            checkAffiliateFileExistence=True, addPicardDictFile=False,
            affiliateFilenameSuffixLs=['nin', 'nhr', 'nsq'],\
            folderName=folderName)

    def registerRefFastaFile(self, refFastaFname=None,
        registerAffiliateFiles=True,
        checkAffiliateFileExistence=True, addPicardDictFile=True,\
        affiliateFilenameSuffixLs=\
            ['fai', 'amb', 'ann', 'bwt', 'pac', 'sa', 'rbwt', 'rpac', 'rsa', \
            'stidx', 'sthash'],
        folderName="reference"):
        """
        suffix here doesn't include ".".
        
        deduce needBWARefIndexJob, needSAMtoolsFastaIndexJob, needPicardFastaDictJob,
            needStampyRefIndexJob from missing suffixes.
        argument folderName
        argument "addPicardDictFile" to offer user option to exclude this file 
            (i.e. in registerBlastNucleotideDatabaseFile).
        
        .dict is output by picard, also required for GATK.
            fai is via "samtools faidx" (index reference). also required for GATK
            amb', 'ann', 'bwt', 'pac', 'sa', 'rbwt', 'rpac', 'rsa' are all bwa index.
            stidx is stampy index.
            sthash is stampy hash.
        stidx (stampy index) and sthash (stampy hash)
        2011-11-11
            if needAffiliatedFiles,
                all other files, with suffix in affiliateFilenameSuffixLs,
                 will be registered (symlinked or copied) as well.
        """
        returnData = PassingData(refFastaFList = [],
            needBWARefIndexJob=False,
            needSAMtoolsFastaIndexJob=False,
            needPicardFastaDictJob=False,
            needStampyRefIndexJob=False,
            needBlastMakeDBJob=False,
            refPicardFastaDictF=None, refSAMtoolsFastaIndexF=None)
        missingSuffixSet = set()
        
        if registerAffiliateFiles:
            refFastaF = File(os.path.join(folderName, os.path.basename(refFastaFname)))
            #use relative path, otherwise, it'll go to absolute path
            # Add it into replica only when needed.
            refFastaF.addPFN(PFN("file://" + refFastaFname, self.input_site_handler))
            if not self.hasFile(refFastaF):
                self.addFile(refFastaF)
            returnData.refFastaFList.append(refFastaF)
            # If it's not needed, assume the index is done and all relevant files
            #  are in absolute path. and no replica transfer
            
            #add extra affiliated files
            suffix2PathToFileLs = {}
            if addPicardDictFile:
                picardDictSuffix = 'dict'
                pathToFile = '%s.%s'%(os.path.splitext(refFastaFname)[0], \
                    picardDictSuffix)
                #remove ".fasta" from refFastaFname
                if checkAffiliateFileExistence and not os.path.isfile(pathToFile):
                    logging.warn(f"{pathToFile} don't exist or not a file. "
                        "Skip registration.")
                    missingSuffixSet.add(picardDictSuffix)
                    #suffix2PathToFileLs.append(pathToFile)
                else:
                    suffix2PathToFileLs[picardDictSuffix] = pathToFile
            for suffix in affiliateFilenameSuffixLs:
                pathToFile = '%s.%s'%(refFastaFname, suffix)
                if checkAffiliateFileExistence and not os.path.isfile(pathToFile):
                    logging.warn(f"{pathToFile} don't exist or not a file. "
                        "Skip registration.")
                    missingSuffixSet.add(suffix)
                    continue
                suffix2PathToFileLs[suffix]= pathToFile
            for suffix, pathToFile in suffix2PathToFileLs.items():
                if checkAffiliateFileExistence and not os.path.isfile(pathToFile):
                    logging.warn(f"{pathToFile} don't exist or not a file. "
                        "Skip registration.")
                    continue
                affiliateF = File(os.path.join(folderName, os.path.basename(pathToFile)))
                #use relative path, otherwise, it'll go to absolute path
                affiliateF.addPFN(PFN("file://" + pathToFile, self.input_site_handler))
                if not self.hasFile(affiliateF):
                    self.addFile(affiliateF)
                returnData.refFastaFList.append(affiliateF)
                
                if suffix=='dict':
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
        if 'nin' in missingSuffixSet or 'nhr' in missingSuffixSet or \
            'nsq' in missingSuffixSet:
            returnData.needBlastMakeDBJob = True
        return returnData

    def registerExecutables(self):
        """
        """
        ParentClass.registerExecutables(self)
        #for OutputVCFSiteStat.py
        self.registerOneExecutable(
            path=os.path.join(self.pymodulePath, "mapper/extractor/tabixRetrieve.sh"),
            name='tabixRetrieve', clusterSizeMultiplier=1)

        # from palos/polymorphism/FindNewRefCoordinatesGivenVCFFolder.py
        self.registerOneExecutable(path=os.path.join(self.pymodulePath, \
            "polymorphism/mapper/LiftOverVCFBasedOnCoordinateMap.py"), \
            name='LiftOverVCFBasedOnCoordinateMap', clusterSizeMultiplier=1)

        self.registerOneExecutable(path=os.path.join(self.pymodulePath, \
            "polymorphism/qc/"+\
            "CalculateLociAndGenomeCoveredAtEachSwitchFrequencyThreshold.py"),
            name='CalculateLociAndGenomeCoveredAtEachSwitchFrequencyThreshold',
            clusterSizeMultiplier=0.01)

        self.registerOneExecutable(path=os.path.join(self.pymodulePath, \
                "mapper/extractor/ExtractFlankingSequenceForVCFLoci.py"), \
            name='ExtractFlankingSequenceForVCFLoci', clusterSizeMultiplier=2)

        self.registerOneExecutable(path=os.path.join(self.pymodulePath, \
            "polymorphism/mapper/FindSNPPositionOnNewRefFromFlankingBlastOutput.py"), \
            name='FindSNPPositionOnNewRefFromFlankingBlastOutput',
            clusterSizeMultiplier=2)

        self.registerOneExecutable(path=os.path.join(self.pymodulePath, \
            "polymorphism/mapper/FindSNPPositionOnNewRefFromFlankingBWAOutput.py"), \
            name='FindSNPPositionOnNewRefFromFlankingBWAOutput',
            clusterSizeMultiplier=1)
        self.registerOneExecutable(path=os.path.join(self.pymodulePath, \
            'db/export/OutputGenomeAnnotation.py'), \
            name='OutputGenomeAnnotation', clusterSizeMultiplier=0.01)
        self.registerOneExecutable(path=os.path.join(self.pymodulePath, \
                'statistics/GenomeMovingAverageStatistics.py'), \
            name='GenomeMovingAverageStatistics', clusterSizeMultiplier=0.1)
        self.registerOneExecutable(path=os.path.join(self.pymodulePath, \
                'mapper/computer/OutputVCFSiteGap.py'), \
            name='OutputVCFSiteGap', clusterSizeMultiplier=1)
        self.registerOneExecutable(path=os.path.join(self.pymodulePath, \
                'plot/PlotLD.py'), \
            name='PlotLD', clusterSizeMultiplier=0)
        
        self.registerOneExecutable(path=os.path.join(self.pymodulePath, 
            'mapper/converter/ConvertBjarniSNPFormat2Yu.py'), \
            name='ConvertBjarniSNPFormat2Yu', clusterSizeMultiplier=1)
        self.registerOneExecutable(path=os.path.join(self.pymodulePath, 
            'mapper/converter/ConvertVCF2BjarniFormat.py'), \
            name='ConvertVCF2BjarniFormat', clusterSizeMultiplier=1)
        self.registerOneExecutable(path=os.path.join(self.pymodulePath, 
            'mapper/converter/ConvertYuSNPFormat2Bjarni.py'), \
            name='ConvertYuSNPFormat2Bjarni', clusterSizeMultiplier=1)

        self.registerOneExecutable(path=os.path.join(self.pymodulePath, 
            'mapper/converter/ConvertYuSNPFormat2EigenStrat.py'), \
            name='ConvertYuSNPFormat2EigenStrat', clusterSizeMultiplier=1)
        self.registerOneExecutable(path=os.path.join(self.pymodulePath, 
            'mapper/converter/ConvertYuSNPFormat2TPED_TFAM.py'), \
            name='ConvertYuSNPFormat2TPED_TFAM', clusterSizeMultiplier=1)
        self.registerOneExecutable(path=os.path.join(self.pymodulePath, 
            'mapper/CalculatePairwiseDistanceOutOfSNPXStrainMatrix.py'),
            name='CalculatePairwiseDistanceOutOfSNPXStrainMatrix',
            clusterSizeMultiplier=0.5)
        #use samtools to extract consensus from bam files
        self.registerOneExecutable(path=os.path.join(self.pymodulePath, 
            'mapper/alignment/ExtractConsensusSequenceFromAlignment.sh'),
            name='ExtractConsensusSequenceFromAlignment',
            clusterSizeMultiplier=0.5)
        #wrapper around psmc's splitfa,
        #  a program that splits fasta files
        self.registerOneExecutable(path=os.path.join(self.pymodulePath,
                "mapper/splitter/splitfa.sh"), \
            name='splitfa', clusterSizeMultiplier=1)
        self.registerOneExecutable(path=os.path.join(self.pymodulePath,
                'mapper/modifier/AppendExtraPedigreeIndividualsToTPED.py'),
            name='AppendExtraPedigreeIndividualsToTPED', clusterSizeMultiplier=1)
        self.registerOneExecutable(path=os.path.join(self.pymodulePath,
                'mapper/converter/ConvertMSOutput2FASTQ.py'),
            name='ConvertMSOutput2FASTQ', clusterSizeMultiplier=1)
        self.registerOneExecutable(path=os.path.join(self.pymodulePath,
                'mapper/extractor/SelectChromosomeSequences.py'),
            name='SelectChromosomeSequences', clusterSizeMultiplier=0.5)
        self.registerOneExecutable(path=os.path.join(self.pymodulePath,
                'plot/PlotGenomeWideData.py'),
            name='PlotGenomeWideData', clusterSizeMultiplier=1)

if __name__ == '__main__':
    from argparse import ArgumentParser
    ap = ArgumentParser()
    ap.add_argument("-i", "--input_path", type=str, required=True,
        help="the path to the input folder or file.")
    ap.add_argument("--inputSuffixList", type=str,
        help='Coma-separated list of input file suffices. Used to exclude input files.'
        'If None, no exclusion. The dot is part of the suffix, .tsv not tsv.'
        'Common zip suffices (.gz, .bz2, .zip, .bz) will be ignored in obtaining the suffix.')
    ap.add_argument("-F", "--pegasusFolderName", type=str,
        help='The path relative to the workflow running root. '
        'This folder will contain pegasus input & output. '
        'It will be created during the pegasus staging process. '
        'It is useful to separate multiple sub-workflows. '
        'If empty or None, everything is in the pegasus root.')
    
    ap.add_argument("-o", "--output_path", type=str, required=True,
        help="The path to the output file that will contain the Pegasus DAG.")
    
    ap.add_argument("--home_path", type=str,
        help="Path to your home folder. Default is ~.")
    ap.add_argument("--javaPath", type=str, default='bin/java',
        help="Path to java. "
        "If relative path, home folder is inserted in the front.")
    ap.add_argument("--pymodulePath", type=str, default="src/pymodule",
        help="Path to the pymodule code folder. "
        "If relative path, home folder is inserted in the front.")
    ap.add_argument("--thisModulePath", type=str,
        help="Path to your own code repo"
        "If relative path, home folder is inserted in the front.")
    ap.add_argument("--plinkPath", type=str, default='bin/plink',
        help="Path to plink, "
        "http://pngu.mgh.harvard.edu/~purcell/plink/index.shtml "
        "If relative path, home folder is inserted in the front.")
    
    ap.add_argument("--tmpDir", type=str, default='/tmp/',
        help='Default: %(default)s. '
        'A local folder for some jobs (MarkDup) to store temp data. '
        '/tmp/ can be too small for high-coverage sequencing.')
    ap.add_argument("--max_walltime", type=int, default=4320,
        help='Default: %(default)s. '
        'Maximum wall time for any job, in minutes. 4320=3 days. '
        'Used in addGenericJob(). Most clusters have upper limit for runtime.')
    ap.add_argument("--jvmVirtualByPhysicalMemoryRatio", type=float,
        default=1.2,
        help='Default: %(default)s. '
        'If a job virtual memory (~1.2X of JVM resident memory) exceeds request, '
        "it will be killed on some clusters. "
        "This will make sure your job has enough memory.")
    
    ap.add_argument("-l", "--site_handler", type=str, required=True,
        help="The name of the computing site where the jobs run and executables are stored. "
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
        "This number decides how many of pegasus jobs should be clustered into one job. "
        "Good if your workflow contains many quick jobs. "
        "It will reduce Pegasus monitor I/O.")
    
    ap.add_argument("--debug", action='store_true',
        help='Toggle debug mode.')
    ap.add_argument("--report", action='store_true',
        help="Toggle verbose mode. Default: %(default)s.")
    ap.add_argument("--needSSHDBTunnel", action='store_true',
        help="If all DB-interacting jobs need a ssh tunnel to "
        "access a database that is inaccessible to computing nodes.")
    args = ap.parse_args()
    instance = AbstractBioinfoWorkflow(
        input_path=args.input_path,
        inputSuffixList=args.inputSuffixList, 
        pegasusFolderName=args.pegasusFolderName,
        output_path=args.output_path,

        home_path=args.home_path,
        javaPath=args.javaPath,
        pymodulePath=args.pymodulePath,
        thisModulePath=args.thisModulePath,
        
        tmpDir=args.tmpDir,
        max_walltime=args.max_walltime,
        jvmVirtualByPhysicalMemoryRatio=args.jvmVirtualByPhysicalMemoryRatio,
        
        site_handler=args.site_handler, 
        input_site_handler=args.input_site_handler,
        cluster_size=args.cluster_size,

        needSSHDBTunnel=args.needSSHDBTunnel,
        debug=args.debug,
        report=args.report)
    instance.setup_run()
    instance.end_run()
