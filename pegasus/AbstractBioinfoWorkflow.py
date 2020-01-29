#!/usr/bin/env python
"""
2020/01/29
    an abstract class for pegasus workflows that work on bioinformatic data
"""
import sys, os, math

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.expanduser('~/script'))


from pymodule import ProcessOptions, getListOutOfStr, PassingData, utils
from pymodule.pegasus import yh_pegasus

from Pegasus.DAX3 import Executable, File, PFN, Link, Job
from AbstractWorkflow import AbstractWorkflow

parentClass = AbstractWorkflow
class AbstractBioinfoWorkflow(parentClass):
    __doc__ = __doc__
    option_default_dict = parentClass.option_default_dict.copy()

    def __init__(self,  **keywords):
        """
        20200129
        """
        parentClass.__init__(self, **keywords)

    def registerPlinkExecutables(self, workflow=None):
        if not workflow:
            workflow = self
        
        self.addExecutableFromPath(path=self.plinkPath, \
            name='plink', clusterSizeMultipler=1)
        self.addExecutableFromPath(path=self.plinkPath, \
            name='plinkNoClustering', clusterSizeMultipler=0)

        #2012.8.10 different plinks so that you can differentiate between different types of plink jobs
        self.addExecutableFromPath(path=self.plinkPath, \
            name='plinkMerge', clusterSizeMultipler=0)
        self.addExecutableFromPath(path=self.plinkPath, \
            name='plinkIBD', clusterSizeMultipler=0)
        self.addExecutableFromPath(path=self.plinkPath, \
            name='plinkConvert', clusterSizeMultipler=1)
        self.addExecutableFromPath(path=self.plinkPath, \
            name='plinkLDPrune', clusterSizeMultipler=1)
        self.addExecutableFromPath(path=self.plinkPath, \
            name='plinkExtract', clusterSizeMultipler=1)
        #2013.07.24
        self.addExecutableFromPath(path=os.path.join(self.pymodulePath, \
            'mapper/modifier/SplitPlinkLMendelFileSNPIDIntoChrPosition.py'), \
            name='SplitPlinkLMendelFileSNPIDIntoChrPosition', clusterSizeMultipler=1)
        
        #2013.07.19
        self.addExecutableFromPath(path=os.path.join(self.pymodulePath, \
            'pedigree/CalculateMendelErrorRateGivenPlinkOutput.py'), \
            name='CalculateMendelErrorRateGivenPlinkOutput', clusterSizeMultipler=1)



    def registerExecutables(self, workflow=None):
        """
        """
        if not workflow:
            workflow = self
        ParentClass.registerExecutables(self, workflow=workflow)

        namespace = self.namespace
        version = self.version
        operatingSystem = self.operatingSystem
        architecture = self.architecture
        clusters_size = self.clusters_size
        site_handler = self.site_handler


        #2013.11.22	#2013.06.25 register tabix
        self.addExecutableFromPath(
            path=self.tabixPath, name='tabix', clusterSizeMultipler=5)

        #2013.11.22 2011.12.21	for OutputVCFSiteStat.py
        self.addExecutableFromPath(
            path=os.path.join(self.pymodulePath, "mapper/extractor/tabixRetrieve.sh"),
            name='tabixRetrieve', clusterSizeMultipler=1)

        #2013.11.22 moved from pymodule/polymorphism/FindNewRefCoordinatesGivenVCFFolderWorkflow.py
        self.addExecutableFromPath(path=os.path.join(self.pymodulePath, \
            "polymorphism/mapper/LiftOverVCFBasedOnCoordinateMap.py"), \
            name='LiftOverVCFBasedOnCoordinateMap', clusterSizeMultipler=1)

        self.addExecutableFromPath(path=os.path.join(workflow.pymodulePath, \
            "polymorphism/qc/CalculateLociAndGenomeCoveredAtEachSwitchFrequencyThreshold.py"), \
            name='CalculateLociAndGenomeCoveredAtEachSwitchFrequencyThreshold', clusterSizeMultipler=0.01)

        self.addExecutableFromPath(path=os.path.join(workflow.pymodulePath, \
                "mapper/extractor/ExtractFlankingSequenceForVCFLoci.py"), \
            name='ExtractFlankingSequenceForVCFLoci', clusterSizeMultipler=2)

        self.addExecutableFromPath(path=os.path.join(workflow.pymodulePath, \
            "polymorphism/mapper/FindSNPPositionOnNewRefFromFlankingBlastOutput.py"), \
            name='FindSNPPositionOnNewRefFromFlankingBlastOutput', clusterSizeMultipler=2)

        self.addExecutableFromPath(path=os.path.join(workflow.pymodulePath, \
            "polymorphism/mapper/FindSNPPositionOnNewRefFromFlankingBWAOutput.py"), \
            name='FindSNPPositionOnNewRefFromFlankingBWAOutput', clusterSizeMultipler=1)


        #2013.08.28
        self.addExecutableFromPath(path=os.path.join(self.pymodulePath, 'Genome/OutputGenomeAnnotation.py'), \
            name='OutputGenomeAnnotation', clusterSizeMultipler=0.01)
        #2013.07.31
        self.addExecutableFromPath(path=os.path.join(self.pymodulePath, 'statistics/GenomeMovingAverageStatistics.py'), \
            name='GenomeMovingAverageStatistics', clusterSizeMultipler=0.1)
        #2013.08.23
        self.addExecutableFromPath(path=os.path.join(self.pymodulePath, 'reducer/ReduceSameChromosomeAlignmentDepthFiles'), \
            name='ReduceSameChromosomeAlignmentDepthFiles', clusterSizeMultipler=0.5)

        #executableList = []
        executableClusterSizeMultiplierList = []	#2012.8.7 each cell is a tuple of (executable, clusterSizeMultipler (0 if u do not need clustering)
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
        ConvertBjarniSNPFormat2Yu.addPFN(PFN("file://" +  os.path.join(self.variationSrcPath, "yhio/ConvertBjarniSNPFormat2Yu.py"), site_handler))
        executableClusterSizeMultiplierList.append((ConvertBjarniSNPFormat2Yu, 1))

        ConvertVCF2BjarniFormat = Executable(namespace=namespace, name="ConvertVCF2BjarniFormat", \
                                            version=version, os=operatingSystem, arch=architecture, installed=True)
        ConvertVCF2BjarniFormat.addPFN(PFN("file://" +  os.path.join(self.pymodulePath, "mapper/converter/ConvertVCF2BjarniFormat.py"), site_handler))
        executableClusterSizeMultiplierList.append((ConvertVCF2BjarniFormat, 1))


        ConvertYuSNPFormat2Bjarni = Executable(namespace=namespace, name="ConvertYuSNPFormat2Bjarni", \
                                            version=version, os=operatingSystem, arch=architecture, installed=True)
        ConvertYuSNPFormat2Bjarni.addPFN(PFN("file://" +  os.path.join(self.variationSrcPath, "yhio/ConvertYuSNPFormat2Bjarni.py"), site_handler))
        executableClusterSizeMultiplierList.append((ConvertYuSNPFormat2Bjarni, 1))

        ConvertYuSNPFormat2EigenStrat = Executable(namespace=namespace, name="ConvertYuSNPFormat2EigenStrat", \
                                            version=version, os=operatingSystem, arch=architecture, installed=True)
        ConvertYuSNPFormat2EigenStrat.addPFN(PFN("file://" +  os.path.join(self.variationSrcPath, "yhio/ConvertYuSNPFormat2EigenStrat.py"), site_handler))
        executableClusterSizeMultiplierList.append((ConvertYuSNPFormat2EigenStrat, 1))

        ConvertYuSNPFormat2TPED_TFAM = Executable(namespace=namespace, name="ConvertYuSNPFormat2TPED_TFAM", \
                                            version=version, os=operatingSystem, arch=architecture, installed=True)
        ConvertYuSNPFormat2TPED_TFAM.addPFN(PFN("file://" +  os.path.join(self.variationSrcPath, "yhio/ConvertYuSNPFormat2TPED_TFAM.py"), site_handler))
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

        self.addExecutables(executableClusterSizeMultiplierList, defaultClustersSize=self.clusters_size)

        self.addExecutableFromPath(path=os.path.join(self.pymodulePath, "plot/PlotVCFtoolsStat.py"), \
                                        name='PlotVCFtoolsStat', clusterSizeMultipler=0)
        
        #2013.07.19
        self.addExecutableFromPath(path=os.path.join(self.pymodulePath, 'mapper/modifier/AppendExtraPedigreeIndividualsToTPED.py'), \
                                        name='AppendExtraPedigreeIndividualsToTPED', clusterSizeMultipler=1)

        self.addExecutableFromPath(path=os.path.join(self.pymodulePath, 'mapper/converter/ConvertMSOutput2FASTQ.py'), \
                                        name='ConvertMSOutput2FASTQ', clusterSizeMultipler=1)

        self.addExecutableFromPath(path=os.path.join(self.pymodulePath, 'mapper/extractor/SelectChromosomeSequences.py'), \
                                        name='SelectChromosomeSequences', clusterSizeMultipler=0.5)

        #2013.2.11 moved from vervet/src/reduce to pymodule/reducer
        self.addExecutableFromPath(path=os.path.join(self.pymodulePath, 'reducer/MergeGenotypeMatrix.py'), \
                                        name='MergeGenotypeMatrix', clusterSizeMultipler=0.2)
        self.addExecutableFromPath(path=os.path.join(self.pymodulePath, 'plot/PlotGenomeWideData.py'), \
                                        name='PlotGenomeWideData', clusterSizeMultipler=1)

        self.bgzipExecutableFile = self.registerOneExecutableAsFile(path=os.path.expanduser("~/bin/bgzip"))	#2013.11.22

if __name__ == '__main__':
    main_class = AbstractBioinfoWorkflow
    po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
    instance = main_class(**po.long_option2value)
    instance.run()