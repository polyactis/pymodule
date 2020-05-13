#!/usr/bin/env python3
"""
2013.1.25 an abstract class for pegasus workflows that work on alignment & VCF files.
"""
import sys, os, math
import copy
from pegaflow.DAX3 import Executable, File, PFN, Link, Job
from palos import ProcessOptions, getListOutOfStr, PassingData, utils
from palos import ngs
from . AbstractAlignmentWorkflow import AbstractAlignmentWorkflow
from . AbstractVCFWorkflow import AbstractVCFWorkflow

ParentClass = AbstractAlignmentWorkflow
class AbstractAlignmentAndVCFWorkflow(ParentClass, AbstractVCFWorkflow):
    __doc__ = __doc__
    option_default_dict = copy.deepcopy(ParentClass.option_default_dict)
    option_default_dict.update({
        ('inputDir', 0, ): ['', 'L', 1, 'input folder that contains vcf or vcf.gz files', ],\
        ('minDepth', 0, float): [0, 'm', 1, 
            'minimum depth for a call to regarded as non-missing', ],\
        ('intervalOverlapSize', 1, int): [300000, 'U', 1, 
            'overlap #bps/#loci between adjacent intervals from one contig/chromosome,'
            'only used for TrioCaller, not for SAMtools/GATK', ],\
        ('intervalSize', 1, int): [5000000, 'Z', 1,
            '#bps/#loci for adjacent intervals from one contig/chromosome (alignment or VCF)', ],\
        })
    def __init__(self,  **keywords):
        """
        2012.1.17
        """
        ParentClass.__init__(self, **keywords)

    registerFilesOfInputDir = AbstractVCFWorkflow.registerFilesOfInputDir

    def setup(self, inputVCFData=None, chr2IntervalDataLs=None, **keywords):
        """
        2013.04.01 derive chr2VCFJobData only when inputVCFData is available
        2013.1.25
        """
        pdata = ParentClass.setup(self, chr2IntervalDataLs=chr2IntervalDataLs,
            inputVCFData=inputVCFData, **keywords)

        #2012.8.26 so that each recalibration will pick up the right vcf
        chr2VCFJobData = {}
        if inputVCFData:
            for jobData in inputVCFData.jobDataLs:
                inputF = jobData.file
                chromosome = self.getChrFromFname(os.path.basename(inputF.name))
                chr2VCFJobData[chromosome] = jobData
        pdata.chr2VCFJobData = chr2VCFJobData
        return pdata

    def registerCustomExecutables(self):

        """
        """
        ParentClass.registerCustomExecutables(self)
        AbstractVCFWorkflow.registerCustomExecutables(self)

    def run(self):
        """
        2011-9-28
        """
        pdata = self.setup_run()

        inputData = self.registerFilesOfInputDir(
            inputDir=self.inputDir,
            input_site_handler=self.input_site_handler, \
            checkEmptyVCFByReading=self.checkEmptyVCFByReading,\
            pegasusFolderName=self.pegasusFolderName)
        if len(inputData.jobDataLs)<=0:
            sys.stderr.write("Error: No VCF files in the input VCF folder %s.\n"%self.inputDir)
            raise
        #adding inputVCFData=... is the key difference from the parent class
        self.addAllJobs(inputVCFData=inputData,
            alignmentDataLs=pdata.alignmentDataLs,
            chr2IntervalDataLs=pdata.chr2IntervalDataLs,
            samtools=self.samtools,
            GenomeAnalysisTKJar=self.GenomeAnalysisTKJar,
            MergeSamFilesJar=self.MergeSamFilesJar, \
            CreateSequenceDictionaryJava=self.CreateSequenceDictionaryJava,
            CreateSequenceDictionaryJar=self.CreateSequenceDictionaryJar,
            BuildBamIndexFilesJava=self.BuildBamIndexFilesJava,
            BuildBamIndexJar=self.BuildBamIndexJar,\
            mv=self.mv, skipDoneAlignment=self.skipDoneAlignment, \
            registerReferenceData=pdata.registerReferenceData,\
            needFastaIndexJob=self.needFastaIndexJob,
            needFastaDictJob=self.needFastaDictJob, \
            data_dir=self.data_dir, no_of_gatk_threads = 1,
            transferOutput=True,
            outputDirPrefix=self.pegasusFolderName)

        self.end_run()

if __name__ == '__main__':
    main_class = AbstractAlignmentAndVCFWorkflow
    po = ProcessOptions(sys.argv, main_class.option_default_dict,
        error_doc=main_class.__doc__)
    instance = main_class(**po.long_option2value)
    instance.run()
