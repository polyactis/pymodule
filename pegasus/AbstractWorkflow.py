#!/usr/bin/env python3
"""
2012.5.23
    a common class for other pymodule workflows
"""
import sys, os, math

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.expanduser('~/script'))

from pymodule import ProcessOptions, getListOutOfStr, PassingData, utils
import yh_pegasus
from pegapy3.Workflow import Workflow
from pegapy3.DAX3 import Executable, File, PFN, Profile, Namespace, Link, Use, Job, Dependency

class AbstractWorkflow(Workflow):
    __doc__ = __doc__
    db_option_dict = {
                    ('drivername', 1,):['postgresql', 'v', 1, 'which type of database? mysql or postgresql', ],\
                    ('hostname', 1, ): ['localhost', 'z', 1, 'hostname of the db server', ],\
                    ('dbname', 1, ): ['', 'd', 1, 'database name', ],\
                    ('schema', 0, ): ['public', 'k', 1, 'database schema name', ],\
                    ('db_user', 1, ): [None, 'u', 1, 'database username', ],\
                    ('db_passwd', 1, ): [None, 'p', 1, 'database password', ],\
                    ('port', 0, ):[None, '', 1, 'database port number'],\
                    ('commit', 0, ):[None, '', 0, 'commit database transaction if there is db transaction'],\
                    ("data_dir", 0, ): ["", 't', 1, 'the base directory where all db-affiliated files are stored. '
                                    'If not given, use the default stored in db.'],\
                    ("local_data_dir", 0, ): ["", 'D', 1, 'this one should contain same files as data_dir but accessible locally. '
                            'If not given, use the default stored in db (db.data_dir). This argument is used to find all input files available.\n '
                            'It should be different from data_dir only when you generate a workflow on one computer and execute it on another which has different data_dir.'],\

                    }
    option_default_dict = {
                        ("pymodulePath", 1, ): ["%s/script/pymodule", '', 1, 'path to the pymodule folder'],\
                        ("variationSrcPath", 1, ): ["%s/script/variation/src", '', 1, 'variation source code folder'],\
                        ("home_path", 1, ): [os.path.expanduser("~"), 'e', 1, 'path to the home directory on the working nodes'],\
                        ("javaPath", 1, ): ["%s/bin/jdk/bin/java", 'J', 1, 'path to java interpreter binary'],\
                        ("plinkPath", 1, ): ["%s/bin/plink", '', 1, 'path to the plink binary, http://pngu.mgh.harvard.edu/~purcell/plink/index.shtml'],\
                        ("pegasusCleanupPath", 1, ): ["%s/bin/pegasus/bin/pegasus-cleanup", '', 1, 'path to pegasus-cleanup executable, it will be registered and run on local universe of condor pool (rather than the vanilla universe)'],\
                        ("pegasusTransferPath", 1, ): ["%s/bin/pegasus/bin/pegasus-transfer", '', 1, 'path to pegasus-transfer executable, it will be registered and run on local universe of condor pool (rather than the vanilla universe)'],\
                        ("site_handler", 1, ): ["condorpool", 'l', 1, 'which site to run the jobs: condorpool, hoffman2'],\
                        ("input_site_handler", 1, ): ["local", 'j', 1, 'which site has all the input files: local, condorpool, hoffman2. '
                            'If site_handler is condorpool, this must be condorpool and files will be symlinked. '
                            'If site_handler is hoffman2, input_site_handler=local induces file transfer and input_site_handler=hoffman2 induces symlink.'],\
                        ('clusters_size', 1, int):[30, 'C', 1, 'For short jobs that will be clustered, how many of them should be clustered int one'],\
                        ('pegasusFolderName', 0, ): ['folder', 'F', 1, 'the folder relative to pegasus workflow root to contain input & output. '
                                'It will be created during the pegasus staging process. It is useful to separate multiple workflows. '
                                'If empty, everything is in the pegasus root.', ],\
                        ('inputSuffixList', 0, ): [None, '', 1, 'coma-separated list of input file suffices. If None, any suffix. '
                            'Suffix include the dot, (i.e. .tsv). Typical zip suffices are excluded (.gz, .bz2, .zip, .bz).'],\
                        ('outputFname', 1, ): [None, 'o', 1, 'xml workflow output file'],\
                        ("tmpDir", 1, ): ["/tmp/", '', 1, 'for MarkDuplicates.jar, etc., default is /tmp/ but sometimes it is too small'],\
                        ('max_walltime', 1, int):[4320, '', 1, 'maximum wall time any job could have, in minutes. 20160=2 weeks.\n'
                            'used in addGenericJob().'],\
                        ('jvmVirtualByPhysicalMemoryRatio', 1, float):[1.0, '', 1, 
                            "if a job's virtual memory (usually 1.2X of JVM resident memory) exceeds request, "
                            "it will be killed on hoffman2. Hence this argument"],\
                        ("thisModulePath", 1, ): ["%s", '', 1, 'path of the module that owns this program. '
                            'used to add executables from this module.'],\
                        ('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
                        ('needSSHDBTunnel', 0, int):[0, 'H', 0, 'DB-interacting jobs need a ssh tunnel (running on cluster behind firewall).'],\
                        ('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']
                        }
                        #('bamListFname', 1, ): ['/tmp/bamFileList.txt', 'L', 1, 'The file contains path to each bam file, one file per line.'],\

    pathToInsertHomePathList = ['javaPath', 'pymodulePath', 'plinkPath', 'variationSrcPath', 'pegasusCleanupPath',\
                            'pegasusTransferPath', "thisModulePath"]

    def __init__(self, inputArgumentLs=None, **keywords):
        """
        2013.06.27 add argumen inputArgumentLs to include everything on the tail of the commandline
        2012.5.23
        """
        Workflow.__init__(self, inputArgumentLs=None, **keywords)
        self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, \
                                                        class_to_have_attr=self)
        #2013.11.24
        self.inputSuffixList = utils.getListOutOfStr(list_in_str=self.inputSuffixList, data_type=str, separator1=',', separator2='-')
        self.inputSuffixSet = set(self.inputSuffixList)
        #2013.06.27
        self.inputArgumentLs = inputArgumentLs
        if self.inputArgumentLs is None:
            self.inputArgumentLs = []

    def constructJobDataFromJob(self, job=None):
        """
        2013.09.05 added vcfFile and tbi_F in structure
        2013.07.18
        """
        if hasattr(job, 'output') and job.output.name.find('.vcf')>=0:
            vcfFile = job.output
            tbi_F = getattr(job, 'tbi_F', None)
        else:
            vcfFile = None
            tbi_F = None
        return PassingData(job=job, jobLs=[job], file=job.output, fileLs=job.outputLs, vcfFile=vcfFile,\
                        tbi_F=tbi_F)

    def registerExecutables(self):
        """
        2012.1.9
        """
        Workflow.registerExecutables(self)
        
        #2013.2.7 convert, an image swissknife program, part of imagemagick
        self.addExecutableFromPath(path="/usr/bin/convert", name='convertImage', clusterSizeMultipler=1)

        #2013.08.23 c++ version of SelectRowsFromMatrix.py
        self.addExecutableFromPath(path=os.path.join(self.pymodulePath, 'mapper/extractor/SelectRowsFromMatrixCC'), \
                                        name='SelectRowsFromMatrixCC', clusterSizeMultipler=1)
        #2012.08.13 SelectRowsFromMatrix is a derivative of AbstractMatrixFileWalker, so use addAbstractMatrixFileWalkerJob()
        self.addExecutableFromPath(path=os.path.join(self.pymodulePath, 'mapper/extractor/SelectRowsFromMatrix.py'), \
                                        name='SelectRowsFromMatrix', clusterSizeMultipler=1)
        self.addExecutableFromPath(path=os.path.join(self.pymodulePath, "mapper/extractor/SelectLineBlockFromFile.py"), 
            name='SelectLineBlockFromFile', clusterSizeMultipler=1)
        self.addExecutableFromPath(path=os.path.join(self.pymodulePath, "plot/AbstractPlot.py"), 
            name='AbstractPlot', clusterSizeMultipler=1)
        self.addExecutableFromPath(path=os.path.join(self.pymodulePath, "plot/PlotYAsBar.py"), 
            name='PlotYAsBar', clusterSizeMultipler=1)
        self.addExecutableFromPath(path=os.path.join(self.pymodulePath, "plot/DrawHistogram.py"), 
            name='DrawHistogram', clusterSizeMultipler=1)

        #2012.8.15 ancestor of SelectRowsFromMatrix,
        self.addExecutableFromPath(path=os.path.join(self.pymodulePath, "yhio/AbstractMatrixFileWalker.py"), 
            name='AbstractMatrixFileWalker', clusterSizeMultipler=1)
        self.addExecutableFromPath(path=self.javaPath, name='java', clusterSizeMultipler=1)
        self.addExecutableFromPath(path=os.path.join(self.pymodulePath, "plot/DrawMatrix.py"), 
            name='DrawMatrix', clusterSizeMultipler=1)
        self.addExecutableFromPath(path=os.path.join(self.pymodulePath, "plot/Draw2DHistogramOfMatrix.py"), 
            name='Draw2DHistogramOfMatrix', clusterSizeMultipler=1)
        # C++ binary
        self.addExecutableFromPath(path=os.path.join(self.pymodulePath, "mapper/CalculateMedianMeanOfInputColumn"), 
            name='CalculateMedianMeanOfInputColumn', clusterSizeMultipler=1)
        self.addExecutableFromPath(path=os.path.join(self.pymodulePath, "statistics/SampleRows.py"), 
            name='SampleRows', clusterSizeMultipler=1)
        #2013.2.11 all reducers
        self.addExecutableFromPath(path=os.path.join(self.pymodulePath, "statistics/EstimateOutliersIn2DData.py"), \
                name='EstimateOutliersIn2DData', clusterSizeMultipler=0)
        self.addExecutableFromPath(path=os.path.join(self.pymodulePath, 'reducer/MergeSameHeaderTablesIntoOne.py'), \
                name='mergeSameHeaderTablesIntoOne', clusterSizeMultipler=0)
        self.addExecutableFromPath(path=os.path.join(self.pymodulePath, 'reducer/MergeSameHeaderTablesIntoOne.py'), \
                name='MergeSameHeaderTablesIntoOne', clusterSizeMultipler=0)
        self.addExecutableFromPath(path=os.path.join(self.pymodulePath, 'reducer/ReduceMatrixByAverageColumnsWithSameKey.py'), \
                name='ReduceMatrixByAverageColumnsWithSameKey', clusterSizeMultipler=0)
        self.addExecutableFromPath(path=os.path.join(self.pymodulePath, 'reducer/ReduceMatrixByChosenColumn.py'), \
                name='ReduceMatrixByChosenColumn', clusterSizeMultipler=0)
        self.addExecutableFromPath(path=os.path.join(self.pymodulePath, 'reducer/ReduceMatrixByMergeColumnsWithSameKey.py'), \
                name='ReduceMatrixByMergeColumnsWithSameKey', clusterSizeMultipler=0)
        self.addExecutableFromPath(path=os.path.join(self.pymodulePath, 'reducer/ReduceMatrixBySumSameKeyColsAndThenDivide.py'), \
                name='ReduceMatrixBySumSameKeyColsAndThenDivide', clusterSizeMultipler=0)
    
    def addSortJob(self, executable=None, executableFile=None, \
                    inputFile=None, outputFile=None, noOfHeaderLines=0, \
                    parentJobLs=None, extraDependentInputLs=None, extraOutputLs=None, transferOutput=False, \
                    extraArguments=None, extraArgumentList=None, sshDBTunnel=None,\
                    job_max_memory=2000, walltime=120, **keywords):
        """
        Examples:

            sortedSNPID2NewCoordinateFile = File(os.path.join(reduceOutputDirJob.output, 'SNPID2NewCoordinates.sorted.tsv.gz'))
            sortSNPID2NewCoordinatesJob = self.addSortJob(inputFile=reduceJob.output, \
                        outputFile=sortedSNPID2NewCoordinateFile, noOfHeaderLines=1, \
                        parentJobLs=[reduceJob], \
                        extraOutputLs=None, transferOutput=False, \
                        extraArgumentList=["-k3,3 -k4,4n"], \
                        sshDBTunnel=None,\
                        job_max_memory=4000, walltime=120)
            # -t$'\t' for sort has to be removed as it won't be passed correctly.
            # the default sort field separator (non-blank to blank) works if no-blank is part of each cell value.
        2014.01.10 use sortHeaderAware executable (from pymodule/shell)
        2013.11.27 use unix sort to sort a file

        """
        if executable is None:
            executable = self.sortHeaderAware
        #if executableFile is None:
        #	executableFile = self.sortExecutableFile
        if extraDependentInputLs is None:
            extraDependentInputLs = []
        if extraArgumentList is None:
            extraArgumentList = []

        extraArgumentList.insert(0, "%s"%(noOfHeaderLines))
        job = self.addGenericJob(executable=executable, \
                    inputFile=inputFile, inputArgumentOption="", \
                    outputFile=outputFile, outputArgumentOption="",
                    parentJobLs=parentJobLs, \
                    extraDependentInputLs=extraDependentInputLs, \
                    extraOutputLs=extraOutputLs, transferOutput=transferOutput, \
                    extraArguments=extraArguments, \
                    extraArgumentList=extraArgumentList, \
                    sshDBTunnel=sshDBTunnel,\
                    job_max_memory=job_max_memory, walltime=walltime)
        return job

    def addStatMergeJob(self, statMergeProgram=None, outputF=None, \
                    parentJobLs=None, extraOutputLs=None,\
                    extraDependentInputLs=None, transferOutput=True, \
                    extraArguments=None, extraArgumentList=None,\
                    key2ObjectForJob=None,\
                    namespace=None, version=None, job_max_memory=1000, **keywords):
        """
        2012.8.10
            use addGenericJob()
        2012.4.3
            make argument namespace, version optional
        2011-11-28
            moved from CalculateVCFStatPipeline.py
        2011-11-17
            add argument extraArguments
        """
        if extraDependentInputLs is None:
            extraDependentInputLs = []

        if extraArgumentList is None:
            extraArgumentList = []
        if key2ObjectForJob is None:
            key2ObjectForJob = {}

        if extraArguments:
            extraArgumentList.append(extraArguments)
        job= self.addGenericJob(executable=statMergeProgram, inputFile=None, outputFile=outputF, \
                parentJobLs=parentJobLs, extraDependentInputLs=extraDependentInputLs, \
                extraOutputLs=extraOutputLs,\
                transferOutput=transferOutput, \
                extraArgumentList=extraArgumentList, key2ObjectForJob=key2ObjectForJob, job_max_memory=job_max_memory, **keywords)
        return job

    def addConvertImageJob(self, inputFile=None, inputArgumentOption=None, \
                    outputFile=None, outputArgumentOption=None, density=None, \
                    resizeDimension=None, \
                    parentJobLs=None, extraDependentInputLs=None, extraOutputLs=None, transferOutput=False, \
                    frontArgumentList=None, extraArguments=None, extraArgumentList=None, job_max_memory=2000,\
                    key2ObjectForJob=None, **keywords):
        """
        2013.2.7 use imagemagick's convert to convert images. examples:
            plotOutputFile = File('%s.eps'%(plotOutputFnamePrefix))
            plotPNGOutputFile = File('%s.png'%(plotOutputFnamePrefix))
            #change dpi to 300
            self.addConvertImageJob(inputFile=plotOutputFile, \
                    outputFile=plotPNGOutputFile, density=300, \
                    resizeDimension=None, \
                    parentJobLs=[psmc_plotJob], extraDependentInputLs=None, extraOutputLs=None, transferOutput=True, \
                    extraArguments=None, frontArgumentList=None, job_max_memory=500)

            #resize by demanding the width = 1800, height will scale accordingly
            self.addConvertImageJob(inputFile=plotOutputFile, \
                    outputFile=plotPNGOutputFile, density=None, \
                    resizeDimension=1800, \
                    parentJobLs=[psmc_plotJob], extraDependentInputLs=None, extraOutputLs=None, transferOutput=True, \
                    extraArguments=None, frontArgumentList=None, job_max_memory=500)
            #resize by demanding the dimension=1800X900
            self.addConvertImageJob(inputFile=plotOutputFile, \
                    outputFile=plotPNGOutputFile, density=None, \
                    resizeDimension='1800X900', \
                    parentJobLs=[psmc_plotJob], extraDependentInputLs=None, extraOutputLs=None, transferOutput=True, \
                    extraArguments=None, frontArgumentList=None, job_max_memory=500)
        """
        if frontArgumentList is None:
            frontArgumentList = []
        if extraOutputLs is None:
            extraOutputLs = []
        if density is not None:
            frontArgumentList.append("-density %s"%(density))
        if resizeDimension is not None:
            frontArgumentList.append("-resize %s"%(resizeDimension))
        #do not pass the inputFileList to addGenericJob() because db arguments need to be added before them.
        job = self.addGenericJob(executable=self.convertImage, inputFile=inputFile, \
                        inputArgumentOption=inputArgumentOption, outputFile=outputFile, \
                        outputArgumentOption=outputArgumentOption, inputFileList=None, parentJobLs=parentJobLs, \
                        extraDependentInputLs=extraDependentInputLs, extraOutputLs=extraOutputLs, \
                        transferOutput=transferOutput, \
                        frontArgumentList=frontArgumentList, extraArguments=extraArguments, extraArgumentList=extraArgumentList,\
                        job_max_memory=job_max_memory, key2ObjectForJob=key2ObjectForJob,\
                        **keywords)
        return job

    def addGenericDBJob(self, executable=None, inputFile=None, inputArgumentOption="-i", \
                    inputFileList=None, argumentForEachFileInInputFileList=None,\
                    outputFile=None, outputArgumentOption="-o", \
                    parentJobLs=None, extraDependentInputLs=None, extraOutputLs=None, transferOutput=False, \
                    extraArguments=None, extraArgumentList=None, job_max_memory=2000,  sshDBTunnel=None, \
                    key2ObjectForJob=None, objectWithDBArguments=None, walltime=None, **keywords):
        """
        2013.3.25 bugfix. moved most post-addGenericDBJob code into addGenericJob()
        2012.10.6
            similar to addGenericJob but these are generic jobs that need db-interacting arguments.

            Argument inputFileList is a list of input files to be added to commandline as the last arguments.
                they would also be added as the job's dependent input.
                Difference from extraDependentInputLs: the latter is purely for dependency purpose, not added as input arguments.
                    So if files have been put in inputFileList, then they shouldn't be in extraDependentInputLs.
        """
        if objectWithDBArguments is None:
            objectWithDBArguments = self
        job = self.addGenericJob(executable=executable, \
                        inputFile=inputFile, inputArgumentOption=inputArgumentOption, \
                        outputFile=outputFile, outputArgumentOption=outputArgumentOption, \
                        inputFileList=inputFileList, argumentForEachFileInInputFileList=argumentForEachFileInInputFileList, \
                        parentJobLs=parentJobLs, \
                        extraDependentInputLs=extraDependentInputLs, extraOutputLs=extraOutputLs, \
                        transferOutput=transferOutput, extraArguments=extraArguments, extraArgumentList=extraArgumentList,\
                        job_max_memory=job_max_memory, sshDBTunnel=sshDBTunnel, key2ObjectForJob=key2ObjectForJob,\
                        objectWithDBArguments=objectWithDBArguments, walltime=walltime,\
                        **keywords)

        #2012.10.6 set the job.input
        if getattr(job, 'input', None) is None and job.inputLs:
            job.input = job.inputLs[0]
        return job

    def addGenericFile2DBJob(self, executable=None, inputFile=None, inputArgumentOption="-i", \
                    inputFileList=None, argumentForEachFileInInputFileList=None,\
                    outputFile=None, outputArgumentOption="-o", \
                    data_dir=None, logFile=None, commit=False,\
                    parentJobLs=None, extraDependentInputLs=None, extraOutputLs=None, transferOutput=False, \
                    extraArguments=None, extraArgumentList=None, job_max_memory=2000,  sshDBTunnel=None, \
                    key2ObjectForJob=None, objectWithDBArguments=None, **keywords):
        """

            job = self.addGenericFile2DBJob(executable=executable, inputFile=None, inputArgumentOption="-i", \
                    outputFile=None, outputArgumentOption="-o", inputFileList=None, \
                    data_dir=None, logFile=logFile, commit=commit,\
                    parentJobLs=parentJobLs, extraDependentInputLs=extraDependentInputLs, \
                    extraOutputLs=None, transferOutput=transferOutput, \
                    extraArguments=extraArguments, extraArgumentList=extraArgumentList, \
                    job_max_memory=job_max_memory,  sshDBTunnel=sshDBTunnel, walltime=walltime,\
                    key2ObjectForJob=None, objectWithDBArguments=self, **keywords)

        2012.12.21 a generic wrapper for jobs that "inserting" data (from file) into database
        """
        if extraArgumentList is None:
            extraArgumentList = []
        if extraOutputLs is None:
            extraOutputLs = []

        if data_dir:
            extraArgumentList.append('--data_dir %s'%(data_dir))
        if commit:
            extraArgumentList.append('--commit')
        if logFile:
            extraArgumentList.extend(["--logFilename", logFile])
            extraOutputLs.append(logFile)
        #do not pass the inputFileList to addGenericJob() because db arguments need to be added before them.
        job = self.addGenericDBJob(executable=executable, inputFile=inputFile, \
                        inputArgumentOption=inputArgumentOption, \
                        inputFileList=inputFileList, argumentForEachFileInInputFileList=argumentForEachFileInInputFileList,\
                        outputFile=outputFile, \
                        outputArgumentOption=outputArgumentOption, \
                        parentJobLs=parentJobLs, \
                        extraDependentInputLs=extraDependentInputLs, extraOutputLs=extraOutputLs, \
                        transferOutput=transferOutput, extraArguments=extraArguments, extraArgumentList=extraArgumentList,\
                        job_max_memory=job_max_memory, sshDBTunnel=sshDBTunnel, key2ObjectForJob=key2ObjectForJob,\
                        objectWithDBArguments=objectWithDBArguments, **keywords)
        return job

    def addGenericJavaJob(self, executable=None, jarFile=None, \
                    inputFile=None, inputArgumentOption=None, \
                    inputFileList=None,argumentForEachFileInInputFileList=None,\
                    outputFile=None, outputArgumentOption=None,\
                    frontArgumentList=None, extraArguments=None, extraArgumentList=None, extraOutputLs=None, \
                    extraDependentInputLs=None, \
                    parentJobLs=None, transferOutput=True, job_max_memory=2000,\
                    key2ObjectForJob=None, no_of_cpus=None, walltime=120, sshDBTunnel=None, **keywords):
        """
        2013.3.23 a generic function to add Java jobs

            fastaDictJob = self.addGenericJavaJob(executable=CreateSequenceDictionaryJava, jarFile=CreateSequenceDictionaryJar, \
                    inputFile=refFastaF, inputArgumentOption="REFERENCE=", \
                    inputFileList=None, argumentForEachFileInInputFileList=None,\
                    outputFile=refFastaDictF, outputArgumentOption="OUTPUT=",\
                    parentJobLs=parentJobLs, transferOutput=transferOutput, job_max_memory=job_max_memory,\
                    frontArgumentList=None, extraArguments=None, extraArgumentList=None, extraOutputLs=None, \
                    extraDependentInputLs=None, no_of_cpus=None, walltime=walltime, sshDBTunnel=None, **keywords)
        """
        if executable is None:
            executable = self.java
        if frontArgumentList is None:
            frontArgumentList = []
        if extraArgumentList is None:
            extraArgumentList = []
        if extraDependentInputLs is None:
            extraDependentInputLs = []
        if extraOutputLs is None:
            extraOutputLs = []


        
        memRequirementObject = self.getJVMMemRequirment(job_max_memory=job_max_memory, minMemory=4000)
        job_max_memory = memRequirementObject.memRequirement
        javaMemRequirement = memRequirementObject.memRequirementInStr

        frontArgumentList = [javaMemRequirement, '-jar', jarFile] + frontArgumentList	#put java stuff in front of other fron arguments
        extraDependentInputLs.append(jarFile)
        job = self.addGenericJob(executable=executable, inputFile=inputFile, \
                    inputArgumentOption=inputArgumentOption,  inputFileList=inputFileList,\
                    argumentForEachFileInInputFileList=argumentForEachFileInInputFileList,\
                    outputFile=outputFile, outputArgumentOption=outputArgumentOption, \
                    parentJob=None, parentJobLs=parentJobLs, extraDependentInputLs=extraDependentInputLs, \
                    extraOutputLs=extraOutputLs, \
                    transferOutput=transferOutput, \
                    frontArgumentList=frontArgumentList, extraArguments=extraArguments, \
                    extraArgumentList=extraArgumentList, job_max_memory=job_max_memory,  sshDBTunnel=sshDBTunnel, \
                    key2ObjectForJob=key2ObjectForJob, no_of_cpus=no_of_cpus, walltime=walltime, **keywords)

        return job

    def setJobOutputFileTransferFlag(self, job=None, transferOutput=False, outputLs=None):
        """
        2012.8.17
            assume all output files in job.outputLs
        """
        if outputLs is None and getattr(job, 'outputLs', None):
            outputLs = job.outputLs

        for output in outputLs:
            job.uses(output, transfer=transferOutput, register=True, link=Link.OUTPUT)

    def addCalculateDepthMeanMedianModeJob(self, executable=None, \
                            inputFile=None, outputFile=None, alignmentID=None, fractionToSample=0.001, \
                            whichColumn=None, maxNumberOfSamplings=1E7, inputStatName=None,\
                            parentJobLs=None, job_max_memory = 500, extraArguments=None, \
                            transferOutput=False, **keywords):
        """
        2013.1.8 moved from vervet.src.alignment.InspectAlignmentPipeline and use addGenericJob()
        2012.6.15 turn maxNumberOfSamplings into integer when passing it to the job
        2012.5.7
            a job to take input of samtoolsDepth
        """
        extraArgumentList = []
        if alignmentID is not None:
            extraArgumentList.append("--alignmentID %s"%(alignmentID))
        if fractionToSample is not None:
            extraArgumentList.append("--fractionToSample %s"%(fractionToSample))
        if whichColumn is not None:
            extraArgumentList.append("--whichColumn %s"%(whichColumn))
        if maxNumberOfSamplings is not None:
            extraArgumentList.append('--maxNumberOfSamplings %d'%(maxNumberOfSamplings))
        if inputStatName is not None:
            extraArgumentList.append("--inputStatName %s"%(inputStatName))
        if extraArguments:
            extraArgumentList.append(extraArguments)
        job= self.addGenericJob(executable=executable, inputFile=inputFile, outputFile=outputFile, \
                parentJobLs=parentJobLs, extraDependentInputLs=None, \
                extraOutputLs=None,\
                transferOutput=transferOutput, \
                extraArgumentList=extraArgumentList, key2ObjectForJob=None, \
                sshDBTunnel=None, job_max_memory=job_max_memory, **keywords)
        return job

    def addDBArgumentsToOneJob(self, job=None, objectWithDBArguments=None):
        """
        2012.8.17
            use long arguments , rather than short ones
        2012.6.5
            tired of adding all these arguments to db-interacting jobs
        """
        if objectWithDBArguments is None:
            objectWithDBArguments = self
        job.addArguments("--drivername", objectWithDBArguments.drivername, "--hostname", objectWithDBArguments.hostname, \
                        "--dbname", objectWithDBArguments.dbname, \
                        "--db_user", objectWithDBArguments.db_user, "--db_passwd %s"%objectWithDBArguments.db_passwd)
        if objectWithDBArguments.schema:
            job.addArguments("--schema", objectWithDBArguments.schema)
        if getattr(objectWithDBArguments, 'port', None):
            job.addArguments("--port=%s"%(objectWithDBArguments.port))
        return job

    def addDBGenomeArgumentsToOneJob(self, job=None, objectWithDBArguments=None):
        """
        2013.07.31 similar to addDBArgumentsToOneJob() but for genome db
        """
        if objectWithDBArguments is None:
            objectWithDBArguments = self
        if self.drivername=='mysql':
            genome_dbname = 'genome'
        else:
            genome_dbname = self.dbname
        job.addArguments("--genome_drivername", objectWithDBArguments.drivername, \
                        "--genome_hostname", objectWithDBArguments.hostname, \
                        "--genome_dbname", genome_dbname, \
                        "--genome_db_user", objectWithDBArguments.db_user, \
                        "--genome_db_passwd %s"%objectWithDBArguments.db_passwd)

        job.addArguments("--genome_schema genome")

        if getattr(objectWithDBArguments, 'port', None):
            job.addArguments("--genome_port=%s"%(objectWithDBArguments.port))
        return job


    def addGzipSubWorkflow(self, inputData=None, transferOutput=True,\
                        outputDirPrefix="", topOutputDirJob=None, report=True, **keywords):
        """
        2012.8.2 bugfix.
        2012.7.19
        """
        if report:
            sys.stderr.write("Adding gzip jobs for %s input job data ... "%(len(inputData.jobDataLs)))

        returnData = PassingData(topOutputDirJob=None)
        returnData.jobDataLs = []
        if inputData:
            if len(inputData.jobDataLs)>0:
                if topOutputDirJob is None:
                    topOutputDir = "%sGzip"%(outputDirPrefix)
                    topOutputDirJob = self.addMkDirJob(outputDir=topOutputDir)

                returnData.topOutputDirJob = topOutputDirJob
                for jobData in inputData.jobDataLs:
                    for inputF in set(jobData.fileLs):	#2013.08.16 do not work on same file
                        inputFBaseName = os.path.basename(inputF.name)
                        outputF = File(os.path.join(topOutputDirJob.output, '%s.gz'%(inputFBaseName)))
                        key2ObjectForJob = {}
                        extraArgumentList = []
                        #make sure set inputArgumentOption&outputArgumentOption to None, \
                        # otherwise addGenericJob will add "-i" and "-o" in front of it
                        job= self.addGenericJob(executable=self.gzip, inputFile=inputF,
                                    inputArgumentOption=None, outputArgumentOption=None,  outputFile=outputF, \
                                    parentJobLs=[topOutputDirJob]+jobData.jobLs, extraDependentInputLs=None, \
                                    extraOutputLs=[],\
                                    transferOutput=transferOutput, \
                                    extraArgumentList=extraArgumentList, key2ObjectForJob=key2ObjectForJob, \
                                    job_max_memory=200, **keywords)
                        """
                        # 2012.8.2 wrong, because -i and -o will be added in front.
                        abstractMapperJob = self.addAbstractMapperLikeJob(workflow, executable=workflow.gzip, \
                                inputF=None, outputF=outputF, \
                                parentJobLs=[topOutputDirJob]+jobData.jobLs, transferOutput=transferOutput, job_max_memory=200,\
                                extraArguments=None, extraDependentInputLs=[inputF], )
                        """
                        returnData.jobDataLs.append(PassingData(jobLs=[job], vcfFile=None, file=outputF,\
                                                fileLs=[outputF]))
        if report:
            sys.stderr.write("no_of_jobs = %s.\n"%(self.no_of_jobs))
        return returnData

    def addAbstractMapperLikeJob(self, executable=None, \
                    inputVCF=None, inputF=None, outputF=None, extraOutputLs=None,\
                    parentJobLs=None, transferOutput=True, job_max_memory=2000,\
                    extraArguments=None, extraArgumentList=None, extraDependentInputLs=None, \
                    sshDBTunnel=None, **keywords):
        """
        2012.10.8 call addGenericJob() instead
        2012.7.19
            moved from AbstractNGSWorkflow to here.
            add argument inputF. inputVCF is not generic enough.
        2012.5.11
        """
        if inputF is None:	#2012.7.19
            inputF = inputVCF
        job= self.addGenericJob(executable=executable, inputFile=inputF, outputFile=outputF, \
                parentJobLs=parentJobLs, extraDependentInputLs=extraDependentInputLs, \
                extraOutputLs=extraOutputLs,\
                transferOutput=transferOutput, \
                extraArguments=extraArguments,\
                extraArgumentList=extraArgumentList, \
                sshDBTunnel=sshDBTunnel, job_max_memory=job_max_memory, **keywords)
        return job

    def addSelectLineBlockFromFileJob(self, executable=None, inputFile=None, outputFile=None,\
                    startLineNumber=None, stopLineNumber=None, parentJobLs=None, extraDependentInputLs=None, \
                    transferOutput=False, \
                    extraArguments=None, job_max_memory=2000, **keywords):
        """
        2012.7.30
        """
        extraArgumentList = ['-s %s'%(startLineNumber), '-t %s'%(stopLineNumber)]
        if extraArguments:
            extraArgumentList.append(extraArguments)

        job= self.addGenericJob(executable=executable, inputFile=inputFile, outputFile=outputFile, \
                        parentJobLs=parentJobLs, extraDependentInputLs=extraDependentInputLs, \
                        extraOutputLs=[],\
                        transferOutput=transferOutput, \
                        extraArgumentList=extraArgumentList, job_max_memory=job_max_memory, **keywords)
        return job
    def getJVMMemRequirment(self, job_max_memory=5000, minMemory=2000, permSizeFraction=0,
                        MaxPermSizeUpperBound=35000):
        """
        20170502 Java 8 does not support PermSize anymore. set permSizeFraction to 0.
        2013.10.27 handle when job_max_memory is None and minMemory is None.
        #2013.06.07 if a job's virtual memory (1.2X=self.jvmVirtualByPhysicalMemoryRatio, of memory request) exceeds request, it'll abort.
            so set memRequirement accordingly.
        2013.05.01 lower permSizeFraction from 0.4 to 0.2
            minimum for MaxPermSize is now minMemory/2
        2013.04.01
            now job_max_memory = MaxPermSize + mxMemory, unless either is below minMemory
            added argument permSizeFraction, MaxPermSizeUpperBound
        2012.8.2
            job_max_memory could be set by user to lower than minMemory.
            but minMemory makes sure it's never too low.
        """
        if job_max_memory is None:
            job_max_memory = 5000
        if minMemory is None:
            minMemory = 2000
        #20170609 permSizeFraction is set to 0 due to newer Java no longer needs it.
        permSizeFraction = 0
        #MaxPermSize_user = int(job_max_memory*permSizeFraction)
        mxMemory_user = int(job_max_memory*(1-permSizeFraction))
        #MaxPermSize= min(MaxPermSizeUpperBound, max(minMemory/2, MaxPermSize_user))
        #PermSize=MaxPermSize*3/4
        mxMemory = max(minMemory, mxMemory_user)
        msMemory = mxMemory*3/4
        #-XX:+UseGCOverheadLimit
            #Use a policy that limits the proportion of the VM's time that is spent in GC before an OutOfMemory error is thrown. (Introduced in 6.)
        #-XX:-UseGCOverheadLimit would disable the policy.
        memRequirementInStr = "-Xms%sm -Xmx%sm"%(msMemory, mxMemory)	# -XX:PermSize=%sm -XX:MaxPermSize=%sm"%\
                    #, PermSize, MaxPermSize)
        
        memRequirement = int(mxMemory*self.jvmVirtualByPhysicalMemoryRatio)
        #2013.06.07 if a job's virtual memory (1.2X of memory request) exceeds request, it'll abort.
        return PassingData(memRequirementInStr=memRequirementInStr, memRequirement=memRequirement)

    def scaleJobWalltimeOrMemoryBasedOnInput(self, realInputVolume=10, baseInputVolume=4, baseJobPropertyValue=120, \
                                            minJobPropertyValue=120, maxJobPropertyValue=1440):
        """
        2013.04.04
            assume it's integer
            walltime is in minutes.
        """
        walltime = min(max(minJobPropertyValue, float(realInputVolume)/float(baseInputVolume)*baseJobPropertyValue), \
                                    maxJobPropertyValue)	#in minutes
        return PassingData(value=int(walltime))

    def addPlotLDJob(self, executable=None, inputFile=None, inputFileList=None, outputFile=None, \
                    outputFnamePrefix=None,
                    whichColumn=None, whichColumnHeader=None, whichColumnPlotLabel=None, title=None, \
                    logY=None, valueForNonPositiveYValue=-1, \
                    missingDataNotation='-nan',\
                    xColumnPlotLabel=None, chrLengthColumnHeader=None, chrColumnHeader=None, \
                    minChrLength=1000000, xColumnHeader=None, pos2ColumnHeader=None, minNoOfTotal=100, maxNoOfTotal=None,\
                    figureDPI=300, formatString='.', ylim_type=2, samplingRate=0.0001,  need_svg=False, logCount=False, \
                    minDist=None, maxDist=None, movingAverageType=2,\
                    parentJobLs=None, \
                    extraDependentInputLs=None, \
                    extraArgumentList=None, extraArguments=None, transferOutput=True,  job_max_memory=2000, **keywords):
        """
        2012.10.25
            expose argument missingDataNotation, minDist, maxDist
        2012.8.18
            use addAbstractPlotJob()
        2012.8.2 moved from vervet/src/CalculateVCFStatPipeline.py
        2012.8.1
        ('outputFname', 1, ): [None, 'o', 1, 'output file for the figure.'],\
            ('minNoOfTotal', 1, int): [100, 'i', 1, 'minimum no of total variants (denominator of inconsistent rate)'],\
            ('title', 1, ): [None, 't', 1, 'title for the figure.'],\
            ('figureDPI', 1, int): [200, 'f', 1, 'dpi for the output figures (png)'],\
            ('formatString', 1, ): ['.', 'a', 1, 'formatString passed to matplotlib plot'],\
            ('ylim_type', 1, int): [1, 'y', 1, 'y-axis limit type, 1: 0 to max. 2: min to max'],\
            ('samplingRate', 1, float): [0.001, 's', 1, 'how often you include the data'],\
            ('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
            ('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.'],\
            ('whichColumn', 0, int): [3, 'w', 1, 'data from this column (index starting from 0) is plotted as y-axis value'],\
            ('whichColumnHeader', 0, ): ["", 'W', 1, 'column label (in the header) for the data to be plotted as y-axis value, substitute whichColumn'],\
            ('logWhichColumn', 0, int): [0, 'g', 0, 'whether to take -log of the whichColumn'],\
            ('whichColumnPlotLabel', 1, ): ['#SNPs in 100kb window', 'D', 1, 'plot label for data of the whichColumn', ],\
            ('chrLengthColumnHeader', 1, ): ['chrLength', 'c', 1, 'label of the chromosome length column', ],\
            ('chrColumnHeader', 1, ): ['CHR', 'C', 1, 'label of the chromosome column', ],\
            ('minChrLength', 1, int): [1000000, 'm', 1, 'minimum chromosome length for one chromosome to be included', ],\
            ('pos1ColumnLabel', 1, ): ['POS1', 'l', 1, 'label of the 1st position column', ],\
            ('pos2ColumnLabel', 1, ): ['POS2', 'p', 1, 'label of the 2nd position column', ],\
            ('posColumnPlotLabel', 1, ): ['distance', 'x', 1, 'x-axis label in  plot', ],\

        """
        if extraArguments is None:
            extraArguments = ""
        if extraArgumentList is None:
            extraArgumentList = []
        if logCount:
            extraArguments += " --logCount "
        if minChrLength is not None:
            extraArguments += " --minChrLength %s "%(minChrLength)
        if chrLengthColumnHeader:
            extraArgumentList.append("--chrLengthColumnHeader %s"%(chrLengthColumnHeader))
        if chrColumnHeader:
            extraArgumentList.append("--chrColumnHeader %s"%(chrColumnHeader))
        if pos2ColumnHeader:
            extraArgumentList.append(' --pos2ColumnHeader %s '%(pos2ColumnHeader))
        if minDist:
            extraArgumentList.append('--minDist %s'%(minDist))
        if maxDist:
            extraArgumentList.append('--maxDist %s'%(maxDist))
        if maxNoOfTotal:
            extraArgumentList.append("--maxNoOfTotal %s"%(maxNoOfTotal))
        if movingAverageType:
            extraArgumentList.append("--movingAverageType %s"%(movingAverageType))

        return self.addAbstractPlotJob(executable=executable, inputFileList=inputFileList, \
                            inputFile=inputFile, outputFile=outputFile, outputFnamePrefix=outputFnamePrefix, whichColumn=whichColumn, \
                            whichColumnHeader=whichColumnHeader, whichColumnPlotLabel=whichColumnPlotLabel, \
                            logY=logY, valueForNonPositiveYValue=valueForNonPositiveYValue, \
                            missingDataNotation=missingDataNotation,\
                            xColumnHeader=xColumnHeader, xColumnPlotLabel=xColumnPlotLabel, title=title, \
                            minNoOfTotal=minNoOfTotal, \
                            figureDPI=figureDPI, formatString=formatString, ylim_type=ylim_type, samplingRate=samplingRate, need_svg=need_svg, \
                            parentJobLs=parentJobLs, extraDependentInputLs=extraDependentInputLs, \
                            extraArgumentList=extraArgumentList, extraArguments=extraArguments, transferOutput=transferOutput, \
                            job_max_memory=job_max_memory, \
                            **keywords)

    def addPlotVCFtoolsStatJob(self, executable=None, inputFileList=None, outputFnamePrefix=None, \
                            whichColumn=None, whichColumnHeader=None, whichColumnPlotLabel=None, need_svg=False, \
                            logY=0, valueForNonPositiveYValue=-1, \
                            xColumnPlotLabel=None, xColumnHeader=None, chrLengthColumnHeader=None, chrColumnHeader=None, \
                            minChrLength=1000000, minNoOfTotal=100,\
                            figureDPI=300, ylim_type=2, samplingRate=0.0001, logCount=False,\
                            tax_id=60711, sequence_type_id=1, chrOrder=None,\
                            parentJobLs=None, \
                            extraDependentInputLs=None, \
                            extraArguments=None, transferOutput=True, job_max_memory=2000, sshDBTunnel=False, **keywords):
        """
        Examples
            outputFnamePrefix = os.path.join(plotOutputDir, 'noOfMendelErrors_along_chromosome')
            self.addPlotVCFtoolsStatJob(executable=workflow.PlotVCFtoolsStat, \
                                inputFileList=[splitPlinkLMendelFileSNPIDIntoChrPositionJob.output], \
                                outputFnamePrefix=outputFnamePrefix, \
                                whichColumn=None, whichColumnHeader="N", whichColumnPlotLabel="noOfMendelErrors", need_svg=False, \
                                logY=0, valueForNonPositiveYValue=-1, \
                                xColumnPlotLabel="genomePosition", chrLengthColumnHeader=None, chrColumnHeader="Chromosome", \
                                minChrLength=100000, xColumnHeader="Start", minNoOfTotal=50,\
                                figureDPI=100, ylim_type=2, samplingRate=1,\
                                tax_id=self.ref_genome_tax_id, sequence_type_id=self.ref_genome_sequence_type_id, chrOrder=1,\
                                parentJobLs=[splitPlinkLMendelFileSNPIDIntoChrPositionJob, plotOutputDirJob], \
                                extraDependentInputLs=None, \
                                extraArguments=None, transferOutput=True, sshDBTunnel=self.needSSHDBTunnel)


        2013.07.26 added argument tax_id, sequence_type_id, chrOrder
        2013.05.27 remove argument positiveLog, rename logWhichColumn to logY
        2012.10.6 use addGenericDBJob() instead of addGenericJob()
        2012.8.31 add argument positiveLog and valueForNonPositiveYValue
        # whichColumnPlotLabel and xColumnPlotLabel should not contain spaces or ( or ). because they will disrupt shell commandline

        2012.8.2 moved from vervet/src/CalculateVCFStatPipeline.py
        2012.8.1

            ('whichColumn', 0, int): [3, 'w', 1, 'data from this column (index starting from 0) is plotted as y-axis value'],\
            ('whichColumnHeader', 0, ): ["", 'W', 1, 'column label (in the header) for the data to be plotted as y-axis value, substitute whichColumn'],\
            ('logY', 0, int): [0, '', 1, 'value 0: nothing; 1: log(), 2: -log(). replacing self.logWhichColumn.'],\
            ('need_svg', 0, ): [0, 'n', 0, 'whether need svg output', ],\
            ('whichColumnPlotLabel', 1, ): ['#SNPs in 100kb window', 'D', 1, 'plot label for data of the whichColumn', ],\
            ('xColumnPlotLabel', 1, ): ['position', 'x', 1, 'x-axis label (posColumn) in manhattan plot', ],\
            ('chrLengthColumnHeader', 1, ): ['chrLength', 'c', 1, 'label of the chromosome length column', ],\
            ('chrColumnHeader', 1, ): ['CHR', 'C', 1, 'label of the chromosome column', ],\
            ('minChrLength', 1, int): [1000000, 'm', 1, 'minimum chromosome length for one chromosome to be included', ],\
            ('xColumnHeader', 1, ): ['BIN_START', 'l', 1, 'label of the position column, BIN_START for binned vcftools output. POS for others.', ],\
            ('outputFnamePrefix', 0, ): [None, 'O', 1, 'output filename prefix (optional).'],\

                ('minNoOfTotal', 1, int): [100, 'i', 1, 'minimum no of total variants (denominator of inconsistent rate)'],\
                ('title', 1, ): [None, 't', 1, 'title for the figure.'],\
                ('figureDPI', 1, int): [200, 'f', 1, 'dpi for the output figures (png)'],\
                ('formatString', 1, ): ['-', '', 1, 'formatString passed to matplotlib plot'],\
                ('ylim_type', 1, int): [1, 'y', 1, 'y-axis limit type, 1: 0 to max. 2: min to max'],\
                ('samplingRate', 1, float): [0.001, 's', 1, 'how often you include the data'],\
        """
        if extraDependentInputLs is None:
            extraDependentInputLs = []
        if inputFileList:
            extraDependentInputLs.extend(inputFileList)
        extraArgumentList = ["--outputFnamePrefix %s"%outputFnamePrefix, '--minNoOfTotal %s'%(minNoOfTotal), \
                            '--figureDPI %s'%(figureDPI), '--ylim_type %s'%(ylim_type), '--samplingRate %s'%(samplingRate), \
                            '--xColumnHeader %s'%(xColumnHeader)]
        extraOutputLs = [File('%s.png'%(outputFnamePrefix)), File('%s_hist.png'%(outputFnamePrefix))]
        if need_svg:
            extraOutputLs.append(File('%s.svg'%(outputFnamePrefix)))
        key2ObjectForJob = {}
        if minChrLength is not None:
            extraArgumentList.append('--minChrLength %s'%(minChrLength))
        if whichColumnHeader:
            extraArgumentList.append("--whichColumnHeader %s"%(whichColumnHeader))
        if whichColumn:
            extraArgumentList.append("--whichColumn %s"%(whichColumn))
        if logY is not None:
            extraArgumentList.append('--logY %s'%(logY))
        if whichColumnPlotLabel:
            extraArgumentList.append("--whichColumnPlotLabel %s"%(whichColumnPlotLabel))
        if xColumnPlotLabel:
            extraArgumentList.append("--xColumnPlotLabel %s"%(xColumnPlotLabel))
        if chrLengthColumnHeader:
            extraArgumentList.append("--chrLengthColumnHeader %s"%(chrLengthColumnHeader))
        if chrColumnHeader:
            extraArgumentList.append("--chrColumnHeader %s"%(chrColumnHeader))
        if logCount:
            extraArgumentList.append("--logCount")
        if valueForNonPositiveYValue:
            extraArgumentList.append("--valueForNonPositiveYValue %s"%(valueForNonPositiveYValue))
        if sequence_type_id:
            extraArgumentList.append("--sequence_type_id %s"%(sequence_type_id))
        if tax_id:
            extraArgumentList.append("--tax_id %s"%(tax_id))
        if chrOrder is not None:
            extraArgumentList.append("--chrOrder %s"%(chrOrder))

        if extraArguments:
            extraArgumentList.append(extraArguments)
        job= self.addGenericDBJob(executable=executable, inputFile=None, outputFile=None, \
                inputFileList=inputFileList, \
                parentJobLs=parentJobLs, extraDependentInputLs=extraDependentInputLs, \
                extraOutputLs=extraOutputLs,\
                transferOutput=transferOutput, \
                extraArgumentList=extraArgumentList, key2ObjectForJob=key2ObjectForJob, job_max_memory=job_max_memory, \
                sshDBTunnel=sshDBTunnel, objectWithDBArguments=self, **keywords)
        return job

    def addPlotGenomeWideDataJob(self, executable=None, inputFileList=None, inputFile=None,\
                            outputFnamePrefix=None, outputFile=None,\
                            whichColumn=None, whichColumnHeader=None, whichColumnPlotLabel=None, \
                            logX=None, logY=None, valueForNonPositiveYValue=-1, \
                            xScaleLog=None, yScaleLog=None,\
                            missingDataNotation='NA',\
                            xColumnPlotLabel=None, xColumnHeader=None, \
                            xtickInterval=None,\
                            drawCentromere=None, chrColumnHeader=None, \
                            minChrLength=100000, minNoOfTotal=100, maxNoOfTotal=None, \
                            figureDPI=300, formatString=".", ylim_type=2, samplingRate=1, logCount=False, need_svg=False,\
                            tax_id=60711, sequence_type_id=1, chrOrder=None,\
                            inputFileFormat=1, outputFileFormat=None,\
                            parentJobLs=None, extraDependentInputLs=None, \
                            extraArguments=None, extraArgumentList=None, \
                            transferOutput=True, job_max_memory=2000, \
                            objectWithDBGenomeArguments=None, sshDBTunnel=False, \
                            **keywords):

        """
        Examples:
            outputFile = File(os.path.join(plotOutputDir, 'noOfMendelErrors_along_chromosome.png'))
            self.addPlotGenomeWideDataJob(inputFileList=None, \
                            inputFile=splitPlinkLMendelFileSNPIDIntoChrPositionJob.output,\
                            outputFile=outputFile,\
                            whichColumn=None, whichColumnHeader="N", whichColumnPlotLabel="noOfMendelErrors", \
                            logX=None, logY=None, valueForNonPositiveYValue=-1, \
                            xScaleLog=None, yScaleLog=None,\
                            missingDataNotation='NA',\
                            xColumnPlotLabel="genomePosition", xColumnHeader="Start", \
                            xtickInterval=20000000,\
                            chrColumnHeader="Chromosome", \
                            minChrLength=None, minNoOfTotal=None, maxNoOfTotal=None, \
                            figureDPI=300, formatString=".", ylim_type=2, samplingRate=1, logCount=False, need_svg=False,\
                            tax_id=self.ref_genome_tax_id, sequence_type_id=self.ref_genome_sequence_type_id, chrOrder=1,\
                            inputFileFormat=1, outputFileFormat=None,\
                            parentJobLs=[splitPlinkLMendelFileSNPIDIntoChrPositionJob, plotOutputDirJob], \
                            extraDependentInputLs=None, \
                            extraArguments=None, extraArgumentList=None, \
                            transferOutput=True, job_max_memory=1000, sshDBTunnel=self.needSSHDBTunnel)

        2013.07.24
        """
        if extraArgumentList is None:
            extraArgumentList=[]
        if executable is None:
            executable = self.PlotGenomeWideData

        if objectWithDBGenomeArguments is None:
            objectWithDBGenomeArguments = self
        """
        #2013.07.31 replaced by addDBGenomeArgumentsToOneJob() in addGenericJob()
        extraArgumentList.extend(['--genome_drivername=%s'%self.drivername,\
                '--genome_hostname=%s'%self.hostname,\
                '--genome_dbname=%s'%(genome_dbname),\
                '--genome_schema=genome',\
                '--genome_db_user=%s'%(self.db_user),\
                '--genome_db_passwd=%s'%(self.db_passwd)])
        """
        if drawCentromere:
            extraArgumentList.append('--drawCentromere')
        if xtickInterval is not None:
            extraArgumentList.append("--xtickInterval %s"%(xtickInterval))
        if chrColumnHeader is not None:
            extraArgumentList.append('--chromosomeHeader %s'%(chrColumnHeader))
        if tax_id:
            extraArgumentList.append("--tax_id %s"%(tax_id))
        if sequence_type_id:
            extraArgumentList.append('--sequence_type_id %s'%(sequence_type_id))
        if chrOrder is not None:
            extraArgumentList.append("--chrOrder %s"%(chrOrder))
        job = self.addAbstractPlotJob(executable=executable, \
            inputFileList=inputFileList, \
            inputFile=inputFile, outputFile=outputFile, \
            outputFnamePrefix=outputFnamePrefix, whichColumn=whichColumn, whichColumnHeader=whichColumnHeader, \
            whichColumnPlotLabel=whichColumnPlotLabel, \
            logX=logX, logY=logY, valueForNonPositiveYValue=valueForNonPositiveYValue, \
            xScaleLog=xScaleLog, yScaleLog=yScaleLog,\
            missingDataNotation=missingDataNotation,\
            xColumnHeader=xColumnHeader, xColumnPlotLabel=xColumnPlotLabel, \
            minNoOfTotal=minNoOfTotal, maxNoOfTotal=maxNoOfTotal,\
            figureDPI=figureDPI, formatString=formatString, ylim_type=ylim_type, \
            samplingRate=samplingRate, need_svg=need_svg, \
            inputFileFormat=inputFileFormat, outputFileFormat=outputFileFormat,\
            parentJobLs=parentJobLs, \
            extraDependentInputLs=extraDependentInputLs, \
            extraArgumentList=extraArgumentList, \
            extraArguments=extraArguments, transferOutput=transferOutput,  job_max_memory=job_max_memory, \
            sshDBTunnel=sshDBTunnel, objectWithDBGenomeArguments=objectWithDBGenomeArguments, \
            **keywords)
        return job

    def addAbstractPlotJob(self, executable=None, inputFileList=None, inputFile=None, outputFile=None, \
                    outputFnamePrefix=None, whichColumn=None, whichColumnHeader=None, whichColumnPlotLabel=None, \
                    logX=None, logY=None, valueForNonPositiveYValue=-1, \
                    xScaleLog=0, yScaleLog=0, \
                    missingDataNotation='NA',\
                    xColumnHeader=None, xColumnPlotLabel=None, title=None, \
                    minNoOfTotal=100, maxNoOfTotal=None,\
                    figureDPI=300, formatString='.', markerSize=None, \
                    ylim_type=2, samplingRate=0.001, legendType=None,\
                    need_svg=False, \
                    inputFileFormat=None, outputFileFormat=None,\
                    parentJob=None, parentJobLs=None, \
                    extraDependentInputLs=None, extraOutputLs=None, \
                    extraArgumentList=None, extraArguments=None, transferOutput=True,  job_max_memory=2000, \
                    sshDBTunnel=False, key2ObjectForJob=None, \
                    objectWithDBArguments=None, **keywords):
        """
        2013.08.28 added argument markerSize
        2013.07.16 added argument legendType
        2012.12.3 added argument title, logX, logY
        2012.10.16 added argument sshDBTunnel, objectWithDBArguments
        2012.8.31 add argument missingDataNotation
        #no spaces or parenthesis or any other shell-vulnerable letters in the x or y axis labels (whichColumnPlotLabel, xColumnPlotLabel)
        2012.8.2 (check AbstractMatrixFileWalker.py or AbstractPlot.py for updated arguments)
            ('outputFname', 0, ): [None, 'o', 1, 'output file for the figure.'],\
            ('minNoOfTotal', 1, int): [100, 'M', 1, 'minimum no of total variants (denominator of inconsistent rate)'],\
            ('title', 0, ): [None, 't', 1, 'title for the figure.'],\
            ('figureDPI', 1, int): [200, 'f', 1, 'dpi for the output figures (png)'],\
            ('formatString', 1, ): ['-', '', 1, 'formatString passed to matplotlib plot'],\
            ('ylim_type', 1, int): [1, 'y', 1, 'y-axis limit type, 1: whatever matplotlib decides. 2: min to max'],\
            ('samplingRate', 1, float): [1, 's', 1, 'how often you include the data, a probability between 0 and 1.'],\
            ('whichColumn', 0, int): [3, 'w', 1, 'data from this column (index starting from 0) is plotted as y-axis value'],\
            ('whichColumnHeader', 0, ): ["", 'W', 1, 'column header for the data to be plotted as y-axis value, substitute whichColumn'],\
            ('whichColumnPlotLabel', 0, ): ['', 'D', 1, 'plot label for data of the whichColumn', ],\
            ('logY', 0, int): [0, '', 1, 'value 0: nothing; 1: log(), 2: -log(). replacing self.logWhichColumn.'],\
            ('valueForNonPositiveYValue', 1, float): [50, '', 1, 'if the whichColumn value is not postive and logWhichColumn is on,\
                    what yValue should be.'],\
            ('xColumnHeader', 1, ): ['', 'l', 1, 'header of the x-axis data column, ' ],\
            ('xColumnPlotLabel', 0, ): ['', 'x', 1, 'x-axis label (posColumn) in manhattan plot', ],\
            ('need_svg', 0, ): [0, 'n', 0, 'whether need svg output', ],\
            ('legendType', 0, int): [0, '', 1, '0: no legend; 1: legend'], \

            inputFileFormat   1: csv-like plain text file; 2: YHPyTables.YHFile; 3: HDF5MatrixFile; . "1"(default)

        """
        if extraOutputLs is None:
            extraOutputLs = []
        if key2ObjectForJob is None:
            key2ObjectForJob = {}
        if executable is None:
            executable = self.AbstractPlot
        if extraDependentInputLs is None:
            extraDependentInputLs = []
        if inputFileList:
            extraDependentInputLs.extend(inputFileList)
        if extraArgumentList is None:
            extraArgumentList = []

        if outputFnamePrefix:
            extraArgumentList.append('--outputFnamePrefix %s'%(outputFnamePrefix))
            if outputFile is None:
                extraOutputLs.append(File('%s.png'%(outputFnamePrefix)))
                if need_svg:
                    extraOutputLs.append(File('%s.svg'%(outputFnamePrefix)))
        if minNoOfTotal is not None:
            extraArgumentList.append('--minNoOfTotal %s'%(minNoOfTotal))
        if maxNoOfTotal:
            extraArgumentList.append("--maxNoOfTotal %s"%(maxNoOfTotal))
        if figureDPI:
            extraArgumentList.append('--figureDPI %s'%(figureDPI))
        if formatString:
            extraArgumentList.append('--formatString %s'%(formatString))
        if markerSize:
            extraArgumentList.append("--markerSize %s"%(markerSize))
        if ylim_type:
            extraArgumentList.append('--ylim_type %s'%(ylim_type))
        if samplingRate is not None:
            extraArgumentList.append('--samplingRate %s'%(samplingRate))
        if legendType!=None:
            extraArgumentList.append("--legendType %s"%(legendType))

        if xColumnHeader:
            extraArgumentList.append('--xColumnHeader %s'%(xColumnHeader))
        if xColumnPlotLabel:
            extraArgumentList.append("--xColumnPlotLabel %s"%(xColumnPlotLabel))
        if whichColumnHeader:
            extraArgumentList.append("--whichColumnHeader %s"%(whichColumnHeader))
        if whichColumn:
            extraArgumentList.append("--whichColumn %s"%(whichColumn))
        if whichColumnPlotLabel:
            extraArgumentList.append("--whichColumnPlotLabel %s"%(whichColumnPlotLabel))
        if title:
            extraArgumentList.append("--title %s"%(title))
        if logX:
            extraArgumentList.append("--logX %s"%(logX))
        if logY:
            extraArgumentList.append('--logY %s'%(logY))
        if xScaleLog:
            extraArgumentList.append("--xScaleLog %s"%(xScaleLog))
        if yScaleLog:
            extraArgumentList.append("--yScaleLog %s"%(yScaleLog))

        if valueForNonPositiveYValue:
            extraArgumentList.append("--valueForNonPositiveYValue %s"%(valueForNonPositiveYValue))
        if inputFileFormat:
            extraArgumentList.append("--inputFileFormat %s"%(inputFileFormat))
        if outputFileFormat:
            extraArgumentList.append("--outputFileFormat %s"%(outputFileFormat))
        if need_svg:
            extraArgumentList.append('--need_svg')
            if not outputFnamePrefix:
                outputFnamePrefix = os.path.splitext(outputFile.name)[0]	#2012.8.20 bugfix.
            extraOutputLs.append(File('%s.svg'%(outputFnamePrefix)))
        if extraArguments:
            extraArgumentList.append(extraArguments)

        job = self.addGenericJob(executable=executable, inputFile=inputFile, outputFile=outputFile, \
                inputFileList = inputFileList,\
                parentJob=parentJob, parentJobLs=parentJobLs, \
                extraDependentInputLs=extraDependentInputLs, \
                extraOutputLs=extraOutputLs, transferOutput=transferOutput, \
                extraArgumentList=extraArgumentList, key2ObjectForJob=key2ObjectForJob, job_max_memory=job_max_memory, \
                sshDBTunnel=sshDBTunnel, objectWithDBArguments=objectWithDBArguments, **keywords)
        return job

    def addAbstractMatrixFileWalkerJob(self, executable=None, inputFileList=None, inputFile=None, outputFile=None, \
                    outputFnamePrefix=None, whichColumn=None, whichColumnHeader=None, \
                    logY=None, valueForNonPositiveYValue=-1, \
                    minNoOfTotal=10,\
                    samplingRate=1, \
                    inputFileFormat=None, outputFileFormat=None,\
                    parentJob=None, parentJobLs=None, extraOutputLs=None, \
                    extraDependentInputLs=None, extraArgumentList=None, \
                    extraArguments=None, transferOutput=True,  job_max_memory=2000, sshDBTunnel=False, \
                    objectWithDBArguments=None, **keywords):
        """
        2012.11.25 more arguments, logY, inputFileFormat, outputFileFormat
        2012.10.16 added argument sshDBTunnel, objectWithDBArguments
        2012.10.15 added extraArgumentList, parentJob
        2012.8.15
            ('outputFname', 0, ): [None, 'o', 1, 'output file for the figure.'],\
            ('minNoOfTotal', 1, int): [100, 'i', 1, 'minimum no of total variants (denominator of inconsistent rate)'],\
            ('samplingRate', 1, float): [1, 's', 1, 'how often you include the data, a probability between 0 and 1.'],\
            ('whichColumn', 0, int): [3, 'w', 1, 'data from this column (index starting from 0) is plotted as y-axis value'],\
            ('whichColumnHeader', 0, ): ["", 'W', 1, 'column header for the data to be plotted as y-axis value, substitute whichColumn'],\
            ('whichColumnPlotLabel', 0, ): ['', 'D', 1, 'plot label for data of the whichColumn', ],\
            ('logWhichColumn', 0, int): [0, 'g', 0, 'whether to take -log of the whichColumn'],\
            ('positiveLog', 0, int): [0, 'p', 0, 'toggle to take log, rather than -log(), \
                only effective when logWhichColumn is toggled. '],\
            ('valueForNonPositiveYValue', 1, float): [50, '', 1, 'if the whichColumn value is not postive and logWhichColumn is on,\
                    what yValue should be.'],\

        """
        return self.addAbstractPlotJob(executable=executable, inputFileList=inputFileList, \
                            inputFile=inputFile, outputFile=outputFile, outputFnamePrefix=outputFnamePrefix, whichColumn=whichColumn, \
                            whichColumnHeader=whichColumnHeader, whichColumnPlotLabel=None, \
                            logY=logY, \
                            valueForNonPositiveYValue=valueForNonPositiveYValue, \
                            missingDataNotation=None,\
                            xColumnHeader=None, xColumnPlotLabel=None, \
                            minNoOfTotal=minNoOfTotal, \
                            figureDPI=None, formatString=None, ylim_type=None, samplingRate=samplingRate, need_svg=False, \
                            parentJob=parentJob, parentJobLs=parentJobLs, \
                            extraOutputLs=extraOutputLs, extraDependentInputLs=extraDependentInputLs, \
                            extraArgumentList=extraArgumentList,\
                            extraArguments=extraArguments, transferOutput=transferOutput, job_max_memory=job_max_memory, \
                            sshDBTunnel=sshDBTunnel, objectWithDBArguments=objectWithDBArguments, **keywords)

    def addAbstractGenomeFileWalkerJob(self, executable=None, inputFileList=None, inputFile=None, outputFile=None, \
                    outputFnamePrefix=None, whichColumn=None, whichColumnHeader=None, \
                    logY=None, valueForNonPositiveYValue=-1, \
                    minNoOfTotal=10,\
                    samplingRate=1, \
                    chrColumnHeader=None, \
                    tax_id=60711, sequence_type_id=1, chrOrder=None,\
                    positionHeader=None,\
                    inputFileFormat=None, outputFileFormat=None,\
                    parentJob=None, parentJobLs=None, \
                    extraDependentInputLs=None, extraArgumentList=None, \
                    extraArguments=None, transferOutput=True,  job_max_memory=2000, sshDBTunnel=False, \
                    objectWithDBGenomeArguments=None, **keywords):
        """
        2013.07.31

        """
        if extraArgumentList is None:
            extraArgumentList=[]

        if objectWithDBGenomeArguments is None:
            objectWithDBGenomeArguments = self
        if chrColumnHeader is not None:
            extraArgumentList.append('--chromosomeHeader %s'%(chrColumnHeader))
        if tax_id:
            extraArgumentList.append("--tax_id %s"%(tax_id))
        if sequence_type_id:
            extraArgumentList.append('--sequence_type_id %s'%(sequence_type_id))
        if chrOrder is not None:
            extraArgumentList.append("--chrOrder %s"%(chrOrder))
        if positionHeader is not None:
            extraArgumentList.append('--positionHeader %s'%(positionHeader))
        return self.addAbstractPlotJob(executable=executable, inputFileList=inputFileList, \
                            inputFile=inputFile, outputFile=outputFile, outputFnamePrefix=outputFnamePrefix, whichColumn=whichColumn, \
                            whichColumnHeader=whichColumnHeader, whichColumnPlotLabel=None, \
                            logY=logY, \
                            valueForNonPositiveYValue=valueForNonPositiveYValue, \
                            missingDataNotation=None,\
                            xColumnHeader=None, xColumnPlotLabel=None, \
                            minNoOfTotal=minNoOfTotal, \
                            figureDPI=None, formatString=None, ylim_type=None, samplingRate=samplingRate, need_svg=False, \
                            parentJob=parentJob, parentJobLs=parentJobLs, extraDependentInputLs=extraDependentInputLs, \
                            extraArgumentList=extraArgumentList,\
                            extraArguments=extraArguments, transferOutput=transferOutput, job_max_memory=job_max_memory, \
                            sshDBTunnel=sshDBTunnel, objectWithDBGenomeArguments=objectWithDBGenomeArguments, \
                            **keywords)


    def addDrawHistogramJob(self, executable=None, inputFileList=None, inputFile=None, outputFile=None, \
                    outputFnamePrefix=None, whichColumn=None, whichColumnHeader=None, whichColumnPlotLabel=None, \
                    xScaleLog=0, yScaleLog=0, \
                    logY=None, valueForNonPositiveYValue=-1, missingDataNotation='NA', title=None, \
                    minNoOfTotal=10,\
                    figureDPI=100, formatString='.', ylim_type=2, samplingRate=0.001, need_svg=False, legendType=None, \
                    logCount=False, inputFileFormat=None, \
                    parentJobLs=None, \
                    extraDependentInputLs=None, \
                    extraArguments=None, transferOutput=True,  job_max_memory=2000, **keywords):
        """
        #no spaces or parenthesis or any other shell-vulnerable letters in the x or y axis labels (whichColumnPlotLabel, xColumnPlotLabel)
        2013.08.15 added argument xScaleLog, yScaleLog, legendType
        2012.8.2
            ('outputFname', 0, ): [None, 'o', 1, 'output file for the figure.'],\
            ('minNoOfTotal', 1, int): [100, 'i', 1, 'minimum no of total variants (denominator of inconsistent rate)'],\
            ('title', 0, ): [None, 't', 1, 'title for the figure.'],\
            ('figureDPI', 1, int): [200, 'f', 1, 'dpi for the output figures (png)'],\
            ('formatString', 1, ): ['-', '', 1, 'formatString passed to matplotlib plot'],\
            ('ylim_type', 1, int): [1, 'y', 1, 'y-axis limit type, 1: 0 to max. 2: min to max'],\
            ('samplingRate', 1, float): [1, 's', 1, 'how often you include the data, a probability between 0 and 1.'],\
            ('whichColumn', 0, int): [3, 'w', 1, 'data from this column (index starting from 0) is plotted as y-axis value'],\
            ('whichColumnHeader', 0, ): ["", 'W', 1, 'column header for the data to be plotted as y-axis value, substitute whichColumn'],\
            ('whichColumnPlotLabel', 0, ): ['', 'D', 1, 'plot label for data of the whichColumn', ],\
            ('logWhichColumn', 0, int): [0, 'g', 0, 'whether to take -log of the whichColumn'],\
            ('positiveLog', 0, int): [0, 'p', 0, 'toggle to take log, rather than -log(), \
                only effective when logWhichColumn is toggled. '],\
            ('valueForNonPositiveYValue', 1, float): [50, '', 1, 'if the whichColumn value is not postive and logWhichColumn is on,\
                    what yValue should be.'],\
            ('need_svg', 0, ): [0, 'n', 0, 'whether need svg output', ],\

        """
        if extraArguments is None:
            extraArguments = ""
        if logCount:
            extraArguments += " --logCount "
        return self.addAbstractPlotJob(executable=executable, inputFileList=inputFileList, \
                            inputFile=inputFile, outputFile=outputFile, outputFnamePrefix=outputFnamePrefix, whichColumn=whichColumn, \
                            whichColumnHeader=whichColumnHeader, whichColumnPlotLabel=whichColumnPlotLabel, \
                            xScaleLog=xScaleLog, yScaleLog=yScaleLog,\
                            logY=logY, valueForNonPositiveYValue=valueForNonPositiveYValue, \
                            missingDataNotation=missingDataNotation,\
                            xColumnHeader=None, xColumnPlotLabel=None, title=title, \
                            minNoOfTotal=minNoOfTotal, \
                            figureDPI=figureDPI, formatString=formatString, ylim_type=ylim_type, \
                            samplingRate=samplingRate, need_svg=need_svg, legendType=legendType, \
                            inputFileFormat=inputFileFormat,\
                            parentJobLs=parentJobLs, extraDependentInputLs=extraDependentInputLs, \
                            extraArguments=extraArguments, transferOutput=transferOutput, job_max_memory=job_max_memory, \
                            **keywords)

    def addDraw2DHistogramOfMatrixJob(self, executable=None, inputFileList=None, inputFile=None, outputFile=None, \
                outputFnamePrefix=None, whichColumn=None, whichColumnHeader=None, whichColumnPlotLabel=None, \
                logX=False, logY=False, logZ=False, valueForNonPositiveYValue=-1, \
                missingDataNotation='NA',\
                xColumnHeader=None, xColumnPlotLabel=None, \
                minNoOfTotal=100,\
                figureDPI=300, formatString='.', samplingRate=0.001, need_svg=False, \
                inputFileFormat=None, outputFileFormat=None,\
                zColumnHeader=None, \
                parentJobLs=None, \
                extraDependentInputLs=None, \
                extraArgumentList=None, extraArguments=None, transferOutput=True,  job_max_memory=2000, **keywords):
        """
        2013.2.8 added argument inputFileFormat
        2013.2.7 executable could be None, default is self.Draw2DHistogramOfMatrix
        2012.11.28 change logX, logY, logZ
        2012.10.7

        """
        if extraArgumentList is None:
            extraArgumentList = []
        if zColumnHeader:
            extraArgumentList.append("--zColumnHeader %s"%(zColumnHeader))
        if logZ:
            extraArgumentList.append("--logZ %s"%(logZ))
        if executable is None:
            executable = self.Draw2DHistogramOfMatrix
        return self.addAbstractPlotJob(executable=executable, inputFileList=inputFileList, \
                            inputFile=inputFile, outputFile=outputFile, outputFnamePrefix=outputFnamePrefix, whichColumn=whichColumn, \
                            whichColumnHeader=whichColumnHeader, whichColumnPlotLabel=whichColumnPlotLabel, \
                            logX=logX, logY=logY, valueForNonPositiveYValue=valueForNonPositiveYValue, \
                            missingDataNotation=missingDataNotation,\
                            xColumnHeader=xColumnHeader, xColumnPlotLabel=xColumnPlotLabel, \
                            minNoOfTotal=minNoOfTotal, \
                            figureDPI=figureDPI, formatString=formatString, ylim_type=None, \
                            samplingRate=samplingRate, need_svg=need_svg, \
                            inputFileFormat=inputFileFormat, outputFileFormat=outputFileFormat,\
                            parentJobLs=parentJobLs, extraDependentInputLs=extraDependentInputLs, \
                            extraArguments=extraArguments, transferOutput=transferOutput, job_max_memory=job_max_memory, \
                            **keywords)

    def setupMoreOutputAccordingToSuffixAndNameTupleList(self, outputFnamePrefix=None, suffixAndNameTupleList=None, extraOutputLs=None, key2ObjectForJob=None):
        """
        2012.8.16
            split from addPlinkJob()
        """
        for suffixNameTuple in suffixAndNameTupleList:
            if len(suffixNameTuple)==1:
                suffix = suffixNameTuple[0]
                name = suffix[1:].replace('.', '_')	#replace dot with underscore. as dot is used to access method/attribute of python object
                # i.e. ".prune.in" is accessible as job.prune_inFile
            elif len(suffixNameTuple)>=2:
                suffix, name = suffixNameTuple[:2]
            outputFile = File('%s%s'%(outputFnamePrefix, suffix))
            extraOutputLs.append(outputFile)
            key2ObjectForJob['%sFile'%(name)] = outputFile

if __name__ == '__main__':
    from argparse import ArgumentParser
    ap = ArgumentParser()
    ap.add_argument("-i", "--input_file", type=str, required=True,
                    help="the path to the input file.")
    ap.add_argument("-o", "--output_file", type=str, required=True,
                    help="the path to the output file.")
    ap.add_argument("-s", "--source_code_dir", type=str, 
            default=os.path.expanduser('~/src/mygit/'), 
            help="the path to the source code dir. (default: %(default)s)")
    args = ap.parse_args()

    main_class = AbstractWorkflow
    instance = AbstractWorkflow(args.arguments, args.output_file)
    instance.run()
