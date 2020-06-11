#!/usr/bin/env python3
"""
2012.5.23
    a common class for other pymodule workflows
"""
import sys, os, math
from palos import ProcessOptions, getListOutOfStr, PassingData, utils
from pegaflow.Workflow import Workflow
from pegaflow.DAX3 import Executable, File, PFN, Profile
from pegaflow.DAX3 import Namespace, Link, Use, Job, Dependency

#path to the source code directory
src_dir = os.path.dirname(os.path.abspath(__file__))

class AbstractWorkflow(Workflow):
    __doc__ = __doc__
    db_option_dict = {
        ('drivername', 1,):['postgresql', 'v', 1, \
            'which type of database? mysql or postgresql', ],\
        ('hostname', 1, ): ['localhost', 'z', 1, 'hostname of the db server', ],\
        ('dbname', 1, ): ['', 'd', 1, 'database name', ],\
        ('schema', 0, ): ['public', 'k', 1, 'database schema name', ],\
        ('db_user', 1, ): [None, 'u', 1, 'database username', ],\
        ('db_passwd', 1, ): [None, 'p', 1, 'database password', ],\
        ('port', 0, ):[None, '', 1, 'database port number'],\
        ('commit', 0, ):[None, '', 0, 
            'commit database transaction if there is db transaction'],
        ("data_dir", 0, ): ["", 't', 1, \
            'the base directory where all db-affiliated files are stored. '
            'If not given, use the default stored in db.'],\
        ("local_data_dir", 0, ): ["", 'D', 1, 
            'A folder contains the same files as data_dir but accessible locally. '
            'If not given, use the default stored in db (db.data_dir). '
            'This argument is used to find all input files available. '
            'Set it differently from data_dir only if you are '
            'generating a workflow on a computer that '
            'has a different path to data_dir than your running site.'],
        }
    option_default_dict = {
        ("home_path", 1, ): [os.path.expanduser("~"), 'e', 1, \
            'path to the home directory on the working nodes'],\
        ("pymodulePath", 1, ): ["src/pymodule", '', 1, 'path to the pymodule folder'],\
        ("thisModulePath", 1, ): ["", '', 1, 
            'path of the module that owns this program. '
            'used to add executables from this module.'],\
        ("javaPath", 1, ): ["bin/java", 'J', 1, 'the path to java binary'],\
        ("site_handler", 1, ): ["condor", 'l', 1, 
            'which computing site to run the jobs. check your Pegasus setup.'],\
        ("input_site_handler", 0, ): ["None", 'j', 1, 
            'the site that stores all the input files: local, condor, hoffman2. '
            'If site_handler is condor, this must be condor and files will be symlinked. '
            'If site_handler is hoffman2, input_site_handler=local induces file transfer '
            'and input_site_handler=hoffman2 induces symlink. '
            'Default is the same as the site_handler.'],\
        ('cluster_size', 1, int):[30, 'C', 1, 
            'For short jobs that will be clustered, how many of them should be clustered int one'],\
        ('pegasusFolderName', 0, ): ['folder', 'F', 1, \
            'The folder relative to pegasus scratch working dir to contain input & output. '
            'It will be created during the pegasus staging process. '
            'It is useful to separate multiple workflows. '
            'If empty, everything is in the pegasus root.', ],\
        ('inputSuffixList', 0, ): [None, '', 1, \
            'coma-separated list of input file suffices. If None, any suffix. '
            'Suffix include the dot, (i.e. .tsv). '
            'Typical zip suffices are excluded (.gz, .bz2, .zip, .bz).'],\
        ('output_path', 1, ): [None, 'o', 1, 'xml workflow output file'],\
        ("tmpDir", 1, ): ["/tmp/", '', 1, \
            'for MarkDuplicates.jar, etc., default is /tmp/ but sometimes it is too small'],\
        ('max_walltime', 1, int):[4320, '', 1, \
            'maximum wall time any job could have, in minutes. 20160=2 weeks.\n'
            'used in addGenericJob().'],\
        ('jvmVirtualByPhysicalMemoryRatio', 1, float):[1.2, '', 1, 
            "if a job's virtual memory (usually 1.2X of JVM resident memory) exceeds request, "
            "it will be killed on hoffman2. Hence this argument"],\
        ('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
        ('needSSHDBTunnel', 0, int):[0, 'H', 0, \
            'DB-interacting jobs need a ssh tunnel (running on cluster behind firewall).'],\
        ('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']
        }
        #("pegasusCleanupPath", 1, ): ["bin/pegasus/bin/pegasus-cleanup", '', 1, \
        #    'path to pegasus-cleanup executable, it will be registered to run in local '
        #    'universe of condor pool (rather than the vanilla universe)'],\
        #("pegasusTransferPath", 1, ): ["bin/pegasus/bin/pegasus-transfer", '', 1, \
        #    'path to pegasus-transfer executable, it will be registered to run in local '
        #     'universe of condor pool (rather than the vanilla universe)'],\
        #('bamListFname', 1, ): ['/tmp/bamFileList.txt', 'L', 1, 
        #     'The file contains path to each bam file, one file per line.'],\

    pathToInsertHomePathList = ['javaPath', 'pymodulePath', 'plinkPath', "thisModulePath"]

    def __init__(self,
        input_path=None,
        inputSuffixList=None, 
        pegasusFolderName='folder',
        output_path=None,

        tmpDir='/tmp/', max_walltime=4320,

        home_path=None,
        pymodulePath="src/pymodule", thisModulePath=None,
        javaPath=None, 
        jvmVirtualByPhysicalMemoryRatio=1.2,
        
        site_handler='condor', input_site_handler='condor',
        cluster_size=30,
        needSSHDBTunnel=False, commit=False,
        debug=False, report=False):
        """
        """
        # Set before ParentClass.__init__()
        if not home_path:
            home_path = os.path.expanduser("~")
        self.pymodulePath = pymodulePath
        self.thisModulePath = thisModulePath
        Workflow.__init__(self,
            input_path=input_path,
            inputSuffixList=inputSuffixList,
            pegasusFolderName=pegasusFolderName,
            output_path=output_path,

            tmpDir=tmpDir, max_walltime=max_walltime, 
            
            home_path=home_path,
            javaPath=javaPath,
            jvmVirtualByPhysicalMemoryRatio=jvmVirtualByPhysicalMemoryRatio,
            
            site_handler=site_handler,
            input_site_handler=input_site_handler,
            cluster_size=cluster_size,

            needSSHDBTunnel=needSSHDBTunnel, commit=commit,
            debug=debug, report=report)

    def constructJobDataFromJob(self, job=None):
        """
        added vcfFile and tbi_F in structure
        """
        if hasattr(job, 'output') and job.output.name.find('.vcf')>=0:
            vcfFile = job.output
            tbi_F = getattr(job, 'tbi_F', None)
        else:
            vcfFile = None
            tbi_F = None
        return PassingData(job=job, jobLs=[job], file=job.output, \
            fileLs=job.outputLs, vcfFile=vcfFile, tbi_F=tbi_F)

    def registerExecutables(self):
        """
        """
        Workflow.registerExecutables(self)
        
        if os.path.isfile("/usr/bin/convert"):
            #convert, an image swissknife program, part of imagemagick
            self.registerOneExecutable(path="/usr/bin/convert",
                name='convertImage', clusterSizeMultiplier=1)
        self.registerOneExecutable(path=os.path.join(self.pymodulePath, 
            "mapper/extractor/SelectLineBlockFromFile.py"), 
            clusterSizeMultiplier=1)
        self.registerOneExecutable(path=os.path.join(self.pymodulePath,
            "plot/AbstractPlot.py"), clusterSizeMultiplier=1)
        self.registerOneExecutable(path=os.path.join(self.pymodulePath,
            "plot/PlotYAsBar.py"), clusterSizeMultiplier=1)
        self.registerOneExecutable(path=os.path.join(self.pymodulePath,
            "plot/DrawHistogram.py"), clusterSizeMultiplier=1)

        self.registerOneExecutable(path=os.path.join(self.pymodulePath,
            "plot/DrawMatrix.py"), clusterSizeMultiplier=1)
        self.registerOneExecutable(path=os.path.join(self.pymodulePath,
            "plot/Draw2DHistogramOfMatrix.py"), clusterSizeMultiplier=1)
        self.registerOneExecutable(path=os.path.join(self.pymodulePath,
            "statistics/SampleRows.py"), 
            clusterSizeMultiplier=1)
        #2013.2.11 all reducers
        self.registerOneExecutable(path=os.path.join(self.pymodulePath,
            "statistics/EstimateOutliersIn2DData.py"), \
            clusterSizeMultiplier=0)
        self.registerOneExecutable(path=os.path.join(self.pymodulePath,
            'reducer/MergeSameHeaderTablesIntoOne.py'), \
            name='mergeSameHeaderTablesIntoOne',
            clusterSizeMultiplier=0)
        self.registerOneExecutable(path=os.path.join(self.pymodulePath,
            'reducer/MergeSameHeaderTablesIntoOne.py'), \
            clusterSizeMultiplier=0)
        self.registerOneExecutable(path=os.path.join(self.pymodulePath,
            'reducer/ReduceMatrixByAverageColumnsWithSameKey.py'), \
            clusterSizeMultiplier=0)
        self.registerOneExecutable(path=os.path.join(self.pymodulePath,
            'reducer/ReduceMatrixByChosenColumn.py'), \
            clusterSizeMultiplier=0)
        self.registerOneExecutable(path=os.path.join(self.pymodulePath,
            'reducer/ReduceMatrixByMergeColumnsWithSameKey.py'), \
            clusterSizeMultiplier=0)
        self.registerOneExecutable(path=os.path.join(self.pymodulePath,
            'reducer/ReduceMatrixBySumSameKeyColsAndThenDivide.py'),
            clusterSizeMultiplier=0)
        """
        # 2013.05.20 DISABLE this
        if self.site_handler=='hcondor' and self.input_site_handler=='hcondor':
            #2013.04.19 to make pegasus cleanup run on local universe of condor pool
            # only enable this on hcondor because
            #	1) its filesystem is very slow and these cleanup & transfer
            #       jobs take forever.
            #	2) workers in vanilla universe expire after certain time.
            #	3) it does not run on ycondor local universe somehow.
            #    pegasus keeps submitting but no condor jobs in the queue.
            # this works because in most of my cases, vanilla universe and
            #    local universe share the same underlying filesystem.
            cleanupExecutable = self.registerOneExecutable(
                path=self.pegasusCleanupPath, name='cleanup',
                clusterSizeMultiplier=0, noVersion=True)
            condorUniverseProfile = Profile(Namespace.CONDOR, key="universe",
                value="local")
            if cleanupExecutable.hasProfile(condorUniverseProfile):
                cleanupExecutable.removeProfile(condorUniverseProfile)
            cleanupExecutable.addProfile(condorUniverseProfile)

            transferExecutable = self.registerOneExecutable(
                path=self.pegasusTransferPath, name='transfer',
                clusterSizeMultiplier=0, noVersion=True)
            condorUniverseProfile = Profile(Namespace.CONDOR, key="universe",
                value="local")
            if transferExecutable.hasProfile(condorUniverseProfile):
                transferExecutable.removeProfile(condorUniverseProfile)
            transferExecutable.addProfile(condorUniverseProfile)
        """

    def addStatMergeJob(self, statMergeProgram=None, outputF=None,
        extraArguments=None, extraArgumentList=None,
        parentJobLs=None, extraOutputLs=None,
        extraDependentInputLs=None, transferOutput=True,
        key2ObjectForJob=None, job_max_memory=1000, **keywords):
        """
        moved from CalculateVCFStatPipeline.py
        """
        if extraDependentInputLs is None:
            extraDependentInputLs = []
        if extraArgumentList is None:
            extraArgumentList = []
        if key2ObjectForJob is None:
            key2ObjectForJob = {}
        if extraArguments:
            extraArgumentList.append(extraArguments)
        job= self.addGenericJob(executable=statMergeProgram,
            inputFile=None, inputArgumentOption='-i',
            outputFile=outputF, outputArgumentOption='-o', 
            extraArgumentList=extraArgumentList,
            parentJobLs=parentJobLs,
            extraDependentInputLs=extraDependentInputLs,
            extraOutputLs=extraOutputLs,
            transferOutput=transferOutput,
            key2ObjectForJob=key2ObjectForJob,
            job_max_memory=job_max_memory, **keywords)
        return job

    def addConvertImageJob(self, inputFile=None, inputArgumentOption=None,
        outputFile=None, outputArgumentOption=None, density=None,
        resizeDimension=None,
        frontArgumentList=None,
        extraArguments=None, extraArgumentList=None,
        parentJobLs=None,
        extraDependentInputLs=None, extraOutputLs=None,
        transferOutput=False,
        job_max_memory=200,
        key2ObjectForJob=None, **keywords):
        """
        2013.2.7 use imagemagick's convert to convert images. examples:
            plotOutputFile = File('%s.eps'%(plotOutputFnamePrefix))
            plotPNGOutputFile = File('%s.png'%(plotOutputFnamePrefix))
            #change dpi to 300
            self.addConvertImageJob(inputFile=plotOutputFile, \
                outputFile=plotPNGOutputFile, density=300, \
                resizeDimension=None, \
                parentJobLs=[psmc_plotJob], extraDependentInputLs=None, 
                extraOutputLs=None, transferOutput=True, \
                extraArguments=None, frontArgumentList=None, job_max_memory=500)

            #resize by demanding the width = 1800, height will scale accordingly
            self.addConvertImageJob(inputFile=plotOutputFile, \
                outputFile=plotPNGOutputFile, density=None, \
                resizeDimension=1800, \
                parentJobLs=[psmc_plotJob], extraDependentInputLs=None, \
                extraOutputLs=None, transferOutput=True, \
                extraArguments=None, frontArgumentList=None, job_max_memory=500)
            #resize by demanding the dimension=1800X900
            self.addConvertImageJob(inputFile=plotOutputFile, \
                outputFile=plotPNGOutputFile, density=None, \
                resizeDimension='1800X900', \
                parentJobLs=[psmc_plotJob], extraDependentInputLs=None, \
                extraOutputLs=None, transferOutput=True, \
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
        #do not pass the inputFileList to addGenericJob() because db arguments
        #  need to be added before them.
        job = self.addGenericJob(executable=self.convertImage,
            inputFile=inputFile,
            inputArgumentOption=inputArgumentOption,
            outputFile=outputFile,
            outputArgumentOption=outputArgumentOption,
            inputFileList=None, 
            frontArgumentList=frontArgumentList,
            extraArguments=extraArguments, 
            extraArgumentList=extraArgumentList,
            parentJobLs=parentJobLs,
            extraDependentInputLs=extraDependentInputLs,
            extraOutputLs=extraOutputLs,
            transferOutput=transferOutput,
            job_max_memory=job_max_memory,
            key2ObjectForJob=key2ObjectForJob,
            **keywords)
        return job
    
    def addCalculateDepthMeanMedianModeJob(self, executable=None,
        inputFile=None, outputFile=None, alignmentID=None,
        fractionToSample=0.001,
        whichColumn=None, maxNumberOfSamplings=1E7,
        inputStatName=None,
        extraArguments=None,
        parentJobLs=None, job_max_memory = 500,
        transferOutput=False, **keywords):
        """
        2012.6.15 turn maxNumberOfSamplings into integer when passing it to the job
        input is output of samtoolsDepth.
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
        job= self.addGenericJob(executable=executable,
            inputFile=inputFile,
            inputArgumentOption='-i',
            outputFile=outputFile,
            outputArgumentOption='-o',
            extraArgumentList=extraArgumentList,
            parentJobLs=parentJobLs,
            transferOutput=transferOutput,
            job_max_memory=job_max_memory, **keywords)
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
        job.addArguments(
            "--genome_drivername", objectWithDBArguments.drivername,
            "--genome_hostname", objectWithDBArguments.hostname, \
            "--genome_dbname", genome_dbname, \
            "--genome_db_user", objectWithDBArguments.db_user, \
            "--genome_db_passwd %s"%objectWithDBArguments.db_passwd)

        job.addArguments("--genome_schema genome")

        if getattr(objectWithDBArguments, 'port', None):
            job.addArguments("--genome_port=%s"%(objectWithDBArguments.port))
        return job


    def addGzipSubWorkflow(self, inputData=None, transferOutput=True,\
        outputDirPrefix="", topOutputDirJob=None, **keywords):
        """
        2012.8.2 bugfix.
        2012.7.19
        """
        sys.stderr.write("Adding gzip jobs for %s input job data ... "%(
            len(inputData.jobDataLs)))
        returnData = PassingData(topOutputDirJob=None)
        returnData.jobDataLs = []
        if inputData:
            if len(inputData.jobDataLs)>0:
                if topOutputDirJob is None:
                    topOutputDir = "%sGzip"%(outputDirPrefix)
                    topOutputDirJob = self.addMkDirJob(outputDir=topOutputDir)

                returnData.topOutputDirJob = topOutputDirJob
                for jobData in inputData.jobDataLs:
                    for inputF in set(jobData.fileLs):
                        #2013.08.16 do not work on same file
                        inputFBaseName = os.path.basename(inputF.name)
                        outputF = File(os.path.join(topOutputDirJob.output, \
                            '%s.gz'%(inputFBaseName)))
                        key2ObjectForJob = {}
                        extraArgumentList = []
                        #make sure set inputArgumentOption&outputArgumentOption to None, \
                        # otherwise addGenericJob will add "-i" and "-o" in front of it
                        job= self.addGenericJob(executable=self.gzip,
                            inputArgumentOption=None,
                            inputFile=inputF,
                            outputArgumentOption=None,
                            outputFile=outputF,
                            extraArgumentList=extraArgumentList,
                            parentJobLs=[topOutputDirJob]+jobData.jobLs,
                            extraDependentInputLs=None,
                            extraOutputLs=None,
                            transferOutput=transferOutput,
                            key2ObjectForJob=key2ObjectForJob,
                            job_max_memory=200, **keywords)
                        returnData.jobDataLs.append(PassingData(jobLs=[job],
                            vcfFile=None, file=outputF, fileLs=[outputF]))
        sys.stderr.write("no_of_jobs = %s.\n"%(self.no_of_jobs))
        return returnData

    def addAbstractMapperLikeJob(self, executable=None, \
        inputVCF=None, inputF=None, outputF=None,
        extraArguments=None, extraArgumentList=None,
        extraOutputLs=None,\
        parentJobLs=None, transferOutput=True, job_max_memory=200,\
        extraDependentInputLs=None, \
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
        job= self.addGenericJob(executable=executable, inputFile=inputF,
            outputFile=outputF, \
            extraArguments=extraArguments,\
            extraArgumentList=extraArgumentList,
            parentJobLs=parentJobLs,
            extraDependentInputLs=extraDependentInputLs,
            extraOutputLs=extraOutputLs,\
            transferOutput=transferOutput, \
            sshDBTunnel=sshDBTunnel, job_max_memory=job_max_memory, **keywords)
        return job

    def addSelectLineBlockFromFileJob(self,
        executable=None,
        inputFile=None, outputFile=None,
        startLineNumber=None, stopLineNumber=None,
        extraArguments=None,
        parentJobLs=None,
        extraDependentInputLs=None,
        transferOutput=False,
        job_max_memory=200, **keywords):
        """
        """
        extraArgumentList = ['-s %s'%(startLineNumber), '-t %s'%(stopLineNumber)]
        if extraArguments:
            extraArgumentList.append(extraArguments)

        job= self.addGenericJob(executable=executable, inputFile=inputFile,
            outputFile=outputFile, \
            extraArgumentList=extraArgumentList,
            parentJobLs=parentJobLs,
            extraDependentInputLs=extraDependentInputLs, \
            extraOutputLs=None,\
            transferOutput=transferOutput, \
            job_max_memory=job_max_memory, \
            **keywords)
        return job

    def addPlotLDJob(self, executable=None,
        inputFile=None, inputFileList=None, \
        outputFile=None, outputFnamePrefix=None,
        whichColumn=None, whichColumnHeader=None, whichColumnPlotLabel=None,
        title=None, logY=None, valueForNonPositiveYValue=-1,
        missingDataNotation='-nan',
        xColumnHeader=None,
        xColumnPlotLabel=None,
        chrLengthColumnHeader=None, chrColumnHeader=None,
        minChrLength=1000000,
        pos2ColumnHeader=None,
        minNoOfTotal=100, maxNoOfTotal=None,
        figureDPI=300, formatString='.', ylim_type=None,
        samplingRate=1, 
        need_svg=False, logCount=False,
        minDist=None, maxDist=None, movingAverageType=2,
        extraArguments=None,
        extraArgumentList=None,
        parentJobLs=None,
        extraDependentInputLs=None,
        transferOutput=True,
        job_max_memory=200, **keywords):
        """
('outputFname', 1, ): [None, 'o', 1, 'output file for the figure.'],\
('minNoOfTotal', 1, int): [100, 'i', 1, 
    'minimum no of total variants (denominator of inconsistent rate)'],\
('title', 1, ): [None, 't', 1, 'title for the figure.'],\
('figureDPI', 1, int): [200, 'f', 1, 'dpi for the output figures (png)'],\
('formatString', 1, ): ['.', 'a', 1, 'formatString passed to matplotlib plot'],\
('ylim_type', 1, int): [1, 'y', 1, 'y-axis limit type, 1: 0 to max. 2: min to max'],\
('samplingRate', 1, float): [0.001, 's', 1, 'how often you include the data'],\
('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.'],\
('whichColumn', 0, int): [3, 'w', 1, 
    'data from this column (index starting from 0) is plotted as y-axis value'],\
('whichColumnHeader', 0, ): ["", 'W', 1, 
    'column label (in the header) for the data to be plotted as y-axis value, substitute whichColumn'],\
('logWhichColumn', 0, int): [0, 'g', 0, 'whether to take -log of the whichColumn'],\
('whichColumnPlotLabel', 1, ): ['#SNPs in 100kb window', 'D', 1, \
    'plot label for data of the whichColumn', ],\
('chrLengthColumnHeader', 1, ): ['chrLength', 'c', 1, 'label of the chromosome length column', ],\
('chrColumnHeader', 1, ): ['CHR', 'C', 1, 'label of the chromosome column', ],\
('minChrLength', 1, int): [1000000, 'm', 1, \
    'minimum chromosome length for one chromosome to be included', ],\
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
            extraArgumentList.append("--chrLengthColumnHeader %s"%\
                (chrLengthColumnHeader))
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

        return self.addAbstractPlotJob(executable=executable,
            inputFileList=inputFileList,
            inputFile=inputFile,
            outputFile=outputFile, outputFnamePrefix=outputFnamePrefix,
            xColumnHeader=xColumnHeader, xColumnPlotLabel=xColumnPlotLabel,
            whichColumn=whichColumn,
            whichColumnHeader=whichColumnHeader,
            whichColumnPlotLabel=whichColumnPlotLabel,
            logY=logY, valueForNonPositiveYValue=valueForNonPositiveYValue,
            missingDataNotation=missingDataNotation,
            title=title,
            minNoOfTotal=minNoOfTotal,
            figureDPI=figureDPI, formatString=formatString,
            ylim_type=ylim_type, 
            samplingRate=samplingRate, need_svg=need_svg,
            extraArguments=extraArguments,
            extraArgumentList=extraArgumentList,
            parentJobLs=parentJobLs,
            extraDependentInputLs=extraDependentInputLs,
            transferOutput=transferOutput,
            job_max_memory=job_max_memory,
            **keywords)

    def addPlotVCFtoolsStatJob(self,
        executable=None,
        inputFileList=None, 
        outputFnamePrefix=None,
        xColumnPlotLabel=None, xColumnHeader=None,
        whichColumn=None, whichColumnHeader=None, whichColumnPlotLabel=None, 
        need_svg=False, logY=0, valueForNonPositiveYValue=-1,
        chrLengthColumnHeader=None, 
        chrColumnHeader=None,
        minChrLength=1000000, minNoOfTotal=100,
        figureDPI=300, ylim_type=None, samplingRate=1,
        logCount=False,
        tax_id=60711, sequence_type_id=1, chrOrder=None,
        extraArguments=None,
        parentJobLs=None,
        extraDependentInputLs=None,
        transferOutput=True, job_max_memory=200,
        sshDBTunnel=False, **keywords):
        """
        Examples

    outputFnamePrefix = os.path.join(plotOutputDir, 'noOfMendelErrors_along_chromosome')
    self.addPlotVCFtoolsStatJob(executable=self.PlotVCFtoolsStat,
        inputFileList=[splitPlinkLMendelFileSNPIDIntoChrPositionJob.output],
        outputFnamePrefix=outputFnamePrefix,
        whichColumn=None, whichColumnHeader="N", whichColumnPlotLabel="noOfMendelErrors",
        need_svg=False,
        logY=0, valueForNonPositiveYValue=-1,
        xColumnPlotLabel="genomePosition",
        chrLengthColumnHeader=None, chrColumnHeader="Chromosome", \
        minChrLength=100000, xColumnHeader="Start", minNoOfTotal=50,\
        figureDPI=100, ylim_type=2, samplingRate=1,\
        tax_id=self.ref_genome_tax_id,
        sequence_type_id=self.ref_genome_sequence_type_id, chrOrder=1,\
        parentJobLs=[splitPlinkLMendelFileSNPIDIntoChrPositionJob, plotOutputDirJob], \
        extraDependentInputLs=None, \
        extraArguments=None, transferOutput=True, sshDBTunnel=self.needSSHDBTunnel)


        remove argument positiveLog, rename logWhichColumn to logY
        # whichColumnPlotLabel and xColumnPlotLabel should not contain spaces or ( or ).
        #  because they will disrupt shell commandline


('whichColumn', 0, int): [3, 'w', 1,
    'data from this column (index starting from 0) is plotted as y-axis value'],\
('whichColumnHeader', 0, ): ["", 'W', 1, 
    'column label (in the header) for the data to be plotted as y-axis value, substitute whichColumn'],\
('logY', 0, int): [0, '', 1, 'value 0: nothing; 1: log(), 2: -log(). replacing self.logWhichColumn.'],\
('need_svg', 0, ): [0, 'n', 0, 'whether need svg output', ],\
('whichColumnPlotLabel', 1, ): ['#SNPs in 100kb window', 'D', 1, 
    'plot label for data of the whichColumn', ],\
('xColumnPlotLabel', 1, ): ['position', 'x', 1, 'x-axis label (posColumn) in manhattan plot', ],\
('chrLengthColumnHeader', 1, ): ['chrLength', 'c', 1, 'label of the chromosome length column', ],\
('chrColumnHeader', 1, ): ['CHR', 'C', 1, 'label of the chromosome column', ],\
('minChrLength', 1, int): [1000000, 'm', 1, 'minimum chromosome length for one chromosome to be included', ],\
('xColumnHeader', 1, ): ['BIN_START', 'l', 1,
    'label of the position column, BIN_START for binned vcftools output. POS for others.', ],\
('outputFnamePrefix', 0, ): [None, 'O', 1, 'output filename prefix (optional).'],\

('minNoOfTotal', 1, int): [100, 'i', 1, 
    'minimum no of total variants (denominator of inconsistent rate)'],\
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
        extraArgumentList = [
            '--xColumnHeader %s'%(xColumnHeader)
            ]
        extraOutputLs = [
            File('%s.png'%(outputFnamePrefix)),
            File('%s_hist.png'%(outputFnamePrefix))
            ]
        key2ObjectForJob = {}
        if ylim_type:
            extraArgumentList.append(f'--ylim_type {ylim_type}')
        if outputFnamePrefix:
            extraArgumentList.append(f'--outputFnamePrefix {outputFnamePrefix}')
        if minNoOfTotal:
            extraArgumentList.append(f'--minNoOfTotal {minNoOfTotal}')
        if figureDPI:
            extraArgumentList.append(f'--figureDPI {figureDPI}')
        if samplingRate:
            extraArgumentList.append(f'--samplingRate {samplingRate}')
        
        if need_svg:
            extraOutputLs.append(File('%s.svg'%(outputFnamePrefix)))
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
            extraArgumentList.append("--chrLengthColumnHeader %s"%\
                (chrLengthColumnHeader))
        if chrColumnHeader:
            extraArgumentList.append("--chrColumnHeader %s"%(chrColumnHeader))
        if logCount:
            extraArgumentList.append("--logCount")
        if valueForNonPositiveYValue:
            extraArgumentList.append("--valueForNonPositiveYValue %s"%\
                (valueForNonPositiveYValue))
        if sequence_type_id:
            extraArgumentList.append("--sequence_type_id %s"%(sequence_type_id))
        if tax_id:
            extraArgumentList.append("--tax_id %s"%(tax_id))
        if chrOrder is not None:
            extraArgumentList.append("--chrOrder %s"%(chrOrder))

        if extraArguments:
            extraArgumentList.append(extraArguments)
        job= self.addDBJob(
            executable=executable,
            inputFile=None,
            outputFile=None,
            inputFileList=inputFileList,
            extraArgumentList=extraArgumentList,
            parentJobLs=parentJobLs,
            extraDependentInputLs=extraDependentInputLs,
            extraOutputLs=extraOutputLs,
            transferOutput=transferOutput,
            key2ObjectForJob=key2ObjectForJob,
            job_max_memory=job_max_memory,
            sshDBTunnel=sshDBTunnel, objectWithDBArguments=self, **keywords)
        return job

    def addPlotGenomeWideDataJob(self, executable=None,
        inputFileList=None,
        inputFile=None,
        outputFnamePrefix=None, outputFile=None,
        xColumnPlotLabel=None, xColumnHeader=None,
        whichColumn=None, whichColumnHeader=None, whichColumnPlotLabel=None,
        logX=None, logY=None, valueForNonPositiveYValue=-1,
        xScaleLog=None, yScaleLog=None,
        missingDataNotation='NA',
        xtickInterval=None,
        drawCentromere=None, chrColumnHeader=None,
        minChrLength=100000, minNoOfTotal=100, maxNoOfTotal=None,
        figureDPI=300, formatString=".", ylim_type=1,
        samplingRate=1,
        logCount=False, need_svg=False,
        tax_id=60711, sequence_type_id=1, chrOrder=None,
        inputFileFormat=1, outputFileFormat=None,
        extraArguments=None,
        extraArgumentList=None,
        parentJobLs=None, extraDependentInputLs=None,
        transferOutput=True, job_max_memory=200, \
        objectWithDBGenomeArguments=None, sshDBTunnel=False, \
        **keywords):

        """
        Examples:
    outputFile = File(os.path.join(plotOutputDir, 'noOfMendelErrors_along_chromosome.png'))
    self.addPlotGenomeWideDataJob(inputFileList=None, \
        inputFile=splitPlinkLMendelFileSNPIDIntoChrPositionJob.output,\
        outputFile=outputFile,\
        whichColumn=None, whichColumnHeader="N",
        whichColumnPlotLabel="noOfMendelErrors", \
        logX=None, logY=None, valueForNonPositiveYValue=-1, \
        xScaleLog=None, yScaleLog=None,\
        missingDataNotation='NA',\
        xColumnPlotLabel="genomePosition", xColumnHeader="Start", \
        xtickInterval=20000000,\
        chrColumnHeader="Chromosome", \
        minChrLength=None, minNoOfTotal=None, maxNoOfTotal=None, \
        figureDPI=300, formatString=".", ylim_type=2, samplingRate=1, 
        logCount=False, need_svg=False,\
        tax_id=self.ref_genome_tax_id, 
        sequence_type_id=self.ref_genome_sequence_type_id, chrOrder=1,\
        inputFileFormat=1, outputFileFormat=None,\
        parentJobLs=[splitPlinkLMendelFileSNPIDIntoChrPositionJob, plotOutputDirJob], \
        extraDependentInputLs=None, \
        extraArguments=None, extraArgumentList=None, \
        transferOutput=True, job_max_memory=1000, sshDBTunnel=self.needSSHDBTunnel)

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
        job = self.addAbstractPlotJob(
            executable=executable,
            inputFileList=inputFileList,
            inputFile=inputFile,
            outputFile=outputFile,
            outputFnamePrefix=outputFnamePrefix,
            xColumnHeader=xColumnHeader, xColumnPlotLabel=xColumnPlotLabel,
            whichColumn=whichColumn,
            whichColumnHeader=whichColumnHeader,
            whichColumnPlotLabel=whichColumnPlotLabel,
            logX=logX, logY=logY,
            valueForNonPositiveYValue=valueForNonPositiveYValue,
            xScaleLog=xScaleLog, yScaleLog=yScaleLog,
            missingDataNotation=missingDataNotation,
            minNoOfTotal=minNoOfTotal, maxNoOfTotal=maxNoOfTotal,
            figureDPI=figureDPI, formatString=formatString,
            ylim_type=ylim_type,
            samplingRate=samplingRate, need_svg=need_svg,
            inputFileFormat=inputFileFormat,
            outputFileFormat=outputFileFormat,
            extraArguments=extraArguments,
            extraArgumentList=extraArgumentList,
            parentJobLs=parentJobLs,
            extraDependentInputLs=extraDependentInputLs,
            transferOutput=transferOutput,
            job_max_memory=job_max_memory,
            sshDBTunnel=sshDBTunnel,
            objectWithDBGenomeArguments=objectWithDBGenomeArguments,
            **keywords)
        return job

    def addAbstractPlotJob(self, executable=None,
        inputFileList=None, inputFile=None,
        outputFile=None,
        outputFnamePrefix=None,
        xColumnHeader=None, xColumnPlotLabel=None,
        whichColumn=None, whichColumnHeader=None,
        whichColumnPlotLabel=None,
        logX=None, logY=None, valueForNonPositiveYValue=-1,
        xScaleLog=0, yScaleLog=0,
        missingDataNotation='NA',
        title=None,
        minNoOfTotal=None, maxNoOfTotal=None,
        figureDPI=300, formatString='.', markerSize=None,
        ylim_type=1, samplingRate=1, legendType=None,
        need_svg=False,
        inputFileFormat=None, outputFileFormat=None,
        extraArguments=None,
        extraArgumentList=None,
        parentJob=None, parentJobLs=None,
        extraDependentInputLs=None, extraOutputLs=None,
        transferOutput=True,
        job_max_memory=200,
        sshDBTunnel=False, key2ObjectForJob=None,
        objectWithDBArguments=None, **keywords):
        """
        #no spaces or parenthesis or any other shell-vulnerable letters
        #  in the x or y axis labels (whichColumnPlotLabel, xColumnPlotLabel)
        Check AbstractMatrixFileWalker.py or AbstractPlot.py for updated arguments:

('outputFname', 0, ): [None, 'o', 1, 'output file for the figure.'],\
('minNoOfTotal', 1, int): [100, 'M', 1, 
    'minimum no of total variants (denominator of inconsistent rate)'],\
('title', 0, ): [None, 't', 1, 'title for the figure.'],\
('figureDPI', 1, int): [200, 'f', 1, 'dpi for the output figures (png)'],\
('formatString', 1, ): ['-', '', 1, 'formatString passed to matplotlib plot'],\
('ylim_type', 1, int): [1, 'y', 1, 'y-axis limit type, 1: whatever matplotlib decides. 2: min to max'],\
('samplingRate', 1, float): [1, 's', 1, 'how often you include the data, a probability between 0 and 1.'],\
('whichColumn', 0, int): [3, 'w', 1, 
    'data from this column (index starting from 0) is plotted as y-axis value'],\
('whichColumnHeader', 0, ): ["", 'W', 1, 
    'column header for the data to be plotted as y-axis value, substitute whichColumn'],\
('whichColumnPlotLabel', 0, ): ['', 'D', 1, 'plot label for data of the whichColumn', ],\
('logY', 0, int): [0, '', 1, 'value 0: nothing; 1: log(), 2: -log(). replacing self.logWhichColumn.'],\
('valueForNonPositiveYValue', 1, float): [50, '', 1, 
    'if the whichColumn value is not postive and logWhichColumn is on,\
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
            extraArgumentList.append("--valueForNonPositiveYValue %s"%\
                (valueForNonPositiveYValue))
        if inputFileFormat:
            extraArgumentList.append("--inputFileFormat %s"%(inputFileFormat))
        if outputFileFormat:
            extraArgumentList.append("--outputFileFormat %s"%(outputFileFormat))
        if need_svg:
            extraArgumentList.append('--need_svg')
            if not outputFnamePrefix:
                outputFnamePrefix = os.path.splitext(outputFile.name)[0]
            extraOutputLs.append(File('%s.svg'%(outputFnamePrefix)))
        if extraArguments:
            extraArgumentList.append(extraArguments)

        job = self.addGenericJob(executable=executable,
            inputArgumentOption='-i',
            inputFile=inputFile, 
            outputArgumentOption='-o',
            outputFile=outputFile,
            inputFileList = inputFileList,
            extraArgumentList=extraArgumentList,
            parentJob=parentJob, parentJobLs=parentJobLs,
            extraDependentInputLs=extraDependentInputLs,
            extraOutputLs=extraOutputLs,
            transferOutput=transferOutput,
            key2ObjectForJob=key2ObjectForJob,
            job_max_memory=job_max_memory,
            sshDBTunnel=sshDBTunnel,
            objectWithDBArguments=objectWithDBArguments,
            **keywords)
        return job

    def addAbstractMatrixFileWalkerJob(self, executable=None,
        inputFileList=None,
        inputFile=None,
        outputFile=None,
        outputFnamePrefix=None,
        whichColumn=None, whichColumnHeader=None,
        logY=None, valueForNonPositiveYValue=-1,
        minNoOfTotal=10,
        samplingRate=1,
        inputFileFormat=None, outputFileFormat=None,
        extraArguments=None,
        extraArgumentList=None,
        parentJob=None, parentJobLs=None, extraOutputLs=None,
        extraDependentInputLs=None,
        transferOutput=True, job_max_memory=200,
        sshDBTunnel=False,
        objectWithDBArguments=None, **keywords):
        """
('outputFname', 0, ): [None, 'o', 1, 'output file for the figure.'],\
('minNoOfTotal', 1, int): [100, 'i', 1, 
    'minimum no of total variants (denominator of inconsistent rate)'],\
('samplingRate', 1, float): [1, 's', 1, 
    'how often you include the data, a probability between 0 and 1.'],\
('whichColumn', 0, int): [3, 'w', 1, 
    'data from this column (index starting from 0) is plotted as y-axis value'],\
('whichColumnHeader', 0, ): ["", 'W', 1, \
    'column header for the data to be plotted as y-axis value, substitute whichColumn'],\
('whichColumnPlotLabel', 0, ): ['', 'D', 1, 'plot label for data of the whichColumn', ],\
('logWhichColumn', 0, int): [0, 'g', 0, 'whether to take -log of the whichColumn'],\
('positiveLog', 0, int): [0, 'p', 0, 'toggle to take log, rather than -log(), \
    only effective when logWhichColumn is toggled. '],\
('valueForNonPositiveYValue', 1, float): [50, '', 1, 
    'if the whichColumn value is not postive and logWhichColumn is on,\
        what yValue should be.'],\

        """
        return self.addAbstractPlotJob(executable=executable, 
            inputFileList=inputFileList,
            inputFile=inputFile,
            outputFile=outputFile, 
            outputFnamePrefix=outputFnamePrefix,
            xColumnHeader=None, xColumnPlotLabel=None,
            whichColumn=whichColumn,
            whichColumnHeader=whichColumnHeader, whichColumnPlotLabel=None,
            logY=logY, \
            valueForNonPositiveYValue=valueForNonPositiveYValue,
            missingDataNotation=None,
            minNoOfTotal=minNoOfTotal,
            figureDPI=None, formatString=None, ylim_type=None,
            samplingRate=samplingRate, need_svg=False,
            extraArguments=extraArguments,
            extraArgumentList=extraArgumentList,
            parentJob=parentJob, parentJobLs=parentJobLs,
            extraOutputLs=extraOutputLs,
            extraDependentInputLs=extraDependentInputLs,
            transferOutput=transferOutput, 
            job_max_memory=job_max_memory,
            sshDBTunnel=sshDBTunnel,
            objectWithDBArguments=objectWithDBArguments, **keywords)

    def addAbstractGenomeFileWalkerJob(self, executable=None,
        inputFileList=None,
        inputFile=None,
        outputFile=None,
        outputFnamePrefix=None,
        whichColumn=None, whichColumnHeader=None,
        logY=None, valueForNonPositiveYValue=-1,
        minNoOfTotal=10,
        samplingRate=1,
        chrColumnHeader=None,
        tax_id=60711, sequence_type_id=1, chrOrder=None,
        positionHeader=None,
        inputFileFormat=None, outputFileFormat=None,
        extraArguments=None,
        extraArgumentList=None,
        parentJob=None, parentJobLs=None,
        extraDependentInputLs=None,
        transferOutput=True, job_max_memory=200,
        sshDBTunnel=False,
        objectWithDBGenomeArguments=None, **keywords):
        """
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
        return self.addAbstractPlotJob(executable=executable,
            inputFileList=inputFileList,
            inputFile=inputFile,
            outputFile=outputFile, 
            outputFnamePrefix=outputFnamePrefix,
            xColumnHeader=None, xColumnPlotLabel=None,
            whichColumn=whichColumn,
            whichColumnHeader=whichColumnHeader, whichColumnPlotLabel=None,
            logY=logY,
            valueForNonPositiveYValue=valueForNonPositiveYValue,
            missingDataNotation=None,
            minNoOfTotal=minNoOfTotal,
            figureDPI=None, formatString=None, ylim_type=None, 
            samplingRate=samplingRate, need_svg=False,
            extraArguments=extraArguments,
            extraArgumentList=extraArgumentList,
            parentJob=parentJob, parentJobLs=parentJobLs, 
            extraDependentInputLs=extraDependentInputLs,
            transferOutput=transferOutput,
            job_max_memory=job_max_memory,
            sshDBTunnel=sshDBTunnel, 
            objectWithDBGenomeArguments=objectWithDBGenomeArguments,
            **keywords)


    def addDrawHistogramJob(self, executable=None,
        inputFileList=None, inputFile=None,
        outputFile=None, outputFnamePrefix=None,
        whichColumn=None, whichColumnHeader=None, 
        whichColumnPlotLabel=None,
        xScaleLog=0, yScaleLog=0,
        logY=None, valueForNonPositiveYValue=-1, missingDataNotation='NA',
        title=None,
        minNoOfTotal=None, figureDPI=300, formatString=None, ylim_type=None,
        samplingRate=1,
        need_svg=False, legendType=None,
        logCount=False, inputFileFormat=None,
        extraArguments=None,
        extraArgumentList=None,
        parentJobLs=None,
        extraDependentInputLs=None,
        transferOutput=True, job_max_memory=200,
        **keywords):
        """
        #no spaces or parenthesis or any other shell-vulnerable letters
        #  in the x or y axis labels (whichColumnPlotLabel, xColumnPlotLabel)
('outputFname', 0, ): [None, 'o', 1, 'output file for the figure.'],\
('minNoOfTotal', 1, int): [100, 'i', 1, 
    'minimum no of total variants (denominator of inconsistent rate)'],\
('title', 0, ): [None, 't', 1, 'title for the figure.'],\
('figureDPI', 1, int): [200, 'f', 1, 'dpi for the output figures (png)'],\
('formatString', 1, ): ['-', '', 1, 'formatString passed to matplotlib plot'],\
('ylim_type', 1, int): [1, 'y', 1, 'y-axis limit type, 1: 0 to max. 2: min to max'],\
('samplingRate', 1, float): [1, 's', 1, 'how often you include the data, a probability between 0 and 1.'],\
('whichColumn', 0, int): [3, 'w', 1, 
    'data from this column (index starting from 0) is plotted as y-axis value'],\
('whichColumnHeader', 0, ): ["", 'W', 1, 
    'column header for the data to be plotted as y-axis value, substitute whichColumn'],\
('whichColumnPlotLabel', 0, ): ['', 'D', 1, 'plot label for data of the whichColumn', ],\
('logWhichColumn', 0, int): [0, 'g', 0, 'whether to take -log of the whichColumn'],\
('positiveLog', 0, int): [0, 'p', 0, 'toggle to take log, rather than -log(), \
    only effective when logWhichColumn is toggled. '],\
('valueForNonPositiveYValue', 1, float): [50, '', 1, 
    'if the whichColumn value is not postive and logWhichColumn is on,\
        what yValue should be.'],\
('need_svg', 0, ): [0, 'n', 0, 'whether need svg output', ],\

        """
        if extraArguments is None:
            extraArguments = ""
        if logCount:
            extraArguments += " --logCount "
        return self.addAbstractPlotJob(executable=executable,
            inputFileList=inputFileList,
            inputFile=inputFile,
            outputFile=outputFile, 
            outputFnamePrefix=outputFnamePrefix,
            xColumnHeader=None, xColumnPlotLabel=None,
            title=title,
            whichColumn=whichColumn,
            whichColumnHeader=whichColumnHeader,
            whichColumnPlotLabel=whichColumnPlotLabel,
            xScaleLog=xScaleLog, yScaleLog=yScaleLog,
            logY=logY,
            valueForNonPositiveYValue=valueForNonPositiveYValue,
            missingDataNotation=missingDataNotation,
            minNoOfTotal=minNoOfTotal,
            figureDPI=figureDPI, formatString=formatString,
            ylim_type=ylim_type,
            samplingRate=samplingRate, need_svg=need_svg,
            legendType=legendType,
            inputFileFormat=inputFileFormat,
            extraArguments=extraArguments,
            extraArgumentList=extraArgumentList,
            parentJobLs=parentJobLs,
            extraDependentInputLs=extraDependentInputLs,
            transferOutput=transferOutput, 
            job_max_memory=job_max_memory,
            **keywords)

    def addDraw2DHistogramOfMatrixJob(self, executable=None,
        inputFileList=None,
        inputFile=None,
        outputFile=None,
        outputFnamePrefix=None,
        xColumnHeader=None, xColumnPlotLabel=None,
        whichColumn=None,
        whichColumnHeader=None, whichColumnPlotLabel=None,
        logX=False, logY=False, logZ=False,
        valueForNonPositiveYValue=-1,
        missingDataNotation='NA',
        minNoOfTotal=100,
        figureDPI=300, formatString='.',
        samplingRate=1,
        need_svg=False,
        inputFileFormat=None, outputFileFormat=None,
        zColumnHeader=None,
        extraArguments=None,
        extraArgumentList=None,
        parentJobLs=None,
        extraDependentInputLs=None,
        transferOutput=True,
        job_max_memory=200, **keywords):
        """
        executable could be None, default is self.Draw2DHistogramOfMatrix

        """
        if extraArgumentList is None:
            extraArgumentList = []
        if zColumnHeader:
            extraArgumentList.append("--zColumnHeader %s"%(zColumnHeader))
        if logZ:
            extraArgumentList.append("--logZ %s"%(logZ))
        if executable is None:
            executable = self.Draw2DHistogramOfMatrix
        return self.addAbstractPlotJob(
            executable=executable, inputFileList=inputFileList,
            inputFile=inputFile,
            outputFile=outputFile,
            outputFnamePrefix=outputFnamePrefix,
            xColumnHeader=xColumnHeader, xColumnPlotLabel=xColumnPlotLabel,
            whichColumn=whichColumn,
            whichColumnHeader=whichColumnHeader,
            whichColumnPlotLabel=whichColumnPlotLabel,
            logX=logX, logY=logY,
            valueForNonPositiveYValue=valueForNonPositiveYValue,
            missingDataNotation=missingDataNotation,
            minNoOfTotal=minNoOfTotal,
            figureDPI=figureDPI, formatString=formatString,
            ylim_type=None,
            samplingRate=samplingRate,
            need_svg=need_svg,
            inputFileFormat=inputFileFormat,
            outputFileFormat=outputFileFormat,
            extraArguments=extraArguments,
            extraArgumentList=extraArgumentList,
            parentJobLs=parentJobLs,
            extraDependentInputLs=extraDependentInputLs,
            transferOutput=transferOutput,
            job_max_memory=job_max_memory,
            **keywords)

    def setupMoreOutputAccordingToSuffixAndNameTupleList(self,
        outputFnamePrefix=None,
        suffixAndNameTupleList=None, extraOutputLs=None, key2ObjectForJob=None):
        """
        split from addPlinkJob()
        """
        for suffixNameTuple in suffixAndNameTupleList:
            if len(suffixNameTuple)==1:
                suffix = suffixNameTuple[0]
                name = suffix[1:].replace('.', '_')
                #replace dot with underscore. as dot is used to access
                #  method/attribute of python object.
                # i.e. ".prune.in" is accessible as job.prune_inFile
            elif len(suffixNameTuple)>=2:
                suffix, name = suffixNameTuple[:2]
            outputFile = File('%s%s'%(outputFnamePrefix, suffix))
            extraOutputLs.append(outputFile)
            key2ObjectForJob['%sFile'%(name)] = outputFile

if __name__ == '__main__':
    from argparse import ArgumentParser
    ap = ArgumentParser()
    ap.add_argument("-i", "--input_path", type=str, required=True,
        help="the path to the input folder or file.")
    ap.add_argument("--inputSuffixList", type=str,
        help='Coma-separated list of input file suffices. Used to exclude input files.'
        'If None, no exclusion. The dot is part of the suffix, .tsv not tsv.'
        'Common zip suffices (.gz, .bz2, .zip, .bz) will be ignored in obtaining the suffix.')
    ap.add_argument("-F", "--pegasusFolderName", default='input',
        help='The path relative to the workflow running root. '
        'This folder will contain pegasus input & output. '
        'It will be created during the pegasus staging process. '
        'It is useful to separate multiple sub-workflows. '
        'If empty or None, everything is in the pegasus root.')
    
    ap.add_argument("-o", "--output_path", type=str, required=True,
        help="The path to the output file that will contain the Pegasus DAG.")
    
    ap.add_argument("--tmpDir", type=str, default='/tmp/',
        help='Default: %(default)s. A local folder for some jobs (MarkDup) to store temp data.'
        '/tmp/ can be too small sometimes.')
    ap.add_argument("--max_walltime", type=int, default=4320,
        help='Default: %(default)s. Maximum wall time for any job, in minutes. 4320=3 days.'
        'Used in addGenericJob(). Most clusters have upper limit for runtime.')
    ap.add_argument("--jvmVirtualByPhysicalMemoryRatio", type=float, default=1.2,
        help='Default: %(default)s. '
        'If a job virtual memory (usually 1.2X of JVM resident memory) exceeds request, '
        "it will be killed on some clusters. "
        "This will make sure your job requests enough memory.")
    
    ap.add_argument("-l", "--site_handler", required=True,
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
    instance = AbstractWorkflow(
        input_path=args.input_path,
        inputSuffixList=args.inputSuffixList, 
        pegasusFolderName=args.pegasusFolderName,
        output_path=args.output_path,
        
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
