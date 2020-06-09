#!/usr/bin/env python3
"""
A generic workflow that map-reduces inputs of one or multiple genomic files (i.e. multi-chromosome, tabix-indexed )
    parent class is AbstractNGSWorkflow.

Examples:
    #2013.11.24
    %s -I ./FindNewRefCoordinates_Method109_vs_3488_BWA_F49.2013.Jul.19T141746/folderReduceGzip/
        -H -C 30 -l hcondor -D ~/NetworkData/vervet/db/
        -o dags/Method109_To_225_LiftOverProbability.xml --inputSuffixList .tsv
        --db_user yh -z localhost
        --ref_ind_seq_id 3280 --ref_genome_tax_id 60711
        --ref_genome_sequence_type_id 9 --ref_genome_version 2
"""
import sys, os, math
__doc__ = __doc__%(sys.argv[0])

import pegaflow
from pegaflow.DAX3 import Executable, File, PFN
from palos import ProcessOptions, PassingData, utils
from palos.io.FastaFile import FastaFile
from palos.pegasus.MapReduceGenomeFileWorkflow import MapReduceGenomeFileWorkflow

ParentClass = MapReduceGenomeFileWorkflow

class TestMapReduceGenomeFileWorkflow(ParentClass):
    __doc__ = __doc__
    option_default_dict = ParentClass.option_default_dict.copy()
    
    # default overlap is 5Mb, overlap is 500K
    option_default_dict[('intervalOverlapSize', 1, int)][0] = 1000000
    option_default_dict[('intervalSize', 1, int)][0] = 5000000
    option_default_dict[('max_walltime', 1, int)][0] = 1300	#under 23 hours
    option_default_dict[('inputSuffixList', 0, )][0] = None	#do not set input suffix
    
    def __init__(self,  **keywords):
        """
        """
        self.needSplitChrIntervalData = True
        ParentClass.__init__(self, **keywords)
        self.needSplitChrIntervalData = True
    
        # output is already gzipped
        self.needGzipPreReduceReturnData = True
        self.needGzipReduceReturnData = False
        
    def connectDB(self):
        """
        """
        ParentClass.connectDB(self)
        
        from palos.db import SunsetDB
        
        db_main = SunsetDB.SunsetDB(drivername=self.drivername, db_user=self.db_user,
            db_passwd=self.db_passwd, hostname=self.hostname, dbname=self.dbname,
            schema=self.schema)
        db_main.setup(create_tables=False)
        self.db_main = db_main
        
        if not self.data_dir:
            self.data_dir = db_main.data_dir
        
        if not self.local_data_dir:
            self.local_data_dir = db_main.data_dir
    
    def registerExecutables(self):
        """
        """
        ParentClass.registerExecutables(self)
        
        self.registerOneExecutable(path=os.path.join(self.pymodulePath,
            "polymorphism/mapper/ComputeLiftOverLocusProbability.py"),\
            name='ComputeLiftOverLocusProbability', \
            clusterSizeMultiplier=1)
    
    def preReduce(self, outputDirPrefix="", passingData=None, transferOutput=True, **keywords):
        """
        2012.9.17
        """
        returnData = ParentClass.preReduce(self, outputDirPrefix=outputDirPrefix,\
            passingData=passingData, transferOutput=transferOutput, **keywords)
        #add a stat merge job and a genome wide plot job
        outputFile = File(os.path.join(self.reduceOutputDirJob.output,
            'locusLiftOverProbability.tsv'))
        self.reduceJob = self.addStatMergeJob(
            statMergeProgram=self.mergeSameHeaderTablesIntoOne, \
            outputF=outputFile, \
            parentJobLs=[self.reduceOutputDirJob],
            extraDependentInputLs=None, transferOutput=False)
        
        sortProbabilityFile = File(os.path.join(self.reduceOutputDirJob.output,
            'locusLiftOverProbability.sorted.tsv'))
        sortProbabilityJob = self.addSortJob(inputFile=self.reduceJob.output,
            outputFile=sortProbabilityFile, \
            parentJobLs=[self.reduceJob], \
            extraOutputLs=None, transferOutput=False, \
            extraArgumentList=["""-k1,1 -k2,3n """], \
            sshDBTunnel=None,\
            job_max_memory=4000, walltime=120)
        #2013.12.3 Tab delimiter syntax (-t$'\t') is removed because it can't be passed correctly.
        #2013.12.3 Tried -t "`/bin/echo -e '\t'`" as well, didn't work either.
        # However since each column field doesn't contain blank,
        #   it is fine to just use the default separator (non-blank to blank).
        
        returnData.jobDataLs.append(self.constructJobDataFromJob(sortProbabilityJob))
        
        outputFile = File(os.path.join(self.plotDirJob.output, 'locusLiftOverProbability.png'))
        self.addPlotGenomeWideDataJob(inputFileList=None,
            inputFile=self.reduceJob.output,
            outputFile=outputFile,\
            whichColumn=None,
            whichColumnHeader="mapPvalue", whichColumnPlotLabel="mapPvalue", \
            logX=None, logY=2, valueForNonPositiveYValue=-1, \
            xScaleLog=None, yScaleLog=None,\
            missingDataNotation='NA',\
            xColumnPlotLabel="genomePosition", xColumnHeader="oldStart", \
            xtickInterval=0,\
            drawCentromere=True, chrColumnHeader="oldChromosome", \
            minChrLength=None, minNoOfTotal=None, maxNoOfTotal=None, \
            figureDPI=100, formatString=".", ylim_type=2,
            samplingRate=1, logCount=False, need_svg=True,\
            tax_id=self.ref_genome_tax_id,
            sequence_type_id=self.ref_genome_sequence_type_id, chrOrder=1,\
            inputFileFormat=1, outputFileFormat=None,\
            parentJobLs=[self.reduceJob], \
            extraDependentInputLs=None, \
            extraArguments=None, extraArgumentList=None, \
            transferOutput=True, job_max_memory=1000, sshDBTunnel=self.needSSHDBTunnel)
        #xtickInterval=0 means no ticks on x-axis.
        
        outputFile = File( os.path.join(self.plotDirJob.output, 'locusLiftOverProbabilityHist.png'))
        #no spaces or parenthesis or any other shell-vulnerable letters in the x
        #   or y axis labels (whichColumnPlotLabel, xColumnPlotLabel)
        self.addDrawHistogramJob(executable=self.DrawHistogram,
            inputFileList=[self.reduceJob.output], \
            outputFile=outputFile, \
            whichColumnHeader="mapPvalue", whichColumnPlotLabel="minusLogLiftOverPvalue", \
            xScaleLog=0, yScaleLog=1, \
            logCount=False, logY=2, valueForNonPositiveYValue=50,\
            minNoOfTotal=10,\
            figureDPI=100, samplingRate=1,legendType=1, \
            parentJobLs=[self.plotDirJob, self.reduceJob], \
            extraDependentInputLs=None, \
            transferOutput=True, job_max_memory=8000)
        
        return returnData
    
    def mapEachInterval(self, inputJobData=None, selectIntervalJobData=None, \
        chromosome=None,intervalData=None,\
        mapEachChromosomeData=None, \
        passingData=None, transferOutput=False, **keywords):
        """
        #. extract flanking sequences from the input Input (ref sequence file => contig ref sequence)
        #. blast them
        #. run FindSNPPositionOnNewRefFromFlankingBlastOutput.py
            #. where hit length match query length, and
            #    no of mismatches <=2 => good => infer new coordinates
        #. output a mapping file between old SNP and new SNP coordinates.
            #. reduce this thing by combining everything
        #. make a new Input file based on the input split Input file
            (replace contig ID , position with the new one's,
                remove the header part regarding chromosomes or replace it)

        """
        returnData = PassingData(no_of_jobs = 0)
        returnData.jobDataLs = []
        passingData.intervalFileBasenamePrefix
        passingData.splitInputFile
        """
        ## 2013.06.19 structures available from passingData, specific to the interval
        passingData.splitInputFile = splitInputFile
        passingData.unitNumber = unitNumber
        passingData.intervalFileBasenamePrefix = '%s_%s_splitInput_u%s'%(
            chromosome, commonPrefix, unitNumber)
        passingData.noOfIndividuals = jobData.file.noOfIndividuals
        passingData.span = self.intervalSize + self.intervalOverlapSize*2
        """
        #add one computing job
        outputFile = File(os.path.join(self.mapDirJob.output,
            "%s.%s.probability.tsv.gz"%(passingData.fileBasenamePrefix,\
            intervalData.interval)))
        locusIntervalDeltaOutputFile = File(os.path.join(self.mapDirJob.output,
            "%s.%s.locusIntervalDelta.tsv.gz"%(passingData.fileBasenamePrefix,
            intervalData.interval)))
        job = self.addAbstractMatrixFileWalkerJob(
            executable=self.ComputeLiftOverLocusProbability, \
            inputFile=selectIntervalJobData.file, outputFile=outputFile, \
            whichColumn=None, whichColumnHeader=None, \
            logY=None, valueForNonPositiveYValue=-1, \
            minNoOfTotal=1, samplingRate=1, \
            inputFileFormat=None, outputFileFormat=None,\
            extraArgumentList=["--locusIntervalDeltaOutputFname", locusIntervalDeltaOutputFile, \
                "--startPosition %s"%(intervalData.start), "--stopPosition %s"%(intervalData.stop)],
            parentJobLs=[selectIntervalJobData.job],
            extraOutputLs=[locusIntervalDeltaOutputFile],\
            transferOutput=transferOutput, job_max_memory=2000, sshDBTunnel=False)
            #For each interval, probabilities are not calculated for loci in
            #  extra segment (from overlapStart to start).
        returnData.jobDataLs.append(self.constructJobDataFromJob(job))
        return returnData

    def reduceEachChromosome(self, chromosome=None, passingData=None,
        mapEachInputDataLs=None, 
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
        #reduce matrix by chosen column and average p-value
        
        outputFile = File(os.path.join(self.reduceEachChromosomeDirJob.output,
            'chr_%s_LocusLiftOverProbability.tsv.gz'%(chromosome)))
        reduceChromosomeJob = self.addStatMergeJob(
            statMergeProgram=self.mergeSameHeaderTablesIntoOne, \
            outputF=outputFile, \
            parentJobLs=[self.reduceEachChromosomeDirJob],extraOutputLs=None, \
            extraDependentInputLs=None, transferOutput=False)
            #extraArgumentList=['--keyColumnLs 0-6 --valueColumnLs 7'],\
        mapEachIntervalDataLs = chromosome2mapEachIntervalDataLs.get(chromosome)
        for mapEachIntervalData in mapEachIntervalDataLs:
            for jobData in mapEachIntervalData.jobDataLs:
                self.addInputToMergeJob(reduceChromosomeJob, parentJobLs=[jobData.job])
            
        #add the reduction job to final stat merge job
        self.addInputToMergeJob(self.reduceJob, parentJobLs=[reduceChromosomeJob])
        
        return returnData
    

if __name__ == '__main__':
    main_class = TestMapReduceGenomeFileWorkflow
    po = ProcessOptions(sys.argv, main_class.option_default_dict,
        error_doc=main_class.__doc__)
    instance = main_class(**po.long_option2value)
    instance.run()	
