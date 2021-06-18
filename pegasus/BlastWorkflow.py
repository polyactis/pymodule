#!/usr/bin/env python3
"""
Examples:
    %s --databaseFname /Network/Data/vervet/db/individual_sequence/524_superContigsMinSize2000.fasta 
        --inputFname ~/script/vervet/data/OphoffMethylation/DMR330K_ProbeSeq.fasta --maxNoOfMismatches 2
        -l condorpool -j condorpool --blockSize 500
        -C 1 -o workflow/BlastDMR330K_ProbeSeqAgainst524.xml
    
    %s --databaseFname ~/NetworkData/vervet/db/individual_sequence/3280_vervet_ref_6.0.3.fasta
        -i ~/script/vervet/data/194SNPData/allSNPFlanks_refAllele.fa
        --maxNoOfMismatches 5 -l hcondor -j hcondor -C 1
        -o dags/Blast/Blast194SNPFlankAgainst3280_5Mismatches.xml
        --blockSize 10 --formatdbPath ~/bin/blast/bin/formatdb
        --blastallPath ~/bin/blast/bin/blastall

"""
import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0])

from pegaflow.api import File
from palos import ProcessOptions, PassingData, utils
from palos.pegasus.AbstractBioinfoWorkflow import AbstractBioinfoWorkflow

ParentClass = AbstractBioinfoWorkflow
class BlastWorkflow(ParentClass):
    __doc__ = __doc__
    option_default_dict = ParentClass.option_default_dict.copy()
    option_default_dict.update({
        ("inputFname", 1, ): ["", 'i', 1, 'the input fasta file'],\
        ("formatdbPath", 1, ): ["/usr/bin/formatdb", 'f', 1,
            'path to formatdb, index fasta database file'],\
        ("blastallPath", 1, ): ["/usr/bin/blastall", 's', 1,
            'path to blastall'],\
        ("blockSize", 1, int): [1000, 'c', 1,
            'how many sequences each blast job handles'],\
        ('databaseFname', 1, ): ['', 'd', 1,
            'filename of the database to blast against, must be indexed', ],\
        ('minNoOfIdentities', 0, int): [None, 'm', 1,
            'minimum number of identities between a query and target', ],\
        ('maxNoOfMismatches', 0, int): [None, 'a', 1,
            'minimum number of mismatches between a query and target', ],\
        ('minIdentityPercentage', 0, float): [None, 'n', 1,
            'minimum percentage of identities between a query and target', ],\
        
        })

    def __init__(self,  **keywords):
        """
        2012.5.23
        """
        ParentClass.__init__(self, **keywords)
    
    def getNoOfSequencesFromFasta(self, inputFastaFname=None):
        """
        2012.5.24
        """
        sys.stderr.write("Getting number of sequences from %s ..."%(inputFastaFname))
        inf = utils.openGzipFile(inputFastaFname)
        no_of_sequences = 0
        for line in inf:
            if line[0]=='>':
                no_of_sequences += 1
        del inf
        sys.stderr.write("%s sequences.\n"%(no_of_sequences))
        return no_of_sequences
    
    def addSplitFastaFileJob(self, executable=None, inputFile=None, outputFnamePrefix=None, \
        noOfSequencesPerSplitFile=1000, filenameSuffix="", noOfTotalSequences=1000000,\
        parentJobLs=[], extraDependentInputLs=[], transferOutput=False, \
        extraArguments=None, job_max_memory=500, **keywords):
        """
        2012.5.24
        """
        noOfSplitFiles = int(math.ceil(noOfTotalSequences/float(noOfSequencesPerSplitFile)))
        suffixLength = len(repr(noOfSplitFiles))
        
        job = self.addGenericJob(executable=executable, inputArgumentOption="-i",
            inputFile=inputFile, extraArgumentList=["--noOfSequences %s"%(noOfSequencesPerSplitFile), \
                "--outputFnamePrefix", outputFnamePrefix,
                '--filenameSuffix %s'%(filenameSuffix),
                '--suffixLength %s'%(suffixLength)],
            parentJobLs=parentJobLs, extraArguments=extraArguments, \
            extraDependentInputLs=extraDependentInputLs,
            job_max_memory=job_max_memory)
        
        for i in range(noOfSplitFiles):	#start from 0
            splitFname = utils.comeUpSplitFilename(outputFnamePrefix=outputFnamePrefix,
                suffixLength=suffixLength, fileOrder=i,\
                filenameSuffix=filenameSuffix)
            splitFile = File(splitFname)
            self.addJobUse(job, file=splitFile, is_input=False, transfer=transferOutput)
        return job
    
    def addBlastWrapperJob(self, executable=None, inputFile=None,
        outputFile=None, outputFnamePrefix=None, databaseFile=None,\
        maxNoOfMismatches=None, minNoOfIdentities=None,
        minIdentityPercentage=None, blastallPath=None, \
        parentJobLs=[], extraDependentInputLs=[], transferOutput=False, \
        extraArguments=None, job_max_memory=2000, **keywords):
        """
        2012.5.24
        """
        extraArgumentList = ['-l %s'%blastallPath, '--databaseFname', databaseFile]
        if outputFnamePrefix:
            extraArgumentList.append("--outputFnamePrefix %s "%(outputFnamePrefix))
        if maxNoOfMismatches:
            extraArgumentList.append('--maxNoOfMismatches %s'%maxNoOfMismatches)
        if minNoOfIdentities:
            extraArgumentList.append("--minNoOfIdentities %s"%(minNoOfIdentities))
        if minIdentityPercentage:
            extraArgumentList.append("--minIdentityPercentage %s"%(minIdentityPercentage))
        if extraArguments:
            extraArgumentList.append(extraArguments)
        return self.addGenericJob(executable=executable, inputFile=inputFile,
            outputFile=outputFile, \
            parentJobLs=parentJobLs, extraDependentInputLs=extraDependentInputLs,
            transferOutput=transferOutput, \
            extraArgumentList=extraArgumentList, job_max_memory=job_max_memory)
        
    def addJobs(self, inputData=None, outputDirPrefix="", ntDatabaseFileList=None, noOfTotalSequences=None, \
        transferOutput=True, makeBlastDBJob=None):
        """
        2012.5.24
        """
        
        sys.stderr.write("Adding blast jobs for %s input ... "%(len(inputData.jobDataLs)))
        no_of_jobs= 0
        
        topOutputDir = "%sBlast"%(outputDirPrefix)
        topOutputDirJob = self.addMkDirJob(outputDir=topOutputDir)
        no_of_jobs += 1
        
        allBlastResultFile = File(os.path.join(topOutputDir, 'blast.tsv'))
        allBlastMergeJob = self.addStatMergeJob(statMergeProgram=self.mergeSameHeaderTablesIntoOne, \
            outputF=allBlastResultFile, transferOutput=transferOutput, parentJobLs=[topOutputDirJob])
        no_of_jobs += 1
        
        ntDatabaseFile = ntDatabaseFileList[0]
        returnData = PassingData()
        returnData.jobDataLs = []
        
        for jobData in inputData.jobDataLs:
            inputF = jobData.output
            outputFnamePrefix = os.path.join(topOutputDir, os.path.splitext(os.path.basename(inputF.name))[0])
            
            splitFastaJob = self.addSplitFastaFileJob(executable=self.SplitFastaFile, inputFile=inputF, outputFnamePrefix=outputFnamePrefix, \
                        noOfSequencesPerSplitFile=self.blockSize, filenameSuffix=".fasta", noOfTotalSequences=noOfTotalSequences,\
                        parentJobLs=jobData.jobLs + [topOutputDirJob], extraDependentInputLs=[], transferOutput=False, \
                        extraArguments=None, job_max_memory=500)
            no_of_jobs += 1
            for splitFastaOutput in splitFastaJob.outputList:
                outputFile = File('%s.tsv'%(splitFastaOutput.name))
                blastJob = self.addBlastWrapperJob(executable=self.BlastWrapper,
                    inputFile=splitFastaOutput, outputFile=outputFile, \
                    outputFnamePrefix=splitFastaOutput.name , databaseFile=ntDatabaseFile,\
                    maxNoOfMismatches=self.maxNoOfMismatches, minNoOfIdentities=self.minNoOfIdentities, \
                    minIdentityPercentage=self.minIdentityPercentage, blastallPath=self.blastallPath, \
                    parentJobLs=[splitFastaJob, makeBlastDBJob], extraDependentInputLs=ntDatabaseFileList, transferOutput=False, \
                    extraArguments=None, job_max_memory=1000)
                
                #add output to some reduce job
                self.addInputToMergeJob(allBlastMergeJob, \
                                inputF=blastJob.output, parentJobLs=[blastJob])
                no_of_jobs += 1
        sys.stderr.write("%s jobs. Done.\n"%(no_of_jobs))
        #include the tfam (outputList[1]) into the fileLs
        returnData.jobDataLs.append(PassingData(jobLs=[allBlastMergeJob], file=allBlastResultFile, \
                                            fileLs=[allBlastResultFile]))
        return returnData
    
                
        
    def registerCustomExecutables(self, **keywords):
        """
        2012.5.23
        """
        self.registerOneExecutable(name="blastall", path=self.blastallPath,
            clusterSizeMultiplier=0.2)
        self.registerOneExecutable(name="formatdb", path=self.formatdbPath,
            clusterSizeMultiplier=0)
        self.registerOneExecutable(name="BlastWrapper",
            path=os.path.join(self.pymodulePath, 'pegasus/mapper/alignment/BlastWrapper.py'),
            clusterSizeMultiplier=0.1)
        self.registerOneExecutable(name="SplitFastaFile",
            path=os.path.join(self.pymodulePath, 'pegasus/mapper/splitter/SplitFastaFile.py'),
            clusterSizeMultiplier=0.1)
    
    def addMakeBlastDBJob(self, executable=None, inputFile=None, \
        parentJobLs=None, extraDependentInputLs=None, transferOutput=False, \
        extraArguments=None, job_max_memory=500, **keywords):
        """
        2012.10.9 use addGenericJob() instead
        2012.5.24
            untested
        """
        extraOutputLs = []
        for suffix in ['.nin', '.nhr', '.nsq']:	#start from 0
            dbIndexFile = File('%s%s'%(inputFile.name, suffix))
            extraOutputLs.append(dbIndexFile)
        # 2013.07.09
        extraOutputLs.append(File("formatdb.log"))
        
        extraArgumentList = ["-p F"]
        job = self.addGenericJob(executable=executable, inputFile=inputFile, outputFile=None, \
                        parentJobLs=parentJobLs, extraDependentInputLs=extraDependentInputLs, \
                        extraOutputLs=extraOutputLs,\
                        transferOutput=transferOutput, \
                        extraArguments=extraArguments, extraArgumentList=extraArgumentList, \
                        key2ObjectForJob=None,\
                        job_max_memory=job_max_memory)
        return job
    
    def run(self):
        """
        2011-7-11
        """
        self.setup_run()
        
        inputData = PassingData(jobDataLs = [])
        inputFile = self.registerOneInputFile(self.inputFname, folderName=self.pegasusFolderName)
        inputData.jobDataLs.append(PassingData(output=inputFile, jobLs=[]))
        noOfTotalSequences= self.getNoOfSequencesFromFasta(inputFastaFname=self.inputFname)
        
        registerReferenceData = self.registerBlastNucleotideDatabaseFile(ntDatabaseFname=self.databaseFname, \
            input_site_handler=self.input_site_handler)
        ntDatabaseFileList = registerReferenceData.refFastaFList
        ntDatabaseFile = ntDatabaseFileList[0]

        if len(ntDatabaseFileList)<4:	#some nt-database index file is missing
            sys.stderr.write("Adding blast-db-making job...")
            makeBlastDBJob = self.addMakeBlastDBJob(executable=self.formatdb,\
                inputFile=ntDatabaseFile, transferOutput=True)
            #add the index files to the ntDatabaseFileList
            ntDatabaseFileList = [ntDatabaseFile] + makeBlastDBJob.outputList
            sys.stderr.write(".\n")
        else:
            makeBlastDBJob = None
        
        self.addJobs(inputData=inputData, outputDirPrefix=self.pegasusFolderName,
            ntDatabaseFileList=ntDatabaseFileList, \
            noOfTotalSequences=noOfTotalSequences, \
            transferOutput=True, makeBlastDBJob=makeBlastDBJob)
        
        self.end_run()

if __name__ == '__main__':
    main_class = BlastWorkflow
    po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
    instance = main_class(**po.long_option2value)
    instance.run()
