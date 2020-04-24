#!/usr/bin/env python3
"""
input_path: 
    if directory, recursively find all fastq/bam files,
    if file, each line is the path to bam.
sample_sheet: (map fastq/bam files to individual information).

this program will
    1. query table individual if the sample is already in db. Create a new entry in db if not.
    2. query table individual_sequence if the sample's sequence is already in db.
        Create one new entry if not.
        Update individual_sequence.path if it's none
    3. If commit, the db actions would be committed
        and files will moved to db-affiliated file-system.
    4. Output the whole unpack workflow to an xml output file.

This program has to be run on the pegasus submission host.
Option "--commit" commits the db transaction.

1. The db commandline arguments will be passed to db-interacting jobs.
    Make sure all computing nodes have access to the db.
2. The workflow must be run on nodes with access to the db and the db-affiliated filesystem.

Examples:
    # run the program on crocea and output a local condor workflow
    %s -i ~/NetworkData/vervet/VRC/ 
        -t /u/home/eeskintmp/polyacti/NetworkData/vervet/db/ 
        --sample_sheet ~/script/vervet/data/VRC_sequencing_20110802.tsv
        -u yh --commit -z dl324b-1.cmb.usc.edu -o /tmp/condorpool.xml

    #2011-8-26	generate a list of all bam file physical paths through find.
    #  (doing this because they are not located on crocea)
    find NetworkData/vervet/raw_sequence/ -name *.bam  > vervet/raw_sequence/bamFileList.txt
    # run the program on the crocea and output a hoffman2 workflow. (with db commit)
    %s
        -i ~/vervet/raw_sequence/bamFileList.txt
        -e /u/home/eeskin/polyacti/
        -t /u/home/eeskin/polyacti/NetworkData/vervet/db/ -u yh
        --sample_sheet xfer.genome.wustl.edu/gxfer3/46019922623327/Vervet_12_4X_README.tsv
        -z dl324b-1.cmb.usc.edu -l hoffman2
        -o unpackAndAdd12_2007Monkeys2DB_hoffman2.xml
        --commit 
    
    # 20120430 run on hcondor, to import McGill 1X data (-y2), (-e) is not necessary
    #   because it's running on hoffman2 and can recognize home folder.
    #   --needSSHDBTunnel means it needs sshTunnel for db-interacting jobs.
    %s
        -i raw_sequence/McGill96_1X/ -z localhost -u yh -j hcondor -l hcondor --commit
        -o dags/AddReads2DB/unpackMcGill96_1X.xml -y2 --needSSHDBTunnel 
        -D NetworkData/vervet/db/ -t NetworkData/vervet/db/
        -e /u/home/eeskin/polyacti
    
    # 20120602 add 18 south-african monkeys RNA read data (-y3),
    #  sequenced by Joe DeYoung's core (from Nam's folder),
    #  later manually changed its tissue id to distinguish them from DNA data (below).
    %s -i SIVpilot/by.Charles.Demultiplexed/
        -z localhost -u yh -j hcondor -l hcondor
        --commit -o dags/AddReads2DB/unpack20SouthAfricaSIVmonkeys.xml -y3 --needSSHDBTunnel
        -D /u/home/eeskin/polyacti/NetworkData/vervet/db/
        -t /u/home/eeskin/polyacti/NetworkData/vervet/db/ -e /u/home/eeskin/polyacti 
        --sample_sheet by.Charles.Demultiplexed/sampleIds.txt
        --minNoOfReads 4000000

    # 20120603 add 24 south-african monkeys DNA read data (-y4),
    #  sequenced by Joe DeYoung's core (from Nam's folder)
    #   --minNoOfReads 4000000
    %s
        -i ~namtran/panasas/Data/HiSeqRaw/Ania/SIVpilot/LowpassWGS/Demultiplexed/
        -z localhost -u yh -j hcondor -l hcondor --commit
        -o dags/unpack20SouthAfricaSIVmonkeysDNA.xml
        -y4 --needSSHDBTunnel -D ～/NetworkData/vervet/db/
        -t ～/NetworkData/vervet/db/ -e ～ --minNoOfReads 4000000
        
    # 2017.04.28 added TCGA sequences (.bam) into db
    %s -i /y/Sunset/tcga/HNSC_TCGA/
        --hostname pdc -u huangyu -j condor -l condor 
        -o dags/unpackTCGAHNSCSamples.xml -y6
        --tissueSourceSiteFname /y/Sunset/tcga/tcga_code_tables/tissueSourceSite.tsv 
        --minNoOfReads 8000000 --dbname pmdb -k sunset --commit

"""
import sys, os
__doc__ = __doc__%(sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0])

import copy, re, csv
import pandas as pd
from pegaflow.DAX3 import File
from palos import ProcessOptions, PassingData, utils
from palos.io.MatrixFile import MatrixFile
from palos.pegasus.AbstractWorkflow import AbstractWorkflow
from palos.db import SunsetDB
import logging

ParentClass=AbstractWorkflow

class ImportIndividualSequence2DB(ParentClass):
    __doc__ = __doc__
    option_default_dict = copy.deepcopy(ParentClass.option_default_dict)
    option_default_dict.update(ParentClass.db_option_dict)
    option_default_dict.update({
        ('input_path', 1, ): ['', 'i', 1, 'If it is a folder, get all .bam/.sam/.fastq files recursively. '
            'If it is a file, every line should be a path to an input file.', ],\
        ('SplitReadFileJarPath', 1, ): ["%s/script/picard/dist/SplitReadFile.jar", '', 1, 
            'path to the SplitReadFile jar', ],\
        ('picard_path', 1, ): ["%s/script/picard.broad/build/libs/picard.jar", '', 1, 
            'path to the new picard jar', ],\
        ('sample_sheet', 0, ): ['', '', 1, 
            'a tsv-format file detailing the corresponding sample ID of each input file.', ],\
        ('minNoOfReads', 1, int): [8000000, '', 1, 'The minimum number of reads in each split fastq file. '
            'The max number of reads in a split file is 2*minNoOfReads.', ],\
        ("sequencer_name", 0, ): ["", '', 1, 'sequencing center of TCGA. parsed from TCGA bacode.'],\
        ("sequence_type_name", 1, ): ["PairedEnd", '', 1,
            'isq.sequence_type_id table column: SequenceType.short_name.'],\
        ("sequence_format", 1, ): ["fastq", 'f', 1, 'Fasta, fastq, etc.'],\
        ("tissueSourceSiteFname", 0, ): ["", '', 1, 'TCGA tissue source site file'],\
        ('inputType', 1, int): [1, 'y', 1, 'input type. 1: TCGA bam files; 2: HCC1187 bam files', ],\
        ('study_id', 1, int): [None, '', 1, 'Field study.id from db, used to group individuals.', ],\
        ('site_id', 1, int): [None, '', 1, 'Field site.id from db, used to group individuals.', ],\
        })
    
    def __init__(self,  **keywords):
        """
        2011-8-3
        """
        # Insert before ParentClass.__init__()
        self.pathToInsertHomePathList.extend(
            ['SplitReadFileJarPath', 'picard_path',]
            )
        ParentClass.__init__(self, **keywords)
    
    def connectDB(self):
        """
        Overwrite the parent function.
        Called in the end of Workflow.__init__().
        Establish the db_main connection for all derivative classes.
        """
        ParentClass.connectDB(self)
        self.db_main = SunsetDB.SunsetDB(drivername=self.drivername, db_user=self.db_user,
                    db_passwd=self.db_passwd, hostname=self.hostname, \
                    dbname=self.dbname, schema=self.schema)
        self.db_main.setup(create_tables=False)

        if not self.data_dir:
            self.data_dir = self.db_main.data_dir
        if not self.local_data_dir:
            self.local_data_dir = self.db_main.data_dir
    
    def readSampleSheet_WUSTLBam(self, sample_sheet):
        """
        20110803 from WUSTL, the sample_sheet looks like this:
    Only the base filename is required, you can provide more as long as the whole path is unique.
Example ("Library" and "Bam Path" are required):

#	FlowCell	Lane	Index Sequence	Library	Common Name	Bam Path	MD5
1	64J6AAAXX	1	VCAC-2007002-1-lib1	African	Green	Monkey	gerald_64J6AAAXX_1.bam	gerald_64J6AAAXX_1.bam.md5
2	64J6AAAXX	2	VCAC-2007006-1-lib1	African	Green	Monkey	gerald_64J6AAAXX_2.bam	gerald_64J6AAAXX_2.bam.md5

        """
        sys.stderr.write("Getting bamBaseFname2SampleID dictionary ...")
        bamBaseFname2SampleID = {}
        reader = csv.reader(open(sample_sheet), delimiter='\t')
        header = reader.next()
        col_name2index = utils.getColName2IndexFromHeader(header, skipEmptyColumn=True)
        sampleIDIndex = col_name2index.get("Library")
        if sampleIDIndex is None:
            sampleIDIndex = col_name2index.get("library")
        bamFnameIndex = col_name2index.get("Bam Path")
        if bamFnameIndex is None:
            bamFnameIndex = col_name2index.get("BAM Path")
        if bamFnameIndex is None:
            bamFnameIndex = col_name2index.get("BAM")
        if bamFnameIndex is None:
            bamFnameIndex = col_name2index.get("bam pathway")
        # i.e. VCAC-2007002-1-lib1
        #sampleIDPattern = re.compile(r'\w+-(\w+)-\d+-\w+')
        # 2012.5.29 i.e. VCAC-VGA00006-AGM0075-lib1 ,
        sampleIDPattern = re.compile(r'\w+-(\w+)-\w+-\w+')
        # VCAC-VZC1014-AGM0055-lib1, VCAC-1996031-VRV0265-lib2a,
        #  VCAC-VKD7-361-VKD7-361-lib1 (VKD7 is to be taken),
        for row in reader:
            code = row[sampleIDIndex]
            pa_search = sampleIDPattern.search(code)
            if pa_search:
                code = pa_search.group(1)
            else:
                sys.stderr.write("Warning: could not parse sample ID from %s. Ignore.\n"%(code))
                continue
            bamFname = row[bamFnameIndex]
            bamBaseFname = os.path.split(bamFname)[1]
            bamBaseFname2SampleID[bamBaseFname] = code
        sys.stderr.write("%s entries.\n"%(len(bamBaseFname2SampleID)))
        return bamBaseFname2SampleID

    def getAllBamFiles(self, inputDir, bamFnameLs=None):
        """
        2011-8-3
            recursively going through the directory to get all bam files
            
        """
        for inputFname in os.listdir(inputDir):
            #get the absolute path
            input_path = os.path.join(inputDir, inputFname)
            if os.path.isfile(input_path):
                prefix, suffix = os.path.splitext(inputFname)
                if suffix=='.bam' or suffix=='.sam':
                    bamFnameLs.append(input_path)
            elif os.path.isdir(input_path):
                self.getAllBamFiles(input_path, bamFnameLs)
        
    def addIndividualSequence(self, db_main=None, 
        code: str=None, name: str=None, 
        tax_id: int =9606,
        study_name: str=None, study_id: int = None, 
        site_id: int =None, 
        tissue_name: str=None, tissue_id: int =None, 
        sequence_batch_id: int =None,
        sequencer_name: str ='HiSeq', sequence_type_name: str ='PairedEnd', 
        sequence_format: str ='fastq', isq_comment: str=None, 
        path_to_original_sequence=None, data_dir: str =None):
        """
        add individual and then add individual_sequence
        """
        individual = db_main.getIndividual(code=code, name=name, tax_id=tax_id, 
            study_name=study_name, study_id=study_id, site_id=site_id)
        individual_sequence = db_main.getIndividualSequence(individual_id=individual.id, 
            sequencer_name=sequencer_name, 
            sequence_type_name=sequence_type_name, sequence_format=sequence_format, 
            path_to_original_sequence=path_to_original_sequence, coverage=None,
            sequence_batch_id=sequence_batch_id, 
            tissue_name=tissue_name, tissue_id=tissue_id,
            filtered=0, version=None, is_contaminated=0, outdated_index=0, 
            data_dir=data_dir, 
            comment=isq_comment)
        
        return individual_sequence
    
    def registerJars(self):
        ParentClass.registerJars(self)
        
        self.registerOneJar(name="SplitReadFileJar", path=self.SplitReadFileJarPath)
        self.registerOneJar(name="PicardJar", path=self.picard_path)

    def registerExecutables(self):
        """
        20170503 No pegasus job clustering for these executables 
            because some of these jobs are certain to fail and
             will affect other good samples' downstream jobs. 
        2012.1.3
        """
        ParentClass.registerExecutables(self)

        self.registerOneExecutable(
            path=os.path.join(self.pymodulePath, 'mapper/splitter/SplitReadFileWrapper.sh'), \
            name='SplitReadFileWrapper', clusterSizeMultiplier=0)
        self.registerOneExecutable(
            path=os.path.join(self.pymodulePath, 'mapper/computer/VerifyFileMD5Sum.py'), \
            name='VerifyFileMD5Sum', clusterSizeMultiplier=0)
        self.registerOneExecutable(
            path=os.path.join(self.pymodulePath, 'db/import/RegisterAndMoveSplitSequenceFiles.py'), \
            name='RegisterAndMoveSplitSequenceFiles', clusterSizeMultiplier=0)
        self.registerOneExecutable(
            path=os.path.join(self.pymodulePath, 'db/import/RegisterIndividualSequence2DB.py'), \
            name='RegisterIndividualSequence2DB', clusterSizeMultiplier=0)
        self.registerOneExecutable(path=self.javaPath, \
            name='SamToFastqJava', clusterSizeMultiplier=0)
        
    def addConvertBamToFastqAndGzipJob(
        self, executable=None, inputF=None, outputFnamePrefix=None, \
        parentJobLs=None, job_max_memory=8000, walltime = 800, 
        extraDependentInputLs=None, transferOutput=False, **keywords):
        """
        2013.04.03 use addGenericJavaJob()
        2012.1.3
            walltime is in minutes (max time allowed on hoffman2 is 24 hours).
            The executable should be convertBamToFastqAndGzip.
            
        """
        if extraDependentInputLs is None:
            extraDependentInputLs = []
        if executable is None:
            executable = self.SamToFastqJava
        output1 = File("%s_1.fastq.gz"%(outputFnamePrefix))
        output2 = File("%s_2.fastq.gz"%(outputFnamePrefix))
        output_unpaired = File("%s_unpaired.fastq.gz"%(outputFnamePrefix))
        extraOutputLs= [output1, output2, output_unpaired]
        extraArgumentList = ["F=", output1, "F2=", output2, "UNPAIRED_FASTQ=", output_unpaired, \
            "VALIDATION_STRINGENCY=LENIENT", "INCLUDE_NON_PF_READS=true"]
        
        job = self.addGenericJavaJob(executable=self.SamToFastqJava, 
            jarFile=self.PicardJar, \
            inputFile=inputF, inputArgumentOption="INPUT=", \
            outputFile=None, outputArgumentOption="",\
            parentJobLs=parentJobLs, transferOutput=transferOutput, \
            frontArgumentList=["SamToFastq"], \
            extraArgumentList=extraArgumentList, extraOutputLs=extraOutputLs, \
            extraDependentInputLs=extraDependentInputLs, \
            job_max_memory=job_max_memory, walltime=walltime, **keywords)
        
        job.output1 = output1
        job.output2 = output2
        job.output_unpaired = output_unpaired
        return job
    
    def addSplitReadFileJob(self, executable=None, \
        inputF=None, outputFnamePrefix=None, outputFnamePrefixTail="",\
        minNoOfReads=5000000, logFile=None, parentJobLs=None, \
        job_max_memory=4000, walltime = 800, \
        extraDependentInputLs=None, transferOutput=False, **keywords):
        """
        Executable is shell/SplitReadFileWrapper.sh,
            which calls "wc -l" to count the number of reads beforehand to 
            derive a proper minNoOfReads, avoiding files with too few reads.
        Run SplitReadFile and generate the output directly into the db-affiliated folders.
        A log file is generated and registered for transfer (so that pegasus won't skip it).
        Walltime is in minutes (max time allowed on hoffman2 is 24 hours).
        """
        if executable is None:
            executable = self.SplitReadFileWrapper
        if extraDependentInputLs is None:
            extraDependentInputLs = []
        frontArgumentList = [self.javaPath, repr(job_max_memory), self.SplitReadFileJar]
        extraArgumentList = [outputFnamePrefix, repr(minNoOfReads)]
        extraDependentInputLs.append(self.SplitReadFileJar)
        if logFile:
            extraArgumentList.append(logFile)
        
        job = self.addGenericJob(executable=executable, 
            inputFile=inputF, inputArgumentOption="", \
            outputFile=None, outputArgumentOption="-o", \
            parentJob=None, parentJobLs=parentJobLs, 
            extraDependentInputLs=extraDependentInputLs, \
            extraOutputLs=[logFile], transferOutput=transferOutput, \
            frontArgumentList=frontArgumentList, \
            extraArgumentList=extraArgumentList, \
            job_max_memory=job_max_memory, walltime=walltime, **keywords)
        return job
    
    def addRegisterIndividualSequence2DBJob(self, executable=None, \
        inputFile=None, individual_id=None, outputFile=None, \
        parentJobLs=None, job_max_memory=100, walltime = 60, commit=0, \
        extraDependentInputLs=None, extraArguments=None, \
        transferOutput=False, sshDBTunnel=1, **keywords):
        """
        20170428
            walltime is in minutes (max time allowed on hoffman2 is 24 hours).
        """
        if extraDependentInputLs is None:
            extraDependentInputLs = []
        if executable is None:
            executable = self.RegisterIndividualSequence2DB
        extraArgumentList =['--individual_id %s'%individual_id]
        for name, value in keywords.items():
            if value is not None and value != "":
                extraArgumentList.append("--%s %s"%(name, value))
        job = self.addGenericFile2DBJob(executable=executable, \
            inputFile=inputFile, inputArgumentOption="-i", \
            outputFile=outputFile, outputArgumentOption="-o", inputFileList=None, \
            data_dir=None, commit=commit,\
            parentJobLs=parentJobLs, extraDependentInputLs=extraDependentInputLs, \
            extraOutputLs=None, transferOutput=transferOutput, \
            extraArguments=extraArguments, extraArgumentList=extraArgumentList, \
            job_max_memory=job_max_memory,  sshDBTunnel=sshDBTunnel, walltime=walltime,\
            key2ObjectForJob=None, objectWithDBArguments=self)
        return job
    
    def addRegisterAndMoveSplitFileJob(self, 
        inputFile=None,  inputDir=None, 
        outputDir:str=None, relativeOutputDir: str=None,
        individual_sequence_id:int=None, \
        individual_sequence_file_raw_id:int=None,
        original_file=None,
        logFile=None, library=None, mate_id=None, \
        parentJobLs=None, 
        job_max_memory:int=100, walltime:int = 60, \
        commit=0, sequence_format:str='fastq',\
        extraDependentInputLs=None, extraArguments=None, \
        transferOutput:bool=False, sshDBTunnel=1):
        """
        walltime is in minutes (max time allowed on hoffman2 is 24 hours).
        """
        if extraDependentInputLs is None:
            extraDependentInputLs = []
        extraArgumentList =['--inputDir', inputDir,
            '--sequence_format', sequence_format]
        if outputDir:
            extraArgumentList.extend(['--outputDir', outputDir])
        if relativeOutputDir:
            extraArgumentList.extend(['--relativeOutputDir', relativeOutputDir])
        if individual_sequence_id:
            extraArgumentList.extend(['--individual_sequence_id', individual_sequence_id])
        if individual_sequence_file_raw_id:
            extraArgumentList.extend(
                ['--individual_sequence_file_raw_id', individual_sequence_file_raw_id])
        if original_file:
            extraArgumentList.extend([
                '--original_file_path', original_file
            ])
            extraDependentInputLs.append(original_file)
        if library:
            extraArgumentList.extend(['--library', library])
        if mate_id:
            extraArgumentList.append('--mate_id %s'%(mate_id))
        
        job = self.addGenericFile2DBJob(executable=self.RegisterAndMoveSplitSequenceFiles, \
            inputFile=inputFile, inputArgumentOption="-i", \
            data_dir=None, logFile=logFile, commit=commit,\
            parentJobLs=parentJobLs, extraDependentInputLs=extraDependentInputLs, \
            extraOutputLs=None, transferOutput=transferOutput, \
            extraArguments=extraArguments, extraArgumentList=extraArgumentList, \
            job_max_memory=job_max_memory, walltime=walltime,\
            sshDBTunnel=sshDBTunnel,
            key2ObjectForJob=None, objectWithDBArguments=self)
        return job
    
    def getInputPathLsFromInput(self, input_path=None, suffixSet=set(['.fastq']), fakeSuffix='.gz'):
        """
        2012.4.30
            this function supercedes self.getAllBamFiles() and it's more flexible.
        """
        input_path_list = []
        if os.path.isdir(input_path):
            self.getFilesWithSuffixFromFolderRecursive(inputFolder=input_path, 
                suffixSet=suffixSet, 
                fakeSuffix=fakeSuffix, return_path_list=input_path_list)
        elif os.path.isfile(input_path):
            inf = open(input_path)
            for line in inf:
                input_path_list.append(line.strip())
            del inf
        else:
            sys.stderr.write("%s is neither a folder nor a file.\n"%(input_path))
            sys.exit(4)
        return input_path_list

    def getSampleID2FastqObjectLsForSouthAfricanDNAFastQ(self, fastq_path_ls=None):
        """        
        In pairs like this:
            VSAA2015_1.fastq.gz
            VSAA2015_2.fastq.gz
        """
        sys.stderr.write("Parsing sample_id2fastq_obj_ls from %s files ..."%(len(fastq_path_ls)))
        sample_id2fastq_obj_ls = {}
        import random
        filename_pattern = re.compile(r'(?P<sample_id>[a-zA-Z0-9]+)_(?P<mate_id>\d).fastq')
        counter = 0
        real_counter = 0
        #McGill's library ID , 7_Index-11, is not unique enough.
        libraryKey2UniqueLibrary = {}
        for fastq_path in fastq_path_ls:
            counter += 1
            search_result = filename_pattern.search(fastq_path)
            if search_result:
                real_counter += 1
                sample_id = search_result.group('sample_id')
                mate_id = search_result.group('mate_id')
                filenameSignature = (sample_id)
                    
                #concoct a unique library ID
                #this combination insures two ends from the same library are grouped together
                libraryKey = (sample_id)
                if libraryKey not in libraryKey2UniqueLibrary:
                    uniqueLibrary = '%s_%s'%(sample_id, repr(random.random())[2:])
                    libraryKey2UniqueLibrary[libraryKey] = uniqueLibrary
                
                uniqueLibrary = libraryKey2UniqueLibrary[libraryKey]
                fastq_obj = PassingData(library=uniqueLibrary, sample_id=sample_id, 
                    mate_id=mate_id, abs_path=fastq_path)
                if sample_id not in sample_id2fastq_obj_ls:
                    sample_id2fastq_obj_ls[sample_id] = []
                sample_id2fastq_obj_ls[sample_id].append(fastq_obj)
            else:
                sys.stderr.write("Error: can't parse sample_id, library, mate_id out of %s.\n"%fastq_path)
                sys.exit(4)
        sys.stderr.write(" %s samples and %s files in the dictionary.\n"%
            (len(sample_id2fastq_obj_ls), real_counter))
        return sample_id2fastq_obj_ls

    def addJobsToImportSouthAfricanDNAFastQ(self, db_main=None, 
            sample_sheet=None, input_path=None, data_dir=None, \
            minNoOfReads=None, commit=None,\
            sequencer_name=None, sequence_type_name=None, sequence_format=None):
        """
        2012.6.2
            input fastq files could be gzipped or not.
            data generated by Joe DeYoung's core, demultiplexed by ICNN.
        """
        fastq_path_ls = self.getInputPathLsFromInput(input_path, suffixSet=set(['.fastq']), 
        fakeSuffix='.gz')
        sample_id2fastq_obj_ls = self.getSampleID2FastqObjectLsForSouthAfricanDNAFastQ(
            fastq_path_ls=fastq_path_ls)
        self.addJobsToSplitAndRegisterFastQ(db_main=db_main, 
            sample_id2fastq_obj_ls=sample_id2fastq_obj_ls, data_dir=data_dir, \
            minNoOfReads=minNoOfReads, commit=commit,\
            sequencer_name=sequencer_name, sequence_type_name=sequence_type_name, 
            sequence_format=sequence_format)

    def readSampleSheetFastQ(self, sample_sheet: str, fastq_path_ls: list):
        """
        20200421
        Only the first two columns are required. The rest is optional.
        But should provide as much as you can.
        Columns of sample_sheet:
            file_path,sample_id,study_name,study_id,site_id,
            sequence_type,sequencer,sequence_batch_id,libarary,
            tissue_name,tissue_id
        """
        logging.info(f"Getting sample_id2data from {sample_sheet} ...")
        filename_pattern = re.compile(r'(?P<other>.*)_(?P<mate_id>\d).(fastq|fq)')
        sample_id2data = {}
        df = pd.read_excel(sample_sheet, header=True)
        no_of_files = 0
        for idx, row in df.iterrows():
            find_path_func = lambda x: x.find(row['file_path']) >= 0
            found_path_ls = list(filter(find_path_func, fastq_path_ls))
            if len(found_path_ls)==1:
                abs_path = found_path_ls[0]
            elif len(found_path_ls)==0:
                logging.warning(f"No file found for {row['file_path']}.")
                continue
            else:
                logging.warning(f">1 files found for {row['file_path']}: "
                    f"{found_path_ls}.")
                continue
            #replace space in sample_id
            sample_id = row['sample_id'].replace(' ', '_')
            if sample_id not in sample_id2data:
                sample_id2data[sample_id] = PassingData(
                    sample_id=sample_id,
                    study_name=row['study_name'], study_id=row['study_id'], 
                    site_id=row['site_id'], fastq_obj_ls=[])
            search_result = filename_pattern.search(row['file_path'])
            mate_id = search_result.group('mate_id')

            fastq_obj = PassingData(
                abs_path = abs_path, 
                library=None, mate_id = mate_id,
                sequence_type=row['sequence_type'],
                sequencer=row['sequencer'], 
                sequence_batch_id=row['sequence_batch_id'],
                tissue_name=row['tissue_name'], tissue_id=row['tissue_id'],
                condition=row['condition'],
                )
            sample_id2data[sample_id].fastq_obj_ls.append(fastq_obj)
            no_of_files +=1
        
        logging.info(f"{len(sample_id2data)} samples with {no_of_files} files.")
        return sample_id2data

    def addJobsToImportFastQ(self, db_main=None, \
        sample_sheet=None, input_path=None, data_dir=None, \
        minNoOfReads=None, commit=None):
        """
        a general importer
        """
        fastq_path_ls = self.getInputPathLsFromInput(
            input_path, suffixSet=set(['.fastq', '.fq']), fakeSuffix='.gz')
        sample_id2data = self.readSampleSheetFastQ(sample_sheet, fastq_path_ls)
        
        logging.info("Adding split-read & register jobs ...")
        filenameKey2PegasusFile = {}
        for sample_id, sample_data in sample_id2data.items():
            individual = db_main.getIndividual(code=sample_id, tax_id=9606, 
                study_name=sample_data.study_name, study_id=sample_data.study_id, 
                site_id=sample_data.site_id)
            fastq_obj_1 = sample_data.fastq_obj_ls[0]
            individual_sequence = db_main.getIndividualSequence(
                individual_id=individual.id, 
                sequencer_name=fastq_obj_1.sequencer,
                sequence_type_name=fastq_obj_1.sequence_type, 
                sequence_format=self.sequence_format, 
                path_to_original_sequence=fastq_obj_1.abs_path, 
                coverage=None,
                sequence_batch_id=fastq_obj_1.sequence_batch_id, 
                tissue_name=fastq_obj_1.tissue_name, 
                tissue_id=fastq_obj_1.tissue_id,
                condition_name=fastq_obj_1.condition,
                filtered=0, version=None, is_contaminated=0, outdated_index=0, 
                data_dir=data_dir)
            
            sequenceAbsDir = os.path.join(data_dir, individual_sequence.path)
            createSequenceAbsDirJob = self.addMkDirJob(outputDir=sequenceAbsDir)
            
            splitOutputDir = f'{individual_sequence.id}'
            #One directory containing split files from both mates is fine,
            # as RegisterAndMoveSplitSequenceFiles could pick up.
            splitOutputDirJob = self.addMkDirJob(outputDir=splitOutputDir)
            for fastq_obj in sample_data.fastq_obj_ls:
                library = fastq_obj.library
                mate_id = fastq_obj.mate_id
                fastq_path = fastq_obj.abs_path
                fastqFile = self.registerOneInputFile(fastq_path)
                fastqFile.sample_id = sample_id
            
                splitFastQFnamePrefix = os.path.join(
                    splitOutputDir, f'{individual_sequence.id}_{mate_id}')
                logFile = File(f'{individual_sequence.id}_{mate_id}.split.log')
                splitReadFileJob1 = self.addSplitReadFileJob(
                    inputF=fastqFile, outputFnamePrefix=splitFastQFnamePrefix, \
                    outputFnamePrefixTail="", minNoOfReads=minNoOfReads, \
                    logFile=logFile, parentJobLs=[splitOutputDirJob], \
                    job_max_memory=4000, walltime = 800, \
                    extraDependentInputLs=None, transferOutput=True)
                
                logFile = File(f'{individual_sequence.id}_{mate_id}.register.log')
                extraArgumentList=['--inputDir', splitOutputDir,
                    '--outputDir', sequenceAbsDir,
                    '--relativeOutputDir', individual_sequence.path,
                    '--sequence_format', self.sequence_format, 
                    '--individual_sequence_id', individual_sequence.id,
                    '--original_file_path', fastqFile,
                    ]
                if library:
                    extraArgumentList.extend(['--library', library])
                if mate_id:
                    extraArgumentList.append(f'--mate_id {mate_id}')
                
                registerJob = self.addGenericFile2DBJob(
                    executable=self.RegisterAndMoveSplitSequenceFiles, \
                    logFile=logFile, commit=commit,\
                    parentJobLs=[splitReadFileJob1, createSequenceAbsDirJob], 
                    extraDependentInputLs=[fastqFile], \
                    extraOutputLs=None, transferOutput=True, \
                    extraArgumentList=extraArgumentList,
                    job_max_memory=100, walltime=60,\
                    sshDBTunnel=self.needneedSSHDBTunnel,
                    objectWithDBArguments=self)
            
        logging.info(f"{self.no_of_jobs} jobs.")
    
    def addJobsToImportUNGCVervetFastQ(self, db_main=None, \
        sample_sheet=None, input_path=None, data_dir=None, \
        minNoOfReads=None, commit=None,\
        sequencer_name=None, sequence_type_name=None, sequence_format=None):
        """
        2013.04.04
        UNGC = UCLA Neuroscience Genomics Core.
        Data generated by Joe DeYoung's core, demultiplexed by ICNN.
        The input fastq files could be gzipped or not.
        """
        sample_id2data = self.readSampleSheetFastQ_UNGC(sample_sheet)
        fastq_path_ls = self.getInputPathLsFromInput(input_path, suffixSet=set(['.fastq']),
            fakeSuffix='.gz')
        sample_id2fastq_obj_ls = self.getSampleID2FastqObjectLsForUNGCFastQ(
            fastq_path_ls=fastq_path_ls, \
            sample_id2data=sample_id2data)
        self.addJobsToSplitAndRegisterFastQ(db_main=db_main, \
            sample_id2fastq_obj_ls=sample_id2fastq_obj_ls, \
            data_dir=data_dir, minNoOfReads=minNoOfReads, commit=commit,\
            sequencer_name=sequencer_name, sequence_type_name=sequence_type_name, 
            sequence_format=sequence_format)

    def readSampleSheetFastQ_UNGC(self, sample_sheet=None):
        """
        2013.04.04
        Format is like this from UNGC  = UCLA Neuroscience Genomics Core:

FCID	Lane	sample ID	sample code	sample name	Index	Description	SampleProject
D1HYNACXX	1	Ilmn Human control pool ( 4plex)	IP1		INDEX IS UNKNOWN	2013-029A
D1HYNACXX	2	UNGC Human Sample 1	S1	AS001A	ATTACTCG	TruSeq DNA PCR Free beta kit	2013-029A
        """
        print(f"Getting sample_id2data from {sample_sheet} ...", flush=True)
        sample_id2data = {}
        reader = MatrixFile(sample_sheet, openMode='r', delimiter=',')
        reader.constructColName2IndexFromHeader()
        sampleIDIndex = reader.getColIndexGivenColHeader("sample ID")
        sampleNameIndex = reader.getColIndexGivenColHeader("sample name")
        libraryIndexIndex = reader.getColIndexGivenColHeader("Index")
        
        for row in reader:
            #2013.04.04 stupid quirks
            sample_id = row[sampleIDIndex].replace(' ', '_')
            sampleName = row[sampleNameIndex]
            libraryIndex = row[libraryIndexIndex]
            if sample_id not in sample_id2data:
                sample_id2data[sample_id] = PassingData(sampleName=sampleName, 
                    libraryIndexList=[])
            if sampleName!=sample_id2data[sample_id].sampleName:
                logging.error("sample_id %s is associated with two different sample names (%s, %s)."\
                    %(sample_id, sampleName, sample_id2data[sample_id].sampleName))
                sys.exit(4)
            sample_id2data[sample_id].libraryIndexList.append(libraryIndex)
        
        print(f"{len(sample_id2data)} entries.", flush=True)
        return sample_id2data
    
    def getSampleID2FastqObjectLsForUNGCFastQ(self, fastq_path_ls=None, 
        sample_id2data=None):
        """
        2013.04.05 bugfix. this ensures two ends from same library have same library.
        2013.04.04
        UNGC  = UCLA Neuroscience Genomics Core
        In pairs like this:
            UNGC_Vervet_Sample_11_1.fastq.gz
            UNGC_Vervet_Sample_11_2.fastq.gz
        """
        sys.stderr.write("Parsing sample_id2fastq_obj_ls from %s files ..."%(len(fastq_path_ls)))
        sample_id2fastq_obj_ls = {}
        import random
        filename_pattern = re.compile(r'(?P<sample_id>[\w\d]+)_(?P<mate_id>\d).fastq')
        counter = 0
        real_counter = 0
        libraryKey2UniqueLibrary = {}
        #McGill's library ID , 7_Index-11, is not unique enough.
        for fastq_path in fastq_path_ls:
            counter += 1
            search_result = filename_pattern.search(fastq_path)
            if search_result:
                real_counter += 1
                sample_id = search_result.group('sample_id')
                mate_id = search_result.group('mate_id')
                if sample_id in sample_id2data:
                    sample_data = sample_id2data.get(sample_id)
                    sample_id = sample_data.sampleName
                    library = '_'.join(sample_data.libraryIndexList)
                    
                    libraryKey = (sample_id, library)
                    if libraryKey not in libraryKey2UniqueLibrary:
                        #this ensures two ends from the same library is marked as the same library.
                        uniqueLibrary = '%s_%s_%s'%(sample_id, library, repr(random.random())[2:])
                        libraryKey2UniqueLibrary[libraryKey] = uniqueLibrary
                    else:
                        print(f"libraryKey {libraryKey} of {fastq_path} is already in "
                        f"libraryKey2UniqueLibrary with unique library = {libraryKey2UniqueLibrary[libraryKey]}.",
                        flush=True)
                    
                    uniqueLibrary = libraryKey2UniqueLibrary[libraryKey]
                    fastq_obj = PassingData(library=uniqueLibrary, sample_id=sample_id, 
                        mate_id=mate_id, abs_path=fastq_path)
                    if sample_id not in sample_id2fastq_obj_ls:
                        sample_id2fastq_obj_ls[sample_id] = []
                    sample_id2fastq_obj_ls[sample_id].append(fastq_obj)
                else:
                    print(f"sample_id {sample_id} not in sample_id2data.", flush=True)
            else:
                print(f"ERROR: can't parse sample_id, mate_id out of {fastq_path}", flush=True)
                raise
        print(f" {len(sample_id2fastq_obj_ls)} samples and {real_counter} files in the dictionary.",
            flush=True)
        return sample_id2fastq_obj_ls
    
    
    def readSampleSheet_SouthAfricanRNAFastQ(self, sample_sheet=None):
        """
        2012.6.2 sample_sheet is tab-delimited, looks like this
            VSAC1012_R      ATCACG  sample2
            VSAF1009_R      ATCACG  sample1
            VSAB2009_RN     TAGCTT  sample2
            VSAB2011_R      TAGCTT  sample1
            VSAB3001_RN     GGCTAC  sample2
            VSAB5004_R      GGCTAC  sample1
            VSAA2015_RN     CTTGTA  sample2
        """
        sys.stderr.write("Getting filenameSignature2SampleID dictionary from %s ..."%(sample_sheet))
        filenameSignature2SampleID = {}
        reader = csv.reader(open(sample_sheet), delimiter='\t')
        sampleIDIndex = 0
        folderNameIndex = 1
        subSampleNameIndex = 2
        sampleIDPattern = re.compile(r'(\w+)_R[\w]*')	# 2012.6.2 VSAA2015_RN or VSAB5004_R
        for row in reader:
            sample_id = row[sampleIDIndex]
            pa_search = sampleIDPattern.search(sample_id)
            if pa_search:
                sample_id = pa_search.group(1)
            else:
                sys.stderr.write("Warning: could not parse sample ID from %s. Ignore.\n"%(sample_id))
                continue
            folderName = row[folderNameIndex]
            subSampleName = row[subSampleNameIndex]
            filenameSignature = (folderName, subSampleName)
            if filenameSignature in filenameSignature2SampleID:
                print("Error: filenameSignature %s already in filenameSignature2SampleID."
                    %(filenameSignature), flush=True)
                sys.exit(3)
            filenameSignature2SampleID[filenameSignature] = sample_id
        sys.stderr.write("%s entries.\n"%(len(filenameSignature2SampleID)))
        return filenameSignature2SampleID
    
    def getSampleID2FastqObjectLsForSouthAfricanRNAFastQ(self, fastq_path_ls=None, 
        filenameSignature2SampleID=None):
        """
        2012.6.2            
            In pairs like this:
            .../GCCAAT/tile_1101_sample1_end1.fastq
            .../GCCAAT/tile_1101_sample1_end2.fastq
        """
        sys.stderr.write("Parsing sample_id2fastq_obj_ls from %s files ..."%(len(fastq_path_ls)))
        sample_id2fastq_obj_ls = {}
        import random
        filename_pattern = re.compile(r'/(?P<folderName>[ACGT]{6})/(?P<library>[\w]+)'
            r'_(?P<subSampleName>sample[12])_end(?P<mate_id>\d).fastq')
        counter = 0
        real_counter = 0
        #McGill's library ID , 7_Index-11, is not unique enough.
        libraryKey2UniqueLibrary = {}
        for fastq_path in fastq_path_ls:
            counter += 1
            search_result = filename_pattern.search(fastq_path)
            if search_result:
                real_counter += 1
                #tile_1101
                library = search_result.group('library')
                folderName = search_result.group('folderName')
                subSampleName = search_result.group('subSampleName')
                mate_id = search_result.group('mate_id')
                filenameSignature = (folderName, subSampleName)
                if filenameSignature in filenameSignature2SampleID:
                    sample_id = filenameSignature2SampleID.get(filenameSignature)
                    
                    #concoct a unique library ID
                    #this combination insures two ends from the same library are grouped together
                    libraryKey = (folderName, library)
                    if libraryKey not in libraryKey2UniqueLibrary:
                        uniqueLibrary = '%s_%s_%s'%(folderName, library, repr(random.random())[2:])
                        libraryKey2UniqueLibrary[libraryKey] = uniqueLibrary
                    
                    uniqueLibrary = libraryKey2UniqueLibrary[libraryKey]
                    fastq_obj = PassingData(library=uniqueLibrary, sample_id=sample_id, 
                        mate_id=mate_id, abs_path=fastq_path)
                    if sample_id not in sample_id2fastq_obj_ls:
                        sample_id2fastq_obj_ls[sample_id] = []
                    sample_id2fastq_obj_ls[sample_id].append(fastq_obj)
                else:
                    sys.stderr.write("%s not in filenameSignature2SampleID.\n"%(filenameSignature))
            else:
                logging.error("Can't parse sample_id, library, mate_id out of %s."%fastq_path)
                sys.exit(4)
        sys.stderr.write(" %s samples and %s files in the dictionary.\n"%(len(sample_id2fastq_obj_ls),
            real_counter))
        return sample_id2fastq_obj_ls
    
    def addJobsToImportSouthAfricanRNAFastQ(self, db_main=None, \
        sample_sheet=None, input_path=None, data_dir=None, \
        minNoOfReads=None, commit=None,\
        sequencer_name=None, sequence_type_name=None, sequence_format=None):
        """
        2012.6.1
        input fastq files could be gzipped or not.
        data generated by Joe DeYoung's core, demultiplexed by ICNN (Charles in particular)
        """
        filenameSignature2SampleID = self.readSampleSheet_SouthAfricanRNAFastQ(sample_sheet)
        fastq_path_ls = self.getInputPathLsFromInput(input_path, \
            suffixSet=set(['.fastq']), fakeSuffix='.gz')
        sample_id2fastq_obj_ls = self.getSampleID2FastqObjectLsForSouthAfricanRNAFastQ(
            fastq_path_ls=fastq_path_ls, \
            filenameSignature2SampleID=filenameSignature2SampleID)
        self.addJobsToSplitAndRegisterFastQ(db_main=db_main, 
            sample_id2fastq_obj_ls=sample_id2fastq_obj_ls, data_dir=data_dir, \
            minNoOfReads=minNoOfReads, commit=commit,\
            sequencer_name=sequencer_name, \
            sequence_type_name=sequence_type_name, sequence_format=sequence_format)		
    
    def getSampleID2FastqObjectLsForMcGillData(self, fastq_path_ls=None):
        """
        2013.1.30
        add fixed monkey ID prefixes (VWP, VGA, etc.) into sampleIDPattern 
        due to new not-all-number monkey IDs.
        
HI.0628.001.D701.VGA00010_R1.fastq.gz  HI.0628.004.D703.VWP00384_R1.fastq.gz  HI.0628.007.D703.VWP10020_R1.fastq.gz
HI.0628.001.D701.VGA00010_R2.fastq.gz  HI.0628.004.D703.VWP00384_R2.fastq.gz  HI.0628.007.D703.VWP10020_R2.fastq.gz
        2012.4.30
        each fastq file looks like 7_Index-11.2006013_R1.fastq.gz, 7_Index-11.2006013_R2.fastq.gz,
            7_Index-10.2005045_replacement_R1.fastq.gz, 7_Index-10.2005045_replacement_R2.fastq.gz
        
        The data from McGill is dated 2012.4.27.
        Each monkey is sequenced at 1X. There are 96 of them. Each library seems to contain 2 samples.
            But each monkey has only 1 library.
        
        8_Index_23.2008126_R1.fastq.gz
        8_Index_23.2008126_R2.fastq.gz
        8_Index_23.2009017_R1.fastq.gz
        8_Index_23.2009017_R2.fastq.gz

        """
        sys.stderr.write("Parsing sample_id2fastq_obj_ls from %s files ..."%(len(fastq_path_ls)))
        sample_id2fastq_obj_ls = {}
        import re, random
        #sampleIDPattern = re.compile(r'(?P<library>[-\w]+)\.(?P<sample_id>\d+)
        # ((_replacement)|(_pool)|())_R(?P<mate_id>\d).fastq.gz')
        sampleIDPattern = re.compile(r'(?P<library>[-\w]+)\.'
            r'(?P<sample_id>((VWP)|(VGA)|(VSA)|())\d+)((_replacement)|(_pool)|())_R(?P<mate_id>\d).fastq.gz')
        counter = 0
        real_counter = 0
        libraryKey2UniqueLibrary = {}
        #McGill's library ID , 7_Index-11, is not unique enough.
        for fastq_path in fastq_path_ls:
            counter += 1
            search_result = sampleIDPattern.search(fastq_path)
            if search_result:
                real_counter += 1
                library = search_result.group('library')
                sample_id = search_result.group('sample_id')
                mate_id = search_result.group('mate_id')
                #concoct a unique library ID
                if library not in libraryKey2UniqueLibrary:
                    libraryKey2UniqueLibrary[library] = '%s_%s'%(library, repr(random.random())[2:])
                uniqueLibrary = libraryKey2UniqueLibrary[library]
                fastq_obj = PassingData(library=uniqueLibrary, \
                    sample_id=sample_id, mate_id=mate_id, abs_path=fastq_path)
                if sample_id not in sample_id2fastq_obj_ls:
                    sample_id2fastq_obj_ls[sample_id] = []
                sample_id2fastq_obj_ls[sample_id].append(fastq_obj)
            else:
                print("Error: can't parse sample_id, library, mate_id out of %s."%
                    fastq_path, flush=True)
                sys.exit(4)
        sys.stderr.write(" %s samples and %s files in the dictionary.\n"%(
            len(sample_id2fastq_obj_ls), real_counter))
        return sample_id2fastq_obj_ls
    
    def addJobsToSplitAndRegisterFastQ(self, db_main=None, 
        sample_id2fastq_obj_ls=None, data_dir=None, \
        minNoOfReads=None, commit=None, \
        sequencer_name=None, sequence_type_name=None, sequence_format=None):
        """
        split out of addJobsToImportMcGillFastQ(), also used in addJobsToImportUNGCVervetFastQ().
        """
        print("Adding split-read & register jobs ...", flush=True)
        filenameKey2PegasusFile = {}
        for sample_id, fastq_obj_ls in sample_id2fastq_obj_ls.items():
            individual_sequence = self.addIndividualSequence(
                db_main=db_main, code=sample_id, \
                site_id=self.site_id, study_id=self.study_id, \
                sequencer_name=sequencer_name, sequence_type_name=sequence_type_name, \
                sequence_format=sequence_format, data_dir=data_dir)
            
            sequenceAbsDir = os.path.join(data_dir, individual_sequence.path)
            createSequenceAbsDirJob = self.addMkDirJob(outputDir=sequenceAbsDir)
            
            splitOutputDir = '%s'%(individual_sequence.id)
            #Same directory containing split files from both mates is fine 
            # as RegisterAndMoveSplitSequenceFiles could pick up.
            splitOutputDirJob = self.addMkDirJob(outputDir=splitOutputDir)
            for fastq_obj in fastq_obj_ls:
                library = fastq_obj.library
                mate_id = fastq_obj.mate_id
                fastq_path = fastq_obj.abs_path
                filenameKey = (library, os.path.basename(fastq_path))
                if filenameKey in filenameKey2PegasusFile:
                    fastqFile = filenameKey2PegasusFile.get(filenameKey)
                    print(f"Error: File {fastq_path} has been registered with sample "
                        f"{fastqFile.sample_id}. Can't happen.", flush=True)
                    sys.exit(3)
                    import pdb
                    pdb.set_trace()
                    continue
                else:
                    fastqFile = self.registerOneInputFile(fastq_path, folderName=library)
                    fastqFile.sample_id = sample_id
                    fastqFile.fastq_obj= fastq_obj
                    filenameKey2PegasusFile[filenameKey] = fastqFile
            
                splitFastQFnamePrefix = os.path.join(splitOutputDir, '%s_%s_%s'%(individual_sequence.id, 
                    library, mate_id))
                logFile = File('%s_%s_%s.split.log'%(individual_sequence.id, library, mate_id))
                splitReadFileJob1 = self.addSplitReadFileJob(
                    inputF=fastqFile, outputFnamePrefix=splitFastQFnamePrefix, \
                    outputFnamePrefixTail="", minNoOfReads=minNoOfReads, \
                    logFile=logFile, parentJobLs=[splitOutputDirJob], \
                    job_max_memory=4000, walltime = 800, \
                    extraDependentInputLs=None, transferOutput=True)
                
                logFile = File('%s_%s_%s.register.log'%(individual_sequence.id, library, mate_id))
                registerJob1 = self.addRegisterAndMoveSplitFileJob(
                    inputDir=splitOutputDir, outputDir=sequenceAbsDir, \
                    relativeOutputDir=individual_sequence.path, logFile=logFile,\
                    individual_sequence_id=individual_sequence.id, \
                    library=library, mate_id=mate_id, \
                    parentJobLs=[splitReadFileJob1, createSequenceAbsDirJob], 
                    job_max_memory=100, walltime = 60, \
                    commit=commit, sequence_format=sequence_format, extraDependentInputLs=None, \
                    transferOutput=True, sshDBTunnel=self.needSSHDBTunnel)
            
        print("%s jobs.\n"%(self.no_of_jobs), flush=True)
        
    def addJobsToImportMcGillFastQ(self, db_main=None, 
        input_path=None, data_dir=None, \
        minNoOfReads=None, commit=None,\
        sequencer_name=None, sequence_type_name=None, sequence_format=None):
        """
        20120430
        Input are .fastq files.
        """
        fastq_path_ls = self.getInputPathLsFromInput(input_path, suffixSet=set(['.fastq']), 
            fakeSuffix='.gz')
        sample_id2fastq_obj_ls = self.getSampleID2FastqObjectLsForMcGillData(fastq_path_ls)
        
        self.addJobsToSplitAndRegisterFastQ(db_main=db_main, \
            sample_id2fastq_obj_ls=sample_id2fastq_obj_ls, data_dir=data_dir, \
            minNoOfReads=minNoOfReads, commit=commit,\
            sequencer_name=sequencer_name, sequence_type_name=sequence_type_name, \
            sequence_format=sequence_format)

    
    def addJobsToImportWUSTLBam(self, db_main=None, sample_sheet=None, 
        input_path=None, data_dir=None, \
        minNoOfReads=None, commit=None,\
        sequencer_name=None, sequence_type_name=None, sequence_format=None):
        """
        20120430
        Input are bam files.
        """
        bamBaseFname2SampleID = self.readSampleSheet_WUSTLBam(sample_sheet)
        bamFnameLs = self.getInputPathLsFromInput(input_path, 
        suffixSet=set(['.bam', '.sam']), fakeSuffix='.gz')
        
        sys.stderr.write("%s total bam files.\n"%(len(bamFnameLs)))
        
        sam2fastqOutputDir = 'sam2fastq'
        sam2fastqOutputDirJob = self.addMkDirJob(outputDir=sam2fastqOutputDir)
        
        for bamFname in bamFnameLs:
            bamBaseFname = os.path.split(bamFname)[1]
            if bamBaseFname not in bamBaseFname2SampleID:
                sys.stderr.write("%s doesn't have code affiliated with.\n"%(bamFname))
                continue
            code = bamBaseFname2SampleID.get(bamBaseFname)
            individual_sequence = self.addIndividualSequence(db_main=db_main, code=code, \
                site_id=self.site_id, study_id=self.study_id,\
                sequencer_name=sequencer_name, sequence_type_name=sequence_type_name, \
                sequence_format=sequence_format, data_dir=data_dir)
            #2012.2.10 stop passing path_to_original_sequence=bamFname to self.addMonkeySequence()
            
            """
            #2012.2.10 temporary, during transition from old records to new ones.
            newISQPath = individual_sequence.constructRelativePathForIndividualSequence()
            newISQPath = '%s_split'%(newISQPath)
            if individual_sequence.path is None or individual_sequence.path !=newISQPath:
                individual_sequence.path = newISQPath
                session.add(individual_sequence)
                session.flush()
            """
            sequenceAbsDir = os.path.join(data_dir, individual_sequence.path)
            createSequenceAbsDirJob = self.addMkDirJob(outputDir=sequenceAbsDir)
            
            bamInputF = self.registerOneInputFile(bamFname)
            
            bamBaseFname = os.path.split(bamFname)[1]
            bamBaseFnamePrefix = os.path.splitext(bamBaseFname)[0]
            library = bamBaseFnamePrefix
            
            outputFnamePrefix = os.path.join(sam2fastqOutputDir, 
                '%s_%s'%(individual_sequence.id, library))
            
            convertBamToFastqAndGzip_job = self.addConvertBamToFastqAndGzipJob(
                executable=self.SamToFastqJava, \
                inputF=bamInputF, outputFnamePrefix=outputFnamePrefix, \
                parentJobLs=[sam2fastqOutputDirJob], job_max_memory=6000, walltime = 800, \
                extraDependentInputLs=None, \
                transferOutput=False)
            
            splitOutputDir = '%s_%s'%(individual_sequence.id, library)
            #same directory containing split files from both mates is fine
            #  as RegisterAndMoveSplitSequenceFiles could pick up.
            splitOutputDirJob = self.addMkDirJob( outputDir=splitOutputDir)
            
            mate_id = 1
            splitFastQFnamePrefix = os.path.join(splitOutputDir, 
                '%s_%s_%s'%(individual_sequence.id, library, mate_id))
            logFile = File('%s_%s_%s.split.log'%(individual_sequence.id, library, mate_id))
            splitReadFileJob1 = self.addSplitReadFileJob(
                executable=self.SplitReadFileWrapper, \
                inputF=convertBamToFastqAndGzip_job.output1, 
                outputFnamePrefix=splitFastQFnamePrefix, \
                outputFnamePrefixTail="", minNoOfReads=minNoOfReads, \
                logFile=logFile, parentJobLs=[convertBamToFastqAndGzip_job, splitOutputDirJob], \
                job_max_memory=6000, walltime = 800, \
                extraDependentInputLs=None, transferOutput=True)
            
            logFile = File('%s_%s_%s.register.log'%(individual_sequence.id, library, mate_id))
            registerJob1 = self.addRegisterAndMoveSplitFileJob(
                executable=self.RegisterAndMoveSplitSequenceFiles, \
                inputDir=splitOutputDir, outputDir=sequenceAbsDir, 
                relativeOutputDir=individual_sequence.path, logFile=logFile,\
                individual_sequence_id=individual_sequence.id, 
                original_file=bamInputF, 
                library=library, mate_id=mate_id, \
                parentJobLs=[splitReadFileJob1, createSequenceAbsDirJob], 
                job_max_memory=100, walltime = 60, \
                commit=commit, sequence_format=sequence_format, extraDependentInputLs=None, \
                transferOutput=True, sshDBTunnel=self.needSSHDBTunnel)
            #handle the 2nd end
            mate_id = 2
            splitFastQFnamePrefix = os.path.join(splitOutputDir, 
                '%s_%s_%s'%(individual_sequence.id, library, mate_id))
            logFile = File('%s_%s_%s.split.log'%(individual_sequence.id, library, mate_id))
            splitReadFileJob2 = self.addSplitReadFileJob(executable=self.SplitReadFileWrapper, \
                inputF=convertBamToFastqAndGzip_job.output2, outputFnamePrefix=splitFastQFnamePrefix, \
                outputFnamePrefixTail="", minNoOfReads=minNoOfReads, \
                logFile=logFile, parentJobLs=[convertBamToFastqAndGzip_job, splitOutputDirJob], \
                job_max_memory=6000, walltime = 800, \
                extraDependentInputLs=None, transferOutput=True)
            
            logFile = File('%s_%s_%s.register.log'%(individual_sequence.id, library, mate_id))
            registerJob1 = self.addRegisterAndMoveSplitFileJob(
                executable=self.RegisterAndMoveSplitSequenceFiles, \
                inputDir=splitOutputDir, outputDir=sequenceAbsDir, 
                relativeOutputDir=individual_sequence.path, logFile=logFile,\
                individual_sequence_id=individual_sequence.id, 
                original_file=bamInputF, 
                library=library, mate_id=mate_id, \
                parentJobLs=[splitReadFileJob2, createSequenceAbsDirJob], 
                job_max_memory=100, walltime = 60, \
                commit=commit, sequence_format=sequence_format, extraDependentInputLs=None, \
                transferOutput=True, sshDBTunnel=self.needSSHDBTunnel)
            
            """
            jobFname = os.path.join(self.jobFileDir, 'job%s.bam2fastq.sh'%(code))
            self.writeQsubJob(jobFname, bamFname, os.path.join(self.data_dir, individual_sequence.path),
                self.vervet_path)
            commandline = 'qsub %s'%(jobFname)
            if self.commit:	#qsub only when db transaction will be committed.
                return_data = runLocalCommand(commandline, report_stderr=True, report_stdout=True)
            """
        sys.stderr.write("%s jobs.\n"%(self.no_of_jobs))

    def saveTCGATissueSourceSiteIntoDB(self, db_main=None, input_path=None):
        """
        20170412
        """
        sys.stderr.write("Parsing TCGA tissue source sites from %s... \n"%input_path)
        reader = MatrixFile(path=input_path, delimiter="\t")
        header = reader.next()
        tssCode2dbEntry = {}
        for row in reader:
            tss_code, source_site_name, study_name, bcr_code = row[:4]
            if tss_code[0]=='0':	#remove the front 0
                tss_code = tss_code[1:]
            study = db_main.getStudy(short_name=study_name)
            site = db_main.getSite(short_name=tss_code, description=source_site_name, \
                region=bcr_code, study_id=study.id)
            tssCode2dbEntry[tss_code] = site
        sys.stderr.write("%s sites.\n"% len(tssCode2dbEntry))
        
        return tssCode2dbEntry

    def getSamplesFromInputDir_TCGA(self, db_main=None, inputDir=None, tax_id=None):
        """
        20170412 Get most info about bams from the .xml file.
        COAD_TCGA/6795b83b-1ca7-4060-a078-6649fa9004bb
            35f008d5-1a15-463e-88e6-955583a33aa6_analysis.xml
            TCGA-AA-3553-10A-01D-1167-02_IlluminaHiSeq-DNASeq_whole.bam
            ...
        """
        print("Getting TCGA samples from %s... "%inputDir, flush=True)
        counter = 0
        noOfSamplesIntoDB = 0
        noOfSamplesWithoutCenter = 0
        
        sample_obj_ls = []
        import xml.etree.ElementTree as ET
        tcga_barcode_re = re.compile(r'TCGA-(?P<tss>[0-9A-Z]{2})-(?P<participant>[0-9A-Z]{4})-'
            r'(?P<tissue>[0-9]{2})[A-Z]-(?P<portion>[0-9]{2}[A-Z])-(?P<plate>[0-9A-Z]{4})(-(?P<center>[0-9]{2})|)')
        for filename in os.listdir(inputDir):
            input_path = os.path.join(inputDir, filename)
            counter += 1
            if os.path.isdir(input_path):
                noOfTargets = 0
                sample_obj = PassingData(uuid=filename)
                for subfilename in os.listdir(input_path):
                    fname_prefix, fname_suffix = os.path.splitext(subfilename)
                    if fname_suffix=='.bam':
                        bamPath = os.path.join(input_path, subfilename)
                        #tcga_barcode = fname_prefix.split('_')[0], difficult to parse
                        # f316e3dff60e757701731d0c0cd94c3d.bam
                        # G92908.TCGA-A6-3809-01A-01D-A46W-08.2.bam
                        # TCGA-AA-3672-01A-01D-0957-02_IlluminaGAII-DNASeq_whole.bam
                        # TCGA-AO-A0JF-01A-11D-A060_130719_SN1120_0270_AC2CVRACXX_s_1_rg.sorted.bam
                        sample_obj.bamPath = bamPath
                        sample_obj.baiPath = "%s.bai"%bamPath
                        noOfTargets+=1
                    elif fname_suffix=='.xml' and fname_prefix[-9:]=="_analysis":
                        #parse to get md5sum
                        tree = ET.parse(os.path.join(input_path, subfilename))
                        root = tree.getroot()
                        description = root.findall("ANALYSIS_SET/ANALYSIS/DESCRIPTION")[0].text
                        #get participant code, tss_id (site.short_name), sample_id (tissue_id),
                        #  center_id (sequencer.id)
                        result = tcga_barcode_re.search(description)
                        if result is None:
                            print("\n  ERROR: not enough (!=7) fields in tcga_barcode of %s/%s: %s."%
                                (input_path, subfilename, description), flush=True)
                            sys.exit(-2)
                        tcga_barcode = result.group()
                        split_row = tcga_barcode.split('-')
                        
                        tss_code = result.group('tss')
                        site = db_main.getSite(short_name=tss_code)
                        
                        participant_code = '-'.join(split_row[:3])
                        tissue_id = int(result.group('tissue'))
                        
                        if result.group('center'):	#some don't have it
                            sequencer_id = int(result.group('center'))
                        else:
                            sequencer_id = None
                            noOfSamplesWithoutCenter += 1
                        db_entry = db_main.getIndividual(code=participant_code, 
                            name=None, sex=None, age=None, \
                            site=site, tax_id=tax_id, study=site.study)
                        sample_obj.db_entry = db_entry
                        sample_obj.tissue_id = tissue_id
                        sample_obj.sequencer_id = sequencer_id
                        sample_obj.tcga_barcode = tcga_barcode
                        
                        #get md5sum
                        fileElemLs = root.findall("ANALYSIS_SET/ANALYSIS/DATA_BLOCK/FILES/FILE")
                        if len(fileElemLs)!=1:
                            print(f"\n  ERROR: number of files with MD5SUM in %s/%s is not one. It is %s."%
                                (input_path, subfilename, len(fileElemLs)), flush=True)
                            sys.exit(-2)
                        filename = fileElemLs[0].get('filename')
                        md5sum = fileElemLs[0].get('checksum')
                        if filename is None or md5sum is None:
                            import pdb
                            pdb.set_trace()
                        if filename[-4:]==".bam":
                            sample_obj.md5sum = md5sum
                        else:
                            print("\n  ERROR: md5sum is for a non-bam file, %s."%(filename),
                                flush=True)
                            sys.exit(-3)
                        noOfTargets+=1
                if noOfTargets==2:
                    sample_obj_ls.append(sample_obj)
                    noOfSamplesIntoDB += 1
        print("%s samples in this folder %s. %s samples into DB. %s samples without center."%(
            counter, inputDir, len(sample_obj_ls), noOfSamplesWithoutCenter), flush=True)
        return sample_obj_ls


    def addJobsToImportTCGABam(self, db_main=None, input_path=None, \
        tax_id=9606, data_dir=None, \
        minNoOfReads=None, commit=None,\
        sequencer_name=None, sequence_type_name=None, sequence_format=None, **keywords):
        """
        20170407
        sequence_type_name is not used in this function. it's determined from tcga barcode.
        input:
            folder of bams
            tissueSourceSite.tsv file
        """
        #parse the tissueSourceSite.tsv file to get all tissue source sites info,
        #  save them into table site, and study_id from db
        self.saveTCGATissueSourceSiteIntoDB(db_main=db_main, 
            input_path=self.tissueSourceSiteFname)
        
        #find all bam in a folder,
        #parse barcode to get participant code, tss_id (site.short_name), 
        # sample_id (tissue_id), center_id (sequencer.id),
        #folder name is uuid
        #uuid is for sequence only, add as isq.comment
        #parse analysis.xml to get md5sum 
        #add all individual, tissue, site, study, sequencer, seq_center into db
        sample_obj_ls = self.getSamplesFromInputDir_TCGA(db_main=db_main, \
            inputDir=input_path, tax_id=tax_id)
        sys.stderr.write("%s total TCGA samples.\n"%(len(sample_obj_ls)))
        
        sam2fastqOutputDir = 'sam2fastq'
        sam2fastqOutputDirJob = self.addMkDirJob(outputDir=sam2fastqOutputDir)
        
        for sample_obj in sample_obj_ls:
            bamInputF = self.registerOneInputFile(sample_obj.bamPath)
            bamFileSize = utils.getFileOrFolderSize(sample_obj.bamPath)
            baiInputF = self.registerOneInputFile(sample_obj.baiPath)
            bamBaseFname = os.path.split(sample_obj.bamPath)[1]
            bamBaseFnamePrefix = os.path.splitext(bamBaseFname)[0]
            library = sample_obj.tcga_barcode
            
            #add VerifyFileMD5Sum, cpu=8 to reduce its IO load
            verifyMD5SumJob = self.addGenericJob(executable=self.VerifyFileMD5Sum, \
                inputFile=bamInputF, \
                inputArgumentOption="-i", \
                outputFile=None, outputArgumentOption="-o", \
                inputFileList=None, argumentForEachFileInInputFileList=None, \
                parentJob=None, parentJobLs=None, extraDependentInputLs=None, \
                extraOutputLs=[bamInputF], transferOutput=False, \
                frontArgumentList=None, \
                extraArguments="--correct_md5sum %s"%sample_obj.md5sum, \
                extraArgumentList=None, \
                job_max_memory=None, sshDBTunnel=None, \
                key2ObjectForJob=None, objectWithDBArguments=None, no_of_cpus=8, walltime=30)
            
            #no_of_cpus=4 to reduce its IO
            #assume a 10GB file needing a 30GB memory with a cap of 150G
            memoryNeeded = min(max(60000, int(bamFileSize/10000000000.0*30000)), 150000)
            outputFnamePrefix = os.path.join(sam2fastqOutputDir, '%s'%(library))
            convertBamToFastqAndGzip_job = self.addConvertBamToFastqAndGzipJob(
                executable=self.SamToFastqJava, \
                inputF=bamInputF, outputFnamePrefix=outputFnamePrefix, \
                parentJobLs=[sam2fastqOutputDirJob, verifyMD5SumJob], \
                job_max_memory=memoryNeeded, no_of_cpus=4, walltime = 800, \
                extraDependentInputLs=[baiInputF], \
                transferOutput=False)
            
            #job to check if each file is empty or not. If empty, exit non-0.
            #  If not, add db isq entry, output isq-id, exit 0.
            #uuid is for sequence only, add as isq.comment
            unpaired_seq_db_idFile = File("%s_unpaired_seq_db_id.txt"%outputFnamePrefix)
            registerUnpairedSequence2DBJob =self.addRegisterIndividualSequence2DBJob(
                inputFile=convertBamToFastqAndGzip_job.output_unpaired, \
                outputFile=unpaired_seq_db_idFile, \
                individual_id=sample_obj.db_entry.id, sequencer_id=sample_obj.sequencer_id,\
                sequence_type_name="SingleRead", sequence_format=sequence_format, \
                copy_original_file=0, tissue_name=None, tissue_id=sample_obj.tissue_id, \
                coverage=None, quality_score_format="Standard", filtered=0,\
                parent_individual_sequence_id=None,\
                read_count=None, no_of_chromosomes=None, \
                sequence_batch_id=None, version=None, \
                is_contaminated=0, outdated_index=0, comment=sample_obj.uuid,\
                data_dir=data_dir,\
                original_sequence_filepath=sample_obj.bamPath, \
                original_sequence_library=library, \
                original_sequence_mate_id=None, \
                original_sequence_md5sum=sample_obj.md5sum, \
                parentJobLs=[convertBamToFastqAndGzip_job], \
                job_max_memory=100, walltime = 60, commit=self.commit, \
                extraDependentInputLs=[bamInputF], extraArguments=None, \
                transferOutput=False, sshDBTunnel=self.needSSHDBTunnel)
            
            splitSROutputDir = 'split_%s_singleRead'%(library)
            splitSROutputDirJob = self.addMkDirJob(outputDir=splitSROutputDir)
            
            splitFastQFnamePrefix = os.path.join(splitSROutputDir, '%s_singleRead'%(library))
            logFile = File('%s_split.log'%(splitFastQFnamePrefix))
            splitSRFileJob = self.addSplitReadFileJob(
                inputF=convertBamToFastqAndGzip_job.output_unpaired, \
                outputFnamePrefix=splitFastQFnamePrefix, \
                outputFnamePrefixTail="", minNoOfReads=minNoOfReads, \
                logFile=logFile, \
                parentJobLs=[registerUnpairedSequence2DBJob, splitSROutputDirJob], \
                job_max_memory=4000, walltime = 800, no_of_cpus=4, \
                extraDependentInputLs=None, transferOutput=True)
            
            #mate_id is set to null for singleRead sequences
            logFile = File('%s_register.log'%(splitFastQFnamePrefix))
            registerSRJob = self.addRegisterAndMoveSplitFileJob(
                inputFile=registerUnpairedSequence2DBJob.output,\
                inputDir=splitSROutputDir, logFile=logFile,\
                library=library, mate_id=None, \
                parentJobLs=[splitSRFileJob], job_max_memory=100, walltime = 60, \
                commit=commit, sequence_format=sequence_format, 
                extraDependentInputLs=None, \
                transferOutput=True, sshDBTunnel=self.needSSHDBTunnel)
            
            #job to check if PE sequence file is empty or not. but only check one mate file.
            pairedEnd_seq_db_idFile = File("%s_pairedEnd_seq_db_id.txt"%(outputFnamePrefix))
            registerPairedEndSequence2DBJob =self.addRegisterIndividualSequence2DBJob(
                inputFile=convertBamToFastqAndGzip_job.output1, \
                outputFile=pairedEnd_seq_db_idFile, \
                individual_id=sample_obj.db_entry.id, sequencer_id=sample_obj.sequencer_id,\
                sequence_type_name="PairedEnd", sequence_format=sequence_format, \
                copy_original_file=0, tissue_name=None, tissue_id=sample_obj.tissue_id, \
                coverage=None, quality_score_format="Standard", filtered=0,\
                parent_individual_sequence_id=None,\
                read_count=None, no_of_chromosomes=None, \
                sequence_batch_id=None, version=None, \
                is_contaminated=0, outdated_index=0, comment=sample_obj.uuid,\
                data_dir=data_dir,\
                original_sequence_filepath=sample_obj.bamPath, \
                original_sequence_library=library, \
                original_sequence_mate_id=None, \
                original_sequence_md5sum=sample_obj.md5sum, \
                parentJobLs=[convertBamToFastqAndGzip_job], \
                job_max_memory=100, walltime = 60, commit=self.commit, \
                extraDependentInputLs=[bamInputF], extraArguments=None, \
                transferOutput=False, sshDBTunnel=self.needSSHDBTunnel)
            
            splitOutputDir = 'split_%s_PE'%(library)
            #same directory containing split files from both mates is fine
            #  as RegisterAndMoveSplitSequenceFiles could pick up.
            splitOutputDirJob = self.addMkDirJob(outputDir=splitOutputDir)
            
            mate_id = 1
            splitFastQFnamePrefix = os.path.join(splitOutputDir, '%s_%s'%(library, mate_id))
            logFile = File('%s_split.log'%(splitFastQFnamePrefix))
            splitReadFileJob1 = self.addSplitReadFileJob(
                inputF=convertBamToFastqAndGzip_job.output1, 
                outputFnamePrefix=splitFastQFnamePrefix, \
                outputFnamePrefixTail="", minNoOfReads=minNoOfReads, \
                logFile=logFile, \
                parentJobLs=[registerPairedEndSequence2DBJob, splitOutputDirJob], \
                job_max_memory=4000, walltime = 800, no_of_cpus=4, \
                extraDependentInputLs=None, transferOutput=True)
            
            logFile = File('%s_register.log'%(splitFastQFnamePrefix))
            registerJob1 = self.addRegisterAndMoveSplitFileJob(
                inputFile=registerPairedEndSequence2DBJob.output,\
                inputDir=splitOutputDir, logFile=logFile,\
                library=library, mate_id=mate_id, \
                parentJobLs=[splitReadFileJob1], job_max_memory=100, walltime = 60, \
                commit=commit, sequence_format=sequence_format, extraDependentInputLs=None, \
                transferOutput=True, sshDBTunnel=self.needSSHDBTunnel)
            #handle the 2nd end
            mate_id = 2
            splitFastQFnamePrefix = os.path.join(splitOutputDir, '%s_%s'%(library, mate_id))
            logFile = File('%s_split.log'%(splitFastQFnamePrefix))
            splitReadFileJob2 = self.addSplitReadFileJob(
                inputF=convertBamToFastqAndGzip_job.output2, 
                outputFnamePrefix=splitFastQFnamePrefix, \
                outputFnamePrefixTail="", minNoOfReads=minNoOfReads, \
                logFile=logFile, \
                parentJobLs=[registerPairedEndSequence2DBJob, splitOutputDirJob], \
                job_max_memory=4000, walltime = 800, no_of_cpus=4, \
                extraDependentInputLs=None, transferOutput=True)
            
            logFile = File('%s_register.log'%(splitFastQFnamePrefix))
            registerJob2 = self.addRegisterAndMoveSplitFileJob(
                inputFile=registerPairedEndSequence2DBJob.output,\
                inputDir=splitOutputDir, logFile=logFile,\
                library=library, mate_id=mate_id, \
                parentJobLs=[splitReadFileJob2], job_max_memory=100, walltime = 60, \
                commit=commit, sequence_format=sequence_format, extraDependentInputLs=None, \
                transferOutput=True, sshDBTunnel=self.needSSHDBTunnel)
        sys.stderr.write("%s jobs.\n"%(self.no_of_jobs))

    def getSamplesFromInputDir_HCC1187(self, db_main=None, inputDir=None, tax_id=None):
        """
        20170607
        inputDir is /y/home/luozhihui/Downloads/
            114G /y/home/luozhihui/Downloads/HCC1187BL_S1.bam
            8.7M /y/home/luozhihui/Downloads/HCC1187BL_S1.bam.bai
            237G /y/home/luozhihui/Downloads/HCC1187C_S1.bam
            8.7M /y/home/luozhihui/Downloads/HCC1187C_S1.bam.bai
        """
        print(f"Getting HCC1187 samples from {inputDir} ... ", flush=True)
        counter = 0
        noOfSamplesIntoDB = 0
        sample_obj_ls = []
        for filename in os.listdir(inputDir):
            input_path = os.path.join(inputDir, filename)
            counter += 1
            if os.path.isfile(input_path):
                fname_prefix, fname_suffix = os.path.splitext(filename)
                if fname_suffix=='.bam':
                    sample_obj = PassingData()
                    code = fname_prefix.split("_")[0]
                    bamPath = input_path
                    sample_obj.bamPath = bamPath
                    sample_obj.baiPath = "%s.bai"%bamPath
                    if code=="HCC1187C":
                        tissue_id = 50	#10 (Blood normal) or 50 (cancer cell line)
                    elif code=="HCC1187BL":
                        tissue_id = 10
                    else:
                        print(f" ERROR: unexpected code {code} in getting tissue_id.",
                            flush=True)
                        import pdb
                        pdb.set_trace()
                        
                    sequencer_id = 10000
                    db_entry = db_main.getIndividual(code=code, name=None, sex=None, age=None, \
                                site_id=self.site_id, tax_id=tax_id, study_id=self.study_id)
                    sample_obj.db_entry = db_entry
                    sample_obj.tissue_id = tissue_id
                    sample_obj.sequencer_id = sequencer_id
                    sample_obj_ls.append(sample_obj)
                    noOfSamplesIntoDB += 1
        print(f"{counter} samples in this folder {inputDir}. "
            f"{len(sample_obj_ls)} samples into DB.",
            flush=True)
        return sample_obj_ls

    def addJobsToImportHCC1187Bam(self, db_main=None, input_path=None, \
        tax_id=9606, data_dir=None, \
        minNoOfReads=None, commit=None,\
        sequencer_name=None, sequence_type_name=None, sequence_format=None, **keywords):
        """
        20170607
        input: A folder of bam files.
        
        - Find all bam in the folder.
        - parse barcode to get participant code, tss_id (site.short_name), 
            sample_id (tissue_id), center_id (sequencer.id),
        - folder name is uuid.
        - uuid is for sequence only, add as isq.comment.
        - parse analysis.xml to get md5sum.
        - add all individual, tissue, site, study, sequencer, seq_center into db.
        """
        sample_obj_ls = self.getSamplesFromInputDir_HCC1187(db_main=db_main, 
            inputDir=input_path, tax_id=tax_id)
        print(f"{len(sample_obj_ls)} TCGA samples in total.", flush=True)
        
        sam2fastqOutputDir = 'sam2fastq'
        sam2fastqOutputDirJob = self.addMkDirJob(outputDir=sam2fastqOutputDir)
        
        for sample_obj in sample_obj_ls:
            bamInputF = self.registerOneInputFile(sample_obj.bamPath)
            bamFileSize = utils.getFileOrFolderSize(sample_obj.bamPath)
            baiInputF = self.registerOneInputFile(sample_obj.baiPath)
            bamBaseFname = os.path.split(sample_obj.bamPath)[1]
            bamBaseFnamePrefix = os.path.splitext(bamBaseFname)[0]
            library = sample_obj.db_entry.code
            
            #no_of_cpus=4 to reduce its IO
            #assume a 10GB file needing a 30GB memory with a cap of 150G
            memoryNeeded = min(max(60000, int(bamFileSize/10000000000.0*30000)), 150000)
            outputFnamePrefix = os.path.join(sam2fastqOutputDir, '%s'%(library))
            convertBamToFastqAndGzip_job = self.addConvertBamToFastqAndGzipJob(
                executable=self.SamToFastqJava, \
                inputF=bamInputF, outputFnamePrefix=outputFnamePrefix, \
                parentJobLs=[sam2fastqOutputDirJob], \
                job_max_memory=memoryNeeded, no_of_cpus=4, walltime = 800, \
                extraDependentInputLs=[baiInputF], \
                transferOutput=False)
            
            #job to check if each file is empty or not. If empty, exit non-0.
            #  If not, add db isq entry, output isq-id, exit 0.
            #uuid is for sequence only, add as isq.comment
            unpaired_seq_db_idFile = File("%s_unpaired_seq_db_id.txt"%outputFnamePrefix)
            registerUnpairedSequence2DBJob =self.addRegisterIndividualSequence2DBJob(
                inputFile=convertBamToFastqAndGzip_job.output_unpaired, \
                outputFile=unpaired_seq_db_idFile, \
                individual_id=sample_obj.db_entry.id, sequencer_id=sample_obj.sequencer_id,\
                sequence_type_name="SingleRead", sequence_format=sequence_format, \
                copy_original_file=0, tissue_name=None, tissue_id=sample_obj.tissue_id, \
                coverage=None, quality_score_format="Standard", filtered=0,\
                parent_individual_sequence_id=None,\
                read_count=None, no_of_chromosomes=None, \
                sequence_batch_id=None, version=None, \
                is_contaminated=0, outdated_index=0, comment=None,\
                data_dir=data_dir,\
                original_sequence_filepath=sample_obj.bamPath, \
                original_sequence_library=library, \
                original_sequence_mate_id=None, \
                original_sequence_md5sum=None, \
                parentJobLs=[convertBamToFastqAndGzip_job], \
                job_max_memory=100, walltime = 60, commit=self.commit, \
                extraDependentInputLs=[bamInputF], extraArguments=None, \
                transferOutput=False, sshDBTunnel=self.needSSHDBTunnel)
            
            splitSROutputDir = 'split_%s_singleRead'%(library)
            splitSROutputDirJob = self.addMkDirJob(outputDir=splitSROutputDir)
            
            splitFastQFnamePrefix = os.path.join(splitSROutputDir, '%s_singleRead'%(library))
            logFile = File('%s_split.log'%(splitFastQFnamePrefix))
            splitSRFileJob = self.addSplitReadFileJob(
                inputF=convertBamToFastqAndGzip_job.output_unpaired, \
                outputFnamePrefix=splitFastQFnamePrefix, \
                outputFnamePrefixTail="", minNoOfReads=minNoOfReads, \
                logFile=logFile, \
                parentJobLs=[registerUnpairedSequence2DBJob, splitSROutputDirJob], \
                job_max_memory=4000, walltime = 800, no_of_cpus=4, \
                extraDependentInputLs=None, transferOutput=True)
            
            #mate_id is set to null for singleRead sequences
            logFile = File('%s_register.log'%(splitFastQFnamePrefix))
            registerSRJob = self.addRegisterAndMoveSplitFileJob(
                inputFile=registerUnpairedSequence2DBJob.output,\
                inputDir=splitSROutputDir, logFile=logFile,\
                library=library, mate_id=None, \
                parentJobLs=[splitSRFileJob], job_max_memory=100, walltime = 60, \
                commit=commit, sequence_format=sequence_format, extraDependentInputLs=None, \
                transferOutput=True, sshDBTunnel=self.needSSHDBTunnel)
            
            #job to check if PE sequence file is empty or not. but only check one mate file.
            pairedEnd_seq_db_idFile = File("%s_pairedEnd_seq_db_id.txt"%(outputFnamePrefix))
            registerPairedEndSequence2DBJob =self.addRegisterIndividualSequence2DBJob(
                inputFile=convertBamToFastqAndGzip_job.output1, \
                outputFile=pairedEnd_seq_db_idFile, \
                individual_id=sample_obj.db_entry.id, sequencer_id=sample_obj.sequencer_id,\
                sequence_type_name="PairedEnd", sequence_format=sequence_format, \
                copy_original_file=0, tissue_name=None, tissue_id=sample_obj.tissue_id, \
                coverage=None, quality_score_format="Standard", filtered=0,\
                parent_individual_sequence_id=None,\
                read_count=None, no_of_chromosomes=None, \
                sequence_batch_id=None, version=None, \
                is_contaminated=0, outdated_index=0, comment=None,\
                data_dir=data_dir,\
                original_sequence_filepath=sample_obj.bamPath, \
                original_sequence_library=library, \
                original_sequence_mate_id=None, \
                original_sequence_md5sum=None, \
                parentJobLs=[convertBamToFastqAndGzip_job], \
                job_max_memory=100, walltime = 60, commit=self.commit, \
                extraDependentInputLs=[bamInputF], extraArguments=None, \
                transferOutput=False, sshDBTunnel=self.needSSHDBTunnel)
    
            splitOutputDir = 'split_%s_PE'%(library)
            #same directory containing split files from both mates is fine as
            #  RegisterAndMoveSplitSequenceFiles could pick up.
            splitOutputDirJob = self.addMkDirJob(outputDir=splitOutputDir)
            
            mate_id = 1
            splitFastQFnamePrefix = os.path.join(splitOutputDir, '%s_%s'%(library, mate_id))
            logFile = File('%s_split.log'%(splitFastQFnamePrefix))
            splitReadFileJob1 = self.addSplitReadFileJob(
                inputF=convertBamToFastqAndGzip_job.output1, outputFnamePrefix=splitFastQFnamePrefix, \
                outputFnamePrefixTail="", minNoOfReads=minNoOfReads, \
                logFile=logFile, \
                parentJobLs=[registerPairedEndSequence2DBJob, splitOutputDirJob], \
                job_max_memory=4000, walltime = 800, no_of_cpus=4, \
                extraDependentInputLs=None, transferOutput=True)
            
            logFile = File('%s_register.log'%(splitFastQFnamePrefix))
            registerJob1 = self.addRegisterAndMoveSplitFileJob(
                inputFile=registerPairedEndSequence2DBJob.output,\
                inputDir=splitOutputDir, logFile=logFile,\
                library=library, mate_id=mate_id, \
                parentJobLs=[splitReadFileJob1], job_max_memory=100, walltime = 60, \
                commit=commit, sequence_format=sequence_format, \
                extraDependentInputLs=None, \
                transferOutput=True, sshDBTunnel=self.needSSHDBTunnel)
            #handle the 2nd end
            mate_id = 2
            splitFastQFnamePrefix = os.path.join(splitOutputDir, '%s_%s'%(library, mate_id))
            logFile = File('%s_split.log'%(splitFastQFnamePrefix))
            splitReadFileJob2 = self.addSplitReadFileJob(
                inputF=convertBamToFastqAndGzip_job.output2, \
                outputFnamePrefix=splitFastQFnamePrefix, \
                outputFnamePrefixTail="", minNoOfReads=minNoOfReads, \
                logFile=logFile, \
                parentJobLs=[registerPairedEndSequence2DBJob, splitOutputDirJob], \
                job_max_memory=4000, walltime = 800, no_of_cpus=4, \
                extraDependentInputLs=None, transferOutput=True)
            
            logFile = File('%s_register.log'%(splitFastQFnamePrefix))
            registerJob2 = self.addRegisterAndMoveSplitFileJob(
                inputFile=registerPairedEndSequence2DBJob.output,\
                inputDir=splitOutputDir, logFile=logFile,\
                library=library, mate_id=mate_id, \
                parentJobLs=[splitReadFileJob2], job_max_memory=100, walltime = 60, \
                commit=commit, sequence_format=sequence_format, \
                extraDependentInputLs=None, \
                transferOutput=True, sshDBTunnel=self.needSSHDBTunnel)
            
        print(f"{self.no_of_jobs} jobs.", flush=True)

    def run(self):
        """
        2011-8-3
        """
        pdata = self.setup_run()
        if self.inputType==1:
            self.addJobsToImportTCGABam(db_main=self.db_main, \
            sample_sheet=self.sample_sheet, \
            input_path=self.input_path, data_dir=self.data_dir, \
            minNoOfReads=self.minNoOfReads, commit=self.commit,\
            sequencer_name=self.sequencer_name, 
            sequence_type_name=self.sequence_type_name, \
            sequence_format=self.sequence_format)
        elif self.inputType==2:
            self.addJobsToImportHCC1187Bam(db_main=self.db_main, \
            sample_sheet=self.sample_sheet, \
            input_path=self.input_path, data_dir=self.data_dir, \
            minNoOfReads=self.minNoOfReads, commit=self.commit,\
            sequencer_name=self.sequencer_name, 
            sequence_type_name=self.sequence_type_name, \
            sequence_format=self.sequence_format)
        # Write the DAX
        self.end_run()
        
if __name__ == '__main__':
    main_class = ImportIndividualSequence2DB
    po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
    instance = main_class(**po.long_option2value)
    instance.run()
