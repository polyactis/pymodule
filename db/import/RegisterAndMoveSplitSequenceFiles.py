#!/usr/bin/env python3
import copy
import logging
import re
import os
import sys
from palos import ProcessOptions, utils
from palos.mapper.AbstractSunsetMapper import AbstractSunsetMapper

ParentClass = AbstractSunsetMapper
class RegisterAndMoveSplitSequenceFiles(ParentClass):
    """
Examples:
    %s --drivername postgresql
        --hostname pdc --dbname pmdb --schema vervet --db_user yh 
        -i isq_id.txt --inputDir 3185_gerald_D14CGACXX_7_GCCAAT
        -o ~/vervet/db/individual_sequence/3185_6059_2002099_GA_0_0
        --relativeOutputDir individual_sequence/3185_6059_2002099_GA_0_0
        --library gerald_D14CGACXX_7_GCCAAT --individual_sequence_id 3185
        --sequence_format fastq --commit
        --mate_id 1 --bamFilePath gerald_D14CGACXX_7_GCCAAT.bam 
        --logFilename  3185_gerald_D14CGACXX_7_GCCAAT_1.register.log

Description:
    Program to register and move/cp split files by picard's
        SplitReadFile.jar to db storage.
    If inputFname is given, it will be parsed to get:
        outputDir, relativeOutputDir, individual_sequence_id, 
        individual_sequence_file_raw_id.
        One variable on each line.
    """
    __doc__ = __doc__%(sys.argv[0])
    parse_args_from_commandline = 0
    option_default_dict = copy.deepcopy(ParentClass.option_default_dict)
    option_default_dict.pop(('outputFname', 0, ))
    option_default_dict.pop(('outputFnamePrefix', 0, ))
    option_default_dict.update({
        ('inputDir', 1, ): ['', '', 1, 
            'input folder that contains split fastq files', ],
        ('original_file_path', 0, ): ['', '', 1, 
            'The original file from which all input files originate '],
        ('outputDir', 0, ): ['', 'o', 1, 
            'output folder to which files from inputDir will be moved.', ],
        ('relativeOutputDir', 0, ): ['', '', 1, 
            'The output folder path relative to db.data_dir. '
            'It should form the last part of outputDir. '
            'If not given, parse from inputFname.', ],
        ('individual_sequence_id', 0, int): [None, '', 1, 
            'individual_sequence.id associated with all files.'],
        ('individual_sequence_file_raw_id', 0, int): [None, '', 1, 
            'individual_sequence_file_raw.id associated with the raw file'],
        ('library', 0, ): [None, 'l', 1, 
            'library name for files in inputDir', ],
        ('parse_args_from_commandline', 0, int): [0, '', 1, 
            'If 1, parse args from the commandline, not from the input file.', ],
        ('mate_id', 0, int): [None, 'm', 1, 
            '1: first end; 2: 2nd end. of paired-end or mate-paired libraries'],
        ("sequence_format", 1, ): ["fastq", 'f', 1, 'fasta, fastq, etc.'],
        ('logFilename', 0, ): [None, 'g', 1, 
            'output file to contain logs. optional.'],
        })

    def __init__(self, **keywords):
        """
        """
        #connectDB(), and setup srcFilenameLs and dstFilenameLs
        ParentClass.__init__(self, inputFnameLs=None, **keywords)
        if self.inputFname and self.parse_args_from_commandline==0:
            self.parseArgumentsFromFile(self.inputFname)
    
    def parseArgumentsFromFile(self, inputFname):
        """
        20190206
        """
        #parse inputFname to get individual_sequence_id & 
        # individual_sequence_file_raw_id and others.
        inputFile = utils.openGzipFile(inputFname)
        input_variable_dict = {}
        for line in inputFile:
            var_name, var_value = line.strip().split(": ")
            input_variable_dict[var_name] = var_value
        inputFile.close()
        
        individual_sequence_id = input_variable_dict.get(
            "individual_sequence_id", self.individual_sequence_id)
        if individual_sequence_id:
            individual_sequence_id = int(individual_sequence_id)
            self.individual_sequence_id = individual_sequence_id

        individual_sequence_file_raw_id = input_variable_dict.get(
            "individual_sequence_file_raw_id", 
            self.individual_sequence_file_raw_id)
        if individual_sequence_file_raw_id:
            individual_sequence_file_raw_id = \
                int(individual_sequence_file_raw_id)
            self.individual_sequence_file_raw_id = \
                individual_sequence_file_raw_id

        self.outputDir = input_variable_dict.get("outputDir", self.outputDir)
        self.relativeOutputDir = input_variable_dict.get("relativeOutputDir",
            self.relativeOutputDir)
        
        relativePathIndex = self.outputDir.find(self.relativeOutputDir)
        noOfCharsInRelativeOutputDir = len(self.relativeOutputDir)
        if self.outputDir[relativePathIndex:relativePathIndex+\
                noOfCharsInRelativeOutputDir]!=self.relativeOutputDir:
            logging.error(f'relativeOutputDir {self.relativeOutputDir} is not'
            f' the last part of outputDir {self.outputDir}.')
            sys.exit(4)
    
    def parseSplitOrderOutOfFilename(self, filename, mate_id=None):
        """
        Parse the split order out of filename. For example:
        For file gerald_81LL0ABXX_4_TTAGGC_2_1.fastq.gz:
            .*_(mate_id)_(split_order).fastq.gz
            mate_id is 2.
            split_order is 1.
        """
        if mate_id:
            middle = mate_id
        else:
            middle = 'singleRead'
        split_order_pattern = re.compile(
            r'(?P<other>.*)_%s_(?P<split_order>\d+).fastq'%(
                middle
            ))
        split_order_search_result = split_order_pattern.search(filename)
        if split_order_search_result:
            split_order = split_order_search_result.group('split_order')
        else:
            split_order = None
        return split_order
    
    def run(self):
        """
        """
        if self.debug:
            import pdb
            pdb.set_trace()
        db_main = self.db_main
        session = db_main.session
        session.begin()
        if self.original_file_path:
            file_raw_db_entry = self.db_main.registerOriginalSequenceFileToDB(
                self.original_file_path, library=self.library,
                individual_sequence_id=self.individual_sequence_id,
                mate_id=self.mate_id, 
                md5sum=None)
            self.individual_sequence_file_raw_id = file_raw_db_entry.id
        
        counter = 0
        real_counter = 0
        for filename in os.listdir(self.inputDir):
            split_order = self.parseSplitOrderOutOfFilename(filename,
                self.mate_id)
            counter += 1
            if split_order:
                #save db entry
                db_entry = db_main.getIndividualSequenceFile(
                    self.individual_sequence_id, library=self.library,
                    mate_id=self.mate_id, \
                    split_order=int(split_order), format=self.sequence_format,
                    filtered=0, parent_individual_sequence_file_id=None, \
                    individual_sequence_file_raw_id=\
                        self.individual_sequence_file_raw_id)
                
                #move the file
                exitCode = db_main.moveFileIntoDBAffiliatedStorage(
                    db_entry=db_entry, filename=filename,
                    inputDir=self.inputDir, outputDir=self.outputDir,
                    relativeOutputDir=self.relativeOutputDir,
                    shellCommand='cp -rL',
                    srcFilenameLs=self.srcFilenameLs,
                    dstFilenameLs=self.dstFilenameLs,
                    constructRelativePathFunction=None)
                if exitCode!=0:
                    logging.error(f"moveFileIntoDBAffiliatedStorage() exits "
                        f"with {exitCode} code.")
                    self.sessionRollback(session)
                    #delete all recorded target files
                    self.cleanUpAndExitOnFailure(exitCode=exitCode)
                real_counter += 1
        
        if self.logFilename:
            outf = open(self.logFilename, 'w')
            outf.write("%s files processed, %s of them added into db.\n"%(\
                counter, real_counter))
            outf.close()
        
        if self.commit:
            try:
                session.commit()
                #delete all source files
                self.rmGivenFiles(filenameLs=self.srcFilenameLs)
            except:
                sys.stderr.write('Except type: %s\n'%repr(sys.exc_info()))
                import traceback
                traceback.print_exc()
                #delete all target files.
                self.cleanUpAndExitOnFailure(exitCode=3)
        else:
            self.sessionRollback(session)
            #delete all target files
            self.cleanUpAndExitOnFailure(exitCode=0)
        
if __name__ == '__main__':
    main_class = RegisterAndMoveSplitSequenceFiles
    po = ProcessOptions(sys.argv, main_class.option_default_dict,\
        error_doc=main_class.__doc__)
    instance = main_class(**po.long_option2value)
    instance.run()
