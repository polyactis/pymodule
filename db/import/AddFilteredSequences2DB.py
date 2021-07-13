#!/usr/bin/env python3
"""
Examples:
    %s -n 130 -u yh 1.fastq.gz 2.fastq.gz
    
    %s 

Description:
    2011-8-18
        follow-up to FilterShortRead.py or picard/dist/FilterRead.jar to add filtered sequences to db
"""

import sys, os, math
__doc__ = __doc__ % (sys.argv[0], sys.argv[0])

from palos import ProcessOptions
from palos.db import SunsetDB
from palos.mapper.AbstractDBJob import AbstractDBJob
from RegisterAndMoveSplitSequenceFiles import RegisterAndMoveSplitSequenceFiles


class AddFilteredSequences2DB(RegisterAndMoveSplitSequenceFiles):
    __doc__ = __doc__
    option_default_dict = AbstractDBJob.option_default_dict.copy()
    option_default_dict.pop(('outputFname', 0, ))
    option_default_dict.pop(('outputFnamePrefix', 0, ))
    option_default_dict.update({
        ('individual_sequence_id', 1, int):
            [None, 'n', 1, 'The individual_sequence id of the input sequence file.', ],\
        ('parent_individual_sequence_file_id', 1, int):
            [None, '', 1, 'ID of the parent of this filtered individual_sequence_file'],\
        ('parse_args_from_commandline', 0, int): [1, '', 1, 
            'If 1, parse args from the commandline, not from the input file.', ],
        ('outputDir', 1, ): ['', '', 1, 'output folder to which split files from inputDir will be moved', ],\
        })
    parse_args_from_commandline = 1

    def __init__(self, **keywords):
        """
        2011-7-11
        """
        RegisterAndMoveSplitSequenceFiles.__init__(self, **keywords)

    def run(self):
        """
        2011-7-11
        """

        if self.debug:
            import pdb
            pdb.set_trace()

        db_main = self.db_main
        session = db_main.session
        session.begin()    # no transaction for input node as there is no data insertion

        individual_sequence = db_main.queryTable(
            SunsetDB.IndividualSequence).get(self.individual_sequence_id)

        if not os.path.isdir(self.outputDir):
            os.makedirs(self.outputDir)

        parent_individual_sequence_file = db_main.queryTable(
            SunsetDB.IndividualSequenceFile).get(
                self.parent_individual_sequence_file_id)
        if not parent_individual_sequence_file:
            sys.stderr.write("parent individual_sequence_file %s doesn't exist in db.\n"%(self.parent_individual_sequence_file_id))
            sys.exit(4)

        db_entry = db_main.copyParentIndividualSequenceFile(
            parent_individual_sequence_file=parent_individual_sequence_file,
            individual_sequence_id=self.individual_sequence_id,
            quality_score_format='Standard', filtered=1)
        if db_entry:
            # move the file
            inputDir, filename = os.path.split(self.inputFname)
            exitCode = db_main.moveFileIntoDBAffiliatedStorage(
                db_entry=db_entry, filename=filename, inputDir=inputDir,
                outputDir=self.outputDir,
                relativeOutputDir=individual_sequence.path,
                shellCommand='cp -rL', srcFilenameLs=self.srcFilenameLs,
                dstFilenameLs=self.dstFilenameLs,
                constructRelativePathFunction=None)
            if exitCode != 0:
                sys.stderr.write("Error: moveNewISQFileIntoDBStorage() exits with code=%s.\n"%(exitCode))
                self.sessionRollback(session)
                # delete all recorded target files
                self.cleanUpAndExitOnFailure(exitCode=exitCode)
        else:
            sys.stderr.write("Error: IndividualSequenceFile db entry is None.\n")
            sys.exit(3)

        outputMsg = "file %s was added into db.\n"%(self.inputFname)
        sys.stderr.write(outputMsg)

        if self.logFilename:
            outf = open(self.logFilename, 'w')
            outf.write(outputMsg)
            outf.close()

        if self.commit:
            try:
                session.commit()
                # delete all source files
                self.rmGivenFiles(filenameLs=self.srcFilenameLs)
            except Exception as e:
                sys.stderr.write('Except type: %s\n'%repr(sys.exc_info()))
                import traceback
                traceback.print_exc()
                # delete all recorded target files
                self.cleanUpAndExitOnFailure(exitCode=3)
        else:
            self.sessionRollback(session)
            # delete all target files
            self.cleanUpAndExitOnFailure(exitCode=0)


if __name__ == '__main__':
    main_class = AddFilteredSequences2DB
    po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
    instance = main_class(**po.long_option2value)
    instance.run()
