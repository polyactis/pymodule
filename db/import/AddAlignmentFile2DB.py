#!/usr/bin/env python3
"""
Examples:
    %s
    -i 2278_634_vs_524_by_2_r4043_sequence_dupMarked.bam
    --baiFile 2278_634_vs_524_by_2_r4043_sequence_dupMarked.bam.bai
    --logFilename 2278_634_vs_524_by_2_r4043_sequence_628C2AAXX_6_2db.log
    --individual_alignment_id 2278
    --data_dir /u/home/eeskin/polyacti/NetworkData/vervet/db/
    --drivername postgresql --hostname localhost --dbname vervetdb
    --schema public --db_user yh
    2278_634_vs_524_by_2_r4043_sequence_628C2AAXX_6_dupMarked.metric

Description:
    Add the alignment file into database
        1. register an alignment entry in db
            a. if individual_alignment_id is not None
            b. if parent_individual_alignment_id is not None:
            c. construct alignment using all other arguments
        2. copy the file (& bai file if it exists and
             other files in the commandline arguments) over.
        3. write the log if instructed so (for workflow purpose)

"""

import sys, os
__doc__ = __doc__%(sys.argv[0])

import copy
import logging
from palos import ProcessOptions, utils
from palos.mapper.AbstractSunsetMapper import AbstractSunsetMapper as ParentClass
from palos.db import SunsetDB

class AddAlignmentFile2DB(ParentClass):
    __doc__ = __doc__
    option_default_dict = copy.deepcopy(ParentClass.option_default_dict)
    #option_default_dict.pop(('inputFname', 0, ))
    option_default_dict.pop(('outputFname', 0, ))
    option_default_dict.pop(('outputFnamePrefix', 0, ))
    option_default_dict.update({
        ('baiFile', 1, ):[None, '', 1,
            'the bai index file'],
        ('individual_alignment_id', 0, int):[None, '', 1, \
            'fetch the db individual_alignment based on this ID'],
        ('individual_sequence_id', 0, int):[None, '', 1, 
            'used to construct individual_alignment'],
        ('ref_sequence_id', 0, int):[None, '', 1, 
            'used to construct individual_alignment'],
        ('alignment_method_id', 0, int):[None, '', 1, 
            'used to construct individual_alignment'],
        ('parent_individual_alignment_id', 0, int):[None, '', 1, 
            'the parent ID of individual_alignment.'
            'if given, an individual_alignment db entry will be created as a '
            'copy of this one.'],
        ('mask_genotype_method_id', 0, int):[None, '', 1, 
            'for alignments coming out of base quality recalibration'],
        ('individual_sequence_file_raw_id', 0, int):[None, '', 1, 
            'for library specific alignment'],
        ('local_realigned', 0, int):[0, '', 1, 
            'value for IndividualAlignment.local_realigned'],
        ('read_group', 0, ):[None, '', 1, 
            'value for IndividualAlignment.read_group. if not given, '
            'it calls IndividualAlignment.getReadGroup()'],
        ('format', 0, ):[None, 'f', 1, 'format for GenotypeFile entry'],
        })
    def __init__(self, inputFnameLs=None, **keywords):
        """
        """
        ParentClass.__init__(self, inputFnameLs=inputFnameLs, **keywords)

    def run(self):
        """
        """
        if self.debug:
            import pdb
            pdb.set_trace()
        session = self.db_main.session

        session.begin()
        if not self.data_dir:
            self.data_dir = self.db_main.data_dir
        data_dir = self.data_dir

        inputFileRealPath = os.path.realpath(self.inputFname)
        logMessage = f"Adding file {self.inputFname} to db."
        print(f'{logMessage}', flush=True)
        self.outputLogMessage(f'{logMessage}\n')

        if os.path.isfile(inputFileRealPath):
            if self.individual_alignment_id:
                ind_aln = self.db_main.queryTable(
                    SunsetDB.IndividualAlignment).get(self.individual_alignment_id)
            elif self.parent_individual_alignment_id:
                ind_aln = self.db_main.copyParentIndividualAlignment(
                    parent_individual_alignment_id=self.parent_individual_alignment_id,
                    mask_genotype_method_id=self.mask_genotype_method_id,\
                    data_dir=self.data_dir,
                    local_realigned=self.local_realigned)
            else:
                #alignment for this library of the individual_sequence
                individual_sequence = self.db_main.queryTable(
                    SunsetDB.IndividualSequence).get(self.individual_sequence_id)
                ind_aln = self.db_main.getAlignment(
                    individual_sequence_id=self.individual_sequence_id,
                    path_to_original_alignment=None,
                    sequencer_name=individual_sequence.sequencer.short_name,
                    sequence_type_name=individual_sequence.sequence_type.short_name,
                    sequence_format=individual_sequence.format,
                    ref_individual_sequence_id=self.ref_sequence_id,
                    alignment_method_id=self.alignment_method_id,
                    alignment_format=self.format,
                    individual_sequence_filtered=individual_sequence.filtered,
                    read_group_added=1,
                    data_dir=data_dir,
                    mask_genotype_method_id=self.mask_genotype_method_id,
                    parent_individual_alignment_id=self.parent_individual_alignment_id,
                    individual_sequence_file_raw_id=self.individual_sequence_file_raw_id,
                    local_realigned=self.local_realigned,
                    read_group=self.read_group)
            needSessionFlush = False
            if not ind_aln.path:
                ind_aln.path = ind_aln.constructRelativePath()
                needSessionFlush = True

            if self.mask_genotype_method_id and \
                    ind_aln.mask_genotype_method_id!=self.mask_genotype_method_id:
                ind_aln.mask_genotype_method_id = self.mask_genotype_method_id
                needSessionFlush = True
            if self.individual_sequence_file_raw_id and \
                    ind_aln.individual_sequence_file_raw_id !=\
                        self.individual_sequence_file_raw_id:
                ind_aln.individual_sequence_file_raw_id = \
                    self.individual_sequence_file_raw_id
                needSessionFlush = True

            if needSessionFlush:
                session.add(ind_aln)
                session.flush()

            try:
                md5sum = utils.get_md5sum(inputFileRealPath)
            except:
                logging.error(f'Except type: {repr(sys.exc_info())}')
                import traceback
                traceback.print_exc()
                self.cleanUpAndExitOnFailure(exitCode=4)

            db_entry = self.db_main.queryTable(SunsetDB.IndividualAlignment).\
                filter_by(md5sum=md5sum).first()
            if db_entry and db_entry.id!=ind_aln.id and \
                db_entry.path and os.path.isfile(os.path.join(data_dir, db_entry.path)):
                logging.warn(f"Another file {db_entry.path} with the identical"
                    f" md5sum {md5sum} as this file {inputFileRealPath}, is "
                    f"already in db.")
                self.sessionRollback(session)
                self.cleanUpAndExitOnFailure(exitCode=3)

            if ind_aln.md5sum is None or \
                    ind_aln.md5sum!=md5sum:
                ind_aln.md5sum = md5sum
                session.add(ind_aln)
                session.flush()
            try:
                #move the file and update the db_entry's path as well
                exitCode = self.db_main.moveFileIntoDBAffiliatedStorage(
                    db_entry=ind_aln,
                    filename=os.path.basename(inputFileRealPath),
                    inputDir=os.path.split(inputFileRealPath)[0], 
                    dstFilename=os.path.join(self.data_dir, ind_aln.path),
                    relativeOutputDir=None,
                    shellCommand='cp -rL',
                    srcFilenameLs=self.srcFilenameLs,
                    dstFilenameLs=self.dstFilenameLs,
                    constructRelativePathFunction=\
                        ind_aln.constructRelativePath)
            except:
                logging.error(f'Except in copying {inputFileRealPath} to '
                    f'db-storage with except info: {repr(sys.exc_info())}.')
                import traceback
                traceback.print_exc()
                self.sessionRollback(session)
                self.cleanUpAndExitOnFailure(exitCode=5)

            if exitCode!=0:
                logging.error(f"moveFileIntoDBAffiliatedStorage() exits with "
                    f"code={exitCode}.")
                self.sessionRollback(session)
                self.cleanUpAndExitOnFailure(exitCode=exitCode)
            try:
                if self.inputFnameLs:
                    #copy other files 
                    for inputFname in self.inputFnameLs:
                        if inputFname!=self.inputFname:
                            # make sure it has not been copied.
                            logMessage = self.db_main.copyFileWithAnotherFilePrefix(
                                inputFname=inputFname,
                                filenameWithPrefix=ind_aln.path,
                                outputDir=self.data_dir,
                                srcFilenameLs=self.srcFilenameLs,
                                dstFilenameLs=self.dstFilenameLs)
                            print(f'{logMessage}', flush=True)
                            self.outputLogMessage(f'{logMessage}\n')

                self.db_main.updateDBEntryPathFileSize(
                    db_entry=ind_aln, data_dir=data_dir)

                ## 2012.7.17 commented out because md5sum is calculated above
                #db_main.updateDBEntryMD5SUM(db_entry=genotypeFile, data_dir=data_dir)
                #copy the bai index file if it exists
                #baiFile = f'{self.inputFname}.bai'
                if not os.path.isfile(self.baiFile):
                    logging.error(f"The bam index file {self.baiFile} does not exist!")
                    self.sessionRollback(session)
                    self.cleanUpAndExitOnFailure(exitCode=5)
                else:
                    srcFilename = self.baiFile
                    dstFilename = os.path.join(self.data_dir, \
                        f'{ind_aln.path}.bai')
                    utils.copyFile(srcFilename=srcFilename, dstFilename=dstFilename)
                    logMessage += f"bai file {srcFilename} has been copied "\
                        f"to {dstFilename}.\n"
                    print(f'{logMessage}', flush=True)
                    self.outputLogMessage(f'{logMessage}\n')
                    self.srcFilenameLs.append(srcFilename)
                    self.dstFilenameLs.append(dstFilename)
            except:
                logging.error(f'Except type: {repr(sys.exc_info())}')
                import traceback
                traceback.print_exc()
                self.sessionRollback(session)
                self.cleanUpAndExitOnFailure(exitCode=5)
        else:
            logMessage = f"{inputFileRealPath} doesn't exist."
            print(f'{logMessage}', flush=True)
            self.outputLogMessage(f'{logMessage}\n')

        if self.commit:
            try:
                session.flush()
                session.commit()
                print(f"Commit() successful.", flush=True)
            except:
                logging.error(f'Except type: {repr(sys.exc_info())}')
                import traceback
                traceback.print_exc()
                self.cleanUpAndExitOnFailure(exitCode=3)
        else:
            #delete all target files but exit gracefully (exit 0)
            self.sessionRollback(session)
            self.cleanUpAndExitOnFailure(exitCode=0)


if __name__ == '__main__':
    main_class = AddAlignmentFile2DB
    po = ProcessOptions(sys.argv, main_class.option_default_dict, \
        error_doc=main_class.__doc__)
    instance = main_class(po.arguments, **po.long_option2value)
    instance.run()
