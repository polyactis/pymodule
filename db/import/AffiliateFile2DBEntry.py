#!/usr/bin/env python3
"""
This program would affiliate a file to an existing database table entry.
    store its path in db,  calculate file size, etc.

Examples:
    %s  -i  OneLibAlignment/6246_depth.tsv.gz
        --logFilename  OneLibAlignment/6246_depth_2db.log
        --db_entry_id 6246 --data_dir /u/home/eeskin/polyacti/NetworkData/vervet/db/
        --drivername postgresql --hostname localhost --dbname vervetdb --db_user yh
        --schema public --tableClassName IndividualAlignment
        --filePathColumnName path_to_depth_file
        --fileSizeColumnName depth_file_size
        --outputFileRelativePath individual_alignment/6246_1517_VEB1010_GA_vs_3488_by_method6_realigned1_reduced0_p4475_m227_depth.tsv.gz

"""

import sys, os, math
__doc__ = __doc__%(sys.argv[0])

import copy
import logging
from palos import ProcessOptions
from palos.mapper.AbstractSunsetMapper import AbstractSunsetMapper as ParentClass
from palos.db import SunsetDB

class AffiliateFile2DBEntry(ParentClass):
    __doc__ = __doc__
    option_default_dict = copy.deepcopy(ParentClass.option_default_dict)
    #option_default_dict.pop(('inputFname', 0, ))
    option_default_dict.pop(('outputFname', 0, ))
    option_default_dict.pop(('outputFnamePrefix', 0, ))
    option_default_dict.update({
        #('inputDir', 1, ): ['', 'i', 1, 
        #   'input folder that contains split fastq files', ],\
        ('db_entry_id', 1, int):[None, '', 1, 
            'ID of the db entry with which the input file is to be affiliated with'],\
        ('tableClassName', 1, ):[None, '', 1, 'table class name'],\
        ('filePathColumnName', 1, ):[None, '', 1, 
            'column name to store file path'],\
        ('fileSizeColumnName', 0, ):[None, '', 1, 
            'column name to store file size'],\
        ('outputFileRelativePath', 1, ):[None, '', 1, 
            'the relative path of the output file, to be attached to self.data_dir'],\
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
        logMessage = "Affiliating file %s to db table %s, entry id=%s."%(
            self.inputFname, self.tableClassName,\
            self.db_entry_id)
        print(f'{logMessage}', flush=True)
        self.outputLogMessage(f'{logMessage}\n')
        
        if os.path.isfile(self.inputFname):
            TableClass = getattr(SunsetDB, self.tableClassName, None)
            db_entry = self.db_main.queryTable(TableClass).get(self.db_entry_id)
            existing_db_entry_file_path = getattr(db_entry, 
                self.filePathColumnName, None)
            if db_entry and existing_db_entry_file_path and \
                    os.path.isfile(os.path.join(data_dir, existing_db_entry_file_path)):
                logging.error(" Another file %s is already in db."%\
                    (existing_db_entry_file_path))
                self.sessionRollback(session)
                self.cleanUpAndExitOnFailure(exitCode=3)
            
            try:
                #copy more files
                logMessage = self.db_main.copyFileWithAnotherFilePrefix(
                    inputFname=self.inputFname,
                    outputFileRelativePath=self.outputFileRelativePath,
                    outputDir=self.data_dir,
                    srcFilenameLs=self.srcFilenameLs,
                    dstFilenameLs=self.dstFilenameLs)
                print(f'{logMessage}', flush=True)
                self.outputLogMessage(f'{logMessage}\n')
                setattr(db_entry, self.filePathColumnName, self.outputFileRelativePath)
                session.add(db_entry)
                session.flush()
                self.db_main.updateDBEntryPathFileSize(db_entry=db_entry,
                    data_dir=self.data_dir, absPath=None,
                    file_path_column_name=self.filePathColumnName,
                    file_size_column_name=self.fileSizeColumnName)
                
                ## 2012.7.17 commented out because md5sum is calculated above
                #db_main.updateDBEntryMD5SUM(db_entry=genotypeFile, data_dir=data_dir)
            except:
                sys.stderr.write('Except type: %s\n'%repr(sys.exc_info()))
                import traceback
                traceback.print_exc()
                self.sessionRollback(session)
                self.cleanUpAndExitOnFailure(exitCode=5)
        else:
            logMessage = "%s doesn't exist."%(inputFileRealPath)
            print(f'{logMessage}', flush=True)
            self.outputLogMessage(f'{logMessage}\n')
        
        if self.commit:
            try:
                session.flush()
                session.commit()
            except:
                sys.stderr.write('Except type: %s\n'%repr(sys.exc_info()))
                import traceback
                traceback.print_exc()
                self.cleanUpAndExitOnFailure(exitCode=3)
        else:
            #delete all target files but exit gracefully (exit 0)
            self.sessionRollback(session)
            self.cleanUpAndExitOnFailure(exitCode=0)

if __name__ == '__main__':
    main_class = AffiliateFile2DBEntry
    po = ProcessOptions(sys.argv, main_class.option_default_dict,
        error_doc=main_class.__doc__)
    instance = main_class(po.arguments, **po.long_option2value)
    instance.run()