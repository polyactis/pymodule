#!/usr/bin/env python3
"""
abstract mapper for Sunset mappers/reducers.
"""
from . AbstractDBJob import AbstractDBJob
import sys, os
import copy
import logging
from palos import ProcessOptions
from palos.db import SunsetDB

class AbstractSunsetMapper(AbstractDBJob):
    __doc__ = __doc__
    option_default_dict = copy.deepcopy(AbstractDBJob.option_default_dict)
    #option_default_dict.pop(('inputFname', 0, ))
    
    def __init__(self, inputFnameLs=None, **keywords):
        """
        """
        #self.connectDB() called within __init__()
        AbstractDBJob.__init__(self, inputFnameLs=inputFnameLs, **keywords)
    
    def connectDB(self):
        """
        """
        self.db_main = SunsetDB.SunsetDB(drivername=self.drivername, 
            hostname=self.hostname,
            dbname=self.dbname, schema=self.schema, port=self.port,
            db_user=self.db_user, db_passwd=self.db_passwd)
        self.db_main.setup(create_tables=False)
    
    def checkIfAlignmentListMatchMethodDBEntry(self, 
        individualAlignmentLs:list=None, methodDBEntry=None, session=None):
        """
        moved from AddVCFFile2DB.py
        """
        #make sure methodDBEntry.individual_alignment_ls is identical to 
        # individualAlignmentLs
        alignmentIDSetInFile = set([alignment.id for alignment in individualAlignmentLs])
        alignmentIDSetInGenotypeMethod = set([alignment.id \
            for alignment in methodDBEntry.individual_alignment_ls])
        if alignmentIDSetInFile!=alignmentIDSetInGenotypeMethod:
            logging.error(f"alignmentIDSetInFile ({repr(alignmentIDSetInFile)}) "
                f"doesn't match alignmentIDSetInFile "
                f"({repr(alignmentIDSetInGenotypeMethod)}).")
            if session:
                session.rollback()
            #delete all target files if there is any
            self.cleanUpAndExitOnFailure(exitCode=2)
            
if __name__ == '__main__':
    main_class = AbstractSunsetMapper
    po = ProcessOptions(sys.argv, main_class.option_default_dict,
        error_doc=main_class.__doc__)
    instance = main_class(**po.long_option2value)
    instance.run()
