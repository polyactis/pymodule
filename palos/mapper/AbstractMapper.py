#!/usr/bin/env python3
"""
Description:
    2012.3.6 an abstract class for mapper
"""

import sys, os, math
import csv
from palos import ProcessOptions, getListOutOfStr, PassingData, utils
import re

class AbstractMapper(object):
    __doc__ = __doc__
    db_option_dict = {
        ('drivername', 1,):['postgresql', 'v', 1,	
            'which type of database? mysql or postgresql. For SQLAlchemy', ],
        ('hostname', 1, ): ['localhost', 'z', 1, 'hostname of the db server', ],
        ('dbname', 1, ): ['vervetdb', 'd', 1, 'database name', ],\
        ('schema', 0, ): [None, '', 1, 'database schema name', ],\
        ('db_user', 1, ): [None, 'u', 1, 'database username', ],\
        ('db_passwd', 1, ): [None, '', 1, 'database password', ],\
        ('port', 0, ):[None, '', 1, 'database port number'],\
        ('commit', 0, int):[0, '', 0, 'commit db transaction'],\
        ('data_dir', 0, ):[None, '', 1, 'path to the db-affliated file storage.'],
        }
    genome_db_option_dict = {
        ('genome_drivername', 1,):['postgresql', '', 1, 
            'which type of database: mysql or postgresql',],
        ('genome_hostname', 1, ): ['uclaOffice', '', 1,
            'hostname of the genome db server', ],
        ('genome_dbname', 1, ): ['vervetdb', '', 1, 'genome database name', ],
        ('genome_schema', 0, ): ['genome', '', 1, 
            'genome database schema name',],
        ('genome_db_user', 1, ): ['yh', '', 1, 'genome database username', ],
        ('genome_db_passwd', 1, ): [None, '', 1, 'genome database password',],
        ('genome_port', 0, ):[None, '', 1, 'genome database port number'],
        }
    option_default_dict = {
        ('inputFname', 0, ): [None, 'i', 1, 'The input file.', ],\
        ('inputDelimiter', 0, ): [None, '', 1, 
            'Delimiter of input file: tab, space, or coma. '
            'If not given, it will guess.', ],\
        ("home_path", 1, ): [os.path.expanduser("~"), 'e', 1, \
            'path to the home directory on the working nodes'],
        ('outputFname', 0, ): [None, 'o', 1, 'output file'],\
        ('outputFnamePrefix', 0, ): [None, 'O', 1, \
            'output filename prefix (optional).'],
        ('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],
        ('report', 0, int):[0, 'r', 0, 
            'toggle report, more verbose stdout/stderr.']
        }

    def __init__(self, inputFnameLs=None, **keywords):
        """
        """
        self.ad = ProcessOptions.process_function_arguments(keywords, \
            self.option_default_dict, error_doc=self.__doc__, \
            class_to_have_attr=self)
        if getattr(self, 'outputFname', None) and hasattr(self,'outputFnamePrefix'):
            self.outputFnamePrefix = os.path.splitext(self.outputFname)[0]
        
        self.inputFnameLs = inputFnameLs
        if self.inputFnameLs is None:
            self.inputFnameLs = []
        
        if getattr(self, 'inputFname', None):
            self.inputFnameLs.insert(0, self.inputFname)
        
        #2013.08.14 setup self.inputFname if self.inputFnameLs is not None
        if self.inputFnameLs and not getattr(self, 'inputFname', None):
            self.inputFname = self.inputFnameLs[0]
        
        if getattr(self, 'inputDelimiter', None):
            if self.inputDelimiter=='tab':
                self.inputDelimiter = '\t'
            if self.inputDelimiter=='space':
                self.inputDelimiter = ' '
            if self.inputDelimiter=='coma':
                self.inputDelimiter = ','
        self.connectDB()
        #2012.8.14 an expression that searches any character in a string
        self.p_char = re.compile(r'[a-df-zA-DF-Z\-]$')
        #no 'e' or 'E', used in scientific number, add '-' and append '$'
    
    def insertHomePath(self, inputPath, home_path):
        """
        2012.8.19 copied from AbstractWorkflow
        """
        if inputPath.find('%s')!=-1:
            inputPath = inputPath%home_path
        return inputPath
    
    def connectDB(self):
        """
        2012.5.11
            place holder.
        """
        pass
        """
        db_main = SunsetDB.SunsetDB(drivername=self.drivername, \
            hostname=self.hostname, database=self.dbname, schema=self.schema,
            port=self.port,
            username=self.db_user, password=self.db_passwd)
        db_main.setup(create_tables=False)
        self.db_main = db_main
        """
    
    def run(self):
        """
        2013.05.01 place holder
        """
        pass

if __name__ == '__main__':
    main_class = AbstractMapper
    po = ProcessOptions(sys.argv, main_class.option_default_dict,
        error_doc=main_class.__doc__)
    instance = main_class(**po.long_option2value)
    instance.run()	