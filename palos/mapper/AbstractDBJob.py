#!/usr/bin/env python3
"""
2012.6.5
    abstract class for db-interacting pegasus jobs
"""
import sys, os
import csv
from palos import ProcessOptions, getListOutOfStr, PassingData, utils
from palos.mapper import AbstractMapper
from . AbstractMapper import AbstractMapper

class AbstractDBJob(AbstractMapper):
    __doc__ = __doc__
    option_default_dict = AbstractMapper.option_default_dict.copy()
    #option_default_dict.pop(('inputFname', 1, ))
    option_default_dict.update({
        ('logFilename', 0, ): [None, '', 1, 'File to hold logs. '
            'Useful if a job is at the end of a pegasus workflow and '
            'has no other output file.'
            'Because otherwise, pegasus optimization will get rid of this job.'
            '(no output file, why need it?)'],\
        ('commit', 0, int):[0, '', 0, 'commit db transaction'],\
        ('sshTunnelCredential', 0, ): ['', '', 1, \
            'a ssh credential to allow machine to access db server. '
            'polyacti@login3, yuhuang@hpc-login2. '
            'If empty or port is empty, no tunnel',],
        })
    option_default_dict.update(AbstractMapper.db_option_dict)
    
    def __init__(self, inputFnameLs=None, **keywords):
        """
        2011-7-11
        """
        #self.connectDB() called within its __init__()
        AbstractMapper.__init__(self, inputFnameLs=inputFnameLs, **keywords)
        
        #2012.7.4 keep track of all the source&destination files,
        #  used by moveNewISQFileIntoDBStorage()
        self.srcFilenameLs = []
        self.dstFilenameLs = []
        #2012.11.18
        if getattr(self, "logFilename", None):
            self.logF = open(self.logFilename, 'w')
        else:
            self.logF = None
    
    def parseInputFile(self, inputFname=None, **keywords):
        """
        2013.08.23
            if a program is adding a file to db-affiliated storage,
             this is used for parsing.
        """
        return PassingData()
    
    def connectDB(self):
        """
        2012.5.11
            place holder.
        """
        pass
    
    def rmGivenFiles(self, filenameLs=[], rmCommand='rm -rf'):
        """
        2012.7.4
            delete all files in filenameLs
        """
        sys.stderr.write("Deleting %s files ...\n"%(len(filenameLs)))
        for filename in filenameLs:
            commandline = '%s %s'%(rmCommand, filename)
            return_data = utils.runLocalCommand(commandline, 
                report_stderr=True, report_stdout=True)
            if return_data.stderr_content:
                sys.stderr.write("commandline %s failed: %s\n"%(commandline,
                    return_data.stderr_content))
        sys.stderr.write(".\n")
    
    def cleanUpAndExitOnFailure(self, exitCode=1):
        """
        2012.7.13 an exit function when the program failed somewhere
        """
        #delete all target files.
        self.rmGivenFiles(filenameLs=self.dstFilenameLs)
        sys.exit(exitCode)
    
    def cleanUpAndExitOnSuccess(self, exitCode=0):
        """
        2012.7.13  an exit function when the program succeeded in the end
        """
        sys.exit(exitCode)
    
    def outputLogMessage(self, logMessage=None, logFilename=None, logF=None):
        """
        2012.11.18
            1. do not close _logF
            2. use self.logF if None available.
            3. append to logFilename if the latter is given.
        2012.7.17
        """
        _logF = None
        if logF is None:
            if logFilename:
                _logF = open(logFilename, 'a')
            elif self.logF:
                _logF = self.logF
        if _logF:
            _logF.write(logMessage)
    
    def closeLogF(self):
        """
        2012.11.18
        """
        if getattr(self, 'logF', None):
            self.logF.close()
    
    def sessionRollback(self, session=None):
        """
        2013.04.05
        wrap try...except around it because sometimes db connection is gone
            halfway through.
        and it causes program to exit without proper cleanup
        """
        try:
            session.rollback()
        except:
            sys.stderr.write('Except type: %s\n'%repr(sys.exc_info()))
            import traceback
            traceback.print_exc()

    def run(self):
        pass


    