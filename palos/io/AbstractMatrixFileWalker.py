#!/usr/bin/env python3
"""
Description:
    Abstract class for programs that walk through a list of matrix-like files.
    Running it will reformat input into tab-delimited tsv matrix file.
    If whichColumn/whichColumnHeader is given, it'll convert it into float or log or -log.
    It combines both input from "-i oneFile.txt" or the trailing standalone arguments.
    It could act as both a mapper (one input) or reducer (multiple input).

"""

import sys, os, math
import random
import logging
from palos import ProcessOptions, getListOutOfStr, PassingData, utils
from palos.mapper.AbstractMapper import AbstractMapper
from palos.io.MatrixFile import MatrixFile
from palos.io.HDF5MatrixFile import HDF5MatrixFile
from palos.io.YHPyTables import YHFile
from palos.ngs.io.VCFFile import VCFFile

class AbstractMatrixFileWalker(AbstractMapper):
    __doc__ = __doc__
    option_default_dict = AbstractMapper.option_default_dict.copy()
    #option_default_dict.update(AbstractMapper.db_option_dict.copy())
    option_default_dict.update({
        ('minNoOfTotal', 1, int): [10, 'M', 1, 'minimum no of data from one file for afterFileFunction() to run'],\
        ('h5TableName', 0, ): [None, '', 1, 'table name in the input file if inputFileFormat is 2'],\
        ('maxNoOfTotal', 0, int): [None, '', 1, 'maximum no of data to sample from one file. if not set, no limit.'],\
        ('samplingRate', 1, float): [1, 's', 1, 'how often you include the data, a probability between 0 and 1.'],\
        ('whichColumn', 0, int): [None, 'w', 1,
            'data from this column (index starting from 0) is plotted as y-axis value'],\
        ('whichColumnHeader', 0, ): [None, 'W', 1, 
            'column header for the data to be plotted as y-axis value, substitute whichColumn'],\
        ('logY', 0, int): [0, '', 1,
            'value 0: nothing; 1: log(); 2: -log(); 3: ln(X); 4: -ln(X). replacing self.logWhichColumn.'],\
        ('valueForNonPositiveYValue', 1, float): [-1, 'v', 1, 
            'default value when log-transformation fails (when value is negative)'],\
        ('missingDataNotation', 0, ): ['NA', '', 1, 
            'coma-separated list of notations for missing data.'
            'missing data will be skipped.'],\
        ('inputFileFormat', 0, int): [1, '', 1, '1: csv-like plain text file; '
            '2: YHPyTables.YHFile; 3: HDF5MatrixFile; 4: VCFFile'],\
        ('outputFileFormat', 0, int): [1, '', 1, '1: csv-like plain text file; '
            '2: YHPyTables.YHFile; 3: HDF5MatrixFile; 4: csv-like matrix, without header'],\
        })
    
    def __init__(self, inputFnameLs=None, **keywords):
        """
        self.missingDataNotation will be processed into a set.
        """
        #self.connectDB() called within its __init__()
        AbstractMapper.__init__(self, inputFnameLs=inputFnameLs, **keywords)
        #if user wants to preserve data in a data structure that is visible throughout reading different files.
        # then use this self.invariantPData.
        self.invariantPData = PassingData(writer=None, headerOutputted=False,
            x_ls = [], y_ls = [], z_ls=[])
        if getattr(self, 'missingDataNotation', None):
            self.missingDataNotation = set(utils.getListOutOfStr(
                self.missingDataNotation, data_type=str, separator2=None))
        else:
            self.missingDataNotation = set()
    
    def connectDB(self):
        """
        split out of __init__() so that derived classes could overwrite this function
        """
        pass
    
    def preFileFunction(self, **keywords):
        """
        """
        if hasattr(self, 'initiatePassingData'):
            return self.initiatePassingData(**keywords)
    
    def initiatePassingData(self, ):
        """
        this function gets called in the beginning of each fileWalker() (for each inputFname)
        """
        pdata = PassingData(x_ls = [], y_ls = [], invariantPData=self.invariantPData)
        #2012.8.16 pass to global data
        self.invariantPData.y_ls = pdata.y_ls
        self.invariantPData.x_ls = pdata.x_ls
        return pdata
    
    def processValue(self, value=None, processType=None,
        valueForNonPositiveValue=None, **keywords):
        """
        processType 0: nothing; 1: log(), 2: -log().
        2012.10.15
            the default value of takeLogarithm depends on (self.logWhichColumn or self.logY)
        """
        if valueForNonPositiveValue is None:
            valueForNonPositiveValue = self.valueForNonPositiveYValue
        value = float(value)
        if processType is not None and processType>=1:
            if value>0:
                if processType==1:
                    value = math.log10(value)
                elif processType==2:
                    value = -math.log10(value)
            else:
                value = valueForNonPositiveValue
        return value
    
    def handleYValue(self, yValue=None, processType=None,
        valueForNonPositiveValue=None, **keywords):
        """
        """
        return self.processValue(value=yValue, processType=processType,
            valueForNonPositiveValue=valueForNonPositiveValue, **keywords)
    
    def _writeHeader(self, header=None, pdata=None, rowDefinition=None):
        """
        called by processHeader() and others (in GenomeMovingAverageStatistics.py)
        """
        if not self.invariantPData.headerOutputted:
            if self.outputFileFormat==1:
                if self.invariantPData.writer and header:
                    self.invariantPData.writer.writerow(header)
            elif getattr(self, 'writer', None) is None and \
                getattr(self.invariantPData, 'writer', None) is None:
                if self.outputFileFormat==2:
                    if not rowDefinition and header:
                        #generate a rowDefinition based on header
                        rowDefinition = []
                        for colID in header:
                            rowDefinition.append((colID, 's2000'))
                    writer = YHFile(self.outputFname, mode='w', rowDefinition=rowDefinition)
                    self.invariantPData.writer = writer
                else:	#for HDF5MatrixFile
                    if not rowDefinition and header:
                        #generate a rowDefinition based on header
                        rowDefinition = []
                        for colID in header:
                            rowDefinition.append((colID, HDF5MatrixFile.varLenStrType))
                    #rowDefinition = [('locus_id','i8'),
                    #   ('chromosome', HDF5MatrixFile.varLenStrType), ('start','i8'), ('stop', 'i8'),
                    #	('score', 'f8'), ('MAC', 'i8'), ('MAF', 'f8')]
                    writer = HDF5MatrixFile(self.outputFname, mode='w', rowDefinition=rowDefinition)
                    self.invariantPData.writer = writer
            else:
                logging.warn("Either self.writer %s, or self.invariantPData.writer %s already exists."%\
                    (getattr(self, 'writer', None), getattr(self.invariantPData, 'writer', None)))
                logging.warn("\t no writer created in processHeader().")
        self.invariantPData.headerOutputted = True
    
    def processHeader(self, header=None, pdata=None, rowDefinition=None):
        """
        only open self.outputFname for write when self.writer and self.invariantPData.writer is not available
        called right after the header of an input file is derived in fileWalker().
        """
        self._writeHeader(header=header, pdata=pdata, rowDefinition=rowDefinition)
    
    
    def processRow(self, row=None, pdata=None):
        """
        return returnValue: 0 if not included or nothing is done on it.
            1 if included or something is carried out on it.
        handles each row in each file, here it replaces the yValue.
        """
        returnValue = 0
        col_name2index = getattr(pdata, 'col_name2index', None)
        y_ls = getattr(pdata, 'y_ls', None)
        if col_name2index and y_ls is not None:
            if self.whichColumnHeader:
                whichColumn = col_name2index.get(self.whichColumnHeader, None)
            elif self.whichColumn:
                whichColumn = self.whichColumn
            else:
                whichColumn = None
            if whichColumn is not None:
                yValue = row[whichColumn]
                if yValue not in self.missingDataNotation:
                    yValue = self.processValue(yValue, processType=self.logY,
                        valueForNonPositiveValue=self.valueForNonPositiveYValue)
                row[whichColumn] = yValue
            if self.invariantPData.writer:
                self.invariantPData.writer.writerow(row)
                returnValue = 1
        return returnValue
    
    def getNumberOfData(self, pdata):
        """
        """
        return len(pdata.y_ls)
        
    def afterFileFunction(self, **keywords):
        """
        """
        if hasattr(self, 'plot'):
            return self.plot(**keywords)
    
    def openOneInputFile(self, inputFname=None):
        """
        2013.09.05 split out of fileWalker() , added VCFFile
        """
        if self.inputFileFormat==2:
            reader = YHFile(inputFname, mode='r', tableName=self.h5TableName)
        elif self.inputFileFormat==3:
            reader = HDF5MatrixFile(inputFname, mode='r')
        elif self.inputFileFormat==4:
            reader = VCFFile(inputFname=inputFname)
        else:
            reader = MatrixFile(inputFname)
        return reader
    
    def fileWalker(self, inputFname=None, preFileFunction=None,
        afterFileFunction=None, processRowFunction=None , run_type=1):
        """
        inputFileFormat==2, support HDF5MatrixFile
        """
        sys.stderr.write("walking through %s ..."%(inputFname))
        counter = 0
        real_counter = 0
        noOfSampled = 0
        if preFileFunction is None:
            preFileFunction = self.preFileFunction
        pdata = preFileFunction()
        if processRowFunction is None:
            processRowFunction = self.processRow
        if afterFileFunction is None:
            afterFileFunction = self.afterFileFunction
        try:
            reader = self.openOneInputFile(inputFname)
            col_name2index = reader.constructColName2IndexFromHeader()
            header = reader.getHeader()
            
            #output the header to the output file if necessary 
            self.processHeader(header=header, pdata=pdata)
            pdata.reader = reader
            pdata.col_name2index = col_name2index
            pdata.filename = inputFname
            
            for row in reader:
                counter += 1
                if self.samplingRate<1 and self.samplingRate>=0:
                    r = random.random()
                    if r>self.samplingRate:
                        continue
                rowReturnValue = processRowFunction(row=row, pdata=pdata)
                if rowReturnValue is None:
                    rowReturnValue = 0
                real_counter += rowReturnValue
                noOfSampled += 1
                if self.maxNoOfTotal and real_counter>self.maxNoOfTotal:
                    break
            if self.getNumberOfData(pdata)>=self.minNoOfTotal:
                afterFileFunction(x_ls=pdata.x_ls, y_ls=pdata.y_ls, pdata=pdata)
            del reader
        except:
            sys.stderr.write('Except type: %s\n'%repr(sys.exc_info()))
            import traceback
            traceback.print_exc()
            sys.exit(3)
        if counter>0:
            fraction = float(noOfSampled)/float(counter)
        else:
            fraction = 0
        sys.stderr.write("%s/%s (%.3f) data sampled. real_counter=%s.\n"%(
            noOfSampled, counter, fraction, real_counter))
    
    def setup(self, **keywords):
        """
        do not open the file if it's a png file.
        run before anything is run.
        """
        writer = None
        if self.outputFileFormat in [1,4]:
            suffix = os.path.splitext(self.outputFname)[1]
            if self.outputFname and suffix!='.png':
                writer = MatrixFile(self.outputFname, mode='w', delimiter='\t')
        else:
            #HDF5MatrixFile
            #can't generate HDF5MatrixFile, because it needs dtypeList
            pass
        #pass it to the invariantPData
        self.invariantPData.writer = writer
        self.writer = writer
    
    def closeFiles(self, **keywords):
        """
        2013.1.27
        """
        if getattr(self.invariantPData, 'writer', None):
            try:
                self.invariantPData.writer.close()
            except:
                sys.stderr.write('Except type: %s\n'%repr(sys.exc_info()))
                import traceback
                traceback.print_exc()
            
            #this would close file too if it's not closed
            del self.invariantPData.writer
        
        if getattr(self, 'writer', None):
            del self.writer
        
    def reduce(self, **keywords):
        """
        run after all files have been walked through
        """
        pass
        
    
    def run(self):
        
        if self.debug:
            import pdb
            pdb.set_trace()
        
        self.setup()
        
        for inputFname in self.inputFnameLs:
            if os.path.isfile(inputFname):
                self.fileWalker(inputFname, afterFileFunction=self.afterFileFunction,
                    run_type=1, processRowFunction=self.processRow)
        
        self.reduce()
        self.closeFiles()
        
            
if __name__ == '__main__':
    main_class = AbstractMatrixFileWalker
    po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
    instance = main_class(po.arguments, **po.long_option2value)
    instance.run()