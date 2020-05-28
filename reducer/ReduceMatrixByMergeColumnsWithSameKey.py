#!/usr/bin/env python3
"""
This program merge lines keyed by keyColumnLs.
    The value columns are merged to a list and outputted by a custom delimiter.
If one input file misses some key combinations,
    those key combinations will have empty data.
    
All input files must have the keys at the same column(s).

Examples:
    #test-merge three identical genotype files
    %s -o /tmp/ccc.tsv --inputDelimiter tab
        /tmp/call_1.tsv /tmp/call_1.tsv /tmp/call_1.tsv
    
"""
import sys, os, math
__doc__ = __doc__%(sys.argv[0])

import copy
import logging
from palos import ProcessOptions, figureOutDelimiter, utils, PassingData
from palos.io.MatrixFile import MatrixFile
from palos.reducer.AbstractReducer import AbstractReducer

class ReduceMatrixByMergeColumnsWithSameKey(AbstractReducer):
    __doc__ = __doc__
    option_default_dict = copy.deepcopy(AbstractReducer.option_default_dict)
    option_default_dict.update({
        ("keyColumnLs", 1, ): [0, 'k', 1, 
            'index(es) of the key in each input file. '
            'Must be same. comma/dash-separated. i.e. 0-2,4 '],
        ("keyHeaderLs", 0, ): [0, '', 1, 
            'header(s) of the key. comma-separated'],
        ('valueColumnLs', 1, ):["1", 'v', 1, 
            'comma/tab-separated list, specifying columns '
            'from which to aggregate total value by key'],
        })

    def __init__(self, inputFnameLs=None, **keywords):
        """
        """
        AbstractReducer.__init__(self, inputFnameLs=inputFnameLs, **keywords)
        if self.keyColumnLs:
            self.keyColumnLs = utils.getListOutOfStr(self.keyColumnLs,
                data_type=int)
        else:
            self.keyColumnLs = []
        if self.keyHeaderLs:
            self.keyHeaderLs = self.keyHeaderLs.split(',')
        else:
            self.keyHeaderLs = []
        
        self.keyColumnSet = set(self.keyColumnLs)
    
    def appendSelectedCellIntoGivenList(self, givenLs=[], inputLs=[],
        indexLs=[]):
        """
        """
        for columnIndex in indexLs:
            if columnIndex<len(inputLs):
                givenLs.append(inputLs[columnIndex])
        return givenLs
    
    def generateKey(self, row, keyColumnLs):
        """
        make sure columnIndex is >=0
        """
        keyLs = []
        for columnIndex in keyColumnLs:
            if columnIndex<len(row) and columnIndex>=0:
                keyLs.append(row[columnIndex])
        key = tuple(keyLs)
        return key
    
    def outputFinalData(self, outputFname, key2dataLs=None, delimiter=None,
        header=None):
        """
        header output is not dependent on key2dataLs anymore 
        """
        writer = MatrixFile(path=outputFname, delimiter=delimiter, mode='w')
        if header and delimiter:
            writer.writerow(header)
        if key2dataLs and delimiter:
            keyLs = sorted(key2dataLs)
            for key in keyLs:
                dataLs = key2dataLs.get(key)
                writer.writerow(list(key) + dataLs)
        writer.close()
    
    def handleNewHeader(self, oldHeader=None, newHeader=None, keyColumnLs=None,
        valueColumnLs=None, keyColumnSet=None):
        """
        """
        originalHeaderLength = len(oldHeader)
        if len(newHeader)==0:	#add the key columns into the new header
            self.appendSelectedCellIntoGivenList(newHeader, oldHeader, keyColumnLs)
        for i in range(originalHeaderLength):
            if i not in keyColumnSet:
                valueColumnLs.append(i)
        self.appendSelectedCellIntoGivenList(newHeader, oldHeader, valueColumnLs)
        return newHeader
    
    def handleValueColumns(self, row, key2dataLs=None, keyColumnLs=[],
        valueColumnLs=[], noOfDataColumnsFromPriorFiles=None,
        visitedKeySet=None):
        """
        """
        key = self.generateKey(row, keyColumnLs)
        if key not in key2dataLs:
            key2dataLs[key] = ['']*noOfDataColumnsFromPriorFiles
        visitedKeySet.add(key)
        
        for columnIndex in valueColumnLs:
            key2dataLs[key].append(row[columnIndex])

    def traverse(self):
        """
        """
        newHeader = []
        key2dataLs = {}
        #key is the keyColumn,
        #  dataLs corresponds to the sum of each column from valueColumnLs 
        noOfDataColumnsFromPriorFiles = 0
        for inputFname in self.inputFnameLs:
            if not os.path.isfile(inputFname):
                if self.exitNonZeroIfAnyInputFileInexistent:
                    logging.error(f'{inputFname} does not exist.')
                    sys.exit(3)
                else:
                    continue
            reader = None
            try:
                inputFile = utils.openGzipFile(inputFname)
                if self.inputDelimiter is None or self.inputDelimiter=='':
                    self.inputDelimiter = figureOutDelimiter(inputFile)
                reader = MatrixFile(file_handle=inputFile,
                    delimiter=self.inputDelimiter)
            except:
                logging.error(f'Except type: {sys.exc_info()}')
                import traceback
                traceback.print_exc()
            
            valueColumnLs = []
            try:
                header = next(reader)
                self.handleNewHeader(header, newHeader, self.keyColumnLs,
                    valueColumnLs, keyColumnSet=self.keyColumnSet)
                if self.noHeader:
                    inputFile.seek(0)
                    reader = MatrixFile(file_handle=inputFile,
                        delimiter=self.inputDelimiter)
            except:
                #in case something wrong (i.e. file is empty)
                logging.error(f'Except type: {sys.exc_info()}')
                import traceback
                traceback.print_exc()
            
            if reader is not None and valueColumnLs:
                visitedKeySet = set()
                for row in reader:
                    try:
                        self.handleValueColumns(row, key2dataLs=key2dataLs,
                            keyColumnLs=self.keyColumnLs,
                            valueColumnLs=valueColumnLs,
                            noOfDataColumnsFromPriorFiles=noOfDataColumnsFromPriorFiles,
                            visitedKeySet=visitedKeySet)
                    except:
                        logging.error(f'Ignore this row: {row}.')
                        logging.error(f'Except type: {sys.exc_info()}')
                        import traceback
                        traceback.print_exc()
                del reader
                #append empty data to keys who are missing in the current file.
                totalKeySet = set(key2dataLs.keys())
                unvisitedKeySet = totalKeySet - visitedKeySet
                for key in unvisitedKeySet:
                    for i in valueColumnLs:
                        key2dataLs[key].append('')
            noOfDataColumnsFromPriorFiles += len(valueColumnLs)
        if self.noHeader:
            newHeader = None
        returnData = PassingData(key2dataLs=key2dataLs,
            delimiter=self.inputDelimiter, header=newHeader)
        return returnData
    
    def run(self):
        
        if self.debug:
            import pdb
            pdb.set_trace()
        
        returnData = self.traverse()
        self.outputFinalData(self.outputFname, key2dataLs=returnData.key2dataLs, 
            delimiter=returnData.delimiter, header=returnData.header)

if __name__ == '__main__':
    main_class = ReduceMatrixByMergeColumnsWithSameKey
    po = ProcessOptions(sys.argv, main_class.option_default_dict,
        error_doc=main_class.__doc__)
    instance = main_class(po.arguments, **po.long_option2value)
    instance.run()
