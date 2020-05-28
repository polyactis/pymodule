#!/usr/bin/env python3
"""
This program sums individual columns from all input files based on keys from the keyColumn.

Examples:
    #testing merge three identical genotype files
    %s -k 0 -v 4,5 -o /tmp/test.tsv trio_inconsistency_summary_hist_homo_het.tsv
        
"""
import sys, os
__doc__ = __doc__%(sys.argv[0])
import logging
from palos import ProcessOptions, figureOutDelimiter, utils, PassingData
from palos.io.MatrixFile import MatrixFile
from ReduceMatrixByMergeColumnsWithSameKey import ReduceMatrixByMergeColumnsWithSameKey

class ReduceMatrixByChosenColumn(ReduceMatrixByMergeColumnsWithSameKey):
    __doc__ = __doc__
    option_default_dict = ReduceMatrixByMergeColumnsWithSameKey.option_default_dict.copy()

    def __init__(self, inputFnameLs, **keywords):
        """
        """
        ReduceMatrixByMergeColumnsWithSameKey.__init__(self, inputFnameLs, **keywords)
        
        if self.valueColumnLs:
            self.valueColumnLs = utils.getListOutOfStr(self.valueColumnLs, data_type=int)
        else:
            self.valueColumnLs = []
    
    def handleNewHeader(self, oldHeader, newHeader, keyColumnLs, valueColumnLs,
        keyColumnSet=None):
        """
        """
        originalHeaderLength = len(oldHeader)
        if len(newHeader)==0:
            self.appendSelectedCellIntoGivenList(newHeader, oldHeader, keyColumnLs)
            self.appendSelectedCellIntoGivenList(newHeader, oldHeader, valueColumnLs)
        return newHeader
    
    def handleValueColumns(self, row, key2dataLs=None, keyColumnLs=[],
        valueColumnLs=[], **keywords):
        """
        initiate key2dataLs[key]=[] and extend it by request
        """
        key = self.generateKey(row, keyColumnLs)
        if key not in key2dataLs:
            key2dataLs[key] = []	#0]*len(valueColumnLs)
        for i in range(len(valueColumnLs)):
            columnIndex = valueColumnLs[i]
            if columnIndex<len(row):
                if len(key2dataLs[key])<=i:
                    #2012.1.17 extend it upon request.
                    key2dataLs[key] += [0]*(i+1-len(key2dataLs[key]))
                value = float(row[columnIndex])
                key2dataLs[key][i] = key2dataLs[key][i] + value
        

    def traverse(self):
        """
        """
        newHeader = []
        key2dataLs = {}
        #key is the keyColumn,
        #  dataLs corresponds to the sum of each column from valueColumnLs 
        delimiter = None
        for inputFname in self.inputFnameLs:
            if not os.path.isfile(inputFname):
                if self.exitNonZeroIfAnyInputFileInexistent:
                    sys.exit(3)
                else:
                    continue
            reader = None
            try:
                inputFile = utils.openGzipFile(inputFname)
                delimiter = figureOutDelimiter(inputFile)
                reader = MatrixFile(file_handle=inputFile, delimiter=delimiter)
            except:
                logging.error(f'Except type: {sys.exc_info()}')
                import traceback
                traceback.print_exc()
            
            try:
                header = next(reader)
                self.handleNewHeader(header, newHeader, self.keyColumnLs,
                    self.valueColumnLs, keyColumnSet=self.keyColumnSet)
                if self.noHeader:
                    inputFile.seek(0)
                    reader = MatrixFile(file_handle=inputFile, delimiter=delimiter)
            except:
                logging.error(f'Except type: {sys.exc_info()}')
                import traceback
                traceback.print_exc()
            
            if reader is not None:
                for row in reader:
                    try:
                        self.handleValueColumns(row, key2dataLs=key2dataLs,
                            keyColumnLs=self.keyColumnLs, valueColumnLs=self.valueColumnLs)
                    except:
                        #in case something wrong (i.e. file is empty)
                        logging.error(f'Ignore this row: {row}.')
                        logging.error(f'Except type: {sys.exc_info()}')
                        import traceback
                        traceback.print_exc()
                del reader
        if self.noHeader:
            newHeader = None
        returnData = PassingData(key2dataLs=key2dataLs, delimiter=delimiter,
            header=newHeader)
        return returnData
    
    def run(self):
        
        if self.debug:
            import pdb
            pdb.set_trace()
        
        returnData = self.traverse()
        
        self.outputFinalData(self.outputFname, returnData.key2dataLs,
            returnData.delimiter, header=returnData.header)

if __name__ == '__main__':
    main_class = ReduceMatrixByChosenColumn
    po = ProcessOptions(sys.argv, main_class.option_default_dict,
        error_doc=main_class.__doc__)
    instance = main_class(po.arguments, **po.long_option2value)
    instance.run()
