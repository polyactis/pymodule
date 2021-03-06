#!/usr/bin/env python3
"""
This program takes mean/median/stdev of values of chosen columns
    over all input files with the same key combinations.

Examples:
    %s -k 0 -v 3 -o output.tsv
        trio_inconsistency.tsv trio_inconsistency1.tsv trio_inconsistency2.tsv

"""
import sys, os, math
__doc__ = __doc__%(sys.argv[0])

from palos import ProcessOptions, figureOutDelimiter, utils, PassingData
import csv, numpy
from ReduceMatrixByChosenColumn import ReduceMatrixByChosenColumn

class ReduceMatrixByAverageColumnsWithSameKey(ReduceMatrixByChosenColumn):
    __doc__ = __doc__
    option_default_dict = ReduceMatrixByChosenColumn.option_default_dict.copy()

    def __init__(self, inputFnameLs, **keywords):
        """
        """
        ReduceMatrixByChosenColumn.__init__(self, inputFnameLs, **keywords)
    
    def handleValueColumns(self, row, key2dataLs=None, keyColumnLs=[],
        valueColumnLs=[], **keywords):
        """
        """
        key = self.generateKey(row, keyColumnLs)
        if key not in key2dataLs:
            key2dataLs[key] = []
            for i in range(len(valueColumnLs)):
                #[[]]*len(valueColumnLs) will fail. every [] is a reference to one [].
                key2dataLs[key].append([])
        for i in range(len(valueColumnLs)):
            columnIndex = valueColumnLs[i]
            if columnIndex<len(row):
                value = float(row[columnIndex])
                key2dataLs[key][i].append(value)
    
    def avgKey2DataLs(self, key2dataLs, no_of_key_columns=1, header=[]):
        """
        1. take mean/median/stdev of every cell in dataLs,
        2. modify newHeader to reflect that
        """
        print(f"Averaging key2dataLs ({len(key2dataLs)} entries ) ...", flush=True)
        newKey2DataLs = {}
        newHeader = []
        keyColHeader = header[:no_of_key_columns]
        valueColHeader = header[no_of_key_columns:]
        newValueColHeader = []
        no_of_value_columns = len(valueColHeader)
        for i in range(no_of_value_columns):
            valueColName = valueColHeader[i]
            newValueColHeader += ['mean_%s'%(valueColName),
                'median_%s'%(valueColName), 'stdev_%s'%(valueColName)]
        
        for key, dataLs in key2dataLs.items():
            if key not in newKey2DataLs:
                newKey2DataLs[key] = []
            no_of_value_columns = len(dataLs)
            for i in range(no_of_value_columns):
                meanValue = numpy.mean(dataLs[i])
                medianValue = numpy.median(dataLs[i])
                stdev = numpy.std(dataLs[i])
                newKey2DataLs[key] += [meanValue, medianValue, stdev]
        print(f"Done.", flush=True)
        return PassingData(key2dataLs= newKey2DataLs, header=keyColHeader + newValueColHeader)
    
    def run(self):
        
        if self.debug:
            import pdb
            pdb.set_trace()
        
        returnData = self.traverse()
        newReturnData = self.avgKey2DataLs(returnData.key2dataLs,
            no_of_key_columns=len(self.keyColumnLs), header=returnData.header)
        self.outputFinalData(self.outputFname, newReturnData.key2dataLs,
            returnData.delimiter, header=newReturnData.header)

if __name__ == '__main__':
    main_class = ReduceMatrixByAverageColumnsWithSameKey
    po = ProcessOptions(sys.argv, main_class.option_default_dict,
        error_doc=main_class.__doc__)
    instance = main_class(po.arguments, **po.long_option2value)
    instance.run()
