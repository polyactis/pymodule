#!/usr/bin/env python3
"""
This program first sums values of chosen columns over all input files
    with the same key combinations.
In the end, it divides values from first two chosen columns
    and appends it to the output as an extra column.

Examples:
    #testing merge three identical genotype files
    %s -k 0 -v 4,5 -o /tmp/test.tsv trio_inconsistency_summary_hist_homo_het.tsv
    
"""
import sys, os, math
__doc__ = __doc__%(sys.argv[0])

from palos import ProcessOptions, PassingData
from ReduceMatrixByChosenColumn import ReduceMatrixByChosenColumn

class ReduceMatrixBySumSameKeyColsAndThenDivide(ReduceMatrixByChosenColumn):
    __doc__ = __doc__
    option_default_dict = ReduceMatrixByChosenColumn.option_default_dict.copy()
    option_default_dict.update({
        ("operatorType", 1, int): [1, 'p', 1, 
            'For the last column, 1: firstColumnValue/2ndColumnValue;\n'
            '2: 1stColumnValue - 2ndColumnValue; 3: divide each column in '
            '--valueColumnLs .. with a fixed value (--fixedValueDenominator ..).'],
        ("fixedValueDenominator", 0, float): [1, '', 1,
            'The denominator used in operatorType 3.'],
        })
    def __init__(self, inputFnameLs, **keywords):
        """
        """
        ReduceMatrixByChosenColumn.__init__(self, inputFnameLs, **keywords)
    
    def divideTwoColumnsInKey2DataLs(self, key2dataLs, no_of_key_columns=1,
        header=[]):
        """
        1. generate a new column by dividing the values from 1st two cells in dataLs,
        2. modify newHeader to reflect that
        """
        print(f"Averaging key2dataLs ({len(key2dataLs)} entries) ...", flush=True)
        if self.noHeader==0:
            keyColHeader = header[:no_of_key_columns]
            valueColHeader = header[no_of_key_columns:]
            if self.operatorType==3:
                for valueSingleColHeader in valueColHeader:
                    header.append("%s_byFixedValue"%(valueSingleColHeader))
            else:
                if len(valueColHeader)>1:
                    if self.operatorType==2:
                        header.append("%s_subtract_%s"%(valueColHeader[0], valueColHeader[1]))
                    else:
                        header.append('%s_by_%s'%(valueColHeader[0], valueColHeader[1]))
            
        for key, dataLs in key2dataLs.items():
            no_of_value_columns = len(dataLs)
            if self.operatorType==3:
                for i in range(no_of_value_columns):
                    #do not iterate over dataLs (for data in dataLs) as dataLs is expanding ... 
                    data = dataLs[i]
                    if self.fixedValueDenominator!=0:
                        ratio = float(data)/self.fixedValueDenominator
                    else:
                        ratio = -1
                    key2dataLs[key].append(ratio)
            else:
                if no_of_value_columns>1:
                    if self.operatorType==2:
                        ratio = float(dataLs[0]) - float(dataLs[1])
                    else:
                        if dataLs[1]!=0:
                            ratio = dataLs[0]/float(dataLs[1])
                        else:
                            ratio = -1
                    key2dataLs[key].append(ratio)
        print("Done.", flush=True)
        return PassingData(key2dataLs=key2dataLs, header=header)
    
    def run(self):
        
        if self.debug:
            import pdb
            pdb.set_trace()
        
        returnData = self.traverse()
        newReturnData = self.divideTwoColumnsInKey2DataLs(
            returnData.key2dataLs, no_of_key_columns=len(self.keyColumnLs),
            header=returnData.header)
        self.outputFinalData(self.outputFname, newReturnData.key2dataLs,
            returnData.delimiter, header=newReturnData.header)

if __name__ == '__main__':
    main_class = ReduceMatrixBySumSameKeyColsAndThenDivide
    po = ProcessOptions(sys.argv, main_class.option_default_dict,
        error_doc=main_class.__doc__)
    instance = main_class(po.arguments, **po.long_option2value)
    instance.run()
