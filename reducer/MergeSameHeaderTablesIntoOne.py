#!/usr/bin/env python3
"""
This program merges all files with the same header into one while retaining the header.
    It can handle gzipped input files.

Examples:
    #testing merge three identical genotype files
    %s -o /tmp/ccc.tsv /tmp/call_1.tsv /tmp/call_1.tsv /tmp/call_1.tsv
    
    %s --exitNonZeroIfAnyInputFileInexistent -o /tmp/ccc.tsv
        /tmp/call_1.tsv /tmp/call_1.tsv /tmp/call_1.tsv
    
"""
import sys, os
__doc__ = __doc__%(sys.argv[0], sys.argv[0])
import logging
import copy
from palos import ProcessOptions, utils
from palos.reducer.AbstractReducer import AbstractReducer

class MergeSameHeaderTablesIntoOne(AbstractReducer):
    __doc__ = __doc__
    option_default_dict = copy.deepcopy(AbstractReducer.option_default_dict)

    def __init__(self, inputFnameLs=None, **keywords):
        """
        """
        AbstractReducer.__init__(self, inputFnameLs=inputFnameLs, **keywords)
    
    def run(self):
        if self.debug:
            import pdb
            pdb.set_trace()
        
        header = None
        outf = utils.openGzipFile(self.outputFname, 'w')
        for inputFname in self.inputFnameLs:
            print(f"File {inputFname} ... ", flush=True)
            if not os.path.isfile(inputFname):
                if self.exitNonZeroIfAnyInputFileInexistent:
                    logging.error(f"{inputFname} doesn't exist.")
                    sys.exit(3)
                else:
                    continue
            inf = utils.openGzipFile(inputFname, 'r')
            if self.noHeader==0:
                #in the case that every input has a common header
                if not header:
                    #if empty string or None, obtain a header
                    try:
                        header = inf.readline()
                        outf.write(header)
                    except:	#in case something wrong (i.e. file is empty)
                        logging.error('Except type: %s'%repr(sys.exc_info()))
                        import traceback
                        traceback.print_exc()
                        print(sys.exc_info())
                else:
                    #skip the header for other input files
                    try:
                        inf.readline()
                    except:
                        #in case something wrong (i.e. file is empty)
                        logging.error('Except type: %s'%repr(sys.exc_info()))
                        import traceback
                        traceback.print_exc()
                        print(sys.exc_info())
            for line in inf:
                isEmpty = self.isInputLineEmpty(line.strip(), inputFile=inf,
                    inputEmptyType=self.inputEmptyType)
                if not isEmpty:
                    outf.write(line)
            print(f"Done.", flush=True)

if __name__ == '__main__':
    main_class = MergeSameHeaderTablesIntoOne
    po = ProcessOptions(sys.argv, main_class.option_default_dict,
        error_doc=main_class.__doc__)
    instance = main_class(po.arguments, **po.long_option2value)
    instance.run()
