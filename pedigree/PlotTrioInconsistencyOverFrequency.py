#!/usr/bin/env python3
"""
2011-9-29
"""

import sys, os, math
import matplotlib; matplotlib.use("Agg")	#to disable pop-up requirement
import csv
from palos import ProcessOptions, PassingData, getColName2IndexFromHeader, figureOutDelimiter
from palos import yh_matplotlib
import pylab, random
from palos.plot.AbstractPlot import AbstractPlot


class PlotTrioInconsistencyOverFrequency(AbstractPlot):
    __doc__ = __doc__
    option_default_dict = AbstractPlot.option_default_dict.copy()
    option_default_dict[('xColumnHeader', 1, )][0] = 'frequency'
    option_default_dict[('xColumnPlotLabel', 0, )][0] = 'frequency'
    option_default_dict[('whichColumnPlotLabel', 0, )][0] = 'inconsistent rate'
    def __init__(self, inputFnameLs=None, **keywords):
        """
        """
        from palos import ProcessOptions
        self.ad = ProcessOptions.process_function_arguments(keywords, 
            self.option_default_dict, error_doc=self.__doc__,
            class_to_have_attr=self)
        self.inputFnameLs = inputFnameLs
    
    @classmethod
    def trioInconsistentRateFileWalker(cls, inputFname, processFunc=None, minNoOfTotal=100, run_type=1):
        """
        only skip except during file opening, not file reading
        """
        try:
            reader = csv.reader(open(inputFname), delimiter=figureOutDelimiter(inputFname))
            header = next(reader)
            col_name2index = getColName2IndexFromHeader(header, skipEmptyColumn=True)
        except:
            sys.stderr.write('Except type: %s\n'%repr(sys.exc_info()))
            import traceback
            traceback.print_exc()
            return
        inconsistent_rate_index = col_name2index.get("inconsistency")
        if run_type==1:
            index_of_x_data = col_name2index.get("stopFrequency")
        elif run_type==2:
            index_of_x_data = col_name2index.get("stop")
        else:
            sys.stderr.write("Unsupported run_type %s in trioInconsistentRateFileWalker().\n"%(run_type))
            sys.exit(3)
        index_of_no_of_total = col_name2index.get("no_of_total")
        inconsistent_rate_ls = []
        x_ls = []
        for row in reader:
            if self.samplingRate<1 and self.samplingRate>=0:
                r = random.random()
                if r>self.samplingRate:
                    continue
            no_of_total = int(float(row[index_of_no_of_total]))
            if no_of_total<=minNoOfTotal:
                continue
            inconsistency = float(row[inconsistent_rate_index])
            inconsistent_rate_ls.append(inconsistency)
            x_data = float(row[index_of_x_data])
            x_ls.append(x_data)
        processFunc(x_ls, inconsistent_rate_ls)
        del reader
    
    def plotXY(self, x_ls, y_ls, ):
        """
        2011-9-30
        """
        pylab.plot(x_ls, y_ls, self.formatString)
        
    
    def run(self):
        
        if self.debug:
            import pdb
            pdb.set_trace()
        
        pylab.clf()
        
        for inputFname in self.inputFnameLs:
            if os.path.isfile(inputFname):
                self.trioInconsistentRateFileWalker(inputFname,
                    processFunc=self.plotXY, minNoOfTotal=self.minNoOfTotal,
                    run_type=1)
        
        if self.title is None:
            title = " %s refs"%(len(self.inputFnameLs))
        else:
            title = self.title
        
        pylab.title(title)
        self.handleXLabel()
        self.handleYLabel()
        
        pylab.savefig(self.outputFname, dpi=self.figureDPI)
        sys.stderr.write("\n")


if __name__ == '__main__':
    main_class = PlotTrioInconsistencyOverFrequency
    po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
    instance = main_class(po.arguments, **po.long_option2value)
    instance.run()
