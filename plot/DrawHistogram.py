#!/usr/bin/env python3
"""
Examples:	
    #Draw histogram of fraction of heterozygotes per individual.
    %s --minNoOfTotal 10 --samplingRate 1 -o /tmp/hetPerMonkey_hist.png
        --whichColumnHeader NoOfHet_by_NoOfTotal
        /tmp/homoHetCountPerSample.tsv

Description:
"""

import sys, os, math
__doc__ = __doc__%(sys.argv[0])

import matplotlib; matplotlib.use("Agg")	#to disable pop-up requirement
from palos import ProcessOptions
from palos.plot import yh_matplotlib
import pylab
from palos.plot.AbstractPlot import AbstractPlot

class DrawHistogram(AbstractPlot):
    __doc__ = __doc__
    option_default_dict = AbstractPlot.option_default_dict
    option_default_dict.pop(('xColumnHeader', 1, ))
    option_default_dict.pop(('xColumnPlotLabel', 0, ))
    option_default_dict.update({
        ('logCount', 0, int): [0, 'l', 0, 
            'whether to take log on the y-axis of the histogram, the raw count. '
            'similar effect to yScaleLog'], \
    })
    def __init__(self, inputFnameLs=None, **keywords):
        """
        """
        AbstractPlot.__init__(self, inputFnameLs=inputFnameLs, **keywords)
    
    def plot(self, x_ls=None, y_ls=None, pdata=None, min_no_of_data_points=10,
        max_no_of_bins=30, min_no_of_bins=10, **kwargs):
        """
        get called by the end of fileWalker() for each inputFname.
        """
        no_of_data_points = len(y_ls)
        if no_of_data_points>=min_no_of_data_points:
            no_of_bins = max(min_no_of_bins, min(max_no_of_bins, no_of_data_points/10))
            n, bins, patches = pylab.hist(y_ls, bins=no_of_bins,
                log=self.logCount, alpha=0.6)
            self.addPlotLegend(plotObject=patches[0],
                legend=os.path.basename(pdata.filename), pdata=pdata)
        
    
    def processRow(self, row=None, pdata=None):
        """
        Skip missing data via self.missingDataNotation.
        handles each row in each file
        """
        col_name2index = getattr(pdata, 'col_name2index', None)
        y_ls = getattr(pdata, 'y_ls', None)
        if col_name2index and y_ls is not None:
            if self.whichColumnHeader:
                whichColumn = col_name2index.get(self.whichColumnHeader, None)
            else:
                whichColumn = self.whichColumn
            
            yValue = row[whichColumn]
            if yValue not in self.missingDataNotation:
                yValue = self.processValue(value=yValue, processType=self.logY)
                y_ls.append(yValue)
    
    def handleXLabel(self,):
        """
        the whichColumnPlotLabel is usually reserved for y-axis, but in histogram, it's on x-axis.
        """
        if getattr(self, 'whichColumnPlotLabel', None):
            pylab.xlabel(self.whichColumnPlotLabel)
        else:
            pylab.xlabel(getattr(self, "whichColumnHeader", ""))
    
    def handleYLabel(self,):
        """
        """
        pass
    
    def handleTitle(self,):
        """
        """
        if self.title:
            pylab.title(self.title)
        else:
            title = yh_matplotlib.constructTitleFromDataSummaryStat(self.invariantPData.y_ls)
            pylab.title(title)

if __name__ == '__main__':
    main_class = DrawHistogram
    po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
    instance = main_class(po.arguments, **po.long_option2value)
    instance.run()
