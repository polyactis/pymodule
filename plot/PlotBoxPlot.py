#!/usr/bin/env python3
"""
Examples:
    %s  -s 1 -o ../plot/AAC_tally_barChart.png -W NumberOfLoci -l AAC -x AAC -i 20 -D NumberOfLoci ./AAC_tally.tsv.gz
    

Description:
    this program draws multiple box plots of yValue. each boxplot's data has the same xValue (on X-axis).
"""

import sys, os, math
__doc__ = __doc__%(sys.argv[0])

import matplotlib; matplotlib.use("Agg")	#to disable pop-up requirement
import csv
import numpy, random, pylab
from palos import ProcessOptions, PassingData
from palos.plot import yh_matplotlib
from AbstractPlot import AbstractPlot

class PlotBoxPlot(AbstractPlot):
    __doc__ = __doc__
#						
    option_default_dict = AbstractPlot.option_default_dict.copy()
    def __init__(self, inputFnameLs=None, **keywords):
        """
        """
        AbstractPlot.__init__(self, inputFnameLs=inputFnameLs, **keywords)
        self.boxPlotXLabelList = []
        
    def initiatePassingData(self, ):
        """
        this function gets called in the beginning of each fileWalker() (for each inputFname)
        """
        pdata = PassingData(xValue2yValueLs={}, x_ls=[], y_ls=[], invariantPData=self.invariantPData)
        self.invariantPData.y_ls = pdata.y_ls
        self.invariantPData.x_ls = pdata.x_ls
        return pdata
    
    def getNumberOfData(self, pdata):
        """
        """
        return len(pdata.xValue2yValueLs)
    
    def handleXLim(self,):
        """
        """
        pylab.xlim(xmin=self.xMin-1, xmax=self.xMax+1)
        
    def handleYLim(self,):
        """
        """
        self._handleYLim()
    
    def plot(self, x_ls=None, y_ls=None, pdata=None,
        color_ls=['r', 'g', 'b', 'y','c', 'k'], width = 0.5, **kwargs):
        """
        get called by the end of fileWalker() for each inputFname.
        """
        xValue2yValueLs = getattr(pdata, 'xValue2yValueLs', None)
        xValueList = []
        y2DList = []
        noOfBoxPlotsBefore = len(self.boxPlotXLabelList)
        for xValue, yValueList in xValue2yValueLs.items():
            xValueList.append(xValue)
            self.boxPlotXLabelList.append(xValue)
            y2DList.append(yValueList)
            yMin = min(yValueList)
            self.setGlobalMinVariable(extremeVariableName='yMin', givenExtremeValue=yMin)
            yMax = max(yValueList)
            self.setGlobalMaxVariable(extremeVariableName='yMax', givenExtremeValue=yMax)
        rects_ls = []
        
        #width = 0.75/len(yValue_list)	#these are for multi-series of bar charts.
        #c = color_ls[i%len(color_ls)]
        #pylab.boxplot(y2DList, positions=xValueList)
        # #this changes the x-location of the box plot
        noOfBoxPlotsNow = len(self.boxPlotXLabelList)
        plotObject = pylab.boxplot(y2DList, positions=range(noOfBoxPlotsBefore, noOfBoxPlotsNow))
        
        #self.addPlotLegend(plotObject=plotObject, 
        # legend=os.path.basename(pdata.filename), pdata=pdata)
        
        self.setGlobalMinVariable(extremeVariableName='xMin', givenExtremeValue=0)
        self.setGlobalMaxVariable(extremeVariableName='xMax', givenExtremeValue=noOfBoxPlotsNow-1)
    
    def processRow(self, row=None, pdata=None):
        """
        handle each row in each file
        """
        col_name2index = getattr(pdata, 'col_name2index', None)
        xValue2yValueLs = getattr(pdata, 'xValue2yValueLs', None)
        y_ls = getattr(pdata, 'y_ls', None)
        if col_name2index and xValue2yValueLs is not None:
            if self.whichColumnHeader:
                whichColumn = col_name2index.get(self.whichColumnHeader, None)
            else:
                whichColumn = self.whichColumn
            x_index = col_name2index.get(self.xColumnHeader, None)
            
            xValue = float(row[x_index])
            xValue = self.processValue(xValue, processType=self.logX)
            yValue = float(row[whichColumn])
            yValue = self.processValue(yValue, processType=self.logY)
            if xValue not in xValue2yValueLs:
                xValue2yValueLs[xValue] = []
            xValue2yValueLs[xValue].append(yValue)
            y_ls.append(yValue)
    
    def saveFigure(self, invariantPData=None, **keywords):
        """
        before saving the figure, mark the xticks
        """
        pylab.xticks(range(len(self.boxPlotXLabelList)), self.boxPlotXLabelList)
        AbstractPlot.saveFigure(self, invariantPData=invariantPData, **keywords)

if __name__ == '__main__':
    main_class = PlotBoxPlot
    po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
    instance = main_class(po.arguments, **po.long_option2value)
    instance.run()
