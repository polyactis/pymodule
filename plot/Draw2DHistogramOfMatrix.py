#!/usr/bin/env python3
"""
The input is tab/coma-delimited, with a header and has at least 3 columns.
The three designated columns must be of float value.
If "-i ..." is given, it is regarded as one of the input files (plus the ones in trailing arguments). 

Examples:	
    %s -i /tmp/Contig315_StKitts_vs_Nevis.tsv --xColumnHeader=StKitts
        --whichColumnHeader=Nevis
        -s 1.0 -o /tmp/Contig315_StKitts_vs_Nevis.2D.png
    
"""

import sys, os, math
__doc__ = __doc__%(sys.argv[0])

import matplotlib; matplotlib.use("Agg")	#to disable pop-up requirement
import logging
from palos import ProcessOptions, PassingData
from palos.plot import yh_matplotlib
import numpy, random
from palos.plot.AbstractPlot import AbstractPlot

class Draw2DHistogramOfMatrix(AbstractPlot):
    __doc__ = __doc__
    option_default_dict = AbstractPlot.option_default_dict.copy()
    for option_key in [('ylim_type', 1, int),]:
        option_default_dict.pop(option_key)
    option_default_dict.update({
        ('zColumnHeader', 0, ): ["", 'z', 1, 'index of the column to be z-axis'],\
        ('logZ', 0, int): [0, '', 1, 'value 0: nothing; 1: log(), 2: -log()'],\
        })
    #('columnForX', 1, ): ["", 'x', 1, 'index of the column to be x-axis'],\
    #('columnForY', 1, ): ["", 'y', 1, 'index of the column to be y-axis'],\
    #('whichColumnPlotLabel', 0, ), ('xColumnHeader', 1, ), ('xColumnPlotLabel', 0, )
    def __init__(self, inputFnameLs, **keywords):
        """
        """
        AbstractPlot.__init__(self, inputFnameLs=inputFnameLs, **keywords)
        
        self.x_value_ls = []
        self.y_value_ls = []
        self.z_value_ls = []
    
    def initiatePassingData(self, ):
        """
        this function gets called in the beginning of each fileWalker() (for each inputFname).
        """
        pdata = PassingData(x_ls = [], y_ls = [], z_ls=[], 
            invariantPData=self.invariantPData)
        return pdata
    
    def processRow(self, row=None, pdata=None):
        """
        """
        col_name2index = getattr(pdata, 'col_name2index', None)
        x_ls = getattr(pdata, 'x_ls', None)
        y_ls = getattr(pdata, 'y_ls', None)
        z_ls = getattr(pdata, 'z_ls', None)
        
        if col_name2index and x_ls is not None and y_ls is not None:
            x_index = col_name2index.get(self.xColumnHeader, None)
            if self.whichColumnHeader:
                y_index = col_name2index.get(self.whichColumnHeader, None)
            else:
                y_index = self.whichColumn
            
            if self.zColumnHeader:
                z_index = col_name2index.get(self.zColumnHeader)
            else:
                z_index = None
            
            xValue = row[x_index]
            yValue = row[y_index]
            if z_index is not None:
                zValue = row[z_index]
            else:
                zValue = None
            
            if yValue not in self.missingDataNotation and xValue not in self.missingDataNotation:
                xValue = self.processValue(value=xValue, processType=self.logX, \
                    valueForNonPositiveValue=self.valueForNonPositiveYValue)
                yValue = self.processValue(yValue, processType=self.logY, \
                    valueForNonPositiveValue=self.valueForNonPositiveYValue)
                if zValue is not None:
                    zValue = self.processValue(zValue, processType=self.logZ, \
                        valueForNonPositiveValue=self.valueForNonPositiveYValue)
                else:
                    zValue = 1	#just a counter
                x_ls.append(xValue)
                y_ls.append(yValue)
                z_ls.append(zValue)
    
    def afterFileFunction(self, x_ls=None, y_ls=None, pdata=None):
        """
        get called by the end of fileWalker() for each inputFname.
        this is a reducer for each file.
        """
        #pylab.plot(x_ls, y_ls, self.formatString)
        self.invariantPData.x_ls.extend(pdata.x_ls)
        self.invariantPData.y_ls.extend(pdata.y_ls)
        self.invariantPData.z_ls.extend(pdata.z_ls)
        
    def saveFigure(self, invariantPData=None):
        """
        executed after all files have been walked through
        """
        
        if self.outputFnamePrefix:
            pngOutputFname = '%s.png'%self.outputFnamePrefix
            svgOutputFname = '%s.svg'%self.outputFnamePrefix
        elif self.outputFname:
            pngOutputFname = self.outputFname
            svgOutputFname = '%s.svg'%(self.outputFname[:-4])
        else:
            logging.error(f"Could not get outputFnamePrefix from "
                f"self.outputFnamePrefix {self.outputFnamePrefix} or "
                f"self.outputFname {self.outputFname}.")
            sys.exit(1)
        
        if self.zColumnHeader:	#take median of the selected Z-value
            colorBarLabelForZ = "median of %s"%(self.zColumnHeader)
            reduce_C_function = numpy.median
        else:	#plain 2-D histogram, Z is the logSum count
            colorBarLabelForZ = 'log10(no of data)'
            reduce_C_function = yh_matplotlib.logSum
        
        #if uniformly distributed over 2D, each hexagon has ~10 points.
        # the maximum gridsize is 30.
        gridsize = min(30, int(math.sqrt(len(invariantPData.x_ls)/10.0)))
        
        xscale = "linear"
        
        yh_matplotlib.drawHexbin(invariantPData.x_ls, invariantPData.y_ls,
            invariantPData.z_ls,
            fig_fname=pngOutputFname, gridsize=gridsize, title=self.handleTitle(),
            xlabel=self.handleXLabel(), ylabel=self.handleYLabel(),
            colorBarLabel=colorBarLabelForZ,
            reduce_C_function=reduce_C_function, dpi=self.figureDPI,
            mincnt=5, marginals=False, xscale=xscale)
            #at least 5 data in one hexagon
        """
        outputFname = '%s_%s_vs_%s_loci_count.png'%(outputFnamePrefix,
            self.columnForX, self.columnForY,)
        yh_matplotlib.drawHexbin(self.x_value_ls, self.y_value_ls, loci_count_ls,
            fig_fname=outputFname, gridsize=gridsize, title=title,
            xlabel=self.columnForX, ylabel=self.columnForY,
            colorBarLabel=colorBarLabelForLociCount,
            reduce_C_function=yh_matplotlib.logSum, dpi=self.figureDPI,
            xscale=xscale)
        """
        
    
    
if __name__ == '__main__':
    main_class = Draw2DHistogramOfMatrix
    po = ProcessOptions(sys.argv, main_class.option_default_dict,
        error_doc=main_class.__doc__)
    instance = main_class(po.arguments, **po.long_option2value)
    instance.run()
