#!/usr/bin/env python3
"""
2013.05.03 a child class of BeagleGenotypeFile. used to describe Beagle likelihood file format, which looks like:
    marker alleleA alleleB 1000_709_1996093_GA_vs_524 1000_709_1996093_GA_vs_524 1000_709_1996093_GA_vs_524 1001_710_1995025_GA_vs_524 1001_710_1995025_GA_vs_524 1001_710_1995025_GA_vs_524 1002_711_2001039_GA_vs_524
    Contig791:1086 C A 1 0 0 0.9997 0.0003 0 0
    Contig791:1649 G C 1 0 0 1 0 0 0
    Contig791:4084 A C 1 0 0 1 0 0 0
    Contig791:4118 A G 1 0 0 1 0 0 0
    Contig791:4143 C A 1 0 0 1 0 0 0
    Contig791:4168 G C 1 0 0 0.9999 0.0001 0 0
    Contig791:4203 C G 1 0 0 1 0 0 0

Example:

    reader = MatrixFile(path='/tmp/input.txt', mode='r')
    reader = MatrixFile('/tmp/input.txt', mode='r')
    reader.constructColName2IndexFromHeader()
    for row in reader:
        row[reader.getColName2IndexFromHeader('KID')]
    
    inf = utils.openGzipFile(path, mode='r')
    reader = MatrixFile(file_handle=inf)
    
    #2013.2.1 writing
    writer = MatrixFile('/tmp/output.txt', mode='w', delimiter='\t')
    writer.writeHeader(...)
    writer.writerow(row)
    writer.close()	

"""

import sys
import os
import copy
from palos.ProcessOptions import  ProcessOptions
from palos.utils import PassingData
from . BeagleGenotypeFile import BeagleGenotypeFile

ParentClass = BeagleGenotypeFile
class BeagleLikelihoodFile(ParentClass):
    __doc__ = __doc__
    option_default_dict = copy.deepcopy(ParentClass.option_default_dict)
    def __init__(self, path=None, **keywords):
        ParentClass.__init__(self, path=path, **keywords)
    
    def constructColName2IndexFromHeader(self):
        """
        First three column is for marker, alleleA, alleleB.
        one sample (diploid) occupies three columns.
            
marker alleleA alleleB 1000_709_1996093_GA_vs_524 1000_709_1996093_GA_vs_524 1000_709_1996093_GA_vs_524 1001_710_1995025_GA_vs_524 1001_710_1995025_GA_vs_524 1001_710_1995025_GA_vs_524 1002_711_2001039_GA_vs_524
Contig791:1086 C A 1 0 0 0.9997 0.0003 0 0
Contig791:1649 G C 1 0 0 1 0 0 0
Contig791:4084 A C 1 0 0 1 0 0 0
        """
        self.header = next(self).genotypeLikelihoodList
        self.col_name2index = {}
        for i in range(len(self.header)):
            sampleID = self.header[i]
            if sampleID not in self.col_name2index:
                self.col_name2index[sampleID] = []
            self.col_name2index[sampleID].append(i)
            #the index corresponds to genotypeLikelihoodList in next().
            #so it starts from 0
        return self.col_name2index
    
    def getLikelihoodListOfOneGenotypeOneSample(self, oneLocus=None, sampleID=None):
        """
        oneLocus is output of next()
        """
        try:
            sampleStartIndex = self.getColIndexGivenColHeader(sampleID)[0]
        except:
            sys.stderr.write('Except type: %s\n'%repr(sys.exc_info()))
            import traceback
            traceback.print_exc()
            sys.stderr.write("  sampleID %s not in beagleLikelihoodFile.\n"%(sampleID))
            raise
        tripleLikelihood = oneLocus.genotypeLikelihoodList[sampleStartIndex:sampleStartIndex+3]
        return tripleLikelihood
    
    def __next__(self):
        try:
            row = next(self.csvFile)
        except:
            raise StopIteration
        if not self.isRealCSV:
            row = row.strip().split()
        markerID, alleleA, alleleB = row[0:3]
        return PassingData(markerID=markerID, alleleA=alleleA, alleleB=alleleB,
            genotypeLikelihoodList=row[3:])

if __name__ == '__main__':
    main_class = BeagleLikelihoodFile
    po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
    instance = main_class(**po.long_option2value)
    instance.run()