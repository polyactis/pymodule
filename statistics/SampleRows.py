#!/usr/bin/env python3
"""
2012.10.15
	Sample rows from a matrix-like file.
	If "-i ..." is given, it is regarded as one of the input files (plus the ones in trailing arguments). 

Examples:
	%s -i /tmp/VRCSamples.tsv --whichColumnHeader=sampleID
		-s 1.0  --sampleSize 5
		--plinkIBDCheckOutputFname PlinkIBDCheck/ibdCheckIBDCheck/LDPrunedMerged_ibdCheck.tsv
		-o /tmp/VRCSamples_sampled.tsv

"""

import sys
__doc__ = __doc__%(sys.argv[0])

import csv
import random
from palos import ProcessOptions, PassingData, figureOutDelimiter
from palos.polymorphism import SNP
from palos.io.AbstractMatrixFileWalker import AbstractMatrixFileWalker

class SampleRows(AbstractMatrixFileWalker):
	__doc__ = __doc__
	option_default_dict = AbstractMatrixFileWalker.option_default_dict.copy()
	option_default_dict.update({
		('sampleSize', 1, int): [None, '', 1, 'number of samples (rows) to be sampled from input, \n\
uniformly if plinkIBDCheckOutputFname is not given.'],\
		('plinkIBDCheckOutputFname', 0, ): [None, '', 1, 'file that contains IBD check result, PI_HAT=relatedness.\n\
at least 3-columns with header: IID1, IID2, PI_HAT. IID1 and IID2 should match the whichColumn (whichColumnHeader) of inputFname.\n\
The sampling will try to avoid sampling close pairs, PI_HAT(i,j)<=maxIBDSharing'],\
		('maxIBDSharing', 1, float): [0.1, '', 1, 'This argument caps the maximum IBD sharing among any pair within the sampled.'],\
		})
	def __init__(self, inputFnameLs=None, **keywords):
		"""
		"""
		AbstractMatrixFileWalker.__init__(self, inputFnameLs=inputFnameLs, **keywords)	#self.connectDB() called within its __init__()
		
	def setup(self, **keywords):
		"""
		2012.10.15
			run before anything is run
		"""
		AbstractMatrixFileWalker.setup(self, **keywords)
		#read in the IBD check result
		if self.plinkIBDCheckOutputFname:
			ibdData = SNP.readAdjacencyListDataIntoMatrix(inputFname=self.plinkIBDCheckOutputFname, rowIDHeader="IID1", colIDHeader="IID2", \
								rowIDIndex=None, colIDIndex=None, \
								dataHeader="PI_HAT", dataIndex=None, hasHeader=True)
		else:
			ibdData = None
		self.ibdData = ibdData
		self.data_matrix = []	#data structure to store all rows during fileWalker()
	
	def processRow(self, row=None, pdata=None):
		"""
		2012.10.7
		"""
		returnValue = 1
		self.data_matrix.append(row)
		self.col_name2index = getattr(pdata, 'col_name2index', None)
		
		return returnValue
	
	def reduce(self, **keywords):
		"""
		2012.10.15
			run after all files have been walked through
		"""
		#sample the data
		noOfTotalRows = len(self.data_matrix)
		real_counter = 0
		if self.sampleSize<noOfTotalRows:
			if self.ibdData:
				col_name2index = getattr(self, 'col_name2index', None)
				sampledRowIndexList = []
				if col_name2index:
					sampledRowIndexSet = set()
					while len(sampledRowIndexSet)<self.sampleSize:
						randomIndex = random.randint(0,noOfTotalRows-1)
						if randomIndex in sampledRowIndexSet:	#must not be in there already
							continue
						row = self.data_matrix[randomIndex]
						if self.whichColumnHeader:
							whichColumn = col_name2index.get(self.whichColumnHeader, None)
						elif self.whichColumn:
							whichColumn = self.whichColumn
						else:
							whichColumn = None
						if whichColumn is not None:
							rowID = row[whichColumn]
							includeInTheSampling = True
							for sampledRowIndex in sampledRowIndexSet:	#not too close to anyone previously sampled
								previousSampledRowID = self.data_matrix[sampledRowIndex][whichColumn]
								#getting the relatedness
								relatedness = self.ibdData.getCellDataGivenRowColID(rowID, previousSampledRowID)
								if relatedness>=self.maxIBDSharing:
									includeInTheSampling = False
							if includeInTheSampling:
								sampledRowIndexSet.add(randomIndex)
					sampledRowIndexList = list(sampledRowIndexSet)
			else:
				sampledRowIndexList = random.sample(range(noOfTotalRows), self.sampleSize)
		else:	#take all
			sampledRowIndexList = range(noOfTotalRows)
		for i in sampledRowIndexList:
			row = self.data_matrix[i]
			self.invariantPData.writer.writerow(row)
			real_counter +=1 
		
		fraction = float(real_counter)/float(noOfTotalRows)
		sys.stderr.write("%s/%s (%.3f) selected.\n"%(real_counter, noOfTotalRows, fraction))
		
		#close the self.invariantPData.writer
		AbstractMatrixFileWalker.reduce(self, **keywords)

if __name__ == '__main__':
	main_class = SampleRows
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(po.arguments, **po.long_option2value)
	instance.run()
