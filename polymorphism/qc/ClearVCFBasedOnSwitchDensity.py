#!/usr/bin/env python3
"""
2013.07.12 program that clears the whole liftover VCF if its switch point density above the maximum.

Examples:
	%s -i folderReduceLiftOverVCF/CAEY.sorted.vcf.gz -o CAEY.sameSite.concordance.tsv
	
"""
import sys, os
__doc__ = __doc__%(sys.argv[0])

from palos import ProcessOptions, PassingData
from palos.io.MatrixFile import MatrixFile
from palos.ngs.io.VCFFile import VCFFile
from palos.mapper.AbstractVCFMapper import AbstractVCFMapper

ParentClass = AbstractVCFMapper
class ClearVCFBasedOnSwitchDensity(ParentClass):
	__doc__ = __doc__
	option_default_dict = ParentClass.option_default_dict.copy()
	option_default_dict.update({
					('switchPointFname', 0, ): ['', '', 1, 'file that has switch points information, output of FindSNPPositionOnNewRefFromFlanking???Output.py ', ],\
					('maxSwitchDensity', 0, float): [0.01, '', 1, 'Maximum switch density (#switches/#loci) for one interval to be included in final variants', ],\
				})
	def __init__(self, inputFnameLs=None, **keywords):
		"""
		"""
		ParentClass.__init__(self, inputFnameLs=inputFnameLs, **keywords)
	
	
	def readInSwitchDensity(self, inputFname=None):
		"""
		2013.07.11
		"""
		sys.stderr.write("Reading in switch density from %s ..."%(inputFname))
		
		reader = MatrixFile(inputFname=inputFname)
		reader.constructColName2IndexFromHeader()
		
		noOfSwitchesPerLocusIndex = reader.getColIndexGivenColHeader("noOfSwitchesPerLocus")
		
		counter = 0
		real_counter = 0
		switchDensity = 0
		for row in reader:
			switchDensity = float(row[noOfSwitchesPerLocusIndex])
			counter += 1
			break
		del reader
		return PassingData(switchDensity=switchDensity)
	
	def run(self):
		if self.debug:
			import pdb
			pdb.set_trace()
		
		outputDir = os.path.split(self.outputFname)[0]
		if outputDir and not os.path.isdir(outputDir):
			os.makedirs(outputDir)
		
		switchDensity = self.readInSwitchDensity(inputFname=self.switchPointFname).switchDensity
		
		reader = VCFFile(inputFname=self.inputFname)
		
		writer = VCFFile(outputFname=self.outputFname, openMode='w')
		writer.metaInfoLs = reader.metaInfoLs
		writer.header = reader.header
		writer.writeMetaAndHeader()
		counter = 0
		real_counter = 0
					
		if switchDensity<=self.maxSwitchDensity:

			for vcfRecord in reader:	#assuming input VCF is sorted
				counter += 1
				real_counter += 1
				writer.writeVCFRecord(vcfRecord)
			
		reader.close()
		writer.close()
		sys.stderr.write("%s (out of %s) records outputted.\n"%(real_counter, counter))
		

if __name__ == '__main__':
	main_class = ClearVCFBasedOnSwitchDensity
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()