#!/usr/bin/env python3
"""
2013.07.12 program that removes any redundant SNPs. By redundant, (chromosome, position) appears >1.

Examples:
	%s -i tmp/testVCF/CAEY.sorted.vcf.gz  -o /tmp/VCF_CAEY.unique.vcf

"""
import sys, os, math
__doc__ = __doc__%(sys.argv[0])


from Bio.Seq import Seq
from palos import ProcessOptions, PassingData
from palos.ngs.io.VCFFile import VCFFile
from . CalculateSameSiteConcordanceInVCF import CalculateSameSiteConcordanceInVCF

ParentClass = CalculateSameSiteConcordanceInVCF
class RemoveRedundantLociFromVCF(ParentClass):
	__doc__ = __doc__
	option_default_dict = ParentClass.option_default_dict.copy()
	option_default_dict.update({
						})
	def __init__(self, inputFnameLs=None, **keywords):
		"""
		"""
		ParentClass.__init__(self, inputFnameLs=inputFnameLs, **keywords)
	
	def run(self):
		if self.debug:
			import pdb
			pdb.set_trace()
		
		
		
		outputDir = os.path.split(self.outputFname)[0]
		if outputDir and not os.path.isdir(outputDir):
			os.makedirs(outputDir)
		
		snp_pos2count =self.readInSNPID2GenotypeVectorLs(self.inputFname, returnType=2).snp_pos2returnData
		
		
		reader = VCFFile(inputFname=self.inputFname)
		
		writer = VCFFile(outputFname=self.outputFname, mode='w')
		writer.metaInfoLs = reader.metaInfoLs
		writer.header = reader.header
		writer.writeMetaAndHeader()
		
		
		counter = 0
		real_counter = 0
		for vcfRecord in reader:	#assuming input VCF is sorted
			counter += 1
			key = (vcfRecord.chromosome, vcfRecord.position)
			frequency = snp_pos2count.get(key)
			if frequency==1:
				writer.writeVCFRecord(vcfRecord)
				real_counter += 1
			
		reader.close()
		writer.close()
		if counter >0:
			fraction=real_counter/float(counter)
		else:
			fraction = 0
		sys.stderr.write("%s (out of %s, %s) snps are unique.\n"%(real_counter, counter, fraction))

if __name__ == '__main__':
	main_class = RemoveRedundantLociFromVCF
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()