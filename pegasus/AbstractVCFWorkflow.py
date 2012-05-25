#!/usr/bin/env python
"""
2012.1.17
	a common class for pegasus workflows that work on VCF variant files
"""
import sys, os, math

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import subprocess, cStringIO
import VervetDB
from pymodule import ProcessOptions, getListOutOfStr, PassingData, yh_pegasus, NextGenSeq
from Pegasus.DAX3 import *
from AbstractNGSWorkflow import AbstractNGSWorkflow

class AbstractVCFWorkflow(AbstractNGSWorkflow):
	__doc__ = __doc__
	option_default_dict = AbstractNGSWorkflow.option_default_dict.copy()
	option_default_dict.update({
						('minDepth', 0, float): [0, 'm', 1, 'minimum depth for a call to regarded as non-missing', ],\
						})
	def __init__(self,  **keywords):
		"""
		2012.1.17
		"""
		AbstractNGSWorkflow.__init__(self, **keywords)
		
	def registerExecutables(self, workflow):
		"""
		"""
		AbstractNGSWorkflow.registerExecutables(self, workflow)
	
	def registerCommonExecutables(self, workflow):
		"""
		"""
		AbstractNGSWorkflow.registerCommonExecutables(self, workflow)
	
	@classmethod
	def registerAllInputFiles(cls, workflow, inputDir, input_site_handler=None, checkEmptyVCFByReading=False, pegasusFolderName=''):
		"""
		2012.5.9
			register the tbi file if it exists
		2012.3.1
			moved from CalculateTrioInconsistencyPipeline.py
		2012.1.9
			the returning data structure is changed to conform to some standard used across several functions
		2011-9-29
			vcf files only
		"""
		sys.stderr.write("Registering input files from %s ..."%(inputDir))
		returnData = PassingData(jobDataLs = [])
		fnameLs = os.listdir(inputDir)
		counter = 0
		real_counter = 0
		previous_reported_real_counter = ''
		for fname in fnameLs:
			counter += 1
			inputFname = os.path.join(inputDir, fname)
			if NextGenSeq.isFileNameVCF(fname, includeIndelVCF=False) and \
					not NextGenSeq.isVCFFileEmpty(inputFname, checkContent=checkEmptyVCFByReading):
				real_counter += 1
				inputBaseFname = os.path.basename(inputFname)
				inputF = File(os.path.join(pegasusFolderName, inputBaseFname))
				inputF.addPFN(PFN("file://" + inputFname, input_site_handler))
				inputF.abspath = inputFname
				workflow.addFile(inputF)
				
				tbi_F_absPath = "%s.tbi"%inputFname
				if os.path.isfile(tbi_F_absPath):	#it exists
					tbi_F = File(os.path.join(pegasusFolderName, "%s.tbi"%inputBaseFname))
					tbi_F.addPFN(PFN("file://" + tbi_F_absPath, input_site_handler))
					tbi_F.abspath = tbi_F_absPath
					workflow.addFile(tbi_F)
				else:
					tbi_F = None
				
				returnData.jobDataLs.append(PassingData(vcfFile=inputF, jobLs=[], tbi_F=tbi_F))
				if real_counter%200==0:
					sys.stderr.write("%s%s"%('\x08'*len(previous_reported_real_counter), real_counter))
					previous_reported_real_counter = repr(real_counter)
		sys.stderr.write("%s non-empty out of %s files.\n"%(len(returnData.jobDataLs), counter))
		return returnData