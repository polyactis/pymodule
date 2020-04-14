#!/usr/bin/env python3
import sys, os, math
import csv, re, copy
sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))
sys.path.insert(0, os.path.join(os.path.expanduser('~/src')))
from pymodule import ProcessOptions, PassingData, utils
from Sunset.mapper.AbstractAccuMapper import AbstractAccuMapper as ParentClass
from pymodule.db import SunsetDB as DBClass

class RegisterIndividualSequence2DB(ParentClass):
	"""
Description:
	program to register individual_sequence. It will check if input file is empty or not. If empty, exit non-0.
		If not, add db isq entry, output isq-id and isq-raw-file-id, exit 0.
Examples:
	%s  --drivername postgresql --hostname 149.142.212.14 --dbname pmdb --schema public --db_user yh 
		--individual_id 3185 -i 3185_gerald_D14CGACXX_7_GCCAAT.fastq.gz -o isq_id_isqf_raw.id.txt
		--sequence_format fastq --commit
		--original_sequence_library gerald_D14CGACXX_7_GCCAAT 
		--original_sequence_mate_id 1
		--original_sequence_filepath gerald_D14CGACXX_7_GCCAAT.bam
	"""
	__doc__ = __doc__%(sys.argv[0])
	option_default_dict = copy.deepcopy(ParentClass.option_default_dict)
	option_default_dict.pop(('outputFnamePrefix', 0, ))
	option_default_dict.update({
						("individual_id", 1, int): [None, '', 1, ''],\
						('tissue_name', 0,): [None, '', 1, 'table tissue.short_name', ],\
						('tissue_id', 0, int): [None, '', 1, 'table tissue.id'],\
						("sequencer_id", 0, int): [None, '', 1, 'isq.sequencer_id'],\
						("sequence_format", 0, ): ["fastq", '', 1, 'fasta, fastq, etc.'],\
						("sequence_type_name", 0, ): ["", '', 1, 'PairedEnd or SingleRead or others in table sequence_type.'],\
						("original_sequence_filepath", 0, ): ["", '', 1, 'path to the original sequence file'],\
						('original_sequence_mate_id', 0, int):[None, '', 1, 'mate_id of original_sequence, only valid if original_sequence contains only one end'],\
						('original_sequence_md5sum', 0, ): ["", '', 1, 'The md5sum for the original_sequence_filepath, not the input file '],\
						('original_sequence_library', 0,): [None, '', 1, 'library name for the original_sequence', ],\
						('copy_original_file', 0, int):[0, '', 1, 'toggle to copy original_sequence_filepath into db-affliated storage'],\
						("quality_score_format", 1, ): ["Standard", '', 1, 'isq'],\
						('coverage', 0, int):[None, '', 1, 'isq'],\
						('filtered', 0, int):[0, '', 1, 'isq'],\
						('parent_individual_sequence_id', 0, int):[None, '', 1, 'isq'],\
						('read_count', 0, int):[None, '', 1, 'isq'],\
						('no_of_chromosomes', 0, int):[None, '', 1, 'isq'],\
						('sequence_batch_id', 0, int):[None, '', 1, 'isq'],\
						('version', 0, int):[None, '', 1, 'isq.version'],\
						('is_contaminated', 0, int):[0, '', 1, 'isq.is_contaminated'],\
						('outdated_index', 0, int):[0, '', 1, 'isq'],\
						('comment', 0, ):[None, '', 1, 'isq.comment'],\
						('data_dir', 0, ):[None, '', 1, 'used to create folder data_dir/isq.path'],\
						})

	def __init__(self,  **keywords):
		"""
		"""
		#connectDB(), and setup srcFilenameLs and dstFilenameLs
		ParentClass.__init__(self, inputFnameLs=None, **keywords)
	
	def run(self):
		"""
		"""
		if self.debug:
			import pdb
			pdb.set_trace()
		#check if inputFname is empty
		inputFile = utils.openGzipFile(self.inputFname)
		char_counter = 0
		for line in inputFile:
			#only need one line
			char_counter += len(line)
			break
		inputFile.close()
		sys.stderr.write("First line character count of %s: %s.\n"%(self.inputFname, char_counter))
		if char_counter==0:
			sys.stderr.write("ERROR: exit due to empty file.\n")
			sys.exit(2)
			
		db_main = self.db_main
		session = db_main.session
		session.begin()
		
		if self.data_dir:
			data_dir = self.data_dir
		else:
			data_dir = db_main.data_dir
		#uuid is for sequence only, add as isq.comment
		individual_sequence = db_main.getIndividualSequence(individual_id=self.individual_id, \
							sequencer_id=self.sequencer_id,\
							sequence_type_name=self.sequence_type_name, \
						sequence_format=self.sequence_format, path_to_original_sequence=self.original_sequence_filepath, \
						copy_original_file=self.copy_original_file,\
						tissue_name=self.tissue_name, tissue_id=self.tissue_id, \
						coverage=self.coverage,\
						subFolder=None, quality_score_format=self.quality_score_format, filtered=self.filtered,\
						parent_individual_sequence_id=self.parent_individual_sequence_id,\
						read_count=self.read_count, no_of_chromosomes=self.no_of_chromosomes, \
						sequence_batch_id=self.sequence_batch_id, version=self.version, data_dir=data_dir,\
						is_contaminated=self.is_contaminated, outdated_index=self.outdated_index, comment=self.comment)
		file_raw_db_entry = None
		if self.original_sequence_filepath:
			file_raw_db_entry = db_main.registerOriginalSequenceFileToDB(self.original_sequence_filepath, 
																	  library=self.original_sequence_library, \
										individual_sequence_id=individual_sequence.id, mate_id=self.original_sequence_mate_id, \
										md5sum=self.original_sequence_md5sum)
		
		#output isq_id to outputFname
		outputDir = os.path.join(data_dir, individual_sequence.path)
		if not os.path.isdir(outputDir):
			os.makedirs(outputDir)
		if self.outputFname:
			outf = open(self.outputFname, 'w')
			outf.write("individual_sequence_id: %s\n"%(individual_sequence.id))
			if file_raw_db_entry:
				outf.write("individual_sequence_file_raw_id: %s\n"%(file_raw_db_entry.id))
			outf.write("outputDir: %s\n"%(outputDir))
			outf.write("relativeOutputDir: %s\n"%(individual_sequence.path))
			outf.close()
		if self.commit:
			session.commit()
		else:
			self.sessionRollback(session)
		
if __name__ == '__main__':
	main_class = RegisterIndividualSequence2DB
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()