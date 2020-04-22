#!/usr/bin/env python3
"""
Description:
	20130731
	abstract walker class for genome files
"""

import sys, os, math
import csv, random, numpy
from palos import ProcessOptions, getListOutOfStr, PassingData, utils
from palos.db import GenomeDB
from palos.io.AbstractMatrixFileWalker import AbstractMatrixFileWalker

ParentClass = AbstractMatrixFileWalker

class AbstractGenomeFileWalker(ParentClass):
	__doc__ = __doc__
	option_default_dict = ParentClass.option_default_dict.copy()
	option_default_dict.update(ParentClass.genome_db_option_dict.copy())
	#option_default_dict.update(ParentClass.db_option_dict.copy())
	
	genome_option_dict = {
		('chromosomeHeader', 1, ): ['chromosome', '', 1, 'header of the column that designates chromosome info' ],\
		('positionHeader', 1, ): ['position', '', 1, 'header of the column that designates position of each locus on chromosome' ],\
		('tax_id', 0, int): [3702, '', 1, 'taxonomy ID of the organism from which '
			'to retrieve the chromosome info', ],\
		('sequence_type_id', 0, int): [1, '', 1, "sequence_type_id (annot_assembly) of the chromosomes, "
			"to retrieve the chromosome info. 1: assembledChromosome, 9: Scaffolds", ],\
		('chrOrder', 0, int): [1, '', 1, 'how to order chromosomes. '
			'1: column genome_order in db; 2: by chromosome size, descending', ],\
		}
	option_default_dict.update(genome_option_dict)
	option_default_dict.update({
					
					
					})
	def __init__(self, inputFnameLs=None, **keywords):
		"""
		"""
		ParentClass.__init__(self, inputFnameLs=inputFnameLs, **keywords)	#self.connectDB() called within its __init__()
	
	def _loadGenomeStructureFromDB(self):
		"""
		2013.07.31
			so that both AbstractGenomeFileWalker and PlotGenomeWideData could call it
		"""
		db_genome = GenomeDB.GenomeDatabase(drivername=self.genome_drivername, username=self.genome_db_user,
						db_passwd=self.genome_db_passwd, hostname=self.genome_hostname, dbname=self.genome_dbname, \
						schema=self.genome_schema)
		db_genome.setup(create_tables=False)
		self.db_genome = db_genome
		
		#chrOrder=1 is to order chromosomes alphabetically
		self.oneGenomeData = db_genome.getOneGenomeData(tax_id=self.tax_id, chr_gap=0, chrOrder=self.chrOrder, \
													sequence_type_id=self.sequence_type_id)
		self.chr_id2cumu_start = self.oneGenomeData.chr_id2cumu_start
		#2013.2.18 get centromere locations
		self.chr_id2centromere = self.oneGenomeData.chr_id2centromere
		
	
	def setup(self, **keywords):
		"""
		2013.07.31
			construct an RBTree dictionary map between windows and their data
		"""
		ParentClass.setup(self, **keywords)
		self._loadGenomeStructureFromDB()
	

if __name__ == '__main__':
	main_class = AbstractGenomeFileWalker
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(po.arguments, **po.long_option2value)
	instance.run()