#!/usr/bin/env python3
""""
2012.8.28 test class for TaxonomyDB.py

Examples:
	%s -u yh -p secret -d vervetdb -k taxonomy

"""
import sys, os, math
__doc__ = __doc__%(sys.argv[0])
import unittest, os, sys, getopt, csv
from palos.db.TaxonomyDB import TaxonomyDB

class TestTaxonomyDB(unittest.TestCase, TaxonomyDB):
	__doc__ = __doc__
	option_default_dict = TaxonomyDB.option_default_dict.copy()
	def __init__(self, testname=None, **keywords):
		"""
		2012.8.28
		"""
		TaxonomyDB.__init__(self, **keywords)
		
		super(TestTaxonomyDB, self).__init__(testname)
		#self.hs = helpspot.HelpSpot(path, user, pword)
		
	def setUp(self):
		print
	
	def test_returnScientificNameGivenTaxID(self):
		self.setup()
		taxID = 60711
		
		#1st just db query
		scientificName = self.getScientificNameGivenTaxID(taxID)
		sys.stderr.write("Scientific name for %s is %s.\n"%(taxID, scientificName))
		self.assertEqual(scientificName, "Chlorocebus sabaeus")
		
		#2nd cached version
		scientificName = self.returnScientificNameGivenTaxID(taxID)
		sys.stderr.write("Scientific name for %s is %s.\n"%(taxID, scientificName))
		self.assertEqual(scientificName, "Chlorocebus sabaeus")

	
	def test_returnTaxIDGivenScientificName(self):
		self.setup()
		scientificName = "Chlorocebus sabaeus"
		#1st just db query
		taxID = self.getTaxIDGivenScientificName(scientificName)
		sys.stderr.write("Taxonomy ID for %s is %s.\n"%(scientificName, taxID))
		self.assertEqual(taxID, 60711)
		
		#2nd cached version
		taxID = self.returnTaxIDGivenScientificName(scientificName)
		sys.stderr.write("Taxonomy ID for %s is %s.\n"%(scientificName, taxID))
		
		self.assertEqual(taxID, 60711)

		
	
if __name__ == '__main__':
	import sys, os, math
	
	TestCaseDict = {1:TestTaxonomyDB,}
	type = 1
	
	from palos import ProcessOptions
	main_class = TestCaseDict[type]
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	suite = unittest.TestSuite()
	suite.addTest(TestCaseDict[type](testname='test_returnScientificNameGivenTaxID', **po.long_option2value))
	suite.addTest(TestCaseDict[type](testname='test_returnTaxIDGivenScientificName', **po.long_option2value))
	#suite.addTest(unittest.makeSuite(TestCaseDict[type]))
	unittest.TextTestRunner(verbosity=2).run(suite)
	#unittest.TextTestRunner().run(suite)
	"""

	instance = main_class(**po.long_option2value)
	instance.setup()
	instance.run()
	"""
