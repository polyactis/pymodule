import os, sys, csv, re

def dict_map(dict, ls, type=1):
	"""
	2008-04-03 copied from annot.codense.common
	10-13-05
		add type 2 to return item itself if mapping is not available
	2006-12-21
		add type 3 to extract a smaller map
	2007-05-14
		bug, "if value" could miss 0
	"""
	if type==3:
		new_list = {}	#it's a dictionary
		for item in ls:
			value = dict.get(item)
			if value is not None:
				new_list[item] = value
	else:
		new_list = []
		for item in ls:
			value = dict.get(item)
			if value is not None:
				new_list.append(value)
			elif type==2:
				new_list.append(item)
	
	return new_list

class PassingData(object):
	"""
	05/09/08
		a class to hold any data structure
	"""
	def __init__(self, **keywords):
		"""
		2008-5-12
			add keyword handling
		"""
		for argument_key, argument_value in keywords.iteritems():
			setattr(self, argument_key, argument_value)
	
	def __str__(self):
		"""
		2010-6-17
			a string-formatting function
		"""
		return_ls = []
		for attribute_name in dir(self):
			if attribute_name.find('__')==0:	#ignore the 
				continue
			return_ls.append("%s = %s"%(attribute_name, getattr(self, attribute_name, None)))
		
		return ", ".join(return_ls)
		
def importNumericArray():
	"""
	2008-07-09
		numarray doesn't have int128
	2008-05-18
		give same numpy types (int, int8 ...) to other numeric modules
	2008-05-18
		add "import array as num"
		should put ImportError in except. but whatever
	2008-05-11
		import whatever available array module
	"""
	try:
		import numpy as num
	except:
		numpy_type2other_ls = ['int', 'int8', 'int16', 'int32', 'int64']
		try:
			import numarray as num
		except:
			import Numeric as num
		for numpy_type in numpy_type2other_ls:	#make sure it has same type names
			numpy_type_in_other = numpy_type[0].upper() + numpy_type[1:]
			setattr(num, numpy_type, getattr(num, numpy_type_in_other))
	return num

def figureOutDelimiter(input_fname, report=0, delimiter_choice_ls = ['\t', ','], use_sniff=False):
	"""
	2008-01-08
		don't use cs.sniff unless the user specifies it. sniff gives you unexpected delimiter when it's a single-column.
	2008-08-28
		nothing weird on hpc-cmb. it's a bug in other code.
		back to 'return None' if input_fname escapes all condition checking.
	2008-08-28
		try 'open(input_fname)' anyway if input_fname escapes all condition checking.
		something weird happened during a mpi job on hpc-cmb. the file is there. but escape the first condition.
	2008-05-25
		now 3 possible types of input_fname
		1. a file name (path)
		2. input_fname is a file object
		3. input_fname is input data, string
		
		for a file object or input file name:
		it could be binary file which doesn't have readline(). have to use this dumb approach due to '\n' might mess up sniff()
	2008-05-21
		csv.Sniffer is handy, use it figure out csv.Sniffer instead.
	2008-05-12
		try tab first
	"""
	if report:
		import sys
		sys.stderr.write("Figuring out delimiter for %s ..."%input_fname)
	cs = csv.Sniffer()
	if isinstance(input_fname, str) and os.path.isfile(input_fname):
		inf = open(input_fname)
	elif isinstance(input_fname, file):	#could be a file object
		inf = input_fname
	elif isinstance(input_fname, str) and not os.path.isfile(input_fname):	#it's the input
		import StringIO
		inf = StringIO.StringIO(input_fname)
	else:
		import sys
		sys.stderr.write("Error: %s is neither a file name nor a file object. But try 'open' anyway.\n"%input_fname)
		return None
	if getattr(inf, 'readline', None) is not None and use_sniff:	#2008-01-08 don't use cs.sniff unless the user specifies it. sniff gives you unexpected delimiter when it's a single-column.
		line = inf.readline()
		delimiter_chosen = cs.sniff(line).delimiter
	else:
		line = inf.read(20000)	##binary file doesn't have readline(). have to use this dumb approach due to '\n' might mess up sniff()
		delimiter_chosen = None
		for delimiter in delimiter_choice_ls:
			delimiter_count = line.count(delimiter)
			if delimiter_count>0:
				delimiter_chosen = delimiter
				break
	del inf
	if report:
		sys.stderr.write("Done.\n")
	return delimiter_chosen

def get_gene_symbol2gene_id_set(curs, tax_id, table='genome.gene_symbol2id', upper_case_gene_symbol=0):
	"""
	2008-07-10 derived from annot.bin.codense.common.get_gene_symbol2gene_id()
	"""
	sys.stderr.write("Getting gene_symbol2gene_id_set...")
	gene_symbol2gene_id_set = {}
	from sets import Set
	curs.execute("select gene_id, gene_symbol from %s where tax_id=%s"%(table, tax_id))
	rows = curs.fetchall()
	for row in rows:
		gene_id, gene_symbol = row
		if upper_case_gene_symbol:
			gene_symbol = gene_symbol.upper()
		if gene_symbol not in gene_symbol2gene_id_set:
			gene_symbol2gene_id_set[gene_symbol] = Set()
		gene_symbol2gene_id_set[gene_symbol].add(gene_id)
	sys.stderr.write(" %s entries. Done.\n"%len(gene_symbol2gene_id_set))
	return gene_symbol2gene_id_set

def get_gene_id2gene_symbol(curs, tax_id, table='genome.gene', upper_case_gene_symbol=0):
	"""
	2008-11-14
		reverse of get_gene_symbol2gene_id_set
	"""
	sys.stderr.write("Getting gene_id2gene_symbol ...")
	gene_id2gene_symbol = {}
	from sets import Set
	curs.execute("select gene_id, gene_symbol from %s where tax_id=%s"%(table, tax_id))
	rows = curs.fetchall()
	for row in rows:
		gene_id, gene_symbol = row
		if upper_case_gene_symbol:
			gene_symbol = gene_symbol.upper()
		if gene_id not in gene_id2gene_symbol:
			gene_id2gene_symbol[gene_id] = gene_symbol
		else:
			sys.stderr.write("Warning: gene_id %s(%s) already in gene_id2gene_symbol with symbol=%s.\n"%(gene_id, gene_symbol, gene_id2gene_symbol[gene_id]))
	sys.stderr.write(" %s entries. Done.\n"%len(gene_id2gene_symbol))
	return gene_id2gene_symbol

class FigureOutTaxID(object):
	__doc__ = "2008-07-29 class to figure out tax_id using postgres database taxonomy schema"
	option_default_dict = {('hostname', 1, ): ['localhost', 'z', 1, 'hostname of the db server', ],\
							('dbname', 1, ): ['graphdb', 'd', 1, 'database name', ],\
							('schema', 1, ): ['taxonomy', 'k', 1, 'database schema name', ],\
							('db_user', 0, ): [None, 'u', 1, 'database username', ],\
							('db_passwd', 0, ): [None, 'p', 1, 'database password', ],\
							}
	def __init__(self,  **keywords):
		"""
		2008-07-29
		"""
		from ProcessOptions import ProcessOptions
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)		
	
	def curs(self):
		from db import db_connect
		conn, curs =  db_connect(self.hostname, self.dbname, self.schema, user=self.db_user, password=self.db_passwd)
		return curs
	
	curs = property(curs)
	
	def scientific_name2tax_id(self):
		scientific_name2tax_id = {}
		curs = self.curs
		curs.execute("SELECT n.name_txt, n.tax_id FROM taxonomy.names n, taxonomy.nodes o where n.name_class='scientific name' and n.tax_id=o.tax_id and o.rank='species'")
		rows = curs.fetchall()
		for row in rows:
			scientific_name, tax_id = row
			scientific_name2tax_id[scientific_name] = tax_id
		return scientific_name2tax_id
	
	scientific_name2tax_id = property(scientific_name2tax_id)
	
	def returnTaXIDGivenScientificName(self, scientific_name):
		return self.scientific_name2tax_id.get(scientific_name)
	
	def returnTaxIDGivenSentence(self, sentence):
		"""
		2008-07-29
		"""
		tax_id_to_return = None
		for scientific_name, tax_id in self.scientific_name2tax_id.iteritems():
			if sentence.find(scientific_name)>=0:
				tax_id_to_return = tax_id
				break
		return tax_id_to_return

def getColName2IndexFromHeader(header, skipEmptyColumn=False):
	"""
	2011-2-11
		add argument skipEmptyColumn
	2008-09-16
		convenient function to read input files with flexible column order.
		One variable doesn't have to be in the same column in different files, as far as the name is same.
	"""
	col_name2index = {}
	for i in range(len(header)):
		column_name = header[i]
		if skipEmptyColumn and not column_name:	#skips empty column
			continue
		col_name2index[column_name] = i
	return col_name2index

def getListOutOfStr(list_in_str, data_type=int, separator1=',', separator2='-'):
	"""
	2009-10-28
		fix a bug,
			run "start_stop_tup = map(int, start_stop_tup)" after making sure start_stop_tup is of length >1.
	2008-10-27
		fix a bug when data_type!=int. it's never right before.
	2008-09-25
		parse a list of a string representation of a list, such as '1,3-7,11'=[1,3,4,5,6,7,11]
		dash-separated representation has to be in integer.
		If all are separated by separator1, it could be in non-int data_type.
		
	"""
	list_to_return = []
	if list_in_str=='' or list_in_str is None:
		return list_to_return
	if type(list_in_str)==int:	#just one integer, put it in and return immediately
		return [list_in_str]
	index_anchor_ls = list_in_str.split(separator1)
	for index_anchor in index_anchor_ls:
		if len(index_anchor)==0:	#nothing there, skip
			continue
		start_stop_tup = index_anchor.split(separator2)
		if len(start_stop_tup)==1:
			list_to_return.append(data_type(start_stop_tup[0]))
		elif len(start_stop_tup)>1:
			start_stop_tup = map(int, start_stop_tup)
			list_to_return += range(start_stop_tup[0], start_stop_tup[1]+1)
	list_to_return = map(data_type, list_to_return)	#2008-10-27
	return list_to_return

def runLocalCommand(commandline, report_stderr=True, report_stdout=False):
		"""
		2008-1-5
			copied from utility/grid_job_mgr/hpc_cmb_pbs.py
		2008-11-07
			command_handler.communicate() is more stable than command_handler.stderr.read()
		2008-11-04
			refactor out of runRemoteCommand()
			
			run a command local (not on the cluster)
		"""
		import subprocess
		import StringIO
		command_handler = subprocess.Popen(commandline, shell=True, \
										stderr=subprocess.PIPE, stdout=subprocess.PIPE)
		#command_handler.wait() #Warning: This will deadlock if the child process generates enough output to a stdout or stderr pipe such that it blocks waiting for the OS pipe buffer to accept more data. Use communicate() to avoid that.
		
		#command_handler.stderr.read() command_handler.stdout.read() also constantly deadlock the whole process.
		
		stderr_content = None
		stdout_content = None
		
		stdout_content, stderr_content = command_handler.communicate()	#to circumvent deadlock caused by command_handler.stderr.read()
		
		output_stdout = None
		output_stderr = None
		if not report_stdout:	#if not reporting, assume the user wanna to have a file handler returned
			output_stdout = StringIO.StringIO(stdout_content)
		if not report_stderr:
			output_stderr = StringIO.StringIO(stderr_content)
		
		if report_stdout:
			sys.stderr.write('stdout of %s: %s \n'%(commandline, stdout_content))
		
		if report_stderr:
			sys.stderr.write('stderr of %s: %s \n'%(commandline, stderr_content))
		
		return_data = PassingData(commandline=commandline, output_stdout=output_stdout, output_stderr=output_stderr,\
								stderr_content=stderr_content, stdout_content=stdout_content)
		return return_data


#2009-2-4 code to link something like AT1G01010.1 or AT1G01040 to NCBI gene id
#refactored out of PutGeneListIntoDB.py
p_acc_ver = re.compile(r'(\w+)\.(\d+)')
p_acc = re.compile(r'(\w+)')	#2008-12-11 only alphanumeric characters in gene_symbol (suzi's file contains weird characters sometimes)	
def getGeneIDSetGivenAccVer(acc_ver, gene_symbol2gene_id_set, noUpperCase=0):
	if not noUpperCase:
		gene_symbol = acc_ver.upper()
	if p_acc_ver.search(gene_symbol):
		gene_symbol, version = p_acc_ver.search(gene_symbol).groups()
	if p_acc.search(gene_symbol):	#2008-12-11 pick out alphanumeric characters
		gene_symbol, = p_acc.search(gene_symbol).groups()
	gene_id_set = gene_symbol2gene_id_set.get(gene_symbol)
	return gene_id_set


def calGreatCircleDistance(lat1, lon1, lat2, lon2, earth_radius=6372.795):
	"""
	2009-4-18
		copied from CreatePopulation.cal_great_circle_distance()
	2007-06-17 copied from 2007-07-11
	http://en.wikipedia.org/wiki/Great-circle_distance
	"""
	import math
	lat1_rad = lat1*math.pi/180
	lon1_rad = lon1*math.pi/180
	lat2_rad = lat2*math.pi/180
	lon2_rad = lon2*math.pi/180
	long_diff = abs(lon1_rad-lon2_rad)
	sin_lat1 = math.sin(lat1_rad)
	cos_lat1 = math.cos(lat1_rad)
	sin_lat2 = math.sin(lat2_rad)
	cos_lat2 = math.cos(lat2_rad)
	spheric_angular_diff = math.atan2(math.sqrt(math.pow(cos_lat2*math.sin(long_diff),2) + math.pow(cos_lat1*sin_lat2-sin_lat1*cos_lat2*math.cos(long_diff),2)), sin_lat1*sin_lat2+cos_lat1*cos_lat2*math.cos(long_diff))
	return earth_radius*spheric_angular_diff
	
def addExtraToFilenamePrefix(filename, extra):
	"""
	2010-4-5
		add the extra bits (string, integer, etc.) before the file name suffix
	"""
	import os
	fname_prefix, fname_suffix = os.path.splitext(filename)
	fname_prefix += '_%s'%extra
	return fname_prefix + fname_suffix

def addExtraLsToFilenamePrefix(filename, extra_ls):
	"""
	2010-5-9
		add a list of the extra bits (string, integer, etc.) before the file name suffix
	"""
	import os
	for extra in extra_ls:
		filename = addExtraToFilenamePrefix(filename, extra)
	return filename

def returnAnyValueIfNothing(string, data_type=int, defaultValue=0):
	"""
	2010-12-15
		used in Transfac.src.GeneASNXML2gene_mapping.return_datetime() in case nothing is returned.
	"""
	if string:
		return data_type(string)
	else:
		return defaultValue


def Denary2Binary(n):
	'''
	2011-2-9
		convert denary integer n to binary string bStr
		
		copied from http://www.daniweb.com/code/snippet216539.html
		
		# convert a decimal (denary, base 10) integer to a binary string (base 2)
		# tested with Python24   vegaseat	6/1/2005
	'''
	bStr = ''
	if n < 0:  raise ValueError, "must be a positive integer"
	if n == 0: return '0'
	while n > 0:
		bStr = str(n % 2) + bStr
		n = n >> 1
	return bStr

def int2bin(n, count=24):
	"""
	2011-2-9
		opposite of Denary2Binary(), same as int(binaryStr, 2)
		
		copied from http://www.daniweb.com/code/snippet216539.html
	
	returns the binary of integer n, using count number of digits
	
	"""
	return "".join([str((n >> y) & 1) for y in range(count-1, -1, -1)])

"""
# this test runs when used as a standalone program, but not as an imported module
# let's say you save this module as den2bin.py and use it in another program
# when you import den2bin the __name__ namespace would now be  den2bin  and the
# test would be ignored
if __name__ == '__main__':
    print Denary2Binary(255)  # 11111111
    
    # convert back to test it
    print int(Denary2Binary(255), 2)  # 255
    
    print
    
    # this version formats the binary
    print int2bin(255, 12)  # 000011111111
    # test it
    print int("000011111111", 2)  # 255
    
    print

    # check the exceptions
    print Denary2Binary(0)
    print Denary2Binary(-5)  # should give a ValueError
"""

def processRegexpString(p_str):
	"""
	2011-4-30
		copied from a pylons controller. used to process regular expressions passed from web client.
	2009-4-3
		if p_str includes '(' or ')', mysql complains: ERROR 1139 (42000): Got error 'parentheses not balanced' from regexp
	"""
	p_str = p_str.replace('(', '\(')	#python re module only needs one '\'
	p_str = p_str.replace(')', '\)')
	p_str_sql = p_str.replace('(', '\\\(')	#mysql needs more
	p_str_sql = p_str_sql.replace(')', '\\\)')
	return PassingData(p_str=p_str, p_str_sql=p_str_sql)

if __name__ == '__main__':
	FigureOutTaxID_ins = FigureOutTaxID()
	print FigureOutTaxID_ins.returnTaxIDGivenSentence('>gi|172045488|ref|NW_001867254.1| Physcomitrella patens subsp. patens PHYPAscaffold_10696, whole genome shotgun sequence')
