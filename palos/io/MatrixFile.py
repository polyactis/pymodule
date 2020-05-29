#!/usr/bin/env python3
"""
A matrix file API like csv.reader, it can deal with any delimiter: coma, tab, space or multi-space.

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
import csv
import os
import re
import sys
from palos.ProcessOptions import  ProcessOptions
from palos import utils, figureOutDelimiter

class MatrixFile(object):
    __doc__ = __doc__
    option_default_dict = {
        ('path', 0, ): [None, 'i', 1, 'Path to the input file.'],\
        ('file_handle', 0, ): [None, '', 1, 'A opened input file handler'],\
        ('mode', 1, ): ['r', '', 1,
            'mode to open the path if file_handle is not presented.'],\
        ('delimiter', 0, ): [None, '', 1, ''],\
        ('header', 0, ): [None, '', 1, 'the header to be in output file, mode=w'],\
        ('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
        ('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']
        }
    def __init__(self, path=None, **keywords):
        self.ad = ProcessOptions.process_function_arguments(keywords,
            self.option_default_dict, error_doc=self.__doc__,
            class_to_have_attr=self)
        if not self.path:
            self.path = path
        
        if self.path and self.file_handle is None:
            self.file_handle = utils.openGzipFile(self.path, mode=self.mode)
        
        #2013.05.03 for easy access
        self.filename = self.path		
        self.csvFile = None
        self.isRealCSV = False
        if self.mode=='r':	#reading mode
            if self.delimiter is None:
                self.delimiter = figureOutDelimiter(self.file_handle)
            
            if self.delimiter=='\t' or self.delimiter==',':
                self.csvFile = csv.reader(self.file_handle, delimiter=self.delimiter)
                self.isRealCSV = True
            else:
                self.csvFile = self.file_handle
                self.isRealCSV = False
        else:	#writing mode
            if not self.delimiter:
                self.delimiter = '\t'
            self.csvFile = csv.writer(self.file_handle, delimiter=self.delimiter)
            self.isRealCSV = True
            #else:
            #	self.csvFile = self.file_handle
            #	self.isRealCSV = False
        self.col_name2index = None
        
        self._row = None	# store the current row being read
        self.headerPattern = re.compile(r'^[a-zA-Z]')
        #default header pattern, line beginned with letter
        self.commentPattern = re.compile(r'^#')	#default, beginned with #
        self.comment_row_list  = []
    
    def _resetInput(self):
        """
        """
        #reset the inf
        self.file_handle.seek(0)
    
    def writeHeader(self,header=None):
        """
        """
        if header is None:
            header = self.header
        if header:
            if self.isRealCSV:
                self.csvFile.writerow(header)
            else:
                self.csvFile.write("%s\n"%(self.delimiter.join(header)))
    
    def constructColName2IndexFromHeader(self):
        """
        """
        self.header = next(self)
        self.col_name2index = utils.getColName2IndexFromHeader(self.header)
        return self.col_name2index
    
    def smartReadHeader(self, headerPattern=None, commentPattern=None):
        """
        Note:
            If an input file does not have a header,
                this function over-reads by one line (stored in self._row)
            so need to process the last self._row before further reading
        Read the header, while ignoring lines fitting the comment pattern.
        Construct col_name2index when a line matching headerPattern is encountered
        
        """
        if headerPattern is None:
            headerPattern = self.headerPattern
        if commentPattern is None:
            commentPattern = self.commentPattern
        row = next(self)
        while commentPattern.search(row[0]):
            #passing all comments
            self.comment_row_list.append(row)
            row = next(self)
        if headerPattern.search(row[0]):
            self.header = row
            self.col_name2index = utils.getColName2IndexFromHeader(self.header)
        else:
            self.col_name2index = None
        return self.col_name2index
    
    
    def getHeader(self):
        """
        """
        return self.header
    
    def getColIndexGivenColHeader(self, colHeader=None):
        """
        """
        if self.col_name2index is None:
            #no header
            return None
        else:
            return self.col_name2index.get(colHeader)
    
    def __iter__(self):
        return self
    
    
    def __next__(self):
        try:
            row = next(self.csvFile)
        except:
            raise StopIteration
        if not self.isRealCSV:
            row = row.strip().split()
        self._row = row
        return row
    
    def close(self):
        del self.csvFile
        self.file_handle.close()
    
    def writerow(self, row=None):
        """
        mimic csv's writerow()
        """
        if row:
            if self.isRealCSV:
                self.csvFile.writerow(row)
            else:
                self.csvFile.write("%s\n"%(self.delimiter.join(row)))
    
    writeRow = writerow
    #same as writerow(), for API consistency.

    def constructDictionary(self, keyColumnIndexList=None,
        valueColumnIndexList=None, keyUniqueInInputFile=False,
        keyDataType=None, valueDataType=None):
        """
        Argument keyUniqueInInputFile
            True when each key appears once and only once in input file.
            Exception would be raise if this is not true.
        If length of keyColumnIndexList is one, 
            key is not a tuple, simply the key.
        If length of valueColumnIndexList is one, 
            value is not a list of tuples, simply a list of values.
        The key is a tuple.
        The value is a list of lists.
                
        i.e.:
        
            alignmentCoverageFile = MatrixFile(path=self.individualAlignmentCoverageFname)
            alignmentCoverageFile.constructColName2IndexFromHeader()
            alignmentReadGroup2coverageLs = alignmentCoverageFile.constructDictionary(
                keyColumnIndexList=[0], valueColumnIndexList=[1])
            alignmentCoverageFile.close()
            
            coverageLs = alignmentReadGroup2coverageLs.get((individualID,))
            return coverageLs[0]
                
            
                
        """
        sys.stderr.write("Constructing a dictionary, keys are column %s, values are column %s. ..."%\
            (repr(keyColumnIndexList), repr(valueColumnIndexList)))
        dc = {}
        counter = 0
        for row in self:
            counter += 1
            keyList = []
            for i in keyColumnIndexList:
                keyData = row[i]
                if keyDataType is not None:
                    keyData = keyDataType(keyData)
                keyList.append(keyData)
            valueList = []
            for i in valueColumnIndexList:
                valueData = row[i]
                if valueDataType is not None:
                    valueData = valueDataType(valueData)
                valueList.append(valueData)
            if len(keyColumnIndexList)>1:
                key = tuple(keyList)
            else:
                key = keyList[0]
            if keyUniqueInInputFile:
                if key in dc:
                    logging.error(f"keyUniqueInInputFile={keyUniqueInInputFile} "
                        f"but this key ({key}) from this row ({row}) is "
                        f"already in dictionary with value={dc.get(key)}.")
                    raise
            else:
                if key not in dc:
                    dc[key] = []
            
            if len(valueColumnIndexList)>1:
                value = valueList
            else:
                value = valueList[0]
            if keyUniqueInInputFile:
                dc[key] = value
            else:
                dc[key].append(value)
        sys.stderr.write("%s unique pairs from %s rows.\n"%(len(dc), counter))
        return dc
    
    
    def run(self):
        """
        """
        
        if self.debug:
            import pdb
            pdb.set_trace()
        

if __name__ == '__main__':
    main_class = MatrixFile
    po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
    instance = main_class(**po.long_option2value)
    instance.run()