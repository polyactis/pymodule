#!/usr/bin/env python3
"""
Five types of figures will be generated.
    qualityHist.png
    quality_per_position.png
    no_of_bases_per_position.png
    diNuc_count.png
    diNuc_quality.png

The db connection is un-used, you can supply any random password.

Examples:
    %s -i gerald_62FGFAAXX_3_1.fastq.gz -u yh -read_sampling_rate 0.0001

"""

import sys, os, copy
__doc__ = __doc__%(sys.argv[0])

import random
#to disable pop-up
import matplotlib; matplotlib.use("Agg")
from palos import PassingData, ProcessOptions, utils
from palos.plot import yh_matplotlib
from palos.mapper.AbstractSunsetMapper import AbstractSunsetMapper

ParentClass = AbstractSunsetMapper
class InspectBaseQuality(ParentClass):
    __doc__ = __doc__
    option_default_dict = copy.deepcopy(ParentClass.option_default_dict)
    option_default_dict.pop(('outputFname', 0, ))
    #option_default_dict.pop(('outputFnamePrefix', 0, ))
    option_default_dict.update({
        ('quality_score_format', 1, ): ['Sanger', 'q', 1, 
            'Can be Sanger (phred+33), or Illumina1.3 (~phred+64). '
            'Illumina1.8+ (after 2011-02) is Sanger.', ],\
        ('read_sampling_rate', 0, float): [0.05, '', 1, 
            'The probability that a read will be included.', ],\
        ('ind_sequence_id', 0, int): [None, '', 1, 
            'The individual_sequence id corresponding to the inputFname, '
            'if not given, no db commit action.', ],\
        })

    def __init__(self,  inputFnameLs=None, **keywords):
        """
        """
        ParentClass.__init__(self, inputFnameLs, **keywords)
        
        if not self.outputFnamePrefix:
            self.outputFnamePrefix = utils.getRealPrefixSuffix(
                self.inputFname)[0]
    
    def getQualityData(self, inputFname, read_sampling_rate=0.05,
        quality_score_format='Sanger'):
        """
        """
        print(f"Getting base quality data from {inputFname} ...", flush=True)
        quality_ls_per_position = []
        quality_ls = []
        no_of_bases_per_position = []
        diNuc2count = {}
        diNuc2quality_ls = {}
        
        inf = utils.openGzipFile(inputFname, 'r')
        counter = 0
        real_counter = 0
        for line in inf:
            if line[0]=='@':
                counter += 1
                coin_toss = random.random()
                base_string = inf.readline().strip()
                inf.readline()
                quality_string = inf.readline().strip()
                if coin_toss<=read_sampling_rate:
                    real_counter += 1
                    read_length = len(base_string)
                    if len(quality_ls_per_position)<read_length:
                        # extend quality_ls_per_position to house more data
                        extraNoOfBases = read_length-len(quality_ls_per_position)
                        for j in range(extraNoOfBases):
                            quality_ls_per_position.append([])
                            no_of_bases_per_position.append(0)
                        
                    for i in range(read_length):
                        base = base_string[i]
                        base_quality = quality_string[i]
                        if quality_score_format=='Illumina1.3':
                            phredScore = utils.converSolexaScoreToPhred(base_quality)
                        else:
                            phredScore = ord(base_quality)-33
                        quality_ls_per_position[i].append(phredScore)
                        quality_ls.append(phredScore)
                        if base!='N':
                            no_of_bases_per_position[i] += 1
                            if i<read_length-1:
                                nextBase = base_string[i+1]
                                if nextBase!='N':
                                    diNuc = base + nextBase
                                    if diNuc not in diNuc2quality_ls:
                                        diNuc2quality_ls[diNuc] = []
                                        diNuc2count[diNuc] = 0
                                    diNuc2quality_ls[diNuc].append(phredScore)
                                    diNuc2count[diNuc] += 1
            if counter%5000==0 and self.report:
                sys.stderr.write("%s%s\t%s"%('\x08'*80, real_counter, counter))
            
            #if baseCount>10000:	#temporary, for testing
            #	break
        del inf
        print(f"{real_counter}/{counter} reads selected.", flush=True)
        return PassingData(
            quality_ls_per_position=quality_ls_per_position, 
            quality_ls=quality_ls, \
            no_of_bases_per_position=no_of_bases_per_position,
            diNuc2quality_ls=diNuc2quality_ls,
            diNuc2count=diNuc2count)
    
    def drawQualityData(self, qualityDataStructure, outputFnamePrefix,
        sequence_id=''):
        """
        """
        print("Making plots on quality data ...", flush=True)
        
        yh_matplotlib.drawHist(qualityDataStructure.quality_ls, 
            title='histogram of phredScore from %s'%(sequence_id),
            xlabel_1D=None,
            outputFname='%s_qualityHist.png'%(outputFnamePrefix),
            min_no_of_data_points=50, needLog=False, dpi=200)
        
        yh_matplotlib.drawBoxPlot(qualityDataStructure.quality_ls_per_position,
            title='quality box plot from %s'%(sequence_id),
            xlabel_1D='base position in read',
            xticks=None,
            outputFname='%s_quality_per_position.png'%(outputFnamePrefix),
            dpi=200)
        
        no_of_bases_per_position = qualityDataStructure.no_of_bases_per_position
        readLength = len(no_of_bases_per_position)
        yh_matplotlib.drawBarChart(range(1, readLength+1),
            no_of_bases_per_position,
            title='no of base calls from %s'%(sequence_id),
            xlabel_1D='base position in read', xticks=None,
            outputFname='%s_no_of_bases_per_position.png'%(outputFnamePrefix),
            bottom=0, needLog=False, dpi=200)
        
        diNuc2count = qualityDataStructure.diNuc2count
        diNuc2quality_ls = qualityDataStructure.diNuc2quality_ls
        
        diNuc_key_ls = sorted(diNuc2count)
        diNuc_count_ls = []
        diNuc_quality_ls_ls = []
        for diNuc in diNuc_key_ls:
            diNuc_count_ls.append(diNuc2count.get(diNuc))
            diNuc_quality_ls_ls.append(diNuc2quality_ls.get(diNuc))
        
        yh_matplotlib.drawBarChart(range(1, len(diNuc_count_ls)+1),
            diNuc_count_ls,
            title='di-nucleotide counts from %s'%(sequence_id),
            xlabel_1D=None, xticks=diNuc_key_ls,
            outputFname='%s_diNuc_count.png'%(outputFnamePrefix),
            bottom=0, needLog=False, dpi=200)
        
        yh_matplotlib.drawBoxPlot(diNuc_quality_ls_ls,
            title='di-Nucleotide quality box plot from %s'%(sequence_id),
            xlabel_1D=None, xticks=diNuc_key_ls,
            outputFname='%s_diNuc_quality.png'%(outputFnamePrefix),
            dpi=200)
        
        print("Done.", flush=True)
    
    def saveDataIntoDB(self, db_main, ind_sequence_id=None):
        """
        """
        pass
    
    def run(self):
        """
        """
        if self.debug:
            import pdb
            pdb.set_trace()
        
        session = self.db_main.session
        session.begin()
        #no transaction for input node as there is no data insertion
        
        qualityDataStructure = self.getQualityData(
            self.inputFname,
            read_sampling_rate=self.read_sampling_rate,
            quality_score_format=self.quality_score_format)
        sequence_id = os.path.split(self.outputFnamePrefix)[1]
        #to be part of title in each figure
        self.drawQualityData(qualityDataStructure, self.outputFnamePrefix,
            sequence_id=sequence_id)
        
        self.saveDataIntoDB(self.db_main, ind_sequence_id=self.ind_sequence_id)
        
        self.outputLogMessage('Inspect Done.\n')

        if self.commit:
            session.flush()
            session.commit()

if __name__ == '__main__':
    main_class = InspectBaseQuality
    po = ProcessOptions(sys.argv, main_class.option_default_dict,
        error_doc=main_class.__doc__)
    instance = main_class(po.arguments, **po.long_option2value)
    instance.run()
