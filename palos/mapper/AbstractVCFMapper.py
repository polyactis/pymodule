#!/usr/bin/env python3
"""
Description:
    2012.1.17 an abstract class for VCF mapper
"""
import sys, os
import copy
from palos import ProcessOptions
from . AbstractMapper import AbstractMapper

class AbstractVCFMapper(AbstractMapper):
    __doc__ = __doc__
    option_default_dict = copy.deepcopy(AbstractMapper.option_default_dict)
    option_default_dict.update({
        ('inputFname', 0, ): ['', 'i', 1, \
            'VCF input file. either plain vcf or gzipped is ok. could be unsorted.',],
        ("chromosome", 0, ): [None, 'c', 1, \
            'chromosome name for these two VCF.'],
        ("chrLength", 1, int): [1, 'l', 1, \
            'length of the reference used for the input VCF file.'],
        ('minDepth', 0, float): [0, 'm', 1, \
            'minimum depth for a call to regarded as non-missing', ],
        })

    def __init__(self,  inputFnameLs=None, **keywords):
        """
        """
        AbstractMapper.__init__(self, inputFnameLs=inputFnameLs, **keywords)

if __name__ == '__main__':
    main_class = AbstractVCFMapper
    po = ProcessOptions(sys.argv, main_class.option_default_dict, 
        error_doc=main_class.__doc__)
    instance = main_class(**po.long_option2value)
    instance.run()