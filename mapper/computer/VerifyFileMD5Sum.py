#!/usr/bin/env python
"""
Examples:
    %s -i ABC.bam --correct_md5sum 486fd350d10cd3d701fa5dee22fca18d

Description:
    2017.04.07
        Verify if the input file's md5sum is equal to the given correct md5sum.
"""

import sys, os, math
__doc__ = __doc__%(sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

from pymodule import ProcessOptions, PassingData, utils
from pymodule.mapper.AbstractMapper import AbstractMapper

class VerifyFileMD5Sum(AbstractMapper):
    __doc__ = __doc__
    option_default_dict = AbstractMapper.option_default_dict.copy()
    option_default_dict.update({
                            ('correct_md5sum', 1, ): ["", 'q', 1, 'The correct md5sum for the input file '],\
                            })
    option_default_dict.pop(('outputFname', 0, ))
    option_default_dict.pop(('outputFnamePrefix', 0, ))
    def __init__(self, inputFnameLs=None, **keywords):
        """
        """
        AbstractMapper.__init__(self, **keywords)
        self.inputFnameLs = inputFnameLs    #2012.3.19 not used
    
    
    def run(self):
        """
        """
        
        if self.debug:
            import pdb
            pdb.set_trace()
        md5sum = utils.get_md5sum(self.inputFname)
        if md5sum==self.correct_md5sum:
            sys.stderr.write("md5sum of %s matches given md5sum %s.\n"%(self.inputFname, self.correct_md5sum))
            sys.exit(0)
        else:
            sys.stderr.write("ERROR: md5sum of %s, %s, does not match given md5sum %s.\n"%(self.inputFname, 
                                                                                           md5sum, self.correct_md5sum))
            sys.exit(2)

    
if __name__ == '__main__':
    main_class = VerifyFileMD5Sum
    po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
    instance = main_class(po.arguments, **po.long_option2value)
    instance.run()