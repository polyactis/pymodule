"""
Palos is a Python3 module developed and used by the yfish group, http://www.yfish.org/.
It contains code related to bioinformatics projects focusing on next-generation sequencing data, 
population genetics, genome-wide association studies, pedigree genetics, etc.
"""
version='0.1.26'
from . ProcessOptions import ProcessOptions, generate_program_doc, process_options, \
    process_function_arguments, turn_option_default_dict2argument_default_dict
from . utils import PassingData, PassingDataList, dict_map, importNumericArray, \
    figureOutDelimiter, get_gene_symbol2gene_id_set, \
    FigureOutTaxID, getColName2IndexFromHeader, getListOutOfStr, runLocalCommand, \
    getGeneIDSetGivenAccVer, calGreatCircleDistance, openGzipFile
from . Genome import GeneModel

from . algorithm import pca_module
from . algorithm.PCA import PCA
from . algorithm.RBTree import RBTree, RBDict, RBTreeIter, RBList, RBNode
from . algorithm.BinarySearchTree import binary_tree
from . io.MatrixFile import MatrixFile
