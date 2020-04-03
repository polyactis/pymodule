"""
2007-10-15. __init__.py for pymodule. pymodule is a concatenation of all common functions/classes.
"""
import os, sys
sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.expanduser('~/script'))

from . ProcessOptions import ProcessOptions, generate_program_doc, process_options, process_function_arguments, turn_option_default_dict2argument_default_dict
from . utils import PassingData, PassingDataList, dict_map, importNumericArray, figureOutDelimiter, get_gene_symbol2gene_id_set, \
	FigureOutTaxID, getColName2IndexFromHeader, getListOutOfStr, runLocalCommand, getGeneIDSetGivenAccVer, calGreatCircleDistance,\
	openGzipFile
from . Genome import GeneModel

#from . plot import yh_matplotlib, yh_matplotlib_artists
#from . plot.yh_matplotlib import assignMatPlotlibHueColorToLs, drawName2FCLegend

#from db import GenomeDB, TaxonomyDB

from . algorithm import pca_module
from . algorithm.PCA import PCA
from . algorithm.RBTree import RBTree, RBDict, RBTreeIter, RBList, RBNode
from . algorithm.BinarySearchTree import binary_tree

#2012.11.14
#from . yhio import SNP
#from . yhio import NextGenSeq, getAssociationLandscapeDataFromHDF5File
#from . yhio.MatrixFile import MatrixFile
#from . yhio.BeagleLikelihoodFile import BeagleLikelihoodFile
#from . yhio.PlinkPedigreeFile import PlinkPedigreeFile
#
#from . yhio.SNP import write_data_matrix, read_data, SNPData, GenomeWideResults, GenomeWideResult, DataObject, \
#	getGenomeWideResultFromFile,\
#	nt2number, number2nt, number2complement, SNPInfo, number2single_char_nt, getGenomeWideResultFromHDF5MatrixFile,\
#	SNPData2RawSnpsData_ls
#from . yhio.TwoSNPData import TwoSNPData, QualityControl
#from . yhio.HDF5MatrixFile import HDF5MatrixFile, addAttributeDictToYHTableInHDF5Group
#from . yhio.YHPyTables import YHTable, YHFile, YHSingleTableFile, castPyTablesRowIntoPassingData
#from . yhio.Association import AssociationTable, AssociationTableFile, LocusMapTable, LocusMapTableFile
#from . yhio.AssociationLandscape import AssociationLandscapeTable, AssociationLandscapeTableFile
#from . yhio.AssociationPeak import AssociationPeakTable, AssociationPeakTableFile
#from . yhio.AssociationLocus import AssociationLocusTable, AssociationLocus2PeakTable, AssociationLocusTableFile
#
#from yhio.BamFile import BamFile
#from . yhio.VCFFile import VCFFile, VCFRecord
#from . yhio import CNV
#from . yhio.CNV import CNVCompare, CNVSegmentBinarySearchTreeKey, get_overlap_ratio, findCorrespondenceBetweenTwoCNVRBDict
#from . yhio.AbstractMatrixFileWalker import AbstractMatrixFileWalker
#from . yhio.latex import outputMatrixInLatexTable, outputFigureInLatex
#from . yhio.FastaFile import FastaFile
