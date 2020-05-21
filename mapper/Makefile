SUBDIRS = extractor
CXXStaticTargets = FindMaxLDBetweenPeakAndEachLocus CalculateColCorBetweenTwoHDF5
CXXLDFLAGS = $(HDF5Lib)
CTARGETS = 

ExtraTargets = CalculateMedianMeanOfInputColumn

#2015.03.13 use dynamic armadillo library in ubuntu. No static armadillo
CalculateMedianMeanOfInputColumn:	%:	%.o
	-$(CXXCOMPILER) $< $(CXXFLAGS) -o $@ $(StaticStartFlags) \
	-lboost_program_options -lboost_iostreams $(StaticEndFlags) -larmadillo

SRCS = CalculateColCorBetweenTwoHDF5.cc FindMaxLDBetweenPeakAndEachLocus.cc \
	CalculateMedianMeanOfInputColumn.cc

include ../Makefile.common
#all:	$(CXXTARGETS) $(CTARGETS) CalculateColCorBetweenTwoHDF5
#overwrite Makefile.common
IncludeDirs = -I`readlink -f ..` -I /usr/include/hdf5/serial/

