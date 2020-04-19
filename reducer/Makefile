SUBDIRS 	= 
CXXTARGETS	= 
CXXLDFLAGS 	= $(BoostLib)
CTARGETS	= 

ExtraTargets = ReduceSameChromosomeAlignmentDepthFiles
SharedLibTargets = ../io/AbstractMatrixFileWalkerCC.so

SRCS = ReduceSameChromosomeAlignmentDepthFiles.cc

ReduceSameChromosomeAlignmentDepthFiles:	%:	%.o $(SharedLibTargets)
	-$(CXXCOMPILER) $< $(SharedLibTargets) $(CXXFLAGS) $(CXXLDFLAGS) -o $@

include ../Makefile.common
#overwrite Makefile.common
IncludeDirs	= -I`readlink -f ..`
