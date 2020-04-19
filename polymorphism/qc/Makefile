SUBDIRS 	= 
CXXTARGETS	=
CXXLDFLAGS 	= $(BoostLib)
CTARGETS	= 

#ExtraTargets	= CheckTwoVCFOverlapCC
SharedLibTargets	= ../../io/AbstractMatrixFileWalkerCC.so 
#CheckTwoVCFOverlapCC:	%:	%.o $(SharedLibTargets)
#	-$(CXXCOMPILER) $< ~/script/vcflib/libvcf.a $(SharedLibTargets) $(CXXFLAGS) -o $@ $(CXXLDFLAGS)
#SRCS	= CheckTwoVCFOverlapCC.cc
include ../../Makefile.common
#all:	$(CXXTARGETS) $(CTARGETS) CalculateColCorBetweenTwoHDF5
#overwrite Makefile.common
IncludeDirs	= -I`readlink -f ../..`
