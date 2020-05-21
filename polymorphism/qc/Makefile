SUBDIRS 	= 
CXXTARGETS	=
CXXLDFLAGS 	= $(BoostLib)
CTARGETS	= 

#ExtraTargets	= CheckTwoVCFOverlapCC
ParentObj = `readlink -f ../../io/AbstractMatrixFileWalkerCC.so`
#CheckTwoVCFOverlapCC:	%:	%.o
#	-$(CXXCOMPILER) $< ~/script/vcflib/libvcf.a $(ParentObj) $(CXXFLAGS) -o $@ $(CXXLDFLAGS)
#SRCS	= CheckTwoVCFOverlapCC.cc
include ../../Makefile.common
#all:	$(CXXTARGETS) $(CTARGETS) CalculateColCorBetweenTwoHDF5
#overwrite Makefile.common
IncludeDirs	= -I`readlink -f ../..`
