
CXXTARGETS	= 
CXXLDFLAGS	= $(BoostLib)
CTARGETS	= 

ExtraTargets	= SelectRowsFromMatrixCC
SharedLibTargets  = ../../io/AbstractMatrixFileWalkerCC.so

SRCS	= SelectRowsFromMatrixCC.cc

SelectRowsFromMatrixCC:	%:	%.o $(SharedLibTargets)
	-$(CXXCOMPILER) $< $(SharedLibTargets) $(CXXFLAGS) $(CXXLDFLAGS) -o $@


include ../../Makefile.common
#overwrite Makefile.common
IncludeDirs	= -I`readlink -f ../..`
