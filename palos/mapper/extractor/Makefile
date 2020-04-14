
CXXTARGETS	= 
CXXLDFLAGS	= $(BoostLib)
CTARGETS	= 

ExtraTargets	= SelectRowsFromMatrixCC
SharedLibTargets	= ~/script/pymodule/io/AbstractMatrixFileWalkerCC.so

SRCS	= SelectRowsFromMatrixCC.cc

SelectRowsFromMatrixCC:	%:	%.o $(SharedLibTargets)
	-$(CXXCOMPILER) $< $(SharedLibTargets) $(CXXFLAGS) $(CXXLDFLAGS) -o $@


include ../../Makefile.common
