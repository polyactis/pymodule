
CXXTARGETS	= 
CXXLDFLAGS	= $(BoostLib)
CTARGETS	= 

ExtraTargets	= SelectRowsFromMatrixCC
ParentObj = `readlink -f ../../io/AbstractMatrixFileWalkerCC.so`

SRCS	= SelectRowsFromMatrixCC.cc

SelectRowsFromMatrixCC:	%:	%.o
	-$(CXXCOMPILER) $< $(ParentObj) $(CXXFLAGS) $(CXXLDFLAGS) -o $@


include ../../Makefile.common
#overwrite Makefile.common
IncludeDirs	= -I`readlink -f ../..`
