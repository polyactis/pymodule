SUBDIRS 	= 
CXXTARGETS	= 
CXXLDFLAGS 	= $(BoostLib)
CTARGETS	= 

ExtraTargets = ReduceSameChromosomeAlignmentDepthFiles
#20200521 relative path (../) is not good.
# The final program can only run in the current folder, 
#   in which the relative path is valid.

ParentObj=`readlink -f ../io/AbstractMatrixFileWalkerCC.so`
SRCS = ReduceSameChromosomeAlignmentDepthFiles.cc
# Do not put $(ParentObj) as a dependency because
#   `readlink ...` is appended as such.
#  readlink itself becomes a dependent.
ReduceSameChromosomeAlignmentDepthFiles: %: %.o
	-$(CXXCOMPILER) $< $(ParentObj) $(CXXFLAGS) \
	$(CXXLDFLAGS) -o $@

include ../Makefile.common
#overwrite Makefile.common
IncludeDirs	= -I`readlink -f ..`
