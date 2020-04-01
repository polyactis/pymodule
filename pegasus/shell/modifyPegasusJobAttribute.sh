#!/bin/bash

source ~/.bash_profile
runTypeDefault=0
jobFilenameSuffixDefault=.sub
if test $# -lt 4; then
	echo "Usage:"
	echo "  $0 pegasusWorkFolder jobFilenamePrefix attributeName attributeNewValue [runType] [jobFilenameSuffix]"
	echo ""
	echo "Note:"
	echo " 1. This program modifies pegasus job scripts (jobFilenamePrefix*.sub, .sub is controlled by argument jobFilenameSuffix) in the working folder." 
	echo " 2. It is only suited for modifying numerical attributes (re pattern: [0-9]*), like changing 82800 in 'TimeToLive>=82800' to 178000. '>=' could any combination of <>= (at least one) and could have optional space in the front or end, i.e. 'TimeToLive = 82800' or 'TimeToLive=   82800'."
	echo " 3. attributeName: attribute of a condor job. i.e. TimeToLive, request_memory, (memory.".
	echo "    To modify memory request, you need to run this program twice to change request_memory and (memory."
	echo " 4. runType: Default runType is $runTypeDefault."
	echo "    0: this program leaves the new content in jobFilename.tmp but do not overwrite the original one."
	echo "    1: this program overwrites (by mv command) the original FILEPATH with FILEPATH.tmp (and FILEPATH.tmp will be gone)."
	echo "    2: this program deletes FILEPATH.tmp left by runType 0."
	echo "    3: same effect as runType 1 but no user-asking, just execute."
	echo " 5. jobFilenameSuffix: Default is $jobFilenameSuffixDefault. Restrict to modify only these files."
	echo
	echo "Examples:"
	echo "	# replace old value TimeToLive with 178000 and overwrite original file"
	echo "	$0 work/LocalRealignmentBQSR_AlnID2828_2847_vsMethod87.2013.Apr.17T230821/ merge_pegasus-AddOrReplaceReadGroupsJava TimeToLive 178000 1"
	echo
	echo "  # dry-run: replace old MaxPermSize in java jobs (original files not affected), .in files."
	echo "	$0 work/LocalRealignmentBQSR_AlnID2828_2847_vsMethod87.2013.Apr.17T230821/ merge_pegasus-AddOrReplaceReadGroupsJava MaxPermSize 4000 0 .in"
	echo
	echo "	# dry-run, change request_memory to 9040 (MB)"
	echo "  $0 work/InspectPopulationMonkeyAlignment_RefSeq3488_AlnMethod6_GATKDOC.2013.Jun.21T105136/  DOCWalkerJava_ID00 request_memory 9040 0"
	echo
	echo "  # dry-run, change memory to 9040 (MB)"
	echo ' 	$0 work/InspectPopulationMonkeyAlignment_RefSeq3488_AlnMethod6_GATKDOC.2013.Jun.21T105136/  DOCWalkerJava_ID00 "(memory" 9040 0'
	echo
	echo
	exit 1
fi

#debug mode
#set -vx
pegasusWorkFolder=$1
jobFilenamePrefix=$2
attributeName=$3
attributeNewValue=$4
runType=$5
jobFilenameSuffix=$6

if [ -z $runType ]; then
	runType=$runTypeDefault
fi

if [ -z $jobFilenameSuffix ]; then
	jobFilenameSuffix=$jobFilenameSuffixDefault
fi

echo attributeName is $attributeName
echo attributeNewValue is $attributeNewValue
echo runType is $runType


affectedFiles=`ls $pegasusWorkFolder/$jobFilenamePrefix*$jobFilenameSuffix`

echo "These files have matches: "
echo  $affectedFiles

counter=0

if test "$runType" = "1"; then
	echo "File overwriting option is on."
	echo -n "Continue? (y/n): "
	read answer
	if [ -z $answer ]; then
		exit 1
	fi
	if [ $answer != 'y' ]; then
		exit 1
	fi
fi

for f in `ls $pegasusWorkFolder/$jobFilenamePrefix*$jobFilenameSuffix`; do
	sed 's/'$attributeName'\([ ]*[<>=][<>=]*[ ]*\)[0-9]*/'$attributeName'\1'$attributeNewValue'/' $f >$f.tmp;
	counter=`echo $counter+1|bc`
	if test "$runType" = "1" || test "$runType" = "3"; then
		mv $f.tmp $f
	elif test "$runType" = "2"; then
		rm $f.tmp
	fi
done

echo "Modified $counter files."
