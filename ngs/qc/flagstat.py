#!/usr/bin/env python
import os,sys

def runflagstatofsamtools(inputFile=None, outputFile=None):
    if outputFile is None:
        outputFile = os.path.basename(inputFile).strip(".bam") + ".stat"
    cmd = "samtools flagstat " + inputFile + " > " + outputFile
    os.system(cmd)

if len(sys.argv) == 2:
    runflagstatofsamtools(inputFile=sys.argv[1])
elif len(sys.argv) ==3:
    runflagstatofsamtools(inputFile=sys.argv[1], outputFile=sys.argv[2])
else :
    sys.stderr.write("argument number is wrong in flagstat...")

