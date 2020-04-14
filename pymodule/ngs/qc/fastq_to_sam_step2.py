#!/usr/bin/env python
import os, sys
import subprocess
import shutil

the_id = sys.argv[1]
outdir = "bamdir_%s" % the_id
if os.path.isfile(outdir):
    sys.stderr.write("Output dir %s is a file. Remove it.\n" % outdir)
    os.remove(outdir)
    os.makedirs(outdir)
elif os.path.isdir(outdir):
    sys.stderr.write("Output dir %s exists, clean it\n" % outdir)
    shutil.rmtree(outdir)
    os.makedirs(outdir)
else:
    os.makedirs(outdir)

picard_path = "/y/home/cl/software/picard.jar"

fastqfile_p1_list = []
fastqfile_p2_list = []
f_name_1 = '/y/Sunset/qualityconversion/fastqfilename/p_%s_1' %the_id
f_name_2 = '/y/Sunset/qualityconversion/fastqfilename/p_%s_2' %the_id


with open(f_name_1,'r') as f_name1, open(f_name_2,'r') as f_name2:
    for f1 in f_name1:
        fastqfile_p1_list.append(f1)
    for f2 in f_name2:
        fastqfile_p2_list.append(f2)

#picard FastqToSam  

for i, p1, p2 in zip(range(len(fastqfile_p1_list)), 
                     fastqfile_p1_list, fastqfile_p2_list):
    outputfile_name = os.path.join(outdir, "b%s.bam" % i)
    p1 = p1.strip()
    p2 = p2.strip()
    cmd = "java -Djava.io.tmpdir=$HOME/java_tmp -jar %s FastqToSam "\
          "F1=%s F2=%s O=%s SM=converted_seq" % (picard_path,
            p1, p2, outputfile_name)
    subprocess.check_call(cmd, stderr=subprocess.STDOUT, shell=True)