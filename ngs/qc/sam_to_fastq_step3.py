#!/usr/bin/env python
import os
import sys, shutil
import subprocess
picard_path = "/y/home/cl/software/picard.jar"

the_id = sys.argv[1]
bamdir = "bamdir_%s" % the_id

fastqfile_p1_list = []
fastqfile_p2_list = []
f_name_1 = '/y/Sunset/qualityconversion/fastqfilename/p_%s_1' %the_id
f_name_2 = '/y/Sunset/qualityconversion/fastqfilename/p_%s_2' %the_id

with open(f_name_1,'r') as f_name1, open(f_name_2,'r') as f_name2:
    for f1 in f_name1:
        fastqfile_p1_list.append(f1)
    for f2 in f_name2:
        fastqfile_p2_list.append(f2)


a = fastqfile_p1_list[0]
outdir = a.strip().split('/')[5]
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

def processfile(fastqfile_list):
    fastqfile_name = []
    for f in fastqfile_list: 
        fastqfile_name.append(f.strip().split('/')[6])
    return fastqfile_name

fastqfile_p1_name = processfile(fastqfile_p1_list)
fastqfile_p2_name = processfile(fastqfile_p2_list)
for i, p1, p2 in zip(range(len(fastqfile_p1_list)), fastqfile_p1_name, fastqfile_p2_name):
    p1 = p1.strip().split('.gz')[0]
    p2 = p2.strip().split('.gz')[0]
    bamfile = os.path.join(bamdir, "b%s.bam" % i)
    output_p1_name = os.path.join(outdir, p1)
    output_p2_name = os.path.join(outdir, p2)
    cmd = "java -Djava.io.tmpdir=$HOME/java_tmp -jar %s SamToFastq I=%s "\
          "FASTQ=%s SECOND_END_FASTQ=%s "\
          "VALIDATION_STRINGENCY=LENIENT" % (picard_path, bamfile,
          output_p1_name, output_p2_name)
    subprocess.check_call(cmd, stderr=subprocess.STDOUT, shell=True)
    cmd2 = 'gzip %s' %output_p1_name
    cmd3 = 'gzip %s' %output_p2_name
    subprocess.check_call(cmd2, stderr=subprocess.STDOUT, shell=True)
    subprocess.check_call(cmd3, stderr=subprocess.STDOUT, shell=True)

