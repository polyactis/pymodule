#!/usr/bin/python3
import os,sys 
import argparse
from pegaflow.api import Workflow
import pegaflow
from palos.db import SunsetDB

parser = argparse.ArgumentParser(description='Convert quality format of fastq files')
parser.add_argument('-i', '--inputfilename', required=True,help='the input file name containing list to be processed')
parser.add_argument('-o','--output_path',default='workflow.yml', help='the output file name, suffix is .dax')
parser.add_argument('-p','--password',  required=True, help='the password used in database connecting')
parser.add_argument('-d','--dbname', default='pmdb', help='database name')
parser.add_argument('-ho','--hostname', default='172.22.99.9', help='hostname of the db server')
parser.add_argument('-s','--schema', default='sunset', help='the password used in database connecting')
parser.add_argument('-u','--db_user', default='cl', help='database username')
args = parser.parse_args()

db_main = SunsetDB.SunsetDB(drivername='postgresql', db_user=args.db_user,db_passwd=args.password, 
                        hostname=args.hostname, dbname=args.dbname, schema=args.schema)
db_main.setup(create_tables=False)
session = db_main.session


wflow = Workflow("pegasus_test")
site_handler = "condor"
picard_path = "/y/home/cl/software/picard.jar"

##############################################################################################
# execute register                                                                           #
##############################################################################################
fastq2sam_exe = pegaflow.registerExecutable(wflow, "fastq_to_sam_step2.py", site_handler)
sam2fastq_exe = pegaflow.registerExecutable(wflow, "sam_to_fastq_step3.py", site_handler)

f_in  = open(args.inputfilename, 'r')
for each_num in f_in:
    the_id = int(each_num)
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
    filenamelst = []
    for f1, f2 in zip(fastqfile_p1_list, fastqfile_p2_list): 
        filenamelst.append(os.path.join(outdir,f1.strip().split('/')[6]))
        filenamelst.append(os.path.join(outdir,f2.strip().split('/')[6]))


    the_id = str(the_id)
    fastq2sam_job = pegaflow.addJob2workflow(wflow, fastq2sam_exe, argv=[the_id],
        job_max_memory=10000, no_of_cpus=1, walltime=40000)
    

    sam2fastq_job = pegaflow.addJob2workflow(wflow, sam2fastq_exe, [the_id], 
        output_file_transfer_list=filenamelst, job_max_memory=10000, no_of_cpus=1,
        walltime=40000)
    wflow.add_dependency(sam2fastq_job, parents=[fastq2sam_job])

wflow.write(args.output_path)
