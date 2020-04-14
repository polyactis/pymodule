#!/usr/bin/python
from pegaflow.DAX3 import *
import os,sys 
sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.expanduser('~/script'))
sys.path.insert(0, os.path.expanduser('~/src'))

import argparse
parser = argparse.ArgumentParser(description='Convert quality format of fastq files')
parser.add_argument('-i', '--inputfilename', required=True,help='the input file name containing list to be processed')
parser.add_argument('-o','--outputfilename',default='workflow.dax', help='the output file name, suffix is .dax')
parser.add_argument('-p','--password',  required=True, help='the password used in database connecting')
parser.add_argument('-d','--dbname', default='pmdb', help='database name')
parser.add_argument('-ho','--hostname', default='172.22.99.9', help='hostname of the db server')
parser.add_argument('-s','--schema', default='sunset', help='the password used in database connecting')
parser.add_argument('-u','--db_user', default='cl', help='database username')
args = parser.parse_args()

from pymodule.db import SunsetDB
db_main = SunsetDB.SunsetDB(drivername='postgresql', db_user=args.db_user,db_passwd=args.password, 
                        hostname=args.hostname, dbname=args.dbname, schema=args.schema)
db_main.setup(create_tables=False)
session = db_main.session

from pymodule.db.SunsetDB import IndividualSequence, IndividualAlignment, IndividualSequenceFile

f_in  = open(args.inputfilename, 'r') 

for each_num in f_in:
    the_id = int(each_num)
    
    fastqfile_p1_list=[]
    fastqfile_p2_list=[]
    recordlist = session.query(IndividualSequenceFile.path, \
    IndividualSequenceFile.split_order).filter(IndividualSequenceFile.\
    individual_sequence_id == the_id).all()
    length = len(recordlist)/2
    adict = {i+1:[] for i in range(length)}
    
    for record in recordlist:
        raw_path = record[0]
        split_order_id = record[1]
        abs_path=os.path.join('/y/Sunset/db',raw_path)
        if not os.path.exists(abs_path):
            abs_path=os.path.join('/y/backup/Sunset_db',raw_path)
            if not os.path.exists(abs_path):
                sys.stderr.write('The path %s is not exist!' %abs_path)
                exit(1)
        adict[split_order_id].append(abs_path)

    f_name_1 = 'p_%s_1' %the_id
    f_name_2 = 'p_%s_2' %the_id

    with open(f_name_1,'w') as f_out1, open(f_name_2,'w') as f_out2:
        for i in range(length):
            i +=1
            for inx in range(2):
                path = adict[i][inx]
                if path.endswith("_1_%s.fastq.gz" % i):
                    f_out1.write("%s\n" % path)
                elif path.endswith("_2_%s.fastq.gz" % i):
                    f_out2.write("%s\n" % path)
                else:
                    sys.stderr.write("Error in getting fastq file\n in path %s" %path)
                    exit(1)
        

    

