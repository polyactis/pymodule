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

workflow_AC = ADAG("pegasus_test")
site_handle = "ycondor"
picard_path = "/y/home/cl/software/picard.jar"

def registerExecutefile(workflow, executeFile):
    architecture = "x86_64"
    operatingSystem = "linux"
    namespace = "pegasus"
    version = "1.0"
    executeName = os.path.basename(executeFile)
    execute = Executable(namespace=namespace, name=executeName, os=operatingSystem, arch=architecture, installed=True)
    execute.addPFN(PFN("file://" + os.path.abspath(executeFile), site_handle))
    workflow.addExecutable(execute)
    return executeName

def addfastqtosamToWorkflow(workflow, excute, argv):
    step = Job(namespace="pegasus", name=excute)
    step.addArguments(*argv)
    workflow.addJob(step)
    return step

def addsamtofastqToWorkflow(workflow, excute,filenamelst, argv):
    step = Job(namespace="pegasus", name=excute)
    step.addArguments(*argv)
    for f in filenamelst:
        step.uses(f, link=Link.OUTPUT)
    workflow.addJob(step)
    return step


def setJobToProperMemoryRequirement(job=None, job_max_memory=500, no_of_cpus=1, walltime=180, sshDBTunnel=0):
    """
    2013.06.22 if job_max_memory is None, then skip setting memory requirement
        if job_max_memory is "" or 0 or "0", then assign 500 (mb) to it.
    2012.8.15
        increase default walltime to 180
    2012.4.16
        add argument sshDBTunnel.
            =1: this job needs a ssh tunnel to access psql db on dl324b-1.
            =anything else: no need for that.
    2011-11-23
        set walltime default to 120 minutes (2 hours)
    2011-11-16
        add more requirements
    2011-11-11
        job_max_memory is in MB.
        walltime is in minutes.
    """
    condorJobRequirementLs = []
    if job_max_memory == "" or job_max_memory == 0 or job_max_memory == "0":
        job_max_memory = 500 
    if job_max_memory is not None:
        job.addProfile(Profile(Namespace.GLOBUS, key="maxmemory", value="%s" % (job_max_memory)))
        job.addProfile(
            Profile(Namespace.CONDOR, key="request_memory", value="%s" % (job_max_memory)))  # for dynamic slots
        condorJobRequirementLs.append("(memory>=%s)" % (job_max_memory))
    # 2012.4.16
    if sshDBTunnel == 1:
        condorJobRequirementLs.append("(sshDBTunnel==%s)" % (sshDBTunnel))  # use ==, not =.

    if no_of_cpus is not None:
        job.addProfile(Profile(Namespace.CONDOR, key="request_cpus", value="%s" % (no_of_cpus)))  # for dynamic slots   

    if walltime is not None:
        # 2013.3.21 scale walltime according to clusters_size
        job.addProfile(Profile(Namespace.GLOBUS, key="maxwalltime", value="%s" % (walltime)))
        # TimeToLive is in seconds
        condorJobRequirementLs.append("(Target.TimeToLive>=%s)" % (int(walltime) * 60))
    # key='requirements' could only be added once for the condor profile
    job.addProfile(Profile(Namespace.CONDOR, key="requirements", value=" && ".join(condorJobRequirementLs)))

##############################################################################################
# execute register                                                                           #
##############################################################################################
fastq2sam_exe = registerExecutefile(workflow_AC, "fastq_to_sam_step2.py")
sam2fastq_exe = registerExecutefile(workflow_AC, "sam_to_fastq_step3.py")


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
    fastq2sam_job = addfastqtosamToWorkflow(workflow_AC, fastq2sam_exe, [the_id])
    setJobToProperMemoryRequirement(job=fastq2sam_job, job_max_memory=10000, no_of_cpus=1, walltime=40000, sshDBTunnel=0)
    

    sam2fastq_job = addsamtofastqToWorkflow(workflow_AC, sam2fastq_exe,filenamelst, [the_id])
    setJobToProperMemoryRequirement(job=sam2fastq_job, job_max_memory=10000, no_of_cpus=1, walltime=40000, sshDBTunnel=0)
    workflow_AC.addDependency(Dependency(parent=fastq2sam_job, child=sam2fastq_job))

with open(args.outputfilename, 'w') as f_out:
    workflow_AC.writeXML(f_out)
