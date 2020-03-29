#!/usr/bin/env python3
"""
An example Pegasus workflow that does not use class.
"""
import copy
import sys, os
from argparse import ArgumentParser
from pegapy3.DAX3 import Executable, File, PFN, Profile, Namespace, Link, ADAG, Use, Job, Dependency

sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))
sys.path.insert(0, os.path.join(os.path.expanduser('~/src')))
site_handler = "ycondor"
version = "1.0"
namespace = "pegasus"

def registerExecutefile(workflow, executeFile):
    architecture = "x86_64"
    operatingSystem = "linux"
    executeName = os.path.basename(executeFile)
    execute = Executable(namespace=namespace, name=executeName,
                         os=operatingSystem, arch=architecture,
                         installed=True, version=version)
    execute.addPFN(PFN("file://" + os.path.abspath(executeFile), site_handler))
    workflow.addExecutable(execute)
    return executeName


def setJobToProperMemoryRequirement(job=None, job_max_memory=500, no_of_cpus=1,
                                    walltime=180, sshDBTunnel=0):
    """
        set walltime default to 120 minutes (2 hours)
        job_max_memory is in MB.
        walltime is in minutes.
    """
    condorJobRequirementLs = []
    if job_max_memory == "" or job_max_memory == 0 or job_max_memory == "0":
        job_max_memory = 500
    if job_max_memory is not None:
        job.addProfile(Profile(Namespace.GLOBUS, key="maxmemory",
                               value="%s" % (job_max_memory)))
        job.addProfile(
            Profile(Namespace.CONDOR, key="request_memory", value="%s"
                    % (job_max_memory)))  # for dynamic slots
        condorJobRequirementLs.append("(memory>=%s)" % (job_max_memory))
    if sshDBTunnel == 1:
        condorJobRequirementLs.append("(sshDBTunnel==%s)" % (sshDBTunnel))
    if no_of_cpus is not None:
        job.addProfile(Profile(Namespace.CONDOR, key="request_cpus", value="%s"
                               % (no_of_cpus)))
    if walltime is not None:
        job.addProfile(Profile(Namespace.GLOBUS, key="maxwalltime",
                               value="%s" % (walltime)))
        condorJobRequirementLs.append("(Target.TimeToLive>=%s)" %
                                      (int(walltime) * 60))
    job.addProfile(Profile(Namespace.CONDOR, key="requirements",
                           value=" && ".join(condorJobRequirementLs)))

def registerOneInputFile(workflow=None, inputFname=None, input_site_handler=None, folderName="", \
                            useAbsolutePathAsPegasusFileName=False,\
                            pegasusFileName=None, checkFileExistence=True):
    """
    Examples:
        pegasusFile = registerOneInputFile(workflow=workflow, inputFname=path, input_site_handler=site_handler, \
                                        folderName=folderName, useAbsolutePathAsPegasusFileName=useAbsolutePathAsPegasusFileName)
    2013.06.29 added argument checkFileExistence
    2013.04.07 raise if inputFname is not a file
    2013.2.14 added argument useAbsolutePathAsPegasusFileName
        This would render the file to be referred as the absolute path on the running computer.
        And pegasus will not seek to symlink or copy/transfer the file.
        set it to True only when you dont want to add the file to the job as INPUT dependency (as it's accessed through abs path).
    2013.1.10 make sure the file is not registed with the workflow already
    2012.3.22
        add abspath attribute to file.
    2012.3.1
        add argument folderName, which will put the file in specific pegasus workflow folder
    2011.12.21
    """
    if not pegasusFileName:
        if useAbsolutePathAsPegasusFileName:
            pegasusFileName = os.path.abspath(inputFname)	#this will stop symlinking/transferring , and also no need to indicate them as file dependency for jobs.
        else:
            pegasusFileName = os.path.join(folderName, os.path.basename(inputFname))
    pegasusFile = File(pegasusFileName)
    if checkFileExistence and not os.path.isfile(inputFname):	#2013.06.29
        sys.stderr.write("Error from registerOneInputFile(): %s does not exist.\n"%(inputFname))
        raise
    pegasusFile.abspath = os.path.abspath(inputFname)
    pegasusFile.absPath = pegasusFile.abspath
    pegasusFile.addPFN(PFN("file://" + pegasusFile.abspath, input_site_handler))
    if not workflow.hasFile(pegasusFile):	#2013.1.10
        workflow.addFile(pegasusFile)
    return pegasusFile


def registerFile(workflow, filepath):
    file = File(os.path.basename(filepath))
    file.addPFN(PFN("file://" + os.path.abspath(filepath), site_handler))
    workflow.addFile(file)
    return file


def addJob2workflow(workflow, excute, file_input_list, file_output_transfer,
                    file_output_notransfer, argv):
    the_job = Job(namespace="pegasus", name=excute, version=version)
    if argv:
        the_job.addArguments(*argv)
    if file_input_list:
        for filename in file_input_list:
            the_job.uses(filename, link=Link.INPUT, transfer=True, register=True)

    if file_output_transfer:
        for filename in file_output_transfer:
            the_job.uses(filename, link=Link.OUTPUT, transfer=True)
    if file_output_notransfer:
        for filename in file_output_notransfer:
            the_job.uses(filename, link=Link.OUTPUT, transfer=False)
    return the_job

if __name__ == '__main__':
    ap = ArgumentParser()
    ap.add_argument("-i", "--input_file", type=str, required=True,
                    help="the path to the input file.")
    ap.add_argument("-o", "--output_file", type=str, required=True,
                    help="the path to the output file.")
    ap.add_argument("-s", "--source_code_dir", type=str, 
            default=os.path.expanduser('~/src/mygit/'), 
            help="the path to the source code dir. (default: %(default)s)")
    args = ap.parse_args()
    workflow_AC = ADAG("pegasus_test")
    input_file = registerOneInputFile(workflow=workflow_AC, inputFname=args.input_file, 
        input_site_handler=site_handler, \
        checkFileExistence=True)
    
    select_gene = registerExecutefile(workflow_AC, os.path.join(args.source_code_dir, "select_gene.py"))
    back_test = registerExecutefile(workflow_AC, os.path.join(args.source_code_dir, "back_test.py"))
    #input_path = "/simm/home/wzx/scratch/num.txt"
    #input_file = registerFile(workflow_AC, input_path)
    
    for i in range(30):
        select_gene_out_fpath = 'select_gene_outfile_de' + str(i) +'.pkl'
        select_gene_outfile = File(select_gene_out_fpath)
        select_gene_job = addJob2workflow(workflow=workflow_AC, excute=select_gene,
                                    file_input_list=[input_file],
                                    file_output_transfer=[select_gene_outfile],
                                    file_output_notransfer=None,
            argv=['-i', input_file, '-n %s'%(i), '-t 1 -o ', select_gene_outfile])
        setJobToProperMemoryRequirement(job=select_gene_job)
        workflow_AC.addJob(select_gene_job)

        backtest_out_fpath = 'backtest_outfile_de' + str(i)+'.csv'
        backtest_outfile = File(backtest_out_fpath)
        back_test_job = addJob2workflow(workflow=workflow_AC, excute=back_test,
                                    file_input_list=[select_gene_outfile],
                                    file_output_transfer=[backtest_outfile],
                                    file_output_notransfer=None,
                                    argv=[select_gene_outfile, backtest_outfile])
        setJobToProperMemoryRequirement(job=back_test_job, job_max_memory=2048)
        workflow_AC.addJob(back_test_job)
        workflow_AC.addDependency(Dependency(parent=select_gene_job, child=back_test_job))
    workflow_AC.writeXML(open(args.output_file, 'w'))
