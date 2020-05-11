#!/usr/bin/env python3
"""
A workflow that counts reads and bases in fastq/fasta files.

Examples:
    # 2012.5.3 run on hoffman2's condorpool, need sshDBTunnel (-H1)
    %s  -i 963-1346 -o dags/ReadCount/read_count_isq_936_1346.xml
        -u yh --commit -z localhost
        --pegasusFolderName readcount --needSSHDBTunnel
        -l hcondor -j hcondor -D NetworkData/vervet/db/
    
    # 2012.3.14 
    %s -i 1-864 -o dags/ReadCount/read_count_isq_1_864.xml
        -u yh -l condor -j condor -z uclaOffice
        --pegasusFolderName readCount --commit
    
"""
import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0])

import subprocess, copy
from pegaflow.DAX3 import File, Link, PFN, Job
from palos import ProcessOptions, getListOutOfStr, PassingData, utils
from palos.ngs.AbstractNGSWorkflow import AbstractNGSWorkflow as ParentClass
from palos.db import SunsetDB
import logging

class CountReadsWorkflow(ParentClass):
    __doc__ = __doc__
    option_default_dict = copy.deepcopy(ParentClass.option_default_dict)
    option_default_dict.update({
        ('ind_seq_id_ls', 1, ): ['', 'i', 1, \
            'a comma/dash-separated list of IndividualSequence.id.'
            'non-fastq entries will be discarded.', ],\
        })

    def __init__(self, ind_seq_id_ls=None, 
        drivername='postgresql', hostname='localhost',
        dbname='', schema='public', port=None,
        db_user=None, db_passwd=None,
        data_dir=None, local_data_dir=None,
        excludeContaminant=False,
        sequence_filtered=None,
        site_id_ls=None,
        country_id_ls=None,
        tax_id_ls="9606",
        sequence_type_id_ls=None,
        sequencer_id_ls=None,
        sequence_batch_id_ls=None,
        version_ls=None,
        site_handler='condor', input_site_handler='condor', cluster_size=30,
        pegasusFolderName='folder', output_path=None,
        tmpDir='/tmp/', max_walltime=4320,
        home_path=None, 
        pymodulePath="src/pymodule",
        needSSHDBTunnel=False, commit=False,
        debug=False, report=False):
        """
        """
        self.ind_seq_id_ls = getListOutOfStr(ind_seq_id_ls, data_type=int)
        ParentClass.__init__(self, 
            inputSuffixList=None, 
            drivername=drivername, hostname=hostname,
            dbname=dbname, schema=schema, port=port,
            db_user=db_user, db_passwd=db_passwd,
            data_dir=data_dir, local_data_dir=local_data_dir,
            ref_ind_seq_id=None,
            samtools_path=None,
            picard_dir=None,
            gatk_path=None,
            gatk2_path=None,
            picard_path=None,
            tabixPath=None,
            vcftoolsPath=None,
            ligateVcfPerlPath=None,
            maxContigID=None,
            minContigID=None,
            contigMaxRankBySize=None,
            contigMinRankBySize=None,
            needFastaIndexJob=False,
            needFastaDictJob=False,
            excludeContaminant=excludeContaminant,
            sequence_filtered=sequence_filtered,
            site_id_ls=site_id_ls,
            country_id_ls=country_id_ls,
            tax_id_ls=tax_id_ls,
            sequence_type_id_ls=sequence_type_id_ls,
            sequencer_id_ls=sequencer_id_ls,
            sequence_batch_id_ls=sequence_batch_id_ls,
            version_ls=version_ls,
            sequence_max_coverage=None,
            sequence_min_coverage=None,
            site_handler=site_handler,
            pegasusFolderName=pegasusFolderName,
            output_path=output_path,
            input_site_handler=input_site_handler,
            cluster_size=cluster_size,
            tmpDir=tmpDir,
            max_walltime=max_walltime, 
            home_path=home_path,
            javaPath=None,
            pymodulePath=pymodulePath,
            plinkPath=None,
            needSSHDBTunnel=needSSHDBTunnel, commit=commit,
            debug=debug, report=report)

    def registerCustomExecutables(self):
        """
        2012.3.14
        """
        ParentClass.registerCustomExecutables(self)
        
        self.registerOneExecutable(
            path=os.path.join(self.pymodulePath,
                'mapper/computer/CountFastqReadBaseCount.py'),
            name='CountFastqReadBaseCount', clusterSizeMultiplier=1)
        
        self.registerOneExecutable(
            path=os.path.join(self.pymodulePath,
                'db/import/PutReadBaseCountIntoDB.py'),
            name='PutReadBaseCountIntoDB', clusterSizeMultiplier=0.2)
        
    
    def registerISQFiles(self, db_main=None, ind_seq_id_ls=[],
        local_data_dir='', pegasusFolderName='', \
        input_site_handler='local'):
        """
        2012.3.14
        """
        print(f"Finding all ISQ-affiliated files of {len(ind_seq_id_ls)} "
            f"ind-seq entries ...", flush=True)
        returnData = PassingData(jobDataLs=[])
        Table = SunsetDB.IndividualSequence
        query = db_main.queryTable(Table).filter(Table.id.in_(ind_seq_id_ls))
        individual_sequence_id_set = set()
        missed_individual_sequence_id_set = set()
        for individual_sequence in query:
            if individual_sequence.individual_sequence_file_ls:
                for individual_sequence_file in individual_sequence.individual_sequence_file_ls:
                    absPath = os.path.join(local_data_dir, 
                        individual_sequence_file.path)
                    if os.path.isfile(absPath):
                        inputF = File(os.path.join(pegasusFolderName,
                            individual_sequence_file.path))
                        inputF.addPFN(PFN("file://" + absPath, input_site_handler))
                        inputF.absPath = absPath
                        self.addFile(inputF)
                        returnData.jobDataLs.append(PassingData(
                            output=inputF, jobLs=[],
                            isq_id=individual_sequence.id,
                            isqf_id=individual_sequence_file.id))
                        individual_sequence_id_set.add(individual_sequence.id)
                    else:
                        missed_individual_sequence_id_set.add(individual_sequence.id)
                        logging.warn("IndividualSequenceFile.id=%s (isq-id=%s) doesn't have "
                            "any affiliated IndividualSequenceFile entries "
                            "while its path %s is not a file."%\
                            (individual_sequence_file.id, individual_sequence.id, absPath))
            elif individual_sequence.path:
                absPath = os.path.join(local_data_dir, individual_sequence.path)
                if os.path.isfile(absPath):
                    inputF = File(os.path.join(pegasusFolderName, individual_sequence.path))
                    inputF.addPFN(PFN("file://" + absPath, input_site_handler))
                    inputF.absPath = absPath
                    self.addFile(inputF)
                    returnData.jobDataLs.append(PassingData(output=inputF, 
                        jobLs=[], isq_id=individual_sequence.id,\
                        isqf_id=None))
                    individual_sequence_id_set.add(individual_sequence.id)
                else:
                    logging.warn("IndividualSequence.id=%s doesn't have any affiliated "
                        "IndividualSequenceFile entries while its path %s is not a file."%\
                        (individual_sequence.id, absPath))
                    missed_individual_sequence_id_set.add(individual_sequence.id)
        
        print(f" {len(returnData.jobDataLs)} files registered for "
            f"{len(individual_sequence_id_set)} individual_sequence entries. "
            f"Missed {len(missed_individual_sequence_id_set)} "
            f"individual-sequence entries.",
            flush=True)
        return returnData
    
    def addPutReadBaseCountIntoDBJob(self, executable=None, inputFileLs=[], \
        logFile=None, commit=False, parentJobLs=[], extraDependentInputLs=[], \
        transferOutput=True, extraArguments=None, \
        job_max_memory=10, sshDBTunnel=1, **keywords):
        """
        20170502 use addData2DBJob()
        2012.5.3
            add argument sshDBTunnel
        2012.3.14
        """
        job = self.addData2DBJob(executable=executable, \
            inputFile=None, inputArgumentOption="-i", \
            outputFile=logFile, outputArgumentOption="--logFilename",
            inputFileList=inputFileLs, \
            data_dir=None, commit=commit,\
            parentJobLs=parentJobLs, extraDependentInputLs=extraDependentInputLs, \
            extraOutputLs=None, transferOutput=transferOutput, \
            extraArguments=extraArguments, extraArgumentList=None, \
            job_max_memory=job_max_memory,  sshDBTunnel=sshDBTunnel,\
            key2ObjectForJob=None, objectWithDBArguments=self, **keywords)
        return job
    
    
    def addCountFastqReadBaseCountJob(self, executable=None, inputFile=None, \
        outputFile=None, isq_id=None, isqf_id=None, \
        parentJobLs=None, extraDependentInputLs=None, transferOutput=True,
        extraArguments=None, \
        job_max_memory=100, **keywords):
        """
        20170503 use addGenericJob()
        2012.3.14
        """
        job = Job(namespace=self.namespace, name=executable.name, version=self.version)
        job.addArguments("--inputFname", inputFile, "--outputFname", outputFile)
        if not extraArguments:
            extraArguments = ""
        if isq_id:
            extraArguments += " --isq_id %s "%(isq_id)
        if isqf_id:
            extraArguments += " --isqf_id %s "%(isqf_id)
        
        job = self.addGenericJob(executable=executable, \
            inputFile=inputFile, \
            inputArgumentOption="-i", \
            outputFile=outputFile, outputArgumentOption="-o", \
            parentJobLs=parentJobLs, \
            extraDependentInputLs=extraDependentInputLs, \
            transferOutput=transferOutput, \
            extraArguments=extraArguments, \
            job_max_memory=job_max_memory, \
            **keywords)
        return job
    
    def addJobs(self, inputData=None, pegasusFolderName="", needSSHDBTunnel=0):
        """
        2012.3.14
        """
        
        sys.stderr.write("Adding read counting jobs on %s input ..."%\
            (len(inputData.jobDataLs)))
        returnJobData = PassingData()
        
        no_of_jobs = 0
        
        topOutputDir = pegasusFolderName
        if topOutputDir:
            topOutputDirJob = self.addMkDirJob(outputDir=topOutputDir)
            no_of_jobs += 1
        else:
            topOutputDirJob = None
        
        finalReduceFile = File(os.path.join(topOutputDir, 'read_base_count.tsv'))
        
        readBaseCountMergeJob = self.addStatMergeJob(
            statMergeProgram=self.mergeSameHeaderTablesIntoOne, \
            outputF=finalReduceFile, transferOutput=True, extraArguments=None,
            parentJobLs=[topOutputDirJob])
        
        logFile = File(os.path.join(topOutputDir, 'PutReadBaseCountIntoDB.log'))
        putCountIntoDBJob = self.addPutReadBaseCountIntoDBJob(
            executable=self.PutReadBaseCountIntoDB,
            inputFileLs=[finalReduceFile],
            logFile=logFile, commit=self.commit,
            parentJobLs=[readBaseCountMergeJob],
            extraDependentInputLs=[],
            transferOutput=True,
            extraArguments=None,
            job_max_memory=10, sshDBTunnel=needSSHDBTunnel)
        no_of_jobs += 2
        for jobData in inputData.jobDataLs:
            #add the read count job
            outputFile = File(os.path.join(topOutputDir, 'read_count_isq_%s_isqf_%s.tsv'%
                (jobData.isq_id, jobData.isqf_id)))
            readCountJob = self.addCountFastqReadBaseCountJob(
                executable=self.CountFastqReadBaseCount, \
                inputFile=jobData.output, outputFile=outputFile,
                isq_id=jobData.isq_id,
                isqf_id=jobData.isqf_id, \
                parentJobLs=jobData.jobLs + [topOutputDirJob],
                extraDependentInputLs=None,
                transferOutput=False, extraArguments=None,
                job_max_memory=10, no_of_cpus=4)
            
            no_of_jobs += 1
            self.addInputToStatMergeJob(statMergeJob=readBaseCountMergeJob, \
                                inputF=readCountJob.output, parentJobLs=[readCountJob])
            
        sys.stderr.write("%s jobs.\n"%(no_of_jobs))
        return putCountIntoDBJob
    
    def setup_run(self):
        """
        2013.04.07 wrap all standard pre-run() related functions into this function.
            setting up for run(), called by run()
        """
        pdata = ParentClass.setup_run(self)
        workflow = pdata.workflow
        
        db_main = self.db_main
        session = db_main.session
        session.begin(subtransactions=True)
        """
        Traceback (most recent call last):
          File "/u/home/eeskin/polyacti/script/vervet/src/db/CountReadsWorkflow.py", line 249, in <module>
            instance.run()
          File "/u/home/eeskin/polyacti/script/vervet/src/db/CountReadsWorkflow.py", line 232, in run
            pdata = self.setup_run()
          File "/u/home/eeskin/polyacti/script/vervet/src/db/CountReadsWorkflow.py", line 217, in setup_run
            session.begin()
          File "/u/home/eeskin/polyacti/lib/python/sqlalchemy/orm/scoping.py", line 139, in do
            return getattr(self.registry(), name)(*args, **kwargs)
          File "/u/home/eeskin/polyacti/lib/python/sqlalchemy/orm/session.py", line 550, in begin
            "A transaction is already begun.  Use subtransactions=True "
        sqlalchemy.exc.InvalidRequestError: A transaction is already begun. 
            Use subtransactions=True to allow subtransactions.
        """
        inputData = self.registerISQFiles(db_main=db_main, ind_seq_id_ls=self.ind_seq_id_ls, \
            local_data_dir=self.local_data_dir, pegasusFolderName=self.pegasusFolderName,\
            input_site_handler=self.input_site_handler)
        
        registerReferenceData = self.getReferenceSequence()
        return PassingData(inputData=inputData,\
            registerReferenceData=registerReferenceData)

    def run(self):
        """
        2011-7-11
        """
        pdata = self.setup_run()
        inputData = pdata.inputData
        self.addJobs(inputData=inputData, pegasusFolderName=self.pegasusFolderName,
            needSSHDBTunnel=self.needSSHDBTunnel)
        self.end_run()

if __name__ == '__main__':
    from argparse import ArgumentParser
    ap = ArgumentParser()
    ap.add_argument('-i', "--ind_seq_id_ls", required=True,
        help='a comma/dash-separated list of IndividualSequence.id.'
            'non-fastq entries will be discarded.')
    ap.add_argument("--drivername", default="postgresql",
        help='Type of database server (default: %(default)s)')
    ap.add_argument("--hostname", default="pdc",
        help='name/IP of database server (default: %(default)s)')
    ap.add_argument("--port", default=None,
        help="database server port (default: %(default)s)")
    ap.add_argument("--dbname", default='pmdb',
        help="database name (default: %(default)s)")
    ap.add_argument('-k', "--schema", default='sunset', 
        help="database schema (default: %(default)s)")
    ap.add_argument("-u", "--db_user", required=True, help="Database user")
    ap.add_argument("-p", "--db_passwd", required=False,
        help="Password of the database user")
    ap.add_argument("-F", "--pegasusFolderName", type=str,
        help='The path relative to the workflow running root. '
        'This folder will contain pegasus input & output. '
        'It will be created during the pegasus staging process. '
        'It is useful to separate multiple sub-workflows. '
        'If empty or None, everything is in the pegasus root.')
    ap.add_argument("-l", "--site_handler", type=str, required=True,
        help="The name of the computing site where the jobs run and "
        "executables are stored. "
        "Check your Pegasus configuration in submit.sh.")
    ap.add_argument("-j", "--input_site_handler", type=str,
        help="It is the name of the site that has all the input files."
        "Possible values can be 'local' or same as site_handler."
        "If not given, it is asssumed to be the same as site_handler and "
        "the input files will be symlinked into the running folder."
        "If input_site_handler=local, the input files will be transferred "
        "to the computing site by pegasus-transfer.")
    ap.add_argument("-C", "--cluster_size", type=int, default=30,
        help="Default: %(default)s. "
        "This number decides how many of pegasus jobs should be clustered "
        "into one job. "
        "Good if your workflow contains many quick jobs. "
        "It will reduce Pegasus monitor I/O.")
    ap.add_argument("-o", "--output_path", type=str, required=True,
        help="The path to the output file that will contain the Pegasus DAG.")
    
    ap.add_argument("--home_path", type=str,
        help="Path to your home folder. Default is ~.")
    ap.add_argument("--pymodulePath", type=str, default="src/pymodule",
        help="Path to the pymodule code folder. "
        "If relative path, home folder is inserted in the front.")
    
    ap.add_argument("--tmpDir", type=str, default='/tmp/',
        help='Default: %(default)s. '
        'A local folder for some jobs (MarkDup) to store temp data. '
        '/tmp/ can be too small for high-coverage sequencing.')
    ap.add_argument("--max_walltime", type=int, default=4320,
        help='Default: %(default)s. '
        'Maximum wall time for any job, in minutes. 4320=3 days. '
        'Used in addGenericJob(). Most clusters have upper limit for runtime.')
    ap.add_argument("--needSSHDBTunnel", action='store_true',
        help="If all DB-interacting jobs need a ssh tunnel to "
        "access a database that is inaccessible to computing nodes.")
    ap.add_argument("-c", "--commit", action='store_true',
        help="Toggle to commit the db transaction (default: %(default)s)")
    ap.add_argument("--debug", action='store_true',
        help='Toggle debug mode.')
    ap.add_argument("--report", action='store_true',
        help="Toggle verbose mode. Default: %(default)s.")

    args = ap.parse_args()
    instance = CountReadsWorkflow(
        ind_seq_id_ls=args.ind_seq_id_ls,

        drivername=args.drivername, 
	    hostname=args.hostname,
        port=args.port,
        dbname=args.dbname, schema=args.schema,
        db_user=args.db_user, db_passwd=args.db_passwd,

        pegasusFolderName=args.pegasusFolderName,
        site_handler=args.site_handler, 
        input_site_handler=args.input_site_handler,
        cluster_size=args.cluster_size,
        output_path=args.output_path,

        home_path=args.home_path,
        pymodulePath=args.pymodulePath,
        
        tmpDir=args.tmpDir,
        max_walltime=args.max_walltime,
        needSSHDBTunnel=args.needSSHDBTunnel,
        commit=args.commit,
        debug=args.debug,
        report=args.report)
    instance.run()
