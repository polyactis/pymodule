#!/usr/bin/env python3

"""
Examples:
    %s -n 130 -u yh 1.fastq.gz 2.fastq.gz

    %s

Description:
    2011-8-16
        script that generates a read-filtering pegasus workflow.
        It uses picard_path/FilterRead.jar and AddFilteredSequences2DB.py.

    1. Be careful with the db connection setting as it'll be passed to the db-registration job.
        Make sure all computing nodes have access to the db.
    2. The workflow has to be run on nodes where they have direct db and db-affiliated file-storage access. 
"""

import sys
import os
__doc__ = __doc__ % (sys.argv[0], sys.argv[0])


from palos import getListOutOfStr
from palos.db import SunsetDB
from palos.ngs.AbstractNGSWorkflow import AbstractNGSWorkflow
from pegaflow.api import File
import getpass


class FilterShortReadPipeline(AbstractNGSWorkflow):
    __doc__ = __doc__

    def __init__(
        self, adapter, adapter2, minimum_length, maximum_length, trim_n,
        quality_cutoff, quality_base, no_of_threads, cutadapt_path,
        drivername='postgresql', hostname='localhost',
        dbname='', schema='public', port=None,
        db_user=None, db_passwd=None,

        data_dir=None, local_data_dir=None,

        ref_ind_seq_id=None,

        ind_seq_id_ls=None,
        excludeContaminant=False,
        sequence_filtered=None,

        samtools_path="bin/samtools",
        picard_dir="script/picard/dist",
        gatk_path="bin/GenomeAnalysisTK1_6_9.jar",
        gatk2_path="bin/GenomeAnalysisTK.jar",
        picard_path="script/picard.broad/build/libs/picard.jar",

        chromosome_type_id=None, 
        ref_genome_tax_id=9606,
        ref_genome_sequence_type_id=1,
        ref_genome_version=15,
        ref_genome_outdated_index=0,

        pegasusFolderName='input',
        site_handler='condor',
        input_site_handler='condor',
        cluster_size=30,
        output_path=None,
        tmpDir='/tmp/',
        max_walltime=4320,
        
        home_path=None,
        javaPath=None,
        pymodulePath="src/pymodule",

        jvmVirtualByPhysicalMemoryRatio=1.2,

        needSSHDBTunnel=False,
        commit=False,
        debug=False,
        report=False):
        """
        """
        self.adapter = adapter
        self.adapter2 = adapter2
        self.minimum_length = minimum_length
        self.maximum_length = maximum_length
        self.trim_n = trim_n
        self.quality_cutoff = quality_cutoff
        self.quality_base = quality_base
        self.no_of_threads = no_of_threads
        self.cutadapt_path = cutadapt_path
        AbstractNGSWorkflow.__init__(
            self,
            drivername=drivername, hostname=hostname,
            dbname=dbname, schema=schema, port=port,
            db_user=db_user, db_passwd=db_passwd,

            data_dir=data_dir, local_data_dir=local_data_dir,

            ref_ind_seq_id=ref_ind_seq_id,

            excludeContaminant=excludeContaminant,
            sequence_filtered=sequence_filtered,

            samtools_path=samtools_path,
            picard_dir=picard_dir,
            gatk_path=gatk_path,
            gatk2_path=gatk2_path,
            picard_path=picard_path,

            contigMaxRankBySize=None,
            contigMinRankBySize=None,

            chromosome_type_id=chromosome_type_id,
            ref_genome_tax_id=ref_genome_tax_id,
            ref_genome_sequence_type_id=ref_genome_sequence_type_id,
            ref_genome_version=ref_genome_version,
            ref_genome_outdated_index=ref_genome_outdated_index,

            pegasusFolderName=pegasusFolderName,
            site_handler=site_handler,
            input_site_handler=input_site_handler,
            cluster_size=cluster_size,
            output_path=output_path,
            tmpDir=tmpDir,
            max_walltime=max_walltime,

            home_path=home_path,
            javaPath=javaPath,
            pymodulePath=pymodulePath,

            jvmVirtualByPhysicalMemoryRatio=jvmVirtualByPhysicalMemoryRatio,

            needSSHDBTunnel=needSSHDBTunnel,
            commit=commit,
            debug=debug, report=report)
        """
        2021-7-12
        """
        self.ind_seq_id_ls = getListOutOfStr(ind_seq_id_ls, data_type=int)

    def registerExecutables(self):
        """
        2021-7-12
        """
        AbstractNGSWorkflow.registerExecutables(self)
        self.registerOneExecutable(
            path=os.path.join(self.pymodulePath,
                              "db/import/AddFilteredSequences2DB.py"),
            name='AddFilteredSequences2DB', clusterSizeMultipler=0.5)
        self.registerOneExecutable(path=self.cutadapt_path, name='cutadapt',
                                   clusterSizeMultipler=0.6)

    def addAddFilteredSequences2DB_job(
         self, executable=None, inputFile=None, individual_sequence_id=None,
         outputDir=None, logFile=None, parent_individual_sequence_file_id=None,
         parentJobLs=None, job_max_memory=100, walltime=60, commit=0,
         extraDependentInputLs=None, transferOutput=False,
         sshDBTunnel=1, **keywords):

        """
        """
        extraArgumentList = [
            f"--individual_sequence_id {individual_sequence_id}",
            f'--outputDir {outputDir}',
            f'--parent_individual_sequence_file_id {parent_individual_sequence_file_id}']
        job = self.addData2DBJob(
            executable=executable, inputFile=inputFile,
            inputArgumentOption="-i", outputFile=None,
            outputArgumentOption="-o", inputFileList=None, data_dir=None,
            logFile=logFile, commit=commit, parentJobLs=parentJobLs,
            extraDependentInputLs=extraDependentInputLs, extraOutputLs=None,
            transferOutput=transferOutput, extraArguments=None,
            extraArgumentList=extraArgumentList, job_max_memory=job_max_memory,
            sshDBTunnel=sshDBTunnel, walltime=walltime, key2ObjectForJob=None,
            objectWithDBArguments=self, **keywords)
        return job

    def addFilterReadJob(
         self, executable=None, parentJobLs=None, extraOutputLs=None,
         job_max_memory=2000, walltime=120, extraDependentInputLs=None,
         transferOutput=False, extraArgumentList=None, no_of_cpus=20):
        """
        """
        job = self.addGenericJob(
            executable=executable, extraOutputLs=extraOutputLs,
            transferOutput=transferOutput, parentJobLs=parentJobLs,
            extraDependentInputLs=extraDependentInputLs,
            extraArgumentList=extraArgumentList, job_max_memory=job_max_memory,
            no_of_cpus=no_of_cpus, walltime=walltime)
        return job

    def getLibrarySplitOrder2DBEntryLs(self, individual_sequence):
        """
        2012.2.10
            generate a dictionary for this program to check if some
            (library, split_order)s have been filtered and recorded in db already
        """
        library_split_order2db_entry_ls = {}
        for individual_sequence_file in individual_sequence.individual_sequence_file_ls:
            key = (individual_sequence_file.library, individual_sequence_file.split_order)
            if key not in library_split_order2db_entry_ls:
                library_split_order2db_entry_ls[key] = []
            library_split_order2db_entry_ls[key].append(individual_sequence_file)
        return library_split_order2db_entry_ls

    def run(self):
        """
        """
        self.setup_run()

        isq_id2LibrarySplitOrder2FileLs = self.db_main.getISQ_ID2LibrarySplitOrder2FileLs(
            self.ind_seq_id_ls, data_dir=self.data_dir,
            filtered=0, ignoreEmptyReadFile=False)
        to_work_ind_seq_id_set = set()
        parent_individual_sequence_file_id_set = set()
        for ind_seq_id, LibrarySplitOrder2FileLs in isq_id2LibrarySplitOrder2FileLs.items():
            parent_individual_sequence = self.db_main.queryTable(
                SunsetDB.IndividualSequence).get(ind_seq_id)
            if parent_individual_sequence is not None and parent_individual_sequence.format == 'fastq':
                """
                check if the child individual_sequence already exists in db or not.
                if it does, what about its files?? if not, go add filtering jobs.
                """
                # 2012.6.8
                individual_sequence = self.db_main.copyParentIndividualSequence(
                    parent_individual_sequence=parent_individual_sequence,
                    parent_individual_sequence_id=ind_seq_id,
                    quality_score_format='Standard', filtered=1,
                    data_dir=self.data_dir)
                library_split_order2filtered_db_entry_ls = self.getLibrarySplitOrder2DBEntryLs(individual_sequence)

                sequenceOutputDirJob = None
                filteredReadOutputDirJob = None
                for key, fileObjLs in LibrarySplitOrder2FileLs.items():
                    if key in library_split_order2filtered_db_entry_ls:
                        sys.stderr.write(
                            "Warning: this pair of filtered individual_sequence_file(s), "
                            f"{repr(key)}, parent_individual_sequence "
                            f"(id={parent_individual_sequence.id}, {parent_individual_sequence.individual.code}), "
                            f"individual_sequence (id={individual_sequence.id}, {individual_sequence.individual.code}) "
                            "are already in db. skip.\n")
                        continue
                    else:
                        if sequenceOutputDirJob is None:
                            sequenceOutputDir = os.path.join(
                                self.data_dir, individual_sequence.path)
                            sequenceOutputDirJob = self.addMkDirJob(
                                outputDir=sequenceOutputDir)
                        if filteredReadOutputDirJob is None:
                            filteredReadOutputDir = os.path.join(
                                os.path.basename(individual_sequence.path))
                            filteredReadOutputDirJob = self.addMkDirJob(
                                outputDir=filteredReadOutputDir)

                    # add filter jobs
                    extraDependentInputLs = []
                    extraOutputLs = []
                    extraArgumentList = [
                        "-a", self.adapter, "-j", self.no_of_threads,
                        "--quality-base", self.quality_base,
                        "-m", self.minimum_length]
                    if self.adapter2 is not None:
                        extraArgumentList.extend(["-A", self.adapter2])
                    if self.maximum_length is not None:
                        extraArgumentList.extend(["-M", self.maximum_length])
                    if self.trim_n:
                        extraArgumentList.append("--trim-n")
                    if self.quality_cutoff is not None:
                        extraArgumentList.extend(["-q", self.quality_cutoff])
                    input_fastq_list = []
                    for i in range(len(fileObjLs)):
                        fileObj = fileObjLs[i]
                        try:  # 2012.7.2
                            inputFile = self.registerOneInputFile(
                                input_path=fileObj.path,
                                folderName='inputIndividualSequenceFile')
                        except Exception as e:
                            import pdb
                            pdb.set_trace()
                        # take the base filename as the output filename. it'll be in scratch/.
                        outputFname = os.path.join(
                            filteredReadOutputDir,
                            os.path.basename(fileObj.path))
                        outputFile = File(outputFname)
                        extraDependentInputLs.append(inputFile)
                        extraOutputLs.append(outputFile)
                        if i == 0:   # 1st mate
                            input_fastq_list.append(inputFile)
                            extraArgumentList.extend(["-o", outputFile])
                        elif i == 1:    # 2nd mate
                            input_fastq_list.append(inputFile)
                            extraArgumentList.extend(["-p", outputFile])
                        else:
                            sys.stderr.write("Error: mate %s appeared in paired-end data (individualSequenceID=%s).\n"%(i+1, ind_seq_id))
                            sys.exit(4)

                    extraArgumentList.extend(input_fastq_list)
                    filterShortRead_job = self.addFilterReadJob(
                        executable=self.cutadapt, extraOutputLs=extraOutputLs,
                        parentJobLs=[filteredReadOutputDirJob],
                        job_max_memory=2000, walltime=120,
                        extraDependentInputLs=extraDependentInputLs,
                        extraArgumentList=extraArgumentList,
                        no_of_cpus=self.no_of_threads,
                        transferOutput=False)
                    for fileObj, outputFile in zip(fileObjLs, extraOutputLs):
                        logFile = File('%s_%s.register.log'%(individual_sequence.id, fileObj.db_entry.id))
                        addFilteredSequences2DB_job = self.addAddFilteredSequences2DB_job(
                            executable=self.AddFilteredSequences2DB,
                            inputFile=outputFile,
                            individual_sequence_id=individual_sequence.id,
                            outputDir=sequenceOutputDir, logFile=logFile,
                            parent_individual_sequence_file_id=fileObj.db_entry.id,
                            parentJobLs=[sequenceOutputDirJob, filterShortRead_job],
                            commit=self.commit, extraDependentInputLs=None,
                            transferOutput=True, sshDBTunnel=self.needSSHDBTunnel)
                    to_work_ind_seq_id_set.add(ind_seq_id)
                    parent_individual_sequence_file_id_set.add(fileObj.db_entry.id)
        sys.stderr.write(
            f"{self.no_of_jobs} jobs, {len(to_work_ind_seq_id_set)} individual_sequence entries, "
            f"{len(parent_individual_sequence_file_id_set)} parent_individual_sequence_file_id s.\n")

        self.end_run()


if __name__ == '__main__':
    from argparse import ArgumentParser
    ap = ArgumentParser()
    ap.add_argument(
        "--drivername", default="postgresql",
        help='Type of database server (default: %(default)s)')
    ap.add_argument(
        '-z', "--hostname", default="pdc",
        help='name/IP of database server (default: %(default)s)')
    ap.add_argument(
        "--port", default=None,
        help="database server port (default: %(default)s)")
    ap.add_argument(
        "--dbname", default='pmdb',
        help="database name (default: %(default)s)")
    ap.add_argument(
        '-k', "--schema", default='sunset', 
        help="database schema (default: %(default)s)")
    ap.add_argument("-u", "--db_user", help="Database user")
    ap.add_argument(
        "-p", "--db_passwd", required=False, help="Password of the database user")

    ap.add_argument(
        '-a', "--ref_ind_seq_id", type=int, required=True,
        help="Select this (IndividualSequence.id) as the reference")

    ap.add_argument(
        '-i', "--ind_seq_id_ls", required=True,
        help='a comma/dash-separated list of IndividualSequence.id.')
    ap.add_argument(
        "--excludeContaminant", action='store_true',
        help='Toggle to exclude sequences from contaminated individuals, '
             '(IndividualSequence.is_contaminated=1).')
    ap.add_argument(
        "--sequence_filtered", type=int, default=None,
        help='Filter with individual_sequence.filtered=THIS_VALUE. '
             'None: everything; 0: unfiltered sequences; 1: filtered sequences.')
    ap.add_argument(
        "--adapter", required=True,
        help="Sequence of an adapter ligated to the 3' end (paired data: of the first read). "
             "The adapter and subsequent bases are trimmed. If a '$' character is appended "
             "('anchoring'), the adapter is only found if it is a suffix of the read.")
    ap.add_argument(
        "-A", "--adapter2",
        help="3' adapter to be removed from second read in a pair.")
    ap.add_argument(
        "-m", "--minimum-length", default=0, type=int,
        help="When trimming paired-end reads, the minimum lengths for R1 and "
        "R2 can be specified separately by separating them with a colon (:). "
        "If the colon syntax is not used, the same minimum length applies to "
        "both reads, as discussed above. Also, one of the values can be "
        "omitted to impose no restrictions. For example, with -m 17:, "
        "the length of R1 must be at least 17, "
        "but the length of R2 is ignored. (default: %(default)s)")
    ap.add_argument(
        "-M", "--maximum-length", default=None,
        help="Maximum lengths can also be specified separately, "
        "see the explanation of -m above. (default: no limit)")
    ap.add_argument(
        "--trim-n", action='store_true', help="Trim N's on ends of reads.")
    ap.add_argument(
        "-q", "--quality-cutoff",
        help="Trim low-quality bases from 5' and/or 3' ends of each read before adapter "
             "removal. Applied to both reads if data is paired. If one value is given, "
             "only the 3' end is trimmed. If two comma-separated cutoffs are given, "
             "the 5' end is trimmed with the first cutoff, the 3' end with the second.")
    ap.add_argument(
        "--quality-base", default=33,
        help="Assume that quality values in FASTQ are encoded as ascii(quality + N). "
             "This needs to be set to 64 for some old Illumina FASTQ files. Default: 33")

    ap.add_argument(
        "-F", "--pegasusFolderName", default='input',
        help='The path relative to the workflow running root. '
             'This folder will contain pegasus input & output. '
             'It will be created during the pegasus staging process. '
             'It is useful to separate multiple sub-workflows. '
             'If empty or None, everything is in the pegasus root.')
    ap.add_argument(
        "-l", "--site_handler", type=str, required=True,
        help="The name of the computing site where the jobs run and "
             "executables are stored. "
             "Check your Pegasus configuration in submit.sh.")
    ap.add_argument(
        "-j", "--input_site_handler", type=str,
        help="It is the name of the site that has all the input files."
             "Possible values can be 'local' or same as site_handler."
             "If not given, it is asssumed to be the same as site_handler and "
             "the input files will be symlinked into the running folder."
             "If input_site_handler=local, the input files will be transferred "
             "to the computing site by pegasus-transfer.")
    ap.add_argument(
        "-C", "--cluster_size", type=int, default=30,
        help="Default: %(default)s. "
        "This number decides how many of pegasus jobs should be clustered "
        "into one job. "
        "Good if your workflow contains many quick jobs. "
        "It will reduce Pegasus monitor I/O.")
    ap.add_argument(
        "-o", "--output_path", type=str, required=True,
        help="The path to the output file that will contain the Pegasus DAG.")

    ap.add_argument("--home_path",
                    help="Path to your home folder. Default is ~.")
    ap.add_argument("--javaPath", default='bin/java',
                    help="Path to java. Default is %(default)s.")
    ap.add_argument("--pymodulePath", type=str, default="src/pymodule",
                    help="Path to the pymodule code folder. "
                    "If relative path, home folder is inserted in the front.")
    ap.add_argument("--cutadapt_path", type=str, default="/y/home/fanxp/.local/bin/cutadapt",
                    help="Path to the binary cutadapt")

    ap.add_argument(
        "--tmpDir", type=str, default='/tmp/',
        help='Default: %(default)s. '
             'A local folder for some jobs (MarkDup) to store temp data. '
             '/tmp/ can be too small for high-coverage sequencing.')
    ap.add_argument(
        "--max_walltime", type=int, default=4320,
        help='Default: %(default)s. '
        'Maximum wall time for any job, in minutes. 4320=3 days. '
        'Used in addGenericJob(). Most clusters have upper limit for runtime.')
    ap.add_argument(
        "--needSSHDBTunnel", action='store_true',
        help="If all DB-interacting jobs need a ssh tunnel to "
             "access a database that is inaccessible to computing nodes.")
    ap.add_argument(
        "-c", "--commit", action='store_true',
        help="Toggle to commit the db transaction (default: %(default)s)")
    ap.add_argument("--debug", action='store_true', help='Toggle debug mode.')
    ap.add_argument("--report", action='store_true',
                    help="Toggle verbose mode. Default: %(default)s.")
    ap.add_argument("--no_of_threads", type=int, default=4,
        help="The number of threads for each job."
            "Default: %(default)s")
    args = ap.parse_args()
    if not args.db_user:
        args.db_user = getpass.getuser()
    if not args.db_passwd:
        args.db_passwd = getpass.getpass(f"Password for {args.db_user}:")

    instance = FilterShortReadPipeline(
        adapter=args.adapter, adapter2=args.adapter2,
        minimum_length=args.minimum_length, maximum_length=args.maximum_length,
        trim_n=args.trim_n, quality_cutoff=args.quality_cutoff,
        quality_base=args.quality_base, no_of_threads=args.no_of_threads,
        cutadapt_path=args.cutadapt_path,
        drivername=args.drivername,
        hostname=args.hostname,
        port=args.port,
        dbname=args.dbname, schema=args.schema,
        db_user=args.db_user, db_passwd=args.db_passwd,

        ref_ind_seq_id=args.ref_ind_seq_id,

        ind_seq_id_ls=args.ind_seq_id_ls,
        excludeContaminant=args.excludeContaminant,
        sequence_filtered=args.sequence_filtered,

        site_handler=args.site_handler,
        input_site_handler=args.input_site_handler,
        cluster_size=args.cluster_size,
        pegasusFolderName=args.pegasusFolderName,
        output_path=args.output_path,
        tmpDir=args.tmpDir,
        max_walltime=args.max_walltime,

        home_path=args.home_path,
        javaPath=args.javaPath,
        pymodulePath=args.pymodulePath,

        needSSHDBTunnel=args.needSSHDBTunnel,
        commit=args.commit,
        debug=args.debug,
        report=args.report)
    instance.run()
