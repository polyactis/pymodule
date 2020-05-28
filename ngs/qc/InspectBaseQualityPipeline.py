#!/usr/bin/env python3
"""
A pegasus workflow to run InspectBaseQuality.py over a list of individual sequences

Examples:
    # on condorpool
    %s -o inspectBaseQuality.xml -u yh -i 1-8,15-130
    
    # use hoffman2 site_handler
    %s -o inspectBaseQuality.xml -u yh -i 1-8,15-130 
        -l hoffman2 --home_path /u/home/eeskin/polyacti
        -t ~/NetworkData/vervet/db

"""
import sys, os, copy
__doc__ = __doc__%(sys.argv[0], sys.argv[0])
import getpass
from palos import ProcessOptions, getListOutOfStr, PassingData, utils
from palos.db import SunsetDB
from palos.ngs.AbstractNGSWorkflow import AbstractNGSWorkflow
from pegaflow.DAX3 import File, PFN

class InspectBaseQualityPipeline(AbstractNGSWorkflow):
    __doc__ = __doc__
    
    def __init__(self,
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
        AbstractNGSWorkflow.__init__(self,
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
        
        self.ind_seq_id_ls = getListOutOfStr(ind_seq_id_ls, data_type=int)
        self.needSplitChrIntervalData = False

    def registerExecutables(self):
        """
        """
        AbstractNGSWorkflow.registerExecutables(self)
        self.registerOneExecutable(
            path=os.path.join(self.pymodulePath, 'ngs/qc/InspectBaseQuality.py'),
            name='InspectBaseQuality', clusterSizeMultiplier=1)
        
    def run(self):
        """
        """
        self.setup_run()
        
        logDir = "Log"
        logDirJob = self.addMkDirJob(outputDir=logDir)
        individualSequenceID2FilePairLs = self.db_main.getIndividualSequenceID2FilePairLs(
            self.ind_seq_id_ls, data_dir=self.data_dir)
        for ind_seq_id, FilePairLs in individualSequenceID2FilePairLs.items():
            individual_sequence = self.db_main.queryTable(
                SunsetDB.IndividualSequence).get(ind_seq_id)
            if individual_sequence is None or individual_sequence.format!='fastq':
                continue
            for filePair in FilePairLs:
                for fileRecord in filePair:
                    relativePath = fileRecord[0]
                    prefix, suffix = utils.getRealPrefixSuffix(
                        os.path.basename(relativePath))
                    if suffix=='.fastq':
                        filepath = os.path.join(self.data_dir, relativePath)
                        #Do not register the input fastq because InspectBaseQuality
                        #  will output directly into self.data_dir.
                        logFile = File(os.path.join(logDir, f'{prefix}.log'))
                        job = self.addDBJob(
                            executable=self.InspectBaseQuality,
                            outputArgumentOption="--logFilename",
                            outputFile=logFile,
                            extraArgumentList=[
                                '-i', filepath,
                                '--read_sampling_rate', '0.01',
                                '--quality_score_format',
                                individual_sequence.quality_score_format,
                                ],
                            parentJobLs=[logDirJob],
                            transferOutput=True,
                            objectWithDBArguments=self,
                            job_max_memory=5000,
                            walltime=120)
        
        self.end_run()
    
if __name__ == '__main__':
    from argparse import ArgumentParser
    ap = ArgumentParser()
    ap.add_argument("--drivername", default="postgresql",
        help='Type of database server (default: %(default)s)')
    ap.add_argument('-z', "--hostname", default="pdc",
        help='name/IP of database server (default: %(default)s)')
    ap.add_argument("--port", default=None,
        help="database server port (default: %(default)s)")
    ap.add_argument("--dbname", default='pmdb',
        help="database name (default: %(default)s)")
    ap.add_argument('-k', "--schema", default='sunset', 
        help="database schema (default: %(default)s)")
    ap.add_argument("-u", "--db_user", help="Database user")
    ap.add_argument("-p", "--db_passwd", required=False,
        help="Password of the database user")

    ap.add_argument('-a', "--ref_ind_seq_id", type=int, required=True,
        help="Select this (IndividualSequence.id) as the reference")

    ap.add_argument('-i', "--ind_seq_id_ls", required=True,
        help='a comma/dash-separated list of IndividualSequence.id.')
    ap.add_argument("--excludeContaminant", action='store_true',
        help='Toggle to exclude sequences from contaminated individuals, '
            '(IndividualSequence.is_contaminated=1).')
    ap.add_argument("--sequence_filtered", type=int, default=None,
        help='Filter with individual_sequence.filtered=THIS_VALUE. '
            'None: everything; 0: unfiltered sequences; 1: filtered sequences.')
    
    ap.add_argument("-F", "--pegasusFolderName", default='input',
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
    
    ap.add_argument("--home_path",
        help="Path to your home folder. Default is ~.")
    ap.add_argument("--javaPath", default='bin/java',
        help="Path to java. Default is %(default)s.")
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
    if not args.db_user:
        args.db_user = getpass.getuser()
    if not args.db_passwd:
        args.db_passwd = getpass.getpass(f"Password for {args.db_user}:")
    
    instance = InspectBaseQualityPipeline(
        drivername = args.drivername,
        hostname = args.hostname,
        port = args.port,
        dbname = args.dbname, schema = args.schema,
        db_user = args.db_user, db_passwd = args.db_passwd,

        ref_ind_seq_id = args.ref_ind_seq_id,

        ind_seq_id_ls = args.ind_seq_id_ls,
        excludeContaminant = args.excludeContaminant,
        sequence_filtered = args.sequence_filtered,

        site_handler = args.site_handler, 
        input_site_handler = args.input_site_handler,
        cluster_size = args.cluster_size,
        pegasusFolderName = args.pegasusFolderName,
        output_path = args.output_path,
        tmpDir = args.tmpDir,
        max_walltime = args.max_walltime,

        home_path = args.home_path,
        javaPath = args.javaPath,
        pymodulePath = args.pymodulePath,
        
        needSSHDBTunnel = args.needSSHDBTunnel,
        commit = args.commit,
        debug = args.debug,
        report=args.report)
    instance.run()