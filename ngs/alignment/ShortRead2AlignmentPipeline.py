#!/usr/bin/env python3
"""
Examples:
    # 2011-8-30 workflow on condor, always commit (--commit)
    %s --ind_seq_id_ls 165-167 -o ShortRead2Alignment_isq_id_165_167_vs_9.xml -u yh -a 9 -l condorpool
        -n1 -z dl324b-1.cmb.usc.edu --commit --needSSHDBTunnel

    # 2011-8-30 a workflow with 454 long-read and short-read PE. need a ref index job (-n1).
    %s --ind_seq_id_ls 165-167 -o ShortRead2Alignment_isq_id_165_167_vs_9.xml -u yh -a 9
        -e /u/home/eeskin/polyacti -l hoffman2 --data_dir NetworkData/vervet/db -n1
        -z dl324b-1.cmb.usc.edu --commit
        --tmpDir /work/ --needSSHDBTunnel

    # 2011-8-30 output a workflow to run alignments on hoffman2's condor pool
    #  (--local_data_dir changes local_data_dir. --data_dir changes data_dir.)
    # 2012.3.20 use /work/ or /u/scratch/p/polyacti/tmp as TMP_DIR for MarkDuplicates.jar (/tmp is too small for 30X genome)
    # 2012.5.4 cluster 4 alignment jobs (before merging) as a unit
    #   (--alignmentJobClustersSizeFraction 0.2), skip done alignment (--skipDoneAlignment)
    # 2012.9.21 add "--needSSHDBTunnel" because AddAlignmentFile2DB need db conneciton
    # 2012.9.21 add "--alignmentPerLibrary" to also get alignment for each library within one individual_sequence
    # 2013.3.15 add "--coreAlignmentJobWallTimeMultiplier 0.5" to reduce wall time for core-alignment (bwa/stampy) jobs by half
    ref=3280; %s --ind_seq_id_ls 632-3230 --sequence_min_coverage 15 --sequence_max_coverage 80
         --site_id_ls 447 --sequence_filtered 1
        --excludeContaminant -a $ref -o dags/ShortRead2AlignmentPipeline_VRCPart1_vs_$ref\_AlnMethod2.xml
        -u yh -l hcondor -j hcondor -z localhost -u yh --commit --tmpDir /work/
        --home_path /u/home/eeskin/polyacti --no_of_aln_threads 1 --skipDoneAlignment
        -D /u/home/eeskin/polyacti/NetworkData/vervet/db/ -t NetworkData/vervet/db/
        --clusters_size 20 --alignment_method_name bwaShortRead
        --coreAlignmentJobWallTimeMultiplier 0.5
        --alignmentJobClustersSizeFraction 0.2
        --needSSHDBTunnel --ref_genome_version 2 --needRefIndexJob --db_passwd secret
        #--alignmentPerLibrary

    # 2011-8-30 a workflow to run on condorpool, no ref index job. Note the site_handler and input_site_handler are both condorpool
    # to enable symlink of input files. need ref index job (--needRefIndexJob).
    # If input_site_handler is "local", pegasus will report error saying it doesn't know how to replica-transfer input files.
    %s --ind_seq_id_ls 176,178-183,207-211
        -o ShortRead2Alignment_8VWP_vs_9_condor_no_refIndex.xml
        -u yh -a 9 -j condorpool -l condorpool --needRefIndexJob -z dl324b-1.cmb.usc.edu -p secret  --commit --needSSHDBTunnel

    # 2011-8-30 a workflow to run on condorpool, no ref index job. Note the site_handler and input_site_handler are both condorpool
    # to enable symlink of input files.
    # If input_site_handler is "local", pegasus will report error saying it doesn't know how to replica-transfer input files.
    %s --ind_seq_id_ls 176,178-183,207-211
        -o ShortRead2Alignment_8VWP_vs_9_condor_no_refIndex.xml
        -u yh -a 9 -j condorpool -l condorpool --needRefIndexJob -z dl324b-1.cmb.usc.edu -p secret  --commit --needSSHDBTunnel

    # 2011-8-30 a workflow to run on uschpc, with ref index job. Note the site_handler and input_site_handler.
    # to enable replica-transfer.
    %s --ind_seq_id_ls 391-397,456,473,493
        -o ShortRead2Alignment_4DeepVRC_6LowCovVRC_392_397_vs_9_uschpc.xml -u yh -a 9
        -j local -l uschpc -n1 -e /home/cmb-03/mn/yuhuang -z 10.8.0.10 -p secret  --commit --needSSHDBTunnel

    # 2011-8-30 a workflow to run on uschpc, Need ref index job (--needRefIndexJob), and 4 threads for each alignment job
    # Note the site_handler, input_site_handler and "--data_dir ..." to enable symlink
    %s --ind_seq_id_ls 391-397,456,473,493
        -o ShortRead2Alignment_4DeepVRC_6LowCovVRC_392_397_vs_9_uschpc.xml -u yh -a 9
        -j uschpc -l uschpc --needRefIndexJob -p secret --commit --no_of_aln_threads 4 --needSSHDBTunnel
        -e /home/cmb-03/mn/yuhuang -z 10.8.0.10
        --data_dir /home/cmb-03/mn/yuhuang/NetworkData/vervet/db/ --javaPath /home/cmb-03/mn/yuhuang/bin/jdk/bin/java

    # 2011-11-16 a workflow to run on uschpc, Need ref index job (--needRefIndexJob), and 4 threads for each alignment job
    # Note the site_handler, input_site_handler. this will stage in all input and output (--notStageOutFinalOutput).
    %s --ind_seq_id_ls 391-397,456,473,493
        -o dags/ShortRead2Alignment/ShortRead2Alignment_4DeepVRC_6LowCovVRC_392_397_vs_9_local2usc.xml -u yh -a 9
        -j local -l uschpc --needRefIndexJob -p secret --commit --no_of_aln_threads 4
        -e /home/cmb-03/mn/yuhuang -z 10.8.0.10
        --javaPath /home/cmb-03/mn/yuhuang/bin/jdk/bin/java
        --needSSHDBTunnel


    #2011-9-13 no ref index job, staging input files from localhost to uschpc, stage output files back to localhost
    # modify the refFastaFile's path in xml manually
    %s --ind_seq_id_ls 1-3 -o ShortRead2Alignment_1_3_vs_524_local2uschpc.xml -u yh -a 524
        -j local -l uschpc --needRefIndexJob -p secret --commit --no_of_aln_threads 4 -e /home/cmb-03/mn/yuhuang -z 10.8.0.10
        --data_dir /Network/Data/vervet/db/
        --needSSHDBTunnel

    # 2011-8-31 output the same workflow above but for condorpool
    %s --ind_seq_id_ls 391-397,456,473,493, -o dags/ShortRead2Alignment/ShortRead2Alignment_4DeepVRC_6LowCovVRC_392_397_vs_9_condorpool.xml
        -u yh -a 9 -j condorpool -l condorpool --needRefIndexJob -z 10.8.0.10  -p secret  --commit --alignmentPerLibrary

    # 2012-4-5 new alignment method, stampy (--alignment_method_name)
    %s --ind_seq_id_ls 167,176,178,182,183,207-211,391-397,456,473,493
        -o dags/ShortRead2Alignment/ShortRead2Alignment_10VWP_4DeepVRC_6LowCovVRC_392_397_vs_508_condorpool.xml
        -u yh -a 508 -j condorpool -l condorpool -n1 -z 10.8.0.10  -p secret  --commit --alignment_method_name stampy

    # 2013.2.28 use the new alignment-method: bwaShortReadHighMismatches
    #double the core alignment (bwa aln) job walltime (=23 hrs) (--coreAlignmentJobWallTimeMultiplier) because it takes much longer
    # set max walltime for any job to be 1 day (--max_walltime 1440)
    ref=3231; %s --ind_seq_id_ls 638 -a $ref
        -o dags/ShortRead2Alignment/ShortRead2AlignmentPipeline_Aethiops_vs_$ref\_AlnMethod5.xml
        -u yh -l hcondor -j hcondor -z localhost -u yh --commit --tmpDir /work/ --home_path /u/home/eeskin/polyacti
        --no_of_aln_threads 1 --skipDoneAlignment -D /u/home/eeskin/polyacti/NetworkData/vervet/db/
        -t /u/home/eeskin/polyacti/NetworkData/vervet/db/ --clusters_size 1 --alignment_method_name bwaShortReadHighMismatches
        --coreAlignmentJobWallTimeMultiplier 2  --needSSHDBTunnel
        --max_walltime 1440

    # 2013.04.04 use the new alignment-method: bwa-mem, get rid of "-q 20" by --additionalArguments " " as mem doesn't support -q.
    # 23*0.1 hrs walltime for the core alignment (bwa mem) jobs (--coreAlignmentJobWallTimeMultiplier 0.1) because it's much faster
    # set max walltime for any job to be 1 day (--max_walltime 1440)
    ref=1;
    %s --ind_seq_id_ls 87 -a $ref --additionalArguments " "
        -o dags/ShortRead2AlignmentPipeline_Aethiops_vs_$ref\_AlnMethod6.xml
        -l hcondor -j hcondor -z pdc -u luozhihui --commit --tmpDir /tmp/
        --no_of_aln_threads 1 --skipDoneAlignment --clusters_size 1 --alignment_method_name mem
        --coreAlignmentJobWallTimeMultiplier 0.1
        --max_walltime 1440
Description:
    2013.04.07
        A program which generates a pegasus workflow dag (xml file) which does the alignment for all available sequences.
        The workflow will stage in (or symlink if site_handler and input_site_handler are same.) every input file.
        It will also stage out every output file.
        Be careful about -R, only toggle it if you know every input individual_sequence_file is not empty.
            Empty read files would fail alignment jobs and thus no final alignment for a few indivdiuals.
        Use "--alignmentJobClustersSizeFraction ..." to cluster the alignment jobs
             if the input read file is small enough (~1Million reads for bwa, ~300K for stampy).
        The arguments related to how many chromosomes/contigs do not matter unless local_realigned=1.
"""
import sys, os
__doc__ = __doc__%(sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0], \
                sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0])

import copy
from pegaflow.DAX3 import Executable, File, PFN, Profile, Namespace
from palos import ProcessOptions, getListOutOfStr, PassingData, utils
from palos.db import SunsetDB
from palos.ngs.AbstractNGSWorkflow import AbstractNGSWorkflow
from . ShortRead2AlignmentWorkflow import ShortRead2AlignmentWorkflow

ParentClass = AbstractNGSWorkflow
class ShortRead2AlignmentPipeline(ParentClass, ShortRead2AlignmentWorkflow):
    __doc__ = __doc__
    option_default_dict = copy.deepcopy(ShortRead2AlignmentWorkflow.option_default_dict)
    option_default_dict.update(ParentClass.option_default_dict.copy())

    option_default_dict.pop(('refSequenceFname', 1, ))
    option_default_dict.update({
        ("alignmentPerLibrary", 0, int): [0, '', 0, 
            'toggle to run alignment for each library of one individual_sequence'],\
        })
    option_default_dict[('local_realigned', 0, int)][0] = 0


    """
    2012.3.29
        default to stage out final output.
        Argument stageOutFinalOutput morphs into notStageOutFinalOutput.
    2011-7-11
    """
    def __init__(self,  **keywords):
        ShortRead2AlignmentWorkflow.__init__(self, **keywords)
        #ParentClass.__init__(self, **keywords)
        if self.ind_seq_id_ls:
            self.ind_seq_id_ls = getListOutOfStr(self.ind_seq_id_ls, data_type=int)

   

if __name__ == '__main__':
    main_class = ShortRead2AlignmentPipeline
    po = ProcessOptions(sys.argv, main_class.option_default_dict,
        error_doc=main_class.__doc__)
    instance = main_class(**po.long_option2value)
    instance.run()
