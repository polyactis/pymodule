#!/usr/bin/env python
"""
run Accurity in the  workflow
    ./downsample.py --outputFname dags/runDownsample_median.xml --data_dir /y/Sunset/db/ \
    --drivername postgresql --hostname pdc --dbname pmdb --db_user luozhihui \
    --db_passwd yfishLab2113 --schema xiandao --ref_ind_seq_id 1 -l ycondor -j ycondor
"""
import sys, os, math
sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))
sys.path.insert(0, os.path.join(os.path.expanduser('~/src')))
import copy
from pymodule import ProcessOptions, getListOutOfStr, PassingData, utils
from pymodule.reducer.AbstractReducer import AbstractReducer
from pegaflow.DAX3 import Executable, PFN, File
from sqlalchemy.sql import text
from . AbstractAccuWorkflow import AbstractAccuWorkflow
from pymodule.db import SunsetDB

parentClass = AbstractAccuWorkflow
class DownsampleWorkflow(parentClass):
    __doc__ = __doc__
    option_default_dict = copy.deepcopy(parentClass.option_default_dict)
    option_default_dict.update(parentClass.db_option_dict.copy())
    option_default_dict.update({
        ("thisModulePath", 1,): ["%s/src/Sunset", '', 1, 'path of the module that owns this program. \
    					used to add executables from this module.'], \
        ("AccurityPath", 1,):["%s/src/Sunset/src_o/main.py", '', 1, 'path of the AccurityPath'], \
        })


    def __init__(self,  **keywords):
        self.pathToInsertHomePathList.extend(["thisModulePath", "AccurityPath"])
        parentClass.__init__(self, **keywords)


    def addDownsamplejob(self, workflow=None, data_dir=None, idDict=None, DownSamplePrefix=None, \
                         downSampleJava=None, downSampleJar=None, transferOutput=False):
        AccurityFolder = "AccurityResult"
        AccurityFolderJob = self.addMkDirJob(executable=workflow.mkdirWrap, outputDir=AccurityFolder)

        sys.stderr.write("Adding downsample jobs for %s individual sequences ..." % (len(idDict)))
        SampleFolder = "%swithSeed1.0" % (DownSamplePrefix)
        SampleFolderJob = self.addMkDirJob(executable=workflow.mkdirWrap, outputDir=SampleFolder)

        alignNormal = self.db_main.queryTable(SunsetDB.IndividualAlignment).get(idDict['normalFile'])
        alignNormalFilePath = os.path.join(data_dir, alignNormal.path)
        inputNormalBamFile = self.registerOneInputFile(inputFname=alignNormalFilePath)
        coverageNormal = int(alignNormal.mean_depth)
        #alignNormalIndiv = self.db_main.queryTable(SunsetDB.IndividualSequence).get(alignNormal.ind_seq_id)
        #coverageNormal = int(alignNormalIndiv.coverage)

        alignTumor = self.db_main.queryTable(SunsetDB.IndividualAlignment).get(idDict['tumorFile'])
        alignTumorFilePath = os.path.join(data_dir, alignTumor.path)
        inputTumorBamFile = self.registerOneInputFile(inputFname=alignTumorFilePath)
        coverageTumor = int(alignTumor.mean_depth)
        #alignTumorIndiv = self.db_main.queryTable(SunsetDB.IndividualSequence).get(alignTumor.ind_seq_id)
        #coverageTumor = int(alignTumorIndiv.coverage)

        job_max_memory = "5000"
        walltime = '600'
        for i in range(1,10):
            jobLs = []
            pair_bam_file_list = []
            probNormal = float(i) / float(coverageNormal)
            probTumor = float(10 - i) / float(coverageTumor)
            outputNormalFile = File(os.path.join(SampleFolder, str(probNormal) + "_normal_downsample.bam"))
            outputTumorFile = File(os.path.join(SampleFolder, str(probTumor) + "_tumor_downsample.bam"))
            mergeJobAndOutputLs = []
            normal_down_sample_job = self.addGenericJavaJob(executable=downSampleJava, jarFile=downSampleJar, \
                                                       inputFile=inputNormalBamFile, inputArgumentOption="INPUT=", \
                                                       inputFileList=None,
                                                       argumentForEachFileInInputFileList=None, \
                                                       outputFile=outputNormalFile, outputArgumentOption="OUTPUT=", \
                                                       parentJobLs=[SampleFolderJob], transferOutput=False,
                                                       job_max_memory=job_max_memory, \
                                                       frontArgumentList=['DownsampleSam'], extraArguments=None,
                                                       extraArgumentList=['PROBABILITY=' + str(probNormal), \
                                                                          'RANDOM_SEED=','1' ,\
                                                                          'STRATEGY=','ConstantMemory', \
                                                                          'VALIDATION_STRINGENCY=','LENIENT'
                                                                          ],
                                                       extraOutputLs=None, \
                                                       extraDependentInputLs=None, no_of_cpus=None, walltime=walltime,
                                                       sshDBTunnel=None)
            mergeJobAndOutputLs.append(PassingData(jobLs=[normal_down_sample_job], file=outputNormalFile))

            tumor_down_sample_job = self.addGenericJavaJob(executable=downSampleJava, jarFile=downSampleJar, \
                                                            inputFile=inputTumorBamFile, inputArgumentOption="INPUT=", \
                                                            inputFileList=None,
                                                            argumentForEachFileInInputFileList=None, \
                                                            outputFile=outputTumorFile, outputArgumentOption="OUTPUT=", \
                                                            parentJobLs=[SampleFolderJob],
                                                            transferOutput=False,
                                                            job_max_memory=job_max_memory, \
                                                            frontArgumentList=['DownsampleSam'], extraArguments=None,
                                                            extraArgumentList=['PROBABILITY=' , str(probTumor), \
                                                                               'RANDOM_SEED=','1', \
                                                                               'STRATEGY=', 'ConstantMemory', \
                                                                               'VALIDATION_STRINGENCY=' ,'LENIENT'
                                                                               ],
                                                            extraOutputLs=None, \
                                                            extraDependentInputLs=None, no_of_cpus=None,
                                                            walltime=walltime,
                                                            sshDBTunnel=None)
            mergeJobAndOutputLs.append(PassingData(jobLs=[tumor_down_sample_job], file=outputTumorFile))

            puritySampleFolder = "puritySample"
            SampleFolderJob = self.addMkDirJob(executable=workflow.mkdirWrap, outputDir=puritySampleFolder)
            purity = str((10-i) * 0.1)
            purityDir = "purity" + str(purity)
            purityFolderJob = self.addMkDirJob(executable=workflow.mkdirWrap, outputDir=os.path.join(puritySampleFolder,purityDir))
            mergedBamFile = File(os.path.join(puritySampleFolder,purityDir, "purity_"+ purity + ".bam"))
            baseCoverage = 4 * 3000000000  # baseline
            minMergeAlignmentWalltime = 240  # in minutes, 4 hours, when coverage is defaultCoverage
            maxMergeAlignmentWalltime = 2980  # in minutes, 2 days
            minMergeAlignmentMaxMemory = 8000  # in MB, when coverage is defaultCoverage
            maxMergeAlignmentMaxMemory = 21000  # in MB

            mergeAlignmentWalltime = self.scaleJobWalltimeOrMemoryBasedOnInput(
                realInputVolume=max(i, 10-i) * 3000000000, \
                baseInputVolume=baseCoverage, baseJobPropertyValue=minMergeAlignmentWalltime, \
                minJobPropertyValue=minMergeAlignmentWalltime, maxJobPropertyValue=maxMergeAlignmentWalltime).value
            mergeAlignmentMaxMemory = self.scaleJobWalltimeOrMemoryBasedOnInput(
                realInputVolume=max(i, 10-i) * 3000000000, \
                baseInputVolume=baseCoverage, baseJobPropertyValue=minMergeAlignmentMaxMemory, \
                minJobPropertyValue=minMergeAlignmentMaxMemory, maxJobPropertyValue=maxMergeAlignmentMaxMemory).value

            MergeJob, bamIndexJob = self.addAlignmentMergeJob(workflow=workflow, \
                                                AlignmentJobAndOutputLs=mergeJobAndOutputLs, \
                                                outputBamFile=mergedBamFile, \
                                                samtools=workflow.samtools, java=workflow.java, \
                                                MergeSamFilesJava=workflow.MergeSamFilesJava,
                                                MergeSamFilesJar=workflow.MergeSamFilesJar, \
                                                BuildBamIndexFilesJava=workflow.IndexMergedBamIndexJava, \
                                                BuildBamIndexJar=workflow.BuildBamIndexJar, \
                                                mv=workflow.mv, namespace=workflow.namespace, \
                                                version=workflow.version, \
                                                transferOutput=transferOutput, \
                                                job_max_memory=mergeAlignmentMaxMemory, \
                                                walltime=mergeAlignmentWalltime, \
                                                parentJobLs=[SampleFolderJob, purityFolderJob])
            normal_part_refer = self.registerOneInputFile(
                inputFname="/y/Sunset/workflow/real_data/downsample/normal_0.2.bam",folderName=os.path.join(puritySampleFolder,purityDir))
            normal_bam_bai = self.registerOneInputFile(
                inputFname="/y/Sunset/workflow/real_data/downsample/normal_0.2.bam.bai",folderName=os.path.join(puritySampleFolder,purityDir))
            pair_bam_file_list.append([mergedBamFile, normal_part_refer])
            AccurityJob = self.doAllAccurityAlignmentJob(workflow=workflow, data_dir=None,  normal_bam_bai=normal_bam_bai,\
                                                         pair_bam_file_list=pair_bam_file_list,\
                                                         outputDirPrefix=None, parentJobLs=[MergeJob, bamIndexJob],\
                                                         AccurityFolder=AccurityFolder, AccurityFolderJob=AccurityFolderJob)

        return jobLs


    def doAllAccurityAlignmentJob(self, workflow=None, data_dir=None, normal_bam_bai=None, pair_bam_file_list = None,\
                               outputDirPrefix=None, parentJobLs=None, AccurityFolder=None, AccurityFolderJob=None):
        sys.stderr.write("Adding Accurity jobs for %s pair individual sequences ..." % (len(pair_bam_file_list)))

        jobLs = []
        for pair_bam in pair_bam_file_list:
            tumor_bam = pair_bam[0]
            tumor_bam_bai = parentJobLs[1].baiFile
            normal_bam = pair_bam[1]

            if tumor_bam is None or normal_bam is None:
                sys.stderr.write("the pair sample bam file is note exist!!!")
                exit(2)
            #tumor_bam_path = os.path.join(data_dir, tumor_bam)
            #tumor_bai_path = tumor_bam_path + ".bai"
            #normal_bam_path = os.path.join(data_dir, normal_bam)
            #normal_bai_path = normal_bam_path + ".bai"
            Accurity_configure_path = os.path.dirname(self.AccurityPath) + "/configure"

            outputList = []
            sample_id = os.path.basename(tumor_bam.name).strip(".bam")
            sample_folder = AccurityFolder + "/" + sample_id

            sample_folder_Job = self.addMkDirJob(executable=workflow.mkdirWrap,outputDir=sample_folder, parentJobLs=parentJobLs.append( AccurityFolderJob))
            outputList.append(File(sample_folder + "/infer.out.tsv"))
            outputList.append(File(sample_folder + "/infer.out.details.tsv"))
            outputList.append(File(sample_folder + "/auto.tsv"))
            outputList.append(File(sample_folder + "/cnv.plot.pdf"))
            outputList.append(File(sample_folder + "/cnv.output.tsv"))
            outputList.append(File(sample_folder + "/rc_ratio_window_count_smoothed.tsv"))
            outputList.append(File(sample_folder + "/rc_ratio_no_of_windows_by_chr.tsv"))
            outputList.append(File(sample_folder + "/cnv.intervel.tsv"))
            outputList.append(File(sample_folder + "/major_allele_fraction_exp_vs_obs.tsv"))
            outputList.append(File(sample_folder + "/peak_bounds.tsv"))
            outputList.append(File(sample_folder + "/rc_logLikelihood.log.tsv"))
            outputList.append(File(sample_folder + "/rc_ratios_of_peaks_based_on_period_from_autocor.tsv"))
            outputList.append(File(sample_folder + "/runTime.log.txt"))

            #tumor_bam_file = self.registerOneInputFile(inputFname=tumor_bam_path)
            #tumor_bai_file = self.registerOneInputFile(inputFname=tumor_bai_path)
            #normal_bam_file = self.registerOneInputFile(inputFname=normal_bam_path)
            #normal_bai_file = self.registerOneInputFile(inputFname=normal_bai_path)
            configure_file = self.registerOneInputFile(inputFname=Accurity_configure_path)
            argumentList = ["-c", configure_file, "-t", tumor_bam, "-n", normal_bam, "-o", sample_folder, "-d", "1", "-l", "4"]
            inputFileList = [tumor_bam, tumor_bam_bai, normal_bam, normal_bam_bai, configure_file]

            job = self.addPurityJobToWorkflow(workflow=self.workflow, executable=self.AccurityExecutableFile,\
                                              argumentList=argumentList, \
                                              inputFileList=inputFileList, outputFileList=outputList, \
                                              parentJobLs=[sample_folder_Job], \
                                              job_max_memory=10000, no_of_cpus=8, walltime=400, sshDBTunnel=0)
            jobLs.append(job)
        return jobLs

    def addPurityJobToWorkflow(self, workflow=None, executable=None, argumentList=None, inputFileList=None,\
                               outputFileList=None, parentJobLs=None, job_max_memory=10000, no_of_cpus=1, \
                               walltime=400, sshDBTunnel=0):
        job = self.addGenericJob(workflow=workflow, executable=executable, inputFile=None, inputArgumentOption=None, \
                          outputFile=None, outputArgumentOption=None, inputFileList=None,
                          argumentForEachFileInInputFileList=None, \
                          parentJob=None, parentJobLs=parentJobLs, extraDependentInputLs=inputFileList,\
                          extraOutputLs=outputFileList, \
                          frontArgumentList=argumentList, extraArguments=None, extraArgumentList=None, \
                          transferOutput=True, sshDBTunnel=sshDBTunnel, \
                          key2ObjectForJob=None, objectWithDBArguments=None, objectWithDBGenomeArguments=None, \
                          no_of_cpus=no_of_cpus, job_max_memory=job_max_memory, walltime=walltime, \
                          max_walltime=None)
        return job


    def run(self):
        IDsample = {"normalFile":380, "tumorFile":381}

        pdata = self.setup_run()
        workflow = pdata.workflow

        self.AccurityExecutableFile = self.addOneExecutableFromPathAndAssignProperClusterSize( \
            path=self.AccurityPath, name='Accuritymainexecutable', clusterSizeMultipler=0)

        self.downSampleJava = self.addOneExecutableFromPathAndAssignProperClusterSize( \
            path=self.javaPath, name='sampleDownexecutable', clusterSizeMultipler=0)

        self.addDownsamplejob(workflow=workflow, data_dir=self.data_dir, idDict=IDsample, DownSamplePrefix=None, \
                         downSampleJava=self.downSampleJava, downSampleJar=self.PicardJar, \
                         transferOutput=False)
        #AccurityJobs = self.doAllAccurityAlignmentJob(workflow=workflow, data_dir=self.data_dir, \
        #                              pair_bam_file_list=entry, outputDirPrefix="")
        #self.addReducerJobtoAccurity(lastJob=AccurityJobs, dir="resultReducer", outputFnameLastStep="resultList.txt")
        self.end_run()


if __name__ == '__main__':
    main_class = DownsampleWorkflow
    po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
    instance = main_class(**po.long_option2value)
    instance.run()
