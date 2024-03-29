#!/usr/bin/env python3
"""
2013.07.09 --ref_ind_seq_id is not required to specify 
    as it infers reference sequence ID from sample IDs in VCF
    (alignment read-group -> reference sequence).
#2012.5.9
    the usual -c (commit) is not here. All DB jobs are run with commit=True.
2012.8.3 if such a workflow with clustering on (several AddVCFFolder2DB jobs crammed into one) fails halfway,
    you can safely re-run it. Already-imported files would be checked and not be imported again (MD5SUM).

Examples:
    #2012.5.11 
    %s -I FilterVCF_VRC_SK_Nevis_FilteredSeq_top1000Contigs/trioCaller_vcftoolsFilter/ 
        -o dags/2DB/AddVCF2DB_FilterVCF_VRC_SK_Nevis_2012.5.6_trioCaller.xml 
        -s ... -u yh -l hcondor -j hcondor  -z localhost 
        -e /u/home/eeskin/polyacti/ -t NetworkData/vervet/db/ -D NetworkData/vervet/db/ 
    
    # 2012.5.10 run on hoffman2 condor, turn on checkEmptyVCFByReading (-E),
    #  require db connection (--needSSHDBTunnel), no clustering (-C1)
    %s -I FilterVCF_VRC_SK_Nevis_FilteredSeq_top1000Contigs/trioCaller_vcftoolsFilter/
        -o dags/AddVCF2DB_FilterVCF_VRC_SK_Nevis_FilteredSeq_trioCaller.xml
        -s ... -E --needSSHDBTunnel -C1
        -l hcondor -j hcondor  -u yh -z localhost 
        -e /u/home/eeskin/polyacti/ 
        -D NetworkData/vervet/db/ -t NetworkData/vervet/db/
    
    # 2012.7.16 add a folder of VCF files to DB without checking zero-loci VCF
    %s -I FilterVCF_VRC_SK_Nevis_FilteredSeq_top1000Contigs.MAC10.MAF.05/trioCaller_vcftoolsFilter/ 
        -o dags/AddVCF2DB_FilterVCF_VRC_SK_Nevis_top1000Contigs.MAC10.MAF.05_trioCaller.xml
        -s ... -l condorpool -j condorpool
        -u yh -z uclaOffice -C1

"""
import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0], sys.argv[0])

import csv, copy
from pegaflow.api import File
from palos import ProcessOptions, getListOutOfStr, PassingData
from palos.db import SunsetDB
from palos.ngs.GenericVCFWorkflow import GenericVCFWorkflow
from palos.ngs.AbstractNGSWorkflow import AbstractNGSWorkflow

class AddVCFFolder2DBWorkflow(GenericVCFWorkflow):
    __doc__ = __doc__
    option_default_dict = copy.deepcopy(AbstractNGSWorkflow.option_default_dict)
    option_default_dict.update({
        ('inputDir', 0, ): ['', 'I', 1, 'input folder that contains vcf or vcf.gz files', ],\
        ('maxContigID', 0, int): [None, 'x', 1, 'if contig ID is beyond this number, '
            'it will not be included. If None or 0, no restriction.', ],\
        ('genotypeMethodShortName', 1, ):[None, 's', 1, 'column short_name of GenotypeMethod table, '
            'will be created if not present in db.'],\
        ('run_type', 1, int): [1, 'y', 1, 'which run_type to run. '],\
        })
    option_default_dict[("thisModulePath", 1, )][0] = '%s/src/Sunset'

    def __init__(self,  **keywords):
        """
        """
        AbstractNGSWorkflow.__init__(self, **keywords)
        self.inputDir = os.path.abspath(self.inputDir)
    
    def registerCustomExecutables(self):
        """
        2011-11-28
        """
        self.registerOneExecutable(path=os.path.join(self.pymodulePath, \
            "db/import/AddGenotypeMethod2DB.py"),
            name='AddGenotypeMethod2DB', clusterSizeMultiplier=1)
        self.registerOneExecutable(path=os.path.join(self.pymodulePath, \
            "db/UpdateGenotypeMethodNoOfLoci.py"),
            name='UpdateGenotypeMethodNoOfLoci', clusterSizeMultiplier=1)

    def addAddGenotypeMethod2DBJob(self, executable=None, inputFile=None, genotypeMethodShortName=None,\
        logFile=None, data_dir=None, commit=False, parentJobLs=None, \
        extraDependentInputLs=None, transferOutput=False, \
        extraArguments=None, job_max_memory=2000, **keywords):
        """
        2012.6.27
        """
        extraArgumentList = ['--genotypeMethodShortName', genotypeMethodShortName]
        if logFile:
            extraArgumentList.extend(["--logFilename", logFile])
        if data_dir:
            extraArgumentList.extend(['--data_dir', data_dir])
        if commit:
            extraArgumentList.append('--commit')
        if extraArguments:
            extraArgumentList.append(extraArguments)
        
        job= self.addGenericJob(
            executable=executable, inputFile=inputFile, outputFile=None, \
            parentJobLs=parentJobLs, extraDependentInputLs=extraDependentInputLs, \
            extraOutputLs=[logFile],\
            transferOutput=transferOutput, \
            extraArgumentList=extraArgumentList, job_max_memory=job_max_memory, **keywords
            )
        self.addDBArgumentsToOneJob(job=job, objectWithDBArguments=self)
        return job
    
    def addUpdateGenotypeMethodNoOfLociJob(
        self, executable=None, genotypeMethodShortName=None, genotypeMethodID=None,\
        logFile=None, data_dir=None, commit=False, parentJobLs=None, 
        extraDependentInputLs=None, transferOutput=False, \
        extraArguments=None, job_max_memory=2000, **keywords):
        """
        2012.6.27
        """
        extraArgumentList = []
        if logFile:
            extraArgumentList.extend(["--logFilename", logFile])
        if data_dir:
            extraArgumentList.extend(['--data_dir', data_dir])
        if commit:
            extraArgumentList.append('--commit')
        if genotypeMethodShortName:
            extraArgumentList.extend(['--genotypeMethodShortName', genotypeMethodShortName, ])
        if genotypeMethodID:
            extraArgumentList.extend(['--genotypeMethodID', genotypeMethodID, ])
        if extraArguments:
            extraArgumentList.append(extraArguments)
        
        job= self.addGenericJob(executable=executable, inputFile=None, outputFile=None, \
            parentJobLs=parentJobLs, extraDependentInputLs=extraDependentInputLs, \
            extraOutputLs=[logFile],\
            transferOutput=transferOutput, \
            extraArgumentList=extraArgumentList, job_max_memory=job_max_memory, **keywords)
        self.addDBArgumentsToOneJob(job=job, objectWithDBArguments=self)
        return job
    
    def addJobs(self, inputData=None, db_main=None, genotypeMethodShortName=None, commit=None,\
            data_dir=None, checkEmptyVCFByReading=False, transferOutput=True,\
            maxContigID=None, outputDirPrefix="", needSSHDBTunnel=False):
        """
        2012.5.9
        """
        sys.stderr.write("Adding VCF2DB jobs for %s vcf files ... "%(len(inputData.jobDataLs)))
        
        
        topOutputDir = "%sVCF2DB"%(outputDirPrefix)
        topOutputDirJob = self.addMkDirJob(outputDir=topOutputDir)
        
        firstVCFFile = inputData.jobDataLs[0].vcfFile
        logFile = File(os.path.join(topOutputDir, 'AddGenotypeMethod2DB.log'))
        addGM2DBJob = self.addAddGenotypeMethod2DBJob(
            executable=self.AddGenotypeMethod2DB, inputFile=firstVCFFile, \
            genotypeMethodShortName=genotypeMethodShortName,\
            logFile=logFile, data_dir=data_dir, commit=commit, parentJobLs=None, 
            extraDependentInputLs=None, transferOutput=True, \
            extraArguments=None, job_max_memory=10, sshDBTunnel=needSSHDBTunnel)
        updateGMlogFile = File(os.path.join(topOutputDir, 'updateGM.log'))
        updateGMNoOfLociJob = self.addUpdateGenotypeMethodNoOfLociJob(
            executable=self.UpdateGenotypeMethodNoOfLoci, \
            genotypeMethodShortName=genotypeMethodShortName,\
            logFile=updateGMlogFile, data_dir=data_dir, commit=commit,
            parentJobLs=[topOutputDirJob], \
            extraDependentInputLs=[], transferOutput=True, \
            extraArguments=None, job_max_memory=20, sshDBTunnel=needSSHDBTunnel)
        
        returnData = PassingData()
        returnData.jobDataLs = []
        for jobData in inputData.jobDataLs:
            inputF = jobData.vcfFile
            if maxContigID:
                contig_id = self.getContigIDFromFname(inputF.name)
                try:
                    contig_id = int(contig_id)
                    if contig_id>maxContigID:	#skip the small contigs
                        continue
                except:
                    sys.stderr.write('Except type: %s\n'%repr(sys.exc_info()))
                    import traceback
                    traceback.print_exc()
            logFile = File(os.path.join(topOutputDir, 'AddVCFFile2DB_%s.log'%(self.getChrFromFname(inputF.name))))
            addVCFJob = self.addAddVCFFile2DBJob(
                executable=self.AddVCFFile2DB, inputFile=inputF, 
                genotypeMethodShortName=genotypeMethodShortName,\
                logFile=logFile, format="VCF", data_dir=data_dir, 
                checkEmptyVCFByReading=checkEmptyVCFByReading, commit=commit, \
                parentJobLs=[addGM2DBJob]+jobData.jobLs, extraDependentInputLs=[], transferOutput=True, \
                extraArguments=None, job_max_memory=1000, sshDBTunnel=needSSHDBTunnel)
            self.add_dependency(updateGMNoOfLociJob, parents=[addVCFJob])
        sys.stderr.write("%s jobs.\n"%(self.no_of_jobs))
        #include the tfam (outputList[1]) into the fileLs
        returnData.jobDataLs.append(PassingData(jobLs=[updateGMNoOfLociJob],
            file=updateGMlogFile, \
            fileLs=[updateGMlogFile]))
        return returnData
    
    def run(self):
        """
        2011-9-28
        """
        
        if self.debug:
            import pdb
            pdb.set_trace()
        
        db_main = self.db_main
        
        # Create a abstract dag
        workflow = self.initiateWorkflow()
        
        self.registerJars(workflow)
        self.registerExecutables(workflow)
        self.registerCustomExecutables(workflow)
        
        inputData = self.registerAllInputFiles(
            inputDir=self.inputDir, input_site_handler=self.input_site_handler, \
            checkEmptyVCFByReading=self.checkEmptyVCFByReading,\
            pegasusFolderName=self.pegasusFolderName)
        if len(inputData.jobDataLs)<=0:
            sys.stderr.write("Warning: Number of VCF files in this folder, %s, <=0.\n"%self.inputDir)
            sys.exit(0)
        
        if self.run_type==1:
            self.addJobs(
                inputData=inputData, db_main=db_main, genotypeMethodShortName=self.genotypeMethodShortName, \
                commit=True,\
                data_dir=self.data_dir, checkEmptyVCFByReading=self.checkEmptyVCFByReading, transferOutput=True,\
                maxContigID=self.maxContigID, outputDirPrefix="", needSSHDBTunnel=self.needSSHDBTunnel)
        else:
            sys.stderr.write("run_type %s not supported.\n"%(self.run_type))
            sys.exit(0)
        # Write the DAX to stdout
        outf = open(self.outputFname, 'w')
        self.writeXML(outf)


if __name__ == '__main__':
    main_class = AddVCFFolder2DBWorkflow
    po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
    instance = main_class(**po.long_option2value)
    instance.run()