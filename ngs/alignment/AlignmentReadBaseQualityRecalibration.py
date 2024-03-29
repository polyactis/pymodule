#!/usr/bin/env python3
"""
Description:
    #2013.04.09  recalibration needs ~500K reads to get accurate estimate.
    #  So adjust the --intervalSize according to coverage.  
    Input: a VCF folder, list of alignment IDs (or use site_id, ref_ind_seq_id to filter )
    It could run local_realignment with or without indel VCF and BQSR with or without SNPVCF. 
    Without a SNPVCF from particular chromosome, it uses a random chromosome. 
    If the interval-span is too short , it use BQSR result from long-intervals.
    
    add alignment file into db (input: parent alignment ID, alignment file path)

Examples:
    # 2012.9.21 run base quality recalibration on VRC alignments (-S 447),
    #  individual_sequence_id from 639-642 (--ind_seq_id_ls ...)
    # filtered sequences (-Q 1), alignment method 2 (-G 2)
    # --contigMaxRankBySize 1000 (top 1000 contigs)
    #  --intervalSize 10000000 (10 million bp for each interval)
    #  --intervalOverlapSize 30000 (30kb overlap between intervals),
    %s --inputDir ~/NetworkData/vervet/db/genotype_file/method_17/
        --ind_seq_id_ls 639-642
        -S 447 -u yh -z localhost --sequence_filtered 1 --alignment_method_id 2
        -a 524 -o dags/BaseQualityRecalibration_VRC447_vsMethod17.xml
        -l hcondor -j hcondor -z localhost -u yh --contigMaxRankBySize 1000 
        --intervalSize 10000000 --intervalOverlapSize 30000
        --indelVCFFolder ...  -e /u/home/eeskin/polyacti
        --local_data_dir NetworkData/vervet/db/
        --data_dir /u/home/eeskin/polyacti/NetworkData/vervet/db/
        --needSSHDBTunnel -J ~/bin/jdk/bin/java --new_mask_genotype_method_id 17
         --commit --skipDoneAlignment

    # 2012.9.18
    %s  --inputDir ~/NetworkData/vervet/db/genotype_file/method_41
        --ind_seq_id_ls 633,634,635,636,637,638 
        --ref_ind_seq_id 524
        -o dags/BaseQualityRecalibration_ISQ633_638_vsMethod41.xml -l hcondor
        -j hcondor -z localhost -u yh --intervalSize 10000000
        --intervalOverlapSize 30000
        -e /u/home/eeskin/polyacti
        --indelVCFFolder ...
        --local_data_dir NetworkData/vervet/db/ --data_dir NetworkData/vervet/db/
        --cluster_size 5 --needSSHDBTunnel -J ~/bin/jdk/bin/java
        --new_mask_genotype_method_id 41
        --commit --skipDoneAlignment
    
    # 2013.3.19 use sequence coverage to filter alignments
    %s  --inputDir ~/NetworkData/vervet/db/genotype_file/method_41
        --sequence_min_coverage 0 --sequence_max_coverage 2
        --ind_seq_id_ls 632-3230
        --ref_ind_seq_id 3280
        -o dags/BaseQualityRecalibration_ISQ632_3230_coverage0_2_vsMethod41.xml
        -l hcondor -j hcondor -z localhost -u yh
        --intervalSize 10000000 --intervalOverlapSize 30000
        -e /u/home/eeskin/polyacti --contigMaxRankBySize 250
        --local_data_dir NetworkData/vervet/db/ --data_dir NetworkData/vervet/db/
        --cluster_size 5 --needSSHDBTunnel -J ~/bin/jdk/bin/java
        --new_mask_genotype_method_id 41
        --indelVCFFolder ~/NetworkData/vervet/db/genotype_file/method_88
        --commit --skipDoneAlignment
        # --ref_genome_version 2 #(optional, as by default,
        #  it gets the outdated_index=0 reference chromosomes from GenomeDB)
        # --ref_genome_outdated_index 1 #to get old reference.
        #  incompatible here as alignment is based on 3280, new ref.
        # --needFastaDictJob --needFastaIndexJob
    
"""
import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0], sys.argv[0])

from pegaflow.api import File
import logging
from palos import ProcessOptions, getListOutOfStr, PassingData, ngs, \
    figureOutDelimiter, getColName2IndexFromHeader, utils
from palos.ngs.io.VCFFile import VCFFile
from ShortRead2Alignment import ShortRead2Alignment

ParentClass = ShortRead2Alignment
class AlignmentReadBaseQualityRecalibration(ParentClass):
    __doc__ = __doc__
    option_default_dict = ParentClass.option_default_dict.copy()
    option_default_dict.update(ParentClass.commonAlignmentWorkflowOptionDict.copy())
    option_default_dict.update(ParentClass.partitionWorkflowOptionDict.copy())
    option_default_dict.update({
        ('new_mask_genotype_method_id', 0, int):[None, '', 1, 
            'which genotype method is used to mask out polymorphic sites '
            'for newly-recalibrated of alignment. '
            'This is different from --mask_genotype_method_id, which is used to filter '
            'input alignments and should be set to None (leave it at default) '
            'for this purpose.'],\
        ('indelVCFFolder', 0, ): [None, '', 1, 
            'folder that contains in-del vcf or vcf.gz files that will '
            'be used for GATK2 indel-realigner, '
            'required if argument local_realigned is non-zero', ],\
        })
    option_default_dict[('intervalSize', 1, int)][0] = 10000000
    option_default_dict[('local_realigned', 0, int)][0] = 0
    option_default_dict[('completedAlignment', 0, int)][0]=1
    """
    option_default_dict.update({
        ('ind_aln_id_ls', 0, ): ['', 'I', 1, 
            'a comma/dash-separated list of IndividualAlignment.id.'
            ' This overrides ind_seq_id_ls.', ],\
        ('inputDir', 0, ): [None, 'i', 1, 
            "folder containing vcf or vcf.gz files that will be used by"
            " CountCovariates to exclude variant sites", ],\
        ("sequence_filtered", 0, int): [None, 'Q', 1, 
            'To filter alignments. None: whatever; 0: unfiltered sequences, '
            '1: filtered sequences'],\
        ("alignment_method_id", 0, int): [None, 'G', 1, 
            'To filter alignments. None: whatever; integer: AlignmentMethod.id'],\
        ('intervalOverlapSize', 1, int): [3000, '', 1, 
            'extension of an interval on each side. overlap size is actually 2X of this number.'
            'Its value should be bigger than maximum read length, to cover '
            'reads that are partially aligned within one interval.', ],\
        ('intervalSize', 1, int): [5000000, '', 1, 
            'size for adjacent intervals from one contig/chromosome', ],\
        ("site_id_ls", 0, ): ["", 'S', 1,
            'comma/dash-separated list of site IDs. individuals must come from these sites.'],\
        ('run_type', 1, int): [1, 'y', 1, '', ],\
        })
    """
    
    def __init__(self,  **keywords):
        """
        """
        ParentClass.__init__(self, **keywords)
        self.chr2IndelVCFJobData = None
        #2013.04.04 mark this variable. setup in setup()
        self.candidateCountCovariatesJob = None
        #2013.04.09 this BQSR count-variates job encompasses one of the top big intervals.
        # replacing equivalent jobs for small intervals
        #  (not accurate if intervals are too small)
        #AlignmentToCallPipeline.__init__(self, **keywords)
        #self.inputDir = os.path.abspath(self.inputDir)
    
    def mapEachAlignment(self, passingData=None, transferOutput=True,
        **keywords):
        """
        2012.9.22
            similar to reduceBeforeEachAlignmentData()
             but for mapping programs that run on one alignment each.
        """
        returnData = PassingData(no_of_jobs = 0)
        returnData.jobDataLs = []
        return returnData
    
    def mapEachChromosome(self, alignmentData=None, chromosome=None,\
        VCFJobData=None, passingData=None, reduceBeforeEachAlignmentData=None,
        transferOutput=True, **keywords):
        """
        2012.9.17
        """
        returnData = PassingData(no_of_jobs = 0)
        returnData.jobDataLs = []
        
        topOutputDirJob = passingData.topOutputDirJob
        
        alignment = alignmentData.alignment
        parentJobLs = alignmentData.jobLs
        bamF = alignmentData.bamF
        baiF = alignmentData.baiF
        bamFnamePrefix = passingData.bamFnamePrefix
        
        
        """
        #2012.9.21 perhaps a downsampling job
        outputFname = os.path.join(topOutputDirJob.output, \
            '%s_%s.bam'%(bamFnamePrefix, overlapFileBasenameSignature))
        outputFile = File(outputFname)
        selectAlignmentJob, bamIndexJob1 = self.addSelectAlignmentJob(
            executable=self.samtools, inputFile=bamF, \
            outputFile=outputFile, region=overlapInterval,
            parentJobLs=[topOutputDirJob] + parentJobLs, \
            extraDependentInputLs=[baiF], transferOutput=False, \
            extraArguments=None, job_max_memory=2000, needBAMIndexJob=True)
        """
        
        """
        #2012.9.21 count covariates job is moved to map()
        recalFile = File(os.path.join(topOutputDirJob.output,
            '%s_%s.recal_data.csv'%(bamFnamePrefix, chromosome)))
        countCovariatesJob = self.addGATKBaseRecalibratorJob(
            GenomeAnalysisTKJar=self.GenomeAnalysisTK2Jar, inputFile=bamF, \
            VCFFile=VCFFile, interval=chromosome, outputFile=recalFile, \
            refFastaFList=passingData.refFastaFList,
            parentJobLs=[topOutputDirJob]+parentJobLs, 
            extraDependentInputLs=[baiF, VCFFile.tbi_F], \
            transferOutput=False, \
            extraArguments=None, job_max_memory=4000)

        self.no_of_jobs += 1
        returnData.countCovariatesJob = countCovariatesJob
        returnData.jobDataLs.append(PassingData(jobLs=[countCovariatesJob],
            file=countCovariatesJob.recalFile, \
            fileLs=[countCovariatesJob.recalFile]))
        """
        
        return returnData
    
    def reduceAfterEachAlignment(self, passingData=None, transferOutput=False,
        data_dir=None, **keywords):
        """
        """
        returnData = PassingData(no_of_jobs = 0)
        returnData.jobDataLs = []
        alignmentJobAndOutputLs = getattr(passingData, 'alignmentJobAndOutputLs', [])
        bamFnamePrefix = passingData.bamFnamePrefix
        topOutputDirJob = passingData.topOutputDirJob
        individual_alignment = passingData.individual_alignment
        reduceOutputDirJob = passingData.reduceOutputDirJob
        
        if len(alignmentJobAndOutputLs)>0:
            #2012.3.29	merge alignment output only when there is something to merge!
            #2013.04.09 create a new child alignment local_realigned =1, etc.
            new_individual_alignment = self.db.copyParentIndividualAlignment(
                parent_individual_alignment_id=individual_alignment.id,\
                mask_genotype_method_id=self.new_mask_genotype_method_id,\
                data_dir=self.data_dir, local_realigned=1)
            
            baseCoverage = 4	#baseline
            actualCoverage = getattr(individual_alignment.individual_sequence,
                'coverage', baseCoverage)
            minMergeAlignmentWalltime = 240
            #in minutes, 4 hours, when coverage is defaultCoverage
            maxMergeAlignmentWalltime = 2880	#in minutes, 2 days
            minMergeAlignmentMaxMemory = 7000
            #in MB, when coverage is defaultCoverage
            maxMergeAlignmentMaxMemory = 12000	#in MB
            
            mergeAlignmentWalltime = self.scaleJobWalltimeOrMemoryBasedOnInput(
                realInputVolume=actualCoverage,
                baseInputVolume=baseCoverage,
                baseJobPropertyValue=minMergeAlignmentWalltime*2,
                minJobPropertyValue=minMergeAlignmentWalltime,
                maxJobPropertyValue=maxMergeAlignmentWalltime).value
            mergeAlignmentMaxMemory = self.scaleJobWalltimeOrMemoryBasedOnInput(
                realInputVolume=actualCoverage,
                baseInputVolume=baseCoverage,
                baseJobPropertyValue=minMergeAlignmentMaxMemory,
                minJobPropertyValue=minMergeAlignmentMaxMemory, 
                maxJobPropertyValue=maxMergeAlignmentMaxMemory).value
            
            # replace read_group with the new one to each alignment job
            newAlignmentJobAndOutputLs = []
            for alignmentJobAndOutput in alignmentJobAndOutputLs:
                # add a AddReadGroup job
                alignmentJob, indexAlignmentJob = alignmentJobAndOutput.jobLs[:2]
                fileBasenamePrefix = os.path.splitext(alignmentJob.output.name)[0]
                outputRGBAM = File("%s.isq_RG.bam"%(fileBasenamePrefix))
                # needBAMIndexJob=False because addAlignmentMergeJob()
                # does not need .bai.
                addRGJob = self.addReadGroupJob(
                    individual_alignment=new_individual_alignment,
                    inputBamFile=alignmentJob.output,
                    outputBamFile=outputRGBAM,
                    needBAMIndexJob=False,
                    parentJobLs=[alignmentJob, indexAlignmentJob],
                    extraDependentInputLs=alignmentJob.outputLs[1:],
                    job_max_memory = 2500, transferOutput=False,
                    walltime=max(180, mergeAlignmentWalltime/20))
                
                newAlignmentJobAndOutputLs.append(PassingData(jobLs=[addRGJob],
                    file=addRGJob.output))
            
            mergedBamFile = File(os.path.join(reduceOutputDirJob.output, \
                '%s_recal.bam'%(bamFnamePrefix)))
            alignmentMergeJob, bamIndexJob = self.addAlignmentMergeJob(
                alignmentJobAndOutputLs=newAlignmentJobAndOutputLs,
                outputBamFile=mergedBamFile,
                needBAMIndexJob=True,
                parentJobLs=[reduceOutputDirJob],
                walltime=mergeAlignmentWalltime,
                job_max_memory=mergeAlignmentMaxMemory,
                transferOutput=False)
            #2012.9.19 add/copy the alignment file to db-affliated storage
            #add the metric file to AddAlignmentFile2DB.py as well
            #  (to be moved into db-affiliated storage)
            logFile = File(os.path.join(reduceOutputDirJob.output,
                '%s_2db.log'%(bamFnamePrefix)))
            alignment2DBJob = self.addAlignmentFile2DBJob(
                executable=self.AddAlignmentFile2DB,
                inputFile=alignmentMergeJob.output,
                baiFile=bamIndexJob.baiFile,
                individual_alignment_id=new_individual_alignment.id,
                mask_genotype_method_id=self.new_mask_genotype_method_id,
                logFile=logFile,
                data_dir=data_dir,
                otherInputFileList=None,
                parentJobLs=[alignmentMergeJob, bamIndexJob],
                transferOutput=transferOutput,
                sshDBTunnel=self.needSSHDBTunnel,
                commit=True,
                job_max_memory=2000,
                walltime=max(180, mergeAlignmentWalltime/2))
            self.no_of_jobs += 1
            returnData.jobDataLs.append(PassingData(jobLs=[alignment2DBJob],
                file=alignment2DBJob.logFile, \
                fileLs=[alignment2DBJob.logFile]))
        return returnData
    
    def addGATKBaseRecalibratorJob(self, executable=None,
        GenomeAnalysisTKJar=None, inputFile=None, \
        VCFFile=None, interval=None, outputFile=None, \
        refFastaFList=[], parentJobLs=None, extraDependentInputLs=None,
        transferOutput=False, \
        extraArguments=None, job_max_memory=2000, no_of_gatk_threads=1,
        **keywords):
        """
        2013.2.14 upgraded to GATK2's "-T BaseRecalibrator"
        http://gatkforums.broadinstitute.org/discussion/44/base-quality-score-recalibration-bqsr
        java -Xmx4g -jar GenomeAnalysisTK.jar -T BaseRecalibrator
                -I my_reads.bam -R resources/Homo_sapiens_assembly18.fasta
            -knownSites bundle/hg18/dbsnp_132.hg18.vcf \
            -knownSites another/optional/setOfSitesToMask.vcf \
            -o recal_data.grp
        2011-12-4
        inputFile is a bam file.
        outputFile is recalFile (csv file).
    
        java -jar ~/script/vervet/bin/GenomeAnalysisTK-1.0.4705/GenomeAnalysisTK.jar -l INFO
            -R ../../../../NCBI/hs_genome.fasta
            -I 454_vs_hg19.3eQTL.minPerBaseAS0.4.minMapQ125.score2.bam
            -T CountCovariates  -cov ReadGroupCovariate  -cov QualityScoreCovariate  
            -cov CycleCovariate    -cov DinucCovariate
            -recalFile  454_vs_hg19.3eQTL.minPerBaseAS0.4.minMapQ125.score2.recal_data.csv
            -B:mask,VCF 454_vs_hg19.3eQTL.minPerBaseAS0.4.minMapQ125.score2.GATK.vcf
        """
        if executable is None:
            executable = self.BaseRecalibratorJava
        if GenomeAnalysisTKJar is None:
            GenomeAnalysisTKJar = self.GenomeAnalysisTK2Jar
        #GATK job
        memRequirementData = self.getJVMMemRequirment(
            job_max_memory=job_max_memory)
        refFastaFile = refFastaFList[0]
        extraArgumentList = [memRequirementData.memRequirementInStr,
            '-jar', GenomeAnalysisTKJar, "-T BaseRecalibrator",\
            "-I", inputFile, "-R", refFastaFile,\
            "-L", interval, self.defaultGATKArguments,
            "-knownSites:vcf", VCFFile,
            "--out", outputFile]
        if extraArguments:
            extraArgumentList.append(extraArguments)
        if extraDependentInputLs is None:
            extraDependentInputLs=[]
        extraDependentInputLs.extend([inputFile, VCFFile, GenomeAnalysisTKJar] + \
            refFastaFList)
        
        job= self.addGenericJob(executable=executable,
            inputFile=None, outputFile=None, \
            parentJobLs=parentJobLs,
            extraDependentInputLs=extraDependentInputLs,
            extraOutputLs=[outputFile],\
            transferOutput=transferOutput, \
            extraArgumentList=extraArgumentList, \
            job_max_memory=memRequirementData.memRequirement, **keywords)
        
        job.recalFile = outputFile
        return job
    
    def addGATKPrintRecalibratedReadsJob(self, executable=None,
        GenomeAnalysisTKJar=None, inputFile=None, \
        recalFile=None, interval=None, outputFile=None, extraOutputLs=None,
        refFastaFList=[], parentJobLs=None, extraDependentInputLs=None,
        transferOutput=False, \
        extraArguments=None, job_max_memory=2000, no_of_gatk_threads=1,
        needBAMIndexJob=False,\
        walltime=None,**keywords):
        """
        2013.04.09 added extraOutputLs
        2013.2.14 upgraded to GATK2's "-T PrintReads"
            java -jar GenomeAnalysisTK.jar -T PrintReads -R reference.fasta
            -I input.bam -BQSR recalibration_report.grp 
                    -o output.bam
        2012.7.27
        java -jar ~/script/vervet/bin/GenomeAnalysisTK-1.0.4705/GenomeAnalysisTK.jar -l INFO 
            -R ../../../../NCBI/hs_genome.fasta
            -I 454_vs_hg19.3eQTL.minPerBaseAS0.4.minMapQ125.score2.bam 
            -T TableRecalibration 
            --out 454_vs_hg19.3eQTL.minPerBaseAS0.4.minMapQ125.score2.recal.bam 
            -recalFile 454_vs_hg19.3eQTL.minPerBaseAS0.4.minMapQ125.score2.recal_data.csv
        """
        if executable is None:
            executable = self.PrintRecalibratedReadsJava
        if GenomeAnalysisTKJar is None:
            GenomeAnalysisTKJar = self.GenomeAnalysisTK2Jar
        if extraOutputLs is None:
            extraOutputLs = []
        if extraDependentInputLs is None:
            extraDependentInputLs=[]
        memRequirementData = self.getJVMMemRequirment(job_max_memory=job_max_memory)
        #MaxPermSize= min(35000, max(1024, job_max_memory*7/9))
        javaMemRequirement = memRequirementData.memRequirementInStr
        job_max_memory = memRequirementData.memRequirement
        refFastaFile = refFastaFList[0]
        extraArgumentList = [javaMemRequirement, '-jar', GenomeAnalysisTKJar,
            "-T PrintReads",
            "-I", inputFile, "-R", refFastaFile,\
            "-L", interval, self.defaultGATKArguments,\
            "-BQSR", recalFile, "--out", outputFile]
        if extraArguments:
            extraArgumentList.append(extraArguments)
        extraDependentInputLs.extend([inputFile, recalFile, GenomeAnalysisTKJar]\
            + refFastaFList)
        extraOutputLs.insert(0, outputFile)
        #2013.04.30 shuld be ahead of any other output
        job= self.addGenericJob(
            executable=executable, inputFile=None,
            outputFile=None, \
            extraArgumentList=extraArgumentList,
            parentJobLs=parentJobLs,
            extraDependentInputLs=extraDependentInputLs,
            extraOutputLs=extraOutputLs,\
            transferOutput=transferOutput, \
            job_max_memory=job_max_memory,
            walltime=walltime,\
            **keywords)
        if needBAMIndexJob:
            # add the index job on the bam file
            bamIndexJob = self.addBAMIndexJob(
                BuildBamIndexFilesJava=self.BuildBamIndexFilesJava, \
                BuildBamIndexJar=self.BuildBamIndexJar, \
                inputBamF=job.output, parentJobLs=[job], \
                transferOutput=transferOutput, job_max_memory=job_max_memory,
                walltime=walltime/2)
        else:
            bamIndexJob = None
        return job, bamIndexJob
    
    def setup(self, inputVCFData=None, chr2IntervalDataLs=None, **keywords):
        """
        2013.04.07
            register the indel vcfs
        """
        pdata = ParentClass.setup(self, inputVCFData=inputVCFData,
            chr2IntervalDataLs=chr2IntervalDataLs, **keywords)
        
        ##2013.04.09 used in case a corresponding VCFFile for one interval
        #  is not available for BQSR
        biggestChromosome = pdata.chrSizeIDList[0][1]
        self.randomSNPVCFJobDataForBQSR = pdata.chr2VCFJobData[biggestChromosome]
        
        self.indelVCFInputData = self.registerFilesOfInputDir(
            inputDir=self.indelVCFFolder, 
            input_site_handler=self.input_site_handler, \
            checkEmptyVCFByReading=self.checkEmptyVCFByReading,\
            pegasusFolderName="%sIndelVCF"%self.pegasusFolderName)
        #2012.8.26 so that each recalibration will pick up the right vcf
        chr2IndelVCFJobData = {}
        if self.indelVCFInputData:
            for jobData in self.indelVCFInputData.jobDataLs:
                inputF = jobData.file
                chromosome = self.getChrFromFname(os.path.basename(inputF.name))
                chr2IndelVCFJobData[chromosome] = jobData
        pdata.chr2IndelVCFJobData = chr2IndelVCFJobData
        self.chr2IndelVCFJobData = chr2IndelVCFJobData
        
        return pdata
    
    def isThisAlignmentComplete(self, individual_alignment=None, data_dir=None):
        """
        2013.04.09 this is more complicated as it tests the local_realigned
            version of individual_alignment is complete or not.
            not individual_alignment itself
        """
        new_individual_alignment = self.db.copyParentIndividualAlignment(
            parent_individual_alignment_id=individual_alignment.id,\
            mask_genotype_method_id=self.new_mask_genotype_method_id,\
            data_dir=self.data_dir, local_realigned=1)
        return self.db.isThisAlignmentComplete(
            individual_alignment=new_individual_alignment, data_dir=data_dir)
    
    def registerCustomExecutables(self):
        """
        2011-11-28
        """
        ParentClass.registerCustomExecutables(self)
        
        #samtools is only used for select alignment, which is very fast,
        #   increase the clustering 
        self.setExecutableClusterSize(executable=self.samtools,
            clusterSizeMultiplier=1)
        
        self.registerOneExecutable(path=self.javaPath, 
            name='BaseRecalibratorJava', clusterSizeMultiplier=1)
        self.registerOneExecutable(path=self.javaPath, 
            name='PrintRecalibratedReadsJava', clusterSizeMultiplier=1)
        self.registerOneExecutable(path=self.javaPath, 
            name='RealignerTargetCreatorJava', clusterSizeMultiplier=1)
        self.registerOneExecutable(path=self.javaPath, 
            name='IndelRealignerJava', clusterSizeMultiplier=1)


if __name__ == '__main__':
    main_class = AlignmentReadBaseQualityRecalibration
    po = ProcessOptions(sys.argv, main_class.option_default_dict,
        error_doc=main_class.__doc__)
    instance = main_class(**po.long_option2value)
    instance.run()
