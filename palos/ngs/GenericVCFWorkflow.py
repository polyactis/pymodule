#!/usr/bin/env python3
"""
Description:
    20120509 Parent of other VCF workflows.

Examples:
    #2012.5.11 convert alignment read group (sample id) into UCLAID
    %s -I FilterVCF_VRC_SK_Nevis_FilteredSeq_top1000Contigs.2012.5.6_trioCaller.2012.5.8T21.42/trioCaller_vcftoolsFilter/ 
        -o dags/SampleIDInUCLAID_FilterVCF_VRC_SK_Nevis_FilteredSeq_top1000Contigs.2012.5.6_trioCaller.2012.5.8.xml 
        -u yh -y4 -l hcondor -j hcondor  -z localhost
        -e /u/home/eeskin/polyacti/
        -t /u/home/eeskin/polyacti/NetworkData/vervet/db/
        -D /u/home/eeskin/polyacti/NetworkData/vervet/db/ 
    
    # 2012.5.10 subset + convert-2-plink.
    # run on hoffman2 condor, minMAC=1 (-n 1), minMAF=0.1 (-f 0.1), maxSNPMissingRate=0 (-L 0)
    #    (turn on checkEmptyVCFByReading, --checkEmptyVCFByReading)
    %s -I FilterVCF_VRC_SK_Nevis_FilteredSeq_top1000Contigs.2012.5.6_trioCaller.2012.5.8T21.42/trioCaller_vcftoolsFilter/
        -o dags/SubsetTo36RNASamplesAndPlink_FilterVCF_VRC_SK_Nevis_FilteredSeq_top1000Contigs.2012.5.6_trioCaller.2012.5.8.xml
        -i ~/script/vervet/data/RNADevelopment_eQTL/36monkeys.phenotypes.txt
        -w ~/script/vervet/data/RNADevelopment_eQTL/36monkeys.inAlignmentReadGroup.tsv
        -n1 -f 0.1 -L 0 -y3 --checkEmptyVCFByReading
        -l hcondor -j hcondor  -u yh -z localhost --needSSHDBTunnel
        -e /u/home/eeskin/polyacti/
        -D /u/home/eeskin/polyacti/NetworkData/vervet/db/
        -t /u/home/eeskin/polyacti/NetworkData/vervet/db/
    
    # 2012.7.16 convert a folder of VCF files into plink, need the db tunnel (--needSSHDBTunnel) for output pedigree in tfam
    # "-V 90 -x 100" are used to restrict contig IDs between 90 and 100.
    %s -I FilterVCF_VRC_SK_Nevis_FilteredSeq_top1000Contigs.MAC10.MAF.05_trioCaller/trioCaller_vcftoolsFilter/ 
        -o dags/ToPlinkFilterVCF_VRC_SK_Nevis_FilteredSeq_top1000Contigs.MAC10.MAF.05.xml
        -y 2 --checkEmptyVCFByReading
        -l condorpool -j condorpool
        -u yh -z uclaOffice  -C 4 --needSSHDBTunnel
        #-V 90 -x 100 
    
    # 2012.7.25 calculate haplotype distance & majority call support stats
    %s -I AlignmentToTrioCall_VRC_FilteredSeq.2012.7.21T0248_VCFWithReplicates/
        -o dags/GetReplicateHaplotypeStat_TrioCall_VRC_FilteredSeq.2012.7.21T0248_VCFWithReplicats.xml
        -y 5 --checkEmptyVCFByReading -l condorpool -j condorpool -u yh -z uclaOffice  -C 1 -a 524
    
    # 2012.8.20 convert method 16 to yu format (-y 6 works for generic VCF, -y 7 adds sample ID conversion first)
    %s -I ~/NetworkData/vervet/db/genotype_file/method_16/
        -o dags/VCF2YuFormat/VCF2YuFormat_Method16.xml
        -y 7 --checkEmptyVCFByReading  -l hcondor -j hcondor  -u yh -z localhost --needSSHDBTunnel -C 2
        -D /u/home/eeskin/polyacti/NetworkData/vervet/db/  -t /u/home/eeskin/polyacti/NetworkData/vervet/db/

    # 2012.8.30 combine multiple VCF into one
    # -s .. is optional. if given, the combined VCF will be added into db.
    %s -I ~/NetworkData/vervet/db/genotype_file/method_10/ -o dags/GenericVCFWorkflow/MultiVCF2OneFile_Method10.xml
        -y 9  -l hcondor -j hcondor -u yh -z localhost --needSSHDBTunnel -C 1
        -s 16HCSAMtoolsMinDP1_2FoldDepth_minMAC8_maxSNPMissing0
        -D /u/home/eeskin/polyacti/NetworkData/vervet/db/ -t /u/home/eeskin/polyacti/NetworkData/vervet/db/
    
"""
import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0], \
    sys.argv[0], sys.argv[0])

import csv, copy
from palos import ProcessOptions, getListOutOfStr, PassingData, \
    figureOutDelimiter, getColName2IndexFromHeader, utils
from palos.ngs.io.VCFFile import VCFFile
from pegaflow.DAX3 import File, Executable, PFN
import pegaflow
#from palos.pegasus.AbstractVCFWorkflow import AbstractVCFWorkflow
from palos.db import SunsetDB
from . AbstractNGSWorkflow import AbstractNGSWorkflow

ParentClass = AbstractNGSWorkflow
class GenericVCFWorkflow(ParentClass):
    __doc__ = __doc__
    option_default_dict = copy.deepcopy(ParentClass.option_default_dict)
    option_default_dict.update({
        ('individualUCLAIDFname', 0, ): [None, 'i', 1, 
            'a file containing individual ucla_id in each row. one column with header UCLAID. ', ],\
        ('vcfSampleIDFname', 0, ): [None, 'w', 1, 
            'a file containing the sample ID (a composite ID including ucla_id) each row. '
            'any samples not in this file will be removed in subset VCF run_type (1, 3)'
            'You can also use individualUCLAIDFname to specify the sample IDs (UCLAID). '
            'Their composite IDs will be inferred from individualUCLAIDFname + first VCF file header.', ],
        ('vcf2Dir', 0, ): ['', '', 1, 
            'the 2nd input folder that contains vcf or vcf.gz files.', ],
        ('vcfSubsetPath', 1, ): ["bin/vcftools/vcf-subset", '', 1, 
            'path to the vcf-subset program'],
        ('run_type', 1, int): [1, 'y', 1, 'which run_type to run. \
            1: subset VCF based on input file containing sample IDs;\
            2: convert to plink format; \
            3: subset + convert-2-plink. MAC & MAF & maxSNPMissingRate applied in the convert-to-plink step.\
            4: ConvertAlignmentReadGroup2UCLAIDInVCF jobs.\
            5: addMergeVCFReplicateHaplotypesJobs to get haplotype distance & majority call support stats.\
            6: VCF2YuFormatJobs, \
            7: ConvertAlignmentReadGroup2UCLAIDInVCF + VCF2YuFormatJobs, \
            8: ConvertAlignmentReadGroup2UCLAIDInVCF + convert to plink format, \
            9: combine all single-chromosome VCF into one. \
            ?: Combine VCF files from two input folder, chr by chr. (not done yet. check CheckTwoVCFOverlapPipeline.py for howto)', ],\
        ("minMAC", 0, int): [None, 'n', 1, 'minimum MinorAlleleCount (by chromosome)'],\
        ("minMAF", 0, float): [None, 'f', 1, 
            'minimum MinorAlleleFrequency (by chromosome)'],\
        ("maxSNPMissingRate", 0, float): [1.0, 'L', 1, 
            'maximum SNP missing rate in one vcf (denominator is #chromosomes)'],\
        ('genotypeMethodShortName', 0, ):[None, 's', 1, 
            'column short_name of GenotypeMethod table,'
            'will be created if not present in db.'
            'for run_type 9, if given the file would be added into db.'],\
        })

    def __init__(self,
        vcfSubsetPath="bin/vcftools/vcf-subset",
        plinkPath="bin/plink",
        **keywords):
        """
        """
        self.vcfSubsetPath = vcfSubsetPath
        self.plinkPath = plinkPath
        # Set before ParentClass.__init__()
        self.pathToInsertHomePathList.extend(['vcfSubsetPath', 'plinkPath'])
        ParentClass.__init__(self, **keywords)
    
    def addPlinkJob(self, executable=None, inputFileList=None,
        parentPlinkJob=None,\
        tpedFile=None, tfamFile=None,\
        pedFile=None, famFile=None, mapFile=None,
        bedFile=None, bimFile=None,
        inputFnamePrefix=None, inputOption='--file', \
        outputFnamePrefix=None, outputOption='--out',\
        makeBED=False, calculateMendelError=False, checkSex=False,
        LDPruneWindowSize=100, LDPruneWindowShiftSize=5,
        LDPruneByPairwiseR2=False, LDPruneMinR2=0.1,\
        LDPruneByRegression=False, LDPruneMinVarianceInflationFactor=2,
        estimatePairwiseGenomeWideIBD=False,
        estimatePairwiseGenomeWideIBDFreqFile=None,
        extractSNPFile=None, recodeOutput=False,
        recodeTransposeOutput=False,
        estimateAlleFrequency=False,
        mergeListFile=None,\
        parentJobLs=None, extraDependentInputLs=None,
        transferOutput=False,
        extraArguments=None, extraArgumentList=None,
        extraOutputLs =None,
        job_max_memory=2000, **keywords):
        """
        i.e.

        bedFnamePrefix = os.path.join(topOutputDir, '%s_bed'%(commonPrefix))
        convertSingleTPED2BEDJob = self.addPlinkJob(executable=self.plink, \
            inputFileList=[],
            tpedFile=modifyTPEDJob.output, tfamFile=tfamJob.tfamFile,\
            outputFnamePrefix=bedFnamePrefix, outputOption='--out',\
            makeBED=True, \
            extraDependentInputLs=None, transferOutput=transferOutput, \
            extraArguments=None, job_max_memory=2000,\
            parentJobLs = convertSingleTPED2BEDParentJobLs)


        convertMergedTPED2BEDJob = self.addPlinkJob(executable=self.plink,
                inputFileList=[tpedFileMergeJob.output, tfamJob.tfamFile], \
                inputFnamePrefix=mergedPlinkFnamePrefix, inputOption='--tfile', \
            outputFnamePrefix=mergedPlinkBEDFnamePrefix, outputOption='--out',\
            makeBED=True, \
            extraDependentInputLs=None, transferOutput=transferOutput, \
            extraArguments=None, job_max_memory=2000, parentJobLs=[mergedOutputDirJob, 
                tpedFileMergeJob, tfamJob])

        mendelFnamePrefix = os.path.join(setupData.mapDirJob.output, '%s'%(commonPrefix))
        if inputJob.output.name[-4:]=='tped':
            inputJob.tfamFile = tfamJob.tfamFile
        plinkMendelJob = self.addPlinkJob(executable=self.plink, \
                parentPlinkJob=inputJob,\
                outputFnamePrefix=mendelFnamePrefix, outputOption='--out',\
                calculateMendelError=True, \
                extraDependentInputLs=None,
                transferOutput=transferOneContigPlinkOutput, \
                extraArguments=None, job_max_memory=2000,\
                parentJobLs =[setupData.mapDirJob, tfamJob]+ jobData.jobLs)

        for plink mendel, LD-prune and other jobs, add extraArguments="--allow-no-sex"
            to include individuals without sex

        2013.07.25 added parentPlinkJob (returned from this function),
            and parse input from that job
        2013.07.24 added argument recodeTransposeOutput (--recode --transpose)
        2012.8.28
            add argument
                estimateAlleFrequency, estimate frequency of input file.
                     "--nonfounders" could be added as well.
                estimatePairwiseGenomeWideIBDFreqFile, is the file from which 
                    IBD check could draw frequency (rather than estimate from founders)

        2012.8.9
            inputFileList is a list of pegasus Files (.ped, .fam, or .tped,
                 .tfam, etc.) or could be supplied individually.

            inputOption could be, "--file" for .ped .map ; "--tfile" for .tped,
                 .tfam; or '--bfile' for .bed, .fam, .bim

            if extractSNPFile or mergeListFile is given, either recodeOutput
                 or makeBED have to be on. otherwise, no output.
            http://pngu.mgh.harvard.edu/~purcell/plink/index.shtml
        """
        if extraDependentInputLs is None:
            extraDependentInputLs = []
        if inputFileList:
            extraDependentInputLs.extend(inputFileList)

        if extraArgumentList is None:
            extraArgumentList = []
        if extraOutputLs is None:
            extraOutputLs = []
        key2ObjectForJob = {}

        if parentPlinkJob:
            if bedFile is None:
                bedFile = getattr(parentPlinkJob, 'bedFile', None)
            if famFile is None:
                famFile = getattr(parentPlinkJob, 'famFile', None)
            if bimFile is None:
                bimFile = getattr(parentPlinkJob, 'bimFile', None)
            if tpedFile is None:
                tpedFile = getattr(parentPlinkJob, 'tpedFile', None)
            if tfamFile is None:
                tfamFile = getattr(parentPlinkJob, 'tfamFile', None)
            if mapFile is None:
                mapFile = getattr(parentPlinkJob, 'mapFile', None)
            if pedFile is None:
                pedFile = getattr(parentPlinkJob, 'pedFile', None)
            if famFile is None:
                famFile = getattr(parentPlinkJob, 'famFile', None)

        if inputOption and inputFnamePrefix:
            extraArgumentList.extend([inputOption, inputFnamePrefix])
        if tpedFile:
            extraDependentInputLs.append(tpedFile)
            extraArgumentList.extend(["--tped", tpedFile])
        if tfamFile:
            extraDependentInputLs.append(tfamFile)
            extraArgumentList.extend(["--tfam", tfamFile])
        if pedFile:
            extraDependentInputLs.append(pedFile)
            extraArgumentList.extend(["--ped", pedFile])
        if famFile:
            extraDependentInputLs.append(famFile)
            extraArgumentList.extend(["--fam", famFile])
        if mapFile:
            extraDependentInputLs.append(mapFile)
            extraArgumentList.extend(["--map", mapFile])
        if bedFile:
            extraDependentInputLs.append(bedFile)
            extraArgumentList.extend(["--bed", bedFile])
        if bimFile:
            extraDependentInputLs.append(bimFile)
            extraArgumentList.extend(["--bim", bimFile])

        if outputFnamePrefix and outputOption:
            extraArgumentList.extend([outputOption, outputFnamePrefix])
        else:
            outputFnamePrefix = 'plink'

        suffixAndNameTupleList = []
            # a list of tuples , in each tuple, 1st element is the suffix.
            #  2nd element is the proper name of the suffix.
            #job.$nameFile will be the way to access the file.
            #if 2nd element (name) is missing, suffix[1:].replace('.', '_')
            #  is the name (dot replaced by _)
        if makeBED:
            extraArgumentList.append('--make-bed')
            suffixAndNameTupleList.extend([['.bed',], ('.fam',), ['.bim',]])
            #, binary map file, is excluded for now
        if calculateMendelError:
            extraArgumentList.append('--mendel')
            suffixAndNameTupleList.extend([('.mendel',), ('.imendel',), 
                ('.fmendel',), ('.lmendel',)])
            #its output is not tab-delimited. rather it's space (multi) delimited.
        if checkSex:
            extraArgumentList.append('--check-sex')
            suffixAndNameTupleList.extend([('.sexcheck',), ('.hh', )])
            #.sexcheck file is accessible as job.sexcheckFile.
            #.hh is heterozygous haplotype genotypes
        if LDPruneByPairwiseR2:
            extraArgumentList.append('--indep-pairwise %s %s %s'%(\
                LDPruneWindowSize, LDPruneWindowShiftSize, LDPruneMinR2))
            suffixAndNameTupleList.extend([('.prune.in',), ('.prune.out',)])
            #".prune.in" is accessible as job.prune_inFile
        if LDPruneByRegression:
            extraArgumentList.append('--indep %s %s %s'%(LDPruneWindowSize,\
                LDPruneWindowShiftSize, LDPruneMinVarianceInflationFactor))
            suffixAndNameTupleList.extend([('.prune.in',), ('.prune.out',)])
            #".prune.in" is accessible as job.prune_inFile
        if estimatePairwiseGenomeWideIBD:
            extraArgumentList.append('--genome')
            suffixAndNameTupleList.extend([('.genome',)])
            #.genome is accessible as job.genomeFile
            if estimatePairwiseGenomeWideIBDFreqFile:
                extraArgumentList.extend(['--read-freq', \
                    estimatePairwiseGenomeWideIBDFreqFile])
                extraDependentInputLs.append(estimatePairwiseGenomeWideIBDFreqFile)
        if extractSNPFile:
            extraArgumentList.extend(['--extract', extractSNPFile])
            extraDependentInputLs.append(extractSNPFile)
        if recodeOutput:
            extraArgumentList.extend(['--recode',])
            suffixAndNameTupleList.extend([('.ped',), ('.map',)])
        if recodeTransposeOutput:
            extraArgumentList.extend(['--recode', "--transpose"])
            suffixAndNameTupleList.extend([('.tped',), ('.tfam',)])
        if estimateAlleFrequency:
            extraArgumentList.append('--freq')
            suffixAndNameTupleList.extend([('.frq',)])

        if mergeListFile:
            extraArgumentList.extend(['--merge-list', mergeListFile])
            extraDependentInputLs.append(mergeListFile)
        if extraArguments:
            extraArgumentList.append(extraArguments)


        self.setupMoreOutputAccordingToSuffixAndNameTupleList(
            outputFnamePrefix=outputFnamePrefix, \
            suffixAndNameTupleList=suffixAndNameTupleList, \
            extraOutputLs=extraOutputLs, key2ObjectForJob=key2ObjectForJob)
        #2013.07.24 add it in the end
        #2012.8.10 left in the folder dying
        logFile = File('%s.log'%(outputFnamePrefix))
        extraOutputLs.append(logFile)

        job= self.addGenericJob(executable=executable, inputFile=None,
            outputFile=None, \
            parentJobLs=parentJobLs,
            extraDependentInputLs=extraDependentInputLs,
            extraOutputLs=extraOutputLs,
            transferOutput=transferOutput,
            extraArgumentList=extraArgumentList,
            key2ObjectForJob=key2ObjectForJob, job_max_memory=job_max_memory,
            **keywords)
        return job
    
    def addVCF2PlinkJobs(self, inputData=None, db_main=None,
        minMAC=None, minMAF=None,
        maxSNPMissingRate=None, transferOutput=True,\
        maxContigID=None, outputDirPrefix="", outputPedigreeAsTFAM=False,
        outputPedigreeAsTFAMInputJobData=None, \
        treatEveryOneIndependent=False,\
        returnMode=3, ModifyTPEDRunType=1, chr_id2cumu_chr_start=None,\
        addUngenotypedDuoParents=False):
        """
    returnMode
        1: only the final merged binary .bed , .fam file and its generation job(s)
        2: only the individual contig/chromosome (whatever in inputDat.jobDataLs)
             binary .bed, .fam files and the associated jobs
        3: 1 & 2 (all individual binary .bed jobs&files + the last merged file/job)
        4: the individual contig/chr non-binary (TPED) job data
         (for Mark mendel error genotype as missing)
        5: 
        
        2013.07.18
            added argument addUngenotypedDuoParents
                for mendel error detection, if an ungenotyped parent in a duo
                 is not present in the genotype file (PED/TPED/BED),
                then plink won't look for its mendel inconsistency 
        2013.02
        2013.1.29 added returnMode 4
        2012.10.22
            change transferOutput of outputPedigreeInTFAMJob to True
        2012.9.13
            add argument treatEveryOneIndependent for
             OutputVRCPedigreeInTFAMGivenOrderFromFile.
        2012.8.20 add outputPedigreeAsTFAMInputJobData, split from input_data.
            outputPedigreeAsTFAMInputJobData.vcfFile must use
             individual_alignment.read_group as sample ID.
            useful in the case of that VCF files have been converted into UCLA IDs.
        2012.8.9
            add argument
                outputPedigreeAsTFAM
                returnMode
            add plink binary 
        2012.7.19 
            add a modifyTPEDJob that modify 2nd column (snp-id) of tped output from default 0 to chr_pos.
                argument ModifyTPEDRunType.
                1: modify snp_id (2nd-column) = chr_phyiscalPosition,\
                2: snp_id=chr_physicalPosition (original data), chr (1st column) = X (chromosome X, for sex check by plink), pos += positionStartBase.,\
                3: snp_id=chr_physicalPosition (original data), chr (1st column) = newChr, pos += positionStartBase
            
            added a GzipSubworkflow in the end to gzip the final merged tped file
            all previous intermediate files are not transferred.
        2012.5.9
        """
        sys.stderr.write("Adding VCF2plink jobs for %s vcf files ... "%(\
            len(inputData.jobDataLs)))
        
        topOutputDir = "%sVCF2Plink"%(outputDirPrefix)
        topOutputDirJob = self.addMkDirJob(outputDir=topOutputDir)
        
        mergedOutputDir = "%sVCF2PlinkMerged"%(outputDirPrefix)
        mergedOutputDirJob = self.addMkDirJob(outputDir=mergedOutputDir)
        
        mergedPlinkFnamePrefix = os.path.join(mergedOutputDir, 'merged')
        mergedTPEDFile = File('%s.tped'%(mergedPlinkFnamePrefix))
        #each input has no header
        tpedFileMergeJob = self.addStatMergeJob(
            statMergeProgram=self.mergeSameHeaderTablesIntoOne, \
            outputF=mergedTPEDFile, transferOutput=False,
            parentJobLs=[mergedOutputDirJob], \
            extraArguments='--noHeader')
        
        if outputPedigreeAsTFAMInputJobData is None:
            outputPedigreeAsTFAMInputJobData = inputData.jobDataLs[0]
        if outputPedigreeAsTFAM and outputPedigreeAsTFAMInputJobData:
            jobData = outputPedigreeAsTFAMInputJobData
            inputF = jobData.vcfFile
            outputFile = File(os.path.join(mergedOutputDir, 'pedigree.tfam'))
            outputPedigreeInTFAMJob = self.addOutputVRCPedigreeInTFAMGivenOrderFromFileJob(
                executable=self.OutputVRCPedigreeInTFAMGivenOrderFromFile, \
                inputFile=inputF, outputFile=outputFile,
                treatEveryOneIndependent=treatEveryOneIndependent,\
                addUngenotypedDuoParents=addUngenotypedDuoParents,\
                parentJobLs=[mergedOutputDirJob]+jobData.jobLs,
                extraDependentInputLs=[], transferOutput=True, \
                extraArguments=None, job_max_memory=2000,
                sshDBTunnel=self.needSSHDBTunnel)
            outputPedigreeInTFAMJob.tfamFile = outputPedigreeInTFAMJob.output
            #so that it looks like a vcf2plinkJob (vcftools job)
        else:
            outputPedigreeInTFAMJob = None
        
        returnData = PassingData()
        returnData.jobDataLs = []
        returnData.tfamJob = None
        #2013.07.25 family file for tped file 
        returnData.famJob = None
        #2013.07.25 family file for bed file
        
        for i in range(len(inputData.jobDataLs)):
            jobData = inputData.jobDataLs[i]
            inputF = jobData.vcfFile
            inputFBaseName = os.path.basename(inputF.name)
            chr_id = self.getChrFromFname(inputFBaseName)
            if maxContigID:
                contig_id = self.getContigIDFromFname(inputFBaseName)
                try:
                    contig_id = int(contig_id)
                    if contig_id>maxContigID:
                        #skip the small contigs
                        continue
                except:
                    sys.stderr.write('Except type: %s\n'%repr(sys.exc_info()))
                    import traceback
                    traceback.print_exc()
            commonPrefix = inputFBaseName.split('.')[0]
            outputFnamePrefix = os.path.join(topOutputDir, '%s'%(commonPrefix))
            if i ==0:
                #need at least one tfam file. 
                transferOneContigPlinkOutput = True
            else:
                transferOneContigPlinkOutput = False
            i += 1
            vcf2plinkJob = self.addFilterJobByvcftools(
                vcftoolsWrapper=self.vcftoolsWrapper,
                inputVCFF=inputF, \
                outputFnamePrefix=outputFnamePrefix, \
                parentJobLs=[topOutputDirJob]+jobData.jobLs, \
                snpMisMatchStatFile=None, \
                minMAC=minMAC, minMAF=minMAF, \
                maxSNPMissingRate=maxSNPMissingRate,\
                extraDependentInputLs=[jobData.tbi_F],
                outputFormat='--plink-tped',
                transferOutput=transferOneContigPlinkOutput)
            #2013.07.19
            if addUngenotypedDuoParents and outputPedigreeInTFAMJob:
                appendExtraIndividualsJob = self.addAbstractMapperLikeJob(
                    executable=self.AppendExtraPedigreeIndividualsToTPED, \
                    inputF=vcf2plinkJob.tpedFile, \
                    outputF=File(os.path.join(topOutputDir, '%s_extraIndividuals.tped'%(commonPrefix))),
                    parentJobLs=[vcf2plinkJob, outputPedigreeInTFAMJob, topOutputDirJob],
                    transferOutput=False, job_max_memory=200,
                    extraArgumentList=["--tfamFname", outputPedigreeInTFAMJob.tfamFile],
                    extraDependentInputLs=[outputPedigreeInTFAMJob.tfamFile])
                modifyTPEDParentJobLs = [appendExtraIndividualsJob]
                modifyTPEDJobInputFile = appendExtraIndividualsJob.output
            else:
                if addUngenotypedDuoParents and outputPedigreeInTFAMJob is None:
                    message = "Warning: addUngenotypedDuoParents (%s) is True but outputPedigreeInTFAMJob (%s, outputPedigreeAsTFAM=%s) is None. so addUngenotypedDuoParents is effectively False."%\
                            (addUngenotypedDuoParents, outputPedigreeInTFAMJob, outputPedigreeAsTFAM)
                    utils.pauseForUserInput(message=message, continueAnswerSet=None, exitType=3)
                modifyTPEDParentJobLs = [vcf2plinkJob]
                modifyTPEDJobInputFile = vcf2plinkJob.tpedFile
            
            #2012.7.20 modify the TPED 2nd column, to become chr_pos (rather than 0)
            modifyTPEDFnamePrefix = os.path.join(topOutputDir, '%s_SNPID_M'%(commonPrefix))
            outputF = File('%s.tped'%(modifyTPEDFnamePrefix))
            modifyTPEDJobExtraArguments = "--run_type %s "%(ModifyTPEDRunType)
            if ModifyTPEDRunType==3 and chr_id2cumu_chr_start:
                newChrID, newCumuStart = chr_id2cumu_chr_start.get(chr_id, (1,0))
                modifyTPEDJobExtraArguments += ' --newChr %s --positionStartBase %s '%(newChrID, newCumuStart)
            modifyTPEDJob = self.addAbstractMapperLikeJob(
                executable=self.ModifyTPED,
                inputF=modifyTPEDJobInputFile, outputF=outputF, \
                parentJobLs=modifyTPEDParentJobLs, transferOutput=False,
                job_max_memory=200,\
                extraArguments=modifyTPEDJobExtraArguments,
                extraDependentInputLs=[])
            
            #add output to the tped merge job
            self.addInputToMergeJob(tpedFileMergeJob, \
                inputF=modifyTPEDJob.output, \
                parentJobLs=[modifyTPEDJob])

            if outputPedigreeInTFAMJob is None:
                tfamJob = vcf2plinkJob
                convertSingleTPED2BEDParentJobLs = [modifyTPEDJob, vcf2plinkJob]
            else:
                tfamJob = outputPedigreeInTFAMJob
                convertSingleTPED2BEDParentJobLs = [modifyTPEDJob, outputPedigreeInTFAMJob]
            
            if returnData.tfamJob is None:
                returnData.tfamJob = tfamJob
            
            if returnMode==4:
                returnData.jobDataLs.append(PassingData(jobLs=[modifyTPEDJob], file=modifyTPEDJob.output, \
                    fileLs=modifyTPEDJob.outputLs))
            elif returnMode==2 or returnMode==3:
                #convert single plink tped file into binary bed file
                #add it to 
                bedFnamePrefix = os.path.join(topOutputDir, '%s_bed'%(commonPrefix))
                convertSingleTPED2BEDJob = self.addPlinkJob(
                    executable=self.plinkConvert, inputFileList=[], 
                    tpedFile=modifyTPEDJob.output, tfamFile=tfamJob.tfamFile,\
                    inputFnamePrefix=None, inputOption=None, \
                    outputFnamePrefix=bedFnamePrefix, outputOption='--out',\
                    makeBED=True, \
                    extraDependentInputLs=None, transferOutput=transferOutput, \
                    extraArguments=None, job_max_memory=2000,\
                    parentJobLs = convertSingleTPED2BEDParentJobLs)
                returnData.jobDataLs.append(PassingData(jobLs=[convertSingleTPED2BEDJob], file=convertSingleTPED2BEDJob.bedFile, \
                                            fileLs=convertSingleTPED2BEDJob.outputLs))
                if returnData.famJob is None:
                    returnData.famJob = convertSingleTPED2BEDJob
        
        if returnMode==1 or returnMode==3:
            #convert merged plain tped file into binary bed files
            mergedPlinkBEDFnamePrefix = os.path.join(mergedOutputDir, 'mergedPlinkBED')
            convertMergedTPED2BEDJob = self.addPlinkJob(
                executable=self.plinkMerge, inputFileList=[], \
                tpedFile=tpedFileMergeJob.output, tfamFile=tfamJob.tfamFile,\
                inputFnamePrefix=None, inputOption=None, \
                outputFnamePrefix=mergedPlinkBEDFnamePrefix,
                outputOption='--out',\
                makeBED=True, \
                extraDependentInputLs=None, transferOutput=transferOutput,
                extraArguments=None, job_max_memory=2000,
                parentJobLs=[mergedOutputDirJob, tpedFileMergeJob, tfamJob])
            returnData.jobDataLs.append(PassingData(jobLs=[convertMergedTPED2BEDJob],
                file=convertMergedTPED2BEDJob.bedFile,
                fileLs=convertMergedTPED2BEDJob.outputLs))
            if returnData.famJob is None:
                returnData.famJob = convertMergedTPED2BEDJob
        ##2012.8.9 gzip workflow is not needed anymore as binary bed is used instead.
        ##2012.7.21 gzip the final output
        gzipInputData = PassingData()
        gzipInputData.jobDataLs = []
        gzipInputData.jobDataLs.append(PassingData(jobLs=[tpedFileMergeJob],
            file=tpedFileMergeJob.output, \
            fileLs=tpedFileMergeJob.outputLs))
        self.addGzipSubWorkflow(inputData=gzipInputData, transferOutput=transferOutput,\
            outputDirPrefix="gzipMergedTPED")
        sys.stderr.write("%s jobs.\n"%(self.no_of_jobs))
        #2013.1.29 return the un-gzipped data so that downstream sub-workflows could work on un-gzipped files
        return returnData
    
    def addVCF2YuFormatJobs(self, inputData=None, transferOutput=True,\
                        maxContigID=None, outputDirPrefix="", \
                        returnMode=1):
        """
        2012.8.20
            argument
                returnMode
                    1=only the final merged file and its generation job(s)
                    2=only the individual contig/chromosome (whatever in inputDat.jobDataLs) converted files and conversion jobs
                    3= 1 & 2 (all individual input binary .bed job&file + the last merging job/file)
        """
        sys.stderr.write("Adding VCF2YuFormat jobs for %s vcf files ... "%(len(inputData.jobDataLs)))
        
        topOutputDir = "%sVCF2BjarniFormat"%(outputDirPrefix)
        topOutputDirJob = self.addMkDirJob(outputDir=topOutputDir)
        
        mergeOutputDir = "%sVCF2YuFormat"%(outputDirPrefix)
        mergeOutputDirJob = self.addMkDirJob(outputDir=mergeOutputDir)
        
        mergeFnamePrefix = os.path.join(mergeOutputDir, 'merged')
        mergeFile = File('%s.csv'%(mergeFnamePrefix))
        #each input has no header
        mergeFileJob = self.addStatMergeJob(statMergeProgram=self.mergeSameHeaderTablesIntoOne, \
                            outputF=mergeFile, transferOutput=False, parentJobLs=[mergeOutputDirJob])
        
        returnData = PassingData()
        returnData.jobDataLs = []
        for i in range(len(inputData.jobDataLs)):
            jobData = inputData.jobDataLs[i]
            inputF = jobData.vcfFile
            inputFBaseName = os.path.basename(inputF.name)
            chr_id = self.getChrFromFname(inputFBaseName)
            if maxContigID:
                contig_id = self.getContigIDFromFname(inputFBaseName)
                try:
                    contig_id = int(contig_id)
                    if contig_id>maxContigID:
                        #skip the small contigs
                        continue
                except:
                    sys.stderr.write('Except type: %s\n'%repr(sys.exc_info()))
                    import traceback
                    traceback.print_exc()
            commonPrefix = inputFBaseName.split('.')[0]
            outputFnamePrefix = os.path.join(topOutputDir, '%s'%(commonPrefix))
            if i ==0:
                #need at least one tfam file. 
                transferOneContigPlinkOutput = True
            else:
                transferOneContigPlinkOutput = False
            i += 1
            bjarniFormatOutput = File('%s.csv'%(outputFnamePrefix))
            vcf2BjarniFormatJob = self.addGenericJob(
                executable=self.ConvertVCF2BjarniFormat,
                inputFile=inputF, inputArgumentOption="-i", \
                outputFile=bjarniFormatOutput, outputArgumentOption="-o", \
                parentJobLs=[topOutputDirJob] + jobData.jobLs,
                extraDependentInputLs=None, extraOutputLs=None, \
                transferOutput=transferOneContigPlinkOutput, \
                extraArguments="--outputDelimiter ,", extraArgumentList=None,
                job_max_memory=2000, sshDBTunnel=None, \
                key2ObjectForJob=None)
            
            #add output to some reduce job
            self.addInputToMergeJob(mergeFileJob, \
                inputF=vcf2BjarniFormatJob.output, \
                parentJobLs=[vcf2BjarniFormatJob])
            
            
            if returnMode==2 or returnMode==3:
                returnData.jobDataLs.append(PassingData(jobLs=[vcf2BjarniFormatJob],
                    file=vcf2BjarniFormatJob.output, \
                    fileLs=vcf2BjarniFormatJob.outputLs))
            

        if returnMode==1 or returnMode==3:
            #conver the merged bjarni format into yu format
#			outputFile = File(os.path.join(mergeOutputDir, 'merge.tsv'))
#			bjar2YuFormatJob = self.addGenericJob(
#               executable=self.ConvertBjarniSNPFormat2Yu,
#               inputFile=mergeFileJob.output, inputArgumentOption="-i", \
#				outputFile=outputFile, outputArgumentOption="-o", \
#				parentJobLs=[mergeFileJob], extraDependentInputLs=None, extraOutputLs=None, \
#				transferOutput=transferOutput, \
#				extraArguments=None, extraArgumentList=None, job_max_memory=8000,
#               sshDBTunnel=None, \
#				key2ObjectForJob=None)
#			no_of_jobs += 1
#			returnData.jobDataLs.append(PassingData(jobLs=[bjar2YuFormatJob],
#               file=bjar2YuFormatJob.output, \
#				fileLs=bjar2YuFormatJob.outputLs))
            pass	#too much memory
        sys.stderr.write("%s jobs. Done.\n"%(self.no_of_jobs))
        return returnData
    
    def addVCFSubsetJob(self, executable=None, vcfSubsetPath=None, 
        sampleIDFile=None,
        inputVCF=None, outputF=None, \
        parentJobLs=None, transferOutput=True, job_max_memory=200,\
        extraArguments=None, extraDependentInputLs=None, **keywords):
        """
        2012.10.10 use addGenericJob
        """
        extraArgumentList = [vcfSubsetPath, sampleIDFile, inputVCF, outputF]
        if extraDependentInputLs is None:
            extraDependentInputLs = []
        extraDependentInputLs.append(sampleIDFile)
        extraDependentInputLs.append(inputVCF)
        extraOutputLs = [outputF]
        job = self.addGenericJob(executable=executable, inputFile=None, \
            outputFile=None, \
            parentJobLs=parentJobLs,
            extraDependentInputLs=extraDependentInputLs,
            extraOutputLs=extraOutputLs, \
            transferOutput=transferOutput, \
            extraArguments=extraArguments, extraArgumentList=extraArgumentList,
            job_max_memory=job_max_memory, sshDBTunnel=None, \
            key2ObjectForJob=None, **keywords)
        return job

    def addVCFSubsetJobs(self, inputData=None, db_main=None,
        sampleIDFile=None, transferOutput=True,\
        refFastaFList=None, GenomeAnalysisTKJar=None,\
        maxContigID=None, outputDirPrefix=""):
        """
        2012.10.5
            add a GATK SelectVariants job to update AC/AF of the final VCF file
            add argument refFastaFList, GenomeAnalysisTKJar
            
        2012.5.9
        """
        if GenomeAnalysisTKJar is None:
            GenomeAnalysisTKJar = self.GenomeAnalysisTKJar
        if refFastaFList is None:
            refFastaFList = self.refFastaFList
        
        sys.stderr.write("Adding vcf-subset jobs for %s vcf files ... "%(len(inputData.jobDataLs)))
        no_of_jobs= 0
        
        
        topOutputDir = "%sVCFSubset"%(outputDirPrefix)
        topOutputDirJob = self.addMkDirJob(outputDir=topOutputDir)
        no_of_jobs += 1
        
        returnData = PassingData()
        returnData.jobDataLs = []
        for jobData in inputData.jobDataLs:
            inputF = jobData.vcfFile
            chr = self.getChrFromFname(inputF.name)
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
            inputFBaseName = os.path.basename(inputF.name)
            commonPrefix = inputFBaseName.split('.')[0]
            outputVCF = File(os.path.join(topOutputDir, '%s.subset.vcf'%(commonPrefix)))
            vcfSubsetJob = self.addVCFSubsetJob(executable=self.vcfSubset,
                vcfSubsetPath=self.vcfSubsetPath,
                sampleIDFile=sampleIDFile,\
                inputVCF=inputF, outputF=outputVCF, \
                parentJobLs=[topOutputDirJob]+jobData.jobLs,
                transferOutput=False, job_max_memory=200,\
                extraArguments=None, extraDependentInputLs=None)
            
            #2012.10.5
            #selectVariants would generate AC, AF so that TrioCaller could read it.
            #samtools uses 'AC1' instead of AC, 'AF1' instead of AF.
            VCF4OutputF = File(os.path.join(topOutputDir, '%s.niceformat.vcf'%commonPrefix))
            vcfConvertJob = self.addSelectVariantsJob(
                SelectVariantsJava=self.SelectVariantsJava, \
                inputF=vcfSubsetJob.output, outputF=VCF4OutputF, \
                refFastaFList=refFastaFList, parentJobLs=[vcfSubsetJob], \
                extraDependentInputLs=[], transferOutput=False, \
                extraArguments=None, job_max_memory=2000, interval=chr)
            
            VCFGzipOutputF = File("%s.gz"%VCF4OutputF.name)
            VCFGzipOutput_tbi_F = File("%s.gz.tbi"%VCF4OutputF.name)
            bgzip_tabix_VCF_job = self.addBGZIP_tabix_Job(
                bgzip_tabix=self.bgzip_tabix, \
                parentJobLs=[vcfConvertJob], inputF=vcfConvertJob.output,
                outputF=VCFGzipOutputF, \
                transferOutput=transferOutput)
            
            returnData.jobDataLs.append(
                PassingData(jobLs=[bgzip_tabix_VCF_job], vcfFile=VCFGzipOutputF, \
                    tbi_F=VCFGzipOutput_tbi_F, \
                    fileLs=[VCFGzipOutputF, VCFGzipOutput_tbi_F]))
            
        sys.stderr.write("%s jobs.\n"%(self.no_of_jobs))
        return returnData
    

    def addSubsetAndVCF2PlinkJobs(self, inputData=None, db_main=None,
        minMAC=None, minMAF=None,\
        maxSNPMissingRate=None, sampleIDFile=None, transferOutput=True,\
        maxContigID=None, outputDirPrefix=""):
        """
        2012.5.9
        """
        vcfSubsetJobData = self.addVCFSubsetJobs(inputData=inputData, db_main=db_main,
            sampleIDFile=sampleIDFile, \
            transferOutput=True, maxContigID=maxContigID, outputDirPrefix="")
        vcf2plinkJobData = self.addVCF2PlinkJobs(inputData=vcfSubsetJobData,
            db_main=db_main, \
            minMAC=minMAC, minMAF=minMAF, maxSNPMissingRate=maxSNPMissingRate,
            transferOutput=transferOutput,\
            outputPedigreeAsTFAM=False, outputPedigreeAsTFAMInputJobData=None,
            maxContigID=maxContigID, outputDirPrefix="")
    
    def addAlignmentReadGroup2UCLAIDJobs(self, inputData=None, db_main=None,
        transferOutput=True,\
        maxContigID=None, outputDirPrefix=""):
        """
        2012.5.9
        """
        sys.stderr.write("Adding alignment read-group -> UCLAID jobs for %s vcf files ... "%(len(inputData.jobDataLs)))
        no_of_jobs= 0
        
        
        topOutputDir = "%sSampleInUCLAID"%(outputDirPrefix)
        topOutputDirJob = self.addMkDirJob(outputDir=topOutputDir)
        no_of_jobs += 1
        
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
            inputFBaseName = os.path.basename(inputF.name)
            commonPrefix = inputFBaseName.split('.')[0]
            outputVCF = File(os.path.join(topOutputDir, '%s.UCLAID.vcf'%(commonPrefix)))
            abstractMapperJob = self.addAbstractMapperLikeJob(
                executable=self.ConvertAlignmentReadGroup2UCLAIDInVCF, \
                inputVCF=inputF, outputF=outputVCF, \
                parentJobLs=[topOutputDirJob]+jobData.jobLs,
                transferOutput=False, job_max_memory=1000,\
                extraArguments=None, extraDependentInputLs=[])
            
            VCFGzipOutputF = File("%s.gz"%outputVCF.name)
            VCFGzipOutput_tbi_F = File("%s.gz.tbi"%outputVCF.name)
            bgzip_tabix_VCF_job = self.addBGZIP_tabix_Job(bgzip_tabix=self.bgzip_tabix, \
                parentJobLs=[abstractMapperJob], inputF=abstractMapperJob.output,
                outputF=VCFGzipOutputF, \
                transferOutput=transferOutput)
            
            returnData.jobDataLs.append(PassingData(jobLs=[bgzip_tabix_VCF_job],
                vcfFile=VCFGzipOutputF, \
                tbi_F=VCFGzipOutput_tbi_F, \
                fileLs=[VCFGzipOutputF, VCFGzipOutput_tbi_F]))
            
        sys.stderr.write("%s jobs.\n"%(self.no_of_jobs))
        return returnData
    
    def addSplitNamVCFJobs(self, inputData=None, db_main=None, transferOutput=True,\
        maxContigID=None, outputDirPrefix=""):
        """
        2012.5.11
            not functional. don't know what to do the fact that SplitNamVCFIntoMultipleSingleChrVCF outputs into a folder
                multiple VCF files (one per chromosome)
        """
        sys.stderr.write("Adding split Nam VCF-file jobs for %s vcf files ... "%(len(inputData.jobDataLs)))
        no_of_jobs= 0
        
        
        topOutputDir = "%sSampleInUCLAID"%(outputDirPrefix)
        topOutputDirJob = self.addMkDirJob(outputDir=topOutputDir)
        no_of_jobs += 1
        
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
            inputFBaseName = os.path.basename(inputF.name)
            commonPrefix = inputFBaseName.split('.')[0]
            outputVCF = File(os.path.join(topOutputDir, '%s.vcf'%(commonPrefix)))
            abstractMapperJob = self.addAbstractMapperLikeJob(executable=self.SplitNamVCFIntoMultipleSingleChrVCF, \
                    inputVCF=inputF, outputF=outputVCF, \
                    parentJobLs=[topOutputDirJob]+jobData.jobLs, transferOutput=False, job_max_memory=200,\
                    extraArguments=None, extraDependentInputLs=[])
            
            VCFGzipOutputF = File("%s.gz"%outputVCF.name)
            VCFGzipOutput_tbi_F = File("%s.gz.tbi"%outputVCF.name)
            bgzip_tabix_VCF_job = self.addBGZIP_tabix_Job(bgzip_tabix=self.bgzip_tabix, \
                    parentJobLs=[abstractMapperJob], inputF=abstractMapperJob.output,
                    outputF=VCFGzipOutputF, \
                    transferOutput=transferOutput)
            
            returnData.jobDataLs.append(PassingData(jobLs=[bgzip_tabix_VCF_job],
                vcfFile=VCFGzipOutputF, \
                tbi_F=VCFGzipOutput_tbi_F, \
                fileLs=[VCFGzipOutputF, VCFGzipOutput_tbi_F]))
            
            no_of_jobs += 2
        sys.stderr.write("%s jobs.\n"%(self.no_of_jobs))
        return returnData
    
    def addMergeVCFReplicateHaplotypesJobs(self, inputData=None, db_main=None,
        transferOutput=True,\
        maxContigID=None, outputDirPrefix="",replicateIndividualTag='copy',
        refFastaFList=None ):
        """
        2012.7.25
            input vcf is output of TrioCaller with replicates.
            this workflow outputs extra debug statistics
                1. replicate haplotype distance to the consensus haplotype
                2. majority support for the consensus haplotype
        """
        sys.stderr.write("Adding MergeVCFReplicateHaplotype jobs for %s vcf files ... "%(len(inputData.jobDataLs)))
        no_of_jobs= 0
        
        
        topOutputDir = "%sMergeVCFReplicateHaplotypeStat"%(outputDirPrefix)
        topOutputDirJob = self.addMkDirJob(outputDir=topOutputDir)
        no_of_jobs += 1
        
        
        haplotypeDistanceMergeFile = File(os.path.join(topOutputDir, 'haplotypeDistanceMerge.tsv'))
        haplotypeDistanceMergeJob = self.addStatMergeJob(self.mergeSameHeaderTablesIntoOne, \
            outputF=haplotypeDistanceMergeFile, transferOutput=False,
            parentJobLs=[topOutputDirJob])
        majoritySupportMergeFile = File(os.path.join(topOutputDir, 'majoritySupportMerge.tsv'))
        majoritySupportMergeJob = self.addStatMergeJob(self.mergeSameHeaderTablesIntoOne, \
            outputF=majoritySupportMergeFile, transferOutput=False,
            parentJobLs=[topOutputDirJob])
        no_of_jobs += 2
        
        returnData = PassingData()
        returnData.jobDataLs = []
        for jobData in inputData.jobDataLs:
            inputF = jobData.vcfFile
            inputFBaseName = os.path.basename(inputF.name)
            commonPrefix = inputFBaseName.split('.')[0]
            outputVCF = File(os.path.join(topOutputDir, '%s.vcf'%(commonPrefix)))
            debugHaplotypeDistanceFile = File(os.path.join(topOutputDir, '%s.haplotypeDistance.tsv'%(commonPrefix)))
            debugMajoritySupportFile = File(os.path.join(topOutputDir, '%s.majoritySupport.tsv'%(commonPrefix)))
            fileSize = utils.getFileOrFolderSize(pegaflow.getAbsPathOutOfFile(inputF))
            memoryRequest = 45000
            memoryRequest = min(42000, max(4000, int(38000*(fileSize/950452059.0))) )
            #extrapolates (33,000Mb memory for a ungzipped VCF file with size=950,452,059)
            #upper bound is 42g. lower bound is 4g.
            #mergeReplicateOutputF = File(os.path.join(trioCallerOutputDirJob.folder, '%s.noReplicate.vcf'%vcfBaseFname))
            #noOfAlignments= len(alignmentDataLs)
            #entireLength = stopPos - startPos + 1	#could be very small for shorter reference contigs
            #memoryRequest = min(42000, max(4000, int(20000*(noOfAlignments/323.0)*(entireLength/2600000.0))) )
                #extrapolates (20000Mb memory for a 323-sample + 2.6Mbase reference length/26K loci)
                #upper bound is 42g. lower bound is 4g.
            mergeVCFReplicateColumnsJob = self.addMergeVCFReplicateGenotypeColumnsJob(
                executable=self.MergeVCFReplicateHaplotypesJava,\
                GenomeAnalysisTKJar=self.GenomeAnalysisTKJar, \
                inputF=inputF, outputF=outputVCF,
                replicateIndividualTag=replicateIndividualTag,
                refFastaFList=refFastaFList,
                debugHaplotypeDistanceFile=debugHaplotypeDistanceFile,
                debugMajoritySupportFile=debugMajoritySupportFile,
                parentJobLs=[topOutputDirJob]+jobData.jobLs,
                extraDependentInputLs=[], transferOutput=False,
                extraArguments=None, job_max_memory=memoryRequest)
            
            #add output to some reduce job
            self.addInputToMergeJob(haplotypeDistanceMergeJob,
                inputF=mergeVCFReplicateColumnsJob.outputLs[1],
                parentJobLs=[mergeVCFReplicateColumnsJob])
            self.addInputToMergeJob(majoritySupportMergeJob,
                inputF=mergeVCFReplicateColumnsJob.outputLs[2],
                parentJobLs=[mergeVCFReplicateColumnsJob])
            no_of_jobs += 1
        sys.stderr.write("%s jobs. Done.\n"%(no_of_jobs))
        
        returnData.jobDataLs.append(PassingData(jobLs=[haplotypeDistanceMergeJob],
            file=haplotypeDistanceMergeFile,
            fileLs=[haplotypeDistanceMergeFile]))
        returnData.jobDataLs.append(PassingData(jobLs=[majoritySupportMergeJob],
            file=majoritySupportMergeFile,
            fileLs=[majoritySupportMergeFile]))
        #2012.7.21 gzip the final output
        newReturnData = self.addGzipSubWorkflow(inputData=returnData, 
            transferOutput=transferOutput,\
            outputDirPrefix="")
        return newReturnData
    
    def generateVCFSampleIDFilenameFromIndividualUCLAIDFname(self,
        db_main=None, individualUCLAIDFname=None,
        vcfSampleIDFname=None, oneSampleVCFFname=None):
        """
            
        """
        sys.stderr.write("Generating vcfSampleIDFname %s from individualUCLAIDFname %s ..."%(\
            vcfSampleIDFname, individualUCLAIDFname))
        
        #first get the set of monkeys to keep from the file
        reader = csv.reader(open(individualUCLAIDFname), 
            delimiter=figureOutDelimiter(individualUCLAIDFname))
        header = next(reader)
        colName2Index = getColName2IndexFromHeader(header)
        UCLAID_col_index = colName2Index.get('UCLAID')
        individualUCLAIDSet = set()
        for row in reader:
            individualUCLAID=row[UCLAID_col_index].strip()
            individualUCLAIDSet.add(individualUCLAID)
        sys.stderr.write(" %s uclaIDs. "%(len(individualUCLAIDSet)))
        del reader
        
        #second, read a sample VCF file and output the samples that have 
        # been in the given set
        writer = csv.writer(open(vcfSampleIDFname, 'w'), delimiter='\t')
        vcfFile = VCFFile(inputFname=oneSampleVCFFname, minDepth=0)
        no_of_samples = 0
        for sample_id in vcfFile.sample_id_ls:
            individual_code = db_main.parseAlignmentReadGroup(sample_id).individual_code
            if individual_code in individualUCLAIDSet:
                no_of_samples += 1
                writer.writerow([sample_id])
        del writer, vcfFile
        sys.stderr.write("%s vcf samples selected.\n"%(no_of_samples))
    
    def addCombineVCFIntoOneJobs(self, inputData=None, data_dir=None,
        maxContigID=None, outputDirPrefix="", genotypeMethodShortName=None,
        needSSHDBTunnel=False,
        transferOutput=True):
        """
        2012.8.30
        """
        sys.stderr.write("Adding jobs to combine all VCF files into one VCF for %s vcf files ... "%(len(inputData.jobDataLs)))
        no_of_jobs= 0
        
        topOutputDir = "%sVCFIntoOne"%(outputDirPrefix)
        topOutputDirJob = self.addMkDirJob(outputDir=topOutputDir)
        no_of_jobs += 1
        
        
        
        #2011-9-22 union of all samtools intervals for one contig
        if genotypeMethodShortName:
            #if this is available, output goes to db-affiliated storage.
            transferUnionOutput = False
        else:
            transferUnionOutput = transferOutput
        unionOutputFname = os.path.join(topOutputDir, 
            '%sVCFIntoOne.vcf.gz'%(len(inputData.jobDataLs)))
        unionOutputF = File(unionOutputFname)
        unionJob = self.addVCFConcatJob(
            concatExecutable=self.concatSamtools,
            parentDirJob=topOutputDirJob, \
            outputF=unionOutputF, transferOutput=transferUnionOutput,
            vcf_job_max_memory=2000)
        
        no_of_jobs += 1
        
        returnData = PassingData()
        returnData.jobDataLs = []
        returnData.jobDataLs.append(PassingData(jobLs=[unionJob],
            file=unionJob.output,
            fileLs=unionJob.outputLs))
        
        for i in range(len(inputData.jobDataLs)):
            jobData = inputData.jobDataLs[i]
            inputF = jobData.vcfFile
            inputFBaseName = os.path.basename(inputF.name)
            chr_id = self.getChrFromFname(inputFBaseName)
            if maxContigID:
                contig_id = self.getContigIDFromFname(inputFBaseName)
                try:
                    contig_id = int(contig_id)
                    if contig_id>maxContigID:
                        #skip the small contigs
                        continue
                except:
                    sys.stderr.write('Except type: %s\n'%repr(sys.exc_info()))
                    import traceback
                    traceback.print_exc()
            commonPrefix = inputFBaseName.split('.')[0]
            outputFnamePrefix = os.path.join(topOutputDir, '%s'%(commonPrefix))
            if i ==0:	#need at least one tfam file. 
                transferOneContigPlinkOutput = True
            else:
                transferOneContigPlinkOutput = False
            i += 1
            
            #add this output to a union job
            self.addInputToMergeJob(unionJob, inputF=inputF,
                parentJobLs=jobData.jobLs, extraDependentInputLs=[jobData.tbi_F])
            
        if genotypeMethodShortName:
            logFile = File(os.path.join(topOutputDir, 'addVCF2DB.log'))
            addVCFJob = self.addAddVCFFile2DBJob(executable=self.AddVCFFile2DB,
                inputFile=unionJob.output, \
                genotypeMethodShortName=genotypeMethodShortName,
                logFile=logFile, format="VCF", data_dir=data_dir,
                checkEmptyVCFByReading=True, commit=True,
                parentJobLs=[unionJob], extraDependentInputLs=[unionJob.tbi_F],
                transferOutput=transferOutput,
                extraArguments=None, job_max_memory=1000, sshDBTunnel=needSSHDBTunnel)
            no_of_jobs += 1
        sys.stderr.write("%s jobs. Done.\n"%(no_of_jobs))
        ##2012.8.9 gzip workflow is not needed anymore as binary bed is used instead.
        ##2012.7.21 gzip the final output
        #newReturnData = self.addGzipSubWorkflow(inputData=returnData,
        #   transferOutput=transferOutput,
        #				outputDirPrefix="")
        return returnData
        
    
    def registerExecutables(self):
        """
        2011-11-28
        """
        ParentClass.registerExecutables(self)
        self.registerOneExecutable(path=os.path.join(self.pymodulePath,
            "mapper/ConvertAlignmentReadGroup2UCLAIDInVCF.py"),
            name='ConvertAlignmentReadGroup2UCLAIDInVCF', clusterSizeMultiplier=1)
        self.registerOneExecutable(path=os.path.join(self.pymodulePath,
            "mapper/SplitNamVCFIntoMultipleSingleChrVCF.py"),
            name='SplitNamVCFIntoMultipleSingleChrVCF', clusterSizeMultiplier=1)
        self.registerOneExecutable(path=os.path.join(self.pymodulePath,
            "mapper/filter/ModifyTPED.py"),
            name='ModifyTPED', clusterSizeMultiplier=1)
        
        self.registerOneExecutable(path=os.path.join(self.pymodulePath,
            "mapper/extractor/vcfSubset.sh"),\
            name='vcfSubset', clusterSizeMultiplier=1)
        if self.vcfSubsetPath:
            self.vcfSubset.vcfSubsetPath = self.vcfSubsetPath
            #vcfSubsetPath is first argument to vcfSubset
            self.vcfSubsetExecutableFile = self.registerOneExecutableAsFile(
                path=self.vcfSubsetPath)
        self.registerOneExecutable(path=self.plinkPath,
            name='plink', clusterSizeMultiplier=1)
        self.registerOneExecutable(path=self.plinkPath,
            name='plinkNoClustering', clusterSizeMultiplier=0)
        self.registerOneExecutable(path=self.plinkPath,
            name='plinkMerge', clusterSizeMultiplier=0)
        self.registerOneExecutable(path=self.plinkPath,
            name='plinkIBD', clusterSizeMultiplier=0)
        self.registerOneExecutable(path=self.plinkPath,
            name='plinkConvert', clusterSizeMultiplier=1)
        self.registerOneExecutable(path=self.plinkPath,
            name='plinkLDPrune', clusterSizeMultiplier=1)
        self.registerOneExecutable(path=self.plinkPath, \
            name='plinkExtract', clusterSizeMultiplier=1)
        self.registerOneExecutable(path=os.path.join(self.pymodulePath,
            'mapper/modifier/SplitPlinkLMendelFileSNPIDIntoChrPosition.py'),
            name='SplitPlinkLMendelFileSNPIDIntoChrPosition', clusterSizeMultiplier=1)
        self.registerOneExecutable(path=os.path.join(self.pymodulePath,
            'pedigree/CalculateMendelErrorRateGivenPlinkOutput.py'),
            name='CalculateMendelErrorRateGivenPlinkOutput', clusterSizeMultiplier=1)

    
    def run(self):
        """
        2011-9-28
        """
        self.setup_run()
        inputData = self.registerAllInputFiles(self.inputDir, 
            input_site_handler=self.input_site_handler, \
            checkEmptyVCFByReading=self.checkEmptyVCFByReading,\
            pegasusFolderName=self.pegasusFolderName,\
            maxContigID=self.maxContigID, \
            minContigID=self.minContigID)
        if len(inputData.jobDataLs)<=0:
            sys.stderr.write("No VCF files in this folder , %s.\n"%self.inputDir)
            sys.exit(0)
        
        if self.individualUCLAIDFname and os.path.isfile(self.individualUCLAIDFname):
            self.generateVCFSampleIDFilenameFromIndividualUCLAIDFname(
                db_main=self.db_main,
                individualUCLAIDFname=self.individualUCLAIDFname,
                vcfSampleIDFname=self.vcfSampleIDFname,
                oneSampleVCFFname=inputData.jobDataLs[0].vcfFile.abspath)
            sampleIDFile = self.registerOneInputFile(self.vcfSampleIDFname)
        elif self.vcfSampleIDFname and os.path.isfile(self.vcfSampleIDFname):
            sampleIDFile = self.registerOneInputFile(self.vcfSampleIDFname)
        else:
            sampleIDFile = None
        
        if self.run_type==1:
            if sampleIDFile is None:
                sys.stderr.write("sampleIDFile is None.\n")
                sys.exit(0)
            self.addVCFSubsetJobs(
                inputData=inputData, db_main=self.db_main, sampleIDFile=sampleIDFile,
                transferOutput=True,\
                refFastaFList=self.refFastaFList,
                GenomeAnalysisTKJar=self.GenomeAnalysisTKJar,\
                maxContigID=self.maxContigID, outputDirPrefix="")
        elif self.run_type==2:
            #2012.8.10 test  ModifyTPEDRunType=3, chr_id2cumu_chr_start=None
            self.addVCF2PlinkJobs(
                inputData=inputData, db_main=self.db_main,
                minMAC=self.minMAC, minMAF=self.minMAF,\
                maxSNPMissingRate=self.maxSNPMissingRate, transferOutput=True,\
                outputPedigreeAsTFAM=False, outputPedigreeAsTFAMInputJobData=None, \
                maxContigID=self.maxContigID, outputDirPrefix="")
        elif self.run_type==3:
            if sampleIDFile is None:
                sys.stderr.write("sampleIDFile is None.\n")
                sys.exit(0)
            self.addSubsetAndVCF2PlinkJobs(
                inputData=inputData, db_main=self.db_main, minMAC=self.minMAC, \
                minMAF=self.minMAF,\
                maxSNPMissingRate=self.maxSNPMissingRate, sampleIDFile=sampleIDFile, transferOutput=True,\
                maxContigID=self.maxContigID, outputDirPrefix="")
        elif self.run_type==4:
            self.addAlignmentReadGroup2UCLAIDJobs(
                inputData=inputData, 
                db_main=self.db_main, transferOutput=True,\
                maxContigID=self.maxContigID, outputDirPrefix="")
        elif self.run_type==5:
            refSequence = SunsetDB.IndividualSequence.get(self.ref_ind_seq_id)
            refFastaFname = os.path.join(self.data_dir, refSequence.path)
            registerReferenceData = self.registerRefFastaFile(refFastaFname,
                registerAffiliateFiles=True,
                checkAffiliateFileExistence=True)
            self.addMergeVCFReplicateHaplotypesJobs(inputData=inputData, 
                db_main=self.db_main, transferOutput=True,\
                maxContigID=self.maxContigID, outputDirPrefix="", 
                replicateIndividualTag='copy', \
                refFastaFList=registerReferenceData.refFastaFList )
        elif self.run_type==6:
            self.addVCF2YuFormatJobs(inputData=inputData, transferOutput=True,\
                maxContigID=self.maxContigID, outputDirPrefix="", \
                returnMode=1)
        elif self.run_type==7:
            #first convert every sample ID from alignment.read_group to simple ucla ID 
            inputData2 = self.addAlignmentReadGroup2UCLAIDJobs(
                inputData=inputData, 
                db_main=self.db_main, \
                transferOutput=False,\
                maxContigID=self.maxContigID, outputDirPrefix="")
            self.addVCF2YuFormatJobs(inputData=inputData2, transferOutput=True,\
                maxContigID=self.maxContigID, outputDirPrefix="", \
                returnMode=1)
        elif self.run_type==8:
            inputData2 = self.addAlignmentReadGroup2UCLAIDJobs(inputData=inputData, 
                db_main=self.db_main, transferOutput=False,\
                maxContigID=self.maxContigID, outputDirPrefix="")
            self.addVCF2PlinkJobs(inputData=inputData2, db_main=self.db_main, \
                minMAC=self.minMAC, minMAF=self.minMAF,\
                maxSNPMissingRate=self.maxSNPMissingRate, transferOutput=True,\
                outputPedigreeAsTFAM=True,
                outputPedigreeAsTFAMInputJobData=inputData.jobDataLs[0],
                maxContigID=self.maxContigID, outputDirPrefix="", returnMode=3)
        elif self.run_type==9:
            self.addCombineVCFIntoOneJobs(inputData=inputData, data_dir=self.data_dir,
                maxContigID=self.maxContigID, outputDirPrefix="", 
                genotypeMethodShortName=self.genotypeMethodShortName, \
                needSSHDBTunnel=self.needSSHDBTunnel, \
                transferOutput=True)
        else:
            sys.stderr.write("run_type %s not supported.\n"%(self.run_type))
            sys.exit(0)
        
        self.end_run()
        


if __name__ == '__main__':
    from argparse import ArgumentParser
    ap = ArgumentParser()
    ap.add_argument("--plinkPath", type=str, default='bin/plink',
        help="Path to plink, "
        "http://pngu.mgh.harvard.edu/~purcell/plink/index.shtml "
        "If relative path, home folder is inserted in the front.")
    
