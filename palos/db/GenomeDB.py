#!/usr/bin/env python3
"""
Examples:
    # 2010-12-15 setup genome schema in vervetdb.
    GenomeDB.py -v postgresql -k genome -d vervetdb -u yh
    
    #setup database in mysql
    #  (if it's already setup, output the pickled geneAnnotation data structure)
    GenomeDB.py -v mysql -z papaya -d genome -k "" -u yh -t 3702
        -o /tmp/geneAnnotationTax3702.pickle

Description:
    This is the genome database ORM.
"""
import sys, os
import logging
from sqlalchemy.engine.url import URL
from sqlalchemy import Unicode, DateTime, String, BigInteger, Integer
from sqlalchemy import UnicodeText, Text, Boolean, Float, Binary, Enum, Table
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import Column
from sqlalchemy import ForeignKey
from sqlalchemy.orm import relationship
from sqlalchemy import UniqueConstraint
from datetime import datetime
from sqlalchemy.schema import MetaData
from sqlalchemy.orm import scoped_session, sessionmaker
from sqlalchemy import UniqueConstraint, create_engine
from sqlalchemy import and_, or_, not_
from sqlalchemy import desc
from palos.utils import PassingData	
from palos.polymorphism.CNV import CNVCompare, CNVSegmentBinarySearchTreeKey
from palos.algorithm.RBTree import RBDict
from palos.Genome import GeneModel
#2010-9-21 although "from Genome import GeneModel" works,
# it causes problem in pickle.load() because Genome is not directly visible outside.
from palos.db import Database, TableClass, get_sequence_segment

Base = declarative_base()
#Set it staticaly because SunsetDB is undefined at this point 
# and it has to be defined after this.
# 20200424 now that schema is set at the engine level.
#  Try not to use _schemaname_.
_schemaname_ = "genome_hg37"

class SequenceType(Base, TableClass):
    """
    2013.3.14 added description
    2008-07-27
        a table storing meta information to be referenced by other tables
    """
    __tablename__ = 'sequence_type'
    __table_args__ = {'schema':_schemaname_}
    id = Column(Integer, primary_key=True)
    short_name = Column(String(223), unique=True)
    description = Column(Text)
    created_by = Column(String(256))
    updated_by = Column(String(256))
    date_created = Column(DateTime, default=datetime.now)
    date_updated = Column(DateTime)

class RawSequence(Base, TableClass):
    """
    A table to store chunks of sequences of entries from AnnotAssembly.
    annot_assembly_gi is not a foreign key anymore.
    annot_assembly_id replaces it.
    If the associated entry in AnnotAssembly is deleted,
        all raw sequences will be deleted as well.
    """
    __tablename__ = 'raw_sequence'
    __table_args__ = {'schema':_schemaname_}
    id = Column(Integer, primary_key=True)
    annot_assembly_id = Column(Integer,
        ForeignKey(f"{_schemaname_}.annot_assembly.id"))
    annot_assembly = relationship('AnnotAssembly')
    annot_assembly_gi = Column(Integer)
    start = Column(Integer)
    stop = Column(Integer)
    sequence = Column(String(10000))	#each fragment is 10kb
    UniqueConstraint('annot_assembly_id', 'start', 'stop', \
        name='raw_sequence_uniq_cnst')

class AnnotAssembly(Base, TableClass):
    """
    2011-8-24
        gi is no longer the primary key.
        a new id is added as primary key.
    2010-12-17
        Column raw_sequence_start_id is no longer a foreign key.
        To avoid circular dependency with RawSequence.
    2008-07-27
        table to store meta info of chromosome sequences
    """
    __tablename__ = 'annot_assembly'
    __table_args__ = {'schema':_schemaname_}
    id = Column(Integer, primary_key=True)
    gi = Column(Integer)
    acc_ver = Column(String(32))
    accession = Column(String(32))
    version = Column(Integer)
    tax_id = Column(Integer)
    chromosome = Column(String(128))
    start = Column(Integer)
    stop = Column(Integer)
    genome_order = Column(Integer)
    orientation = Column(String(1))
    sequence = Column(String(10000))
    raw_sequence_start_id = Column(Integer)
    original_path = Column(Text)
    sequence_type_id = Column(Integer,
        ForeignKey(f"{_schemaname_}.sequence_type.id"))
    sequence_type = relationship('SequenceType')
    chromosome_type_id = Column(Integer, 
        ForeignKey(f"{_schemaname_}.chromosome_type.id"))
    chromosome_type = relationship('ChromosomeType')
    genome_annotation_list = relationship('GenomeAnnotation',
        back_populates='annot_assembly',cascade='all,delete')
    #non-zero means outdated. to allow multiple outdated alignments
    outdated_index = Column(Integer, default=0)
    comment = Column(Text)
    created_by = Column(String(256))
    updated_by = Column(String(256))
    date_created = Column(DateTime, default=datetime.now)
    date_updated = Column(DateTime)
    UniqueConstraint(
        'accession', 'version','tax_id','chromosome', 'start', 'stop',
        'orientation', 'sequence_type_id', 'chromosome_type_id',
        'outdated_index')


class ChromosomeType(Base, TableClass):
    """
    store autosome, X, Y, mitochondrial
    """
    __tablename__ = 'chromosome_type'
    __table_args__ = {'schema':_schemaname_}
    id = Column(Integer, primary_key=True)    
    short_name = Column(String(223), unique=True)
    comment = Column(Text)
    created_by = Column(String(256))
    updated_by = Column(String(256))
    date_created = Column(DateTime, default=datetime.now)
    date_updated = Column(DateTime)

class GenomeAnnotation(Base, TableClass):
    """
    2013.2.17
    """
    __tablename__ = 'genome_annotation'
    __table_args__ = {'schema':_schemaname_}
    id = Column(Integer, primary_key=True)
    strand = Column(String(4))
    start = Column(Integer)
    stop = Column(Integer)
    description = Column(Text)
    comment = Column(Text)
    annot_assembly_id = Column(Integer,
        ForeignKey(f"{_schemaname_}.annot_assembly.id"))
    annot_assembly = relationship('AnnotAssembly')
    genome_annotation_type_id = Column(Integer, \
        ForeignKey(f"{_schemaname_}.genome_annotation_type.id"))
    genome_annotation_type = relationship('GenomeAnnotationType')
    created_by = Column(String(256))
    updated_by = Column(String(256))
    date_created = Column(DateTime, default=datetime.now)
    date_updated = Column(DateTime)
    UniqueConstraint('strand', 'start', 'stop', 'annot_assembly_id',\
        'genome_annotation_type_id')


class GenomeAnnotationType(Base, TableClass):
    """
    2013.2.17
    table to store centromere, telemere, euchromatin, heterochromatin
    """
    __tablename__ = 'genome_annotation_type'
    __table_args__ = {'schema':_schemaname_}
    id = Column(Integer, primary_key=True)
    short_name = Column(String(223), unique=True)
    description = Column(Text)
    comment = Column(Text)
    created_by = Column(String(256))
    updated_by = Column(String(256))
    date_created = Column(DateTime, default=datetime.now)
    date_updated = Column(DateTime)

class EntrezgeneType(Base, TableClass):
    """
    2008-07-28
        store the entrez gene types
    """
    __tablename__ = 'entrezgene_type'
    __table_args__ = {'schema':_schemaname_}
    id = Column(Integer, primary_key=True)
    type = Column(String(223), unique=True)
    created_by = Column(String(256))
    updated_by = Column(String(256))
    date_created = Column(DateTime, default=datetime.now)
    date_updated = Column(DateTime)

class GeneCommentaryType(Base, TableClass):
    """
    2008-07-28
    """
    __tablename__ = 'gene_commentary_type'
    __table_args__ = {'schema':_schemaname_}
    id = Column(Integer, primary_key=True)
    type = Column(String(223), unique=True)
    created_by = Column(String(256))
    updated_by = Column(String(256))
    date_created = Column(DateTime, default=datetime.now)
    date_updated = Column(DateTime)

class GeneCommentary(Base, TableClass):
    """
    2010-12-15
        gene linked to table Gene
    2008-07-28
        store different mRNAs/Peptides from the same gene or mRNA
    """
    __tablename__ = 'gene_commentary'
    __table_args__ = {'schema':_schemaname_}
    id = Column(Integer, primary_key=True)
    accession = Column(String(32))
    version = Column(Integer)
    gi = Column(Integer)
    gene_id = Column(Integer, ForeignKey(f"{_schemaname_}.gene.id"))
    gene = relationship('Gene')
    gene_commentary_id = Column(Integer,
        ForeignKey(f"{_schemaname_}.gene_commentary.id"))
    gene_commentary = relationship('GeneCommentary')
    #gene_commentary_ls = relationship('GeneCommentary',
    #   back_populates='gene_commentary', cascade='all,delete')
    start = Column(Integer)
    stop = Column(Integer)
    gene_commentary_type_id = Column(Integer, 
        ForeignKey(f"{_schemaname_}.gene_commentary_type.id"))
    gene_commentary_type = relationship('GeneCommentaryType')
    label = Column(Text)
    text = Column(Text)
    comment = Column(Text)
    gene_segment_ls = relationship('GeneSegment',
        back_populates='gene_commentary', cascade='all,delete')
    created_by = Column(String(256))
    updated_by = Column(String(256))
    date_created = Column(DateTime, default=datetime.now)
    date_updated = Column(DateTime)
    UniqueConstraint('gene_id', 'start', 'stop',
        'gene_commentary_type_id', 'gene_commentary_id')
    
    def getSequence(self, box_ls):
        """
2012.5.16
    it uses annot_assembly_id preferentially, if it is null,
        then uses genomic_gi.
2012.3.26
    get_sequence_segment()'s API has been changed.
    But this way of calling get_sequence_segment() is still outdated.
        valid only for banyan's mysql genome db
        because annot_assembly_id was not used.
        check GenomeDatabase.getSequenceSegment()
2009-01-03
        """
        seq = ''
        curs = __metadata__.bind
        #or self.table.bind or self.table.metadata.bind
        for box in box_ls:
            genomic_start, genomic_stop = box[:2]
            if self.gene.annot_assembly_id:
                annot_assembly_id = self.gene.annot_assembly_id
                gi = None
            elif self.gene.genomic_gi:
                gi = self.gene.genomic_gi
                annot_assembly_id = None
            else:
                gi = None
                annot_assembly_id = None
            schema = self._descriptor.table_options.get('schema')
            seq += get_sequence_segment(curs, gi=gi, start=genomic_start,
                stop=genomic_stop,
                annot_assembly_id=annot_assembly_id,
                annot_assembly_table='%s.%s'%(schema, AnnotAssembly.__tablename__),
                raw_sequence_table='%s.%s'%(schema, RawSequence.__tablename__))
        if self.gene.strand=='-1' and seq:	#need to reverse complement
            from Bio.Seq import Seq
            from Bio.Alphabet import IUPAC
            seq1 = Seq(seq, IUPAC.unambiguous_dna)
            seq2 = seq1.reverse_complement()
            seq = seq2.data
        return seq
    
    def getCDSsequence(self):
        """
        2009-01-03 implemented
        2008-09-22 implement later
        """
        if not getattr(self, 'protein_box_ls', None):
            self.construct_protein_box()
        
        self.cds_sequence = self.getSequence(self.protein_box_ls)
        #self.mrna_sequence = self.getSequence(self.mrna_box_ls)
    
    def construct_mrna_box(self):
        """
        2012.5.15
            GeneSegment.gene_commentary_type could be either exon or mRNA.
            add gene_segment.id into the mrna_box
            
        2010-8-18
            update to fetch exon/mRNA only because introns are now added in GeneSegment.
        2008-09-22
        """
        #either exon or mRNA type
        gene_commentary_type_id_ls = []
        gene_commentary_type = GeneCommentaryType.query.filter_by(type='exon').first()
        if gene_commentary_type:
            gene_commentary_type_id_ls.append(gene_commentary_type.id)
        gene_commentary_type = GeneCommentaryType.query.filter_by(type='mRNA').first()
        if gene_commentary_type:
            gene_commentary_type_id_ls.append(gene_commentary_type.id)
        
        gene_segment_ls = GeneSegment.query.\
            filter(GeneSegment.gene_commentary_type_id.in_(
                gene_commentary_type_id_ls)).\
            filter_by(gene_commentary_id=self.id)
        self.mrna_box_ls = []
        for gene_segment in gene_segment_ls:
            self.mrna_box_ls.append([gene_segment.start, gene_segment.stop,
                gene_segment.id])
        self.mrna_box_ls.sort()
        return self.mrna_box_ls
    
    def construct_protein_box(self):
        """
        2012.5.15
            GeneSegment.gene_commentary_type could be either CDS or peptide.
            add gene_segment.id into the protein_box
        2008-09-22
            find the corresponding protein gene_commentary if it's available.
            and construct a protein_box_ls
        """
        self.protein_box_ls = []
        if len(self.gene_commentary_ls)==1:
            gene_commentary = self.gene_commentary_ls[0]
        elif len(self.gene_commentary_ls)>1:
            logging.warn(f'More than 1 gene_commentary_ls for this '
                f'commentary id={self.id}, gene_id={self.gene_id}.')
            gene_commentary = self.gene_commentary_ls[0]
        else:
            gene_commentary = None
        if gene_commentary:
            #either exon or mRNA type
            gene_commentary_type_id_ls = []
            gene_commentary_type = GeneCommentaryType.query.filter_by(type='CDS').first()
            if gene_commentary_type:
                gene_commentary_type_id_ls.append(gene_commentary_type.id)
            gene_commentary_type = GeneCommentaryType.query.\
                filter_by(type='peptide').first()
            if gene_commentary_type:
                gene_commentary_type_id_ls.append(gene_commentary_type.id)
            
            gene_segment_ls = GeneSegment.query.\
                filter(GeneSegment.gene_commentary_type_id.in_(
                    gene_commentary_type_id_ls)).\
                filter_by(gene_commentary_id=gene_commentary.id)
            for gene_segment in gene_segment_ls:
                self.protein_box_ls.append([gene_segment.start, \
                    gene_segment.stop, gene_segment.id])
            self.protein_label = gene_commentary.label
            self.protein_comment = gene_commentary.comment
            self.protein_text = gene_commentary.text
        else:
            self.protein_label = None
            self.protein_comment = None
            self.protein_text = None
        self.protein_box_ls.sort()
        return self.protein_box_ls
    
    def constructAnnotatedBox(self):
        """
2010-8-18
    deal with the db (GeneSegment, etc.) populated by TAIRGeneXML2GenomeDB.py,
        not GeneASNXML2gene_mapping.py.
    GeneSegment store, UTR, CDS, exon, intron in their commentary_type already.
    
    #each entry is a tuple, (start, stop, box_type, is_translated, gene_segment.id, protein_box_index)
        """
        #each entry is a tuple, (start, stop, box_type, is_translated,
        #  gene_segment.id, protein_box_index)
        box_ls = []
        if len(self.gene_commentary_ls)==1:	#it's translated into protein.
            protein_commentary = self.gene_commentary_ls[0]
        elif len(self.gene_commentary_ls)>1:
            logging.warn(f'More than 1 gene_commentary_ls for this '
                f'commentary id={self.id}, gene_id={self.gene_id}.')
            protein_commentary = self.gene_commentary_ls[0]
        else:
            protein_commentary = None
        
        query = GeneSegment.query.filter_by(gene_commentary_id=self.id)
        if protein_commentary:
            #restrict the gene segment to intron only. 
            # exons will be covered by CDS and UTR.
            gene_commentary_type = GeneCommentaryType.query.\
                filter_by(type='intron').first()
            if gene_commentary_type:
                query = query.filter_by(
                    gene_commentary_type_id=gene_commentary_type.id)
        
        for gene_segment in query:
            box_ls.append([gene_segment.start, gene_segment.stop,\
                gene_segment.gene_commentary_type.type, 0,
                gene_segment.id, None])
        if protein_commentary:
            for gene_segment in protein_commentary.gene_segment_ls:
                if gene_segment.gene_commentary_type.type.find('UTR')!=-1:
                    is_translated = 0
                else:
                    is_translated = 1
                
                box_ls.append([gene_segment.start, gene_segment.stop,\
                    gene_segment.gene_commentary_type.type, \
                    is_translated, gene_segment.id,  None])
        box_ls.sort()
        return box_ls
        
    def construct_annotated_box(self):
        """
2012.5.15
    add gene_segment.id(or None) and cumulativeWithinCDSUTRAndIntronLen into box_ls
2008-10-01
    fix a bug that coordinates of a whole untranslated mrna block are
        replaced by that of a protein block
2008-10-01
    fix a bug that a whole untranslated mrna block got totally omitted.
2008-09-22
    combine mrna_box_ls and protein_box_ls to partition the whole gene into finer segments.
    box_ls = []	#each entry is a tuple,
        (start, stop, box_type, is_translated, protein_box_index,
            gene_segment.id, detailed_box_type, \
            cumulativeWithinCDSUTRAndIntronLen, exon_number)
    box_type = 'intron' or 'exon'. is_translated = 0 or 1. if it's translated,
        protein_box_index is index in protein_box_ls.
        """
        #each entry is a tuple, (start, stop, box_type, is_translated, protein_box_index)
        self.box_ls = []
        cumulativeWithinCDSUTRAndIntronLen=0
        #the cumulative length of all withinCDS UTRs and introns from the 5'\
        # till the designated CDS
        CDS_5_end_pos = -1	#5' end position of the whole CDS sequence
        if not hasattr(self, 'mrna_box_ls'):
            self.construct_mrna_box()
        
        if not hasattr(self, 'protein_box_ls'):
            self.construct_protein_box()
        
        no_of_mrna_boxes = len(self.mrna_box_ls)
        no_of_prot_boxes = len(self.protein_box_ls)
        if no_of_prot_boxes==0:
            # 2012.5.15 no protein, it's a non-coding gene
            default_detailed_box_type = 'exon'
        else:
            if self.gene.strand=='+1':
                default_detailed_box_type = '3UTR'
                #this is used after all proteins boxes
            else:
                default_detailed_box_type = '5UTR'
        j = 0	#index in protein_box_ls
        no_of_introns = 0
        exon_number = 0
        for i in range(no_of_mrna_boxes):
            mrna_start, mrna_stop, mrna_gene_segment_id = self.mrna_box_ls[i][:3]
            exon_number += 1
            if i>0:	#add intron if this is not the first exon
                intron_start = self.mrna_box_ls[i-1][1]+1	#stop of previous exon + 1
                intron_stop = mrna_start-1	#start of current exon + 1
                self.box_ls.append((intron_start, intron_stop, 'intron', 0, None, None, \
                    'intron', cumulativeWithinCDSUTRAndIntronLen, None))
                    #intron, no gene_segment
                no_of_introns += 1
                if CDS_5_end_pos!=-1:	#CDS has already begun
                    nonCDSLength = abs(intron_stop-intron_start)+1
                    cumulativeWithinCDSUTRAndIntronLen += nonCDSLength
            if j<no_of_prot_boxes:
                prot_start, prot_stop, prot_gene_segment_id = self.protein_box_ls[j][:3]
                if prot_start>=mrna_start and prot_stop<=mrna_stop:
                    if prot_start>mrna_start:	#one untranslated exon
                        if self.gene.strand=='+1':
                            detailed_box_type = '5UTR'
                        else:
                            detailed_box_type = '3UTR'
                        self.box_ls.append((mrna_start, prot_start-1, 'exon', 0, None, \
                            mrna_gene_segment_id, detailed_box_type, \
                            cumulativeWithinCDSUTRAndIntronLen, exon_number))
                            #UTR before CDS
                        if CDS_5_end_pos!=-1:	#CDS has already begun
                            nonCDSLength = abs(prot_start - mrna_start)
                            cumulativeWithinCDSUTRAndIntronLen += nonCDSLength
                    self.box_ls.append((prot_start, prot_stop, 'exon', 1, j, \
                        prot_gene_segment_id, 'CDS', 
                        cumulativeWithinCDSUTRAndIntronLen, exon_number))
                        #CDS
                    if CDS_5_end_pos==-1:	#first CDS
                        CDS_5_end_pos = prot_start
                    if prot_stop<mrna_stop:	#one more untranslated exon
                        if self.gene.strand=='+1':
                            detailed_box_type = '3UTR'
                        else:
                            detailed_box_type = '5UTR'
                        self.box_ls.append((prot_stop+1, mrna_stop, 'exon', 0, None, \
                            mrna_gene_segment_id, detailed_box_type, \
                            cumulativeWithinCDSUTRAndIntronLen, exon_number))
                        #UTR after CDS
                        if CDS_5_end_pos!=-1:	#CDS has already begun
                            nonCDSLength = abs(mrna_stop-prot_stop)
                            cumulativeWithinCDSUTRAndIntronLen += nonCDSLength
                    j += 1	#push protein box index up
                elif prot_stop<mrna_start:
                    #not supposed to happen
                    logging.error(f"protein box: [{prot_start}, {prot_stop}] of"
                        f" gene-id={self.gene_id} is ahead of mrna box "
                        f"[{mrna_start}, {mrna_stop}].")
                elif prot_start>mrna_stop:
                    if self.gene.strand=='+1':
                        detailed_box_type = '5UTR'
                    else:
                        detailed_box_type = '3UTR'
                    self.box_ls.append((mrna_start, mrna_stop, 'exon', 0, None, \
                        mrna_gene_segment_id, detailed_box_type, \
                        cumulativeWithinCDSUTRAndIntronLen, exon_number))
                    #2008-10-1 UTR
                    if CDS_5_end_pos!=-1:	#CDS has already begun
                            nonCDSLength = abs(mrna_stop-mrna_start+1)
                            cumulativeWithinCDSUTRAndIntronLen += nonCDSLength
                elif prot_start<=mrna_stop and prot_stop>mrna_stop:
                    #not supposed to happen
                    logging.error(f"protein box: [{prot_start}, {prot_stop}] of"
                        f" gene-id={self.gene_id} overlaps partially with "
                        f"mRNA box [{mrna_start}, {mrna_stop}].")
            else:	#passing all protein boxes
                self.box_ls.append((mrna_start, mrna_stop, 'exon', 0, None, \
                    mrna_gene_segment_id, default_detailed_box_type, \
                    cumulativeWithinCDSUTRAndIntronLen, exon_number))
                if CDS_5_end_pos!=-1:	#CDS has already begun
                    nonCDSLength = abs(mrna_stop-mrna_start+1)
                    cumulativeWithinCDSUTRAndIntronLen += nonCDSLength
                #UTR after all CDSs or exon
        self.CDS_5_end_pos = CDS_5_end_pos
        self.no_of_introns = no_of_introns
        return self.box_ls
    
    def outputSequenceInUpperLowerCase(self, outf):
        """
        2012.6.18
            get the box_ls from self.gene_segment_ls.
            If it's mRNA commentary, this function outputs pre-splicing mRNA.
            If it's protein commentary, this function outputs CDS sequence with introns.
        2012.6.11
            Segments present in GeneSegment are outputted in upper case.
            All inter-between segments are outputted in lower case.
        """
        """
        #2012.6.18 commented out
        import copy
        #use copy to avoid modifying self.mrna_box_ls and self.protein_box_ls
        box_ls = copy.deepcopy(self.construct_mrna_box())
        
        if not box_ls:	#this is a protein gene commentary.
            box_ls = copy.deepcopy(self.construct_protein_box())
        """
        box_ls = []
        for geneSegment in self.gene_segment_ls:
            box_ls.append([geneSegment.start, geneSegment.stop])
        box_ls.sort()
        
        outf.write(">GeneCommentary%s_fromGene%s\n"%(self.id, self.gene.id))
        if self.start<box_ls[0][0]:	#need to output the beginning in lower case.
            firstBoxStart = box_ls[0][0]
            sequence = self.getSequence([[self.start, firstBoxStart-1]])
            outf.write("%s"%(sequence.lower()))
        
        for i in range(len(box_ls)):
            box = box_ls[i]
            sequence = self.getSequence([[box[0], box[1]]])
            outf.write("%s"%(sequence.upper()))
            if i<len(box_ls)-1:
                #output the inter-box sequences in lower case
                next_box = box_ls[i+1]
                sequence = self.getSequence([[box[1]+1, next_box[0]-1]])
                outf.write("%s"%(sequence.lower()))
        if self.stop>box_ls[-1][1]:
            lastBoxStop = box_ls[-1][1]
            sequence = self.getSequence([[lastBoxStop+1, self.stop]])
            outf.write("%s"%(sequence.lower()))
        outf.write("\n")
    
    def getProteinSequence(self,):
        """
        2012.6.18
            assuming this is a mRNA gene commentary. won't work on its derivative protein commentary.
            get the protein_box_ls
            get all nucleotide sequence
            translate it into protein sequence
            return a fastq record.
        """
        self.getCDSsequence()
        #self.cds_sequence
        from Bio.Seq import Seq
        from Bio.Alphabet import IUPAC
        cds_seq = Seq(self.cds_sequence, IUPAC.unambiguous_dna)
        peptide = cds_seq.translate()
        return peptide
        
    
class GeneSegment(Base, TableClass):
    """
    2008-07-28
        table to store position info of segments (exon, intron, CDS, UTR ...)
         of gene-products (mRNA, peptide) in table GeneProduct
    """
    __tablename__ = 'gene_segment'
    __table_args__ = {'schema':_schemaname_}
    id = Column(Integer, primary_key=True)
    gene_commentary_id = Column(Integer,
        ForeignKey(f"{_schemaname_}.gene_commentary.id"))
    gene_commentary = relationship('GeneCommentary')
    gi = Column(Integer)
    start = Column(Integer)
    stop = Column(Integer)
    gene_commentary_type_id = Column(Integer,
        ForeignKey(f"{_schemaname_}.gene_commentary_type.id"))
    gene_commentary_type = relationship('GeneCommentaryType')
    created_by = Column(String(256))
    updated_by = Column(String(256))
    date_created = Column(DateTime, default=datetime.now)
    date_updated = Column(DateTime)
    UniqueConstraint('gene_commentary_id', 'start', 'stop',
        'gene_commentary_type_id')

class Gene(Base, TableClass):
    """
    2012.4.26
        change the unique constraint to include annot_assembly_id,
         entrezgene_type_id and remove type_of_gene.
    2010-12-15
        move all EntrezgeneMapping-exclusive elements into Gene.
        EntrezgeneMapping will disappear.
    2008-07-27
        table to store meta info of genes
    """
    __tablename__ = 'gene'
    __table_args__ = {'schema':_schemaname_}
    id = Column(Integer, primary_key=True)
    tax_id = Column(Integer)
    ncbi_gene_id = Column(Integer, unique=True)
    gene_symbol = Column(String(128))
    locustag = Column(String(128))
    synonyms = Column(Text)
    dbxrefs = Column(Text)
    chromosome = Column(String(115))
    strand = Column(String(4))
    start = Column(Integer)
    stop = Column(Integer)
    map_location = Column(String(256))
    description = Column(Text)
    type_of_gene = Column(String(128))
    symbol_from_nomenclature_authority = Column(String(256))
    full_name_from_nomenclature_authority = Column(Text)
    nomenclature_status = Column(String(64))
    other_designations = Column(Text)
    modification_date = Column(DateTime, default=datetime.now)
    genomic_accession = Column(String(32))
    genomic_version = Column(Integer)
    genomic_gi = Column(Integer)
    annot_assembly_id = Column(Integer,
        ForeignKey(f"{_schemaname_}.annot_assembly.id"))
    annot_assembly = relationship('AnnotAssembly')
    entrezgene_type_id = Column(Integer,
        ForeignKey(f"{_schemaname_}.entrezgene_type.id"))
    entrezgene_type = relationship('EntrezgeneType')
    comment = Column(Text)
    gene_commentary_ls = relationship('GeneCommentary',
        back_populates='gene',cascade='all,delete')
    created_by = Column(String(256))
    updated_by = Column(String(256))
    date_created = Column(DateTime, default=datetime.now)
    date_updated = Column(DateTime)
    UniqueConstraint('tax_id', 'locustag', 'chromosome', 'strand', 'start',
        'stop', 'entrezgene_type_id', 'annot_assembly_id')

class Gene2go(Base, TableClass):
    """
    2010-12-15
        add go_qualifier as part of the unique constraint
        "NOT" means the gene is not associated with specified GO.
    2008-07-27
        table to store mapping between gene and GO
    """
    __tablename__ = 'gene2go'
    __table_args__ = {'schema':_schemaname_}
    id = Column(Integer, primary_key=True)
    tax_id = Column(Integer)
    gene_id = Column(Integer,
        ForeignKey(f"{_schemaname_}.gene.id"))
    gene = relationship('Gene')
    go_id = Column(String(32))
    evidence = Column(String(3))
    go_qualifier = Column(String(64))
    go_description = Column(Text)
    pubmed_ids = Column(Text)
    category = Column(String(16))
    created_by = Column(String(256))
    updated_by = Column(String(256))
    date_created = Column(DateTime, default=datetime.now)
    date_updated = Column(DateTime)
    UniqueConstraint('tax_id', 'gene_id', 'go_id', 'evidence', 'category',
        'go_qualifier')

class Gene2Family(Base, TableClass):
    """
    2010-8-19 a table recording which family TE belongs to
    """
    __tablename__ = 'gene2family'
    __table_args__ = {'schema':_schemaname_}
    id = Column(Integer, primary_key=True)
    gene_id = Column(Integer,
        ForeignKey(f"{_schemaname_}.gene.id"))
    gene = relationship('Gene')
    family_id = Column(Integer,
        ForeignKey(f"{_schemaname_}.gene_family.id"))
    family = relationship('GeneFamily')
    created_by = Column(String(256))
    updated_by = Column(String(256))
    date_created = Column(DateTime, default=datetime.now)
    date_updated = Column(DateTime)
    UniqueConstraint('gene_id', 'family_id')


class GeneFamily(Base, TableClass):
    """
    2010-8-19
        family names and etc.
    """
    __tablename__ = 'gene_family'
    __table_args__ = {'schema':_schemaname_}
    id = Column(Integer, primary_key=True)
    short_name = Column(String(212), unique=True)
    parent_family_id = Column(Integer,
        ForeignKey(f"{_schemaname_}.gene_family.id"))
    parent_family = relationship('GeneFamily')
    family_type = Column(String(256))
    description = Column(Text)
    created_by = Column(String(256))
    updated_by = Column(String(256))
    date_created = Column(DateTime, default=datetime.now)
    date_updated = Column(DateTime)

class README(Base, TableClass):
    __tablename__ = 'readme'
    __table_args__ = {'schema':_schemaname_}
    id = Column(Integer, primary_key=True)
    title = Column(Text)
    description = Column(Text)
    created_by = Column(String(256))
    updated_by = Column(String(256))
    date_created = Column(DateTime, default=datetime.now)
    date_updated = Column(DateTime)

class Gene_symbol2id(Base, TableClass):
    """
    2010-6-21
        add created_by, updated_by, date_created, date_updated
    2008-07-27
        a derived table from Gene in order to map all available
         gene names (symbol, locustag, synonym) to gene-id
    """
    __tablename__ = 'gene_symbol2id'
    __table_args__ = {'schema':_schemaname_}
    id = Column(Integer, primary_key=True)
    tax_id = Column(Integer)
    gene_symbol = Column(String(256))
    gene_id = Column(Integer,
        ForeignKey(f"{_schemaname_}.gene.id"))
    gene = relationship('Gene')
    symbol_type = Column(String(64))
    created_by = Column(String(256))
    updated_by = Column(String(256))
    date_created = Column(DateTime, default=datetime.now)
    date_updated = Column(DateTime)
    UniqueConstraint('gene_id', 'gene_symbol', 'symbol_type')

def getEntrezgeneAnnotatedAnchor(db, tax_id):
    """
    2011-1-25
        row.gene_id -> row.id
    2008-08-13
        similar to annot.bin.codense.common.get_entrezgene_annotated_anchor, 
        but use elixir db interface
    """
    sys.stderr.write("Getting entrezgene_annotated_anchor ...")
    chromosome2anchor_gene_tuple_ls = {}
    gene_id2coord = {}
    offset_index = 0
    block_size = 5000
    rows = db.queryTable(Gene).filter_by(tax_id=tax_id).offset(offset_index).\
        limit(block_size)
    while rows.count()!=0:
        for row in rows:
            genomic_gi = row.genomic_gi
            chromosome = row.chromosome
            gene_id = row.id
            strand = row.strand
            start = row.start
            stop = row.stop
            if strand=='1' or strand=='+1' or strand=='+':
                strand = '+'
            elif strand=='-1' or strand=='-':
                strand = '-'
            else:
                strand = strand
            
            if chromosome not in chromosome2anchor_gene_tuple_ls:
                chromosome2anchor_gene_tuple_ls[chromosome] = []
            
            chromosome2anchor_gene_tuple_ls[chromosome].append((start, gene_id))
            gene_id2coord[gene_id] = (start, stop, strand, genomic_gi)
        
        rows = db.queryTable(Gene).filter_by(tax_id=tax_id).\
            offset(offset_index).limit(block_size)
    for chromosome in chromosome2anchor_gene_tuple_ls:
        #sort the list
        chromosome2anchor_gene_tuple_ls[chromosome].sort()
    sys.stderr.write("Done.\n")
    return chromosome2anchor_gene_tuple_ls, gene_id2coord

class OneGenomeData(PassingData):
    """
    chrOrder (order of chromosomes)
        1: order by column genome_order, ascendingly
        2: order by size, descendingly
        3: order by column genome_order, ascendingly
    
    2011-11-28 add two arguments
        self.chrOrder = None
        self.sequence_type_id = None
    2011-3-25
        a structure to wrap data related to one genome
         (identified by tax_id)
        
        All these data used to be handled by class GenomeDatabase.
            Previous solution is not good because GenomeDatabase
            is an umbrella of genome data from different species.
    """
    def __init__(self, db_genome=None, tax_id=3702, maxPseudoChrSize=3000000000, **keywords):
        """
        2012.8.13
            add argument maxPseudoChrSize, 
            which determines the max pseudo chromosome size (for self.chr_id2cumu_chr_start)
        #2011-3-25
        """
        #PassingData.__init__(self)
        self.db_genome = db_genome
        self._chr_id2size = None
        self._chr_id2cumu_size = None
        self._chr_id2cumu_start = None
        self._chr_id2cumu_chr_start = None
        #2012.8.13 the key is chromosome ID.
        # the value is (pseudoChromosomeID, cumuStart).
        # This data structure is used to re-organize the old chromosomes into new/pseudo chromosomes.
        # The max length of each pseudo chromosome is <=self.maxPseudoChrSize. 
        
        #2013.2.17
        self._chr_id2centromere = None
        self._chr_id2annot_assembly = None
        
        self.tax_id = tax_id
        self.maxPseudoChrSize = maxPseudoChrSize
        
        self.chr_gap = None
        self._chr_id_ls = None
        #2013.07.29 changed from chr_id_ls to _chr_id_ls 
        
        self.chrOrder = None
        self.sequence_type_id = None
        self._cumuSpan2ChrRBDict = None
        
        PassingData.__init__(self, **keywords)
        #keywords could contain chr_gap, chrOrder
    
    @property
    def chr_id_ls(self):
        """
        2013.07.29 return or construct a list of chromosome IDs according to user's chrOrder
        """
        if self._chr_id_ls is None:
            self.chr_id_ls = (self.tax_id, self.chrOrder, )
        return self._chr_id_ls
    
    @chr_id_ls.setter
    def chr_id_ls(self, argument_ls=[]):
        """
        #2014.01.12 added outdated_index=0 to filter AnnotAssembly
        2013.07.29
        2014.01.12 set outdated_index to 0
        """
        if argument_ls and len(argument_ls)>0:
            tax_id = argument_ls[0]
            if len(argument_ls)>1:
                chrOrder = argument_ls[1]
            else:
                chrOrder = self.chrOrder
        else:
            tax_id = self.tax_id
            chrOrder = self.chrOrder
        sys.stderr.write("Getting chr_id_ls with tax_id=%s, chrOrder=%s ..."%(tax_id, chrOrder))
        
        query = self.db_genome.queryTable(AnnotAssembly).filter_by(tax_id=tax_id).\
            filter_by(start=1).filter_by(outdated_index=0)
        if self.chrOrder==2:
            query = query.order_by(AnnotAssembly.stop)
        else:	#1 or 3
            query = query.order_by(AnnotAssembly.genome_order)
        if self.sequence_type_id:
            query = query.filter_by(sequence_type_id=self.sequence_type_id)
        
        self._chr_id_ls = []
        
        for row in query:
            self._chr_id_ls.append(row.chromosome)
        sys.stderr.write("%s chromosomes.\n"%(len(self._chr_id_ls)))
    
    @property
    def chr_id2size(self):
        """
        2011-3-13
        """
        if self._chr_id2size is None:
            self.chr_id2size = (self.tax_id, )
        return self._chr_id2size
    
    @chr_id2size.setter
    def chr_id2size(self, argument_ls=[]):
        """
        #2014.01.12 added outdated_index=0 to filter AnnotAssembly
        #2011-11-21
            order the chromosomes according to self.chrOrder
            select sequences based on self.sequence_type_id
        2011-3-12
            modified from get_chr_id2size() of variation/src/common.py
            keywords could include, tax_id=3702.
        2008-10-07 curs could be elixirdb.metadata.bind
        2007-10-12
        """
        if argument_ls and len(argument_ls)>0:
            tax_id = argument_ls[0]
        else:
            tax_id = self.tax_id
        sys.stderr.write("Getting chr_id2size for tax_id %s ..."%(tax_id))
        
        #query = AnnotAssembly.query.filter_by(tax_id=tax_id).filter_by(start=1)
        """
        if self.chrOrder==2:
            orderByString = 'order by stop desc'
        else:
            orderByString = ""
        if self.sequence_type_id:
            extraCondition = " and sequence_type_id=%s "%self.sequence_type_id
        else:
            extraCondition = ""
        
        query = self.db_genome.metadata.bind.execute(
            "select id, chromosome, stop from genome.%s where tax_id=%s and start=1 %s %s"%\
            (AnnotAssembly.table.name, tax_id, extraCondition, orderByString))
        """
        query = self.db_genome.queryTable(AnnotAssembly).filter_by(tax_id=tax_id).\
            filter_by(start=1).filter_by(outdated_index=0)
        if self.chrOrder==2:
            query = query.order_by(AnnotAssembly.stop)
        else:	#1 or 3
            query = query.order_by(AnnotAssembly.genome_order)
        if self.sequence_type_id:
            query = query.filter_by(sequence_type_id=self.sequence_type_id)
        
        self._chr_id2size = {}
        self._chr_id2annot_assembly = {}
        
        for row in query:
            self._chr_id2annot_assembly[row.chromosome] = row
            chr_id = row.chromosome
            size = row.stop
            self._chr_id2size[chr_id] = size
        #self._chr_id2size = chr_id2size
        sys.stderr.write("%s chromosomes.\n"%(len(self._chr_id2size)))
    
    @property
    def chr_id2cumu_size(self):
        """
        2011-3-13
        """
        if self._chr_id2cumu_size is None:
            self.chr_id2cumu_size = (self.tax_id, self.chr_gap)
        return self._chr_id2cumu_size
    
    @chr_id2cumu_size.setter
    def chr_id2cumu_size(self, argument_list):
        """
        2011-3-12
            modified from get_chr_id2cumu_size() of variation/src/common.py
            cumu_size of one chr = sum of length of
             (all prior chromosomes + current chromosome) + all gaps between them
        2008-02-04
            add _chr_id_ls
            turn chr_id all into 'str' form
        2008-02-01
            add chr_gap, copied from variation.src.misc
        2007-10-16
        """
        tax_id, chr_gap = argument_list[:2]
        sys.stderr.write("Getting chr_id2cumu_size for %s ..."%tax_id)
        if self._chr_id2size is None:
            self.chr_id2size = (tax_id,)
        #if chr_gap not specified, take the 1/5th of the average chromosome size as its value
        if chr_gap==None:
            chr_size_ls = self.chr_id2size.values()
            chr_gap = int(sum(chr_size_ls)/(5.0*len(chr_size_ls)))
        
        if self._chr_id_ls is None:
            self.chr_id_ls = (tax_id, self.chrOrder, )
        
        #_chr_id_ls = sorted(self.chr_id2size)
        #no more chromosome 0
        first_chr = self.chr_id_ls[0] 
        #chr_id_ls might not be continuous integers. so dictionary is better
        self._chr_id2cumu_size = {first_chr:self.chr_id2size[first_chr]}
        for i in range(1,len(self.chr_id_ls)):
            chr_id = self.chr_id_ls[i]
            prev_chr_id = self.chr_id_ls[i-1]
            self._chr_id2cumu_size[chr_id] = \
                self._chr_id2cumu_size[prev_chr_id] + chr_gap + self.chr_id2size[chr_id]
        sys.stderr.write("%s chromosomes.\n"%(len(self._chr_id2cumu_size)))
    
    @property
    def chr_id2cumu_start(self):
        """
        2011-3-13
        """
        if self._chr_id2cumu_start is None:
            self.chr_id2cumu_start = (self.tax_id, self.chr_gap)
        return self._chr_id2cumu_start
    
    @chr_id2cumu_start.setter
    def chr_id2cumu_start(self, argument_list):
        """
        chrOrder
            =1: alphabetical order
            =2: by size , descendingly
        
        2011-11-21
            deal with chrOrder
                =1: alphabetical order
                =2: by size , descendingly
        2011-4-22
            cumu_start is now 0-based, which makes it easy to generate adjusted coordinates.
                new_start = cumu_start + start.
        2011-3-15
            Treat one genome of multiple chromsomes as one continuous chromosome.
            
            For one chr, cumu_start = sum of length of all prior chromosomes + all gaps between them.
            
            cumu_start is 0-based.
            
            The difference between chr_id2cumu_size and chr_id2cumu_start is 
                that the former includes the length of
                the current chromosome and the gap before it.
        """
        tax_id, chr_gap = argument_list[:2]
        sys.stderr.write("Getting chr_id2cumu_start for %s, ..."%tax_id)
        if self._chr_id2size is None:
            self.chr_id2size = [tax_id]
        #if chr_gap not specified, take the 1/5th of the average chromosome size as its value
        if chr_gap==None:
            chr_size_ls = self.chr_id2size.values()
            chr_gap = int(sum(chr_size_ls)/(5.0*len(chr_size_ls)))
        
        if self._chr_id_ls is None:
            self.chr_id_ls = (tax_id, self.chrOrder, )
        
        #chr_id_ls = sorted(self.chr_id2size)
        """
        if self.chrOrder==2:
            size_chr_id_ls = [(value, key) for key, value in self.chr_id2size.items()]
            size_chr_id_ls.sort()
            size_chr_id_ls.reverse()
            chr_id_ls = [row[1] for row in size_chr_id_ls]
        """
        first_chr = self.chr_id_ls[0] 
        self._chr_id2cumu_start = {first_chr:0}
        #chr_id_ls might not be continuous integers. so dictionary is better
        #start from 0.
        for i in range(1, len(self.chr_id_ls)):
            chr_id = self.chr_id_ls[i]
            prev_chr_id = self.chr_id_ls[i-1]
            self._chr_id2cumu_start[chr_id] = \
                self._chr_id2cumu_start[prev_chr_id] + chr_gap + self.chr_id2size[prev_chr_id]
        sys.stderr.write("%s chromosomes.\n"%(len(self._chr_id2cumu_start)))
    
    @property
    def chr_id2cumu_chr_start(self):
        """
        2012.8.13
        """
        if self._chr_id2cumu_chr_start is None:
            self.chr_id2cumu_chr_start = (self.maxPseudoChrSize, self.tax_id, self.chr_gap)
        return self._chr_id2cumu_chr_start
    
    @chr_id2cumu_chr_start.setter
    def chr_id2cumu_chr_start(self, argument_list):
        """
        2012.8.13
        this structure is slightly different from chr_id2cumu_start.
            chr_id2cumu_start is for grouping all chromosomes into one giant chromosomes.
            chr_id2cumu_chr_start is for groupping all chromosomes into several chromosomes,
                size of each is <= self.maxPseudoChrSize. 
        """
        
        maxPseudoChrSize, tax_id, chr_gap = argument_list[:3]
        sys.stderr.write("Getting chr_id2cumu_chr_start for tax-id=%s, "
            " chr-gap=%s, maxPseudoChrSize=%s ..."%(
                tax_id, chr_gap, maxPseudoChrSize))
        if self._chr_id2size is None:
            self.chr_id2size = (tax_id, )
        #if chr_gap not specified, take the 1/5th of the average chromosome size as its value
        if chr_gap==None:
            chr_size_ls = self.chr_id2size.values()
            chr_gap = int(sum(chr_size_ls)/(5.0*len(chr_size_ls)))
        
        if self._chr_id_ls is None:
            self.chr_id_ls = (tax_id, self.chrOrder, )
        """
        chr_id_ls = sorted(self.chr_id2size)
        if self.chrOrder==2:
            size_chr_id_ls = [(value, key) for key, value in self.chr_id2size.items()]
            size_chr_id_ls.sort()
            size_chr_id_ls.reverse()
            chr_id_ls = [row[1] for row in size_chr_id_ls]
        """
        first_chr = self.chr_id_ls[0]
        newChrID = 1
        self._chr_id2cumu_chr_start = {first_chr:(newChrID, 0)}
            #chr_id_ls might not be continuous integers. so dictionary is better
            #start from 0.
        for i in range(1, len(self.chr_id_ls)):
            chr_id = self.chr_id_ls[i]
            prev_chr_id = self.chr_id_ls[i-1]
            prev_new_chr, prev_cumu_start = self._chr_id2cumu_chr_start[prev_chr_id]
            new_cumu_start = prev_cumu_start + chr_gap + self.chr_id2size[prev_chr_id]
            if new_cumu_start>=self.maxPseudoChrSize:
                newChrID += 1	#start a new fake chromosome
                new_cumu_start =0	#reset
            else:
                newChrID = prev_new_chr
            self._chr_id2cumu_chr_start[chr_id] = (newChrID, new_cumu_start)
        sys.stderr.write("%s chromosomes, last newChrID=%s.\n"%(len(self._chr_id2cumu_chr_start), newChrID))
    
    @property
    def chr_id2centromere(self):
        """
        2013.2.17
        """
        if self._chr_id2centromere is None:
            self.chr_id2centromere = (self.maxPseudoChrSize, self.tax_id, self.chr_gap)
        return self._chr_id2centromere
    
    @chr_id2centromere.setter
    def chr_id2centromere(self, argument_list):
        """
        2013.2.17
        """
        
        maxPseudoChrSize, tax_id, chr_gap = argument_list[:3]
        sys.stderr.write("Getting chr_id2centromere for tax-id=%s, chr-gap=%s, maxPseudoChrSize=%s ..."%\
                        (tax_id, chr_gap, maxPseudoChrSize))
        if self._chr_id2cumu_start is None:
            self.chr_id2cumu_start = (self.tax_id, self.chr_gap)
        
        self._chr_id2centromere ={}
        for chr_id, annot_assembly in self._chr_id2annot_assembly.items():
            for genome_annotation in annot_assembly.genome_annotation_list:
                if genome_annotation.genome_annotation_type_id==1:
                    self._chr_id2centromere[chr_id] = genome_annotation
        sys.stderr.write("%s centromeres.\n"%(len(self._chr_id2centromere)))
    
    @property
    def cumuSpan2ChrRBDict(self):
        """
        2011-3-25
        """
        if self._cumuSpan2ChrRBDict is None:
            self.cumuSpan2ChrRBDict = (self.tax_id, self.chr_gap)
        return self._cumuSpan2ChrRBDict
    
    @cumuSpan2ChrRBDict.setter
    def cumuSpan2ChrRBDict(self, argument_list):
        """
        2011-4-22
            adjust because chr_id2cumu_start is now 0-based.
            the position in cumuSpan2ChrRBDict is 1-based.
        2011-3-25
            turn it into a setter of cumuSpan2ChrRBDict
            usage example:
                tax_id = 3702
                chr_gap = 0
                self.cumuSpan2ChrRBDict = (tax_id, chr_gap)
            
        2011-3-16
            Treat one genome of multiple chromsomes as one continuous chromosome (uber-chomosome).
            This function creates a RB dictionary, which holds a map between 
                uberchromosome coordinate and individual chromosome coordinate.
        """
        tax_id, chr_gap = argument_list[:2]
        sys.stderr.write("Creating cumuSpan2ChrRBDict for tax_id=%s, chr_gap=%s \n"%(
            tax_id, chr_gap))
        self._cumuSpan2ChrRBDict = RBDict()
        if self._chr_id2size is None:
            self.chr_id2size = [tax_id,]
        if self._chr_id2cumu_start is None:
            self.chr_id2cumu_start = (tax_id, chr_gap)
        for chr_id, cumu_start in self.chr_id2cumu_start.items():
            chr_size = self.chr_id2size.get(chr_id)
            span_ls=[cumu_start+1, cumu_start+chr_size]
            segmentKey = CNVSegmentBinarySearchTreeKey(
                chromosome=0,
                span_ls=span_ls,
                min_reciprocal_overlap=0.00000000000001,)
                #2010-8-17 overlapping keys are regarded as separate instances
                #  as long as they are not identical.
            if segmentKey not in self._cumuSpan2ChrRBDict:
                self._cumuSpan2ChrRBDict[segmentKey] = [chr_id, 1, chr_size]
            else:
                logging.error("%s of chr %s is already in cumuSpan2ChrRBDict."%\
                    (segmentKey, chr_id))
        sys.stderr.write("%s chromosomes.\n"%(len(self._cumuSpan2ChrRBDict)))

class GenomeDatabase(Database):
    __doc__ = __doc__
    option_default_dict = Database.option_default_dict.copy()
    option_default_dict[('drivername', 1,)][0] = 'postgresql'
    option_default_dict.update({
        ('geneAnnotationPickleFname', 0, ): ['', 'o', 1, 
            'The optional output file to contain a pickled object with three attributes: \n'
            '1. gene_id2model \n'
            '2. chr_id2gene_id_ls \n'
            '3. geneSpanRBDict. \n'
            'OPTIONAL: only used when you run GenomeDB.py as a standalone program'
            ],
        ('tax_id', 0, int): [3702, 't', 1, 'Taxonomy ID for geneAnnotationPickleFname.'],\
        })
    def __init__(self, **keywords):
        """
        wrap _chr_id2size, _chr_id2cumu_size, _chr_id2cumu_start etc. into OneGenomeData class.
        """
        Database.__init__(self, **keywords)
        #2011-3-25
        self.tax_id2genomeData = {}
    
    
    def get_gene_id2model(self, tax_id=3702, needSequence=False):
        """
        argument needSequence: whether CDS_sequence,
            mrna_sequence is fetched for each gene_commentary
        2010-10-11
            bug fix: several genes could have the same segmentKey in geneSpanRBDict.
        2009-1-3
            add cds_sequence, mrna_sequence to new_gene_commentary
        2008-10-01
            sort EntrezgeneMapping by chromosome, start
        2008-10-1
            similar to variation.src.GenomeBrowser.get_gene_id2model()
             but adapts to the new genome db schema
            construct a data structure to hold whole-genome annotation of genes
        """
        sys.stderr.write("Getting gene_id2model and chr_id2gene_id_ls ...\n")
        gene_id2model = {}
        chr_id2gene_id_ls = {}
        geneSpanRBDict = RBDict()
        i = 0
        block_size = 50000
        query = self.queryTable(Gene).filter_by(tax_id=tax_id).\
            order_by(Gene.chromosome).order_by(Gene.start)
        rows = query.offset(i).limit(block_size)
        while rows.count()!=0:
            for row in rows:
                chromosome = row.chromosome
                gene_id = row.id
                
                if chromosome not in chr_id2gene_id_ls:
                    chr_id2gene_id_ls[chromosome] = []
                chr_id2gene_id_ls[chromosome].append(gene_id)
                if gene_id not in gene_id2model:
                    gene_id2model[gene_id] = GeneModel(gene_id=gene_id, 
                        ncbi_gene_id=row.ncbi_gene_id, chromosome=chromosome,
                        gene_symbol=row.gene_symbol,\
                        locustag=row.locustag, map_location=row.map_location,\
                        type_of_gene=row.entrezgene_type.type,
                        type_id=row.entrezgene_type_id,\
                        start=row.start, stop=row.stop, strand=row.strand,
                        tax_id=row.tax_id,
                        description =row.description)
                        #2010-8-19 add description
                    
                    gene_model = gene_id2model[gene_id]
                    if gene_model.chromosome and gene_model.start and gene_model.stop:
                        segmentKey = CNVSegmentBinarySearchTreeKey(
                            chromosome=gene_model.chromosome, \
                            span_ls=[gene_model.start, gene_model.stop], \
                            min_reciprocal_overlap=1, strand=gene_model.strand,
                            gene_id=gene_model.gene_id,
                            gene_start=gene_model.start,
                            gene_stop=gene_model.stop)
                        #2010-8-17 any overlap short of being identical is tolerated.
                        if segmentKey not in geneSpanRBDict:
                            geneSpanRBDict[segmentKey] = []
                        geneSpanRBDict[segmentKey].append(gene_id)
                
                    for gene_commentary in row.gene_commentary_ls:
                        if not gene_commentary.gene_commentary_id:
                            #ignore gene_commentary that are derived from other gene_commentary_ls.
                            #  they'll be handled within the parental gene_commentary.
                            #gene_commentary.constructAnnotatedBox()
                            gene_commentary.construct_annotated_box()
                            #2010-9-21 it sets protein_label and etc.
                            #  constructAnnotatedBox() doesn't do this.
                            if needSequence:
                                gene_commentary.getCDSsequence()
                            #pass everything to a database-independent object, easy for pickling
                            new_gene_commentary = GeneModel(gene_id=gene_id, \
                                gene_commentary_id=gene_commentary.id, \
                                start=gene_commentary.start, stop=gene_commentary.stop, \
                                label=gene_commentary.label, comment=gene_commentary.comment, \
                                text=gene_commentary.text,\
                                gene_commentary_type=gene_commentary.gene_commentary_type.type,\
                                protein_label=gene_commentary.protein_label,\
                                protein_comment=gene_commentary.protein_comment,\
                                protein_text=gene_commentary.protein_text,\
                                mrna_box_ls=gene_commentary.mrna_box_ls,\
                                protein_box_ls=gene_commentary.protein_box_ls,\
                                box_ls=gene_commentary.box_ls,
                                cds_sequence = getattr(gene_commentary, 'cds_sequence', None),
                                mrna_sequence = getattr(gene_commentary, 'mrna_sequence', None))
                            gene_id2model[gene_id].gene_commentary_ls.append(new_gene_commentary)
                else:
                    logging.error(f"gene {gene_id} already exists in gene_id2model.")
                i += 1
            if self.report:
                sys.stderr.write("%s\t%s\t"%('\x08'*80, i))
            rows = query.offset(i).limit(block_size)
        
        sys.stderr.write("Done.\n")
        return gene_id2model, chr_id2gene_id_ls, geneSpanRBDict
    
    def createGenomeRBDict(self, tax_id=3702, max_distance=20000, debug=False):
        """
        2012.5.15
            use gene_commentary.construct_annotated_box() instead of
                gene_commentary.constructAnnotatedBox()
            one more step toward deprecating get_gene_id2model()
        2012.3.19
            add an attribute to geneSpanRBDict:
                genomeRBDict.genePadding (=max_distance)
        2011-3-24
            stop casting row.chromosome into integer (just string type)
            add ncbi_gene_id to oneGeneData
        2011-3-10
            copied from FindCNVContext.py
        2011-1-27
            require Gene.start, Gene.stop, Gene.chromosome not null.
            Gene.id is the new gene_id (was Gene.gene_id).
        2010-10-3
            bug fixed: (chr, start, stop) is not unique.
            There are genes with the same coordinates.
        """
        sys.stderr.write("Creating a RBDict for all genes from organism %s ... \n"%tax_id)
        genomeRBDict = RBDict()
        genomeRBDict.genePadding = max_distance
        query = self.queryTable(Gene).filter_by(tax_id=tax_id).filter(Gene.start!=None).\
            filter(Gene.stop!=None).filter(Gene.chromosome!=None)
        counter = 0
        real_counter = 0
        for row in query:
            #try:	# convert to integer except when "C" or "M"/mitochondria is encountered.
            #	chromosome = int(row.chromosome)
            # #integer chromosomes should be converted as CNV.chromosome is integer.
            #except:
            chromosome = row.chromosome
            segmentKey = CNVSegmentBinarySearchTreeKey(chromosome=chromosome, \
                span_ls=[max(1, row.start - max_distance), row.stop + max_distance], \
                min_reciprocal_overlap=1,)
                #2010-8-17 overlapping keys are regarded as separate instances as long as they are not identical.
            if segmentKey not in genomeRBDict:
                genomeRBDict[segmentKey] = []
            oneGeneData = PassingData(strand = row.strand, gene_id = row.id, \
                gene_start = row.start, \
                gene_stop = row.stop, geneCommentaryRBDictLs=[],\
                ncbi_gene_id=row.ncbi_gene_id, type_of_gene=row.type_of_gene)
            counter += 1
            for gene_commentary in row.gene_commentary_ls:
                if not gene_commentary.gene_commentary_id:
                    # ignore gene_commentary that are derived from other gene_commentary_ls. 
                    # they'll be handled within the parental gene_commentary.
                    geneCommentaryRBDict = RBDict()
                    geneCommentaryRBDict.gene_commentary_id = gene_commentary.id
                    box_ls = gene_commentary.construct_annotated_box()
                    geneCommentaryRBDict.box_ls = box_ls
                    geneCommentaryRBDict.protein_box_ls= gene_commentary.protein_box_ls
                    gene_commentary.getCDSsequence()
                    geneCommentaryRBDict.cds_sequence = gene_commentary.cds_sequence
                    geneCommentaryRBDict.gene_commentary_type_name = gene_commentary.gene_commentary_type.type
                    geneCommentaryRBDict.start = gene_commentary.start
                    geneCommentaryRBDict.stop = gene_commentary.stop
                    geneCommentaryRBDict.CDS_5_end_pos = gene_commentary.CDS_5_end_pos
                    geneCommentaryRBDict.no_of_introns = gene_commentary.no_of_introns
                    
                    #box_ls = gene_commentary.constructAnnotatedBox()
                    #2012.5.15 constructAnnotatedBox()
                    #  requires GeneSegment to have UTR,exon,CDS,intron already in place
                    #box_ls = gene_commentary.box_ls
                    no_of_boxes = len(box_ls)
                    
                    numberPorter = PassingData(cds_number = 0,\
                                            intron_number = 0,\
                                            utr_number = 0,\
                                            exon_number = 0)
                    for i in range(no_of_boxes):
                        if row.strand == "-1":	#reverse
                            box = box_ls[-i-1]
                        else:
                            box = box_ls[i]
                        start, stop, box_type, is_translated, protein_box_index, \
                            gene_segment_id, detailed_box_type, \
                            cumulativeWithinCDSUTRAndIntronLen, exon_number = box[:9]
                        
                        geneSegmentKey = CNVSegmentBinarySearchTreeKey(
                            chromosome=chromosome, 
                            span_ls=[start, stop], \
                            min_reciprocal_overlap=1, label=detailed_box_type,
                            is_translated=is_translated,
                            protein_box_index=protein_box_index,utr_number=None,\
                            cds_number=None, intron_number=None, exon_number=None, \
                            gene_segment_id=gene_segment_id, 
                            cumulativeWithinCDSUTRAndIntronLen=cumulativeWithinCDSUTRAndIntronLen)
                        #2010-8-17 overlapping keys are regarded as
                        #  separate instances as long as they are not identical.
                        
                        if detailed_box_type=='3UTR' or detailed_box_type=='5UTR' or detailed_box_type=='UTR':
                            numberPorter.utr_number += 1
                            geneSegmentKey.utr_number = numberPorter.utr_number
                        elif detailed_box_type=='CDS':
                            numberPorter.cds_number += 1
                            geneSegmentKey.cds_number = numberPorter.cds_number
                        elif detailed_box_type=='intron':
                            numberPorter.intron_number += 1
                            geneSegmentKey.intron_number = numberPorter.intron_number
                        if box_type=='exon':
                            #detailed_box_type could be exon only for non-coding genes. 
                            geneSegmentKey.exon_number = exon_number
                        
                        geneCommentaryRBDict[geneSegmentKey] = None
                        real_counter += 1
                    oneGeneData.geneCommentaryRBDictLs.append(geneCommentaryRBDict)
            genomeRBDict[segmentKey].append(oneGeneData)
            if counter%1000==0:
                sys.stderr.write("%s%s\t%s"%('\x08'*100, counter, real_counter))
                if debug:
                    break
        sys.stderr.write("%s%s\t%s\n"%('\x08'*100, counter, real_counter))
        sys.stderr.write("%s Done.\n"%(str(genomeRBDict)))
        return genomeRBDict
    
    def dealWithGenomeRBDict(self, genomeRBDictPickleFname, tax_id=3702, \
        max_distance=20000, debug=None):
        """
        2011-3-10
        """
        sys.stderr.write("Dealing with genomeRBDict ...")
        import pickle
        if genomeRBDictPickleFname:
            if os.path.isfile(genomeRBDictPickleFname):
                #if this file is already there, suggest to un-pickle it.
                picklef = open(genomeRBDictPickleFname)
                genomeRBDict = pickle.load(picklef)
                del picklef
            else:
                #if the file doesn't exist, but the filename is given,
                #  pickle snps_context_wrapper into it
                genomeRBDict = self.createGenomeRBDict(tax_id=tax_id, \
                    max_distance=max_distance, debug=debug)
                #2008-09-07 pickle the snps_context_wrapper object
                picklef = open(genomeRBDictPickleFname, 'w')
                pickle.dump(genomeRBDict, picklef, -1)
                picklef.close()
        else:
            genomeRBDict = self.createGenomeRBDict(tax_id=tax_id, \
                max_distance=max_distance, debug=debug)
        sys.stderr.write("%s unique genomic spans. Done.\n"%(len(genomeRBDict)))
        return genomeRBDict
    
    def getOneGenomeData(self, tax_id=3702, chr_gap=0, chrOrder=1, \
        sequence_type_id=None, \
        maxPseudoChrSize=3000000000, **keywords):
        """
        chrOrder (order of chromosomes)
            1: order by column genome_order, ascendingly
            2: order by size, descendingly
            3: order by column genome_order, ascendingly
        
        2012.8.13
            add argument maxPseudoChrSize, keywords
        2011-11-21
            add argument chrOrder (order of chromosomes)
                1: whatever order in database
                2: order by size, descendingly
        2011-3-25
            API to get & set OneGenomeData for one taxonomy.
        """
        if tax_id not in self.tax_id2genomeData:
            oneGenomeData = OneGenomeData(db_genome=self, tax_id=tax_id, \
                chr_gap=chr_gap, chrOrder=chrOrder,\
                sequence_type_id=sequence_type_id,
                maxPseudoChrSize=maxPseudoChrSize, **keywords)
            self.tax_id2genomeData[tax_id] = oneGenomeData
        return self.tax_id2genomeData.get(tax_id)
    
    def outputGenomeSequence(self, tax_id=None, sequence_type_id=1, \
        fastaTitlePrefix='chr', outputDir=None, chunkSize=70):
        """
        #2014.01.12 added outdated_index=0 to filter AnnotAssembly
        2011-7-7
at this moment, The output is with symap in mind.

i.e.:
outputDir = os.path.expanduser("~/script/variation/bin/symap_3_3/data/pseudo/vervet177BAC/sequence/pseudo/")
db_genome.outputGenomeSequence(tax_id=60711, sequence_type_id=10,
    fastaTitlePrefix='BAC', outputDir=outputDir, chunkSize=70)
        """
        sys.stderr.write("Outputting genome sequences from tax %s, seq-type %s to %s ...\n"\
            %(tax_id, sequence_type_id, outputDir))
        query = self.queryTable(AnnotAssembly).filter_by(tax_id=tax_id).\
            filter_by(sequence_type_id=sequence_type_id).\
            filter_by(outdated_index=0).order_by(AnnotAssembly.id)
        for row in query:
            fastaTitle = '%s%s'%(fastaTitlePrefix, row.id)
            outputFname = os.path.join(outputDir, '%s.seq'%(fastaTitle))
            outf = open(outputFname, 'w')
            outf.write(">%s\n"%(fastaTitle))
            sub_query = self.queryTable(RawSequence).filter_by(annot_assembly_id=row.id).\
                order_by(RawSequence.start)
            for seq_row in sub_query:
                sequence = seq_row.sequence
                while len(sequence)>0:
                    outf.write("%s\n"%(sequence[:chunkSize]))
                    sequence = sequence[chunkSize:]
            del outf
        sys.stderr.write("Done.\n")
    
    def getSequenceSegment(self, tax_id=None, chromosome=None, start=None,
        stop=None, schema='genome'):
        """
        #2014.01.12 added outdated_index=0 to filter AnnotAssembly
        2011-10-24
            start and stop are 1-based.
        """
        chromosome = str(chromosome)
        query = self.queryTable(AnnotAssembly).filter_by(tax_id=tax_id).\
            filter_by(chromosome=chromosome).filter_by(outdated_index=0).\
            order_by(AnnotAssembly.id)
        annot_assembly_id_ls = []
        for row in query:
            annot_assembly_id_ls.append(row.id)
        if len(annot_assembly_id_ls)>1:
            logging.warn("more than 1 entries from tax_id %s, chromosome %s. take first one."\
                %(tax_id, chromosome))
        elif len(annot_assembly_id_ls)==0:
            logging.warn("No entry available from tax_id %s, chromosome %s."%(\
                tax_id, chromosome))
            return ""
        annot_assembly_id = annot_assembly_id_ls[0]
        curs = __metadata__.bind
        #or self.table.bind or self.table.metadata.bind
        seq = get_sequence_segment(curs, gi=None, start=start, stop=stop, \
            annot_assembly_id=annot_assembly_id,\
            annot_assembly_table='%s.%s'%(schema, AnnotAssembly.__tablename__), \
            raw_sequence_table='%s.%s'%(schema, RawSequence.__tablename__))
        return seq
    
    def getSequenceType(self, short_name=None, entry_id=None, **keywords):
        """
        2013.3.14
        """
        db_entry = self.checkIfEntryInTable(TableClass=SequenceType, 
            short_name=short_name, entry_id=entry_id)
        if not db_entry:
            db_entry = SequenceType(short_name=short_name)
            self.session.add(db_entry)
            self.session.flush()
        return db_entry
    
    def getTopNumberOfChomosomes(self, contigMaxRankBySize=100, \
        contigMinRankBySize=1, tax_id=60711, sequence_type_id=9,\
        version=None, chromosome_type_id=0, chromosome=None, outdated_index=0):
        """
        #2014.01.12 added outdated_index=0 to filter AnnotAssembly
        2013.3.19 added argument outdated_index
        2013.3.14 added argument version, chromosome
        2013.2.15 added argument chromosome_type_id
        chromosome_type_id: what type of chromosomes to be included,
             same as table genome.chromosome_type.
            0: all, 1: autosomes, 2: X, 3:Y, 4: mitochondrial
        2011-11-7
            copied from vervet/src/AlignmentToCallPipeline.py
        2011-11-6
            rename argument topNumberOfContigs to contigMaxRankBySize
            add argument contigMinRankBySize
        2011-9-13
            return chr2size instead of a set of ref names
        2011-7-12
            get all the top contigs
        """
        no_of_contigs_to_fetch = contigMaxRankBySize-contigMinRankBySize+1
        print(f"Getting {no_of_contigs_to_fetch} chromosomes with rank "
            f"(by size) between {contigMinRankBySize} and {contigMaxRankBySize}:",
            flush=True)
        chr2size = {}
        query = self.queryTable(AnnotAssembly).filter_by(tax_id=tax_id).\
            filter_by(sequence_type_id=sequence_type_id).\
            filter_by(outdated_index=0)
        if chromosome_type_id:
            query = query.filter_by(chromosome_type_id=chromosome_type_id)
        if version is not None:
            query = query.filter_by(version=version)
        if chromosome:
            query = query.filter_by(chromosome=chromosome)
        if outdated_index is not None:
            query = query.filter_by(outdated_index=outdated_index)
        #order by size from big to small
        query = query.order_by(desc(AnnotAssembly.stop))
        counter = 0
        for row in query:
            counter += 1
            if counter>=contigMinRankBySize and counter<=contigMaxRankBySize:
                chr2size[row.chromosome] = row.stop
            if len(chr2size)>=no_of_contigs_to_fetch:
                break
        print(f"{len(chr2size)} chromosomes.", flush=True)
        return chr2size
    
    def getChromosomeType(self, short_name=None, id=None, comment=None,
        **keywords):
        """
        2013.3.14
        """
        TableClass=ChromosomeType
        db_entry = self.checkIfEntryInTable(TableClass=TableClass, \
            short_name=short_name, entry_id=id)
        if not db_entry:
            db_entry = TableClass(short_name=short_name, comment=comment)
            self.session.add(db_entry)
            self.session.flush()
        return db_entry
    
    def checkAnnotAssembly(self, id=None, accession = None, \
        version=None, tax_id=None, chromosome =None, \
        start=None, stop =None, orientation=None, \
        sequence_type_id=None, chromosome_type_id=None,\
        outdated_index=0):
        """
        2013.3.15 added argument outdated_index
        2013.3.14
            all arguments part of unique key
        """
        if id is not None:
            db_entry = self.checkIfEntryInTable(
                TableClass=AnnotAssembly, entry_id=id)
            if db_entry:
                return db_entry
        
        query = self.queryTable(AnnotAssembly)
        if accession:
            query = query.filter_by(accession=accession)
        if version:
            query = query.filter_by(version=version)
        if tax_id:
            query = query.filter_by(tax_id=tax_id)
        if chromosome:
            query = query.filter_by(chromosome=chromosome)
        if start:
            query = query.filter_by(start=start)
        if stop:
            query = query.filter_by(stop=stop)
        if sequence_type_id:
            query = query.filter_by(sequence_type_id=sequence_type_id)
        if chromosome_type_id:
            query = query.filter_by(chromosome_type_id=chromosome_type_id)
        if outdated_index:
            query = query.filter_by(outdated_index=outdated_index)
        db_entry = query.first()
        return db_entry
    
    def getAnnotAssembly(self, id=None, gi=None, acc_ver=None,
        accession = None, version =None, tax_id=None, chromosome =None,
        start =1, stop =None, orientation=None, sequence = None,\
        raw_sequence_start_id=None, original_path=None, \
        sequence_type_id=None, sequence_type_name= None, \
        chromosome_type_id=None, chromosome_type_name=None, \
        outdated_index=0,\
        comment=None, **keywords):
        """
        2013.3.15 added argument outdated_index
        2013.3.14 more fully-fledged
        2012.4.26
            return db entry of a chromosome
        """
        if not sequence_type_id:
            if sequence_type_name:
                sequence_type = self.getSequenceType(short_name=sequence_type_name)
                sequence_type_id=sequence_type.id
        if not chromosome_type_id:
            if chromosome_type_name:
                chromosome_type = self.getChromosomeType(short_name=chromosome_type_name)
                chromosome_type_id=chromosome_type.id
        db_entry=self.checkAnnotAssembly(id=id, accession=accession, \
            version=version, tax_id=tax_id, \
            chromosome=chromosome, start=start, stop=stop,
            orientation=orientation,
            sequence_type_id=sequence_type_id,
            chromosome_type_id=chromosome_type_id,
            outdated_index=outdated_index)
        if not db_entry:
            db_entry = AnnotAssembly(gi =gi, acc_ver=acc_ver, accession=accession,
                version=version, tax_id = tax_id, chromosome = chromosome,
                start = start, stop =stop, orientation=orientation,
                sequence=sequence, raw_sequence_start_id = raw_sequence_start_id,
                original_path=original_path,
                sequence_type_id = sequence_type_id,
                chromosome_type_id=chromosome_type_id,
                outdated_index=outdated_index,\
                comment = comment)
            self.session.add(db_entry)
            self.session.flush()
        return db_entry
    
    def checkGenomeAnnotation(self, entry_id=None, strand=None, start=None,
        stop =None, annot_assembly_id=None,
        genome_annotation_type_id=None, description=None, comment=None):
        """
        2013.08.01
        """
        if entry_id is not None:
            db_entry = self.checkIfEntryInTable(
                TableClass=GenomeAnnotation, entry_id=entry_id)
            if db_entry:
                return db_entry
        
        query = self.queryTable(GenomeAnnotation)
        if strand:
            query = query.filter_by(strand=strand)
        if start:
            query = query.filter_by(start=start)
        if stop:
            query = query.filter_by(stop=stop)
        if annot_assembly_id:
            query = query.filter_by(annot_assembly_id=annot_assembly_id)
        if genome_annotation_type_id:
            query = query.filter_by(
                genome_annotation_type_id=genome_annotation_type_id)
        db_entry = query.first()
        return db_entry
    
    def getGenomeAnnotation(self, entry_id=None, strand=None,
        start=None, stop =None,
        annot_assembly_id=None, \
        tax_id=None, chromosome=None, version=None, \
        sequence_type_name=None, sequence_type_id=None, \
        genome_annotation_type_id=None, genome_annotation_type_name=None,
        description=None, comment=None, **keywords):
        """
        2013.8.1
        """
        if not annot_assembly_id:
            if tax_id and chromosome:
                annot_assembly = self.getAnnotAssembly(id=None, gi=None,
                    acc_ver=None, accession=None, version=version,
                    tax_id=tax_id, chromosome=chromosome, start=None, stop=None,
                    orientation=None, sequence=None, raw_sequence_start_id=None,
                    original_path=None, sequence_type_id=sequence_type_id,
                    sequence_type_name=sequence_type_name,
                    chromosome_type_id=None,
                    chromosome_type_name=None, outdated_index=None,
                    comment=None)
                annot_assembly_id = annot_assembly.id
        if not genome_annotation_type_id:
            if genome_annotation_type_name:
                genome_annotation_type = self.getGenomeAnnotationType(
                    short_name=genome_annotation_type_name, description=None,
                    comment=None)
                genome_annotation_type_id = genome_annotation_type.id
        db_entry=self.checkGenomeAnnotation(entry_id=entry_id, strand=strand,
            start=start, stop=stop,
            annot_assembly_id=annot_assembly_id, \
            genome_annotation_type_id=genome_annotation_type_id, \
            description=description, comment=comment)
        if not db_entry:
            db_entry = GenomeAnnotation(strand = strand, start = start,\
                stop =stop, annot_assembly_id=annot_assembly_id, \
                genome_annotation_type_id = genome_annotation_type_id, \
                description=description, comment = comment)
            self.session.add(db_entry)
            self.session.flush()
        return db_entry
    
    def checkGenomeAnnotationType(self, entry_id=None, short_name=None):
        """
        2013.08.01
        """
        if id is not None:
            db_entry = self.checkIfEntryInTable(
                TableClass=GenomeAnnotationType, entry_id=entry_id)
            if db_entry:
                return db_entry
        
        query = self.queryTable(GenomeAnnotationType)
        if short_name:
            query = query.filter_by(short_name=short_name)
        db_entry = query.first()
        return db_entry
    
    def getGenomeAnnotationType(self, entry_id=None, short_name=None, \
        description=None, comment=None, **keywords):
        """
        2013.8.1
        """
        db_entry=self.checkGenomeAnnotationType(entry_id=entry_id,
            short_name=short_name)
        if not db_entry:
            db_entry = GenomeAnnotationType(short_name = short_name, \
                description=description, comment = comment)
            self.session.add(db_entry)
            self.session.flush()
        return db_entry
    
    def getGeneCommentary(self, gene_id=None, start=None, stop=None,
        gene_commentary_type_id=None, gene_commentary_id=None,\
        accession=None, version=None, gi=None, label=None, text=None,
        comment=None,
        gene=None, gene_commentary=None):
        """
        2012.4.26
        """
        if gene_id is None and gene.id is not None:
            gene_id = gene.id
        if gene_commentary_id is None and gene_commentary is not None:
            gene_commentary_id = gene_commentary.id
        
        if gene_id and start and stop and gene_commentary_type_id:
            query = self.queryTable(GeneCommentary)
            query = query.filter_by(gene_id=gene_id).filter_by(start=start).\
                filter_by(stop=stop).\
                filter_by(gene_commentary_type_id=gene_commentary_type_id)
            if gene_commentary_id:
                query = query.filter_by(gene_commentary_id=gene_commentary_id)
            geneCommentary = query.first()
        else:
            geneCommentary = None
        if not geneCommentary:
            geneCommentary = GeneCommentary(start=start, stop=stop, \
                gene_id=gene_id, gene_commentary_type_id=gene_commentary_type_id,
                gene_commentary_id=gene_commentary_id, accession=accession,
                version=version,\
                gi=gi, label=label,text=text, comment=comment)
            geneCommentary.gene = gene
            geneCommentary.gene_commentary = gene_commentary
            self.session.add(geneCommentary)
            self.session.flush()
        return geneCommentary
    
    def getGeneCommentaryType(self, gene_commentary_type_id=None,
        commentary_type=None):
        """
        2012.5.15
            moved from transfac/src/GeneASNXML2gene_mapping.py
        """
        if gene_commentary_type_id:
            gene_commentary_type = self.checkIfEntryInTable(
                TableClass=GeneCommentaryType, entry_id=gene_commentary_type_id)
        elif commentary_type:
            gene_commentary_type = self.queryTable(GeneCommentaryType).\
                filter_by(type=commentary_type).first()
        else:
            gene_commentary_type = None
        if not gene_commentary_type:
            if self.debug:
                logging.warn("\t Gene-commentary_type %s not in db yet."%\
                    commentary_type)
            gene_commentary_type = GeneCommentaryType(type=commentary_type)
            if gene_commentary_type_id:
                gene_commentary_type.id = gene_commentary_type_id
            self.session.add(gene_commentary_type)
            self.session.flush()
        return gene_commentary_type
    
    def getGene(self, annot_assembly_id=None, locustag=None, tax_id=None, 
        chromosome=None, strand=None,\
        start=None, stop=None, type_of_gene=None, entrezgene_type_id=None, 
        chromosome_sequence_type_id=None,
        annot_assembly=None, ncbi_gene_id=None, gene_symbol=None, 
        synonyms=None, dbxrefs=None, map_location=None, description=None, 
        genomic_accession=None, genomic_version=None, genomic_gi=None, \
        symbol_from_nomenclature_authority=None,
        full_name_from_nomenclature_authority=None,\
        nomenclature_status=None, modification_date=None):
        """
        2012.4.26
using_table_options(UniqueConstraint('tax_id', 'locustag', 'chromosome', 
    'strand', 'start', 'stop', 'type_of_gene'))
        """
        if annot_assembly_id is None and annot_assembly:
            annot_assembly_id = annot_assembly.id
        if annot_assembly_id is None and chromosome and tax_id  and chromosome_sequence_type_id:
            annot_assembly = self.getAnnotAssembly(tax_id=tax_id, \
                sequence_type_id=chromosome_sequence_type_id, chromosome=chromosome)
            annot_assembly_id = annot_assembly.id
        
        if annot_assembly_id and start and stop and tax_id and chromosome and entrezgene_type_id:
            query = self.queryTable(Gene).filter_by(annot_assembly_id=annot_assembly.id).\
                filter_by(start=start).filter_by(stop=stop).\
                filter_by(tax_id=tax_id).filter_by(chromosome=chromosome).\
                filter_by(entrezgene_type_id=entrezgene_type_id)
            if type_of_gene:
                query = query.filter_by(type_of_gene=type_of_gene)
            gene = query.first()
        else:
            gene = None
        if gene is None:
            gene = Gene(annot_assembly_id=annot_assembly_id, tax_id=tax_id, 
                chromosome=chromosome, locustag=locustag, strand=strand,
                start=start, stop=stop, type_of_gene=type_of_gene, 
                entrezgene_type_id=entrezgene_type_id,\
                ncbi_gene_id=ncbi_gene_id, gene_symbol=gene_symbol, 
                synonyms=synonyms, dbxrefs=dbxrefs, map_location=map_location,
                description=description, genomic_accession=genomic_accession, 
                genomic_version=genomic_version, genomic_gi=genomic_gi, \
                symbol_from_nomenclature_authority=symbol_from_nomenclature_authority,
                full_name_from_nomenclature_authority=full_name_from_nomenclature_authority,
                nomenclature_status=nomenclature_status,
                modification_date=modification_date)
            gene.annot_assembly = annot_assembly
            self.session.add(gene)
            self.session.flush()
        return gene
    
    def constructGeneSegments(self, box_ls=[], gene_commentary=None, gi=None,
        commentary_type=None):
        """
        2012.5.15
            add argument commentary_type to stop
                replicating gene_commentary.gene_commentary_type
        2012.4.26
        """
        gene_segment_ls = []
        if commentary_type:
            gene_commentary_type = self.getGeneCommentaryType(
                commentary_type=commentary_type)
        else:
            gene_commentary_type = gene_commentary.gene_commentary_type
        for box in box_ls:
            start, stop = box[:2]
            
            gene_segment = GeneSegment(start=start, stop=stop, gi=gi,
                gene_commentary_type=gene_commentary_type)
            gene_segment.gene_commentary = gene_commentary
            self.session.add(gene_segment)
            gene_segment_ls.append(gene_segment)
        self.session.flush()
        return gene_segment_ls
    
    def run(self):
        """
        """
        if self.debug:
            import pdb
            pdb.set_trace()
        
        if self.geneAnnotationPickleFname and self.tax_id:
            import pickle
            pickle_fname = os.path.expanduser(self.geneAnnotationPickleFname)
            #2012.3.26 stopped using '~/at_gene_model_pickelf',
            if os.path.isfile(pickle_fname):
                #2011-1-20 check if the pickled file already exists or not
                logging.error("File %s already exists, no gene model pickle output."\
                    %(pickle_fname))
                sys.exit(2)
            else:
                #2008-10-01	get gene model and pickle it into a file
                gene_id2model, chr_id2gene_id_ls, geneSpanRBDict = \
                    self.get_gene_id2model(tax_id=self.tax_id)
                gene_annotation = PassingData()
                gene_annotation.gene_id2model = gene_id2model
                gene_annotation.chr_id2gene_id_ls = chr_id2gene_id_ls
                gene_annotation.geneSpanRBDict = geneSpanRBDict
                picklef = open(os.path.expanduser(pickle_fname), 'w')
                pickle.dump(gene_annotation, picklef, -1)
                picklef.close()

def get_entrezgene_annotated_anchor(curs, tax_id,
    entrezgene_mapping_table='genome.gene',\
    annot_assembly_table='genome.annot_assembly'):
    """
    2011-1-25
        moved from annot.bin.common
        change e.gene_id to e.id
    2010-8-1 make sure chromosome is not null
    
    2008-08-12
        deal with genes with no strand info
        only start of the gene gets into chromosome2anchor_gene_tuple_ls. not stop.
    12-11-05
    
    12-17-05
        for TFBindingSiteParse.py
    """
    sys.stderr.write("Getting entrezgene_annotated_anchor ...")
    chromosome2anchor_gene_tuple_ls = {}
    gene_id2coord = {}
    curs.execute(f"select e.genomic_gi, a.chromosome, e.id, e.strand, "
        f"e.start, e.stop from {entrezgene_mapping_table} e, "
        f"{annot_assembly_table} a where e.genomic_gi = a.gi and "
        f"e.tax_id={tax_id} and a.chromosome is not null")
        #2010-8-1 make sure chromosome is not null
    
    #curs.execute("fetch 5000 from eaa_crs")
    rows = curs.fetchall()
    counter = 0#2010-8-1
    #while rows:
    for row in rows:
        genomic_gi, chromosome, gene_id, strand, start, stop = row
        gene_id = int(gene_id)
        if strand=='1' or strand=='+1' or strand=='+':
            strand = '+'
        elif strand=='-1' or strand=='-':
            strand = '-'
        else:
            strand = strand
        
        if chromosome not in chromosome2anchor_gene_tuple_ls:
            chromosome2anchor_gene_tuple_ls[chromosome] = []
        
        chromosome2anchor_gene_tuple_ls[chromosome].append((start, gene_id))
        gene_id2coord[gene_id] = (start, stop, strand, genomic_gi)
        #curs.execute("fetch 5000 from eaa_crs")
        #rows = curs.fetchall()
        counter += 1
    for chromosome in chromosome2anchor_gene_tuple_ls:	#sort the list
        chromosome2anchor_gene_tuple_ls[chromosome].sort()
    #curs.execute("close eaa_crs")
    sys.stderr.write("%s genes (including AS-isoforms). Done.\n"%counter)
    return chromosome2anchor_gene_tuple_ls, gene_id2coord

if __name__ == '__main__':
    import sys, os
    from palos import ProcessOptions
    main_class = GenomeDatabase
    po = ProcessOptions(sys.argv, main_class.option_default_dict, \
        error_doc=main_class.__doc__)
    
    instance = main_class(**po.long_option2value)
    instance.setup()
    instance.run()
    """
    if instance.debug:
        import pdb
        pdb.set_trace()
    
    import pickle
    #2011-1-20 check if the pickled file already exists or not
    pickle_fname = po.arguments[0]
    #2012.3.26 stopped using '~/at_gene_model_pickelf', instead use the 1st argument
    if os.path.isfile(os.path.expanduser(pickle_fname)):
        logging.error(f"File {pickle_fname} already exists, "
            "no gene model pickle output."
    else:
        #2008-10-01	get gene model and pickle it into a file
        gene_id2model, chr_id2gene_id_ls, geneSpanRBDict = \
            instance.get_gene_id2model(tax_id=3702)
        gene_annotation = PassingData()
        gene_annotation.gene_id2model = gene_id2model
        gene_annotation.chr_id2gene_id_ls = chr_id2gene_id_ls
        gene_annotation.geneSpanRBDict = geneSpanRBDict
        picklef = open(os.path.expanduser(pickle_fname), 'w')
        pickle.dump(gene_annotation, picklef, -1)
        picklef.close()
    """
    import sqlalchemy as sql
    #print dir(Gene)
    #print Gene.table.c.keys()
    #results = instance.session.query(Gene)
    #results = instance.session.execute(sql.select([Gene.table]))
    #print dir(results)
    #print dir(Gene.query)
    
    #import pdb
    #pdb.set_trace()
    i = 0
    block_size = 10
    rows = Gene.query.offset(i).limit(block_size)
    print(dir(rows))
    while rows.count()!=0:
        print(rows.count())
        for row in rows:
            i += 1
            print(row.id, row.ncbi_gene_id, row.gene_symbol)
        if i>=5*block_size:
            break
        rows = Gene.query.offset(i).limit(block_size)