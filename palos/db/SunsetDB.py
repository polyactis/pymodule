#!/usr/bin/env python3
"""
This is the ORM architecture of SunsetDB, build on top of sqlalchemy.

Examples:
	#setup database in postgresql
	%s -v postgresql -u yh -z pdc -d pmdb -k sunset
	
	#setup database in mysql
	%s -u yh -z papaya.usc.edu
	
"""
import sys, os, copy
__doc__ = __doc__%(sys.argv[0], sys.argv[0])

import hashlib
from sqlalchemy import Unicode, DateTime, String, BigInteger, Integer, UnicodeText, Text, Boolean, Float, Binary, Enum, Table
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import Column
from sqlalchemy import ForeignKey
from sqlalchemy.orm import relationship
from sqlalchemy import UniqueConstraint
from sqlalchemy import and_, or_, not_
from sqlalchemy.types import LargeBinary

from datetime import datetime
from palos import ProcessOptions, utils, PassingData
from palos import ngs
from palos.db import Database, TableClass, AbstractTableWithFilename
from palos.utils import runLocalCommand
from palos.utils import returnZeroFunc
from palos.ngs.io.VCFFile import VCFFile

Base = declarative_base()
#20190111 have to set it staticaly because SunsetDB is undefined at this point and it has to be defined after this.
_schemaname_ = "sunset"

	
user2group_table = Table('user2group',Base.metadata,
	Column('group_id',Integer,ForeignKey(_schemaname_+'.acl_group.id')),
	Column('user_id',Integer,ForeignKey(_schemaname_+'.acl_user.id'))
	)
	
group2phenotype_method_table = Table('group2phenotype_method',Base.metadata,
	Column('group_id',Integer,ForeignKey(_schemaname_+'.acl_group.id')),
	Column('phenotype_method_id',Integer,ForeignKey(_schemaname_+'.phenotype_method.id'))
	)

individual2group_table = Table('individual2group',Base.metadata,
	Column('group_id',Integer,ForeignKey(_schemaname_+'.acl_group.id')),
	Column('individual_id',Integer,ForeignKey(_schemaname_+'.individual.id'))
	)

individual2user_table = Table('individual2user',Base.metadata,
	Column('user_id',Integer,ForeignKey(_schemaname_+'.acl_user.id')),
	Column('individual_id',Integer,ForeignKey(_schemaname_+'.individual.id'))
	)

individual2batch_table = Table('individual2batch',Base.metadata,
	Column('individual_id',Integer,ForeignKey(_schemaname_+'.individual.id')),
	Column('sequence_batch_id',Integer,ForeignKey(_schemaname_+'.sequence_batch.id'))
	)

genotype_method2individual_alignment_table = Table('genotype_method2individual_alignment',Base.metadata,
	Column('individual_alignment_id',Integer,ForeignKey(_schemaname_+'.individual_alignment.id')),
	Column('genotype_method_id',Integer,ForeignKey(_schemaname_+'.genotype_method.id'))
	)

alignment_depth_interval_method2individual_alignment_table = Table('alignment_depth_interval_method2individual_alignment',Base.metadata,
	Column('alignment_depth_interval_method_id',Integer,ForeignKey(_schemaname_+'.alignment_depth_interval_method.id')),
	Column('individual_alignment_id',Integer,ForeignKey(_schemaname_+'.individual_alignment.id'))

	)

user2phenotype_method_table = Table('user2phenotype_method',Base.metadata,
	Column('phenotype_method_id',Integer,ForeignKey(_schemaname_+'.phenotype_method.id')),
	Column('user_id',Integer,ForeignKey(_schemaname_+'.acl_user.id'))
	)

class AnalysisMethod(Base, TableClass):
	"""
	2011-4-5
		record the analysis method used in like ScoreMethod or others.
	"""	
	__tablename__ = 'analysis_method'
	__table_args__ = {'schema':_schemaname_}
        
	id = Column(Integer,primary_key=True)
	short_name = Column(String(120))
	method_description = Column(String(8000))
	smaller_score_more_significant = Column(Integer)
	created_by = Column(String(200))
	updated_by = Column(String(200))
	date_created = Column(DateTime, default=datetime.now())
	date_updated = Column(DateTime)
	#using_options(tablename='analysis_method', metadata=__metadata__, session=__session__)
	#using_table_options(mysql_engine='InnoDB')
	
	score_method_list_in_analysis_method = relationship('ScoreMethod',back_populates='analysis_method',cascade='all,delete')

class README(Base, TableClass):
	__tablename__ = 'readme'
	__table_args__ = {'schema':_schemaname_}
        
	id = Column(Integer,primary_key=True)
	title = Column(String(2000))
	description = Column(Text)
	created_by = Column(String(128))
	updated_by = Column(String(128))
	date_created = Column(DateTime, default=datetime.now())
	date_updated = Column(DateTime)
	#using_options(tablename='readme', metadata=__metadata__, session=__session__)
	#using_table_options(mysql_engine='InnoDB')

class Family(Base, TableClass):
	__tablename__ = 'family'
	__table_args__ = {'schema':_schemaname_}
        
	id = Column(Integer,primary_key=True)
	short_name = Column(String(256))
	description = Column(Text)
	created_by = Column(String(128))
	updated_by = Column(String(128))
	date_created = Column(DateTime, default=datetime.now())
	date_updated = Column(DateTime)
	#using_options(tablename='family', metadata=__metadata__, session=__session__)
	#using_table_options(mysql_engine='InnoDB')
	
	individual_list_in_family = relationship('Individual',back_populates="family",cascade='all,delete')

class Country(Base):
	"""
	2011-3-1
	"""
	__tablename__ = 'country'
	__table_args__ = {'schema':_schemaname_}
	
	id = Column(Integer,primary_key=True)
	name = Column(String(100), unique=True)
	abbr = Column(String(10))
	capital = Column(Text)
	latitude = Column(Float)
	longitude = Column(Float)
	created_by = Column(String(128))
	updated_by = Column(String(128))
	date_created = Column(DateTime, default=datetime.now())
	date_updated = Column(DateTime)
	#using_options(tablename='country', metadata=__metadata__, session=__session__)
	#using_table_options(mysql_engine='InnoDB')
	
	site_list_in_country = relationship('Site',back_populates='country',cascade='all,delete')

class Site(Base, TableClass):
	"""
	20170411 add study_id
	2012.12.6 add unique constraint on (latitude, longitude)
	2011-4-29
		add column altitude
	2011-3-1
	"""
	__tablename__ = 'site'
	__table_args__ = {'schema':_schemaname_}
        
	id = Column(Integer,primary_key=True)
	short_name = Column(String(256))
	description = Column(Text)
	latitude = Column(Float)
	longitude = Column(Float)
	altitude = Column(Float)	#2011-4-29
	#study = ManyToOne('%s.Study'%(__name__), colname='study_id', ondelete='CASCADE', onupdate='CASCADE')	#2017.04.11
	study_id = Column(Integer,ForeignKey(_schemaname_+".study.id"))
	city = Column(String(100))
	stateprovince = Column(String(100))
	region = Column(String(100))
	zippostal = Column(String(20))
	#country = ManyToOne("%s.Country"%(__name__), colname='country_id', ondelete='CASCADE', onupdate='CASCADE')
	country_id = Column(Integer,ForeignKey(_schemaname_+".country.id"))
	created_by = Column(String(128))
	updated_by = Column(String(128))
	date_created = Column(DateTime, default=datetime.now())
	date_updated = Column(DateTime)
	#using_options(tablename='site', metadata=__metadata__, session=__session__)
	#using_table_options(mysql_engine='InnoDB')
	#using_table_options(UniqueConstraint('short_name', 'latitude', 'longitude', 'city', 'stateprovince', 'country_id'))
	#using_table_options(UniqueConstraint('short_name', 'latitude', 'longitude'))
	UniqueConstraint('short_name', 'latitude', 'longitude', 'city', 'stateprovince', 'country_id',name='site_sllcsc')
	UniqueConstraint('short_name', 'latitude', 'longitude',name='site_sll')
        
	study = relationship('Study',back_populates='site_list_in_study')
	country = relationship('Country',back_populates='site_list_in_country')
	individual_list_in_site = relationship('Individual',back_populates="site",cascade='all,delete')

class Group(Base):
	"""
	2011-3-1
	"""
	__tablename__ = 'acl_group'
	__table_args__ = {'schema':_schemaname_}
        
	id = Column(Integer,primary_key=True)
	#name = Column(Unicode(512),required=True)
	name = Column(Unicode(512))
	#user_ls = ManyToMany('%s.User'%(__name__),tablename='user2group', local_colname='group_id')
	#phenotype_method_ls = ManyToMany("%s.PhenotypeMethod"%(__name__), tablename='group2phenotype_method', local_colname='group_id')
	#individual_ls = ManyToMany("%s.Individual"%(__name__),tablename='individual2group', local_colname='group_id')
	#using_table_options(mysql_engine='InnoDB', useexisting=True)
	#using_options(tablename='acl_group', metadata=__metadata__, session=__session__)
	#group is preserved keyword in postgresql (mysql likely)
	
	user_ls = relationship('User',secondary=user2group_table,back_populates='group_ls')
	phenotype_method_ls = relationship('PhenotypeMethod',secondary=group2phenotype_method_table,back_populates='group_ls')
	individual_ls = relationship('Individual',secondary=individual2group_table,back_populates='group_ls')
	
	def __repr__(self):
		return (u'<Group: name=%s>' % self.name).encode('utf-8')

class User(Base):
	"""
	2011-3-1
	"""
	__tablename__ = 'acl_user'
	__table_args__ = {'schema':_schemaname_}
        
	id = Column(Integer,primary_key=True)
	title = Column(String(4))
	realname = Column(Unicode(512))
	email = Column(String(100))
	username = Column(String(10))
	#_password = Column(String(40), colname='password', synonym='password')
	_password = Column(String(40))
	organisation = Column(Unicode(100))
	#isAdmin = Column(postgresql.Enum(('Y','N'), name=is_admin_enum_type), default='N', required=True,)
	#isAdmin = Column(Enum("Y","N", name="is_admin_enum_type"), default='N', required=True,)
	isAdmin = Column(Enum("Y","N", name="is_admin_enum_type"), default='N')
	#group_ls = ManyToMany('%s.Group'%(__name__), tablename='user2group', local_colname='user_id')
	#phenotype_method_ls = ManyToMany("%s.PhenotypeMethod"%(__name__),tablename='user2phenotype_method', local_colname='user_id')
	#individual_ls = ManyToMany("%s.Individual"%(__name__),tablename='individual2user', local_colname='user_id')
	created_by = Column(String(128))
	updated_by = Column(String(128))
	date_created = Column(DateTime, default=datetime.now())
	date_updated = Column(DateTime)
	#using_options(tablename='acl_user', metadata=__metadata__, session=__session__)
	#user is preserved keyword in postgresql (mysql likely)
	#using_table_options(mysql_engine='InnoDB')
	#using_table_options(UniqueConstraint('realname', 'username'))
	UniqueConstraint('realname', 'username',name='acl_user_realname_username')
	
	individual_list_in_acl_user = relationship('Individual',back_populates="collector",cascade='all,delete')
	phenotype_method_list_in_user = relationship('PhenotypeMethod',back_populates='collector',cascade='all,delete')
	group_ls = relationship('Group',secondary=user2group_table,back_populates='user_ls')
	phenotype_method_ls = relationship('PhenotypeMethod',secondary=user2phenotype_method_table,back_populates='user_ls')
	individual_ls = relationship('Individual',secondary=individual2user_table,back_populates='user_ls')
	
	def validate_password(self, password):
		"""
        Check the password against existing credentials.
        
        :param password: the password that was provided by the user to
            try and authenticate. This is the clear text version that we will
            need to match against the hashed one in the database.
        :type password: unicode object.
        :return: Whether the password is valid.
        :rtype: bool
        
        """
		hashed_pass = hashlib.sha1()
		hashed_pass.update(self._password[:8] + password)
		return self._password[8:] == hashed_pass.hexdigest()

	def _set_password(self, password):
		"""
		encrypts password on the fly using the encryption
		algo defined in the configuration
		"""
		self._password = self.__encrypt_password(password)

	def _get_password(self):
		"""
		returns password
		"""
		return self._password

	password = property(_get_password,_set_password)

	@classmethod
	def __encrypt_password(cls, password):
		"""
		Hash the given password with the specified algorithm. Valid values
		for algorithm are 'md5' and 'sha1'. All other algorithm values will
		be essentially a no-op.
		"""
		hashed_password = password
		
		if isinstance(password, str):
			password_8bit = password.encode('UTF-8')
		else:
			password_8bit = password
		
		salt = hashlib.sha1()
		salt.update(os.urandom(60))
		salt_text = salt.hexdigest()
		hash = hashlib.sha1()
		hash.update(salt_text[:8] + password_8bit)
		hashed_password = salt_text[:8] + hash.hexdigest()
		print(f"{'*'*20} {salt_text[:8]} {hashed_password[:8]}")
		
		if not isinstance(hashed_password, str):
			hashed_password = hashed_password.decode('UTF-8')

		return hashed_password

class Individual(Base, TableClass):
	"""
	2013.3.15 renamed is_contaminated to outdated_index
	2013.3.13 added study column, 
	2012.12.6 get rid of latitude, longitude, altitude, it's now rolled into site.
	2012.9.27 add is_contaminated
	2012.7.5 added column sequence_batch_ls	
	2012.6.19 add column target_coverage
	2012.2.7
		make code unique
	2011-9-8
		add column comment
	2011-4-8
		change (firstname, lastname) to one field, name.
	2011-3-1
		add tax_id, collector, site, group_ls, user_ls, latitude, longitude
	"""
	__tablename__ = 'individual'
	__table_args__ = {'schema':_schemaname_}
	
	id = Column(Integer,primary_key=True)
	#family = ManyToOne('%s.Family'%(__name__), colname='family_id', ondelete='CASCADE', onupdate='CASCADE')
	family_id = Column(Integer,ForeignKey(_schemaname_+'.family.id'))
	code = Column(String(256), unique=True)
	name = Column(String(256))
	sex = Column(String(256))
	birthdate = Column(DateTime)
	birthplace = Column(String(256))
	access = Column(Enum("public", "restricted", name="access_enum_type"), default='restricted')
	tax_id =Column(Integer)	#2011-3-1
	age = Column(Integer)	#2011-4-28
	collection_date = Column(DateTime)	#2011-4-28
	#collector = ManyToOne("%s.User"%(__name__), colname='collector_id', ondelete='CASCADE', onupdate='CASCADE')	#2011-3-1
	collector_id = Column(Integer,ForeignKey(_schemaname_+".acl_user.id"))
	#site = ManyToOne("%s.Site"%(__name__), colname='site_id', ondelete='CASCADE', onupdate='CASCADE')	#2011-3-1
	site_id = Column(Integer,ForeignKey(_schemaname_+".site.id"))
	#group_ls = ManyToMany('%s.Group'%(__name__),tablename='individual2group', local_colname='individual_id', remote_colname='group_id')	#2011-3-1
	#user_ls = ManyToMany('%s.User'%(__name__), tablename='individual2user', local_colname='individual_id', remote_colname='user_id')	#2011-3-1
	comment = Column(String(4096))
	target_coverage = Column(Integer)	#2012.6.19
	outdated_index = Column(Integer, default=0)	#2012.9.27 is_contaminated, 2013.3.15 changed. any non-zero means outdated. 
	#sequence_batch_ls = ManyToMany('%s.SequenceBatch'%(__name__), tablename='individual2batch', local_colname='individual_id')	#2012.7.5
	#study = ManyToOne('%s.Study'%(__name__), colname='study_id', ondelete='CASCADE', onupdate='CASCADE')	#2013.03.13
	study_id = Column(Integer,ForeignKey(_schemaname_+".study.id"))
	#ManyToOne('SequenceBatch', colname='sequence_batch_id', ondelete='CASCADE', onupdate='CASCADE')	#2012.7.5
	created_by = Column(String(128))
	updated_by = Column(String(128))
	date_created = Column(DateTime, default=datetime.now())
	date_updated = Column(DateTime)
	#using_options(tablename='individual', metadata=__metadata__, session=__session__)
	#using_table_options(mysql_engine='InnoDB')
	#using_table_options(UniqueConstraint('family_id', 'code', 'tax_id'))
	UniqueConstraint('family_id', 'code', 'tax_id',name='individual_family_id_code_tax_id')
	
	family = relationship('Family',back_populates="individual_list_in_family")
	collector = relationship('User',back_populates="individual_list_in_acl_user")
	site = relationship('Site',back_populates="individual_list_in_site")
	study = relationship('Study',back_populates="individual_list_in_study")
	ind_seq_list_in_individual = relationship('IndividualSequence',back_populates='individual',cascade='all,delete')
	ind2ind_list_in_individual1 = relationship('Ind2Ind',back_populates='individual1',foreign_keys='Ind2Ind.individual1_id')
	ind2ind_list_in_individual2 = relationship('Ind2Ind',back_populates='individual2',foreign_keys='Ind2Ind.individual2_id')
	genotype_list_in_individual = relationship('Genotype',back_populates='individual',cascade='all,delete')
	phenotype_list_in_individual = relationship('Phenotype',back_populates='individual',cascade='all,delete')
	group_ls = relationship('Group',secondary=individual2group_table, back_populates='individual_ls')
	user_ls = relationship('User',secondary=individual2user_table, back_populates='individual_ls')
	sequence_batch_ls = relationship('SequenceBatch', secondary=individual2batch_table, back_populates='individual_ls')
	
	@classmethod
	def getIndividualsForACL(cls, user = None, db=None):
		"""
		2011-3-1
			get all individuals that could be accessed by this user
		"""
		TableClass = Individual
		query = db.queryTable(Individual)
		clause = and_(TableClass.access == 'public')
		if user is not None:
			if user.isAdmin == 'Y':
				return query
			clause = or_(clause,TableClass.collector == user, TableClass.user_ls.any(User.id == user.id),
						TableClass.group_ls.any(Group.id.in_([group.id for group in user.group_ls])))
		query = query.filter(clause)
		return query
	
	def checkACL(self,user):
		"""
		2011-3-1
			check if the user could access this individual
		"""
		if self.access == 'public':
			return True
		if user is None:
			return False
		if user.isAdmin == 'Y':
			return True
		if self.collector == user: 
			return True
		if user in self.user_ls:
			return True
		if [group in self.group_ls for group in user.group_ls]: 
			return True
		return False
	
	def codeSexInNumber(self):
		"""
		2012.9.4 if sex is empty or None, return 0
		2011.12.4
			represent male as 1. represent female as 2.
		"""
		if self.sex:
			if self.sex[0]=='M':
				return 1
			else:
				return 2
		else:
			return 0
	
	def getCurrentAge(self):
		"""
		2012.2.22
			get the most current age of the monkey
		"""
		if self.birthdate:
			from datetime import datetime
			return datetime.now().year - self.birthdate.year
		elif self.collection_date and self.age:
			extraYears = datetime.now().year - self.collection_date.year
			if self.age:
				ageInYears = self.age
			else:
				ageInYears = 0
			return ageInYears + extraYears
		else:
			return None

class Ind2Ind(Base, TableClass):
	"""
	2013.3.13 added study column
	2011-5-5
		add a unique constraint
	"""
	__tablename__ = 'ind2ind'
	__table_args__ = {'schema':_schemaname_}
	
	id = Column(Integer,primary_key=True)
	#individual1 = ManyToOne('%s.Individual'%(__name__), colname='individual1_id', ondelete='CASCADE', onupdate='CASCADE')
	individual1_id = Column(Integer,ForeignKey(_schemaname_+'.individual.id'))
	#individual2 = ManyToOne('%s.Individual'%(__name__), colname='individual2_id', ondelete='CASCADE', onupdate='CASCADE')
	individual2_id = Column(Integer,ForeignKey(_schemaname_+'.individual.id'))
	#relationship_type = ManyToOne('%s.RelationshipType'%(__name__), colname='relationship_type_id', ondelete='CASCADE', onupdate='CASCADE')
	relationship_type_id = Column(Integer,ForeignKey(_schemaname_+".relationship_type.id"))
	#study = ManyToOne('%s.Study'%(__name__), colname='study_id', ondelete='CASCADE', onupdate='CASCADE')	#2013.03.13
	study_id = Column(Integer,ForeignKey(_schemaname_+'.study.id'))
	created_by = Column(String(128))
	updated_by = Column(String(128))
	date_created = Column(DateTime, default=datetime.now())
	date_updated = Column(DateTime)
	#using_options(tablename='ind2ind')
	#using_table_options(mysql_engine='InnoDB')
	#using_table_options(UniqueConstraint('individual1_id', 'individual2_id', 'relationship_type_id', 'study_id'))

	UniqueConstraint('individual1_id', 'individual2_id', 'relationship_type_id', 'study_id',name='ind2ind_iirs')
	
	relationship_type = relationship('RelationshipType',back_populates='ind2ind_list_in_relationship_type')
	study = relationship('Study', back_populates='ind2ind_list_in_study')
	individual1 = relationship('Individual',back_populates='ind2ind_list_in_individual1',foreign_keys='Ind2Ind.individual1_id')
	individual2 = relationship('Individual',back_populates='ind2ind_list_in_individual2',foreign_keys='Ind2Ind.individual2_id')

class RelationshipType(Base, TableClass):
	__tablename__ = 'relationship_type'
	__table_args__ = {'schema':_schemaname_}
	
	id = Column(Integer,primary_key=True)
	short_name = Column(String(256), unique=True)
	description = Column(Text)
	created_by = Column(String(128))
	updated_by = Column(String(128))
	date_created = Column(DateTime, default=datetime.now())
	date_updated = Column(DateTime)
	#using_options(tablename='relationship_type', metadata=__metadata__, session=__session__)
	#using_table_options(mysql_engine='InnoDB')

	ind2ind_list_in_relationship_type = relationship('Ind2Ind',back_populates='relationship_type',cascade='all,delete')

class AlignmentMethod(Base):
	"""
	2011-9-1
		add column command
	2011-3-3
	"""
	__tablename__ = 'alignment_method'
	__table_args__ = {'schema':_schemaname_}
	
	id = Column(Integer,primary_key=True)
	short_name = Column(String(256), unique=True)
	command = Column(String(256))	#sub-command of bwa
	description = Column(Text)
	#individual_alignment_ls = OneToMany("%s.IndividualAlignment"%(__name__))
	created_by = Column(String(128))
	updated_by = Column(String(128))
	date_created = Column(DateTime, default=datetime.now())
	date_updated = Column(DateTime)
	#using_options(tablename='alignment_method', metadata=__metadata__, session=__session__)
	#using_table_options(mysql_engine='InnoDB')
	
	individual_alignment_ls = relationship('IndividualAlignment',back_populates='alignment_method',cascade='all,delete')

class IndividualAlignment(Base, AbstractTableWithFilename):
	"""
	2013.08.08 added column path_to_depth_file, depth_file_size
	2013.04.11 added column reduce_reads
	2013.04.01 added columns local_realigned, version
	2012.9.21 rename IndividualAlignment.ind_sequence to individual_sequence
		add file_size
	#2012.9.19 to distinguish alignments from different libraries/lanes/batches
		add individual_sequence_file_raw
	2012.7.26 add column parent_individual_alignment, mask_genotype_method
	2012.7.14 add md5sum
	2012.6.13
		add column outdated_index to accommodate old (incomplete or bad alignments) 
	2012.4.2
		add a couple of statistics, perc_reads_mapped, ... perc_mapped_to_diff_chrs.
	2011-11-28
		add pass_qc_read_base_count which counts the number of bases in all pass-QC reads (base quality>=20, read mapping quality >=30).
	2011-9-19
		add mean_depth, read_group_added
	2011-8-3
		add column median_depth, mode_depth
	2011-3-3
	"""
	__tablename__ = 'individual_alignment'
	__table_args__ = {'schema':_schemaname_}
	
	id = Column(Integer,primary_key=True)
	#individual_sequence = ManyToOne('%s.IndividualSequence'%__name__, colname='ind_seq_id', ondelete='CASCADE', onupdate='CASCADE')
	ind_seq_id = Column(Integer,ForeignKey(_schemaname_+'.individual_sequence.id'))
	#ref_sequence = ManyToOne('%s.IndividualSequence'%__name__, colname='ref_ind_seq_id', ondelete='CASCADE', onupdate='CASCADE')
	ref_ind_seq_id = Column(Integer,ForeignKey(_schemaname_+'.individual_sequence.id'))
	#alignment_method = ManyToOne('%s.AlignmentMethod'%__name__, colname='alignment_method_id', ondelete='CASCADE', onupdate='CASCADE')
	alignment_method_id = Column(Integer,ForeignKey(_schemaname_+'.alignment_method.id'))
	#genotype_method_ls = ManyToMany("%s.GenotypeMethod"%__name__,tablename='genotype_method2individual_alignment', local_colname='individual_alignment_id')
	path = Column(Text)
	path_to_depth_file = Column(Text)	#2013.08.08
	depth_file_size = Column(BigInteger)	#2013.08.08
	format = Column(String(512))
	median_depth = Column(Float)	#2011-8-2
	mode_depth = Column(Float)	#2011-8-2
	mean_depth = Column(Float)	#2011-9-12
	pass_qc_read_base_count = Column(BigInteger)	#2011-11-28	QC = (base quality>=20, read mapping quality >=30). 
	read_group_added = Column(Integer, default=0)	# 2011-9-15 0=No, 1=Yes
	perc_reads_mapped = Column(Float)	#2012.4.2
	perc_secondary = Column(Float)	#20170603
	perc_supplementary = Column(Float)	#20170603
	perc_duplicates = Column(Float)	#2012.4.2
	perc_paired = Column(Float)	#2012.4.2
	perc_properly_paired = Column(Float)	#2012.4.2
	perc_both_mates_mapped = Column(Float)	#2012.4.2
	perc_singletons = Column(Float)	#2012.4.2
	perc_mapped_to_diff_chrs = Column(Float)	#2012.4.2
	perc_mapq5_mapped_to_diff_chrs = Column(Float)	#2012.4.2
	total_no_of_reads = Column(BigInteger)	#2012.4.2
	local_realigned = Column(Integer, default=0)	#2013.04.01
	reduce_reads = Column(Integer, default=0)	#2013.04.11
	outdated_index = Column(Integer, default=0)	#2012.6.13 any non-zero means outdated. to allow multiple outdated alignments
	md5sum = Column(Text, unique=True)
	file_size = Column(BigInteger)	#2012.9.21
	read_group = Column(Text)	#2013.04.08 record read_group here so that if getReadGroup() changes. it'll be fine.
	#2012.7.26 the parent individual_alignment
	#parent_individual_alignment = ManyToOne('%s.IndividualAlignment'%__name__, colname='parent_individual_alignment_id', ondelete='CASCADE', onupdate='CASCADE')
	parent_individual_alignment_id = Column(Integer,ForeignKey(_schemaname_+'.individual_alignment.id', 
											ondelete='CASCADE',onupdate='CASCADE'))
	#2012.7.26 mask loci of the alignment out for read-recalibration 
	#mask_genotype_method = ManyToOne('%s.GenotypeMethod'%__name__, colname='mask_genotype_method_id', ondelete='CASCADE', onupdate='CASCADE')
	mask_genotype_method_id = Column(Integer,ForeignKey(_schemaname_+'.genotype_method.id'))
	#2012.9.19 to distinguish alignments from different libraries/lanes/batches
	#individual_sequence_file_raw = ManyToOne('%s.IndividualSequenceFileRaw'%__name__, colname='individual_sequence_file_raw_id', \
	#								ondelete='CASCADE', onupdate='CASCADE')
	individual_sequence_file_raw_id = Column(Integer,ForeignKey(_schemaname_+'.individual_sequence_file_raw.id'))
	created_by = Column(String(128))
	updated_by = Column(String(128))
	date_created = Column(DateTime, default=datetime.now())
	date_updated = Column(DateTime)
	#using_options(tablename='individual_alignment', metadata=__metadata__, session=__session__)
	#using_table_options(mysql_engine='InnoDB')
	#using_table_options(UniqueConstraint('ind_seq_id', 'ref_ind_seq_id', 'alignment_method_id', 'outdated_index', \
	#							'parent_individual_alignment_id', 'mask_genotype_method_id',\
	#							'individual_sequence_file_raw_id', 'local_realigned', 'reduce_reads'))
	UniqueConstraint('ind_seq_id', 'ref_ind_seq_id', 'alignment_method_id', 'outdated_index', \
					'parent_individual_alignment_id', 'mask_genotype_method_id',\
					'individual_sequence_file_raw_id', 'local_realigned', 'reduce_reads',\
					name='individual_alignment_iraopmilr')
	
	individual_sequence = relationship('IndividualSequence',back_populates='individual_alignment_list_in_individual_sequence',
									foreign_keys='IndividualAlignment.ind_seq_id')
	ref_sequence = relationship('IndividualSequence',back_populates='individual_alignment_list_with_this_reference',
							foreign_keys='IndividualAlignment.ref_ind_seq_id')
	alignment_method = relationship('AlignmentMethod',back_populates='individual_alignment_ls')
	mask_genotype_method = relationship('GenotypeMethod',back_populates='individual_alignment_list_in_genotype_method')
	individual_sequence_file_raw = relationship('IndividualSequenceFileRaw',back_populates='ind_alig_list_in_ind_seq_file_raw')
	individual_alignment_consensus_seq_list = relationship('IndividualAlignmentConsensusSequence',
											back_populates='individual_alignment',cascade='all,delete')
	parent_individual_alignment = relationship('IndividualAlignment',
											foreign_keys='IndividualAlignment.parent_individual_alignment_id')
	#child_individual_alignment = relationship('IndividualAlignment', back_populates='parent_individual_alignment')
	genotype_method_ls = relationship('GenotypeMethod',secondary=genotype_method2individual_alignment_table, 
									back_populates='individual_alignment_ls')
	alignment_depth_interval_method_ls = relationship('AlignmentDepthIntervalMethod', 
										secondary=alignment_depth_interval_method2individual_alignment_table,
										back_populates="individual_alignment_ls")
	def getReadGroup(self):
		"""
		2013.04.08 revert it back to its old form (don't change it EVER) as alignment files has it embedded.
		2013.04.01 adding alignment_method_id
		2013.3.18 bugfix: sequencer is now a db entry.
		2012.9.19
			add three more optional additions to the read group
		2011-11-02
			read group is this alignment's unique identifier in a sam/bam file.
		"""
		if self.read_group:
			return self.read_group
		
		sequencer = self.individual_sequence.sequencer
		read_group = '%s_%s_%s_%s_vs_%s'%(self.id, self.ind_seq_id, \
							self.individual_sequence.individual.code, \
							sequencer.short_name, self.ref_ind_seq_id)
		return read_group

	def constructBaseNamePrefix(self):
		"""
		2013.04.08 for constructRelativePath()
		"""
		read_group = self.getReadGroup()
		
		read_group = '%s_by_method%s_realigned%s_reduced%s'%(read_group,
						self.alignment_method_id, self.local_realigned, self.reduce_reads)
		if self.parent_individual_alignment_id:
			read_group += '_p%s'%(self.parent_individual_alignment_id)
		if self.mask_genotype_method_id:
			read_group += '_m%s'%(self.mask_genotype_method_id)
		if self.individual_sequence_file_raw_id:
			read_group += '_r%s'%(self.individual_sequence_file_raw_id)
		return read_group

	def getCompositeID(self):
		"""
		2012.1.25
			ID to be used to identify members of trios. almost same as getReadGroup()
		"""
		compositeID = '%s_%s_%s_%s'%(self.id, self.ind_seq_id, self.individual_sequence.individual.id, \
								self.individual_sequence.individual.code)
		return compositeID
	
	folderName = 'individual_alignment'
	def constructRelativePath(self, subFolder='individual_alignment', **keywords):
		"""
		"""
		if subFolder is None:
			subFolder = self.folderName
		#'/' must not be put in front of the relative path.
		# otherwise, os.path.join(self.data_dir, dst_relative_path) will only take the path of dst_relative_path.
		baseNamePrefix = self.constructBaseNamePrefix()
		pathPrefix = '%s/%s'%(subFolder, baseNamePrefix)
		dst_relative_path = '%s.%s'%(pathPrefix, self.format)
		
		return dst_relative_path
	
	def constructPSMCPlotLabel(self, **keywords):
		"""
		2013.2.16 moved from PSMCOnAlignmentWorkflow.py
		2013.2.12
		"""
		return '%s_%s%dMD'%(self.individual_sequence.individual.code, self.individual_sequence.individual.sex, self.median_depth)
	
	
class IndividualAlignmentConsensusSequence(Base, AbstractTableWithFilename):
	"""
	2013.2.5 input for psmc, extracted from individual alignment.
		The extraction takes so long. so use this to store them.
	"""
	__tablename__ = 'individual_alignment_consensus_sequence'
	__table_args__ = {'schema':_schemaname_}
	
	id = Column(Integer, primary_key=True)
	#individual_alignment = ManyToOne('%s.IndividualAlignment'%__name__, colname='individual_alignment_id', \
	#								ondelete='CASCADE', onupdate='CASCADE')
	individual_alignment_id = Column(Integer,ForeignKey(_schemaname_+'.individual_alignment.id'))
	path = Column(Text, unique=True)
	format = Column(String(512))
	minDP = Column(Integer)
	maxDP = Column(Integer)
	minBaseQ = Column(Integer, default=20)
	minMapQ = Column(Integer, default=30)
	minRMSMapQ = Column(Integer, default=10)	#root mean squared mapping quality of reads covering the locus
	minDistanceToIndel = Column(Integer, default=5)
	
	no_of_chromosomes = Column(Integer)
	no_of_bases = Column(BigInteger)
	md5sum = Column(Text, unique=True)
	file_size = Column(BigInteger)	#2012.9.21
	original_path = Column(Text)
	comment = Column(Text)
	created_by = Column(String(128))
	updated_by = Column(String(128))
	date_created = Column(DateTime, default=datetime.now())
	date_updated = Column(DateTime)
	#using_options(tablename='individual_alignment_consensus_sequence', metadata=__metadata__, session=__session__)
	#using_table_options(mysql_engine='InnoDB')
	#using_table_options(UniqueConstraint('individual_alignment_id', 'minDP', 'maxDP', 'minBaseQ',\
	#									'minMapQ', 'minRMSMapQ', 'minDistanceToIndel', 'no_of_chromosomes'))
	
	UniqueConstraint('individual_alignment_id', 'minDP', 'maxDP', 'minBaseQ', \
					'minMapQ', 'minRMSMapQ', 'minDistanceToIndel', 'no_of_chromosomes', \
					name='individual_alignment_consensus_sequence_immmmmmn')
	
	individual_alignment = relationship('IndividualAlignment', back_populates='individual_alignment_consensus_seq_list')
	
	folderName = 'individual_alignment_consensus_sequence'
	
	def constructRelativePath(self, data_dir=None, subFolder=None, sourceFilename=None, **keywords):
		"""
		2013.2.5
		"""
		if not subFolder:
			subFolder = self.folderName
		outputDirRelativePath = subFolder
		if data_dir and outputDirRelativePath:
			#make sure the final output folder is created. 
			outputDirAbsPath = os.path.join(data_dir, outputDirRelativePath)
			if not os.path.isdir(outputDirAbsPath):
				os.makedirs(outputDirAbsPath)
		
		filename_part_ls = []
		if self.id:
			filename_part_ls.append(self.id)
		if self.individual_alignment_id:
			filename_part_ls.append('aln%s'%(self.individual_alignment_id))
		if self.minDP is not None:
			filename_part_ls.append('minDP%s'%(self.minDP))
		if self.maxDP is not None:
			filename_part_ls.append("maxDP%s"%(self.maxDP))
		if self.no_of_chromosomes is not None:
			filename_part_ls.append("%sChromosomes"%(self.no_of_chromosomes))
		
		filename_part_ls = map(str, filename_part_ls)
		fileRelativePath = os.path.join(outputDirRelativePath, '%s.fastq.gz'%('_'.join(filename_part_ls)))
		return fileRelativePath
	
class IndividualSequence(Base, AbstractTableWithFilename):
	"""
	2013.3.13 added/changed sequence_batch, sequencer, sequence_type, no_of_chromosomes
	2012.2.27
		add column read_count
	2012.1.26
		add column individual_sequence_file_raw_ls
	2011-8-30
		add column chromosome
	2011-8-18
		add field original_path, quality_score_format, parent_individual_sequence, filtered
	2011-8-5
		change type of base_count to BigInteger
	2011-8-3
		add column base_count
	2011-5-8
		add column coverage
	2011-3-3
	"""
	__tablename__ = 'individual_sequence'
	__table_args__ = {'schema':_schemaname_}
	
	id = Column(Integer, primary_key=True)
	#individual = ManyToOne('%s.Individual'%(__name__), colname='individual_id', ondelete='CASCADE', onupdate='CASCADE')
	individual_id = Column(Integer, ForeignKey(_schemaname_ + '.individual.id'))
	#sequencer = ManyToOne('%s.Sequencer'%(__name__), colname='sequencer_id', ondelete='CASCADE', onupdate='CASCADE')
	sequencer_id = Column(Integer, ForeignKey(_schemaname_ + '.sequencer.id'))
		# 454, GA, Sanger
	#sequence_type = ManyToOne('%s.SequenceType'%(__name__), colname='sequence_type_id', ondelete='CASCADE', onupdate='CASCADE')
	sequence_type_id = Column(Integer, ForeignKey(_schemaname_ + '.sequence_type.id'))
		#genome, contig, SR (single-end read) or PE ...
	no_of_chromosomes = Column(Integer)	#1,2,4,5,X,Y,etc
	#tissue  = ManyToOne('%s.Tissue'%(__name__), colname='tissue_id', ondelete='CASCADE', onupdate='CASCADE')	#2011-5-9
	tissue_id = Column(Integer, ForeignKey(_schemaname_ + '.tissue.id'))
	coverage = Column(Float)	#2011-5-8
	read_count = Column(BigInteger)	#2012.2.27
	base_count = Column(BigInteger)	#2011-8-2
	path = Column(Text, unique=True)	#storage folder path
	format = Column(String(512))	#fasta, fastq
	original_path = Column(Text)	#the path to the original file
	quality_score_format = Column(String(512))	#Standard=Phred+33 (=Sanger), Illumina=Phred+64 (roughly, check pymodule/utils for exact formula)
	# Illumina1.8+ (after 2011-02) is Standard.
	#parent_individual_sequence = ManyToOne('%s.IndividualSequence'%(__name__), colname='parent_individual_sequence_id', ondelete='SET NULL', onupdate='CASCADE')
	parent_individual_sequence_id = Column(Integer, ForeignKey(_schemaname_ + '.individual_sequence.id', 
															ondelete='SET NULL', onupdate='CASCADE'))
	filtered = Column(Integer, default=0)	#0 means not. 1 means yes.
	#individual_sequence_file_ls = OneToMany("%s.IndividualSequenceFile"%(__name__))
	#individual_sequence_file_raw_ls = OneToMany("%s.IndividualSequenceFileRaw"%(__name__))
	#sequence_batch = ManyToOne('%s.SequenceBatch'%(__name__), colname='sequence_batch_id', ondelete='CASCADE', onupdate='CASCADE')	#2013.03.13
	sequence_batch_id = Column(Integer, ForeignKey(_schemaname_ + '.sequence_batch.id'))
	is_contaminated = Column(Integer, default=0)	#2013.3.15 field to mark whether it's contaminated or not.
	outdated_index = Column(Integer, default=0)	#2013.3.15 any non-zero means outdated. to allow multiple outdated alignments
	version = Column(Integer, default=1)
	comment = Column(Text)
	created_by = Column(String(128))
	updated_by = Column(String(128))
	date_created = Column(DateTime, default=datetime.now())
	date_updated = Column(DateTime)
	#using_options(tablename='individual_sequence', metadata=__metadata__, session=__session__)
	#using_table_options(mysql_engine='InnoDB')
	#using_table_options(UniqueConstraint('individual_id', 'sequencer_id', 'sequence_type_id', 'tissue_id',\
	#	'filtered','no_of_chromosomes','format', 'parent_individual_sequence_id', 'sequence_batch_id', 'version', \
	#	'is_contaminated', 'outdated_index'))
	UniqueConstraint('individual_id', 'sequencer_id', 'sequence_type_id', 'tissue_id', \
					'filtered', 'no_of_chromosomes', 'format', 'parent_individual_sequence_id', 'sequence_batch_id', 'version', \
					'is_contaminated', 'outdated_index', name='individual_sequence_isstfnfpsvio')
	
	individual_alignment_list_in_individual_sequence = relationship('IndividualAlignment', back_populates='individual_sequence', 
													cascade='all,delete', foreign_keys='IndividualAlignment.ind_seq_id')
	individual_alignment_list_with_this_reference = relationship('IndividualAlignment', back_populates='ref_sequence', 
													cascade='all,delete', foreign_keys='IndividualAlignment.ref_ind_seq_id')
	individual = relationship('Individual', back_populates='ind_seq_list_in_individual')
	sequencer = relationship('Sequencer', back_populates='ind_seq_in_sequencer')
	sequence_type = relationship('SequenceType', back_populates='ind_seq_in_sequence_type')
	tissue = relationship('Tissue', back_populates='ind_seq_in_tissue')
	parent_individual_sequence = relationship('IndividualSequence', foreign_keys='IndividualSequence.parent_individual_sequence_id')
	#child_individual_sequence = relationship('IndividualSequence', back_populates='parent_individual_sequence')
	
	individual_sequence_file_ls = relationship('IndividualSequenceFile', back_populates='individual_sequence', 
											cascade='all,delete')
	individual_sequence_file_raw_ls = relationship('IndividualSequenceFileRaw', back_populates='individual_sequence', 
												cascade='all,delete')
	sequence_batch = relationship('SequenceBatch', back_populates='ind_seq_list_in_sequence_batch')
	ind_seq2_list_in_ind1_seq = relationship('IndividualSequence2Sequence', back_populates='individual1_sequence', 
											cascade='all,delete', foreign_keys='IndividualSequence2Sequence.individual1_sequence_id')
	ind_seq2_list_in_ind2_seq = relationship('IndividualSequence2Sequence', back_populates='individual2_sequence', 
											cascade='all,delete', foreign_keys='IndividualSequence2Sequence.individual2_sequence_id')
	locus_list_in_ref_sequence = relationship('Locus', back_populates='ref_sequence', cascade='all,delete')
	genotype_method_list_in_ind_seq = relationship('GenotypeMethod', back_populates='ref_sequence', cascade='all,delete')
	alignemnt_depth_interval_method_list_in_ind_seq = relationship('AlignmentDepthIntervalMethod', back_populates='ref_sequence', 
																cascade='all,delete')
	
	def constructRelativePath(self, subFolder='individual_sequence', **keywords):
		"""
		2012.7.13 link to constructRelativePathForIndividualSequence()
		"""
		return self.constructRelativePathForIndividualSequence(subFolder=subFolder)
	
	folderName='individual_sequence'
	def constructRelativePathForIndividualSequence(self, subFolder=None):
		"""
		2013.2.16 use filename_part_ls
		2012.2.10
			add "split" in the end of the path
		2011-8-3
			called by getIndividualSequence() and other outside programs
		"""
		#'/' must not be put in front of the relative path.
		# otherwise, os.path.join(self.data_dir, dst_relative_path) will only take the path of dst_relative_path.
		
		if subFolder is None:
			subFolder = self.folderName
		filename_part_ls = []
		if self.id:
			filename_part_ls.append(self.id)
		if self.individual_id:
			filename_part_ls.append('indID%s'%(self.individual_id))
		if self.individual_id:
			filename_part_ls.append('code%s'%(self.individual.code))
		if self.sequencer_id:
			filename_part_ls.append("sequencer%s"%(self.sequencer_id))
		if self.sequence_type_id:
			filename_part_ls.append("seqType%s"%(self.sequence_type_id))
		if self.tissue:
			filename_part_ls.append("tissueID%s"%(self.tissue.id))
		if self.filtered is not None:
			filename_part_ls.append("filtered%s"%(self.filtered))
		if self.no_of_chromosomes is not None:
			filename_part_ls.append("%sChromosomes"%(self.no_of_chromosomes))
		if self.sequence_batch_id:
			filename_part_ls.append("batch%s"%(self.sequence_batch_id))
		if self.version is not None:
			filename_part_ls.append("version%s"%(self.version))
		
		filename_part_ls = map(str, filename_part_ls)
		
		dst_relative_path = '%s/%s'%(subFolder, '_'.join(filename_part_ls))
		return dst_relative_path
	
	def constructPSMCPlotLabel(self, **keywords):
		"""
		2013.2.16
		"""
		return '%s_%sISQ%s'%(self.individual.code, self.individual.sex, self.id)
	

class SequenceBatch(Base):
	"""
	2012.7.5 table to store batch names so that IndividualSequence could refer to them.
	"""
	__tablename__ = 'sequence_batch'
	__table_args__ = {'schema':_schemaname_}
	
	id = Column(Integer, primary_key=True)
	short_name = Column(String(256), unique=True)
	description = Column(Text)
	coverage = Column(Integer)
	#individual_ls = ManyToMany('%s.Individual'%(__name__), tablename='individual2batch', \
	#						local_colname='sequence_batch_id')	#2012.7.5
	created_by = Column(String(128))
	updated_by = Column(String(128))
	date_created = Column(DateTime, default=datetime.now())
	date_updated = Column(DateTime)
	#using_options(tablename='sequence_batch', metadata=__metadata__, session=__session__)
	#using_table_options(mysql_engine='InnoDB')
	
	ind_seq_list_in_sequence_batch = relationship('IndividualSequence', back_populates='sequence_batch', cascade='all,delete')
	individual_ls = relationship('Individual',secondary=individual2batch_table,back_populates='sequence_batch_ls')


class Study(Base):
	"""
	2013.3.12 table used to group individuals
		table Individual and Ind2Ind refers to this table.
	"""
	__tablename__ = 'study'
	__table_args__ = {'schema': _schemaname_}
	
	id = Column(Integer, primary_key=True)
	short_name = Column(String(256), unique=True)
	description = Column(Text)
	created_by = Column(String(128))
	updated_by = Column(String(128))
	date_created = Column(DateTime, default=datetime.now())
	date_updated = Column(DateTime)
	#using_options(tablename='study', metadata=__metadata__, session=__session__)
	#using_table_options(mysql_engine='InnoDB')
	
	site_list_in_study = relationship('Site', back_populates='study', cascade='all,delete')
	individual_list_in_study = relationship('Individual', back_populates='study', cascade='all,delete')
	ind2ind_list_in_study = relationship('Ind2Ind', back_populates='study', cascade='all,delete')


class IndividualSequenceFile(Base, AbstractTableWithFilename):
	"""
	#2013.04.30 added file_size
	2012.7.14 add md5sum
	2012.2.27
		add column read_count
	2012.2.10
		add column parent_individual_sequence_file_id
	2012.1.26
	2011-10-31
		a table recording a number of files (split fastq mostly) affiliated with one IndividualSequence entry
	"""
	__tablename__ = 'individual_sequence_file'
	__table_args__ = {'schema':_schemaname_}
	
	id = Column(Integer, primary_key=True)
	#individual_sequence = ManyToOne('%s.IndividualSequence'%(__name__), colname='individual_sequence_id', ondelete='CASCADE', onupdate='CASCADE')
	individual_sequence_id = Column(Integer, ForeignKey(_schemaname_ + '.individual_sequence.id'))
	#individual_sequence_file_raw = ManyToOne('%s.IndividualSequenceFileRaw'%(__name__), colname='individual_sequence_file_raw_id', \
	#										ondelete='CASCADE', onupdate='CASCADE')
	individual_sequence_file_raw_id = Column(Integer, ForeignKey(_schemaname_ + '.individual_sequence_file_raw.id'))
	library = Column(Text)	#id for the preparation library
	split_order = Column(Integer)	# the number that designates the order of this split fastq file within the large file
	mate_id = Column(Integer)	# id of the mate pair. 1 = 1st end. 2 = 2nd end. null = single-end.
	read_count = Column(BigInteger)	#2012.2.27
	base_count = Column(BigInteger)	#2011-8-2
	path = Column(Text, unique=True)	#path to the actual file
	format = Column(String(512))	#fasta, fastq
	quality_score_format = Column(String(512), default='Standard')
		#Standard=Phred+33 (=Sanger), Illumina=Phred+64 (roughly, check pymodule/utils for exact formula)
		# Illumina1.8+ (after 2011-02) is Standard.
	filtered = Column(Integer, default=0)	#0 means not. 1 means yes.
	md5sum = Column(Text, unique=True)
	#parent_individual_sequence_file = ManyToOne('%s.IndividualSequenceFile'%(__name__), colname='parent_individual_sequence_file_id', \
	#							ondelete='CASCADE', onupdate='CASCADE')
	parent_individual_sequence_file_id = Column(Integer, ForeignKey(_schemaname_ + '.individual_sequence_file.id', ondelete='cascade', onupdate='cascade'))
	file_size = Column(BigInteger)	#2013.04.30
	created_by = Column(String(128))
	updated_by = Column(String(128))
	date_created = Column(DateTime, default=datetime.now())
	date_updated = Column(DateTime)
	#using_options(tablename='individual_sequence_file', metadata=__metadata__, session=__session__)
	#using_table_options(mysql_engine='InnoDB')
	#using_table_options(UniqueConstraint('individual_sequence_id', 'library', 'split_order', 'mate_id', 'filtered',\
	#									'parent_individual_sequence_file_id'))
	UniqueConstraint('individual_sequence_id', 'library', 'split_order', 'mate_id', 'filtered', \
					'parent_individual_sequence_file_id', \
					name='individual_sequence_file_ilsmfp')
	
	individual_sequence = relationship('IndividualSequence', back_populates='individual_sequence_file_ls')
	individual_sequence_file_raw = relationship('IndividualSequenceFileRaw', back_populates='individual_sequence_file_ls')
	parent_individual_sequence_file = relationship('IndividualSequenceFile')
	
	
	def constructRelativePath(self, subFolder='individual_sequence', sourceFilename="", **keywords):
		"""
		2012.7.13 
		"""
		folderRelativePath = self.individual_sequence.constructRelativePath(subFolder=subFolder)
		relativePath = os.path.join(folderRelativePath, '%s_%s'%(self.id, sourceFilename))
		return relativePath
	


class IndividualSequenceFileRaw(Base, AbstractTableWithFilename):
	"""
	2013.04.03 added column original_path, to take the meaning of "path"
		"path" is now reserved for db-affiliated file.
	2012.7.12 add column file_size
	2012.4.30
		add column mate_id
	2012.2.27
		add column read_count
	2012.1.26
		this table is used to store the bam files from WUSTL. The bam files hold either single-end or paired-end reads
			in one file.
		Technically, this table is parent of IndividualSequenceFile.
	"""
	__tablename__ = 'individual_sequence_file_raw'
	__table_args__ = {'schema':_schemaname_}
	
	id = Column(Integer, primary_key=True)
	#individual_sequence = ManyToOne('%s.IndividualSequence'%(__name__), colname='individual_sequence_id', ondelete='CASCADE', onupdate='CASCADE')
	individual_sequence_id = Column(Integer, ForeignKey(_schemaname_ + '.individual_sequence.id'))
	#individual_sequence_file_ls = OneToMany("%s.IndividualSequenceFile"%(__name__))
	library = Column(Text)	#id for the preparation library
	mate_id = Column(Integer)	#2012.4.30 to handle fastq raw sequence files.
	read_count = Column(BigInteger)	#2012.2.27
	base_count = Column(BigInteger)
	path = Column(Text)	#path to the file
	original_path = Column(Text)	#path to the original file
	md5sum = Column(Text, unique=True)	#used to identify each raw file
	quality_score_format = Column(String(512), default='Standard')
		#Standard=Phred+33 (=Sanger), Illumina=Phred+64 (roughly, check pymodule/utils for exact formula)
		# Illumina1.8+ (after 2011-02) is Standard.
	file_size = Column(BigInteger)	#2012.7.12
	created_by = Column(String(128))
	updated_by = Column(String(128))
	date_created = Column(DateTime, default=datetime.now())
	date_updated = Column(DateTime)
	#using_options(tablename='individual_sequence_file_raw', metadata=__metadata__, session=__session__)
	#using_table_options(mysql_engine='InnoDB')
	#using_table_options(UniqueConstraint('library', 'md5sum', 'mate_id'))
	UniqueConstraint('library', 'md5sum', 'mate_id', name='individual_sequence_file_raw_library_md5sum_mate_id')
	
	ind_alig_list_in_ind_seq_file_raw = relationship('IndividualAlignment', back_populates='individual_sequence_file_raw', cascade='all,delete')
	individual_sequence = relationship('IndividualSequence', back_populates='individual_sequence_file_raw_ls')
	individual_sequence_file_ls = relationship('IndividualSequenceFile', back_populates='individual_sequence_file_raw', cascade='all,delete')

class Tissue(Base, TableClass):
	"""
	2011-5-9
	"""
	__tablename__ = 'tissue'
	__table_args__ = {'schema':_schemaname_}
	
	id = Column(Integer, primary_key=True)
	short_name = Column(String(512), unique=True)
	description = Column(Text)
	created_by = Column(String(128))
	updated_by = Column(String(128))
	date_created = Column(DateTime, default=datetime.now())
	date_updated = Column(DateTime)
	#using_options(tablename='tissue', metadata=__metadata__, session=__session__)
	#using_table_options(mysql_engine='InnoDB')

	ind_seq_in_tissue = relationship('IndividualSequence', back_populates='tissue', cascade='all,delete')

class IndividualSequence2Sequence(Base, TableClass):
	"""
	2011-3-3
	"""
	__tablename__ = 'individual_seq2seq_map'
	__table_args__ = {'schema':_schemaname_}
	
	id = Column(Integer, primary_key=True)
	#individual1_sequence = ManyToOne('%s.IndividualSequence'%(__name__), 
	#								colname='individual1_sequence_id', ondelete='CASCADE', onupdate='CASCADE')
	individual1_sequence_id = Column(Integer, ForeignKey(_schemaname_ + '.individual_sequence.id'))
	individual1_chr = Column(String(512))
	individual1_start = Column(Integer)
	individual1_stop = Column(Integer)
	#individual2_sequence = ManyToOne('%s.IndividualSequence'%(__name__), 
	#								colname='individual2_sequence_id', ondelete='CASCADE', onupdate='CASCADE')
	individual2_sequence_id = Column(Integer, ForeignKey(_schemaname_ + '.individual_sequence.id'))
	individual2_chr = Column(String(512))
	individual2_start = Column(Integer)
	individual2_stop = Column(Integer)
	created_by = Column(String(128))
	updated_by = Column(String(128))
	date_created = Column(DateTime, default=datetime.now())
	date_updated = Column(DateTime)
	#using_options(tablename='individual_seq2seq_map')
	#using_table_options(mysql_engine='InnoDB')
	#using_table_options(UniqueConstraint('individual1_sequence_id', 'individual1_chr', 'individual1_start', 'individual1_stop',\
	#									'individual2_sequence_id', 'individual2_chr', 'individual2_start', 'individual2_stop'))
	UniqueConstraint('individual1_sequence_id', 'individual1_chr', 'individual1_start', 'individual1_stop', \
					'individual2_sequence_id', 'individual2_chr', 'individual2_start', 'individual2_stop', \
					name='individual_seq2seq_map_iiiiiiii')
	individual1_sequence = relationship('IndividualSequence', back_populates='ind_seq2_list_in_ind1_seq', foreign_keys='IndividualSequence2Sequence.individual1_sequence_id')
	individual2_sequence = relationship('IndividualSequence', back_populates='ind_seq2_list_in_ind2_seq', foreign_keys='IndividualSequence2Sequence.individual2_sequence_id')

class Locus(Base, TableClass):
	"""
	2012.5.2
		make ref_seq_id, ref_ind_seq_id, locus_type_id required (not null)
	2012.4.29
		change locus_method_ls to locus_type
		include locus_type_id into the unique key constraint
	2011-4-5
		add table locus_method
	2011-2-3
		add ref_allele, alt_allele
	"""
	__tablename__ = 'locus'
	__table_args__ = {'schema':_schemaname_}
	
	id = Column(Integer, primary_key=True)
	chromosome = Column(String(512))
	start = Column(Integer)
	stop = Column(Integer)
	#ref_seq = ManyToOne('%s.AlleleSequence'%(__name__), colname='ref_seq_id', required=True, ondelete='CASCADE', onupdate='CASCADE')
	ref_seq_id = Column(Integer, ForeignKey(_schemaname_ + '.allele_sequence.id'))
	#alt_seq = ManyToOne('%s.AlleleSequence'%(__name__), colname='alt_seq_id', ondelete='CASCADE', onupdate='CASCADE')
	alt_seq_id = Column(Integer, ForeignKey(_schemaname_ + '.allele_sequence.id'))
	#ref_sequence = ManyToOne('%s.IndividualSequence'%(__name__), colname='ref_ind_seq_id', required=True, ondelete='CASCADE', onupdate='CASCADE')
	ref_ind_seq_id = Column(Integer, ForeignKey(_schemaname_ + '.individual_sequence.id'))
	#locus_type = ManyToOne('%s.LocusType'%(__name__), colname='locus_type_id', required=True, ondelete='CASCADE', onupdate='CASCADE')
	locus_type_id = Column(Integer, ForeignKey(_schemaname_ + '.locus_type.id'))
	## which study, or SNPs/ indels 
	
	#locus_method_ls = ManyToMany('LocusMethod',tablename='locus2locus_method', local_colname='locus_id', \
	#							remote_colname='locus_method_id')
	created_by = Column(String(128))
	updated_by = Column(String(128))
	date_created = Column(DateTime, default=datetime.now())
	date_updated = Column(DateTime)
	#using_options(tablename='locus', metadata=__metadata__, session=__session__)
	#using_table_options(mysql_engine='InnoDB')
	#using_table_options(UniqueConstraint('chromosome', 'start', 'stop', 'ref_ind_seq_id', 'locus_type_id'))
	UniqueConstraint('chromosome', 'start', 'stop', 'ref_ind_seq_id', 'locus_type_id', \
					name='locus_chromosome_start_stop_ref_ind_seq_id_locus_type_id')
	ref_seq = relationship('AlleleSequence', back_populates='locus_list_in_ref_seq', foreign_keys='Locus.ref_seq_id')
	alt_seq = relationship('AlleleSequence', back_populates='locus_list_in_alt_seq', foreign_keys='Locus.alt_seq_id')
	ref_sequence = relationship('IndividualSequence', back_populates='locus_list_in_ref_sequence')
	locus_type = relationship('LocusType', back_populates='locus_list_in_locus_type')
	locus_annotation_list_in_locus = relationship('LocusAnnotation', back_populates='locus', cascade='all,delete')
	locus_context_list_in_locus = relationship('LocusContext', back_populates='locus', cascade='all,delete')
	genotype_list_in_locus = relationship('Genotype', back_populates='locus', cascade='all,delete', foreign_keys='Genotype.locus_id')
	genotype_list_in_target_locus = relationship('Genotype', back_populates='target_locus', cascade='all,delete',
											foreign_keys='Genotype.target_locus_id')
	locus_score_list_in_locus = relationship('LocusScore',back_populates='locus',cascade='all,delete')
	
class LocusScore(Base, TableClass):
	"""
	2011-4-5
		score of locus
	"""
	__tablename__ = 'locus_score'
	__table_args__ = {'schema':_schemaname_}
	
	id = Column(Integer, primary_key=True)
	#locus = ManyToOne('%s.Locus'%(__name__), colname='locus_id', ondelete='CASCADE', onupdate='CASCADE')
	locus_id = Column(Integer,ForeignKey(_schemaname_+'.locus.id'))
	#score_method = ManyToOne('%s.ScoreMethod'%(__name__), colname='score_method_id', ondelete='CASCADE', onupdate='CASCADE')
	score_method_id = Column(Integer, ForeignKey(_schemaname_ + '.score_method.id'))
	score = Column(Float)
	rank = Column(Integer)
	#object = Column(LargeBinary(134217728), deferred=True)	#a python dictionary to store other attributes
	object = Column(LargeBinary(134217728))
	created_by = Column(String(128))
	updated_by = Column(String(128))
	date_created = Column(DateTime, default=datetime.now())
	date_updated = Column(DateTime)
	#using_options(tablename='locus_score', metadata=__metadata__, session=__session__)
	#using_table_options(mysql_engine='InnoDB')
	#using_table_options(UniqueConstraint('locus_id', 'score_method_id'))
	UniqueConstraint('locus_id', 'score_method_id', name='locus_score_locus_id_score_method_id')

	score_method = relationship('ScoreMethod', back_populates='locus_score_list_in_score_method')
	locus = relationship('Locus',back_populates='locus_score_list_in_locus')

class LocusAnnotation(Base):
	"""
	2012.5.18 add field which_codon
	2012.5.14 add fields for structural variation
	2012.4.25
		copied from SNPAnnotation of Stock_250kDB
	2009-01-05
		table to store annotation of SNPs, like synonymous, non-synonymous, ...
		information finer than SnpsContext
	"""
	__tablename__ = 'locus_annotation'
	__table_args__ = {'schema':_schemaname_}
	
	id = Column(Integer, primary_key=True)
	#locus = ManyToOne('%s.Locus'%__name__, colname='locus_id', ondelete='CASCADE', onupdate='CASCADE')
	locus_id = Column(Integer, ForeignKey(_schemaname_ + '.locus.id'))
	#locus_context = ManyToOne('%s.LocusContext'%__name__, colname='locus_context_id', ondelete='CASCADE', onupdate='CASCADE')
	locus_context_id = Column(Integer, ForeignKey(_schemaname_ + '.locus_context.id'))
	gene_id = Column(Integer)
	gene_commentary_id = Column(Integer)
	gene_segment_id = Column(Integer)
	#locus_annotation_type = ManyToOne('%s.LocusAnnotationType'%__name__, colname='locus_annotation_type_id', ondelete='CASCADE', onupdate='CASCADE')
	locus_annotation_type_id = Column(Integer, ForeignKey(_schemaname_ + '.locus_annotation_type.id'))
	which_exon_or_intron = Column(Integer)
	pos_within_codon = Column(Integer)	#which position in the tri-nucleotide codon, this locus is at. for synonymou/non-syn nucleotide changes.
	which_codon = Column(Integer)	#2012.5.18 which AA this locus affects if synonymous, or non-synonymous, 
	label = Column(Text)	#what type of gene segment it is, exon, cds, intron, 3UTR, 5UTR
	utr_number = Column(Integer)	#2012.5.14 
	cds_number = Column(Integer)
	intron_number = Column(Integer)
	exon_number = Column(Integer)
	overlap_length = Column(Integer)	#for structural variation
	overlap_fraction_in_gene = Column(Float)		#for structural variation
	overlap_fraction_in_locus = Column(Float)	#for structural variation
	comment = Column(String(512))
	created_by = Column(String(200))
	updated_by = Column(String(200))
	date_created = Column(DateTime, default=datetime.now())
	date_updated = Column(DateTime)
	#using_options(tablename='locus_annotation', metadata=__metadata__, session=__session__)
	#using_table_options(mysql_engine='InnoDB')
	#using_table_options(UniqueConstraint('locus_id', 'gene_commentary_id', 'locus_annotation_type_id'))
	UniqueConstraint('locus_id', 'gene_commentary_id', 'locus_annotation_type_id', name='locus_annotation_lgl')
	
	locus = relationship('Locus', back_populates='locus_annotation_list_in_locus')
	locus_context = relationship('LocusContext', back_populates='locus_annotation_list_in_locus_context')
	locus_annotation_type = relationship('LocusAnnotationType', back_populates='locus_annotation_list_in_locus_annotation_type')

class LocusAnnotationType(Base):
	"""
	2012.4.25
		copied from SNPAnnotationType of Stock_250kDB
	2009-01-05
		table to store types of SNP annotation, like synonymous, non-synonymous, ...
	"""
	__tablename__ = 'locus_annotation_type'
	__table_args__ = {'schema':_schemaname_}
	
	id = Column(Integer, primary_key=True)
	short_name = Column(String(256), unique=True)
	description = Column(String(8192))
	created_by = Column(String(200))
	updated_by = Column(String(200))
	date_created = Column(DateTime, default=datetime.now())
	date_updated = Column(DateTime)
	#using_options(tablename='locus_annotation_type', metadata=__metadata__, session=__session__)
	#using_table_options(mysql_engine='InnoDB')
	
	locus_annotation_list_in_locus_annotation_type = relationship('LocusAnnotation', back_populates='locus_annotation_type', cascade='all,delete')

class LocusContext(Base):
	"""
	2012.5.14
		add overlap_length, overlap_fraction_in_locus, overlap_fraction_in_gene
	2012.4.25
		copied from SnpsContext of Stock_250kDB
	"""
	__tablename__ = 'locus_context'
	__table_args__ = {'schema':_schemaname_}
	
	id = Column(Integer, primary_key=True)
	#locus = ManyToOne('%s.Locus'%__name__, colname='locus_id', ondelete='CASCADE', onupdate='CASCADE')
	locus_id = Column(Integer, ForeignKey(_schemaname_ + '.locus.id'))
	disp_pos = Column(Float)	#[0,1) is for within gene fraction, <=-1 is upstream. >=1 is downstream
	gene_id = Column(Integer)
	gene_strand = Column(String(2))
	#left_or_right = Column(String(200), deferred=True)
	left_or_right = Column(String(200))
	#disp_pos_comment = Column(String(2000), deferred=True)
	disp_pos_comment = Column(String(2000))
	overlap_length = Column(Integer)	#for structural variation
	overlap_fraction_in_gene = Column(Float)		#for structural variation
	overlap_fraction_in_locus = Column(Float)	#for structural variation
	#created_by = Column(String(200), deferred=True)
	created_by = Column(String(200))
	#updated_by = Column(String(200), deferred=True)
	updated_by = Column(String(200))
	#date_created = Column(DateTime, default=datetime.now, deferred=True)
	date_created = Column(DateTime, default=datetime.now())
	#date_updated = Column(DateTime, deferred=True)
	date_updated = Column(DateTime)
	#using_options(tablename='locus_context', metadata=__metadata__, session=__session__)
	#using_table_options(mysql_engine='InnoDB')
	#using_table_options(UniqueConstraint('locus_id', 'gene_id'))
	UniqueConstraint('locus_id', 'gene_id', name='locus_context_locus_id_gene_id')
	
	locus_annotation_list_in_locus_context = relationship('LocusAnnotation', back_populates='locus_context', cascade='all,delete')
	locus = relationship('Locus', back_populates='locus_context_list_in_locus')

class ScoreMethod(Base, TableClass):
	"""
	2011-4-5
	"""
	__tablename__ = 'score_method'
	__table_args__ = {'schema':_schemaname_}
	
	id = Column(Integer, primary_key=True)
	short_name = Column(String(256), unique=True)
	filename = Column(String(1000), unique=True)
	original_filename = Column(Text)
	description = Column(Text)
	min_maf = Column(Float)
	#genotype_method = ManyToOne('%s.GenotypeMethod'%__name__, colname='genotype_method_id', ondelete='CASCADE', onupdate='CASCADE')
	genotype_method_id = Column(Integer, ForeignKey(_schemaname_ + '.genotype_method.id'))
	#analysis_method = ManyToOne('%s.AnalysisMethod'%__name__, colname='analysis_method_id', ondelete='CASCADE', onupdate='CASCADE')
	analysis_method_id = Column(Integer, ForeignKey(_schemaname_ + '.analysis_method.id'))
	#phenotype_method = ManyToOne('%s.PhenotypeMethod'%__name__, colname='phenotype_method_id', ondelete='CASCADE', onupdate='CASCADE')
	phenotype_method_id = Column(Integer, ForeignKey(_schemaname_ + '.phenotype_method.id'))
	#transformation_method = ManyToOne('%s.TransformationMethod'%__name__, colname='transformation_method_id', ondelete='CASCADE', onupdate='CASCADE')
	transformation_method_id = Column(Integer, ForeignKey(_schemaname_ + '.transformation_method.id'))
	#score_method_type = ManyToOne('%s.ScoreMethodType'%__name__, colname='score_method_type_id', ondelete='CASCADE', onupdate='CASCADE')
	score_method_type_id = Column(Integer, ForeignKey(_schemaname_ + '.score_method_type.id'))
	created_by = Column(String(128))
	updated_by = Column(String(128))
	date_created = Column(DateTime, default=datetime.now())
	date_updated = Column(DateTime)
	#using_options(tablename='score_method', metadata=__metadata__, session=__session__)
	#using_table_options(mysql_engine='InnoDB')
	#using_table_options(UniqueConstraint('genotype_method_id', 'analysis_method_id', 'phenotype_method_id', \
	#							'score_method_type_id', 'transformation_method_id'))
	UniqueConstraint('genotype_method_id', 'analysis_method_id', 'phenotype_method_id', \
					'score_method_type_id', 'transformation_method_id', \
					name='score_method_gapst')
	locus_score_list_in_score_method = relationship('LocusScore', back_populates='score_method', cascade='all,delete')
	genotype_method = relationship('GenotypeMethod', back_populates='score_method_list_in_genotype_method')
	analysis_method = relationship('AnalysisMethod', back_populates='score_method_list_in_analysis_method')
	phenotype_method = relationship('PhenotypeMethod', back_populates='score_method_list_in_phenotype_method')
	transformation_method = relationship('TransformationMethod', back_populates='score_method_list_in_transformation_method')
	score_method_type = relationship('ScoreMethodType', back_populates='score_method_list_in_score_method_type')

class TransformationMethod(Base):
	"""
	2011-4-5
	"""
	__tablename__ = 'transformation_method'
	__table_args__ = {'schema':_schemaname_}
	
	id = Column(Integer, primary_key=True)
	name = Column(String(30))
	description = Column(Text)
	formular = Column(String(100))
	function = Column(String(20))
	#using_options(tablename='transformation_method', metadata=__metadata__, session=__session__)
	#using_table_options(mysql_engine='InnoDB')

	score_method_list_in_transformation_method = relationship('ScoreMethod', back_populates='transformation_method', cascade='all,delete')
	
class ScoreMethodType(Base):
	"""
	2011-4-5
	"""
	__tablename__ = 'score_method_type'
	__table_args__ = {'schema':_schemaname_}
	
	id = Column(Integer, primary_key=True)
	short_name = Column(String(30), unique=True)
	method_description = Column(String(8000))
	data_description = Column(String(8000))
	created_by = Column(String(200))
	updated_by = Column(String(200))
	date_created = Column(DateTime, default=datetime.now())
	date_updated = Column(DateTime)
	#using_options(tablename='score_method_type', metadata=__metadata__, session=__session__)
	#using_table_options(mysql_engine='InnoDB')

	score_method_list_in_score_method_type = relationship('ScoreMethod', back_populates='score_method_type', cascade='all,delete')

class Sequencer(Base):
	"""
	"""
	__tablename__ = 'sequencer'
	__table_args__ = {'schema':_schemaname_}
	
	id = Column(Integer, primary_key=True)
	short_name = Column(String(128), unique=True)
	description = Column(String(8000))
	#seq_center = ManyToOne('%s.SeqCenter'%(__name__), colname='seq_center_id', ondelete='CASCADE', onupdate='CASCADE')
	seq_center_id = Column(Integer, ForeignKey(_schemaname_ + '.seq_center.id'))
	center_type= Column(String(200)) 
	created_by = Column(String(200))
	updated_by = Column(String(200))
	date_created = Column(DateTime, default=datetime.now())
	date_updated = Column(DateTime)
	#using_options(tablename='sequencer', metadata=__metadata__, session=__session__)
	#using_table_options(mysql_engine='InnoDB')
	
	ind_seq_in_sequencer = relationship('IndividualSequence', back_populates='sequencer', cascade='all,delete')
	seq_center = relationship('SeqCenter', back_populates='sequencer_ls')

class SeqCenter(Base):
	"""
	20170407 sequencing center
	"""
	__tablename__ = 'seq_center'
	__table_args__ = {'schema':_schemaname_}
	
	id = Column(Integer, primary_key=True)
	short_name = Column(String(256), unique=True)
	description = Column(String(8000))
	center_short_name = Column(String(200))
	#sequencer_ls = OneToMany("%s.Sequencer"%(__name__))
	created_by = Column(String(200))
	updated_by = Column(String(200))
	date_created = Column(DateTime, default=datetime.now())
	date_updated = Column(DateTime)
	#using_options(tablename='seq_center', metadata=__metadata__, session=__session__)
	#using_table_options(mysql_engine='InnoDB')
	
	sequencer_ls = relationship('Sequencer', back_populates='seq_center', cascade='all,delete')

class SequenceType(Base):
	"""
	2013.3.13
		PE
		genome
		scaffold
		BAC
		SR
		simulationPE
	"""
	__tablename__ = 'sequence_type'
	__table_args__ = {'schema':_schemaname_}
	
	id = Column(Integer, primary_key=True)
	short_name = Column(String(30), unique=True)
	description = Column(String(8000))
	read_length_mean = Column(BigInteger)
	paired_end = Column(Integer, default=0)
	insert_size_mean = Column(Integer)
	insert_size_variance = Column(Float)
	per_base_error_rate = Column(Float)
	created_by = Column(String(200))
	updated_by = Column(String(200))
	date_created = Column(DateTime, default=datetime.now())
	date_updated = Column(DateTime)
	#using_options(tablename='sequence_type', metadata=__metadata__, session=__session__)
	#using_table_options(mysql_engine='InnoDB')

	ind_seq_in_sequence_type = relationship('IndividualSequence', back_populates='sequence_type', cascade='all,delete')

class LocusType(Base, TableClass):
	"""
	2012.4.29
		renamed to LocusType
	2011-4-5
		to mark different sets of loci
	"""
	__tablename__ = 'locus_type'
	__table_args__ = {'schema':_schemaname_}
	
	id = Column(Integer, primary_key=True)
	short_name = Column(String(256), unique=True)
	description = Column(Text)
	created_by = Column(String(128))
	updated_by = Column(String(128))
	date_created = Column(DateTime, default=datetime.now())
	date_updated = Column(DateTime)
	#using_options(tablename='locus_type', metadata=__metadata__, session=__session__)
	#using_table_options(mysql_engine='InnoDB')

	locus_list_in_locus_type = relationship('Locus', back_populates='locus_type', cascade='all,delete')

class AlleleType(Base, TableClass):
	__tablename__ = 'allele_type'
	__table_args__ = {'schema':_schemaname_}
	
	id = Column(Integer, primary_key=True)
	short_name = Column(String(256), unique=True)
	description = Column(Text)
	created_by = Column(String(128))
	updated_by = Column(String(128))
	date_created = Column(DateTime, default=datetime.now())
	date_updated = Column(DateTime)
	#using_options(tablename='allele_type', metadata=__metadata__, session=__session__)
	#using_table_options(mysql_engine='InnoDB')
	
	genotype_list_in_allele_type = relationship('Genotype', back_populates='allele_type', cascade='all,delete')

class AlleleSequence(Base, TableClass):
	"""
	2011-2-4
		to store the base(s) of an allele
	"""
	__tablename__ = 'allele_sequence'
	__table_args__ = {'schema':_schemaname_}
	
	id = Column(Integer, primary_key=True)
	sequence = Column(Text)
	sequence_length = Column(Integer)
	comment = Column(Text)
	created_by = Column(String(128))
	updated_by = Column(String(128))
	date_created = Column(DateTime, default=datetime.now())
	date_updated = Column(DateTime)
	#using_options(tablename='allele_sequence', metadata=__metadata__, session=__session__)
	#using_table_options(mysql_engine='InnoDB')
	#using_table_options(UniqueConstraint('sequence'))
	UniqueConstraint('sequence', name='allele_sequence_sequence')
	
	locus_list_in_ref_seq = relationship('Locus', back_populates='ref_seq', cascade='all,delete', foreign_keys='Locus.ref_seq_id')
	locus_list_in_alt_seq = relationship('Locus', back_populates='alt_seq', cascade='all,delete', foreign_keys='Locus.alt_seq_id')
	genotype_list_in_allele_sequence = relationship('Genotype', back_populates='allele_sequence', cascade='all,delete')
        
class Genotype(Base, TableClass):
	"""
	2011-2-3
		
	"""
	__tablename__ = 'genotype'
	__table_args__ = {'schema':_schemaname_}
	
	id = Column(Integer, primary_key=True)
	#individual = ManyToOne('%s.Individual'%(__name__), colname='individual_id', ondelete='CASCADE', onupdate='CASCADE')
	individual_id = Column(Integer, ForeignKey(_schemaname_ + '.individual.id'))
	#locus = ManyToOne('%s.Locus'%(__name__), colname='locus_id', ondelete='CASCADE', onupdate='CASCADE')
	locus_id = Column(Integer, ForeignKey(_schemaname_ + '.locus.id'))	
	chromosome_copy = Column(Integer, default=0)
	#on which chromosome copy (multi-ploid), 0=unknown (for un-phased), 1 =1st chromosome, so on
	
	#allele_type = ManyToOne('%s.AlleleType'%(__name__), colname='allele_type_id', ondelete='CASCADE', onupdate='CASCADE')
	allele_type_id = Column(Integer, ForeignKey(_schemaname_ + '.allele_type.id'))
	# SNP/MNP/Indel/inversion
	
	#allele_sequence = ManyToOne('%s.AlleleSequence'%(__name__), colname='allele_sequence_id', ondelete='CASCADE', onupdate='CASCADE')
	allele_sequence_id = Column(Integer, ForeignKey(_schemaname_ + '.allele_sequence.id'))
	#the actual allele
	allele_sequence_length = Column(Integer)
	score = Column(Float)
	#target_locus = ManyToOne('%s.Locus'%(__name__), colname='target_locus_id', ondelete='CASCADE', onupdate='CASCADE')	#for translocated allele only
	target_locus_id = Column(Integer, ForeignKey(_schemaname_ + '.locus.id'))
	#genotype_method = ManyToOne('%s.GenotypeMethod'%(__name__), colname='genotype_method_id', ondelete='CASCADE', onupdate='CASCADE')
	genotype_method_id = Column(Integer, ForeignKey(_schemaname_ + '.genotype_method.id'))
	comment = Column(Text)
	created_by = Column(String(128))
	updated_by = Column(String(128))
	date_created = Column(DateTime, default=datetime.now())
	date_updated = Column(DateTime)
	#using_options(tablename='genotype', metadata=__metadata__, session=__session__)
	#using_table_options(mysql_engine='InnoDB')
	#using_table_options(UniqueConstraint('individual_id', 'locus_id', 'chromosome_copy'))
	UniqueConstraint('individual_id', 'locus_id', 'chromosome_copy', name='genotype_individual_id_locus_id_chromosome_copy')
	
	individual = relationship('Individual', back_populates='genotype_list_in_individual')
	locus = relationship('Locus', back_populates='genotype_list_in_locus', foreign_keys='Genotype.locus_id')
	allele_type = relationship('AlleleType', back_populates='genotype_list_in_allele_type')
	allele_sequence = relationship('AlleleSequence', back_populates='genotype_list_in_allele_sequence')
	target_locus = relationship('Locus', back_populates='genotype_list_in_target_locus', foreign_keys='Genotype.target_locus_id')
	genotype_method = relationship('GenotypeMethod', back_populates='genotype_list_in_genotype_method')



class GenotypeMethod(Base, AbstractTableWithFilename):
	"""
	2013.3.6 add column is_phased
	2012.9.6 add column min_neighbor_distance
	2012.7.12 add meta-columns: min_depth ... min_MAF
	2012.7.11
		get rid of all file-related columns, moved to GenotypeFile.
	2011-2-4
		file format:
						locus1	locus2
			individual1	allele1/allele2
			individual2
	"""
	__tablename__ = 'genotype_method'
	__table_args__ = {'schema':_schemaname_}
	
	id = Column(Integer, primary_key=True)
	short_name = Column(String(256), unique=True)
	description = Column(Text)
	path = Column(Text, unique=True)
	#individual_alignment_ls = ManyToMany("%s.IndividualAlignment"%(__name__), tablename='genotype_method2individual_alignment', \
	#						local_colname='genotype_method_id')
	#ref_sequence = ManyToOne('%s.IndividualSequence'%(__name__), colname='ref_ind_seq_id', ondelete='CASCADE', onupdate='CASCADE')
	ref_ind_seq_id = Column(Integer, ForeignKey(_schemaname_ + '.individual_sequence.id'))
	#parent = ManyToOne("%s.GenotypeMethod"%(__name__), colname='parent_id', ondelete='CASCADE', onupdate='CASCADE')
	parent_id = Column(Integer, ForeignKey(_schemaname_ + '.genotype_method.id', ondelete='CASCADE', onupdate='CASCADE'))
	#genotype_file_ls = OneToMany('%s.GenotypeFile'%(__name__))	#2012.7.17
	no_of_individuals = Column(Integer)
	no_of_loci = Column(BigInteger)
	min_depth = Column(Float)
	max_depth = Column(Float)
	max_missing_rate = Column(Float)
	min_maf = Column(Float)
	min_neighbor_distance = Column(Integer)
	is_phased = Column(Integer, default=0)	#2013.3.6. unphased = 0, phased = 1
	created_by = Column(String(128))
	updated_by = Column(String(128))
	date_created = Column(DateTime, default=datetime.now())
	date_updated = Column(DateTime)
	#using_options(tablename='genotype_method', metadata=__metadata__, session=__session__)
	#using_table_options(mysql_engine='InnoDB')
	#using_table_options(UniqueConstraint('short_name', 'ref_ind_seq_id'))
	UniqueConstraint('short_name', 'ref_ind_seq_id', name='genotype_method_short_name_ref_ind_seq_id')
	
	individual_alignment_list_in_genotype_method = relationship('IndividualAlignment', back_populates='mask_genotype_method', cascade='all,delete')
	score_method_list_in_genotype_method = relationship('ScoreMethod', back_populates='genotype_method', cascade='all,delete')
	genotype_list_in_genotype_method = relationship('Genotype', back_populates='genotype_method', cascade='all,delete')
	ref_sequence = relationship('IndividualSequence', back_populates='genotype_method_list_in_ind_seq')
	parent = relationship('GenotypeMethod')
	genotype_file_ls = relationship('GenotypeFile', back_populates='genotype_method', cascade='all,delete')
	individual_alignment_ls = relationship('IndividualAlignment', secondary=genotype_method2individual_alignment_table,
										back_populates='genotype_method_ls')
	
	folderName = 'genotype_file'
	def constructRelativePath(self, subFolder='genotype_file', **keywords):
		"""
		2012.7.12
		"""
		#'/' must not be put in front of the relative path.
		# otherwise, os.path.join(self.data_dir, dst_relative_path) will only take the path of dst_relative_path.
		dst_relative_path = '%s/method_%s'%(subFolder, self.id)
		return dst_relative_path
	
	def getFileSize(self, data_dir=None):
		"""
		2012.7.12
		"""
		if data_dir is None:
			data_dir = SunsetDB.get_data_dir()
		return utils.getFileOrFolderSize(os.path.join(data_dir, self.path))
	
class GenotypeFile(Base, AbstractTableWithFilename):
	"""
	2012.8.30 add column no_of_chromosomes
	2012.7.11
		modify it so that this holds files, each of which corresponds to GenotypeMethod.
	2011-2-4
		file format:
			locus.id	allele_order	allele_type	seq.id	score	target_locus
	"""
	__tablename__ = 'genotype_file'
	__table_args__ = {'schema':_schemaname_}
	
	id = Column(Integer, primary_key=True)
	path = Column(Text, unique=True)
	original_path = Column(Text)
	md5sum = Column(Text)	# unique=True
	file_size = Column(BigInteger)	#2012.7.12
	chromosome = Column(Text)
	no_of_chromosomes = Column(BigInteger, default=1)	#2012.8.30 BigInteger in case some huge number of contigs
	no_of_individuals = Column(Integer)
	no_of_loci = Column(BigInteger)
	format = Column(Text)
	#genotype_method = ManyToOne('%s.GenotypeMethod'%(__name__), colname='genotype_method_id', ondelete='CASCADE', onupdate='CASCADE')
	genotype_method_id = Column(Integer, ForeignKey(_schemaname_ + '.genotype_method.id'))
	comment = Column(String(4096))
	created_by = Column(String(128))
	updated_by = Column(String(128))
	date_created = Column(DateTime, default=datetime.now())
	date_updated = Column(DateTime)
	#using_options(tablename='genotype_file', metadata=__metadata__, session=__session__)
	#using_table_options(mysql_engine='InnoDB')
	#using_table_options(UniqueConstraint('genotype_method_id', 'chromosome', 'format', 'no_of_chromosomes'))
	UniqueConstraint('genotype_method_id', 'chromosome', 'format', 'no_of_chromosomes', name='genotype_file_gcfn')
	
	genotype_method = relationship('GenotypeMethod', back_populates='genotype_file_ls')
	
	folderName = 'genotype_file'
	def constructRelativePath(self, subFolder='genotype_file', sourceFilename="", **keywords):
		"""
		2012.8.30
			name differently when the no_of_chromosomes is >1
		2012.7.12
			path relative to self.data_dir
		"""
		folderRelativePath = self.genotype_method.constructRelativePath(subFolder=subFolder)
		#'/' must not be put in front of the relative path.
		# otherwise, os.path.join(self.data_dir, dst_relative_path) will only take the path of dst_relative_path.
		folderRelativePath = folderRelativePath.lstrip('/')
		if self.no_of_chromosomes<=1:
			dst_relative_path = os.path.join(folderRelativePath, '%s_chr_%s_%s'%(self.id, self.chromosome, sourceFilename))
		else:	#files with more than 1 chromosome
			# 2012.8.30
			# it'll be something like genotype_file/method_6_$id_$format_#chromosomes_sourceFilename
			# do not put it in genotype_file/method_6/ folder.
			folderRelativePath= folderRelativePath.rstrip('/')
			dst_relative_path = '%s_%s_%schromosomes_%s'%(folderRelativePath, self.id, self.no_of_chromosomes, sourceFilename)
		return dst_relative_path
	
	def getFileSize(self, data_dir=None):
		"""
		2012.7.12
		"""
		if data_dir is None:
			data_dir = SunsetDB.get_data_dir()
		return utils.getFileOrFolderSize(os.path.join(data_dir, self.path))


class AlignmentDepthIntervalMethod(Base, AbstractTableWithFilename):
	"""
	2013.08.16
	"""
	__tablename__ = 'alignment_depth_interval_method'
	__table_args__ = {'schema':_schemaname_}
	
	id = Column(Integer, primary_key=True)
	short_name = Column(String(256), unique=True)
	description = Column(Text)
	path = Column(Text, unique=True)
	#individual_alignment_ls = ManyToMany("%s.IndividualAlignment"%(__name__), tablename='alignment_depth_interval_method2individual_alignment', \
	#						local_colname='alignment_depth_interval_method_id')
	#parent = ManyToOne("%s.AlignmentDepthIntervalMethod"%(__name__), colname='parent_id', ondelete='CASCADE', onupdate='CASCADE')
	parent_id = Column(Integer, ForeignKey(_schemaname_ + '.alignment_depth_interval_method.id', ondelete='CASCADE', onupdate='CASCADE'))
	#alignment_depth_interval_file_ls = OneToMany('%s.AlignmentDepthIntervalFile'%(__name__))
	#ref_sequence = ManyToOne('%s.IndividualSequence'%(__name__), colname='ref_ind_seq_id', ondelete='CASCADE', onupdate='CASCADE')
	ref_ind_seq_id = Column(Integer, ForeignKey(_schemaname_ + '.individual_sequence.id'))
	no_of_alignments = Column(Integer)
	no_of_intervals = Column(BigInteger)
	sum_median_depth = Column(Float)	#across all alignments
	sum_mean_depth = Column(Float)	#across all alignments
	min_segment_length = Column(Integer)	#2013.09.02 this is a parameter of segmentation program
	created_by = Column(String(128))
	updated_by = Column(String(128))
	date_created = Column(DateTime, default=datetime.now())
	date_updated = Column(DateTime)
	#using_options(tablename='alignment_depth_interval_method', metadata=__metadata__, session=__session__)
	#using_table_options(mysql_engine='InnoDB')
	#using_table_options(UniqueConstraint('short_name', 'parent_id'))
	
	parent = relationship('AlignmentDepthIntervalMethod')
	alignment_depth_interval_file_ls = relationship('AlignmentDepthIntervalFile', back_populates='alignment_depth_interval_method', 
												cascade='all,delete')
	ref_sequence = relationship('IndividualSequence', back_populates='alignemnt_depth_interval_method_list_in_ind_seq')
	individual_alignment_ls = relationship('IndividualAlignment', 
										secondary=alignment_depth_interval_method2individual_alignment_table,
										back_populates="alignment_depth_interval_method_ls")
	
	folderName="alignment_depth_interval_file"
	def constructRelativePath(self, subFolder=None, **keywords):
		"""
		2013.08.16
		"""
		if not subFolder:
			subFolder =  self.folderName
		#'/' must not be put in front of the relative path.
		# otherwise, os.path.join(self.data_dir, dst_relative_path) will only take the path of dst_relative_path.
		dst_relative_path = '%s/method_%s'%(subFolder, self.id)
		return dst_relative_path
	
class AlignmentDepthIntervalFile(Base, AbstractTableWithFilename):
	"""
	2013.08.16
	"""
	__tablename__ = 'alignment_depth_interval_file'
	__table_args__ = (
		UniqueConstraint('alignment_depth_interval_method_id', 'chromosome', 'format', 'no_of_chromosomes', \
		                 name='alignment_depth_interval_file_acfn'),
		{'schema':_schemaname_},
	)
	
	id = Column(Integer, primary_key=True)
	path = Column(Text, unique=True)
	original_path = Column(Text)
	md5sum = Column(Text)	# unique=True
	file_size = Column(BigInteger)	#2012.7.12
	chromosome = Column(Text)
	chromosome_size = Column(BigInteger)
	no_of_chromosomes = Column(BigInteger, default=1)	#2012.8.30 BigInteger in case some huge number of contigs
	no_of_intervals = Column(BigInteger)
	min_interval_value = Column(Float)
	max_interval_value = Column(Float)
	mean_interval_value = Column(Float)
	median_interval_value = Column(Float)
	min_interval_length = Column(Integer)	#2013.08.30 length of the shortest interval in file
	max_interval_length = Column(Integer)
	median_interval_length = Column(Integer)
	format = Column(Text)	# BED or other
	#alignment_depth_interval_method = ManyToOne('%s.AlignmentDepthIntervalMethod'%(__name__), \
	#								colname='alignment_depth_interval_method_id', ondelete='CASCADE', onupdate='CASCADE')
	alignment_depth_interval_method_id = Column(Integer, ForeignKey(_schemaname_ + '.alignment_depth_interval_method.id'))
	comment = Column(String(4096))
	created_by = Column(String(128))
	updated_by = Column(String(128))
	date_created = Column(DateTime, default=datetime.now())
	date_updated = Column(DateTime)
	#using_options(tablename='alignment_depth_interval_file', metadata=__metadata__, session=__session__)
	#using_table_options(mysql_engine='InnoDB')
	#using_table_options(UniqueConstraint('alignment_depth_interval_method_id', 'chromosome', 'format', 'no_of_chromosomes'))
	
	alignment_depth_interval_method = relationship('AlignmentDepthIntervalMethod', back_populates='alignment_depth_interval_file_ls')

	folderName=AlignmentDepthIntervalMethod.folderName
	def constructRelativePath(self, subFolder=None, sourceFilename="", **keywords):
		"""
		2013.8.16
		"""
		if not subFolder:
			subFolder =  self.folderName
		folderRelativePath = self.alignment_depth_interval_method.constructRelativePath(subFolder=subFolder)
		#'/' must not be put in front of the relative path.
		# otherwise, os.path.join(self.data_dir, dst_relative_path) will only take the path of dst_relative_path.
		folderRelativePath = folderRelativePath.lstrip('/')
		
		filename_part_ls = []
		if self.id:
			filename_part_ls.append(self.id)
		
		if self.no_of_chromosomes<=1:
			filename_part_ls.append('chr_%s'%(self.chromosome))
		else:
			filename_part_ls.append('%schromosomes'%(self.no_of_chromosomes))
		
		filename_part_ls = map(str, filename_part_ls)
		fileRelativePath = os.path.join(folderRelativePath, '%s_%s'%('_'.join(filename_part_ls), sourceFilename))
		return fileRelativePath
	
	def getFileSize(self, db=None):
		"""
		2013.8.16
		"""
		return utils.getFileOrFolderSize(os.path.join(db.data_dir, self.path))

class Phenotype(Base, TableClass):
	__tablename__ = 'phenotype'
	__table_args__ = {'schema':_schemaname_}
	
	id = Column(Integer, primary_key=True)
	#individual = ManyToOne('%s.Individual'%(__name__), colname='individual_id', ondelete='CASCADE', onupdate='CASCADE')
	individual_id = Column(Integer, ForeignKey(_schemaname_ + '.individual.id'))
	#phenotype_method = ManyToOne('%s.PhenotypeMethod'%(__name__), colname='phenotype_method_id', ondelete='CASCADE', onupdate='CASCADE')
	phenotype_method_id = Column(Integer, ForeignKey(_schemaname_ + '.phenotype_method.id'))
	value = Column(Float)
	comment = Column(Text)
	replicate = Column(Integer, default=1)
	created_by = Column(String(128))
	updated_by = Column(String(128))
	date_created = Column(DateTime, default=datetime.now())
	date_updated = Column(DateTime)
	#using_options(tablename='phenotype', metadata=__metadata__, session=__session__)
	#using_table_options(mysql_engine='InnoDB')
	#using_table_options(UniqueConstraint('individual_id', 'replicate', 'phenotype_method_id'))
	UniqueConstraint('individual_id', 'replicate', 'phenotype_method_id', name='phenotype_individual_id_replicate_phenotype_method_id')
	
	individual = relationship('Individual', back_populates='phenotype_list_in_individual')
	phenotype_method = relationship('PhenotypeMethod', back_populates='phenotype_list_in_phenotype_method')

class BiologyCategory(Base):
	#2011-5-4
	__tablename__ = 'biology_category'
	__table_args__ = {'schema':_schemaname_}
	
	id = Column(Integer, primary_key=True)
	short_name = Column(String(256), unique=True)
	description = Column(Text)
	#phenotype_method_ls = OneToMany("%s.PhenotypeMethod"%(__name__))
	created_by = Column(String(128))
	updated_by = Column(String(128))
	date_created = Column(DateTime, default=datetime.now())
	date_updated = Column(DateTime)
	#using_options(tablename='biology_category', metadata=__metadata__, session=__session__)
	#using_table_options(mysql_engine='InnoDB')

	phenotype_method_ls = relationship('PhenotypeMethod', back_populates='biology_category', cascade='all,delete')

class PhenotypeMethod(Base, TableClass):
	"""
	2011-5-4
		add biology_category, phenotype_scoring
	2011-3-1
		add user_ls/group_ls information
	"""
	__tablename__ = 'phenotype_method'
	__table_args__ = {'schema':_schemaname_}
	
	id = Column(Integer, primary_key=True)
	short_name = Column(String(256), unique=True)
	description = Column(Text)
	phenotype_scoring = Column(Text)
	#biology_category = ManyToOne("%s.BiologyCategory"%(__name__), colname='biology_category_id', \
	#							ondelete='CASCADE', onupdate='CASCADE')
	biology_category_id = Column(Integer, ForeignKey(_schemaname_ + '.biology_category.id'))
	#collector = ManyToOne("%s.User"%(__name__),colname='collector_id')
	collector_id = Column(Integer, ForeignKey(_schemaname_ + '.acl_user.id'))
	access = Column(Enum("public", "restricted", name="access_enum_type"), default='restricted')
	#group_ls = ManyToMany('%s.Group'%(__name__),tablename='group2phenotype_method', local_colname='phenotype_method_id', remote_colname='group_id')
	#user_ls = ManyToMany('%s.User'%(__name__),tablename='user2phenotype_method', local_colname='phenotype_method_id', remote_colname='user_id')
	created_by = Column(String(128))
	updated_by = Column(String(128))
	date_created = Column(DateTime, default=datetime.now())
	date_updated = Column(DateTime)
	#using_options(tablename='phenotype_method', metadata=__metadata__, session=__session__)
	#using_table_options(mysql_engine='InnoDB')
	
	score_method_list_in_phenotype_method = relationship('ScoreMethod', back_populates='phenotype_method', cascade='all,delete')
	phenotype_list_in_phenotype_method = relationship('Phenotype', back_populates='phenotype_method', cascade='all,delete')
	biology_category = relationship('BiologyCategory', back_populates='phenotype_method_ls')
	collector = relationship('User', back_populates='phenotype_method_list_in_user')
	group_ls = relationship('Group',secondary=group2phenotype_method_table,back_populates='phenotype_method_ls')
	user_ls = relationship('User',secondary=user2phenotype_method_table,back_populates='phenotype_method_ls')
	
	@classmethod
	def getPhenotypesForACL(cls, biology_category_id = None,user = None):
		#
		TableClass = PhenotypeMethod
		query = TableClass.query
		"""
		if biology_category_id is not None:
			query = query.filter(TableClass.biology_category_id==biology_category_id)
		"""
		clause = and_(TableClass.access == 'public')
		if user is not None:
			if user.isAdmin == 'Y':
				return query
			clause = or_(clause, TableClass.collector == user, TableClass.user_ls.any(User.id == user.id),
						TableClass.group_ls.any(Group.id.in_([group.id for group in user.group_ls])))
		query = query.filter(clause)
		return query
	
	def checkACL(self,user):
		if self.access == 'public':
			return True
		if user is None:
			return False
		if user.isAdmin == 'Y':
			return True
		if self.collector == user: 
			return True
		if user in self.user_ls:
			return True
		if [group in self.group_ls for group in user.group_ls]: 
			return True
		return False

class SunsetDB(Database):
	__doc__ = __doc__
	option_default_dict = copy.deepcopy(Database.option_default_dict)
	option_default_dict.update({
						('drivername', 1,):['postgresql', 'v', 1, 'which type of database? mysql or postgresql', ],\
						('hostname', 1, ):['localhost', 'z', 1, 'hostname of the db server', ],\
						('dbname', 1, ):['pmdb', 'd', 1, '',],\
						('schema', 0, ): [None, 'k', 1, 'database schema name', ],\
						})
	def __init__(self, **keywords):
		"""
		2008-10-06
			add option 'pool_recycle' to recycle connection. MySQL typically close connections after 8 hours.
			__metadata__.bind = create_engine(self._url, pool_recycle=self.pool_recycle)
		2008-07-09
		"""
		Database.__init__(self, **keywords)
		self.dbID2monkeyDBEntry = {}	#2012.9.6
		self.READMEClass = README	#2012.12.18 required to figure out data_dir
	
	def isThisAlignmentComplete(self, individual_alignment=None, data_dir=None, returnFalseIfInexitentFile=False):
		"""
		2013.05.04 added argument returnFalseIfInexitentFile
		2013.04.29 whether os.path.isfile(alignmentAbsPath) is true or not is not mandatory anymore.
			Will give a warning though. 
		2013.04.10 bugfix. individual_alignment.path could be None.
		2013.03.28
		"""
		if data_dir is None:
			data_dir = self.data_dir
		if not individual_alignment.path:
			return False
		alignmentAbsPath= os.path.join(data_dir, individual_alignment.path)
		#2012.3.29	check if the alignment exists or not. if it already exists, no alignment jobs.
		if individual_alignment.file_size is not None and individual_alignment.file_size>0:
			if not os.path.isfile(alignmentAbsPath):
				if returnFalseIfInexitentFile:	#2013.05.04
					sys.stderr.write("Warning: Skip alignment file for id=%s, %s does not exist but file_size is recorded in DB.\n"%(individual_alignment.id, alignmentAbsPath))
					return False
				else:
					sys.stderr.write("Warning: Alignment file for id=%s, %s does not exist but file_size is recorded in DB.\n"%(individual_alignment.id, alignmentAbsPath))
			return True
		else:
			return False

	def getUniqueAlleleSequence(self, sequence):
		"""
		2011-2-11
		"""
		db_entry = self.queryTable(AlleleSequence).filter_by(sequence=sequence).first()
		if not db_entry:
			db_entry = AlleleSequence(sequence=sequence, sequence_length=len(sequence))
			self.session.add(db_entry)
			self.session.flush()
		return db_entry
	
	def getAlignmentMethod(self, alignment_method_name=''):
		"""
		2011-5-6
		"""
		db_entry = self.queryTable(AlignmentMethod).filter_by(short_name=alignment_method_name).first()
		if not db_entry:
			db_entry = AlignmentMethod(short_name=alignment_method_name)
			self.session.add(db_entry)
			self.session.flush()
		return db_entry
	
	def getAlignment(self, individual_code=None, individual_id=None, individual=None, individual_sequence_id=None, \
					path_to_original_alignment=None, sequencer_name='GA', sequencer_id=None, \
					sequence_type_name=None, sequence_type_id=None, sequence_format='fastq', \
					ref_individual_sequence_id=10, \
					alignment_method_name='bwa-short-read', alignment_method_id=None, alignment_method=None,\
					alignment_format='bam', subFolder='individual_alignment', \
					createSymbolicLink=False, individual_sequence_filtered=0, read_group_added=None, data_dir=None, \
					outdated_index=0, mask_genotype_method_id=None, parent_individual_alignment_id=None,\
					individual_sequence_file_raw_id=None, md5sum=None, local_realigned=0, read_group=None,\
					reduce_reads=0):
		"""
		i.e.
			individual_alignment = self.db_vervet.getAlignment(individual_sequence_id=self.individual_sequence_id,\
										path_to_original_alignment=None, sequencer=individual_sequence.sequencer,\
										sequence_type=individual_sequence.sequence_type, sequence_format=individual_sequence.format, \
										ref_individual_sequence_id=self.ref_sequence_id, \
										alignment_method_id=self.alignment_method_id, alignment_format=self.format,\
										individual_sequence_filtered=individual_sequence.filtered, read_group_added=1,
										data_dir=data_dir, \
										mask_genotype_method_id=self.mask_genotype_method_id, \
										parent_individual_alignment_id=self.parent_individual_alignment_id,\
										individual_sequence_file_raw_id=self.individual_sequence_file_raw_id,\
										local_realigned=self.local_realigned, read_group=self.read_group)
										
			oneLibraryAlignmentEntry = db_vervet.getAlignment(individual_code=individual_sequence.individual.code, \
												individual_sequence_id=individual_sequence.id,\
									path_to_original_alignment=None, sequencer_id=individual_sequence.sequencer_id, \
									sequence_type_id=individual_sequence.sequence_type_id, sequence_format=individual_sequence.format, \
									ref_individual_sequence_id=refSequence.id, \
									alignment_method_name=alignment_method.short_name, alignment_format=alignment_format,\
									individual_sequence_filtered=individual_sequence.filtered, read_group_added=1,
									data_dir=data_dir, individual_sequence_file_raw_id=minIsqFileRawID,\
									local_realigned=self.local_realigned)
		
		2013.04.11 added argument reduce_reads
		2013.04.09 added argument read_group and process it
		2013.04.05 added argument local_realigned
		2012.9.19 add argument individual_sequence_file_raw_id, individual_id, individual,
			alignment_method_id, alignment_method, md5sum
		2012.7.26 add argument mask_genotype_method_id & parent_individual_alignment_id
		2012.6.13 add argument outdated_index
		2012.2.24
			add argument data_dir
		2011-9-15
			add argument read_group_added
		2011-8-30
			add argument individual_sequence_id
			use constructRelativePathForIndividualAlignment() to come up path
		2011-7-11
			add argument createSymbolicLink. default to False
				if True, create a symbolic from source file to target, instead of cp
		2011-5-6
			subFolder is the name of the folder in self.data_dir that is used to hold the alignment files.
		"""
		if data_dir is None:
			data_dir = self.data_dir
		if individual_code:
			individual = self.getIndividual(individual_code)
		elif individual_id:
			individual = self.queryTable(Individual).get(individual_id)
		
		if individual_sequence_id:
			individual_sequence = self.queryTable(IndividualSequence).get(individual_sequence_id)
			if individual is None:
				individual = individual_sequence.individual
		elif individual:
			individual_sequence = self.getIndividualSequence(individual_id=individual.id, sequencer_name=sequencer_name, \
								sequencer_id=sequencer_id, sequence_type_id=sequence_type_id,\
								sequence_type_name=sequence_type_name,\
								sequence_format=sequence_format, filtered=individual_sequence_filtered)
		else:
			sys.stderr.write("Error: not able to get individual_sequence cuz individual_sequence_id=%s; individual_code=%s; individual_id=%s .\n"%\
							(individual_sequence_id, individual_code, individual_id))
			return None
		
		if alignment_method_name:
			alignment_method = self.getAlignmentMethod(alignment_method_name=alignment_method_name)
		elif alignment_method_id:
			alignment_method = self.queryTable(AlignmentMethod).get(alignment_method_id)
		elif alignment_method is None:
			sys.stderr.write("Error: not able to get alignment_method cuz alignment_method_name=%s and alignment_method_id=%s.\n"%\
							(alignment_method_name, alignment_method_id))
			return None
		
		
		db_entry = self.checkIndividualAlignment(individual_code=None, individual_id=individual.id, individual_sequence_id=individual_sequence_id, \
					path_to_original_alignment=path_to_original_alignment, sequencer_name=sequencer_name, sequencer_id=sequencer_id, \
					sequence_type_name=sequence_type_name, sequence_type_id=sequence_type_id, sequence_format=sequence_format, \
					ref_individual_sequence_id=ref_individual_sequence_id, \
					alignment_method_name=alignment_method_name, alignment_method_id=alignment_method_id, alignment_method=alignment_method,\
					alignment_format=alignment_format, subFolder=subFolder, \
					createSymbolicLink=createSymbolicLink, individual_sequence_filtered=individual_sequence_filtered, \
					read_group_added=read_group_added, data_dir=data_dir, \
					outdated_index=outdated_index, mask_genotype_method_id=mask_genotype_method_id, \
					parent_individual_alignment_id=parent_individual_alignment_id,\
					individual_sequence_file_raw_id=individual_sequence_file_raw_id, md5sum=md5sum, \
					local_realigned=local_realigned, reduce_reads=reduce_reads)
		if not db_entry:
			db_entry = IndividualAlignment(ind_seq_id=individual_sequence.id, ref_ind_seq_id=ref_individual_sequence_id,\
							alignment_method_id=alignment_method.id, format=alignment_format, read_group_added=read_group_added,\
							outdated_index=outdated_index, mask_genotype_method_id=mask_genotype_method_id,\
							parent_individual_alignment_id=parent_individual_alignment_id,\
							individual_sequence_file_raw_id=individual_sequence_file_raw_id, md5sum=md5sum,\
							local_realigned=local_realigned, read_group=read_group, reduce_reads=reduce_reads)
			self.session.add(db_entry)
			self.session.flush()
			
			#copy the file over
			
			if path_to_original_alignment and (os.path.isfile(path_to_original_alignment) or os.path.isdir(path_to_original_alignment)):
				#'/' must not be put in front of the relative path.
				# otherwise, os.path.join(data_dir, dst_relative_path) will only take the path of dst_relative_path.
				dst_relative_path = db_entry.constructRelativePath(subFolder=subFolder)
				#update its path in db to the relative path
				db_entry.path = dst_relative_path
				
				dst_pathname = os.path.join(data_dir, dst_relative_path)
				dst_dir = os.path.join(data_dir, subFolder)
				if not os.path.isdir(dst_dir):	#the upper directory has to be created at this moment.
					commandline = 'mkdir %s'%(dst_dir)
					return_data = runLocalCommand(commandline, report_stderr=True, report_stdout=True)
				if createSymbolicLink:
					commandline = 'ln -s %s %s'%(path_to_original_alignment, dst_pathname)
				else:
					commandline = 'cp -r %s %s'%(path_to_original_alignment, dst_pathname)
				return_data = runLocalCommand(commandline, report_stderr=True, report_stdout=True)
				
				#2011-7-11 copy the samtools index file as well
				indexFname = '%s.bai'%(path_to_original_alignment)
				if os.path.isfile(indexFname):
					dstIndexPathname = '%s.bai'%(dst_pathname)
					if createSymbolicLink:
						commandline = 'ln -s %s %s'%(indexFname, dstIndexPathname)
					else:
						commandline = 'cp -r %s %s'%(indexFname, dstIndexPathname)
					return_data = runLocalCommand(commandline, report_stderr=True, report_stdout=True)
				
				#since the .path has been updated, so add & flush
				self.session.add(db_entry)
				self.session.flush()
		## 2013.04.09 adding the read_group
		if db_entry.read_group is None:	#2013.04.09
			if not read_group:
				read_group = db_entry.getReadGroup()
			db_entry.read_group = read_group
			self.session.add(db_entry)
			self.session.flush()
		return db_entry
	
	
	def getAlignments(self, ref_ind_seq_id=None, ind_seq_id_ls=None, ind_aln_id_ls=None, alignment_method_id=None, \
					data_dir=None, sequence_type_name=None, sequence_type_id=None, outdated_index=0, mask_genotype_method_id=None, \
					parent_individual_alignment_id=None, individual_sequence_file_raw_id_type=1,\
					country_id_ls=None, tax_id_ls=None, excludeAlignmentWithoutLocalFile=True, local_realigned=1,\
					reduce_reads=None, completedAlignment=None, completeAlignmentCheckFunction=None):
		"""
		2013.05.04 added argument completeAlignmentCheckFunction
			if mask_genotype_method_id is 0 or '0', then this requires alignment.mask_genotype_method_id to be null. 
			if mask_genotype_method_id is None or '', then it's not checked .
			ditto for parent_individual_alignment_id
		2013.05.03 added argument completedAlignment
		2013.04.11 added argument reduce_reads
		2013.04.05 added argument local_realigned
		2012.11.29 added argument excludeAlignmentWithoutLocalFile, (exclude an alignment if it does not exist in local storage)
		2012.9.23 add argument country_id_ls and tax_id_ls
		2012.9.20 rename aln_method_id to alignment_method_id
		2012.9.19 & 2012.9.22 add argument individual_sequence_file_raw_id_type:
			1: only all-library-fused libraries, (individual_sequence_file_raw_id is null)
			2: only library-specific alignments, (individual_sequence_file_raw_id is non-null)
			3 or else: both all-library-fused and library-specific alignments, 
				(no filter based on individual_sequence_file_raw_id)
		2012.7.26 added argument mask_genotype_method_id & parent_individual_alignment_id
		2012.6.13 add argument outdated_index
		2012.4.13
			moved from AlignmentToCallPipeline.py
		2012.4.5
			select alignment using AND between all input arguments
		2011-11-27
			add argument sequence_type
		2011-9-16
			order each alignment by id. It is important because this is the order that gatk&samtools take input bams.
			#Read group in each bam is beginned by alignment.id. GATK would arrange bams in the order of read groups.
			# while samtools doesn't do that and vcf-isect could combine two vcfs with columns in different order.
		2011-9-13
			add argument data_dir, to filter out alignments that don't exist on file storage
		2011-8-31
			add argument aln_method_id
		2011-7-12
		
		"""
		if completedAlignment is not None:
			if completedAlignment==0 or completedAlignment=='0' or completedAlignment=='':
				completedAlignment = False
			elif completedAlignment!=False:
				completedAlignment = True
		
		sys.stderr.write("Get alignments to from local_realigned=%s, reduce_reads=%s \n\
	mask_genotype_method_id=%s, parent_individual_alignment_id=%s, completedAlignment=%s ..."%\
							( \
							local_realigned, reduce_reads, \
							mask_genotype_method_id, parent_individual_alignment_id, completedAlignment,\
							)
						)
		if data_dir is None:
			data_dir = self.data_dir
		alignmentLs = []
		TableClass = IndividualAlignment
		query = self.queryTable(TableClass)
		if ind_aln_id_ls:
			sys.stderr.write("Adding filter via %s alignment IDs ... \n"%(len(ind_aln_id_ls)) )
			query = query.filter(TableClass.id.in_(ind_aln_id_ls))
		if ind_seq_id_ls:
			sys.stderr.write("Adding filter via %s sequence IDs  ... \n"%(len(ind_seq_id_ls)))
			query = query.filter(TableClass.ind_seq_id.in_(ind_seq_id_ls))
		if ref_ind_seq_id:
			sys.stderr.write("Adding filter reference sequence ID %s ... \n"%(ref_ind_seq_id))
			query = query.filter_by(ref_ind_seq_id=ref_ind_seq_id)
		if alignment_method_id:
			sys.stderr.write("Adding filter alignment_method_id=%s ... \n"%(alignment_method_id))
			query = query.filter_by(alignment_method_id=alignment_method_id)
		
		if country_id_ls:	#2012.9.23 not sure whether it'll work
			sys.stderr.write("Adding filter %s countries: %s ... \n"%(getattr(country_id_ls, '__len__', returnZeroFunc)(), \
																	repr(country_id_ls)))
			#2013.04.09 five-hierarchy query, old way does not work
			individual_id_ls_query = self.queryTable(Individual).filter(Individual.site.has(Site.country_id.in_(country_id_ls)))
			individual_id_ls = [row.id for row in individual_id_ls_query]
			query = query.filter(TableClass.individual_sequence.has(IndividualSequence.individual_id.in_(individual_id_ls)))
		if tax_id_ls:
			sys.stderr.write("Adding filter %s taxonomies: %s ... \n"%(getattr(tax_id_ls, '__len__', returnZeroFunc)(),\
													repr(tax_id_ls)))
			#2013.04.09 4-hierarchy query, old way does not work
			ind_seq_id_ls_query = self.queryTable(IndividualSequence).filter(IndividualSequence.individual.has(Individual.tax_id.in_(tax_id_ls)))
			ind_seq_id_ls = [row.id for row in ind_seq_id_ls_query]
			query = query.filter(TableClass.ind_seq_id.in_(ind_seq_id_ls))
		
		if not ind_aln_id_ls and not ind_seq_id_ls and not ref_ind_seq_id:
			sys.stderr.write("Both ind_seq_id_ls and ind_aln_id_ls are empty and ref_ind_seq_id is None. no alignment to be fetched.\n")
			sys.exit(3)
		
		if individual_sequence_file_raw_id_type==1:	#only all-library-fused alignments
			sys.stderr.write("Adding filter individual_sequence_file_raw_id_type=%s ... \n"%(individual_sequence_file_raw_id_type))
			query = query.filter(TableClass.individual_sequence_file_raw_id==None)
		elif individual_sequence_file_raw_id_type==2:	#only library-specific alignments
			sys.stderr.write("Adding filter individual_sequence_file_raw_id_type=%s ... \n"%(individual_sequence_file_raw_id_type))
			query = query.filter(TableClass.individual_sequence_file_raw_id!=None)
		else:
			#2012.9.22 do nothing = include both 
			pass
		
		if local_realigned is not None:	#2013.04.09
			sys.stderr.write("Adding filter local_realigned=%s ... \n"%(local_realigned))
			query = query.filter_by(local_realigned=local_realigned)
		if reduce_reads is not None:	#2013.04.11
			sys.stderr.write("Adding filter reduce_reads=%s ... \n"%(reduce_reads))
			query = query.filter_by(reduce_reads=reduce_reads)
		#order by TableClass.id is important because this is the order that gatk&samtools take input bams.
		#Read group in each bam is beginned by alignment.id. GATK would arrange bams in the order of read groups.
		# while samtools doesn't do that and vcf-isect could combine two vcfs with columns in different order.
		if outdated_index is not None:	#2013.04.11
			sys.stderr.write("Adding filter outdated_index=%s ... \n"%(outdated_index))
			query = query.filter_by(outdated_index=outdated_index)
		if mask_genotype_method_id is not None and mask_genotype_method_id!='':	#2013.05.03
			sys.stderr.write("Adding filter mask_genotype_method_id=%s ... \n"%(mask_genotype_method_id))
			if mask_genotype_method_id==0 or mask_genotype_method_id=='0':	#only null
				query = query.filter_by(mask_genotype_method_id=None)
			else:
				query = query.filter_by(mask_genotype_method_id=mask_genotype_method_id)
		if parent_individual_alignment_id is not None and parent_individual_alignment_id!='':
			sys.stderr.write("Adding filter parent_individual_alignment_id=%s ... \n"%(parent_individual_alignment_id))
			if parent_individual_alignment_id==0 or parent_individual_alignment_id=='0':
				query = query.filter_by(parent_individual_alignment_id=None)
			else:
				query = query.filter_by(parent_individual_alignment_id=parent_individual_alignment_id)
		
		query = query.order_by(TableClass.id)
		if sequence_type_name:
			sequence_type = self.getSequenceType(short_name=sequence_type_name)
			sequence_type_id = sequence_type.id
		if sequence_type_id is not None:
			sys.stderr.write("Adding filter sequence_type_id=%s ... \n"%(sequence_type_id))
		
		for row in query:
			if sequence_type_id is not None and row.individual_sequence.sequence_type_id!=sequence_type_id:
				continue
			if row.path:	#it's not None
				if completedAlignment is not None:
					if completeAlignmentCheckFunction is None:
						isAlignmentCompleted = self.isThisAlignmentComplete(individual_alignment=row, data_dir=data_dir)
					else:
						isAlignmentCompleted = completeAlignmentCheckFunction(individual_alignment=row, data_dir=data_dir)
					if completedAlignment==isAlignmentCompleted:
						pass
					else:
						sys.stderr.write("SunsetDB.getAlignments() warning: completeness (%s) of alignment %s, (file=%s, read_group=%s, file_size=%s) does not match given(%s). Skip.\n"%\
										(isAlignmentCompleted, row.id, row.path, row.getReadGroup(), row.file_size, completedAlignment))
						continue
				abs_path = os.path.join(data_dir, row.path)
				if excludeAlignmentWithoutLocalFile:
					if os.path.isfile(abs_path):
						alignmentLs.append(row)
					else:
						sys.stderr.write("SunsetDB.getAlignments() Warning: file %s does not exist. so skipping this alignment.\n"%(abs_path))
				else:
					alignmentLs.append(row)
				
		sys.stderr.write("%s alignments Done.\n"%(len(alignmentLs)))
		return alignmentLs
	
	def getProperAlignmentGivenIndividualID(self, individual_id=None, ref_ind_seq_id=524, alignment_method_id=2):
		"""
		2012.11.26
			Definition of proper alignment:
				1. not outdated
				2. from filtered reads.
				3. ref_ind_seq_id is the most recent reference (524 now)
				4. alignment_method_id is the consensus one(=2).
			
			query the view_alignment_with_country.
		"""
		query_string = "select * from view_alignment_with_country"
		where_condition_ls = ["filtered=1 and outdated_index=0 and ref_ind_seq_id=%s and alignment_method_id=%s "%\
							(ref_ind_seq_id, alignment_method_id)]
		if individual_id:
			where_condition_ls.append("individual_id=%s"%(individual_id))
		query_string = "%s where %s "%(query_string, " and ".join(where_condition_ls))
		query = self.metadata.bind.execute(query_string)
		return query.fetchone()
	
	def getAlignmentsFromVCFSampleIDList(self, sampleIDList=None):
		"""
		2012.8.15
			each sample ID is constructed through IndividualAlignment.getReadGroup()
			use vcfFile.getSampleIDList() to generate sampleIDList. vcfFile is pymodule.VCFFile.
				vcfFile.sample_id_ls is not good because its index 0 is "ref".
		"""
		no_of_samples = 0
		if sampleIDList:
			no_of_samples = len(sampleIDList)
		sys.stderr.write("Getting alignments from %s samples ...\n"%(no_of_samples))
		alignmentLs = []
		if sampleIDList:
			for read_group in sampleIDList:
				individualAlignment = self.parseAlignmentReadGroup(read_group).individualAlignment
				individualAlignment.sampleID = read_group
				alignmentLs.append(individualAlignment)
		sys.stderr.write("%s alignments.\n"%(len(alignmentLs)))
		return alignmentLs
	
	def getAlignmentsFromAlignmentIDList(self, alignmentIDList=None):
		"""
		2013.08.23
		"""
		no_of_samples = 0
		if alignmentIDList:
			no_of_samples = len(alignmentIDList)
		sys.stderr.write("Getting alignments from %s samples ...\n"%(no_of_samples))
		alignmentLs = []
		if alignmentIDList:
			for alignmentID in  alignmentIDList:
				individualAlignment = self.checkIfEntryInTable(TableClass=IndividualAlignment, id=int(alignmentID))
				alignmentLs.append(individualAlignment)
		sys.stderr.write("%s alignments.\n"%(len(alignmentLs)))
		return alignmentLs
	
	def getAlignmentsFromVCFFile(self, inputFname=None):
		"""
		2013.1.2. moved from db/OutputVRCPedigreeInTFAMGivenOrderFromFile.py
			inputFname is a VCFFile containing genotypes from alignments of SunsetDB
		"""
		vcfFile = VCFFile(inputFname=inputFname)
		alignmentLs = self.getAlignmentsFromVCFSampleIDList(vcfFile.getSampleIDList())
		#vcfFile.sample_id_ls is not good because its index 0 is "ref"
		vcfFile.close()
		return alignmentLs
	
	def getMonkeyID2ProperAlignment(self, ref_ind_seq_id=524, alignment_method_id=2, idType=1):
		"""
		2012.11.30
			Definition of proper alignment:
				1. not outdated
				2. from filtered reads.
				3. ref_ind_seq_id is the most recent reference (524 now)
				4. alignment_method_id is the consensus one(=2).
			
			query the view_alignment_with_country.
			
			idType
				1: individual.id
				2: individual.code
		"""
		query_string = "select * from view_alignment_with_country"
		where_condition_ls = ["filtered=1 and outdated_index=0 and ref_ind_seq_id=%s and alignment_method_id=%s "%\
							(ref_ind_seq_id, alignment_method_id)]
		query_string = "%s where %s "%(query_string, " and ".join(where_condition_ls))
		query = self.metadata.bind.execute(query_string)
		monkeyID2ProperAlignment = {}
		for row in query:
			if idType ==1:
				monkeyID = row.individual_id
			else:
				monkeyID = row.code
			if monkeyID not in monkeyID2ProperAlignment:
				monkeyID2ProperAlignment[monkeyID] = row
			else:
				sys.stderr.write("Warning: monkey %s has >1 proper alignment. Only used the 1st one. (ref_ind_seq_id=%s, alignment_method_id=%s)\n"%\
								(ref_ind_seq_id, alignment_method_id))
		return monkeyID2ProperAlignment
	
	def filterAlignments(self, data_dir=None, alignmentLs=None, min_coverage=None, max_coverage=None, \
						individual_site_id=None, sequence_filtered=None,\
						individual_site_id_set=None, mask_genotype_method_id=None, parent_individual_alignment_id=None,\
						country_id_set=None, tax_id_set=None, excludeContaminant=False, is_contaminated=None, \
						excludeTissueIDSet=set([6]),\
						local_realigned=1, reduce_reads=None, completedAlignment=None, \
						alignment_method_id=None, outdated_index=None,\
						completeAlignmentCheckFunction=None, report=True):
		"""
		2013.10.04 added alignment_method_id, outdated_index
		2013.07.03 added argument is_contaminated (whether to fetch contaminated samples or not)
		2013.05.04 added argument completeAlignmentCheckFunction
			if mask_genotype_method_id is 0 or '0', then this requires alignment.mask_genotype_method_id to be null. 
			if mask_genotype_method_id is None or '', then it's not checked .
		2013.05.03 added argument completedAlignment
		2013.04.11 added argument reduce_reads
		2013.04.05 added argument local_realigned
		2013.3.15 use individual_sequence.is_contaminated, instead of individual_sequence.individual.is_contaminated
		2013.3.15 added min_coverage
		2012.10.2 add argument excludeTissueIDSet, default=6 (RNASASamples). Most alignments have either tissue_id=5 (ACD-blood) or null.
		2012.9.27 add argument excludeContaminant
		2012.9.22 add argument country_id_set, tax_id_set
		2012.7.26 added argument mask_genotype_method_id & parent_individual_alignment_id
		2012.5.8
			bugfix, individual_site_id_set could be none. so no_of_sites is un-defined.
		2012.4.13
			moved from AlignmentToCallPipeline.py
			become classmethod
			add argument individual_site_id_set
		2012.4.2
			add argument sequence_filtered
		2011-11-22
			447 in "individual_site_id=447" is VRC.
		"""
		if completedAlignment is not None:
			if completedAlignment==0 or completedAlignment=='0' or completedAlignment=='':
				completedAlignment = False
			elif completedAlignment!=False:
				completedAlignment = True
		if report:
			sys.stderr.write("Filter %s alignments to select individual_sequence, %s<=coverage <=%s & site-id=%s & \n\
	sequence_filtered=%s & from %s sites & %s countries & %s taxonomies & local_realigned=%s, reduce_reads=%s \n\
	excludeContaminant=%s, is_contaminated=%s, excludeTissueIDSet=%s, \n\
	mask_genotype_method_id=%s, parent_individual_alignment_id=%s, completedAlignment=%s, \n\
	alignment_method_id=%s, outdated_index=%s  ..."%\
							(len(alignmentLs), min_coverage, max_coverage, individual_site_id, sequence_filtered, \
							getattr(individual_site_id_set, '__len__', returnZeroFunc)(),\
							getattr(country_id_set, '__len__', returnZeroFunc)(),\
							getattr(tax_id_set, '__len__', returnZeroFunc)(),\
							local_realigned, reduce_reads, excludeContaminant, is_contaminated, repr(excludeTissueIDSet),\
							mask_genotype_method_id, parent_individual_alignment_id, completedAlignment,\
							alignment_method_id, outdated_index,\
							)
						)
		newAlignmentLs = []
		for alignment in alignmentLs:
			if not alignment:
				continue
			if min_coverage is not None and alignment.individual_sequence.coverage<min_coverage:
				continue
			if max_coverage is not None and alignment.individual_sequence.coverage>max_coverage:
				continue
			if individual_site_id is not None and alignment.individual_sequence.individual.site_id!=individual_site_id:
				continue
			if sequence_filtered is not None and alignment.individual_sequence.filtered!=sequence_filtered:
				continue
			if individual_site_id_set and alignment.individual_sequence.individual.site_id not in individual_site_id_set:
				#2012.4.13
				continue
			if mask_genotype_method_id is not None and mask_genotype_method_id!='':
				if mask_genotype_method_id==0 or mask_genotype_method_id=='0':	#require mask_genotype_method_id to be null
					if alignment.mask_genotype_method_id is not None:
						continue
				elif alignment.mask_genotype_method_id!=mask_genotype_method_id:
					continue
			if parent_individual_alignment_id is not None and alignment.parent_individual_alignment_id!=parent_individual_alignment_id:
				continue
			if local_realigned is not None and alignment.local_realigned!=local_realigned:
				continue
			if reduce_reads is not None and alignment.reduce_reads!=reduce_reads:
				continue
			if country_id_set:
				if alignment.individual_sequence.individual.site is None:
					sys.stderr.write("Warning: alignment (id=%s, path=%s, %s) has no site.\n"%(alignment.id, alignment.path,\
																			alignment.individual_sequence.individual.code))
					continue
				elif (alignment.individual_sequence.individual.site is None or \
									alignment.individual_sequence.individual.site.country_id not in country_id_set):
					continue
			if tax_id_set:
				if alignment.individual_sequence.individual.tax_id is None:
					sys.stderr.write("Warning: alignment (id=%s, path=%s, %s) has no tax_id.\n"%(alignment.id, alignment.path,\
																	alignment.individual_sequence.individual.code))
					continue
				elif alignment.individual_sequence.individual.tax_id not in tax_id_set:
					continue
			if excludeContaminant and alignment.individual_sequence.is_contaminated:	#2012.9.27
				continue
			if is_contaminated is not None and alignment.individual_sequence.is_contaminated!=is_contaminated: #2013.07.03
				continue
			if excludeTissueIDSet and alignment.individual_sequence.tissue_id in excludeTissueIDSet:	#2012.10.2
				continue
			if completedAlignment is not None:
				if completeAlignmentCheckFunction is None:
					isAlignmentCompleted = self.isThisAlignmentComplete(individual_alignment=alignment, data_dir=data_dir)
				else:
					isAlignmentCompleted = completeAlignmentCheckFunction(individual_alignment=alignment, data_dir=data_dir)
				if completedAlignment==isAlignmentCompleted:
					pass
				else:
					sys.stderr.write("SunsetDB.getAlignments() warning: completeness (%s) of alignment %s, (file=%s, read_group=%s) does not match given(%s). Skip.\n"%\
									(isAlignmentCompleted, alignment.id, alignment.path, alignment.getReadGroup(), completedAlignment))
					continue
			if alignment_method_id is not None and alignment.alignment_method_id!=alignment_method_id:
				continue
			if outdated_index is not None and alignment.outdated_index!=outdated_index:
				continue
			newAlignmentLs.append(alignment)
		if report:
			sys.stderr.write(" kept %s alignments. Done.\n"%(len(newAlignmentLs)))
		return newAlignmentLs
	
	@classmethod
	def getCumulativeAlignmentMedianDepth(cls, alignmentLs=[], defaultSampleAlignmentDepth=10):
		"""
		2012.8.7
		"""
		sys.stderr.write("Getting cumulative median depth of %s alignments ..."%(len(alignmentLs)))
		cumulativeDepth = 0
		for alignment in alignmentLs:
			if alignment and alignment.median_depth is not None:
				medianDepth = alignment.median_depth
			else:
				medianDepth = defaultSampleAlignmentDepth
			cumulativeDepth += medianDepth
		sys.stderr.write("=%s.\n"%(cumulativeDepth))
		return cumulativeDepth
	
	def getAlleleType(self, allele_type_name):
		"""
		2011-2-11
		"""
		db_entry = self.queryTable(AlleleType).filter_by(short_name=allele_type_name).first()
		if not db_entry:
			db_entry = AlleleType(short_name=allele_type_name)
			self.session.add(db_entry)
			self.session.flush()
		return db_entry
	
	def getIndividual(self, code=None, name=None, sex=None, age=None, latitude=None, longitude=None, \
					altitude=None, site=None, site_name=None, site_id=None,\
					city=None, stateprovince=None, country_name=None, \
					collection_date=None, collector=None, \
					tax_id=None, birthdate=None, comment=None, study=None, study_id=None):
		"""
		Examples:
			individual = db_vervet.getIndividual(code=code)
			individual = db_vervet.getIndividual(code=code, site=None, site_name='YBK', country_name='Gambia')
			
		2012.12.6 in the case of creating a new individual, and site is None (to  be created on the fly),
			site_name and country_name must be not null.
			latitude, longitude, altitude is no longer part of table Individual.
		2011-9-8
			add argument comment
		2011-5-6
			add birthdate
		2011-5-5
			take first letter of sex and upper case.
		2011-2-11
			code can't be None
		"""
		db_entry = None
		if not db_entry and code is not None:
			query = self.queryTable(Individual).filter_by(code=code)
			if tax_id:
				query = query.filter_by(tax_id=tax_id)
			db_entry = query.first()
		if not db_entry:
			if sex is None or sex=='?' or sex=='':	#2011-4-29
				sex = None
			elif len(sex)>=1:	#2011-5-5 take the first letter
				sex = sex[0].upper()
			if not site and site_name and country_name:
				site = self.getSite(description=site_name, city=city, \
								stateprovince=stateprovince, country_name=country_name,\
								latitude=latitude, longitude=longitude, altitude=altitude)
			db_entry = Individual(code=code, name=name, sex=sex, age=age, site=site, site_id=site_id,\
						collection_date=collection_date, collector=collector,\
						tax_id=tax_id, birthdate=birthdate,\
						comment=comment, study=study, study_id=study_id)
			self.session.add(db_entry)
			self.session.flush()
		return db_entry
	
	def checkIndividualSequence(self, individual_id=None, sequencer_id=None, sequencer_name=None, \
							sequence_type_name=None, sequence_type_id=None, \
						sequence_format=None, tissue_name=None, tissue_id=None, \
						filtered=0,\
						parent_individual_sequence_id=None,\
						no_of_chromosomes=None, sequence_batch_id=None, version=None, \
						is_contaminated=0, outdated_index=0, returnFirstEntry=True):
		"""
		2013.04.03 bugfix
		2013.3.15
		"""
		if sequencer_id is None and sequencer_name:
			sequencer = self.getSequencer(short_name=sequencer_name)
			if sequencer:
				sequencer_id=sequencer.id
		if sequence_type_id is None and sequence_type_name:
			sequence_type = self.getSequenceType(short_name=sequence_type_name)
			if sequence_type:
				sequence_type_id=sequence_type.id
		
		query = self.queryTable(IndividualSequence).filter_by(individual_id=individual_id).filter_by(filtered=filtered)
		query = query.filter_by(sequencer_id=sequencer_id)
		query = query.filter_by(sequence_type_id=sequence_type_id)
		if sequence_format:
			query = query.filter_by(format=sequence_format)
		query = query.filter_by(is_contaminated=is_contaminated).filter_by(outdated_index=outdated_index)
		
		if tissue_name:
			tissue = self.getTissue(short_name=tissue_name)
			query = query.filter_by(tissue_id=tissue.id)
		elif tissue_id:
			query = query.filter_by(tissue_id=tissue_id)
		else:
			query = query.filter_by(tissue_id=None)
		if parent_individual_sequence_id:
			query = query.filter_by(parent_individual_sequence_id=parent_individual_sequence_id)
		else:
			query = query.filter_by(parent_individual_sequence_id=None)
		
		if sequence_batch_id is not None:
			query = query.filter_by(sequence_batch_id=sequence_batch_id)
		else:
			query = query.filter_by(sequence_batch_id=None)
		
		if version is not None:	#default is 1. so if argument is None, don't query it
			query = query.filter_by(version=version)
		
		if no_of_chromosomes is not None:
			query = query.filter_by(no_of_chromosomes=no_of_chromosomes)
		else:
			query = query.filter_by(no_of_chromosomes=None)
		
		no_of_entries = query.count()
		if no_of_entries>1:
			sys.stderr.write("Error, >1 entries (%s) returned for IndividualSequence \
		(individual_id=%s, filtered=%s, \
		sequencer_id=%s, sequence_type_id=%s, sequence_format=%s, is_contaminated=%s, outdated_index=%s,\
		tissue_id=%s, parent_individual_sequence_id=%s, sequence_batch_id=%s, version=%s, no_of_chromosomes=%s).\n"%\
							(no_of_entries, individual_id, filtered, sequencer_id, sequence_type_id,\
							sequence_format, is_contaminated, outdated_index,\
							tissue_id, parent_individual_sequence_id, sequence_batch_id, version,\
							no_of_chromosomes))
			raise
		if returnFirstEntry:
			db_entry = query.first()
			return db_entry
		else:
			return query
	
	def checkIndividualAlignment(self, individual_code=None, individual_id=None, individual=None, individual_sequence_id=None, \
					path_to_original_alignment=None, sequencer_name='GA', sequencer_id=None, \
					sequence_type_name=None, sequence_type_id=None, sequence_format='fastq', \
					ref_individual_sequence_id=10, \
					alignment_method_name='bwa-short-read', alignment_method_id=None, alignment_method=None,\
					alignment_format='bam', subFolder='individual_alignment', \
					createSymbolicLink=False, individual_sequence_filtered=0, read_group_added=None, data_dir=None, \
					outdated_index=0, mask_genotype_method_id=None, parent_individual_alignment_id=None,\
					individual_sequence_file_raw_id=None, md5sum=None, local_realigned=0, reduce_reads=None):
		"""
		2013.04.11 added argument reduce_reads
		2013.04.05 split out of getAlignment()
		"""
		if data_dir is None:
			data_dir = self.data_dir
		if individual_code:
			individual = self.getIndividual(individual_code)
		elif individual_id:
			individual = self.queryTable(Individual).get(individual_id)
		
		if individual_sequence_id:
			individual_sequence = self.queryTable(IndividualSequence).get(individual_sequence_id)
			if individual is None:
				individual = individual_sequence.individual
		elif individual:
			individual_sequence = self.getIndividualSequence(individual_id=individual.id, sequencer_name=sequencer_name, \
								sequencer_id=sequencer_id, sequence_type_id=sequence_type_id,\
								sequence_type_name=sequence_type_name,\
								sequence_format=sequence_format, filtered=individual_sequence_filtered)
		else:
			sys.stderr.write("Error: not able to get individual_sequence cuz individual_sequence_id=%s; individual_code=%s; individual_id=%s .\n"%\
							(individual_sequence_id, individual_code, individual_id))
			return None
		
		if alignment_method_name:
			alignment_method = self.getAlignmentMethod(alignment_method_name=alignment_method_name)
		elif alignment_method_id:
			alignment_method = self.queryTable(AlignmentMethod).get(alignment_method_id)
		elif alignment_method is None:
			sys.stderr.write("Error: not able to get alignment_method cuz alignment_method_name=%s and alignment_method_id=%s.\n"%\
							(alignment_method_name, alignment_method_id))
			return None
		
		query = self.queryTable(IndividualAlignment).filter_by(ind_seq_id=individual_sequence.id).\
				filter_by(ref_ind_seq_id=ref_individual_sequence_id).\
				filter_by(alignment_method_id=alignment_method.id).filter_by(outdated_index=outdated_index).\
				filter_by(mask_genotype_method_id=mask_genotype_method_id).\
				filter_by(parent_individual_alignment_id=parent_individual_alignment_id).\
				filter_by(individual_sequence_file_raw_id=individual_sequence_file_raw_id)
		
		if alignment_format:
			query = query.filter_by(format=alignment_format)
		if local_realigned is not None:
			query = query.filter_by(local_realigned=local_realigned)
		if reduce_reads is not None:
			query = query.filter_by(reduce_reads=reduce_reads)
		
		no_of_entries = query.count()
		if no_of_entries>1:
			sys.stderr.write("Error, >1 entries (%s) returned for IndividualAlignment \
		(individual_id=%s, ref_individual_sequence_id=%s, \
		sequencer_id=%s, sequence_type_id=%s, sequence_format=%s, filtered=%s, outdated_index=%s,\
		alignment_method_id=%s,mask_genotype_method_id =%s, parent_individual_alignment_id=%s, individual_sequence_file_raw_id=%s,\
		alignment_format=%s, local_realigned=%s, reduce_reads=%s).\n"%\
							(no_of_entries, individual_id, ref_individual_sequence_id, sequencer_id, sequence_type_id,\
							sequence_format, individual_sequence_filtered, outdated_index,\
							alignment_method_id, mask_genotype_method_id, parent_individual_alignment_id, individual_sequence_file_raw_id,\
							alignment_format, local_realigned, reduce_reads))
			raise
		db_entry = query.first()
		return db_entry
	
	def getIndividualSequence(self, individual_id=None, sequencer_id=None, sequencer_name=None, \
							sequence_type_name=None, sequence_type_id=None, \
						sequence_format=None, path_to_original_sequence=None, copy_original_file=False,\
						tissue_name=None, tissue_id=None, \
						coverage=None,\
						subFolder=None, quality_score_format="Standard", filtered=0,\
						parent_individual_sequence_id=None,\
						read_count=None, no_of_chromosomes=None, sequence_batch_id=None, version=None, data_dir=None,\
						is_contaminated=0, outdated_index=0, comment=None):
		"""
		2013.3.15 added is_contaminated, outdated_index
		2013.3.13 read_count, no_of_chromosomes, sequencer_id, sequence_type_id, sequence_batch_id, version
		2012.6.3
			columns that are None become part of the db query to see if entry is in db already
		2012.2.24
			add argument data_dir
		2012.2.10
			path_to_original_sequence is only given when you want to copy the file to db storage.
			add argument parent_individual_sequence_id
		2011-8-30
			add argument filtered
		2011-8-18
			add argument quality_score_format, default to "Standard"
		2011-8-3
			the path field is now considered a folder (rather than a file).
		2011-5-7
			subFolder is the name of the folder in self.data_dir that is used to hold the sequence files.
		"""
		if not data_dir:
			data_dir = self.data_dir
		if not subFolder:
			subFolder = IndividualSequence.folderName
		
		if sequencer_id is None and sequencer_name:
			sequencer = self.getSequencer(short_name=sequencer_name)
			if sequencer:
				sequencer_id=sequencer.id
		if sequence_type_id is None and sequence_type_name:
			sequence_type = self.getSequenceType(short_name=sequence_type_name)
			if sequence_type:
				sequence_type_id=sequence_type.id
		
		db_entry = self.checkIndividualSequence(individual_id=individual_id, sequencer_id=sequencer_id, \
						sequencer_name=sequencer_name, sequence_type_name=sequence_type_name, \
						sequence_type_id=sequence_type_id, sequence_format=sequence_format, tissue_name=tissue_name,\
						tissue_id=tissue_id, filtered=filtered, \
						parent_individual_sequence_id=parent_individual_sequence_id, \
						no_of_chromosomes=no_of_chromosomes, sequence_batch_id=sequence_batch_id, \
						version=version, is_contaminated=is_contaminated, outdated_index=outdated_index)
		if not db_entry:
			tissue = self.getTissue(db_entry_id=tissue_id, short_name=tissue_name)
			db_entry = IndividualSequence(individual_id=individual_id, sequencer_id=sequencer_id, \
										sequence_type_id=sequence_type_id,\
									format=sequence_format, original_path=path_to_original_sequence,\
									tissue=tissue, coverage=coverage, \
									quality_score_format=quality_score_format, filtered=filtered,\
									parent_individual_sequence_id=parent_individual_sequence_id, \
									read_count=read_count, no_of_chromosomes=no_of_chromosomes,\
									sequence_batch_id=sequence_batch_id, version=version, \
									is_contaminated=is_contaminated, outdated_index=outdated_index,\
									comment=comment)
			#to make db_entry.id valid
			self.session.add(db_entry)
			self.session.flush()
			
			dst_relative_path = db_entry.constructRelativePathForIndividualSequence(subFolder=subFolder)
			
			#update its path in db to the relative path
			db_entry.path = dst_relative_path
			
			dst_abs_path = os.path.join(data_dir, dst_relative_path)
			if copy_original_file and path_to_original_sequence and (os.path.isfile(path_to_original_sequence) or os.path.isdir(path_to_original_sequence)):
				dst_dir = os.path.join(data_dir, subFolder)
				if not os.path.isdir(dst_dir):	#the upper directory has to be created at this moment.
					commandline = 'mkdir -p %s'%(dst_dir)
					return_data = runLocalCommand(commandline, report_stderr=True, report_stdout=True)
				if not os.path.isdir(dst_abs_path):	#2011-8-3 create the directory to host all sequences.
					commandline = 'mkdir -p %s'%(dst_abs_path)
					return_data = runLocalCommand(commandline, report_stderr=True, report_stdout=True)
				commandline = 'cp -r %s %s'%(path_to_original_sequence, dst_abs_path)
				return_data = runLocalCommand(commandline, report_stderr=True, report_stdout=True)
			self.session.add(db_entry)
			self.session.flush()
		return db_entry

	def getIndividualDBEntryViaDBID(self, db_id=None):
		"""
		2012.9.6
			has a cache inside
		"""
		if db_id in self.dbID2monkeyDBEntry:
			return self.dbID2monkeyDBEntry.get(db_id)
		else:
			monkey = self.queryTable(Individual).get(db_id)
			self.dbID2monkeyDBEntry[db_id] = monkey
			return monkey
	
	def copyParentIndividualSequence(self, parent_individual_sequence=None, parent_individual_sequence_id=None,\
									subFolder='individual_sequence', quality_score_format='Standard', filtered=1,\
									data_dir=None):
		"""
		2012.6.10
			call getIndividualSequence to construct individual_sequence, rather than construct it from here.
		2012.2.10
			constructRelativePathForIndividualSequence is now moved to Table IndividualSequence
			add argument filtered
		2011-9-1
			add quality_score_format
		2011-8-11
			make a copy of parent_individual_sequence_id individual_sequence entity
		"""
		if parent_individual_sequence is None:
			parent_individual_sequence = self.queryTable(IndividualSequence).get(parent_individual_sequence_id)
		
		pis = parent_individual_sequence
		
		individual = self.queryTable(Individual).get(pis.individual_id)
		
		individual_sequence = self.getIndividualSequence(individual_id=parent_individual_sequence.individual_id, \
						sequencer_id=parent_individual_sequence.sequencer.id, \
						sequence_type_id=parent_individual_sequence.sequence_type.id,\
						sequence_format=parent_individual_sequence.format, path_to_original_sequence=None, \
						tissue_id=getattr(parent_individual_sequence.tissue, 'id', None), coverage=None,\
						quality_score_format=quality_score_format, filtered=filtered,\
						parent_individual_sequence_id=parent_individual_sequence.id, \
						data_dir=data_dir, subFolder=subFolder)
		return individual_sequence
	
	def copyParentIndividualSequenceFile(self, parent_individual_sequence_file=None, parent_individual_sequence_file_id=None,\
									individual_sequence_id=None,\
									quality_score_format='Standard', filtered=1):
		"""
		2012.2.14
			call self.getIndividualSequenceFile() instead
		2012.2.10
		"""
		if parent_individual_sequence_file is None:
			if parent_individual_sequence_file_id is not None:
				parent_individual_sequence_file = self.queryTable(IndividualSequenceFile).get(parent_individual_sequence_file_id)
			else:
				sys.stderr.write("Warning: parent individual_sequence_id is None. No IndividualSequenceFile instance to be created.\n")
				return None
		
		
		parent = parent_individual_sequence_file
		
		
		db_entry = self.getIndividualSequenceFile(individual_sequence_id, library=parent.library, mate_id=parent.mate_id, \
									split_order=parent.split_order, format=parent.format,\
									filtered=filtered, parent_individual_sequence_file_id=parent.id, \
									individual_sequence_file_raw_id=parent.individual_sequence_file_raw_id,\
									quality_score_format=quality_score_format)
		return db_entry
	
	def copyParentIndividualAlignment(self, parent_individual_alignment=None, \
									parent_individual_alignment_id=None, mask_genotype_method=None, \
									mask_genotype_method_id=None, data_dir=None, local_realigned=0, read_group=None,\
									reduce_reads=None):
		"""
		Examples:
			in AlignmentReduceReadsWorkflow.py
			new_individual_alignment = self.db.copyParentIndividualAlignment(parent_individual_alignment_id=individual_alignment.id,\
										data_dir=self.data_dir, local_realigned=individual_alignment.local_realigned,\
										reduce_reads=1)
		
		2013.04.11 added argument reduce_reads
		2013.04.09 added argument read_group
		2013.04.05 added argument local_realigned
		2012.7.28
		"""
		if parent_individual_alignment is None and parent_individual_alignment_id:
			parent_individual_alignment = self.queryTable(IndividualAlignment).get(parent_individual_alignment_id)
		if mask_genotype_method is None and mask_genotype_method_id:
			mask_genotype_method = self.queryTable(GenotypeMethod).get(mask_genotype_method_id)
		
		individual_sequence = parent_individual_alignment.individual_sequence
		ref_sequence = parent_individual_alignment.ref_sequence
		alignment_method = parent_individual_alignment.alignment_method
		individual_alignment = self.getAlignment(individual_code=individual_sequence.individual.code, \
								individual_sequence_id=individual_sequence.id,\
					path_to_original_alignment=None, sequencer_id=individual_sequence.sequencer_id, \
					sequence_type_id=individual_sequence.sequence_type_id, sequence_format=individual_sequence.format, \
					ref_individual_sequence_id=ref_sequence.id, \
					alignment_method_name=alignment_method.short_name, alignment_format=parent_individual_alignment.format,\
					individual_sequence_filtered=individual_sequence.filtered, read_group_added=1,
					data_dir=data_dir, outdated_index=0, mask_genotype_method_id=getattr(mask_genotype_method, 'id', None),\
					parent_individual_alignment_id=parent_individual_alignment.id,\
					individual_sequence_file_raw_id=parent_individual_alignment.individual_sequence_file_raw_id,\
					local_realigned=local_realigned, read_group=read_group, reduce_reads=reduce_reads)
					#read-group addition is part of pipeline
		#if not individual_alignment.path:
		#	individual_alignment.path = individual_alignment.constructRelativePath()
		#	session.add(individual_alignment)
		#	session.flush()
		return individual_alignment
	
	def getIndividualSequenceFileRaw(self, individual_sequence_id=None, library=None, md5sum=None, path=None, \
									original_path=None, mate_id=None, file_size=None):
		"""
		2013.04.03 argument original_path
		2012.7.12 add argument file_size
		2012.4.30
			add argument mate_id
		2012.2.14
		"""
		#query first
		query = self.queryTable(IndividualSequenceFileRaw).filter_by(individual_sequence_id=individual_sequence_id)
		if library:
			query = query.filter_by(library=library)
		if md5sum:
			query = query.filter_by(md5sum=md5sum)
		if path:
			path = os.path.realpath(path)	#get the realpath
			query = query.filter_by(path=path)
		if original_path:
			original_path = os.path.realpath(original_path)	#get the realpath
			query = query.filter_by(original_path=original_path)
		
		if mate_id:
			query = query.filter_by(mate_id=mate_id)
		db_entry = query.first()
		if not db_entry:
			if file_size is None:	#2012.7.12
				if path and os.path.isfile(path):
					file_size = utils.getFileOrFolderSize(path)
				elif original_path and os.path.isfile(original_path):
					file_size = utils.getFileOrFolderSize(path)
			db_entry = IndividualSequenceFileRaw(individual_sequence_id=individual_sequence_id, library=library, md5sum=md5sum, \
										path=path, mate_id=mate_id, file_size=file_size, original_path=original_path)
			self.session.add(db_entry)
			self.session.flush()
		return db_entry
	
	def registerOriginalSequenceFileToDB(self, original_sequence_filepath=None, library=None, 
										 individual_sequence_id=None, mate_id=None, md5sum=None):
		"""
		20190206 moved from RegisterIndividualSequence2DB.py
		20130403 original_sequence_filepath passed to original_path
		2012.4.30
			add mate_id
		2012.1.27
			1. run md5sum
			2. check if it already exists in db
				if not, add it into db
				if yes, exit the program. no more work.
		"""
		sys.stderr.write("Register the original sequence file %s into db ..."%(original_sequence_filepath))
		if not md5sum:
			md5sum = utils.get_md5sum(original_sequence_filepath)
		db_entry = self.queryTable(IndividualSequenceFileRaw).filter_by(md5sum=md5sum).first()
		if db_entry:
			sys.stderr.write("Warning: another file %s with the identical md5sum %s (library=%s) as this file %s is already in db.\n"%\
								(db_entry.path, md5sum, library, original_sequence_filepath))
		else:
			file_size = utils.getFileOrFolderSize(original_sequence_filepath)
			db_entry = self.getIndividualSequenceFileRaw(individual_sequence_id=individual_sequence_id, library=library, md5sum=md5sum,
											original_path=original_sequence_filepath, mate_id=mate_id, file_size=file_size)
		mate_id2split_order_ls = {}
		for individual_sequence_file in db_entry.individual_sequence_file_ls:
			mate_id = individual_sequence_file.mate_id
			if mate_id is None:	#sequence entries without mate_id are just from one mate.
				mate_id = 1
			if mate_id not in mate_id2split_order_ls:
				mate_id2split_order_ls[mate_id] = []
			mate_id2split_order_ls[mate_id].append(individual_sequence_file.split_order)
		if len(mate_id2split_order_ls)>2:
			sys.stderr.write("Error: db sequence files spawned from file %s (md5sum=%s) form %s(>2) mates.\
				Unless this bam file contains reads from >2 mates, reads from this file should be stored in db already.\n"%\
				(db_entry.path, md5sum, len(mate_id2split_order_ls)))
			sys.exit(4)
		return db_entry
	
	def getIndividualSequenceFile(self, individual_sequence_id, library=None, mate_id=None, split_order=None, format=None,\
							filtered=0, parent_individual_sequence_file_id=None, individual_sequence_file_raw_id=None,\
							quality_score_format='Standard'):
		"""
		2012.5.2
			add all filter_by() straight to self.queryTable(IndividualSequenceFile)
		2012.2.14
		"""
		#query first
		query = self.queryTable(IndividualSequenceFile).filter_by(individual_sequence_id=individual_sequence_id).filter_by(library=library).\
					filter_by(mate_id=mate_id).filter_by(split_order=split_order).filter_by(format=format).filter_by(filtered=filtered).\
					filter_by(parent_individual_sequence_file_id=parent_individual_sequence_file_id).\
					filter_by(individual_sequence_file_raw_id=individual_sequence_file_raw_id)
		db_entry = query.first()
		if not db_entry:
			db_entry = IndividualSequenceFile(individual_sequence_id=individual_sequence_id,\
							library=library, mate_id=mate_id, split_order=split_order,\
							format=format, filtered=filtered, \
							parent_individual_sequence_file_id=parent_individual_sequence_file_id,\
							individual_sequence_file_raw_id=individual_sequence_file_raw_id,\
							quality_score_format=quality_score_format)
			self.session.add(db_entry)
			self.session.flush()
		return db_entry
	
	def getRelationshipType(self, relationship_type_name=None):
		"""
		2011-5-5
			fetch (and add) a RelationshipType entry
		"""
		query = self.queryTable(RelationshipType).filter_by(short_name=relationship_type_name)
		db_entry = query.first()
		if not db_entry:
			db_entry = RelationshipType(short_name=relationship_type_name)
			self.session.add(db_entry)
			self.session.flush()
		return db_entry
	
	def checkSpecificRelationOfIndividual2(self, individual2=None, relationship_type_name=None):
		"""
		2012.1.23
			check to see if specific relationship of individual2 exists in Ind2Ind already.
			This is to detect errors where one individual has >1 father/mother in the database.
		"""
		relationship_type = self.getRelationshipType(relationship_type_name)
		query = self.queryTable(Ind2Ind).filter_by(individual2_id=individual2.id).filter_by(relationship_type_id=relationship_type.id)
		db_entry = query.first()
		if not db_entry:
			return None
		else:
			return db_entry
	
	def checkIndividualAlignmentConsensusSequence(self, individual_alignment_id=None, minDP=None, maxDP=None, minBaseQ=None, minMapQ=None,\
						minRMSMapQ=None, minDistanceToIndel=None, no_of_chromosomes=None, **keywords):
		"""
		2013.2.8 check whether one IndividualAlignmentConsensusSequence is in db or not.
		"""
		
		query = self.queryTable(IndividualAlignmentConsensusSequence).filter_by(individual_alignment_id=individual_alignment_id).\
			filter_by(minDP=minDP).filter_by(maxDP=maxDP)
		if minBaseQ is not None:
			query = query.filter_by(minBaseQ=minBaseQ)
		if minMapQ is not None:
			query = query.filter_by(minMapQ=minMapQ)
		if minRMSMapQ:
			query = query.filter_by(minRMSMapQ=minRMSMapQ)
		if minDistanceToIndel:
			query = query.filter_by(minDistanceToIndel=minDistanceToIndel)
		if no_of_chromosomes:
			query = query.filter_by(no_of_chromosomes=no_of_chromosomes)
		
		db_entry = query.first()
		if db_entry:
			return db_entry
		else:
			return None
	
	def getIndividualAlignmentConsensusSequence(self, individual_alignment_id=None, format=None, minDP=None, maxDP=None, minBaseQ=None, \
											minMapQ=None,\
						minRMSMapQ=None, minDistanceToIndel=None, no_of_chromosomes=None,no_of_bases=None, \
						original_path=None, data_dir=None, **keywords):
		"""
		2013.2.8 get one IndividualAlignmentConsensusSequence from db
		"""
		db_entry = self.checkIndividualAlignmentConsensusSequence(individual_alignment_id=individual_alignment_id, minDP=minDP, \
									maxDP=maxDP, minBaseQ=minBaseQ, minMapQ=minMapQ,\
									minRMSMapQ=minRMSMapQ, minDistanceToIndel=minDistanceToIndel, no_of_chromosomes=no_of_chromosomes)
		
		if not db_entry:
			if original_path:
				original_path = os.path.abspath(original_path)
			db_entry = IndividualAlignmentConsensusSequence(individual_alignment_id=individual_alignment_id, \
								format=format,\
								minDP=minDP, maxDP=maxDP,\
								minBaseQ=minBaseQ, minMapQ=minMapQ, minRMSMapQ=minRMSMapQ, minDistanceToIndel=minDistanceToIndel,\
								no_of_chromosomes=no_of_chromosomes, no_of_bases=no_of_bases,
								original_path=original_path,  **keywords)
			self.session.add(db_entry)
			self.session.flush()
			if db_entry.path and db_entry.file_size is None:
				self.updateDBEntryPathFileSize(db_entry=db_entry, data_dir=data_dir)
			if db_entry.path and db_entry.md5sum is None:
				self.updateDBEntryMD5SUM(db_entry=db_entry, data_dir=data_dir)
		return db_entry
	
	def getInd2Ind(self, individual1=None, individual2=None, relationship_type_name=None):
		"""
		2011-5-5
			fetch (and add) a Ind2Ind entry
		"""
		relationship_type = self.getRelationshipType(relationship_type_name)
		query = self.queryTable(Ind2Ind).filter_by(individual1_id=individual1.id).filter_by(individual2_id=individual2.id).\
			filter_by(relationship_type_id=relationship_type.id)
		db_entry = query.first()
		if not db_entry:
			db_entry = Ind2Ind(individual1=individual1, individual2=individual2, relationship_type=relationship_type)
			self.session.add(db_entry)
			self.session.flush()
		return db_entry
	
	def checkAlignmentDepthIntervalMethod(self, db_entry_id=None, short_name=None, \
					ref_ind_seq_id=None, parent_id=None, \
					min_segment_length=None,**keywords):
		"""
		2013.09.02
		"""
		if db_entry_id:
			db_entry = self.checkIfEntryInTable(TableClass=AlignmentDepthIntervalMethod, id=db_entry_id)
		else:
			db_entry = None
		if not db_entry:
			query = self.queryTable(AlignmentDepthIntervalMethod).filter_by(short_name=short_name)
			if ref_ind_seq_id:
				query = query.filter_by(ref_ind_seq_id=ref_ind_seq_id)
			if parent_id:
				query = query.filter_by(parent_id=parent_id)
			if min_segment_length:
				query = query.filter_by(min_segment_length=min_segment_length)
			db_entry = query.first()
		if db_entry:
			return db_entry
		else:
			return None
	
	def getAlignmentDepthIntervalMethod(self, short_name=None, description=None, ref_ind_seq_id=None, \
						individualAlignmentLs=None, parent_db_entry=None, parent_id=None, \
						no_of_alignments=None, no_of_intervals=None, sum_median_depth=None, sum_mean_depth=None,\
						min_segment_length=None,\
						data_dir=None):
		"""
		2013.09.02 added min_segment_length
		2013.08.23
			examples:
				alignmentDepthIntervalMethod = self.db_vervet.getAlignmentDepthIntervalMethod(short_name="725VRCAlignmentDepthInterval", \
														individualAlignmentLs=individualAlignmentLs)
		"""
		if not short_name:
			sys.stderr.write("Error: short_name (%s) is empty.\n"%(short_name))
			return None
		
		if ref_ind_seq_id is None and individualAlignmentLs:
			firstAlignment = individualAlignmentLs[0]
			ref_ind_seq_id =firstAlignment.ref_ind_seq_id
		
		if parent_db_entry and parent_id is None:
			parent_id = parent_db_entry.id
		
		db_entry = self.checkAlignmentDepthIntervalMethod(db_entry_id=None, short_name=short_name, \
					ref_ind_seq_id=ref_ind_seq_id, parent_id=parent_id, \
					min_segment_length=min_segment_length)
		if not db_entry:
			if individualAlignmentLs:
				if no_of_alignments is None:
					no_of_alignments = len(individualAlignmentLs)
				if sum_median_depth is None:
					sum_median_depth = sum([db_entry.median_depth for db_entry in individualAlignmentLs])
				if sum_mean_depth is None:
					sum_mean_depth = sum([db_entry.mean_depth for db_entry in individualAlignmentLs])
			
			db_entry = AlignmentDepthIntervalMethod(short_name=short_name, description=description, ref_ind_seq_id=ref_ind_seq_id,\
								parent_id=parent_id, no_of_alignments=no_of_alignments, no_of_intervals=no_of_intervals, \
								sum_median_depth=sum_median_depth, sum_mean_depth=sum_mean_depth, min_segment_length=min_segment_length)
			db_entry.individual_alignment_ls = individualAlignmentLs
			self.session.add(db_entry)
			self.session.flush()
			if not data_dir:
				data_dir = self.data_dir
			
			db_entry.path = db_entry.constructRelativePath()
			folderAbsPath = os.path.join(data_dir, db_entry.path)
			
			if not os.path.isdir(folderAbsPath):	#the upper directory has to be created at this moment.
				commandline = 'mkdir -p %s'%(folderAbsPath)
				return_data = runLocalCommand(commandline, report_stderr=True, report_stdout=True)
			self.session.add(db_entry)
			self.session.flush()
		return db_entry
	
	def checkAlignmentDepthIntervalFile(self, db_entry_id=None, alignment_depth_interval_method_id=None, \
					chromosome=None, no_of_chromosomes=None, format=None,\
					**keywords):
		"""
		2013.08.23
		"""
		if db_entry_id:
			db_entry = self.checkIfEntryInTable(TableClass=AlignmentDepthIntervalFile, id=db_entry_id)
		else:
			db_entry = None
		if not db_entry:
			query = self.queryTable(AlignmentDepthIntervalFile).filter_by(alignment_depth_interval_method_id=alignment_depth_interval_method_id)
			if format:
				query = query.filter_by(format=format)
			if no_of_chromosomes:
				query = query.filter_by(no_of_chromosomes=no_of_chromosomes)
			if chromosome:
				query = query.filter_by(chromosome=chromosome)
			
			db_entry = query.first()
		if db_entry:
			return db_entry
		else:
			return None
	
	def getAlignmentDepthIntervalFile(self, alignment_depth_interval_method=None, alignment_depth_interval_method_id=None,  \
					path=None, file_size=None, \
					chromosome=None, chromosome_size=None, no_of_chromosomes=None, no_of_intervals=None,\
					format="tsv",\
					mean_interval_value=None, median_interval_value=None, \
					min_interval_value=None, max_interval_value=None, \
					min_interval_length=None, max_interval_length=None, median_interval_length=None,\
					md5sum=None, original_path=None, data_dir=None):
		"""
		2013.08.30 added min_interval_length, max_interval_length, median_interval_length
		2013.8.23
			add argument no_of_chromosomes
		"""
		if alignment_depth_interval_method_id is None and alignment_depth_interval_method.id:
			alignment_depth_interval_method_id = alignment_depth_interval_method.id
			
		
		db_entry = self.checkAlignmentDepthIntervalFile(db_entry_id=None, \
									alignment_depth_interval_method_id=alignment_depth_interval_method_id, \
									chromosome=chromosome, \
									no_of_chromosomes=no_of_chromosomes, format=format)
		if original_path:
			original_path = os.path.realpath(original_path)
		
		if not db_entry:
			db_entry = AlignmentDepthIntervalFile(path=path, \
							alignment_depth_interval_method_id=alignment_depth_interval_method_id, \
							chromosome=chromosome, \
							chromosome_size=chromosome_size, no_of_chromosomes=no_of_chromosomes, no_of_intervals=no_of_intervals,\
							file_size=file_size, format=format,\
							mean_interval_value=mean_interval_value, median_interval_value=median_interval_value, \
							min_interval_value=min_interval_value, max_interval_value=max_interval_value,\
							min_interval_length=min_interval_length, max_interval_length=max_interval_length, \
							median_interval_length=median_interval_length,\
							md5sum=md5sum, original_path=original_path)
			if not data_dir:
				data_dir = self.data_dir
			if db_entry.file_size is None and db_entry.path:
				db_entry.file_size = db_entry.getFileSize(data_dir=data_dir)
			self.session.add(db_entry)
			self.session.flush()
			
			folderAbsPath = os.path.join(data_dir, db_entry.alignment_depth_interval_method.constructRelativePath())
			
			if not os.path.isdir(folderAbsPath):	#the upper directory has to be created at this moment.
				commandline = 'mkdir -p %s'%(folderAbsPath)
				return_data = runLocalCommand(commandline, report_stderr=True, report_stdout=True)
		return db_entry
	
	def getGenotypeMethod(self, short_name=None, description=None, ref_ind_seq_id=None, \
						individualAlignmentLs=None, parent_genotype_method=None, parent_id=None, \
						min_depth=None, max_depth=None, max_missing_rate=None, min_maf=None,\
						no_of_individuals=None, no_of_loci=None, is_phased=0, data_dir=None, \
						subFolder='genotype_file'):
		"""
		2012.7.12
			examples:
				genotypeMethod = self.db_vervet.getGenotypeMethod(short_name=self.genotypeMethodShortName, \
														individualAlignmentLs=individualAlignmentLs)
		"""
		if not short_name:
			sys.stderr.write("Error: short_name (%s) is empty.\n"%(short_name))
			return None
		
		if ref_ind_seq_id is None and individualAlignmentLs:
			firstAlignment = individualAlignmentLs[0]
			ref_ind_seq_id =firstAlignment.ref_ind_seq_id
		
		if parent_genotype_method and parent_id is None:
			parent_id = parent_genotype_method.id
		
		query = self.queryTable(GenotypeMethod).filter_by(short_name=short_name)
		
		if ref_ind_seq_id:
			query = query.filter_by(ref_ind_seq_id=ref_ind_seq_id)
		if parent_id:
			query = query.filter_by(parent_id=parent_id)
		
		
		db_entry = query.first()
		if not db_entry:
			db_entry = GenotypeMethod(short_name=short_name, description=description, ref_ind_seq_id=ref_ind_seq_id,\
								parent_id=parent_id, min_depth=min_depth, max_depth=max_depth, \
								max_missing_rate=max_missing_rate, min_maf=min_maf,\
								no_of_individuals=no_of_individuals, no_of_loci=no_of_loci, is_phased=is_phased)
			db_entry.individual_alignment_ls = individualAlignmentLs
			self.session.add(db_entry)
			self.session.flush()
			if not data_dir:
				data_dir = self.data_dir
			
			db_entry.path = db_entry.constructRelativePath(subFolder=subFolder)
			folderAbsPath = os.path.join(data_dir, db_entry.path)
			
			if not os.path.isdir(folderAbsPath):	#the upper directory has to be created at this moment.
				commandline = 'mkdir -p %s'%(folderAbsPath)
				return_data = runLocalCommand(commandline, report_stderr=True, report_stdout=True)
			self.session.add(db_entry)
			self.session.flush()
		return db_entry
		
	def getGenotypeFile(self, genotype_method=None, genotype_method_id=None, chromosome=None, format=None, path=None, file_size=None, \
					no_of_individuals=None, no_of_loci=None, md5sum=None, original_path=None, data_dir=None, subFolder='genotype_file',\
					no_of_chromosomes=1):
		"""
		2012.8.30
			add argument no_of_chromosomes
		2012.7.13
		"""
		if genotype_method_id is None and genotype_method.id:
			genotype_method_id = genotype_method.id
			
		query = self.queryTable(GenotypeFile).filter_by(genotype_method_id=genotype_method_id)

		if format:
			query = query.filter_by(format=format)
		if md5sum:	#2012.8.6 if format is not given, then use md5sum as the sole unique identifier 
			query = query.filter_by(md5sum=md5sum)
		if no_of_chromosomes:
			query = query.filter_by(no_of_chromosomes=no_of_chromosomes)
		if chromosome:
			query = query.filter_by(chromosome=chromosome)
		
		db_entry = query.first()
		if original_path:
			original_path = os.path.realpath(original_path)
		
		if not db_entry:
			
			db_entry = GenotypeFile(format=format, path=path, no_of_individuals=no_of_individuals, no_of_loci=no_of_loci, \
							chromosome=chromosome, md5sum=md5sum, file_size=file_size, \
							original_path=original_path, genotype_method_id=genotype_method_id, no_of_chromosomes=no_of_chromosomes)
			if not data_dir:
				data_dir = self.data_dir
			if db_entry.file_size is None and db_entry.path:
				db_entry.file_size = db_entry.getFileSize(data_dir=data_dir)
			self.session.add(db_entry)
			self.session.flush()
			
			folderAbsPath = os.path.join(data_dir, db_entry.genotype_method.constructRelativePath(subFolder=subFolder))
			
			if not os.path.isdir(folderAbsPath):	#the upper directory has to be created at this moment.
				commandline = 'mkdir -p %s'%(folderAbsPath)
				return_data = runLocalCommand(commandline, report_stderr=True, report_stdout=True)
			
		return db_entry
	
	def getLocus(self, chr=None, start=None, stop=None, ref_seq=None, alt_seq=None, ref_ind_seq_id=None, \
				locus_type_id=None, locus_type=None, ref_seq_id=None, alt_seq_id=None):
		"""
		2012.5.2
			add ref_ind_seq_id & locus_type_id, locus_type
		2011-2-11
		
		"""
		if ref_seq_id is None and ref_seq is not None:
			ref_seq_id = ref_seq.id
		if alt_seq_id is None and alt_seq is not None:
			alt_seq_id = alt_seq.id
		if locus_type_id is None and locus_type is not None:
			locus_type_id = locus_type.id
		
		query = self.queryTable(Locus).filter_by(chromosome=chr).filter_by(start=start).filter_by(stop=stop).filter_by(ref_seq_id=ref_seq_id).\
					filter_by(alt_seq_id=alt_seq_id).filter_by(locus_type_id=locus_type_id).\
					filter_by(ref_ind_seq_id=ref_ind_seq_id)
					
		db_entry = query.first()
		if not db_entry:
			db_entry = Locus(chromosome=chr, start=start, stop=stop, ref_ind_seq_id=ref_ind_seq_id, locus_type_id=locus_type_id,\
							ref_seq_id=ref_seq_id, alt_seq_id=alt_seq_id)
			if locus_type:
				db_entry.locus_type = locus_type
			if ref_seq:
				db_entry.ref_seq = ref_seq
			if alt_seq:
				db_entry.alt_seq = alt_seq
			self.session.add(db_entry)
			self.session.flush()
		return db_entry
	
	def getLocusContext(self, locus_id=None, gene_id=None, gene_strand=None, disp_pos=None, overlap_length=None,\
					overlap_fraction_in_locus=None, overlap_fraction_in_gene=None):
		"""
		2012.5.14
		"""
		locus_context = self.queryTable(LocusContext).filter_by(locus_id=locus_id).filter_by(gene_id=gene_id).first()
		
		if not locus_context:
			locus_context = LocusContext(locus_id=locus_id, gene_id = gene_id, \
								gene_strand=gene_strand, disp_pos=disp_pos, \
								overlap_length=overlap_length, \
								overlap_fraction_in_locus=overlap_fraction_in_locus, \
								overlap_fraction_in_gene=overlap_fraction_in_gene)
			self.session.add(locus_context)
		return locus_context
	
	def getLocusAnnotation(self, locus_id=None, locus_context_id=None, locus_context=None, gene_id=None, gene_commentary_id=None, \
						gene_segment_id=None, locus_annotation_type =None, locus_annotation_type_id=None,\
						which_exon_or_intron = None, pos_within_codon=None, which_codon=None, \
						label=None, utr_number=None, cds_number=None, intron_number=None,\
						exon_number=None, overlap_length=None, overlap_fraction_in_locus=None, overlap_fraction_in_gene=None,\
						comment=None):
		"""
		2012.5.14
		"""
		if locus_context_id is None and locus_context:
			locus_context_id = locus_context.id
		if locus_annotation_type_id is None and locus_annotation_type:
			locus_annotation_type_id = locus_annotation_type.id
		locus_annotation = self.queryTable(LocusAnnotation).filter_by(locus_id=locus_id).filter_by(locus_context_id=locus_context.id).\
								filter_by(gene_id=gene_id).\
								filter_by(gene_commentary_id= gene_commentary_id).\
								filter_by(gene_segment_id= gene_segment_id).\
								filter_by(locus_annotation_type_id=locus_annotation_type_id).first()
		if not locus_annotation:
			locus_annotation = LocusAnnotation(locus_id=locus_id, locus_context_id=locus_context_id, gene_id=gene_id,\
							gene_commentary_id = gene_commentary_id, \
							gene_segment_id=gene_segment_id, \
							locus_annotation_type=locus_annotation_type, locus_annotation_type_id=locus_annotation_type_id,\
							which_exon_or_intron=which_exon_or_intron, pos_within_codon=pos_within_codon, \
							which_codon=which_codon, label=label, \
							utr_number = utr_number, cds_number = cds_number, \
							intron_number = intron_number, exon_number = exon_number,\
							overlap_length=overlap_length, \
							overlap_fraction_in_locus=overlap_fraction_in_locus, overlap_fraction_in_gene=overlap_fraction_in_gene,\
							comment=comment)
			locus_annotation.locus_context = locus_context
			self.session.add(locus_annotation)
		return locus_annotation
	
	def getLocusAnnotationType(self, locus_annotation_type_short_name=None):
		"""
		2012.5.16
		
		"""
		ty = self.queryTable(LocusAnnotation).filter_by(short_name=locus_annotation_type_short_name).first()
		if not ty:
			ty = LocusAnnotationType(short_name=locus_annotation_type_short_name)
			self.session.add(ty)
			self.session.flush()
		return ty
	
	def getSequencer(self, short_name=None, seq_center_id=None):
		"""
		2013.3.13
		
		"""
		db_entry = self.queryTable(Sequencer).filter_by(short_name=short_name).first()
		if not db_entry and short_name is not None:
			db_entry = Sequencer(short_name=short_name, seq_center_id=seq_center_id)
			self.session.add(db_entry)
			self.session.flush()
		return db_entry
	
	def getSequenceType(self, short_name=None, id=None, read_length_mean=None, paired_end=0, \
					insert_size_mean=None, insert_size_variance=None, per_base_error_rate=None, **keywords):
		"""
		2013.3.13	#add column id
		"""
		db_entry = self.checkIfEntryInTable(TableClass=SequenceType, short_name=short_name, id=id)
		if not db_entry:
			db_entry = SequenceType(short_name=short_name, read_length_mean=read_length_mean,\
								paired_end=paired_end, insert_size_mean=insert_size_mean, \
								insert_size_variance=insert_size_variance, per_base_error_rate=per_base_error_rate)
			self.session.add(db_entry)
			self.session.flush()
		return db_entry
	
	def getAlleleSequence(self, sequence=None, comment=None):
		"""
		2012.5.2
		
		"""
		db_entry = self.queryTable(AlleleSequence).filter_by(sequence=sequence).first()
		if not db_entry:
			db_entry = AlleleSequence(sequence=sequence, comment=comment)
			self.session.add(db_entry)
			self.session.flush()
		return db_entry
	
	def getSequenceBatch(self, short_name=None, description=None):
		"""
		2012.7.5
		"""
		db_entry = self.queryTable(SequenceBatch).filter_by(short_name=short_name).first()
		if not db_entry:
			db_entry = SequenceBatch(short_name=short_name, description=description)
			self.session.add(db_entry)
			self.session.flush()
		return db_entry
	
	def getSite(self, short_name=None, description=None, city=None, stateprovince=None, region=None, country_name=None, \
			latitude=None, longitude=None, altitude=None, study_id=None):
		"""
		20170412 added study_id
		20121206 added argument short_name
		20110428
		"""
		if country_name:
			country = self.getCountry(country_name=country_name)
		else:
			country = None
		if short_name :
			short_name = short_name
		elif description:
			short_name = description
		elif city:
			short_name = city
		else:
			short_name = None
		query = self.queryTable(Site)
		if country:
			query = query.filter_by(country_id=country.id)
		if short_name:
			query = query.filter_by(short_name=short_name)
		if city:
			query = query.filter_by(city=city)
		if stateprovince:
			query = query.filter_by(stateprovince=stateprovince)
		if region:
			query = query.filter_by(region=region)
		if latitude:
			query = query.filter_by(latitude=latitude)
		if longitude:
			query = query.filter_by(longitude=longitude)
		if altitude:
			query = query.filter_by(altitude=altitude)
		if study_id:
			query = query.filter_by(study_id=study_id)
		db_entry = query.first()
		if not db_entry:
			db_entry = Site(short_name=short_name, description=description, city=city, stateprovince=stateprovince, \
						region=region, country=country, \
						latitude=latitude, longitude=longitude, altitude=altitude, study_id=study_id)
			self.session.add(db_entry)
			self.session.flush()
		return db_entry
	
	def getStudy(self, short_name=None, description=None):
		"""
		20170412
		"""
		db_entry = self.queryTable(Study).filter_by(short_name=short_name).first()
		if not db_entry:
			db_entry = Study(short_name=short_name, description=description)
			self.session.add(db_entry)
			self.session.flush()
		return db_entry
	
	def getTissue(self, db_entry_id=None, short_name=None, description=None):
		"""
		2011-5-9
			fetch a tissue entry (create one if none existent)
		"""
		db_entry=None
		if db_entry_id:
			db_entry = self.queryTable(Tissue).get(db_entry_id)
		
		if not db_entry and short_name:
			db_entry = self.queryTable(Tissue).filter_by(short_name=short_name).first()
		if not db_entry and (db_entry_id or short_name):
			db_entry = Tissue(id=db_entry_id, short_name=short_name, description=description)
			self.session.add(db_entry)
			self.session.flush()
		return db_entry
	
	def getCountry(self, country_name=None, capital=None, abbr=None, latitude=None, longitude=None):
		"""
		2011-4-28
		
		"""
		db_entry = self.queryTable(Country).filter_by(name=country_name).first()
		if not db_entry:
			db_entry = Country(name=country_name, capital=capital, abbr=abbr, latitude=latitude, longitude=longitude)
			self.session.add(db_entry)
			self.session.flush()
		return db_entry
	
	def getUser(self, name=None, username=None):
		"""
		2011-4-28
		"""
		query = self.queryTable(User).filter_by(realname=name)
		if username:
			query = query.filter_by(username=username)
		db_entry = query.first()
		if not db_entry:
			db_entry = User(realname=name, username=username, _password="123456")
			self.session.add(db_entry)
			self.session.flush()
		return db_entry
	
	def getPhenotype(self, phenotype_name=None, value=None, replicate=None, individual_id=None, comment=None, collector_name=None):
		"""
		2011-4-28
		"""
		phenotype_method = self.getPhenotypeMethod(phenotype_name=phenotype_name, collector_name=collector_name)
		if replicate is None:
			replicate = 1
		db_entry = self.queryTable(Phenotype).filter_by(phenotype_method_id=phenotype_method.id).filter_by(individual_id=individual_id).\
			filter_by(replicate=replicate).first()
		if not db_entry:
			db_entry = Phenotype(phenotype_method=phenotype_method, value=value, replicate=replicate,\
								individual_id=individual_id, comment=comment)
			self.session.add(db_entry)
			self.session.flush()
		return db_entry
		
		
	def getPhenotypeMethod(self, phenotype_name=None, collector_name=None):
		"""
		2011-4-28
		"""
		db_entry = self.queryTable(PhenotypeMethod).filter_by(short_name=phenotype_name).first()
		if not db_entry:
			if collector_name:
				collector = self.getUser(name=collector_name)
			else:
				collector = None
			db_entry = PhenotypeMethod(short_name=phenotype_name, collector=collector)
			self.session.add(db_entry)
			self.session.flush()
		return db_entry
	
	def getIndividualSequenceID2FilePairLs(self, individualSequenceIDList=None, data_dir=None, needPair=True, checkOldPath=False):
		"""
		2013.3.14 replace SR, PE with individual_sequence.sequence_type
		2012.2.10
			add argument checkOldPath.
				True: find files in IndividualSequence.path[:-6], (=path without the trailing '_split').
					This is the old format.
				False: find files in IndividualSequence.path 
		2011-8-30
			filename in individualSequenceID2FilePairLs is path relative to data_dir
		2011-8-28
			add argument needPair
			copied from MpiBaseCount.py
		2011-8-5
		"""
		sys.stderr.write("Getting individualSequenceID2FilePairLs ...")
		individualSequenceID2FilePairLs = {}
		if not data_dir:
			data_dir = self.data_dir
		for individualSequenceID in individualSequenceIDList:
			individual_sequence = self.queryTable(IndividualSequence).get(individualSequenceID)
			if individual_sequence and individual_sequence.path:
				if checkOldPath:
					path = individual_sequence.path[:-6]
				else:
					path = individual_sequence.path
				abs_path = os.path.join(data_dir, path)
				if individual_sequence.id not in individualSequenceID2FilePairLs:
					individualSequenceID2FilePairLs[individual_sequence.id] = []
				if os.path.isfile(abs_path):
					fileRecord = [path, individual_sequence.format, individual_sequence.sequence_type, individual_sequence.sequencer]
						#"SR" means it's single-end
					individualSequenceID2FilePairLs[individual_sequence.id].append([fileRecord])
				elif os.path.isdir(abs_path):	#it's a folder, sometimes it's nothing there
					if individual_sequence.sequence_type.paired_end==1:
						isPE = True
					else:
						isPE = False
					pairedEndPrefix2FileLs = ngs.getPEInputFiles(abs_path, isPE=isPE)
					for pairedEndPrefix, fileLs in pairedEndPrefix2FileLs.items():
						if isPE and len(fileLs)==2 and fileLs[0] and fileLs[1]:	#PE
							filename = os.path.join(path, fileLs[0])	#take one file only
							fileRecord = [filename, individual_sequence.format, individual_sequence.sequence_type, individual_sequence.sequencer]
							#"PE" means it's paired-end
							filename2 = os.path.join(path, fileLs[1])	#take one file only
							fileRecord2 = [filename2, individual_sequence.format, individual_sequence.sequence_type, individual_sequence.sequencer]
							#"PE" means it's paired-end
							individualSequenceID2FilePairLs[individual_sequence.id].append([fileRecord, fileRecord2])	#"PE" means it's paired-end
						else:
							for filename in fileLs:	#usually should be only one file
								if filename:
									filename = os.path.join(path, filename)
									fileRecord = [filename, individual_sequence.format, individual_sequence.sequence_type, individual_sequence.sequencer]
									#"SR" means it's single-end
									individualSequenceID2FilePairLs[individual_sequence.id].append([fileRecord])
		sys.stderr.write("%s individual sequences. Done.\n"%(len(individualSequenceID2FilePairLs)))
		return individualSequenceID2FilePairLs
	
	def getISQ_ID2LibrarySplitOrder2FileLs(self, individualSequenceIDList=None, data_dir=None, filtered=None, \
										ignoreEmptyReadFile=True, is_contaminated=0, outdated_index=0):
		"""
		2013.04.05 added argument is_contaminated, outdated_index
		2012.3.19
			add argument, ignoreEmptyReadFile
		2012.2.24
			filtered=None means "no filtering based on this field.".
			
		2012.2.10
			If for one (isq_id, librarySplitOrder), there is only one mate (single-end).
				The isq_id2LibrarySplitOrder2FileLs only stores one file object (FileLs is of length 1).
			Length of FileLs is commesurate with the number of ends.
		"""
		sys.stderr.write("Getting isq_id2LibrarySplitOrder2FileLs for %s isq entries ..."%(len(individualSequenceIDList)))
		isq_id2LibrarySplitOrder2FileLs = {}
		if not data_dir:
			data_dir = self.data_dir
		counter = 0
		
		for individualSequenceID in individualSequenceIDList:
			individual_sequence = self.queryTable(IndividualSequence).get(individualSequenceID)
			if not individual_sequence:	#not present in db, ignore
				continue
			if is_contaminated is not None:
				if individual_sequence.is_contaminated!=is_contaminated:
					sys.stderr.write(" individual_sequence %s 's is_contaminated=%s (!=%s). ignore.\n"%\
									(individual_sequence.id, individual_sequence.is_contaminated, is_contaminated))
					continue
			if outdated_index is not None:
				if individual_sequence.outdated_index!=outdated_index:
					sys.stderr.write(" individual_sequence %s 's outdated_index=%s (!=%s). ignore.\n"%\
									(individual_sequence.id, individual_sequence.outdated_index, outdated_index))
					continue
			for individual_sequence_file in individual_sequence.individual_sequence_file_ls:
				path = os.path.join(data_dir, individual_sequence_file.path)
				if filtered is not None and individual_sequence_file.filtered!=filtered:	#skip entries that don't matched the filtered argument
					continue
				if ignoreEmptyReadFile:	#2012.3.19	ignore empty read files.
					if individual_sequence_file.read_count is None:	#calculate it on the fly
						baseCountData = ngs.getReadBaseCount(path, onlyForEmptyCheck=True)
						read_count = baseCountData.read_count
					else:
						read_count = individual_sequence_file.read_count
					if read_count==0:
						continue
				if individualSequenceID not in isq_id2LibrarySplitOrder2FileLs:
					isq_id2LibrarySplitOrder2FileLs[individualSequenceID] = {}
				counter += 1
				LibrarySplitOrder2FileLs = isq_id2LibrarySplitOrder2FileLs[individualSequenceID]
				library = individual_sequence_file.library
				split_order = individual_sequence_file.split_order
				mate_id = individual_sequence_file.mate_id
				if mate_id is None:
					mate_id = 1
				key = (library, split_order)
				if key not in LibrarySplitOrder2FileLs:
					LibrarySplitOrder2FileLs[key] = []
				if len(LibrarySplitOrder2FileLs[key])<mate_id:
					for i in range(mate_id-len(LibrarySplitOrder2FileLs[key])):	#expand the list to match the number of mates
						LibrarySplitOrder2FileLs[key].append(None)
				isq_file_obj = PassingData(db_entry=individual_sequence_file, path=path)
				LibrarySplitOrder2FileLs[key][mate_id-1] = isq_file_obj
		
		sys.stderr.write("%s individual sequence files from %s isq entries.\n"%(counter, len(isq_id2LibrarySplitOrder2FileLs)))
		return isq_id2LibrarySplitOrder2FileLs
	
	def filterIndividualSequenceList(self, individual_sequence_list=None, min_coverage=None, max_coverage=None, individual_site_id=None, \
						individual_site_id_set=None, individual_id_set=None, sequence_type_id_set=None,\
						sequencer_id_set=None, sequence_filtered=None,\
						sequence_batch_id_set=None, parent_individual_sequence_id_set=None, \
						version_set=None,\
						country_id_set=None, tax_id_set=None, excludeContaminant=False, excludeTissueIDSet=set([6]),\
						outdated_index=0,\
						report=True):
		"""
		2013.3.15 added argument outdated_index
		2013.3.15 use individual_sequence.is_contaminated, instead of individual_sequence.individual.is_contaminated
		2013.3.15
		"""
		if report:
			sys.stderr.write("Filter %s individual_sequence entries to select %s=<coverage <=%s & site-id=%s & sequence_filtered=%s & from %s sites, \n\
	%s countries, %s taxonomies, %s individuals, %s sequence types, %s sequencers, %s sequence batches, %s parent isqs, \n\
	%s versions; excludeContaminant=%s, outdated_index=%s, %s tissues to be excluded,..."%\
							(len(individual_sequence_list), min_coverage, max_coverage, individual_site_id, sequence_filtered, \
							getattr(individual_site_id_set, '__len__', returnZeroFunc)(),\
							getattr(country_id_set, '__len__', returnZeroFunc)(),\
							getattr(tax_id_set, '__len__', returnZeroFunc)(),\
							getattr(individual_id_set, '__len__', returnZeroFunc)(),\
							getattr(sequence_type_id_set, '__len__', returnZeroFunc)(),\
							getattr(sequencer_id_set, '__len__', returnZeroFunc)(),\
							getattr(sequence_batch_id_set, '__len__', returnZeroFunc)(),\
							getattr(parent_individual_sequence_id_set, '__len__', returnZeroFunc)(),\
							getattr(version_set, '__len__', returnZeroFunc)(),\
							excludeContaminant, outdated_index,\
							getattr(excludeTissueIDSet, '__len__', returnZeroFunc)(),\
							)
						)
		listToReturn = []
		cumulative_coverage = 0
		for individual_sequence in individual_sequence_list:
			if not individual_sequence:
				continue
			if min_coverage is not None and individual_sequence.coverage<min_coverage:
				continue
			if max_coverage is not None and individual_sequence.coverage>max_coverage:
				continue
			if individual_site_id is not None and individual_sequence.individual.site_id!=individual_site_id:
				continue
			if sequence_filtered is not None and individual_sequence.filtered!=sequence_filtered:
				continue
			if individual_site_id_set and individual_sequence.individual.site_id not in individual_site_id_set:
				#2012.4.13
				continue
			if individual_id_set and individual_sequence.individual_id not in individual_id_set:
				continue
			if sequence_type_id_set and individual_sequence.sequence_type_id not in sequence_type_id_set:
				continue
			if sequencer_id_set and individual_sequence.sequencer_id not in sequencer_id_set:
				continue
			if sequence_batch_id_set and individual_sequence.sequence_batch_id not in sequence_batch_id_set:
				continue
			if parent_individual_sequence_id_set and \
				individual_sequence.parent_individual_sequence_id not in parent_individual_sequence_id_set:
				continue
			if version_set and individual_sequence.version not in version_set:
				continue
			if country_id_set:
				if individual_sequence.individual.site is None:
					sys.stderr.write("Warning: isq (id=%s, path=%s, %s) has no site.\n"%(individual_sequence.id, individual_sequence.path,\
																			individual_sequence.individual.code))
					continue
				elif (individual_sequence.individual.site is None or \
									individual_sequence.individual.site.country_id not in country_id_set):
					continue
			if tax_id_set:
				if individual_sequence.individual.tax_id is None:
					sys.stderr.write("Warning: alignment (id=%s, path=%s, %s) has no tax_id.\n"%(individual_sequence.id, individual_sequence.path,\
																	individual_sequence.individual.code))
					continue
				elif individual_sequence.individual.tax_id not in tax_id_set:
					continue
			if excludeContaminant and individual_sequence.is_contaminated:	#2012.9.27
				continue
			if outdated_index is not None and individual_sequence.outdated_index!=outdated_index:	#2013.3.15
				continue
			if excludeTissueIDSet and individual_sequence.tissue_id in excludeTissueIDSet:	#2012.10.2
				continue
			listToReturn.append(individual_sequence)
			if individual_sequence.coverage is not None:
				cumulative_coverage += individual_sequence.coverage
		if report:
			sys.stderr.write(" kept %s individual_sequence entries, cumulative_coverage=%s.\n"%\
							(len(listToReturn), cumulative_coverage))
		return listToReturn
	
	def getISQDBEntryLsForAlignment(self, individualSequenceIDList=None, data_dir=None, filtered=None, ignoreEmptyReadFile=True,\
								outdated_index=0):
		"""
		2013.3.15 added argument outdated_index
		2012.9.19 similar to getISQ_ID2LibrarySplitOrder2FileLs()
			isqLs
				library2Data
					isqFileRawID2Index
					isqFileRawDBEntryLs
					splitOrder2Index
					fileObjectPairLs
		"""
		sys.stderr.write("Getting isqLs given %s isq-IDs ..."%(len(individualSequenceIDList)))
		isqLs = []
		if not data_dir:
			data_dir = self.data_dir
		counter = 0
		cumulative_coverage = 0
		for individualSequenceID in individualSequenceIDList:
			individual_sequence = self.queryTable(IndividualSequence).get(individualSequenceID)
			if not individual_sequence:	#not present in db, ignore
				continue
			if outdated_index is not None and individual_sequence.outdated_index!=outdated_index:	#2013.3.15
				continue
			library2Data = {}	#2012.9.19
			for individual_sequence_file in individual_sequence.individual_sequence_file_ls:
				path = os.path.join(data_dir, individual_sequence_file.path)
				if filtered is not None and individual_sequence_file.filtered!=filtered:	#skip entries that don't matched the filtered argument
					continue
				if ignoreEmptyReadFile:	#2012.3.19	ignore empty read files.
					if individual_sequence_file.read_count is None:	#calculate it on the fly
						baseCountData = ngs.getReadBaseCount(path, onlyForEmptyCheck=True)
						read_count = baseCountData.read_count
					else:
						read_count = individual_sequence_file.read_count
					if read_count==0:
						continue
				isqFileRawID = individual_sequence_file.individual_sequence_file_raw_id
				
				
				counter += 1
				library = individual_sequence_file.library
				splitOrder = individual_sequence_file.split_order
				mate_id = individual_sequence_file.mate_id
				if library not in library2Data:
					library2Data[library] = PassingData(isqFileRawID2Index = {}, isqFileRawDBEntryLs=[], \
											splitOrder2Index={}, fileObjectPairLs=[])
				#for convenience
				isqFileRawID2Index = library2Data[library].isqFileRawID2Index
				isqFileRawDBEntryLs = library2Data[library].isqFileRawDBEntryLs
				splitOrder2Index = library2Data[library].splitOrder2Index
				fileObjectPairLs = library2Data[library].fileObjectPairLs
				if isqFileRawID not in isqFileRawID2Index:
					isqFileRawID2Index[isqFileRawID] = len(isqFileRawID2Index)
					isqFileRawDBEntryLs.append(individual_sequence_file.individual_sequence_file_raw)
				
				if splitOrder not in splitOrder2Index:
					splitOrder2Index[splitOrder] = len(splitOrder2Index)
					fileObjectPairLs.append([])
				fileObjectPairIndex = splitOrder2Index[splitOrder]
				if mate_id is None:
					mate_id = 1
				noOfFileObjects = len(fileObjectPairLs[fileObjectPairIndex])
				if noOfFileObjects<mate_id:
					for i in range(mate_id-noOfFileObjects):	#expand the list to match the number of mates
						fileObjectPairLs[fileObjectPairIndex].append(None)
				isq_file_obj = PassingData(db_entry=individual_sequence_file, path=path)
				fileObjectPairLs[fileObjectPairIndex][mate_id-1] = isq_file_obj
			#add it to the individual_sequence db entry
			individual_sequence.library2Data = library2Data
			isqLs.append(individual_sequence)
			if individual_sequence.coverage is not None:
				cumulative_coverage += individual_sequence.coverage
		sys.stderr.write("%s individual sequence files from %s isq entries. cumulative_coverage=%s.\n"%\
						(counter, len(isqLs), cumulative_coverage))
		return isqLs
	
	
	def constructPedigreeGraphOutOfAlignments(self, alignmentLs=None, useAlignmentIDAsNodeID=False):
		"""
		2013.06.13 added argument useAlignmentIDAsNodeID
			default is False, which means individual.id is used as node id.
			If True, then alignment.id is used as node id.
				If multiple alignments for one individual, only id of the first alignment will be used.
		2012.4.20
			Warning if one individual appears in the graph >1 times. 
		2011-12-15
			copied from AlignmentToTrioCallPipeline.py
			
		construct a directed graph (edge: from parent to child) of which nodes are all from alignmentLs.
		"""
		sys.stderr.write("Construct pedigree out of %s alignments, useAlignmentIDAsNodeID=%s... "%\
						(len(alignmentLs), useAlignmentIDAsNodeID))
		from palos.algorithm import graph
		DG=graph.DiGraphWrapper()
		
		individual_id2alignmentLs = {}
		for alignment in alignmentLs:
			individual_id = alignment.individual_sequence.individual_id
			if individual_id not in individual_id2alignmentLs:
				individual_id2alignmentLs[individual_id] = [alignment]
			else:
				individual_id2alignmentLs[individual_id].append(alignment)
				sys.stderr.write("Warning: individual_id %s appears in >%s alignments.\n"%(individual_id, \
																	len(individual_id2alignmentLs[individual_id])))
				continue
			if useAlignmentIDAsNodeID:
				node_id = alignment.id
			else:
				node_id = individual_id
			if DG.has_node(node_id):
				sys.stderr.write("Warning: individual_id %s already in the graph.\n"%(individual_id))
			else:
				DG.add_node(node_id)
		
		for row in self.queryTable(Ind2Ind):
			if row.individual1_id in individual_id2alignmentLs and row.individual2_id in individual_id2alignmentLs:
				if useAlignmentIDAsNodeID:
					node1_id = individual_id2alignmentLs[row.individual1_id][0].id
					node2_id = individual_id2alignmentLs[row.individual2_id][0].id
				else:
					node1_id = row.individual1_id
					node2_id = row.individual2_id
				DG.add_edge(node1_id, node2_id)
		import networkx
		sys.stderr.write("%s nodes (%s alignments). %s edges. %s connected components.\n"%(DG.number_of_nodes(), \
							len(alignmentLs), DG.number_of_edges(), \
							networkx.number_connected_components(DG.to_undirected())))
		return PassingData(DG=DG, individual_id2alignmentLs=individual_id2alignmentLs)
	
	def constructPedigree(self, directionType=1, individualIDSet=None):
		"""
		2013.08.09 added argument individualIDSet
		2012.9.6 add argument directionType
			1: parent -> child as edge
			2: child -> parent as edge
			3: undirected
		2012.1.23
		"""
		sys.stderr.write("Constructing pedigree from db directionType=%s..."%(directionType))
		from palos.algorithm import graph
		if directionType==3:
			DG = graph.GraphWrapper()
		else:
			DG = graph.DiGraphWrapper()
		
		for row in self.queryTable(Ind2Ind):
			if individualIDSet and (row.individual1_id not in individualIDSet or row.individual2_id not in individualIDSet):
				#2013.08.09 ignore nodes that are not in individualIDSet if the latter is not None
				continue
			if directionType==2:
				DG.add_edge(row.individual2_id, row.individual1_id)
			else:
				DG.add_edge(row.individual1_id, row.individual2_id)
		
		"""
		#2012.6.19 initialization, but found unnecessary. increasePredecessorNoOfDescendantCount() can do it itself.
		for node in DG:
			DG.node[node]["noOfDescendant"] = 0	#initialize this
			DG.node[node]["noOfChildren"] = len(DG.successors(node))	#
		"""
		import networkx
		sys.stderr.write("%s nodes. %s edges. %s connected components.\n"%(DG.number_of_nodes(), DG.number_of_edges(), \
					networkx.number_connected_components(DG.to_undirected())))
		return DG
	
	def accumulatePredecessorAttributeRecursively(self, DG=None, source=None, attributeName='noOfDescendant'):
		"""
		2012.11.30
			bugfix. use copy.deepcopy around attributeInitValue or attributeIncrementValue everytime.
		2012.6.19
			attributeName=
				noOfDescendant: wrong, don't use it
				descendantList: treated as a set
				descendantSet: same as descendantList
			a recursive function to accumulate the attribute of the predecessors of the source node
				by absorbing the attribute of the source node.
			The algorithm starts from the bottom of the pedigree. called in calculateCumulativeAttributeForEachNodeInPedigree() 
		"""
		import copy
		predecessorList = DG.predecessors(source)
		if attributeName =='noOfDescendant':
			attributeInitValue = 0
			attributeIncrementValue = 1
		elif attributeName =='descendantList' or attributeName=='descendantSet':
			attributeInitValue = set()
			attributeIncrementValue = set([source])
		else:	#same as noOfDescendant
			attributeInitValue = 0
			attributeIncrementValue = 1
		if attributeName not in DG.node[source]:	#never accessed
			DG.node[source][attributeName] = copy.deepcopy(attributeInitValue)	#avoid "by-reference"
		for node in predecessorList:
			if attributeName not in DG.node[node]:
				DG.node[node][attributeName] = copy.deepcopy(attributeInitValue)
			if attributeName=='descendantList' or attributeName=='descendantSet':
				DG.node[node][attributeName] |= copy.deepcopy(DG.node[source][attributeName])
				DG.node[node][attributeName] |= copy.deepcopy(attributeIncrementValue)	#source itself
			else:
				DG.node[node][attributeName] += copy.deepcopy(DG.node[source][attributeName])
				DG.node[node][attributeName] += copy.deepcopy(attributeIncrementValue)	#source itself
			self.accumulatePredecessorAttributeRecursively(DG, source=node, attributeName=attributeName)
	
	def calculateCumulativeAttributeForEachNodeInPedigree(self, DG=None, attributeName='noOfDescendant'):
		"""
		2012.6.19
			the algorithm starts from the bottom (leaves) of the graph and do a reverse BFS
		"""
		leafNodeList = []
		for node in DG:
			if len(DG.successors(node))==0:
				leafNodeList.append(node)
		sys.stderr.write("%s leaf nodes.\n"%(len(leafNodeList)))
		for leafNode in leafNodeList:
			self.accumulatePredecessorAttributeRecursively(DG, source=leafNode, attributeName=attributeName)
		
	
	def getNoOfDescendant(self, individual_id=None, attributeName = 'noOfDescendant'):
		"""
		2012.6.19
			this function is wrong. don't use it. use getDescendantSet() to get a unique set of descendants and get its length.
			It's wrong because in BFS upwards traversal, if one parent has >1 children, its descendant count would be added to
				that of grandparent for >1 times, which is bad.
		"""
		if not getattr(self, 'pedigreeDG', None):
			self.pedigreeDG = self.constructPedigree()
		
		if individual_id not in self.pedigreeDG.node:
			return None
		
		firstNode = self.pedigreeDG.nodes()[0]
		nodeProperty = self.pedigreeDG.node[firstNode]
		if attributeName not in nodeProperty:
			self.calculateCumulativeAttributeForEachNodeInPedigree(self.pedigreeDG, attributeName=attributeName)
		return self.pedigreeDG.node[individual_id][attributeName]
	
	def getDescendantSet(self, individual_id=None, attributeName = 'descendantList'):
		"""
		2012.6.19
			descendantList is treated as descendantSet.
			similar to getNoOfDescendant() (which is wrong) but return a set of descendant nodes
		"""
		if not getattr(self, 'pedigreeDG', None):
			self.pedigreeDG = self.constructPedigree()
		
		if individual_id not in self.pedigreeDG.node:
			return None
		
		firstNode = self.pedigreeDG.nodes()[0]
		nodeProperty = self.pedigreeDG.node[firstNode]
		if attributeName not in nodeProperty:	#not run yet
			self.calculateCumulativeAttributeForEachNodeInPedigree(self.pedigreeDG, attributeName=attributeName)
		return self.pedigreeDG.node[individual_id][attributeName]
	
	def getAncestorSet(self, individual_id=None, attributeName = 'descendantList'):
		"""
		2012.11.30
			use a reverse pedigree DG
		2012.6.19
			descendantList is treated as descendantSet.
			similar to getNoOfDescendant() (which is wrong) but return a set of descendant nodes
		"""
		if not getattr(self, 'pedigreeReverseDG', None):
			self.pedigreeReverseDG = self.constructPedigree(directionType=2)
		
		if individual_id not in self.pedigreeReverseDG.node:
			return None
		
		firstNode = self.pedigreeReverseDG.nodes()[0]
		nodeProperty = self.pedigreeReverseDG.node[firstNode]
		if attributeName not in nodeProperty:	#not run yet
			self.calculateCumulativeAttributeForEachNodeInPedigree(self.pedigreeReverseDG, attributeName=attributeName)
		return self.pedigreeReverseDG.node[individual_id][attributeName]
	
	def getPedigreeSplitStructure(self, pedigreeGraph=None, removeFamilyFromGraph=True):
		"""
		2013.06.13
			This function splits a pedigree into list of trios, duos, singletons.
			argument removeFamilyFromGraph decides whether
				1. there will be overlap (removeFamilyFromGraph=False) among the trios/duos/singletons
				2. or not (removeFamilyFromGraph=True).
				"Overlap" means when nuclear families or half-siblings exist in the pedigree,
					whether parents should be replicated in order for them to be associated with each kid during split.
			one example in BeagleAndTrioCallerOnVCFWorkflow.py
		"""
		# work on the pedigree graph to figure out if singleton, trio, duo file will exist.
		#figure out the singletons, duos, trios, using the function in AlignmentToTrioCall
		pedigreeSplitStructure = PassingData(familySize2familyLs={})
				
		#find trios first
		pedigreeSplitStructure.familySize2familyLs[3] = self.findFamilyFromPedigreeGivenSize(DG=pedigreeGraph, familySize=3, \
																							removeFamilyFromGraph=removeFamilyFromGraph)
		#find duos
		pedigreeSplitStructure.familySize2familyLs[2] = self.findFamilyFromPedigreeGivenSize(DG=pedigreeGraph, familySize=2, \
																							removeFamilyFromGraph=removeFamilyFromGraph)
		#find singletons (familySize=1 => noOfIncomingEdges=0, noOfIncomingEdges=0 => will not be parents of others)
		pedigreeSplitStructure.familySize2familyLs[1] = self.findFamilyFromPedigreeGivenSize(DG=pedigreeGraph, familySize=1, \
																							removeFamilyFromGraph=removeFamilyFromGraph, \
																							noOfOutgoingEdges=0)
		return pedigreeSplitStructure
	
	def findFamilyFromPedigreeGivenSize(self, DG=None, familySize=3, removeFamilyFromGraph=True, noOfOutgoingEdges=None):
		"""
		2012.4.2
			add noOfOutgoingEdges to find true singletons.
		2011.12.15
			copied from AlignmentToTrioCallPipeline.py
		2011-12.4
			DG is a directed graph.
			It views the graph from the offspring point of view because it's designed to find max number of trios (for trioCaller).
			Nuclear families (size>3) could not be handled by trioCaller and not considered in this algorithm.
			
			1. sort all nodes by their out degree ascendingly
				nodes with higher out degree would be checked later in the process to minimize the disruption of dependent trios.
				(if this offspring happens to be parent in other trios).
			2. if a node's incoming degree is same as familySize-1, then a family of given size is found.
		"""
		sys.stderr.write("Finding families of size %s .."%(familySize))
		familyLs = []	#each element of this list is either [singleton] or [parent, child] or [father, mother, child]
		allNodes = DG.nodes()
		out_degree_node_ls = []
		for node in allNodes:
			out_degree_node_ls.append([DG.out_degree(node), node])
		out_degree_node_ls.sort()	#nodes with low out degree would be considered first.
		
		for out_degree, node in out_degree_node_ls:
			if DG.has_node(node):	#it could have been removed as trios/duos were found
				#edges = DG.in_edges(node)
				noOfIncomingEdges = DG.in_degree(node)
				if noOfOutgoingEdges is not None and out_degree!=noOfOutgoingEdges:
					continue
				#noOfIncomingEdges = len(edges)
				if noOfIncomingEdges==familySize-1:	# a trio is found
					parents = DG.predecessors(node)
					family = parents + [node]
					familyLs.append(family)
					if removeFamilyFromGraph:
						DG.remove_nodes_from(family)
				elif noOfIncomingEdges>familySize-1:
					sys.stderr.write("Warning: noOfIncomingEdges for node %s is %s (family size =%s).\n"%(node, \
																noOfIncomingEdges, noOfIncomingEdges+1))
		sys.stderr.write(" found %s.\n"%(len(familyLs)))
		return familyLs
	
	@classmethod
	def parseAlignmentCompositeID(cls, compositeID):
		"""
		2012.5.4
			composite ID is used to designate members in trio, similar to read groups used in bam files but different.
			It's generated via IndividualAlignment.getCompositeID()
		
				compositeID = '%s_%s_%s_%s'%(self.id, self.ind_seq_id, self.individual_sequence.individual.id, \
								self.individual_sequence.individual.code)
		"""
		split_id_ls = compositeID.split('_')
		individual_alignment_id = int(split_id_ls[0])
		ind_seq_id = int(split_id_ls[1])
		individual_id = int(split_id_ls[2])
		individual_code = split_id_ls[3]
		return PassingData(individual_alignment_id=individual_alignment_id, individual_sequence_id=ind_seq_id,\
						individual_id=individual_id, individual_code=individual_code)
	
	@classmethod
	def parseAlignmentReadGroupWithoutDB(cls, read_group):
		"""
		2012.5.11
			read-group is used to designate samples in alignment files and also in the ensuing multi-sample VCF files,
				similar to compositeID used in trio-representations but different.
			It's generated via IndividualAlignment.getReadGroup()
				read_group = '%s_%s_%s_%s_vs_%s'%(self.id, self.ind_seq_id, self.individual_sequence.individual.code, \
								sequencer, self.ref_ind_seq_id)
								
		"""
		individual_alignment_id = None
		ind_seq_id = None
		individual_code = None
		sequencer = None
		ref_ind_seq_id = None
		individual_id = None
		try:
			split_id_ls = read_group.split('_')
			individual_alignment_id = int(split_id_ls[0])
			ind_seq_id = int(split_id_ls[1])
			individual_code = split_id_ls[2]
			sequencer = split_id_ls[3]
			ref_ind_seq_id = int(split_id_ls[5])	#index 4 is "vs"
			individual_id = None
		except:
			sys.stderr.write('Except type: %s\n'%repr(sys.exc_info()))
			import traceback
			traceback.print_exc()
			
			
		return PassingData(individual_alignment_id=individual_alignment_id, individual_sequence_id=ind_seq_id,\
						individual_id=individual_id, individual_code=individual_code, ref_ind_seq_id=ref_ind_seq_id,\
						sequencer=sequencer, individualAlignment=None)
	
	def parseAlignmentReadGroup(self, read_group):
		"""
		2013.3.18 bugfix: sequencer is now a db entry.
		2012.10.8 no more try except
		2012.5.10
			read-group is used to designate samples in alignment files and also in the ensuing multi-sample VCF files,
				similar to compositeID used in trio-representations but different.
			It's generated via IndividualAlignment.getReadGroup()
				read_group = '%s_%s_%s_%s_vs_%s'%(self.id, self.ind_seq_id, self.individual_sequence.individual.code, \
								sequencer, self.ref_ind_seq_id)
								
		"""
		split_id_ls = read_group.split('_')
		try:
			individual_alignment_id = int(split_id_ls[0])
		except:	#only catch the string error, not catching the db error,
			sys.stderr.write('Except type: %s\n'%repr(sys.exc_info()))
			import traceback
			traceback.print_exc()
			individualAlignment = None
			individual_alignment_id = None
			ind_seq_id = None
			individual_code = None
			sequencer = None
			ref_ind_seq_id = None
			individual_id = None
			return PassingData(individual_alignment_id=individual_alignment_id, individual_sequence_id=ind_seq_id,\
						individual_id=individual_id, individual_code=individual_code, ref_ind_seq_id=ref_ind_seq_id,\
						sequencer=sequencer, individualAlignment=individualAlignment)
		
		individualAlignment = self.queryTable(IndividualAlignment).get(individual_alignment_id)
		ind_seq_id = individualAlignment.individual_sequence.id
		individual_code = individualAlignment.individual_sequence.individual.code
		sequencer = individualAlignment.individual_sequence.sequencer
		ref_ind_seq_id = individualAlignment.ref_sequence.id
		individual_id = individualAlignment.individual_sequence.individual.id
		"""
		except:
			sys.stderr.write('Except type: %s\n'%repr(sys.exc_info()))
			import traceback
			traceback.print_exc()
			individualAlignment = None
			individual_alignment_id = None
			ind_seq_id = None
			individual_code = None
			sequencer = None
			ref_ind_seq_id = None
			individual_id = None
		"""
		return PassingData(individual_alignment_id=individual_alignment_id, individual_sequence_id=ind_seq_id,\
						individual_id=individual_id, individual_code=individual_code, ref_ind_seq_id=ref_ind_seq_id,\
						sequencer=sequencer.short_name, individualAlignment=individualAlignment)
	
	def parseDBEntryGivenDBAffiliatedFilename(self, filename=None, TableClass=None):
		"""
		2012.8.15
			almost all db-affliated filenames are beginned with the db_entry's primary ID, separated by _.
			this is a generic function, to be called by others with different TableClass
		"""
		baseFilename = os.path.basename(filename)
		try:
			dbID = int(baseFilename.split('_')[0])
			db_entry = self.queryTable(TableClass).get(dbID)
		except:
			sys.stderr.write("Warning: could not parse dbID from %s.\n"%(baseFilename))
			db_entry = None
		return db_entry
	
	def parseGenotypeFileGivenDBAffiliatedFilename(self, filename=None):
		"""
		2012.8.15
			call parseDBEntryGivenDBAffiliatedFilename
		"""
		return self.parseDBEntryGivenDBAffiliatedFilename(filename=filename, TableClass=GenotypeFile)
	
	def parseIndividualAlignmentGivenDBAffiliatedFilename(self, filename=None):
		"""
		2012.8.15
			call parseDBEntryGivenDBAffiliatedFilename
		"""
		return self.parseDBEntryGivenDBAffiliatedFilename(filename=filename, TableClass=IndividualAlignment)
	
	def findISQFoldersNotInDB(self, data_dir=None, subFolder='individual_sequence',):
		"""
		2012.7.11 scan the individual_sequence folder to see which subfolder is not in individual_sequence.path
			i.e.
				orphanPathList = db_vervet.findISQFoldersNotInDB(data_dir=None, subFolder='individual_sequence',)
				orphanPathList = db_vervet.findISQFoldersNotInDB(data_dir='/u/home/eeskin2/polyacti/NetworkData/vervet/db/', \
														subFolder='individual_sequence',)
		"""
		sys.stderr.write("Finding individual_sequence folders that are not in db (orphaned)... ")
		dbPath2dbEntry = {}
		for row in self.queryTable(IndividualSequence):
			if row.path and row.path not in dbPath2dbEntry:
				dbPath2dbEntry[row.path] = row
		sys.stderr.write("%s items in db.\n"%(len(dbPath2dbEntry)))
		
		if not data_dir:
			data_dir = self.data_dir
		
		topFolder = os.path.join(data_dir, subFolder)
		orphanPathList = []
		for item in os.listdir(topFolder):
			itemAbsPath = os.path.join(topFolder, item)
			if os.path.isdir(itemAbsPath):
				isqPath = os.path.join(subFolder, item)
				if isqPath not in dbPath2dbEntry:
					orphanPathList.append(itemAbsPath)
		sys.stderr.write(" %s found.\n"%(len(orphanPathList)))
		return orphanPathList
	
	def findAlignmentFileNotInDB(self, data_dir=None, subFolder='individual_alignment',):
		"""
		2012.7.11 scan the individual_alignment folder to see which file is not in individual_alignment.path
			i.e.
				orphanPathList = db_vervet.findAlignmentFileNotInDB(data_dir=None, subFolder='individual_alignment',)
		"""
		sys.stderr.write("Finding individual_alignment folders that are not in db (orphaned)... ")
		dbPath2dbEntry = {}
		for row in self.queryTable(IndividualAlignment):
			if row.path and row.path not in dbPath2dbEntry:
				dbPath2dbEntry[row.path] = row
		
		sys.stderr.write("%s items in db.\n"%(len(dbPath2dbEntry)))
		
		if not data_dir:
			data_dir = self.data_dir
		
		topFolder = os.path.join(data_dir, subFolder)
		orphanPathList = []
		for item in os.listdir(topFolder):
			itemAbsPath = os.path.join(topFolder, item)
			
			if os.path.isfile(itemAbsPath) and utils.getRealPrefixSuffixOfFilenameWithVariableSuffix(itemAbsPath)[1]=='.bam':
				#make sure the suffix is bam. loads of .bai files are not recorded in db.
				isqPath = os.path.join(subFolder, item)
				if isqPath not in dbPath2dbEntry:
					orphanPathList.append(itemAbsPath)
		sys.stderr.write(" %s found.\n"%(len(orphanPathList)))
		return orphanPathList
		
	def updateGenotypeMethodNoOfLoci(self, db_entry=None, format=None, no_of_chromosomes=None):
		"""
		2012.8.30
		2012.8.15
			add argument format and aggregate the number of loci according to format
		2012.7.17
		"""
		no_of_loci = 0
		formatAndNoOfChromosomes2NoOfLoci = {}
		query = self.queryTable(GenotypeFile).filter_by(genotype_method_id=db_entry.id)
		if format:
			query = query.filter_by(format=format)
		if no_of_chromosomes:	#2012.8.30
			query = query.filter_by(no_of_chromosomes = no_of_chromosomes)
		
		for genotype_file in query:
			key = (genotype_file.format, genotype_file.no_of_chromosomes)
			if key not in formatAndNoOfChromosomes2NoOfLoci:
				formatAndNoOfChromosomes2NoOfLoci[key] = 0
			formatAndNoOfChromosomes2NoOfLoci[key] += genotype_file.no_of_loci
		noOfLociList = formatAndNoOfChromosomes2NoOfLoci.values()
		#do some checking here to make sure every number in noOfLociList agrees with each other
		db_entry.no_of_loci = max(noOfLociList)	#take the maximum
		self.session.add(db_entry)
		self.session.flush()
		
	def updateOneEntryByReducingItsAffiliatedTables(self, db_entry=None, affiliate_table_query=None, \
									affiliate_table_field_name='no_of_loci', \
									db_entry_field_name='no_of_loci', ):
		"""
		2013.08.23
		"""
		reduce_value = 0
		for affiliate_table_entry in affiliate_table_query:
			reduce_value += getattr(affiliate_table_entry, affiliate_table_field_name, 0)
		setattr(db_entry, db_entry_field_name, reduce_value)
		self.session.add(db_entry)
		self.session.flush()

	def updateAlignmentDepthIntervalMethodNoOfIntervals(self, db_entry=None, format=None, no_of_chromosomes=None):
		"""
		2013.08.23
		"""
		query = self.queryTable(AlignmentDepthIntervalFile).filter_by(alignment_depth_interval_method_id=db_entry.id)
		if format:
			query = query.filter_by(format=format)
		if no_of_chromosomes:	#2012.8.30
			query = query.filter_by(no_of_chromosomes = no_of_chromosomes)
		self.updateOneEntryByReducingItsAffiliatedTables(db_entry=db_entry, affiliate_table_query=query,\
							affiliate_table_field_name='no_of_intervals', db_entry_field_name='no_of_intervals')
		
	
if __name__ == '__main__':
	main_class = SunsetDB
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	schema = po.long_option2value["schema"]
	if schema != _schemaname_:
		sys.stderr.write("ERROR: Schema given by commandline argument, %s, is different from the static one, %s.\n"%
			(schema, _schemaname_))
		sys.exit(2)
	instance = main_class(**po.long_option2value)
	instance.setup(create_tables=True, Base=Base)
	sys.stderr.write("data_dir is %s.\n"%(instance.data_dir))

	import pdb
	pdb.set_trace()
