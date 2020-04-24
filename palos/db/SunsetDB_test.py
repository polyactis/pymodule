#!/usr/bin/env python3
"""
Examples:
    #setup database in postgresql
    %s -v postgresql -u yh -z pdc -d pmdb -k sunset
    
    #setup database in mysql
    %s -u yh -z papaya.usc.edu
    
Description:
    2017.03.17
    This is the ORM architecture of SunsetDB.
"""
import sys, os, math, copy
__doc__ = __doc__%(sys.argv[0], sys.argv[0])

from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import Column, Integer, String, DateTime, func, Float, Text
from sqlalchemy import create_engine
from sqlalchemy import ForeignKey
from sqlalchemy.orm import relationship
from sqlalchemy import UniqueConstraint, create_engine
from sqlalchemy import and_, or_, not_
from sqlalchemy.engine.url import URL
from sqlalchemy.orm import scoped_session, sessionmaker
from sqlalchemy.schema import ThreadLocalMetaData, MetaData
from sqlalchemy.types import LargeBinary

from datetime import datetime

from palos import ProcessOptions, PassingData
from palos.db import Database, TableClass


Base = declarative_base()

_schemaname_ = "sunset"

class Country(Base, TableClass):
    """
    2017-4-17
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

    site_list = relationship("Site",back_populates="country",cascade='all,delete')
    
class Site(Base, TableClass):
    """
    20170411 add study_id
    2012.12.6 add unique constraint on (latitude, longitude)
    2011-4-29
        add column altitude
    2011-3-1
    2017-4-17
        use new sqlalchemy syntax
    """
    __tablename__ = 'site'
    __table_args__ = {'schema':_schemaname_}

    id = Column(Integer,primary_key=True)
    short_name = Column(String(256))
    description = Column(Text)
    latitude = Column(Float)
    longitude = Column(Float)
    altitude = Column(Float)    #2011-4-29
    study_id = Column(Integer,ForeignKey(_schemaname_+".study.id"))  #2017-04-17
    #study = ManyToOne('%s.Study'%(__name__), colname='study_id', ondelete='CASCADE', onupdate='CASCADE')    #2017.04.11
    city = Column(String(100))
    stateprovince = Column(String(100))
    region = Column(String(100))
    zippostal = Column(String(20))
    country_id = Column(Integer,ForeignKey(_schemaname_+".country.id"))  #2017-04-17
    #country = ManyToOne("%s.Country"%(__name__), colname='country_id', ondelete='CASCADE', onupdate='CASCADE')
    created_by = Column(String(128))
    updated_by = Column(String(128))
    date_created = Column(DateTime, default=datetime.now())
    date_updated = Column(DateTime)
    UniqueConstraint('short_name', 'latitude', 'longitude', 'city', 'stateprovince', 'country_id',name='site_sllcsc')
    UniqueConstraint('short_name', 'latitude', 'longitude',name='site_short_name_latitude_longtude')
    
    study = relationship("Study",back_populates="site_list")
    country = relationship("Country",back_populates="site_list")

class Study(Base,TableClass):
    """
    2013.3.12 table used to group individuals
        table Individual and Ind2Ind refers to this table.
    2017.4.17
        use new sqlalchemy syntax
    """
    __tablename__ = 'study'
    __table_args__ = {'schema':_schemaname_}
    
    id = Column(Integer,primary_key=True)
    short_name = Column(String(256), unique=True)
    description = Column(Text)
    created_by = Column(String(128))
    updated_by = Column(String(128))
    date_created = Column(DateTime, default=datetime.now())
    date_updated = Column(DateTime)
    
    site_list = relationship("Site",back_populates="study",cascade='all,delete')
    

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

class SunsetDB(Database):
    __doc__ = __doc__
    option_default_dict = copy.deepcopy(Database.option_default_dict)
    def __init__(self, **keywords):
        Database.__init__(self, **keywords)
        ##self.setup_engine()    #2012.12.18 it needs __metadata__, __session__ from each db-definition file to be ready. can't be run here.
        self.READMEClass = README    #2012.12.18 required to figure out data_dir
        self.dbID2indEntry = {}

    def getCountry(self, country_name=None, capital=None, abbr=None, latitude=None, longitude=None):
        """
        2011-4-28
        
        """
        db_entry = self.session.query(Country).filter_by(name=country_name).first()
        if not db_entry:
            db_entry = Country(name=country_name, capital=capital, abbr=abbr, latitude=latitude, longitude=longitude)
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
        query = self.session.query(Site)
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
        db_entry = self.session.query(Study).filter_by(short_name=short_name).first()
        if not db_entry:
            db_entry = Study(short_name=short_name, description=description)
            self.session.add(db_entry)
            self.session.flush()
        return db_entry
    
if __name__ == '__main__':
    main_class = SunsetDB
    po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
    instance = main_class(**po.long_option2value)
    #instance.setup(create_tables=True)
    sys.stderr.write("data_dir is %s.\n"%(instance.data_dir))
    import pdb
    pdb.set_trace()