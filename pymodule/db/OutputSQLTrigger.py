#!/usr/bin/env python3
"""
Examples:
	#for postgres tables
	OutputSQLTrigger.py  -i readme
	
	for mysql tables
	OutputSQLTrigger.py  -i gene,readme -y 2

Description:
	Output triggers for standard fields, created_by, updated_by, date_created and date_updated
	
"""

import sys, os, math
sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

mysql_trigger = """
DELIMITER |
CREATE TRIGGER before_insert_%s BEFORE INSERT ON %s
  FOR EACH ROW BEGIN
        if NEW.created_by is null then
               set NEW.created_by = USER();
        end if;
        if NEW.date_created is null then
               set NEW.date_created = CURRENT_TIMESTAMP();
        end if;
  END;
|

CREATE TRIGGER before_update_%s BEFORE UPDATE ON %s
  FOR EACH ROW BEGIN
        set NEW.updated_by = USER();
        set NEW.date_updated = CURRENT_TIMESTAMP();
  END;
|

DELIMITER ;
"""

postgres_function = """
create function set_default_created_by() returns trigger as '
	begin
	if new.created_by is null then
		new.created_by := current_user;
	end if;
	if new.date_created is null then
		new.date_created := current_timestamp;
	end if;
	return new;
	end;
	'
	language plpgsql;

create function set_default_updated_by() returns trigger as '
	begin
	new.updated_by := current_user;
	new.date_updated := current_timestamp;
	return new;
	end;
	'
	language plpgsql;
"""

postgres_trigger = """
create trigger insert_%s before insert on %s
	for each row execute procedure set_default_created_by();

create trigger update_%s before update on %s
	for each row execute procedure set_default_updated_by();
"""
        
class OutputSQLTrigger(object):
	__doc__ = __doc__
	option_default_dict = {('type', 1, int):[1, 'y', 1, 'which type of database? mysql (2) or postgres (1)', ],\
						('table_names', 1, ): [None, 'i', 1, 'coma-separated names of the table', ],\
						('schema', 0, ): [None, 'k', 1, 'set the search_path for psql. only for postgresql.', ],\
						('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
						('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
	def __init__(self,  **keywords):
		"""
		2008-07-27
		"""
		from pymodule import ProcessOptions
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
		
		self.table_names = self.table_names.split(',')
		self.printTriggerDict = {1: self.printPostGresTrigger,\
								2: self.printMySQLTrigger}
		
	def printMySQLTrigger(self, table_name):
		"""
		"""
		print(mysql_trigger%(table_name, table_name, table_name, table_name))
		
	def printPostGresTrigger(self, table_name):
		print(postgres_trigger%(table_name, table_name, table_name, table_name))
	
	def run(self):
		"""
		2008-07-27
		"""
		if self.type==1 and self.schema:
			print("set search_path to %s;"%(self.schema))
		if self.type==1:
			print(postgres_function)
		for table_name in self.table_names:
			self.printTriggerDict[self.type](table_name)

if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = OutputSQLTrigger
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.run()
