#!/usr/bin/env python
"""
Example of embedding matplotlib in an application and interacting with a treeview to store data.  Double click on an entry to update plot data

2015.04.16 input matrix file should have a header.
	By default first two columns are of string type, thereafter are of float type.
	
	A user can also append something like "|float" to each header to specify type.
	For example:
		name|string	position|int	chr|str	value|float
	
	Numeric value can not be empty. By default -1000 (changeable on the fly) is treated as NA in plotting.
	
	float/number/double/numeric all means float.
	str/string all means string.
"""
import __init__	#used to know the path to this file itself
import os, sys, pygtk
pygtk.require('2.0')
import gtk, gtk.glade, gobject
from gtk import gdk
import gnome
import gnome.ui
import math, random, re, copy
import matplotlib
matplotlib.use('GTKAgg')  # or 'GTK'
from matplotlib.backends.backend_gtkagg import FigureCanvasGTKAgg as FigureCanvas

from matplotlib.figure import Figure
from matplotlib.backends.backend_gtkagg import NavigationToolbar2GTKAgg as NavigationToolbar

import yh_gnome, csv, traceback
from pymodule.yhio.SNP import SNPData, read_data
from pymodule.utils import figureOutDelimiter
from pymodule import MatrixFile
from pymodule.plot import yh_matplotlib
from variation.src.qc.FilterStrainSNPMatrix import FilterStrainSNPMatrix

class ValuePreProcessor(object):
	"""
	2014.07.31
	"""
	def __init__(self, na=None, scalar=None, addition=None, logScale=None, errorColumnIndex=None):
		self.na = na
		self.scalar = scalar
		self.addition = addition
		self.logScale = logScale
		self.errorColumnIndex = errorColumnIndex
		

class DataMatrixGuiXYProbe(gtk.Window):
	"""
	2009-3-13
		migrated from QCVisualize.py. now become a standalone program and able to read data from a file and plot ...
		QCVisualize.py inherits from here
	2008-02-05
		embed it into a bigger gnome app, add more buttons, and change the __init__()
	2008-01-01
		class to visualize the results from QualityControl.py
	"""
	def __init__(self, plot_title='', id_is_strain=1, header=None, strain_acc_list=None, category_list=None, data_matrix=None):
		"""
		2008-01-10
			use a paned window to wrap the scrolledwindow and the canvas
			so that the relative size of canvas to the scrolledwindow could be adjusted by the user.
		"""
		prog = gnome.program_init('DataMatrixGuiXYProbe', '0.1')	#this must be called before any initialization for gnome app
		
		program_path = os.path.dirname(__init__.__file__)	#sys.argv[0])
		xml = gtk.glade.XML(os.path.join(program_path, 'DataMatrixGuiXYProbe.glade'))
		xml.signal_autoconnect(self)
		self.app1 = xml.get_widget("app1")
		self.app1.connect("delete_event", gtk.main_quit)
		self.app1.set_default_size(800, 1000)
		self.app1.set_title(plot_title)
		
		self.plot_title = plot_title
		self.id_is_strain = id_is_strain
		self.header = header
		self.strain_acc_list = strain_acc_list
		self.category_list = category_list
		self.data_matrix = data_matrix
		
		self.inputDataHeaders = None
		
		self.columnTypes = None
		self.columnHeaders = None
		self.columnEditableFlagList = None
		self.list2D = None
		
		self.vbox1 = xml.get_widget("vbox1")
		self.treeview_matrix = xml.get_widget("treeview_matrix")
		
		# matplotlib stuff
		fig = Figure(figsize=(8,8))
		self.canvas = FigureCanvas(fig)  # a gtk.DrawingArea
		self._idClick = self.canvas.mpl_connect('button_press_event', self.onUserClickCanvas)
		self.vpaned1 = xml.get_widget("vpaned1")
		self.vpaned1.add2(self.canvas)
		
		#vbox.pack_start(self.canvas, True, True)
		self.ax = fig.add_subplot(111)
		self.treeview_matrix.connect('row-activated', self.plot_row)
		
		toolbar = NavigationToolbar(self.canvas, self.app1)
		self.vbox1.pack_start(toolbar, False, False)
		
		self.checkbutton_label_dot = xml.get_widget('checkbutton_label_dot')
		self.entry_dot_label_column = xml.get_widget('entry_dot_label_column')
		self.dataLabelColumnIndexAndSeparatorList = None
		self.dataLabelNumericItemIndexList = None 
		
		self.entry_x_na = xml.get_widget('entry_x_na')
		self.entry_y_na = xml.get_widget('entry_y_na')
		
		self.entry_multiply_x = xml.get_widget('entry_multiply_x')
		self.entry_multiply_y = xml.get_widget('entry_multiply_y')
		self.entry_add_x = xml.get_widget('entry_add_x')
		self.entry_add_y = xml.get_widget('entry_add_y')

		self.entry_x_error = xml.get_widget("entry_x_error")
		self.entry_y_error = xml.get_widget("entry_y_error")
		self.checkbutton_logX = xml.get_widget('checkbutton_logX')	#2014.06.09
		self.checkbutton_logY = xml.get_widget('checkbutton_logY')
		self.checkButtonPlotOnlySelected = xml.get_widget('checkButtonPlotOnlySelected')	#2015.04.16
		self.entry_x_column = xml.get_widget('entry_x_column')
		self.entry_y_column = xml.get_widget('entry_y_column')
		
		self.entry_filters = xml.get_widget("entry_filters")
		#self.checkbutton_histLogX = xml.get_widget('checkbutton_histLogX')	#2014.06.09
		#self.checkbutton_histLogY = xml.get_widget('checkbutton_histLogY')
		
		self.entry_hist_column = xml.get_widget('entry_hist_column')
		self.entry_no_of_bins = xml.get_widget('entry_no_of_bins')	#2009-5-20
		self.entry_plot_title = xml.get_widget('entry_plot_title')
		self.entry_plot_title.set_text(self.plot_title)
		
		self.filechooserdialog_save = xml.get_widget("filechooserdialog_save")
		self.filechooserdialog_save.connect("delete_event", yh_gnome.subwindow_hide)
		
		self.entry_sampling_probability = xml.get_widget("entry_sampling_probability")
		self.filechooserdialog_open = xml.get_widget("filechooserdialog_open")
		self.filechooserdialog_open.connect("delete_event", yh_gnome.subwindow_hide)
		
		self.app1_appbar1 = xml.get_widget('app1_appbar1')
		self.app1_appbar1.push('Status Message.')	#import gnome.ui has to be executed.
		
		self.treeview_matrix.connect('cursor-changed', self.update_no_of_selected, self.app1_appbar1)
		self.app1.show_all()
		
		self.xValuePreProcessor = None
		self.yValuePreProcessor = None
		
		self.x_error_column_index = None
		self.y_error_column_index = None
		
		self.typeName2PythonType = {
				"str": str,
				"string": str,
				"numeric":float,
				"number":float,
				"double":float,
				"float":float,
				"integer":int,
				"int":int}
		
		#self.add_events(gdk.BUTTON_PRESS_MASK|gdk.KEY_PRESS_MASK|gdk.KEY_RELEASE_MASK)
	
	def parseDataLabelColumns(self, inputText):
		"""
		2015.04.16
		"""
		splitP= re.compile(r'([,/|\.\-_=\?:;"\'^%$@+])')	#any single character included could be used as splitter
		self.dataLabelColumnIndexAndSeparatorList = splitP.split(inputText)
		self.dataLabelNumericItemIndexList = []
		for i in range(len(self.dataLabelColumnIndexAndSeparatorList)):
			if not splitP.match(self.dataLabelColumnIndexAndSeparatorList[i]):	#it's a column index
				self.dataLabelColumnIndexAndSeparatorList[i] = int(self.dataLabelColumnIndexAndSeparatorList[i])
				self.dataLabelNumericItemIndexList.append(i)
		return self.dataLabelColumnIndexAndSeparatorList, self.dataLabelNumericItemIndexList 
	
	def getDataPointLabel(self, dataRow):
		"""
		2015.04.16
		"""
		dataLabelInList = copy.deepcopy(self.dataLabelColumnIndexAndSeparatorList)
		for i in self.dataLabelNumericItemIndexList:
			dataLabelInList[i] = str(dataRow[self.dataLabelColumnIndexAndSeparatorList[i]])
		return ''.join(dataLabelInList)
	
	def onUserClickCanvas(self, event):
		"""
		2016.06.07 bugfix. zero(0) generates false condition while only empty string or None should generate false.
		2009-3-13
			use (x_lim[1]-x_lim[0])/200. as the resolution for a dot to be called identical to a data point.
			similar for the y_data
		2009-3-13
			deal with checkbutton_label_dot, entry_dot_label_column, entry_x_column, entry_y_column
		2008-01-01
			derived from on_click_row() of QualityControl.py
			reaction when user clicked in the plot
		"""
		# get the x and y coords, flip y from top to bottom
		x, y = event.x, event.y
		toLabelDataPoint = self.checkbutton_label_dot.get_active()
		
		self.parseDataLabelColumns(self.entry_dot_label_column.get_text())
		
		x_column = int(self.entry_x_column.get_text())
		y_column = int(self.entry_y_column.get_text())
		x_lim = self.ax.get_xlim()
		x_grain_size = abs(x_lim[1]-x_lim[0])/200.
		y_lim = self.ax.get_ylim()
		y_grain_size = abs(y_lim[1]-y_lim[0])/200.
		if event.button==1:
			if event.inaxes is not None:
				sys.stderr.write("x_grain_size: %s, y_grain_size: %s\n"%(x_grain_size, y_grain_size))
				sys.stderr.write("data coords: %s, %s\n"%(event.xdata, event.ydata))
				if self.list2D is None:
					return
				for row in self.list2D:
					if row[x_column] is not None and row[x_column]!="" and row[y_column] is not None and row[y_column]!="" :	#not empty
						#2016.06.07
						try:
							x_data = float(row[x_column])
							y_data = float(row[y_column])
							x = self.processDataValue(x_data, self.xValuePreProcessor)
							if x is None:
								continue
							y = self.processDataValue(y_data, self.yValuePreProcessor)
							if y is None:
								continue
							if abs(x-event.xdata)<x_grain_size and abs(y-event.ydata)<y_grain_size:
								dataLabel = self.getDataPointLabel(row)
								if toLabelDataPoint:
									self.ax.text(event.xdata, event.ydata, dataLabel, size=8)
									self.canvas.draw()
								sys.stderr.write("%s: %s, %s: %s, raw xy=(%s, %s), scaled xy=(%s,%s), dataLabel: %s.\n"%\
												(self.columnHeaders[0], row[0], self.columnHeaders[1], row[1], 
												row[x_column], row[y_column], x,y, dataLabel))
						except:
							sys.stderr.write("Column %s, %s of row (%s), could not be converted to float. skip.\n"%\
											(x_column, y_column, repr(row)))
	
	def setUserDataPreprocessingFlags(self):
		"""
		2014.07.25
		"""
		
		self.xValuePreProcessor = ValuePreProcessor(na = self.entry_x_na.get_text())
		self.yValuePreProcessor = ValuePreProcessor(na = self.entry_y_na.get_text())
		
		if self.entry_multiply_x.get_text():
			self.xValuePreProcessor.scalar = float(self.entry_multiply_x.get_text())
		if self.entry_multiply_y.get_text():
			self.yValuePreProcessor.scalar = float(self.entry_multiply_y.get_text())

		if self.entry_add_x.get_text():
			self.xValuePreProcessor.addition = float(self.entry_add_x.get_text())
		if self.entry_add_y.get_text():
			self.yValuePreProcessor.addition = float(self.entry_add_y.get_text())
		
		if self.entry_x_error.get_text():
			self.xValuePreProcessor.errorColumnIndex = int(self.entry_x_error.get_text())
		if self.entry_y_error.get_text():
			self.yValuePreProcessor.errorColumnIndex = int(self.entry_y_error.get_text())
		
		if self.checkbutton_logX.get_active():
			self.xValuePreProcessor.logScale = True
		if self.checkbutton_logY.get_active():
			self.yValuePreProcessor.logScale = True
	
	def processDataValue(self, value=None, valuePreProcessor=None):
		"""
		2014.07.31
		"""
		
		if valuePreProcessor.na is not None and (value==valuePreProcessor.na or float(value)==float(valuePreProcessor.na)):
			return None
		value = float(value)
		if valuePreProcessor.scalar is not None:
			value = value*valuePreProcessor.scalar
		if valuePreProcessor.addition is not None:
			value = value + valuePreProcessor.addition
		return value
	
	def decorateAxisLabel(self, label=None, valuePreProcessor=None):
		"""
		2014.07.31
		"""
		if valuePreProcessor.scalar is not None:
			label = "%s*%s"%(valuePreProcessor.scalar, label)
		if valuePreProcessor.addition:
			label = "%s+%s"%(label, valuePreProcessor.addition)
		return label
	
	def filterDataRow(self, dataRow):
		"""
		2015.4.16
			unfinished
		"""
		logicSplitP= re.compile(r'(AND|OR)')
		equationSplitP = re.compile(r'>=|>|=|<|<=')
		filtersText = self.entry_filters.get_text()
		logicSplitP.split(filtersText)
		
		self.dataLabelColumnIndexAndSeparatorList = splitP.split(inputText)
		self.dataLabelNumericItemIndexList = []
		for i in range(len(self.dataLabelColumnIndexAndSeparatorList)):
			if not splitP.match(self.dataLabelColumnIndexAndSeparatorList[i]):	#it's a column index
				self.dataLabelColumnIndexAndSeparatorList[i] = int(self.dataLabelColumnIndexAndSeparatorList[i])
				self.dataLabelNumericItemIndexList.append(i)

	
	def plotXY(self, ax, canvas, liststore, plot_title='', 
			chosen_index_ls=[]):
		"""
		2015.01.28 add summary stats to title
		2014.04.29 add error bars
		2009-3-13
			rename plot_NA_mismatch_rate to plotXY()
		2008-02-05
			chosen_index => chosen_index_ls
		2007-12-14
		"""
		x_column = int(self.entry_x_column.get_text())
		y_column = int(self.entry_y_column.get_text())
		self.setUserDataPreprocessingFlags()   
		
		plot_title = self.entry_plot_title.get_text()
		
		min_x = 1
		min_y = 1
		max_x = 0
		max_y = 0
		
		x_ls = []
		x_error_ls = []
		y_ls = []
		y_error_ls = []
		
		x_chosen_ls = []
		x_chosen_error_ls = []
		y_chosen_ls = []
		y_chosen_error_ls = []
		
		noOfValuesDiscarded = 0
		chosen_index_set = set(chosen_index_ls)
		for i in range(len(liststore)):
			row = liststore[i]
			x = row[x_column]
			y = row[y_column]
			if x=='' or y=='':	#2015.04.16 bugfix 2013.07.12 skip if empty cells
				continue
				noOfValuesDiscarded += 1
			x = self.processDataValue(x, self.xValuePreProcessor)
			if x is None:
				continue
				noOfValuesDiscarded += 1
			y = self.processDataValue(y, self.yValuePreProcessor)
			if y is None:
				continue
				noOfValuesDiscarded += 1
			#self.filterDataRow()
			if self.xValuePreProcessor.errorColumnIndex is not None:
				x_error = row[self.xValuePreProcessor.errorColumnIndex]
			else:
				x_error = 0
			if self.yValuePreProcessor.errorColumnIndex is not None:
				y_error = row[self.yValuePreProcessor.errorColumnIndex]
			else:
				y_error = 0

			if x<min_x:
				min_x = x
			if x>max_x:
				max_x = x
			if y<min_y:
				min_y = y
			if y>max_y:
				max_y = y
			
			if i in chosen_index_set:
				x_chosen_ls.append(x)
				y_chosen_ls.append(y)
				x_chosen_error_ls.append(x_error)
				y_chosen_error_ls.append(y_error)
			else:
				x_ls.append(x)
				y_ls.append(y)
				x_error_ls.append(x_error)
				y_error_ls.append(y_error)
		
		sys.stderr.write("WARNING: %s values were discarded.\n"%(noOfValuesDiscarded))
		ax.clear()
		if self.xValuePreProcessor.logScale:
			ax.set_xscale('log')
		if self.yValuePreProcessor.logScale:
			ax.set_yscale('log')
		
		#plot unselected data only if 
		if not self.checkButtonPlotOnlySelected.get_active():
			if self.x_error_column_index is not None and self.y_error_column_index is not None:
				ax.errorbar(x_ls, y_ls, xerr=x_error_ls, yerr=y_error_ls, ecolor='g', fmt='o')
			else:
				ax.plot(x_ls, y_ls, 'o', alpha=0.7, linewidth=0, linestyle="None")
		
		
		"""
		#diagonal line give a rough feeling about the notion, more NA, worse calling
		diagonal_start = min(min_x, min_y)-0.1
		diagonal_end = max(max_x, max_x)+0.1
		ax.plot([diagonal_start, diagonal_end],[diagonal_start, diagonal_end])
		"""
		if x_chosen_ls and y_chosen_ls:	#highlight
			titleWithStats = "Highlighted data\n" + yh_matplotlib.constructTitleFromTwoDataSummaryStat(x_chosen_ls, y_chosen_ls)
			
			ax.plot(x_chosen_ls, y_chosen_ls, 'o', alpha=0.6,  c='r')
			if self.x_error_column_index is not None and self.y_error_column_index is not None:
				ax.errorbar(x_chosen_ls, y_chosen_ls, xerr=x_chosen_error_ls, yerr=y_chosen_error_ls, ecolor='r', color='r', fmt='o')
		else:	#take all data
			titleWithStats = yh_matplotlib.constructTitleFromTwoDataSummaryStat(x_ls+x_chosen_ls, y_ls+y_chosen_ls)
		if plot_title:
			ax.set_title("%s %s"%(plot_title, titleWithStats))
		else:
			ax.set_title(titleWithStats)
		
		xlabel = "(%s)"%self.columnHeaders[x_column]
		xlabel = self.decorateAxisLabel(xlabel, self.xValuePreProcessor)
		ax.set_xlabel(xlabel)
		ylabel = "(%s)"%self.columnHeaders[y_column]
		ylabel = self.decorateAxisLabel(ylabel, self.yValuePreProcessor)
		ax.set_ylabel(ylabel)
		canvas.draw()
	
	def plot_row(self, treeview, path, view_column):
		if self._idClick==None:
			self._idClick = self.canvas.mpl_connect('button_press_event', self.onUserClickCanvas)
		self.plotXY(self.ax, self.canvas, self.liststore, self.plot_title, path)
	
	def setupColumns(self, treeview):
		"""
		2009-3-13
		"""
		if not getattr(self, 'columnHeaders', None):
			sys.stderr.write("Nothing in columns yet.\n")
			return
		self.liststore = gtk.ListStore(*self.columnTypes)
		#self.add_columns(self.treeview_matrix)
		yh_gnome.create_columns(self.treeview_matrix, self.columnHeaders, self.columnEditableFlagList, self.liststore)
		yh_gnome.fill_treeview(self.treeview_matrix, self.liststore, self.list2D, reorderable=True)
		self.treeselection = self.treeview_matrix.get_selection()
	
	def on_button_PlotXY_clicked(self, widget, data=None):
		"""
		2008-02-12
		to update the no_of_selected rows (have to double click a row to change a cursor if it's multiple selection)
		2008-02-05
		"""
		if self._idClick==None:
			self._idClick = self.canvas.mpl_connect('button_press_event', self.onUserClickCanvas)
		pathlist_strains1 = []
		self.treeselection.selected_foreach(yh_gnome.foreach_cb, pathlist_strains1)
		index_ls = []
		for path in pathlist_strains1:
			index_ls.append(path[0])
		self.app1_appbar1.push("%s rows selected."%len(pathlist_strains1))
		self.plotXY(self.ax, self.canvas, self.liststore, self.plot_title, index_ls)
	
	def on_button_UnSelectAll_clicked(self, widget, data=None):
		"""
		2015.04.16
		"""
		self.treeselection = self.treeview_matrix.get_selection()
		self.treeselection.unselect_all()

		
	def on_button_save_clicked(self, widget, data=None):
		"""
		2008-02-05
		"""
		self.filechooserdialog_save.show_all()
	
	def on_button_filechooserdialog_cancel_ok_clicked(self, widget, data=None):
		"""
		2008-02-05
		"""
		self.filechooserdialog_save.hide()
	
	def on_button_filechooserdialog_save_ok_clicked(self, widget, data=None):
		"""
		2008-02-12
		to update the no_of_selected rows (have to double click a row to change a cursor if it's multiple selection)
		2008-02-05
		"""
		output_fname = self.filechooserdialog_save.get_filename()
		self.filechooserdialog_save.hide()
		pathlist_strains1 = []
		self.treeselection.selected_foreach(yh_gnome.foreach_cb, pathlist_strains1)
		self.app1_appbar1.push("%s rows selected."%len(pathlist_strains1))
		if self.header and self.strain_acc_list and self.category_list and self.data_matrix:
			selected_index_set = set()
			for path in pathlist_strains1:
				row = self.liststore[path[0]]
				id = row[0]
				index_in_data_matrix = row[-1]
				selected_index_set.add(index_in_data_matrix)
				if self.id_is_strain:
					id = id[1:-1].split(',')	#id is a tuple of (ecotypeid,duplicate)
					self.strain_acc_list[index_in_data_matrix] = id[0].strip()	#remove extra space
					self.category_list[index_in_data_matrix] = id[1].strip()
				#else:
				#	self.header[index_in_data_matrix+2] = id
			FilterStrainSNPMatrix_instance = FilterStrainSNPMatrix()
			if self.id_is_strain:
				rows_to_be_tossed_out = set(range(len(self.strain_acc_list))) - selected_index_set
				FilterStrainSNPMatrix_instance.write_data_matrix(self.data_matrix, output_fname, self.header, self.strain_acc_list, self.category_list,\
								rows_to_be_tossed_out, cols_to_be_tossed_out=set(), nt_alphabet=0)
			else:
				cols_to_be_tossed_out = set(range(len(self.header)-2)) - selected_index_set
				FilterStrainSNPMatrix_instance.write_data_matrix(self.data_matrix, output_fname, self.header, self.strain_acc_list, self.category_list,\
								rows_to_be_tossed_out=set(), cols_to_be_tossed_out=cols_to_be_tossed_out, nt_alphabet=0)
	
	def show_all(self):
		"""
		2008-02-05
			preserve the old interface. in order not to change anything in plot_col_NA_mismatch_rate() and plot_row_NA_mismatch_rate() of QualityControl.py
		"""
		self.app1.show_all()
	
	def on_button_PlotHistogram_clicked(self, widget, data=None):
		"""
		2016.04.14 bugfix: skip "" /empty cells
		2015.01.28 add summary stats to title
		2009-5-20
			get the number of bins from entry_no_of_bins 
		2009-3-13
			draw histogram of specific hist_column
		2008-02-06
		"""
		if not getattr(self, 'columnHeaders', None):
			sys.stderr.write("Nothing in columns yet.\n")
			return
		self.setUserDataPreprocessingFlags()
		
		self.ax.clear()
		self.canvas.mpl_disconnect(self._idClick)	#drop the signal handler
		self._idClick = None	#reset the _idClick
		hist_ls = []
		hist_column = int(self.entry_hist_column.get_text())
		noOfValuesDiscarded = 0
		for i in range(len(self.liststore)):
			x = self.liststore[i][hist_column]
			if x=='':	#2016.04.14 bugfix skip if empty cells  #if not x:
				noOfValuesDiscarded += 1
				continue
			x = self.processDataValue(x, self.xValuePreProcessor)
			if x is None:
				continue
				noOfValuesDiscarded += 1
			if self.xValuePreProcessor.logScale:
				if x>0:
					x = math.log10(x)
				else:
					sys.stderr.write("x value %s, not good for log10.\n"%(x))
					continue
					noOfValuesDiscarded += 1
			hist_ls.append(x)
		sys.stderr.write("WARNING: %s values were discarded.\n"%(noOfValuesDiscarded))
		
		title = "%s %s %s"%(self.plot_title, self.columnHeaders[hist_column],
				yh_matplotlib.constructTitleFromDataSummaryStat(hist_ls))
		self.ax.set_title(title);	#"Histogram of %s %s"%(self.plot_title, self.columnHeaders[hist_column]))
		no_of_bins = int(self.entry_no_of_bins.get_text())
		
		#if self.x_logScale:
		#	self.ax.set_xscale('log')
		if self.yValuePreProcessor.logScale:
			self.ax.set_yscale('log')
		
		xlabel = "(%s)"%self.columnHeaders[hist_column]
		xlabel = self.decorateAxisLabel(xlabel, self.xValuePreProcessor)
		if self.xValuePreProcessor.logScale:
			xlabel = "log10(%s)"%(xlabel)
		self.ax.set_xlabel(xlabel)
		self.ax.hist(hist_ls, no_of_bins)
		self.canvas.draw()
	
	def update_no_of_selected(self, treeview, app1_appbar1):
		"""
		2008-02-12
			to update the no_of_selected rows (have to double click a row to change a cursor if it's multiple selection)
		"""
		pathlist_strains1 = []
		self.treeselection.selected_foreach(yh_gnome.foreach_cb, pathlist_strains1)
		app1_appbar1.push("%s rows selected."%len(pathlist_strains1))
		return True
	
	def parseDataHeader(self, dataHeaders=None):
		"""
		2016.04.15 inserted a first column to denote order of data
		2015.04.16
		"""
		no_of_cols = len(dataHeaders)
		self.columnHeaders = ['']*(no_of_cols+1)
		self.columnHeaders[0] = "0 Order"
		self.columnTypes = [str]*(no_of_cols+1)
		self.columnTypes[0] = int
		self.columnEditableFlagList = [False]*(no_of_cols+1)
		
		for i in range(no_of_cols):
			header = dataHeaders[i]
			tmp_ls = header.split('|')
			columnHeader = tmp_ls[0]
			if len(tmp_ls)>1:
				columnTypeName = tmp_ls[1]
			elif i<2:	#by default first two columns are of string type
				columnTypeName = "string"
			else:	#by default columns after first two are of numeric type
				columnTypeName = 'number'
			column_type = self.typeName2PythonType.get(columnTypeName, str)
			self.columnHeaders[i+1] = '%s %s'%(i+1, columnHeader)
			self.columnTypes[i+1] = column_type
			if column_type==str:
				self.columnEditableFlagList[i+1] = True
		
	
	def readInDataToPlot(self, input_fname, sampling_probability=1.0):
		"""
		2016.04.15 report summary of data and inserted a first column to denote order of data
		2015.04.16 use parseDataHeader()
			convert each column data according to self.columnTypes
		2015.01.23 added argument sampling_probability to sub-sample data
		2013.07.11 use MatrixFile to read in the file
		2009-5-20
			add the column index into the column header for easy picking
		2009-3-13
			wrap the float conversion part into try...except to report what goes wrong
		2009-3-13
		"""
		if sampling_probability>1 or sampling_probability<0:
			sampling_probability=1.0
		reader = MatrixFile(inputFname=input_fname)
		self.inputDataHeaders = reader.next()
		self.parseDataHeader(self.inputDataHeaders)
		
		self.list2D = []
		dimOfRawData = [0,0]
		dimOfList2D = [0,0]
		
		for row in reader:
			dimOfRawData[0] += 1
			dimOfRawData[1] = len(row)
			
			if sampling_probability>0 and sampling_probability<1:
				if random.random()>sampling_probability:	#skip
					continue
			new_row = ['']*(len(row)+1)	#first column is the order of data
			for i in range(len(row)):
				try:
					new_row[i+1] = self.columnTypes[i+1](row[i])
				except:
					sys.stderr.write("Error in converting column %s data %s to type %s.\n"%(i, row[i], self.columnTypes[i+1]))
					sys.stderr.write('Except type: %s\n'%repr(sys.exc_info()))
					traceback.print_exc()
			
			dimOfList2D[0] += 1
			dimOfList2D[1] = len(new_row)
			new_row[0] = dimOfList2D[0]	#order of this data
			self.list2D.append(new_row)
			
		reader.close()
		sys.stderr.write("Dimension of raw data: %s. Dimension of displayed data: %s."%(dimOfRawData, dimOfList2D))
		self.setupColumns(self.treeview_matrix)
		#update status to reflect the input filename
		self.app1.set_title(os.path.basename(input_fname))
		self.app1_appbar1.push(input_fname)
		self.plotXY(self.ax, self.canvas, self.liststore, self.plot_title)
		
	
	def readInRawMatrixData(self, input_fname):
		"""
		2009-3-13
		"""
		delimiter = figureOutDelimiter(input_fname)
		self.header, self.strain_acc_list, self.category_list, self.data_matrix = read_data(input_fname, delimiter=delimiter)
		
	def on_imagemenuitem_open_activate(self, widget, data=None):
		"""
		2009-3-13
		"""
		self.filechooserdialog_open.show_all()
	
	def on_button_fileopen_cancel_clicked(self, widget, data=None):
		"""
		2015.01.23
		"""
		self.filechooserdialog_open.hide()
		
	
	def on_button_fileopen_ok_clicked(self, widget, data=None):
		"""
		2009-3-13
		"""
		input_fname = self.filechooserdialog_open.get_filename()
		sampling_probability = float(self.entry_sampling_probability.get_text())
		self.filechooserdialog_open.hide()
		self.readInDataToPlot(input_fname, sampling_probability)
	
	def on_entry_plot_title_change(self, widget, data=None):
		"""
		2009-3-13
			upon any change in the entry_plot_title
		"""
		self.plot_title = self.entry_plot_title.get_text()
	
if __name__ == '__main__':
	prog = gnome.program_init('DataMatrixGuiXYProbe', '0.1')
	instance = DataMatrixGuiXYProbe()
	gtk.main()
