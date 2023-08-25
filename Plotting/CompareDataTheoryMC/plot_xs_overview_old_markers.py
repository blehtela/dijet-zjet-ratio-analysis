#!/usr/bin/env python2
#-*- coding:utf-8 -*-


#####################################################
#													#
#	XS overview plots :								#
#	coparing data, theory, MC						#
#	Taking input from the ROOT file					#
#	and from fastNLO tables							#
#													#
# 	This is the original version of the script		#
#	In newer version: only lines for MC, theory		#
#													#
#	Created by B. Schillinger, 21.02.2020			#
#	Last modified: B.Schillinger, 21.02.2020		#
#													#
#####################################################

import argparse
import glob, os, sys
import string
import timeit
import math
import matplotlib as mpl
mpl.use('Cairo')			#Cairo offline backend --> for e.g. pdf or png output
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
from matplotlib.ticker import (FormatStrFormatter, LogFormatter, NullFormatter, ScalarFormatter, AutoMinorLocator, MultipleLocator)
from matplotlib.lines import Line2D
from matplotlib import cm
import numpy as np

import fastnlo
import ROOT
from ROOT import TFile, TH1D, TH1F
from ROOT import gROOT

#dictionary --> ybys bin to color
_ybysbin_col = {'yb0ys0':'forestgreen', 'yb0ys1':'mediumblue', 'yb0ys2':'orange', 'yb1ys0':'firebrick', 'yb1ys1':'deepskyblue', 'yb2ys0':'mediumpurple'}
_ybysbin_marker = {'yb0ys0':"o", 'yb0ys1':"^", 'yb0ys2':"s", 'yb1ys0':"d", 'yb1ys1':"P", 'yb2ys0':"v"}
_ybysbin_label = {'yb0ys0':r'$0 \leq y_b < 1$	$0 \leq y^{\ast} < 1$', 'yb0ys1':r'$0 \leq y_b < 1$	$1 \leq y^{\ast} < 2$', 'yb0ys2':r'$0 \leq y_b < 1$	$2 \leq y^{\ast} < 2.4$', 'yb1ys0':r'$1 \leq y_b < 2$	$0 \leq y^{\ast} < 1$', 'yb1ys1':r'$1 \leq y_b < 2$	$1 \leq y^{\ast} < 2$', 'yb2ys0':r'$2 \leq y_b < 2.4$	$0 \leq y^{\ast} < 1$'}

_ybysbin_xind_up = {'yb0ys0':62, 'yb0ys1':54, 'yb0ys2':40, 'yb1ys0':50, 'yb1ys1':43, 'yb2ys0':37}




#----------------------------------------------------------------------------------------------------------------#
#												  function definitions											 # 
#----------------------------------------------------------------------------------------------------------------#
#reuse functions that I wrote for dijetplot_fig7.py
#function to read histogram from rootfile and store values in arrays
def read_rootfile(rootfile, objname, verbose=False):
	#histo = TH1D(gROOT.FindObject(objname)) #is a histogram

	print("[plot_xs_overview.py]: Reading %s from %s" %(objname, rootfile))
	print("--------------------------------------------------------------------------------")
	histo = rootfile.Get(objname)
	#rootfile.ls()
	nbins = histo.GetNbinsX()
	if(verbose==True): print "nbins: ", nbins
	entrieslist=[]
	errorlist=[]		#need errorbars for datapoints
	x_axis_list=[]	  	#list of bincenters --> will become x-axis (only read once [with the ratio])
	low_edges_list=[]
	up_edges_list=[]
	for j in range(1, nbins+1):
		entrieslist.append(histo.GetBinContent(j))
		if(verbose==True): print "entry %s is %s" %(j, entrieslist[-1])

		#if("ratio" in objname):
		x_axis_list.append(histo.GetBinCenter(j))
		low_edges_list.append(histo.GetXaxis().GetBinLowEdge(j))
		up_edges_list.append(histo.GetXaxis().GetBinUpEdge(j))
		errorlist.append(histo.GetBinError(j))
			
	histentries = np.array(entrieslist)
	entries_err = np.array(errorlist)
	x_axis = np.array(x_axis_list)
	low_bb = np.array(low_edges_list)
	up_bb = np.array(up_edges_list)

	if(verbose==True):
		print "histentries: \n", histentries, "\n"
		print "x_axis: \n", x_axis, "\n"
		print "lower bin bounds: \n", low_bb, "\n"
		print "upper bin bounds: \n", up_bb, "\n"

	return histentries, entries_err, x_axis, low_bb, up_bb
	

#function to read XS from fastnlo table
def read_fnlotable(fnlotable, pdfset, order, scale_var_type):
	fnlo = fastnlo.fastNLOLHAPDF(fnlotable, pdfset)	#takes automatically pdfmember 0 ?
	if(order==0):	#only LO desired
		fnlo.SetContributionON(fastnlo.kFixedOrder, 0, True)
		fnlo.SetContributionON(fastnlo.kFixedOrder, 1, False)
	elif(order==1):	#switch LO and NLO on
		fnlo.SetContributionON(fastnlo.kFixedOrder, 0, True)
		fnlo.SetContributionON(fastnlo.kFixedOrder, 1, True)
	
	fnlo.CalcCrossSection()
	xs_ord = fnlo.GetCrossSection()
	xs_array = np.array(xs_ord)

	#convert from pb to fb (multiply by 1000)
	xs_array = 1000*xs_array


	#should also calculate uncertainties --> PDF unc and Scale unc!!!
	#Calculate Scale Uncertainty (6P variation)
	rel_scale_unc = np.array(fnlo.GetScaleUncertaintyVec(scale_var_type)) 	#for chosen order
	###########
	# structure of rel_scale_unc:
	# rel_scale_unc[0,:] means XS
	# rel_scale_unc[1,:] means rel. unc upwards
	# rel_scale_unc[2,:] means rel. unc downwards
	##########
	
	#calculating absolute scale uncertainty in chosen order
	abs_scale_unc = np.empty([2, len(xs_array)])
	#abs_scale_unc = np.empty([2, xs_array.size])
	abs_scale_unc[0] = np.multiply(xs_array, rel_scale_unc[2])				#absolute scale uncertainty downwards
	abs_scale_unc[1] = np.multiply(xs_array, rel_scale_unc[1])				#absolute scale uncertainty upwards


	#Calculate PDF Uncertainty for chosen pdf
	rel_pdf_unc = np.array(fnlo.GetPDFUncertaintyVec(fastnlo.kLHAPDF6))		#for chosen order
	###########
	# structure of rel_pdf_unc:
	# rel_pdf_unc[0,:] means XS
	# rel_pdf_unc[1,:] means rel. unc upwards
	# rel_pdf_unc[2,:] means rel. unc downwards
	##########

	#calculating absolute PDF uncertainty in chosen order
	abs_pdf_unc = np.empty([2, len(xs_array)])
	abs_pdf_unc[0] = np.multiply(xs_array, rel_pdf_unc[2])					#absolute pdf uncertainty downwards
	abs_pdf_unc[1] = np.multiply(xs_array, rel_pdf_unc[1])					#absolute pdf uncertainty upwards
	



	#return xs_array, entries_err, x_axis	#x_axis should be the same as for the root hists
	return xs_array, abs_scale_unc, abs_pdf_unc



#function to plot a xs in the xs plot
#for individual ybys bins
def plot_xs(bname, ax, x_axis, xs_values, xs_err, m_color, m_size, m_alpha, process, labeling, legbool=False):
##low_bb and up_bb are used for drawing the "x-error" == binwidth

	if(legbool==True):
		labeling = labeling
	else:
		labeling = None

	#plot xs as errorbar plot
	#originally markersize was set to 10.
	ax.errorbar(x_axis, xs_values, xerr=0.0, yerr=xs_err, elinewidth=1.2, linewidth=1.0, marker=_ybysbin_marker[bname], ms=m_size, color=m_color, fillstyle='none', fmt='.', label=labeling)
	ax.grid(True, which="major", axis="y", color="k", linestyle="dotted", linewidth=1)
	#plt.grid(True, which="major", linestyle="dotted", color='k')

	#ax.set_xlim(40, 1200)
	#ax.set_ylim(y_min+0.0000001, y_max+0.1)

	ax.set_xlabel(r'$\mathrm{p_{T,avg}}$ /GeV', fontsize=12, horizontalalignment='right')
	ax.xaxis.set_label_coords(1.00, -0.06)
	#plt.xlabel(r'$\mathrm{p_{T,avg}}$ /GeV', horizontalalignment='right', loc='right')
	ax.set_ylabel(r'XS Data, Theory and MC in fb/GeV', fontsize=12, horizontalalignment='center', verticalalignment='top', y=0.5, rotation=90, labelpad=24)



#-----------------------------------------------------------------------------------------------#
#											  main-function										# 
#-----------------------------------------------------------------------------------------------#

def main():
	# start timer
	# measuring wall clock time - not necessary, just for fun.
	start_time = timeit.default_timer()



		#------------------------------------------------------------------------------#
		#								  parser									   # 
		#------------------------------------------------------------------------------#
		#here: not a function for the parsing job, but part of main-function
		#defining arguments and options

	parser = argparse.ArgumentParser(epilog='', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	# Positional arguments
	parser.add_argument('data_rootfile', type=TFile, nargs='?',
						help='Required argument! File with data results that shall be plotted. (filename glob)')
	parser.add_argument('mc_rootfile', type=TFile, nargs='?',
						help='Required argument! File with MC results that shall be plotted. (filename glob)')
	#parser.add_argument('fnlo_table', nargs='?',
	#					help='Required argument! FastNLO table to be evaluated.')
	parser.add_argument('fnlo_tables', type=str, nargs=6,	#nargs='+'
						help='Required argument! FastNLO tables for the six ybys bins to be evaluated.')

	parser.add_argument('process', type=str, nargs='?', choices=["zjet", "dijet"],
						help="Required argument! Specifies which analysis is considered: zjet for Z+jet analysis; dijet for Dijet analysis.") 

	# Optional arguments
	parser.add_argument('--alphamarker', '-a', type=float, default=1.0, nargs='?', action='store',
						help='If set to True: draw marker in summary plot with given input alpha.')

	parser.add_argument('--pdfset', '-p', default='CT14nlo', nargs='?',
						help='PDF set for fnlo table evaluation. Default is CT14nlo.')

	parser.add_argument('--runperiod', '-r', type=str, nargs='?', default=None, action='store',
						help='Specifies which run period has been used: A, B, C, D')
	parser.add_argument('--generator', '-g', type=str, nargs='?', default=None, action='store',
						help='Specifies which generator has been used: herwig7, pythia8, HT_madgraph_pythia8')
	parser.add_argument('--asymmetric', default=False, action='store_true',
						help='If chosen, use asymmetric (6P) scale variations; otherwise use symmetric ones (2P).')

	parser.add_argument('--ybysbin', '-b', type=str, nargs='?', default=None,
						help="Required argument! Specifies yboost ystar bin name.")

	parser.add_argument('--verbose', '-v', default=False, action='store_true',
						help="Increase verbosity if chosen.")


	if len(sys.argv)==1:
		parser.print_help(sys.stderr)
		sys.exit(1)

	# Parse arguments
	args = vars(parser.parse_args())

	# input filename
	mc_rootfile = args['mc_rootfile']
	data_rootfile = args['data_rootfile']
	#fnlo_table = args['fnlo_table']
	fnlo_tables = args['fnlo_tables']	#this is a list
	process = args['process']			#specifies whether looking at dijet or zjet
	ybys = args['ybysbin']				#name of ybys bin
	m_alpha = args['alphamarker']
	pdfset = args['pdfset']
	runperiod = args['runperiod']	#optional argument, but good for labeling --> to be improved
	generator = args['generator']	#optional argument, pythia or herwig
	verbose = args['verbose']

	#scale variation choice for uncertainty:
	if args['asymmetric']:
		scale_var_type = fastnlo.kAsymmetricSixPoint
		variation_type = 'Scale uncertainty (6P)'
	else:
		scale_var_type = fastnlo.kSymmetricTwoPoint
		variation_type = 'Scale uncertainty (2P)'


	#rootfile = TFile(args['rootfile'])
	print('\n')
	print("[plot_xs_overview.py]: Taking MC input from %s" %mc_rootfile)
	print("[plot_xs_overview.py]: Taking data input from %s" %data_rootfile)
	print("[plot_xs_overview.py]: Taking theory input from %s" %fnlo_tables)


	#prepare the lists (will later be arrays)
	if(ybys==None):
		bin_names = np.array(["yb0ys0", "yb0ys1", "yb0ys2", "yb1ys0", "yb1ys1", "yb2ys0"])
	else:
		bin_names = np.array([ybys])



	#-------------------------------------------------------#
	#	preparing the plot for overview of all 6 ybys bin	#
	#-------------------------------------------------------#
	#for the summary (all bins) of the XS plots
	gs_xs = gridspec.GridSpec(3,3)
	fig_xs = plt.figure(figsize=(8,7), dpi=300)
	#ax1 = plt.subplot(gs[:-1,:])
	ax_xs = plt.subplot(xscale="log", yscale="log")
	#patches = []



	#looping through the ybys bins

	for i in range(bin_names.size): #should usually be 6
		bname = bin_names[i]
		if(verbose==True): print "[plot_xs_overview.py]: Now in bin ", bname
		#reading the data and mc points in the current ybys bin from the rootfile
		#read (statistical?) errors on the histogram entries at the same time
		#furthermore: read x-axis and lower and upper bin bounds
		data, data_err, data_xax, data_low_bb, data_up_bb = read_rootfile(data_rootfile, "CrossSections_byBinwidth/%s/sumhisto_%s_Run%s_%s"%(bname, process, runperiod, bname))
		mc, mc_err, mc_xax, mc_low_bb, mc_up_bb = read_rootfile(mc_rootfile, "CrossSections_byBinwidth/%s/sumhisto_%s_%s_%s"%(bname, process, generator, bname))
		#fnloxs = read_fnlotable(fnlo_table)
		if(verbose==True): print "[plot_xs_overview.py]: Reading fnlo table: ", fnlo_tables[i]
		fnloxs_LO, scale_unc_LO, pdf_unc_LO = read_fnlotable(fnlo_tables[i], pdfset, 0, scale_var_type)
		fnloxs_NLO, scale_unc_NLO, pdf_unc_NLO = read_fnlotable(fnlo_tables[i], pdfset, 1, scale_var_type)


		#as long as we do not calculate pdf and scale uncs
		#TO DO: include ScaleUnc and PDFUnc calculation in read_fnlotable() function!!!
		nentries = fnloxs_LO.size
		fnloLO_err = np.full(nentries,0)
		fnloNLO_err = np.full(nentries,0)


		#check the rel unc by giving output
		if(verbose==True):
			print "Scale unc NLO down: \n", scale_unc_NLO[0]
			print "Scale unc NLO up: \n", scale_unc_NLO[1]

		#calculate error using scale unc and pdf unc (quadratic sum)
		#fnlo_relerrNLO = np.sqrt(pow(rel_scale_unc_NLO, 2), pow(rel_pdf_unc_NLO, 2))
		#fnlo_errNLO = np.multiply(fnloxs_NLO, fnlo_relerrNLO)
		fnlo_errLO_down = np.sqrt(pow(scale_unc_LO[0,:],2)+pow(pdf_unc_LO[0,:],2))		#LO down
		fnlo_errLO_up = np.sqrt(pow(scale_unc_LO[1,:],2)+pow(pdf_unc_LO[1,:],2))		#LO up
		fnlo_errNLO_down = np.sqrt(pow(scale_unc_NLO[0,:],2)+pow(pdf_unc_NLO[0,:],2))	#NLO down
		fnlo_errNLO_up = np.sqrt(pow(scale_unc_NLO[1,:],2)+pow(pdf_unc_NLO[1,:],2))		#NLO up


		#fnloLO_err = zip(fnlo_errLO_down, fnlo_errLO_up)	#would make tuples
		#fnloNLO_err = zip(fnlo_errNLO_down, fnlo_errNLO_up)

		fnloLO_err = np.zeros([2, nentries])
		fnloLO_err[0] = fnlo_errLO_down
		fnloLO_err[1] = fnlo_errLO_up

		fnloNLO_err = np.zeros([2, nentries])
		fnloNLO_err[0] = fnlo_errNLO_down
		fnloNLO_err[1] = fnlo_errNLO_up
		
		

		#convert from pb to fb (multiply by 1000) --> now already handled in read_fnlotable()
		#fnloxs_LO = 1000*fnloxs_LO
		#fnloxs_NLO = 1000*fnloxs_NLO


		#remove the lower and higher bins that are not present in fastnlo tables
		#get upper limit from dictionary _ybysbin_xind_up
		#data_xax = data_xax[7:62]	#bin 7 is included, bin 62 excluded

		if(process=="zjet"):
			data_xax = data_xax[7:_ybysbin_xind_up[bname]]
			data = data[7:_ybysbin_xind_up[bname]]
			data_err = data_err[7:_ybysbin_xind_up[bname]]
			mc_err = mc_err[7:_ybysbin_xind_up[bname]]
			mc = mc[7:_ybysbin_xind_up[bname]]

		else:
			#-----------------------------------------------------------------------------------#
			#	Remove the low-pTavg bins, that are empty due to dijet trigger (ptavg < 56GeV)	#
			#	Additionally: restrictions from fnlo tables (see above in original version)		#
			#	--> numpy array slicing															#
			#-----------------------------------------------------------------------------------#
			data_xax = data_xax[11:_ybysbin_xind_up[bname]]
			data = data[11:_ybysbin_xind_up[bname]]
			data_err = data_err[11:_ybysbin_xind_up[bname]]
			mc_err = mc_err[11:_ybysbin_xind_up[bname]]
			mc = mc[11:_ybysbin_xind_up[bname]]
			fnloxs_LO = fnloxs_LO[4:]
			fnloxs_NLO = fnloxs_NLO[4:]
			fnloLO_err = fnloLO_err[:,4:]
			fnloNLO_err = fnloNLO_err[:,4:]



		#-----------------------------------------------#
		#		plotting the XS of current ybys bin		#
		#-----------------------------------------------#
	#	plot_xs(bname,ax_xs, data_xax, data, data_err, _ybysbin_col[bname], m_alpha, process, "Data Run%s in %s"%(runperiod, _ybysbin_label[bname]), True) #data
		plot_xs(bname, ax_xs, data_xax, mc, mc_err, 'saddlebrown', 4, m_alpha, process, "MC %s in %s"%(generator, bname), False)#mc
		plot_xs(bname, ax_xs, data_xax, fnloxs_LO, fnloLO_err, 'gray', 4, m_alpha, process, "LO Table %s in %s"%(pdfset, bname), False)#fnlo LO
		plot_xs(bname, ax_xs, data_xax, fnloxs_NLO, fnloNLO_err, 'black', 4, m_alpha, process, "NLO Table %s in %s"%(pdfset, bname), False)#fnlo LO
		plot_xs(bname,ax_xs, data_xax, data, data_err, _ybysbin_col[bname], 6, m_alpha, process, "Data Run%s in %s"%(runperiod, _ybysbin_label[bname]), True) #data

		#ax_xs.set_title("%s XS in %s (Data Run%s, MC %s, FixedOrder Theory with %s)"%(process, bname, runperiod, generator, pdfset))
		ax_xs.set_title("%s XS overview (Data Run%s, MC %s, FixedOrder Theory with %s)"%(process, runperiod, generator, pdfset), fontsize=11)







	#save the figures: should be only one in total
	#-------------------#
	# legend and save 	#
	#-------------------#
	#create extra legend for MC and theory
	legend_elements = [Line2D([0],[0],color='saddlebrown', lw=4, label="MC %s in each ybys bin"%(generator)), Line2D([0],[0],color='gray', lw=4, label="LO Table %s in each ybys bin"%(pdfset)), Line2D([0],[0],color='black', lw=4, label="NLO Table %s in each ybys bin"%(pdfset))]
	
	#save plots for current ybys bin
	#xs plot:
	#plt.legend()
	ax_xs.legend(fontsize=10, numpoints=1, loc='lower left')
	first_legend = ax_xs.get_legend()
	ax_xs.add_artist(first_legend)
	second_legend = plt.legend(handles=legend_elements, loc='upper right')


	xlocs = ax_xs.get_xticks(minor=True)
	xlabels = ax_xs.get_xticklabels(minor=True)
	print "xlocs: "
	print xlocs
	print "xlabels: "
	print xlabels

	#xlocs --> 48 minor label locations
	#xlabels --> 48 minor labels

	x_minticks = [20, 30, 40, 50, 60, 70, 80, 90, 200, 300, 400, 500, 600, 700, 800, 900, 2000, 3000, 4000]
	x_minticklabels = ['','30','','','60', '', '', '', '200', '', '', '500','', '700', '', '', '2000', '3000']
	x_majticks = [100, 1000]
	x_majticklabels = ['100', '1000']
	ax_xs.set_xticks(x_majticks, minor=False)	#set major ticks loc
	ax_xs.set_xticks(x_minticks, minor=True)	#set minor ticks loc

	ax_xs.set_xticklabels(x_majticklabels, minor=False)
	ax_xs.set_xticklabels(x_minticklabels, minor=True)

	plt.tight_layout()
	fig_xs.savefig("%s_xs_overview_theory_Run%s_%s_NotShifted.png"%(process, runperiod, generator))


	# stop timer
	stop_time = timeit.default_timer()
	timediff = stop_time-start_time
	print("[plot_xs_overview.py]: Elapsed time: %s sec = %s min" %(timediff, round(timediff/60., 2)))



########################################################
## function to make better uncertainty plots (shaded)
## could maybe be replaced by .flatten() ??
def steppify_bin(arr, isx=False):
	""" 
	Produce stepped array of arr, needed for example for stepped fill_betweens. 
	Pass all x bin edges to produce stepped x arr and all y bincontents to produce 
	stepped bincontents representation 
	steppify_bin([1,2,3], True)  
	-> [1,2,2,3] 
	steppify_bin([5,6]) 
	-> [5,5,6,6] 
	""" 
	if isx:
		newarr = np.array(zip(arr[:-1], arr[1:])).ravel()
	else:
		newarr = np.array(zip(arr, arr)).ravel()

	return newarr

###########################################################



if __name__ == "__main__":
	main()








