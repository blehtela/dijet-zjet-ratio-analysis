#!/usr/bin/env python2
#-*- coding:utf-8 -*-

#############################################################
#														   	#
# 	Plotting the Dijet / Z+jet ratios						#
# 	Taking input from analysis carried out within ROOT		#
# 	This is based on "ratios.py", but only plotting a		#
#	summary plot containing all the ybys bins.				#
#															#
#	This is an updated version of ratios_root.py			#															
#	It contains the following changes:						#
#	> Only display the ratio-range available in data		#
#	> Changed axis ranges accordingly						#
#	> Plus some smaller adjustments							#
#															#
#	--> Removed unused and commented-out lines				#
#	--> Only create a summary plot							#
#													   		#
# 	Created by B. Schillinger, 29.08.2019				 	#
# 	Last modified by B. Schillinger, 24.03.2020		   		#
#													   		#
#############################################################

import argparse
import glob, os, sys
import string
import timeit
import matplotlib as mpl
mpl.use('Cairo')			#Cairo offline backend --> for e.g. pdf or png output
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
from matplotlib.ticker import (FormatStrFormatter, LogFormatter, NullFormatter, ScalarFormatter, AutoMinorLocator, MultipleLocator)
from matplotlib import cm
import numpy as np

#import fastnlo
import ROOT
from ROOT import TFile, TH1D, TH1F
from ROOT import gROOT

#dictionary --> ybys bin to color
_ybysbin_col = {'yb0ys0':'forestgreen', 'yb0ys1':'mediumblue', 'yb0ys2':'orange', 'yb1ys0':'firebrick', 'yb1ys1':'deepskyblue', 'yb2ys0':'mediumpurple'}
_ybysbin_marker = {'yb0ys0':"o", 'yb0ys1':"^", 'yb0ys2':"s", 'yb1ys0':"d", 'yb1ys1':"P", 'yb2ys0':"v"}
_ybysbin_label = {'yb0ys0':r'$0 \leq y_b < 1$	$0 \leq y^{\ast} < 1$', 'yb0ys1':r'$0 \leq y_b < 1$	$1 \leq y^{\ast} < 2$', 'yb0ys2':r'$0 \leq y_b < 1$	$2 \leq y^{\ast} < 2.4$', 'yb1ys0':r'$1 \leq y_b < 2$	$0 \leq y^{\ast} < 1$', 'yb1ys1':r'$1 \leq y_b < 2$	$1 \leq y^{\ast} < 2$', 'yb2ys0':r'$2 \leq y_b < 2.4$	$0 \leq y^{\ast} < 1$'}

_ybysbin_xind_up_zjetdata = {'yb0ys0':40, 'yb0ys1':30, 'yb0ys2':20, 'yb1ys0':37, 'yb1ys1':30, 'yb2ys0':26}	#constraints coming from rootfiles (zjet data) --> applied on rootfiles

#----------------------------------------------------------------------------------------------------------------#
#												  function definitions											 # 
#----------------------------------------------------------------------------------------------------------------#
#reuse functions that I wrote for dijetplot_fig7.py
#function to read histogram from rootfile and store values in arrays
def read_rootfile(rootfile, objname):
	#histo = TH1D(gROOT.FindObject(objname)) #is a histogram

	print("[ratios.py]: Reading %s from %s" %(objname, rootfile))
	print("--------------------------------------------------------------------------------")
	histo = rootfile.Get(objname)
	#rootfile.ls()
	nbins = histo.GetNbinsX()
	print "nbins: ", nbins
	entrieslist=[]
	errorlist=[]		#need errorbars for ratio datapoints
	x_axis_list=[]	  #list of bincenters --> will become x-axis (only read once [with the ratio])
	low_edges_list=[]
	up_edges_list=[]
	for j in range(1, nbins+1):
		entrieslist.append(histo.GetBinContent(j))
		print "entry %s is %s" %(j, entrieslist[-1])

		if("ratio" in objname):
			x_axis_list.append(histo.GetBinCenter(j))
			low_edges_list.append(histo.GetXaxis().GetBinLowEdge(j))
			up_edges_list.append(histo.GetXaxis().GetBinUpEdge(j))
			errorlist.append(histo.GetBinError(j))
			
	histentries = np.array(entrieslist)
	ratio_err = np.array(errorlist)
	x_axis = np.array(x_axis_list)
	low_bb = np.array(low_edges_list)
	up_bb = np.array(up_edges_list)

	print "histentries: \n", histentries, "\n"
	print "x_axis: \n", x_axis, "\n"
	print "lower bin bounds: \n", low_bb, "\n"
	print "upper bin bounds: \n", up_bb, "\n"

	return histentries, ratio_err, x_axis, low_bb, up_bb
	


#function for plotting the points for the ratio data/theory with errorbars
#def plot_ratio(bname, ax, ax_sum, x_axis, ratio, binerr, low_bb, up_bb, m_alpha, runperiod, generator): ##low_bb and up_bb are used for drawing the "x-error" == binwidth
def plot_ratio(bname, ax_sum, x_axis, ratio, binerr, low_bb, up_bb, m_alpha, runperiod, generator): ##low_bb and up_bb are used for drawing the "x-error" == binwidth

	x_min= low_bb[0]
	x_max = up_bb[-1]
	print("xmin: %s, xmax: %s" %(x_min, x_max))

	y_min = min(ratio)
	y_max = max(ratio)
	print("ymin: %s, ymax: %s" %(y_min, y_max))

	#for x-error-bars
	xerr_low = np.subtract(x_axis, low_bb)
	xerr_up = np.subtract(up_bb, x_axis)
	#zip(xerr_low, xerr_up) #-->only assigns [(low1, up1), (low2, up2), ...]

	xerrors = np.array([xerr_low, xerr_up])

	#plotting (originally marker size was at 12...)
	ax_sum.errorbar(x_axis, ratio, xerr=xerrors, yerr=binerr, elinewidth=1.2, linewidth=0.0, marker=_ybysbin_marker[bname], ms=8, color=_ybysbin_col[bname], fillstyle='none', fmt='.', label=_ybysbin_label[bname], alpha=m_alpha)


#-----------------------------------------------------------------------------------------------#
#											  main-function										# 
#-----------------------------------------------------------------------------------------------#

def main():
	# start timer
	# measuring wall clock time - not necessary, just for fun.
	start_time = timeit.default_timer()

	# getting the parsed arguments
	parser = argparse.ArgumentParser(epilog='', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	# Positional arguments
	parser.add_argument('rootfile', type=TFile, nargs='?',
						help='Required argument! File with results that shall be plotted. (filename glob)')
	parser.add_argument('datatype', type=str, nargs='?',
						help='Required argument! Specify if file is DATA ("data") or MONTE CARLO ("mc") or THEORY "fnlo".')

	# Optional arguments
	parser.add_argument('--alphamarker', '-a', type=float, default=1.0, nargs='?', action='store',
						help='If set to True: draw marker in summary plot with given input alpha.')

	parser.add_argument('--runperiod', '-r', type=str, nargs='?', default=None, action='store',
						help='Specifies which run period has been used: A, B, C, D')

	parser.add_argument('--generator', '-g', type=str, nargs='?', default=None, action='store',
						help='Specifies which generator has been used: herwig7, pythia8, HT_madgraph_pythia8')

	if len(sys.argv)==1:
		parser.print_help(sys.stderr)
		sys.exit(1)


	# Parse arguments
	args = vars(parser.parse_args())

	#args = args_parsing() # -- do not use this, impractical!

	# input filename
	rootfile = args['rootfile']
	m_alpha = args['alphamarker']
	datatype = args['datatype']
	runperiod = args['runperiod']	#optional argument, but good for labeling --> to be improved
	generator = args['generator']	#optional argument, pythia or herwig
	#rootfile = TFile(args['rootfile'])
	print('\n')
	print("[ratios.py]: Taking input from %s" %rootfile)
	print("[ratios.py]: Looking at %s" %datatype)

	if(generator=="pythia8"):
		generatorname="Pythia8"
	elif(generator=="herwig7"):
		generatorname="Herwig7"


	#prepare the lists (will later be arrays)
	bin_names = np.array(["yb0ys0", "yb0ys1", "yb0ys2", "yb1ys0", "yb1ys1", "yb2ys0"])

	#for the summary (all bins)
	gs = gridspec.GridSpec(3,3)
	fig_tot = plt.figure(figsize=(8,7), dpi=300)
	#ax1 = plt.subplot(gs[:-1,:])
	ax1_sum = plt.subplot(xscale="log", yscale="log")
	#patches = []

	#adjustments that are only done once for the summary plot
	#ax1_sum.axhline(y=1.0, xmin=0, xmax=1, linewidth=0.4, color=_ybysbin_col[bname], linestyle='dashed', dashes=(5,10))

 
	for j in range(0, 6):
		#current bin name:
		bname = bin_names[j]

		#individual plots:
		gs = gridspec.GridSpec(3,3)
   		fig = plt.figure(figsize=(7,7), dpi=300)
   		#ax1 = plt.subplot(gs[:-1,:])
		ax1 = plt.subplot()
 
		#reading the ratio in the current ybys bin from the rootfile
		#read (statistical?) errors on the ratio at the same time
		#furthermore: read x-axis and lower and upper bin bounds
		ratio, err, xax, low_bb, up_bb = read_rootfile(rootfile, "Standard/%s/ratio_%s"%(bname, bname))

		#Do some slicing (lower bins are present in the ratio, but empty due to dijet trigger etc.)
		#CHANGED to more sophisticated slicing (see also plot_xs_ratios.py script etc.)
		#why was this set to [6:62] before?
		ratio = ratio[11:_ybysbin_xind_up_zjetdata[bname]]
		err = err[11:_ybysbin_xind_up_zjetdata[bname]]
		xax = xax[11:_ybysbin_xind_up_zjetdata[bname]]
		low_bb = low_bb[11:_ybysbin_xind_up_zjetdata[bname]]
		up_bb = up_bb[11:_ybysbin_xind_up_zjetdata[bname]]
	
		print "xax: ", xax
		print "Ratio in %s"%bname
		print ratio

		#plot the ratio for this bin
		plot_ratio(bname, ax1_sum, xax, ratio, err, low_bb, up_bb, m_alpha, runperiod, generator)

	#take it from ratios_fnlo.py, where i use my standard scheme
	#---------------------------#
	#	legend, title and save	#
	#---------------------------#
	#ax_ratio.grid(True, which="both", axis="y", color="k", linestyle="dotted", linewidth=0.8, alpha=0.2)
	ax1_sum.grid(True, which="minor", axis="y", color="k", linestyle="dotted", linewidth=0.8, alpha=0.2)
	ax1_sum.grid(True, which="major", axis="y", color="k", linestyle="dotted", linewidth=0.8, alpha=0.6)
	ax1_sum.grid(True, which="both", axis="x", color="k", linestyle="dotted", linewidth=0.8, alpha=0.2)
	ax1_sum.set_xlabel(r'$\mathrm{p_{T,avg}}$ /GeV', fontsize=13, horizontalalignment='right')
	ax1_sum.xaxis.set_label_coords(1.00, -0.06)
	ax1_sum.set_ylabel(r'XS Ratios: Dijet / Zjet', fontsize=13, horizontalalignment='center', verticalalignment='top', y=0.5, rotation=90, labelpad=24)

	#Changed the ticks! --> for old version see ratios_root.py
	x_minticks = [50, 60, 70, 80, 90, 200, 300, 400, 500, 600, 700, 800, 900, 2000]
	x_minticklabels = ['','60', '', '', '', '200', '', '', '500','', '700', '', '', '2000']
	x_majticks = [100, 1000]
	x_majticklabels = ['100', '1000']
	ax1_sum.set_xticks(x_majticks, minor=False)	#set major ticks loc
	ax1_sum.set_xticks(x_minticks, minor=True)		#set minor ticks loc
	ax1_sum.set_xticklabels(x_majticklabels, minor=False, fontsize=12)	#set major ticks label
	ax1_sum.set_xticklabels(x_minticklabels, minor=True, fontsize=12)	#set minor titcks label

	#ytick labels
	for tick in ax1_sum.yaxis.get_major_ticks():
		tick.label.set_fontsize(12)


	#test for suitable y-axis range
	ax1_sum.set_ylim(1e3, 4*1e6)
	
	ax1_sum.legend(fontsize=12, numpoints=1, loc='upper right')
	#ax1_sum.set_title("Cross Section Ratios: Dijet over Z+jet (%s)"%datatype, fontsize=12)

	plt.tight_layout()
	if(datatype=="data"):
		#ax1_sum.set_title("Cross Section Ratios in Data:\n Dijet over Z+jet (Run %s)"%runperiod, fontsize=12)
		#ax1_sum.set_title("Cross Section Ratios: Dijet over Z+jet\n CMS Data (2018 Run %s)"%runperiod, fontsize=14)
		ax1_sum.set_title("Cross Section Ratios: Dijet over Z+jet\n CMS Data (2018 Run %s, 31.93 fb$^{-1}$)"%runperiod, fontsize=14) #version with lumi
		fig_tot.savefig("ratio_summary_data_Run%s.png"%runperiod)
	elif(datatype=="mc"):
		#ax1_sum.set_title("Cross Section Ratios in Monte Carlo:\n Dijet over Z+jet (%s)"%generator, fontsize=12)
		ax1_sum.set_title("Cross Section Ratios: Dijet over Z+jet\n Monte Carlo Simulation (%s without/with MadGraph)"%generatorname, fontsize=14)
		fig_tot.savefig("ratio_summary_mc_%s.png"%generator)

	#now zoom in:
	'''
	ax1_sum.set_ylim(1e3, 1e7)	#some points won't be displayed
	ax1_sum.text(0.02, 0.04, "zoomed in y-axis", fontsize=12, fontstyle='italic', ha='left', va='bottom', transform=ax1_sum.transAxes)
	plt.tight_layout()
	if(datatype=="data"):
		fig_tot.savefig("ratio_summary_data_Run%s_zoomed.png"%runperiod)
	elif(datatype=="mc"):
		fig_tot.savefig("ratio_summary_mc_%s_zoomed.png"%generator)
	'''

	#test with even smaller x-range
	ax1_sum.set_xlim(50,1200)
	plt.tight_layout()
	if(datatype=="data"):
		fig_tot.savefig("ratio_summary_data_Run%s_smallxrange_50to1200.png"%runperiod)
	elif(datatype=="mc"):
		fig_tot.savefig("ratio_summary_mc_%s_smallxrange_50to1200.png"%generator)
	
	


	# stop timer
	stop_time = timeit.default_timer()
	timediff = stop_time-start_time
	print("[ratios_root.py]: Elapsed time: %s sec = %s min" %(timediff, round(timediff/60., 2)))



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
