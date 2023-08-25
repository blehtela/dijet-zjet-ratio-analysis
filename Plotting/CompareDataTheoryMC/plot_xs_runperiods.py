#!/usr/bin/env python2
#-*- coding:utf-8 -*-


#####################################################
#													#
#	XS plot of all runperiods:						#
#	coparing data A, B, C, D.						#
#	Taking input from the ROOT file					#
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
def read_rootfile(rootfile, objname):
	#histo = TH1D(gROOT.FindObject(objname)) #is a histogram

	print("[plots_data_mc.py]: Reading %s from %s" %(objname, rootfile))
	print("--------------------------------------------------------------------------------")
	histo = rootfile.Get(objname)
	#rootfile.ls()
	nbins = histo.GetNbinsX()
	print "nbins: ", nbins
	entrieslist=[]
	errorlist=[]		#need errorbars for datapoints
	x_axis_list=[]	  	#list of bincenters --> will become x-axis (only read once [with the ratio])
	low_edges_list=[]
	up_edges_list=[]
	for j in range(1, nbins+1):
		entrieslist.append(histo.GetBinContent(j))
		print "entry %s is %s" %(j, entrieslist[-1])

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

	print "histentries: \n", histentries, "\n"
	print "x_axis: \n", x_axis, "\n"
	print "lower bin bounds: \n", low_bb, "\n"
	print "upper bin bounds: \n", up_bb, "\n"

	return histentries, entries_err, x_axis, low_bb, up_bb
	

#function to plot a xs in the xs plot
#for individual ybys bins
def plot_xs(bname, ax, x_axis, xs_values, xs_err, m_color, m_alpha, process, labeling, legbool=False):
##low_bb and up_bb are used for drawing the "x-error" == binwidth

	#plot xs as errorbar plot
	#originally markersize was set to 10.
	ax.errorbar(x_axis, xs_values, xerr=0.0, yerr=xs_err, elinewidth=1.2, linewidth=1.0, marker=_ybysbin_marker[bname], ms=6, color=m_color, fillstyle='none', fmt='.', label=labeling)
	ax.grid(True, which="major", axis="y", color="k", linestyle="dotted", linewidth=1)
	#plt.grid(True, which="major", linestyle="dotted", color='k')

	#ax.set_xlim(40, 1200)
	#ax.set_ylim(y_min+0.0000001, y_max+0.1)

	ax.set_xlabel(r'$\mathrm{p_{T,avg}}$ /GeV', fontsize=12, horizontalalignment='right')
	ax.xaxis.set_label_coords(1.00, -0.06)
	#plt.xlabel(r'$\mathrm{p_{T,avg}}$ /GeV', horizontalalignment='right', loc='right')
	ax.set_ylabel(r'XS in Data in fb/GeV', fontsize=12, horizontalalignment='center', verticalalignment='top', y=0.5, rotation=90, labelpad=24)


#function to plot the data/mc ratio for either dijet or zjet
#for individual ybys bins
def plot_ratio(bname, ax, x_axis, ratio, ratio_err, m_color, m_alpha, process, labeling):
##low_bb and up_bb are used for drawing the "x-error" == binwidth

	#plot the ratio as errorbar
	ax.errorbar(x_axis, ratio, xerr=0.0, yerr=ratio_err, elinewidth=1.2, linewidth=1.0, marker=_ybysbin_marker[bname], ms=6, color=m_color, fillstyle='none', fmt='.', label=labeling)

	#plt.yscale("log")
	plt.xscale("log")
	ax.grid(True, which="both", axis="y", color="k", linestyle="dotted", linewidth=1, alpha=0.8)
	plt.grid(True, which="major", linestyle="dotted", color='k', alpha=0.7)

	#ax.set_xlim(40, 1200)
	#ax.set_ylim(y_min+0.0000001, y_max+0.1)

	ax.set_xlabel(r'$\mathrm{p_{T,avg}}$ /GeV', fontsize=12, horizontalalignment='right')
	ax.xaxis.set_label_coords(1.00, -0.06)
	ax.set_ylabel(r'Data XS Ratios', fontsize=12, horizontalalignment='center', verticalalignment='top', y=0.5, rotation=90, labelpad=24)

	ax.set_title(r'%s in %s (%s)'%(process,bname, _ybysbin_label[bname]))
	ax.legend(fontsize=10, numpoints=1, loc='upper left')




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
	parser.add_argument('data_rootfiles', type=TFile, nargs=4,
						help='Required argument! Files with XS for Run A, B, C, D (in this order!) that shall be plotted.')
	parser.add_argument('process', type=str, nargs='?', choices=["zjet", "dijet"],
						help="Required argument! Specifies which analysis is considered: zjet for Z+jet analysis; dijet for Dijet analysis.") 

	# Optional arguments
	parser.add_argument('--alphamarker', '-a', type=float, default=1.0, nargs='?', action='store',
						help='If set to True: draw marker in summary plot with given input alpha.')
	parser.add_argument('--ybysbin', '-b', type=str, nargs='?', default=None,
						help="Specifies yboost ystar bin. If None, plot all 6 bins.")
	parser.add_argument('--xsplot', '-x', default=False, action='store_true', 
						help="When chosen, plot also xs-overview.")




	if len(sys.argv)==1:
		parser.print_help(sys.stderr)
		sys.exit(1)

	# Parse arguments
	args = vars(parser.parse_args())

	# input filename
	data_rootfiles = args['data_rootfiles']
	process = args['process']		#specifies whether looking at dijet or zjet
	ybys = args['ybysbin']			#name of ybys bin
	xsplot = args['xsplot']			#bool whether xs shall be plotted
	m_alpha = args['alphamarker']

	print('\n')
	print("[plot_xs_runperiods.py]: Taking data input from \n %s" %data_rootfiles)


	#prepare the lists (will later be arrays)
	#will loop through the ybys bins if bname==None
	if(ybys==None):	#no ybys bin specified --> will loop through all 6
		bin_names = np.array(["yb0ys0", "yb0ys1", "yb0ys2", "yb1ys0", "yb1ys1", "yb2ys0"])
	else:
		bin_names = np.array([ybys])

	#for(int i=0; i!=bin_names.size; ++i):
	for i in range(6):	#loops from 0 to 5
		bname = bin_names[i]
		#preparing the plot for current ybys bin
		#for the summary (all runperiods) of the XS plots
		gs_xs = gridspec.GridSpec(3,3)
		fig_xs = plt.figure(figsize=(7,7), dpi=300)
		ax_xs = plt.subplot(xscale="log", yscale="log")

		#the ratio axes
		gs_ratio = gridspec.GridSpec(3,3)
		fig_ratio = plt.figure(figsize=(7,7), dpi=300)
		ax_ratio = plt.subplot(xscale="log", yscale="linear")



		#reading the data points in the current ybys bin from the rootfile
		#read (statistical?) errors on the histogram entries at the same time
		#furthermore: read x-axis and lower and upper bin bounds
		dataA, dataA_err, dataA_xax, dataA_low_bb, dataA_up_bb = read_rootfile(data_rootfiles[0], "CrossSections_byBinwidth/%s/sumhisto_%s_RunA_%s"%(bname, process, bname))
		dataB, dataB_err, dataB_xax, dataB_low_bb, dataB_up_bb = read_rootfile(data_rootfiles[1], "CrossSections_byBinwidth/%s/sumhisto_%s_RunB_%s"%(bname, process, bname))
		dataC, dataC_err, dataC_xax, dataC_low_bb, dataC_up_bb = read_rootfile(data_rootfiles[2], "CrossSections_byBinwidth/%s/sumhisto_%s_RunC_%s"%(bname, process, bname))
		dataD, dataD_err, dataD_xax, dataD_low_bb, dataD_up_bb = read_rootfile(data_rootfiles[3], "CrossSections_byBinwidth/%s/sumhisto_%s_RunD_%s"%(bname, process, bname))


		print("data A: \n %s"%dataA)
		print("data B: \n %s"%dataB)
		print("data C: \n %s"%dataC)
		print("data D: \n %s"%dataD)


		#-----------------------#
		# calculate the ratios	#
		#-----------------------#
		ratio_AtoD = np.divide(dataA, dataD)
		ratio_BtoD = np.divide(dataB, dataD)
		ratio_CtoD = np.divide(dataC, dataD)
		ratio_DtoD = np.divide(dataD, dataD)

		#calculate the relative xs errors:
		dataA_relerr = np.divide(dataA_err, dataA)
		dataB_relerr = np.divide(dataB_err, dataB)
		dataC_relerr = np.divide(dataC_err, dataC)
		dataD_relerr = np.divide(dataD_err, dataD)

		#calculate the ratio relative error
		ratioA_relerr = np.sqrt(pow(dataA_relerr,2)+pow(dataD_relerr,2))	#np.sqrt()-->element wise sqrt
		ratioB_relerr = np.sqrt(pow(dataB_relerr,2)+pow(dataD_relerr,2))
		ratioC_relerr = np.sqrt(pow(dataC_relerr,2)+pow(dataD_relerr,2))
		ratioD_relerr = np.sqrt(pow(dataD_relerr,2)+pow(dataD_relerr,2))

		#calculate the absolute error:
		ratioA_err = np.multiply(ratio_AtoD, ratioA_relerr)
		ratioB_err = np.multiply(ratio_BtoD, ratioB_relerr)
		ratioC_err = np.multiply(ratio_CtoD, ratioC_relerr)

		ratioD_err = np.full(dataD.size, 0)	#set to zero as this is just the normalisation
		#ratioD_err = np.zeros


		#-----------------------#
		#		plotting		#
		#-----------------------#
		#creating the plots for the chosen ybys bin

		#plotting the three cross sections (take always same x-axis, should not differ)
		if(xsplot==True):
			plot_xs(bname, ax_xs, dataA_xax, dataA, dataA_err, _ybysbin_col[bname], m_alpha, process, "Data Run A in %s"%(_ybysbin_label[bname]), True)	#runperiod A
			plot_xs(bname, ax_xs, dataA_xax, dataB, dataB_err, _ybysbin_col[bname], m_alpha, process, "Data Run B in %s"%(_ybysbin_label[bname]), True)	#runperiod B
			plot_xs(bname, ax_xs, dataA_xax, dataC, dataC_err, _ybysbin_col[bname], m_alpha, process, "Data Run C in %s"%(_ybysbin_label[bname]), True)	#runperiod C
			plot_xs(bname, ax_xs, dataA_xax, dataD, dataD_err, _ybysbin_col[bname], m_alpha, process, "Data Run D in %s"%(_ybysbin_label[bname]), True)	#runperiod D

			ax_xs.set_title("%s XS runperiod overview (%s)"%(process, bin_names[i]))
			ystart,  yend = ax_xs.get_ylim()
			#ax_xs.set_ylim((ystart, yend))
			#ax_xs.set_yticks(np.arange(math.floor(ystart), math.ceil(yend), step=0.1))
			#ax_xs.set_ylim((0.0, 2.0))
			#ax_xs.set_yticks(np.arange(0.0, 2.0, step=0.1))
			#xstart, xend = ax_xs.get_xlim()	#get the xrange of the xs plot, set the same for the ratio plot
			#ax_xs.set_xlim((xstart, xend))

			#save plots for current ybys bin
			#xs plot:
			#plt.legend()
			ax_xs.legend(fontsize=14, numpoints=1, loc='lower left')
			plt.tight_layout()
			fig_xs.savefig("%s_xs_runperiods_%s.png"%(process, bname))



		#plot the ratios
		plot_ratio(bname, ax_ratio, dataA_xax, ratio_AtoD, ratioA_err, 'gold', m_alpha, process, "RunA/RunD")			#A over D
		plot_ratio(bname, ax_ratio, dataA_xax, ratio_BtoD, ratioB_err, 'orangered', m_alpha, process, "RunB/RunD")		#B over D
		plot_ratio(bname, ax_ratio, dataA_xax, ratio_CtoD, ratioC_err, 'darkturquoise', m_alpha, process, "RunC/RunD")	#C over D
		plot_ratio(bname, ax_ratio, dataA_xax, ratio_DtoD, ratioD_err, 'yellowgreen', m_alpha, process, "RunD/RunD")	#D over D

		ax_ratio.set_ylim((0.0,2.0))
		ax_ratio.axhline(1.0, linestyle='--', color='k')



		#save the figure: 
		#-------------------#
		# legend and save 	#
		#-------------------#
		#save plots for current ybys bin
		ax_ratio.legend(fontsize=12, numpoints=1, loc='upper left')
		plt.tight_layout()
		fig_ratio.savefig("%s_xs_runperiods_ratios_%s.png"%(process, bname))




	# stop timer
	stop_time = timeit.default_timer()
	timediff = stop_time-start_time
	print("[plot_xs_runperiods.py]: Elapsed time: %s sec = %s min" %(timediff, round(timediff/60., 2)))



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




