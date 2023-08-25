#!/usr/bin/env python2
#-*- coding:utf-8 -*-


#####################################################
#													#
#	XS plots and ratio:								#
#	coparing data, theory, MC						#
#	Taking input from the ROOT file					#
#	and from fastNLO tables							#
#													#
#	Created by B. Schillinger, 09.02.2020			#
#	Last modified: B.Schillinger, 09.02.2020		#
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
	

#function to read XS from fastnlo table
def read_fnlotable(fnlotable, pdfset, order):
	fnlo = fastnlo.fastNLOLHAPDF(fnlotable, pdfset)
	if(order==0):	#only LO desired
		fnlo.SetContributionON(fastnlo.kFixedOrder, 0, True)
		fnlo.SetContributionON(fastnlo.kFixedOrder, 1, False)
	elif(order==1):	#switch LO and NLO on
		fnlo.SetContributionON(fastnlo.kFixedOrder, 0, True)
		fnlo.SetContributionON(fastnlo.kFixedOrder, 1, True)
	
	fnlo.CalcCrossSection()
	xs_ord = fnlo.GetCrossSection()
	xs_array = np.array(xs_ord)

	#should also calculate uncertainties --> PDF unc and Scale unc!!!
	#return xs_array, entries_err, x_axis	#x_axis should be the same as for the root hists
	return xs_array



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
	ax.set_ylabel(r'XS Data, Theory and MC in fb/GeV', fontsize=12, horizontalalignment='center', verticalalignment='top', y=0.5, rotation=90, labelpad=24)


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
	ax.set_ylabel(r'Ratios', fontsize=12, horizontalalignment='center', verticalalignment='top', y=0.5, rotation=90, labelpad=24)

	ax.set_title(r'%s in %s (%s)'%(process,bname, _ybysbin_label[bname]))
	ax.legend(fontsize=10, numpoints=1, loc='upper left')


	#add it to the summary plot
	#ax_ratio_sum.errorbar(x_axis, ratio, xerr=0.0, yerr=ratio_err, elinewidth=1.2, linewidth=1.0, marker=_ybysbin_marker[bname], ms=6, color=_ybysbin_col[bname], fillstyle='none', fmt='.', label=_ybysbin_label[bname], alpha=m_alpha)
	#does the following work?
	#ax_ratio_sum.grid(True, which="both", axis="y", color="k", linestyle="dotted", linewidth=1)




#function to plot XS for data and mc in upper part of plot
#and ratio data/mc in lower part of plot
#for individual ybys bins
def plot_combi(bname, ax_xs, ax_ratio, x_axis, data, data_err, ratio, ratio_err, m_color, m_alpha, process, labeling):

	#call the XS plotting function for the upper part of the plot
	plot_xs(bname, ax_xs, x_axis, data, data_err, m_color, m_alpha, process, labeling, False)
	ax_xs.set_xlabel("")
	ax_xs.set_yscale("log")

	#call the ratio plotting function for the lower part of the plot
	plot_ratio(bname, ax_ratio, x_axis, ratio, ratio_err, m_color, m_alpha, process, labeling)
	ax_ratio.set_title("")
	ax_ratio.set_yscale("linear")
	#ax_ratio.set_ylim(-1.5, 2.0)	#test
	ax_ratio.set_ylim(ymin=0.0, ymax=2.0)	#test
	#ax_ratio.get_legend().remove()

	xstart, xend = ax_xs.get_xlim()	#get the xrange of the xs plot, set the same for the ratio plot
	ax_ratio.set_xlim((xstart, xend))



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
	parser.add_argument('mc_rootfile', type=TFile, nargs='?',
						help='Required argument! File with MC results that shall be plotted. (filename glob)')
	parser.add_argument('data_rootfile', type=TFile, nargs='?',
						help='Required argument! File with data results that shall be plotted. (filename glob)')
	parser.add_argument('fnlo_table', nargs='?',
						help='Required argument! FastNLO table to be evaluated.')
	parser.add_argument('process', type=str, nargs='?', choices=["zjet", "dijet"],
						help="Required argument! Specifies which analysis is considered: zjet for Z+jet analysis; dijet for Dijet analysis.") 
	parser.add_argument('ybysbin', type=str, nargs='?',
						help="Required argument! Specifies yboost ystar bin name.")

	# Optional arguments
	parser.add_argument('--alphamarker', '-a', type=float, default=1.0, nargs='?', action='store',
						help='If set to True: draw marker in summary plot with given input alpha.')

	parser.add_argument('--pdfset', '-p', default='CT14nlo', nargs='?',
						help='PDF set for fnlo table evaluation. Default is CT14nlo.')

	parser.add_argument('--runperiod', '-r', type=str, nargs='?', default=None, action='store',
						help='Specifies which run period has been used: A, B, C, D')
	parser.add_argument('--generator', '-g', type=str, nargs='?', default=None, action='store',
						help='Specifies which generator has been used: herwig7, pythia8')


	#parser.add_argument('--alphamarker', '-a', default=False, action='store_true'
	#					help='If set to True: draw marker in summary plot with alpha=0.9.')

	if len(sys.argv)==1:
		parser.print_help(sys.stderr)
		sys.exit(1)

	# Parse arguments
	args = vars(parser.parse_args())

	# input filename
	mc_rootfile = args['mc_rootfile']
	data_rootfile = args['data_rootfile']
	fnlo_table = args['fnlo_table']
	process = args['process']		#specifies whether looking at dijet or zjet
	bname = args['ybysbin']			#name of ybys bin
	m_alpha = args['alphamarker']
	pdfset = args['pdfset']
	runperiod = args['runperiod']	#optional argument, but good for labeling --> to be improved
	generator = args['generator']	#optional argument, pythia or herwig
	#rootfile = TFile(args['rootfile'])
	print('\n')
	print("[plots_data_theory_mc.py]: Taking MC input from %s" %mc_rootfile)
	print("[plots_data_theory_mc.py]: Taking data input from %s" %data_rootfile)
	print("[plots_data_theory_mc.py]: Taking theory input from %s" %fnlo_table)


	#prepare the lists (will later be arrays)
	bin_names = np.array(["yb0ys0", "yb0ys1", "yb0ys2", "yb1ys0", "yb1ys1", "yb2ys0"])


	#---------------------------------------------------#
	#	preparing the three plots for chosen ybys bin	#
	#---------------------------------------------------#
	#for the summary (all bins) of the XS plots
	gs_xs = gridspec.GridSpec(3,3)
	fig_xs = plt.figure(figsize=(7,7), dpi=300)
	#ax1 = plt.subplot(gs[:-1,:])
	ax_xs = plt.subplot(xscale="log", yscale="log")
	#patches = []

	#for the summary (all bins) of the data/mc ratio plots
	gs_ratio = gridspec.GridSpec(3,3)
	fig_ratio = plt.figure(figsize=(7,7), dpi=300)
	ax_ratio = plt.subplot(xscale="log",yscale="linear")

	#for the summary of the combined XS and ratio plots
	gs_combi = gridspec.GridSpec(3,3)
	fig_combi = plt.figure(figsize=(7,7), dpi=300)
	ax_combi_xs = plt.subplot(gs_combi[:-1,:], xscale="log", yscale="log")			#upper part of plot
	ax_combi_ratio = plt.subplot(gs_combi[-1:,:], xscale="log", yscale="linear")	#lower part of plot

	ax_combi_ratio.grid(True, which="both", axis="y", color="k", linestyle="dotted", linewidth=1)

	
	#reading the data and mc points in the current ybys bin from the rootfile
	#read (statistical?) errors on the histogram entries at the same time
	#furthermore: read x-axis and lower and upper bin bounds
	if(process=="dijet"):
		data, data_err, data_xax, data_low_bb, data_up_bb = read_rootfile(data_rootfile, "CrossSections_byBinwidth/%s/sumhisto_dijet_Run%s_%s"%(bname, runperiod, bname))
		mc, mc_err, mc_xax, mc_low_bb, mc_up_bb = read_rootfile(mc_rootfile, "CrossSections_byBinwidth/%s/sumhisto_dijet_%s_%s"%(bname, generator, bname))
		fnloxs_LO = read_fnlotable(fnlo_table, pdfset, 0)
		fnloxs_NLO = read_fnlotable(fnlo_table, pdfset, 1)
		print("data: \n %s"%data)

	elif(process=="zjet"):
		data, data_err, data_xax, data_low_bb, data_up_bb = read_rootfile(data_rootfile, "CrossSections_byBinwidth/%s/sumhisto_zjet_Run%s_%s"%(bname, runperiod, bname))
		mc, mc_err, mc_xax, mc_low_bb, mc_up_bb = read_rootfile(mc_rootfile, "CrossSections_byBinwidth/%s/sumhisto_zjet_%s_%s"%(bname, generator, bname))
		#fnloxs = read_fnlotable(fnlo_table)
		fnloxs_LO = read_fnlotable(fnlo_table, pdfset, 0)
		fnloxs_NLO = read_fnlotable(fnlo_table, pdfset, 1)


	#as long as we do not calculate pdf and scale uncs
	nentries = fnloxs_LO.size
	fnloLO_err = np.full(nentries,0)
	fnloNLO_err = np.full(nentries,0)

	#convert from pb to fb (multiply by 1000)
	fnloxs_LO = 1000*fnloxs_LO
	fnloxs_NLO = 1000*fnloxs_NLO


	#remove the lower and higher bins that are not present in fastnlo tables
	#get upper limit from dictionary _ybysbin_xind_up
	#data_xax = data_xax[7:62]	#bin 7 is included, bin 62 excluded
	data_xax = data_xax[7:_ybysbin_xind_up[bname]]
	data = data[7:_ybysbin_xind_up[bname]]
	data_err = data_err[7:_ybysbin_xind_up[bname]]
	mc_err = mc_err[7:_ybysbin_xind_up[bname]]
	mc = mc[7:_ybysbin_xind_up[bname]]

	#CALCULATE THE RATIO DATA/MC! ARRAYS! --> how to handle the errors?
	print("data again: \n %s"%data)
	ratio_datamc = np.divide(data, mc)
	ratio_datafnloNLO = np.divide(data, fnloxs_NLO)
	ratio_datafnloLO = np.divide(data, fnloxs_LO)

	#try simple error propagation for absolute error of ratio:
	# ratio_err = Ratio*Sqrt((data_err/data)^2 + (mc_err/mc)^2)
	# calculate relative errors:
	data_relerr = np.divide(data_err, data)	# relative error
	mc_relerr = np.divide(mc_err, mc)		# relative error

	#calculate ratio relative error:
	ratio_relerr = np.sqrt(pow(data_relerr,2)+pow(mc_relerr,2))	#numpy.sqrt() returns sqrt for array elementwise

	#calulate ratio absolute error:
	ratio_err = np.multiply(ratio_datamc, ratio_relerr)	#elementwise numpy multiplication of arrays



	#-----------------------#
	#		plotting		#
	#-----------------------#
	#creating the plots for the chosen ybys bin

	#plotting the three cross sections
	plot_xs(bname, ax_xs, data_xax, data, data_err, _ybysbin_col[bname], m_alpha, process, "Data Run%s in %s"%(runperiod, _ybysbin_label[bname]), True)#data
	plot_xs(bname, ax_xs, data_xax, mc, mc_err, 'saddlebrown', m_alpha, process, "MC %s in %s"%(generator, bname), True)#mc
	plot_xs(bname, ax_xs, data_xax, fnloxs_LO, fnloLO_err, 'gray', m_alpha, process, "LO Table %s"%(pdfset), True)#fnlo LO
	plot_xs(bname, ax_xs, data_xax, fnloxs_NLO, fnloLO_err, 'black', m_alpha, process, "NLO Table %s"%(pdfset), True)#fnlo LO

	ax_xs.set_title("%s XS in %s (Data Run%s, MC %s, FixedOrder Theory with %s)"%(process, bname, runperiod, generator, pdfset))


	#plotting the ratios
	plot_ratio(bname, ax_ratio, data_xax, ratio_datamc, ratio_err, _ybysbin_col[bname], m_alpha, process, "data (Run%s) / MC (%s)"%(runperiod, generator))#data/mc
	plot_ratio(bname, ax_ratio, data_xax, ratio_datafnloLO, fnloLO_err, 'gray', m_alpha, process, "data (Run%s) / LO theo."%runperiod)#data/LO
	plot_ratio(bname, ax_ratio, data_xax, ratio_datafnloNLO, fnloLO_err, 'black', m_alpha, process, "data (Run%s) / NLO theo."%runperiod)#data/NLO

	ax_ratio.set_title("%s XS ratio in %s (Data Run%s, MC %s, FixedOrder Theory with %s)"%(process, bname, runperiod, generator, pdfset))

	ystart,  yend = ax_ratio.get_ylim()
	#ax_ratio.set_ylim((ystart, yend))
	#ax_ratio.set_yticks(np.arange(math.floor(ystart), math.ceil(yend), step=0.1))
	ax_ratio.set_ylim((0.0, 2.0))
	ax_ratio.set_yticks(np.arange(0.0, 2.0, step=0.1))
	#xstart, xend = ax_xs.get_xlim()	#get the xrange of the xs plot, set the same for the ratio plot
	#ax_ratio.set_xlim((xstart, xend))



	print "data/LO"
	print ratio_datafnloLO
	print "data/NLO"
	print ratio_datafnloNLO


	#plot the combi plots
#def plot_combi(bname, ax_xs, ax_ratio, x_axis, data, data_err, ratio, ratio_err, col, m_alpha, process, labeling):
	#plot_combi(bname, ax_combi_xs, ax_combi_ratio, data_xax, data_

	#do the same for the combination axes... should write better function for this
	plot_xs(bname, ax_combi_xs, data_xax, data, data_err, _ybysbin_col[bname], m_alpha, process, "Data Run%s in %s"%(runperiod, _ybysbin_label[bname]), True)#data
	plot_xs(bname, ax_combi_xs, data_xax, mc, mc_err, 'saddlebrown', m_alpha, process, "MC %s in %s"%(generator, bname), True)#mc
	plot_xs(bname, ax_combi_xs, data_xax, fnloxs_LO, fnloLO_err, 'gray', m_alpha, process, "LO Table %s"%(pdfset), True)#fnlo LO
	plot_xs(bname, ax_combi_xs, data_xax, fnloxs_NLO, fnloLO_err, 'black', m_alpha, process, "NLO Table %s"%(pdfset), True)#fnlo LO

	plot_ratio(bname, ax_combi_ratio, data_xax, ratio_datamc, ratio_err, _ybysbin_col[bname], m_alpha, process, "data (Run%s) / MC (%s)"%(runperiod, generator))#data/mc
	plot_ratio(bname, ax_combi_ratio, data_xax, ratio_datafnloLO, fnloLO_err, 'gray', m_alpha, process, "data (Run%s) / LO theo."%runperiod)#data/LO
	plot_ratio(bname, ax_combi_ratio, data_xax, ratio_datafnloNLO, fnloLO_err, 'black', m_alpha, process, "data (Run%s) / NLO theo."%runperiod)#data/NLO

	ax_combi_xs.set_title("%s in %s (Data Run%s, MC %s, FixedOrder Theory with %s)"%(process, bname, runperiod, generator, pdfset))
	ax_combi_ratio.set_ylim((0.0, 2.0))


	#save the figures: should be three in total
	#-------------------#
	# legend and save 	#
	#-------------------#
	#save plots for current ybys bin
	#xs plot:
	#plt.legend()
	ax_xs.legend(fontsize=14, numpoints=1, loc='lower left')
	ax_combi_xs.legend(fontsize=12, numpoints=1, loc='lower left')
	ax_combi_ratio.legend(fontsize=9, numpoints=1, loc='upper left')
	plt.tight_layout()
	fig_xs.savefig("%s_xs_%s_theory_Run%s_%s.png"%(process, bname, runperiod, generator))
	fig_ratio.savefig("%s_ratio_%s_theory_Run%s_%s.png"%(process, bname, runperiod, generator))
	fig_combi.savefig("%s_combi_%s_theory_Run%s_%s.png"%(process, bname, runperiod, generator))


	# stop timer
	stop_time = timeit.default_timer()
	timediff = stop_time-start_time
	print("[plots_data_theory_mc.py]: Elapsed time: %s sec = %s min" %(timediff, round(timediff/60., 2)))



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


