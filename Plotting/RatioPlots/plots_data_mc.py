#!/usr/bin/env python2
#-*- coding:utf-8 -*-


#####################################################
#													#
#	Various plots of data and MC					#
#	Taking input from the ratio ROOT file			#
#													#
#	Created by B. Schillinger, 29.11.2019			#
#													#
#####################################################

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
	

#function to plot data and mc for either dijet or zjet into same plot
#for individual ybys bins
#def plot_xs(bname, ax, ax_xs_sum, x_axis, data, data_err, data_low_bb, data_up_bb, mc, mc_err, m_alpha, runperiod, generator):
def plot_xs(bname, ax, ax_xs_sum, x_axis, data, data_err, mc, mc_err, m_alpha, process, runperiod, generator):
##low_bb and up_bb are used for drawing the "x-error" == binwidth

	'''
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
	#print "xerr_low: \n", xerr_low
	#print "xerr_up: \n", xerr_up

	xerrors = np.array([xerr_low, xerr_up])
	'''


	#plot data and mc as errorbar
	ax.errorbar(x_axis, data, xerr=0.0, yerr=data_err, elinewidth=1.2, linewidth=1.0, marker=_ybysbin_marker[bname], ms=10, color=_ybysbin_col[bname], fillstyle='none', fmt='.', label=_ybysbin_label[bname])
	ax.errorbar(x_axis, mc, xerr=0.0, yerr=mc_err, elinewidth=1.2, linewidth=1.0, marker=_ybysbin_marker[bname], ms=10, color='k', fillstyle='none', fmt='.', label=_ybysbin_label[bname])
	plt.yscale("log")
	plt.xscale("log")

	plt.grid(True, which="major", linestyle="dotted", color='k')


	#ax.set_xlim(40, 1200)
	#ax.set_ylim(y_min+0.0000001, y_max+0.1)

	ax.set_xlabel(r'$\mathrm{p_{T,avg}}$ /GeV', fontsize=12, horizontalalignment='right')
	ax.xaxis.set_label_coords(1.00, -0.06)
	#plt.xlabel(r'$\mathrm{p_{T,avg}}$ /GeV', horizontalalignment='right', loc='right')
	ax.set_ylabel(r'XS Data and MC', fontsize=12, horizontalalignment='center', verticalalignment='top', y=0.5, rotation=90, labelpad=24)

	if((runperiod!=None)and(generator!=None)):
		ax.set_title(r'%s in %s (Run %s, %s)'%(process,bname,runperiod,generator))
	else:
		ax.set_title(r'%s in %s'%(process,bname))

	ax.legend(fontsize=14, numpoints=1, loc='upper left')

	'''
	if(runperiod!=None):
		ax.text(0.02, 1.01, "2018 data Run %s"%runperiod, fontsize=13, fontstyle='italic', ha='left', va='bottom', transform=ax.transAxes)          #text: Run period
	else:
		ax.text(0.02, 1.01, "2018 monte carlo", fontsize=13, fontstyle='italic', ha='left', va='bottom', transform=ax.transAxes)          			#text: monte carlo
	'''

	#add it to the summary plot
	ax_xs_sum.errorbar(x_axis, data, xerr=0.0, yerr=data_err, elinewidth=1.2, linewidth=1.0, marker=_ybysbin_marker[bname], ms=10, color=_ybysbin_col[bname], fillstyle='none', fmt='.', label=_ybysbin_label[bname], alpha=m_alpha)
	ax_xs_sum.errorbar(x_axis, mc, xerr=0.0, yerr=mc_err, elinewidth=1.2, linewidth=1.0, marker=_ybysbin_marker[bname], ms=10, color='k', fillstyle='none', fmt='.', label=_ybysbin_label[bname], alpha=m_alpha)




#function to plot the data/mc ratio for either dijet or zjet
#for individual ybys bins
def plot_ratio(bname, ax, ax_ratio_sum, x_axis, ratio, ratio_err, m_alpha, process, runperiod, generator):
##low_bb and up_bb are used for drawing the "x-error" == binwidth

	#plot the ratio as errorbar
	ax.errorbar(x_axis, data, xerr=0.0, yerr=ratio_err, elinewidth=1.2, linewidth=1.0, marker=_ybysbin_marker[bname], ms=10, color=_ybysbin_col[bname], fillstyle='none', fmt='.', label=_ybysbin_label[bname])
	plt.yscale("log")
	plt.xscale("log")

	plt.grid(True, which="major", linestyle="dotted", color='k')

	#ax.set_xlim(40, 1200)
	#ax.set_ylim(y_min+0.0000001, y_max+0.1)

	ax.set_xlabel(r'$\mathrm{p_{T,avg}}$ /GeV', fontsize=12, horizontalalignment='right')
	ax.xaxis.set_label_coords(1.00, -0.06)
	ax.set_ylabel(r'Ratio Data/MC', fontsize=12, horizontalalignment='center', verticalalignment='top', y=0.5, rotation=90, labelpad=24)

	if((runperiod!=None)and(generator!=None)):
		ax.set_title(r'%s in %s (Run %s, %s)'%(process,bname,runperiod,generator))
	else:
		ax.set_title(r'%s in %s'%(process,bname))

	ax.legend(fontsize=14, numpoints=1, loc='upper left')


	#add it to the summary plot
	ax_ratio_sum.errorbar(x_axis, data, xerr=0.0, yerr=ratio_err, elinewidth=1.2, linewidth=1.0, marker=_ybysbin_marker[bname], ms=10, color=_ybysbin_col[bname], fillstyle='none', fmt='.', label=_ybysbin_label[bname], alpha=m_alpha)



#function to plot XS for data and mc in upper part of plot
#and ratio data/mc in lower part of plot
#for individual ybys bins
def plot_combi():
	return 0






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
	parser.add_argument('process', type=str, nargs='?', choices=["zjet", "dijet"],
						help="Required argument! Specifies which analysis is considered: zjet for Z+jet analysis; dijet for Dijet analysis.") 

	# Optional arguments
	parser.add_argument('--alphamarker', '-a', type=float, default=1.0, nargs='?', action='store',
						help='If set to True: draw marker in summary plot with given input alpha.')

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
	process = args['process']		#specifies whether looking at dijet or zjet
	m_alpha = args['alphamarker']
	runperiod = args['runperiod']	#optional argument, but good for labeling --> to be improved
	generator = args['generator']	#optional argument, pythia or herwig
	#rootfile = TFile(args['rootfile'])
	print('\n')
	print("[plots_data_mc.py]: Taking MC input from %s" %mc_rootfile)
	print("[plots_data_mc.py]: Taking data input from %s" %data_rootfile)


	#prepare the lists (will later be arrays)
	bin_names = np.array(["yb0ys0", "yb0ys1", "yb0ys2", "yb1ys0", "yb1ys1", "yb2ys0"])


	#---------------------------------------#
	#	preparing the three summary plots	#
	#---------------------------------------#
	#for the summary (all bins) of the XS plots
	gs_xs = gridspec.GridSpec(3,3)
	fig_xs = plt.figure(figsize=(7,7), dpi=300)
	#ax1 = plt.subplot(gs[:-1,:])
	ax_xs_sum = plt.subplot(xscale="log", yscale="log")
	#patches = []

	#for the summary (all bins) of the data/mc ratio plots
	gs_ratio = gridspec.GridSpec(3,3)
	fig_ratio = plt.figure(figsize=(7,7), dpi=300)
	ax_ratio_sum = plt.subplot()

	#for the summary of the combined XS and ratio plots
	gs_combi = gridspec.GridSpec(3,3)
	fig_combi = plt.figure(figsize=(7,7), dpi=300)
	ax_combi_xs = plt.subplot(gs_combi[:-1,:])			#upper part of plot
	ax_combi_ratio = plt.subplot(gs_combi[-1:,:])			#lower part of plot




	#adjustments that are only done once for the summary plot
	#ax1_sum.axhline(y=1.0, xmin=0, xmax=1, linewidth=0.4, color=_ybysbin_col[bname], linestyle='dashed', dashes=(5,10))
	plt.yscale("log")
	plt.xscale("log")
	plt.grid(True, which="major", linestyle="dotted", color='k')

	#ax_xs_sum.yscale("log")
	#ax_xs_sum.xscale("log")

	'''
	if(datatype=="data"):
		ax1_sum.set_xlim(40, 1200)
		ax1_sum.set_ylim(0.0000001, 0.1)
	elif(datatype=="mc"):
		ax1_sum.set_xlim(10, 2000)
		ax1_sum.set_ylim(1.0, 1000000.0)
	'''

	#---------------------------------------#
	#		plotting for each ybys bin		#
	#---------------------------------------#
	#loop through bins and create one plot for each ybys bin, containing MC and data
	for j in range(0, 6):
		#current bin name:
		bname = bin_names[j]

		#for the individual (ybys bin) XS plot
		fig_xs_bin = plt.figure(figsize=(7,7), dpi=300)
		ax_xs_bin = plt.subplot()

		#for the individual data/mc ratio plot
		fig_ratio_bin = plt.figure(figsize=(7,7), dpi=300)
		ax_ratio_bin = plt.subplot()

		#for the individual combined XS and ratio plot
		gs_combi_bin = gridspec.GridSpec(3,3)
		fig_combi_bin = plt.figure(figsize=(7,7), dpi=300)
		ax_combi_xs_bin = plt.subplot(gs_combi_bin[:-1,:])			#upper part of plot
		ax_combi_ratio_bin = plt.subplot(gs_combi_bin[-1:,:])			#lower part of plot

 
		#reading the data and mc points in the current ybys bin from the rootfile
		#read (statistical?) errors on the histogram entries at the same time
		#furthermore: read x-axis and lower and upper bin bounds
		if(process=="dijet"):
			data, data_err, data_xax, data_low_bb, data_up_bb = read_rootfile(data_rootfile, "Standard/%s/dijet_%s"%(bname, bname))
			mc, mc_err, mc_xax, mc_low_bb, mc_up_bb = read_rootfile(mc_rootfile, "Standard/%s/dijet_%s"%(bname, bname))

		elif(process=="zjet"):
			data, data_err, data_xax, data_low_bb, data_up_bb = read_rootfile(data_rootfile, "Standard/%s/zjet_%s"%(bname, bname))
			mc, mc_err, mc_xax, mc_low_bb, mc_up_bb = read_rootfile(mc_rootfile, "Standard/%s/zjet_%s"%(bname, bname))


		#CALCULATE THE RATIO DATA/MC! ARRAYS! --> how to handle the errors?
		cur_ratio = np.divide(data, mc)



		#-----------------------#
		#		plotting		#
		#-----------------------#
		#creating the plots for the current ybys bin

		#MC and data in same plot
		plot_xs(bname, ax_xs_bin, ax_xs_sum, data_xax, data, data_err, mc, mc_err, m_alpha, process, runperiod, generator)

		#ratio data/mc
		#plot_ratio()

		#combination: xs and ratio
		#plot_combi()
	


		#-------------------#
		# legend and save 	#
		#-------------------#
		#save plots for current ybys bin
		#xs plot:
		#plt.legend()
		ax_xs_bin.legend(fontsize=14, numpoints=1, loc='upper left')
		if((runperiod!=None)and(generator!=None)):
			fig_xs_bin.savefig("%s_xs_data_Run%s_mc_%s.png"%(bname,runperiod,generator))
		else:
			fig_xs_bin.savefig("%s_xs_data_mc.png"%bname)

		#ratio plot:
		ax_ratio_bin.legend(fontsize=14, numpoints=1, loc='upper left')
		if((runperiod!=None)and(generator!=None)):
			fig_ratio_bin.savefig("%s_ratio_data_Run%s_mc_%s.png"%(bname,runperiod,generator))
		else:
			fig_ratio_bin.savefig("%s_ratio_data_mc.png")

		#combination plot
		ax_combi_xs_bin.legend(fontsize=14, numpoints=1, loc='upper left')	#only one legend for upper plot
		if((runperiod!=None)and(generator!=None)):
			fig_combi_bin.savefig("%s_combi_data_Run%s_mc_%s.png"%(bname,runperiod,generator))
		else:
			fig_combi_bin.savefig("%s_combi_data_mc.png"%bname)


	
	#---------------------------------------#
	#	finalising the three summary plots	#
	#---------------------------------------#

	# CROSS SECTION SUMMARY
	#adjustments that are only done once for the summary plot
	ax_xs_sum.set_xlabel(r'$\mathrm{p_{T,avg}}$ /GeV', fontsize=12, horizontalalignment='right')
	ax_xs_sum.xaxis.set_label_coords(1.00, -0.06)
	#plt.xlabel(r'$\mathrm{p_{T,avg}}$ /GeV', horizontalalignment='right', loc='right')
	ax_xs_sum.set_ylabel(r'XS Data and MC', fontsize=12, horizontalalignment='center', verticalalignment='top', y=0.5, rotation=90, labelpad=24)
	if((runperiod!=None)and(generator!=None)):
		ax_xs_sum.set_title(r'Cross Section Summary (Run %s, %s)'%(runperiod,generator), fontsize=14, loc='right', horizontalalignment='right')
	else:
		ax_xs_sum.set_title(r'Cross Section Summary', fontsize=14, loc='right', horizontalalignment='right')

	'''
	#set x-labels of summary plot manually --> more labels, scalar formatting
	locs, labels = plt.xticks()
	#plt.xticks((60, 100, 200, 400, 600, 1000), ('60', '100', '200', '400', '600', '1000'))
	print("xticks: locs, labels: ", locs, labels)
	locs, labels = plt.xticks()
	print("xticks: locs, labels: ", locs, labels)

	#taking care of the tick labels (still to be optimized.)
	#set minor ticks from 1 to 1000 --> 2, 3, 4, 5, 6, 7, 8, 9, 20, 30, 40, 50, 60, 70, 80, 90, etc.
	ax1_sum.set_xticklabels(('', '', '', '', '', '', '', '', '', '', r'$40$','','60', '', '', '', '200', '', '400', '', '600', '', '', ''), minor=True)
	#set major ticks from 1 to 1000 --> 1, 10, 100, 1000
	ax1_sum.set_xticklabels(('1','10','100','1000'), minor=False)
	'''	
	
	#ax1_sum.text(0.06, 0.92, "CMS", fontsize=26, weight='bold', ha='left', va='bottom', transform=ax1.transAxes)                        #text: CMS-logo
	#ax1_sum.text(0.06, 0.86, "private work", fontsize=20, fontstyle='italic', ha='left', va='bottom', transform=ax1.transAxes)          #text: 8TeV info
	###ax1_sum.text(0.02, 1.01, "CMS", fontsize=17, weight='bold', ha='left', va='bottom', transform=ax1.transAxes)                        #text: CMS-logo
	#ax1_sum.text(0.16, 1.01, "Preliminary", fontsize=13, fontstyle='italic', ha='left', va='bottom', transform=ax1.transAxes)          #text: 8TeV info
	###ax1_sum.text(0.16, 1.01, "private work", fontsize=13, fontstyle='italic', ha='left', va='bottom', transform=ax1.transAxes)          #text: private work

	#for now: just give run period as info --> later add more and in a nicer way
	'''
	if(runperiod!=None):
		ax1_sum.text(0.16, 1.01, "2018 data Run %s"%runperiod, fontsize=13, fontstyle='italic', ha='left', va='bottom', transform=ax1.transAxes)		#text: Run period
	else:
		ax1_sum.text(0.16, 1.01, "2018 monte carlo", fontsize=13, fontstyle='italic', ha='left', va='bottom', transform=ax1.transAxes)          		#text: monte carlo
	'''

	plt.legend()
	ax_xs_sum.legend(fontsize=14, numpoints=1, loc='upper left')

	if((generator!=None)and(runperiod!=None)):
		fig_xs.savefig("%s_xs_overview_mc_%s_data_%s.png"%(process,generator,runperiod))
	else:
		fig_xs.savefig("%s_xs_overview_mc_data.png"%(process,generator,runperiod))


	#RATIO SUMMARY

	#COMBIPLOT SUMMARY


	# stop timer
	stop_time = timeit.default_timer()
	timediff = stop_time-start_time
	print("[ratios.py]: Elapsed time: %s sec = %s min" %(timediff, round(timediff/60., 2)))



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
