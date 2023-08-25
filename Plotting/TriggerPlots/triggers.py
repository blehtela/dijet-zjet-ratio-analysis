#!/usr/bin/env python2
#-*- coding:utf-8 -*-

#################################################################
#															   	#
# Plotting the different steps of the trigger fits				#
# Taking input from dijet analysis carried out within ROOT		#
#													   			#
# Created by B. Schillinger, 03.09.2019				 			#
# Last modified by B. Schillinger, 03.09.2019		   			#
#													   			#
#################################################################

#--------------------- UNDER CONSTRUCTION -----------------------------------------------------#
# ---> Rewriting my ratio plotting script ratios.py in a way that it plots my trigger studies! #

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
_ybysbin_label = {'yb0ys0':r'$0 \leq y_b < 1$	$0 \leq y^{\ast} < 1$', 'yb0ys1':r'$0 \leq y_b < 1$	$1 \leq y^{\ast} < 2$', 'yb0ys2':r'$0 \leq y_b < 1$	$2 \leq y^{\ast} < 3$', 'yb1ys0':r'$1 \leq y_b < 2$	$0 \leq y^{\ast} < 1$', 'yb1ys1':r'$1 \leq y_b < 2$	$1 \leq y^{\ast} < 2$', 'yb2ys0':r'$2 \leq y_b < 3$	$0 \leq y^{\ast} < 1$'}


#----------------------------------------------------------------------------------------------------------#
#												  parser												   # 
#----------------------------------------------------------------------------------------------------------#
#function for the parsing job
#defining arguments and options
def args_parsing():
	parser = argparse.ArgumentParser(epilog='', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	# Positional arguments
	parser.add_argument('rootfile', type=TFile, nargs='?',
						help='Required argument! File with results that shall be plotted. (filename glob)')
	# Optional arguments

	# Parse arguments
	args = vars(parser.parse_args())
	return args


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
def plot_ratio(bname, ax, ax_sum, x_axis, ratio, binerr, low_bb, up_bb): ##low_bb and up_bb are used for drawing the "x-error" == binwidth

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

	ax.axhline(y=1.0, xmin=0, xmax=1, linewidth=0.4, color=_ybysbin_col[bname], linestyle='dashed', dashes=(5,10))
	ax.errorbar(x_axis, ratio, xerr=xerrors, yerr=binerr, elinewidth=1.2, linewidth=0.0, marker=_ybysbin_marker[bname], ms=12, color=_ybysbin_col[bname], fillstyle='none', fmt='.', label=_ybysbin_label[bname])	#label=bname
	plt.yscale("log")
	plt.xscale("log")

	plt.grid(True, which="major", linestyle="dotted", color='k')

	ax.set_xlim(40, 1200)
	ax.set_ylim(y_min+0.0000001, y_max+0.1)

	ax.set_xlabel(r'$\mathrm{p_{T,avg}}$ /GeV', fontsize=12, horizontalalignment='right')
	ax.xaxis.set_label_coords(1.00, -0.06)
	#plt.xlabel(r'$\mathrm{p_{T,avg}}$ /GeV', horizontalalignment='right', loc='right')
	ax.set_ylabel(r'XS ratios: Zjet / Dijet', fontsize=12, horizontalalignment='center', verticalalignment='top', y=0.5, rotation=90, labelpad=24)
	ax.set_title(r'%s'%bname)

	ax.legend(fontsize=14, numpoints=1, loc='upper left')

	#add it to the total plot
	ax_sum.errorbar(x_axis, ratio, xerr=xerrors, yerr=binerr, elinewidth=1.2, linewidth=0.0, marker=_ybysbin_marker[bname], ms=12, color=_ybysbin_col[bname], fillstyle='none', fmt='.', label=_ybysbin_label[bname])




###### THE FOLLOWING FUNCTION IS NOT USED!! REMAINS FROM OLD IMPLEMENTATION!!! ########################
#function for plotting the experimental and theoretical uncertainty as "steppified band"
def plot_unc_filled(ax, patches, x_axis, binbounds, low_unc, up_unc, linecol, fillcol, uncname):

	#ax.semilogx(x_axis, up_unc, '.', ms=4, ls='dashed', linewidth=0.6, color=linecol, label='')
	#plot upper and lower line
	
	print "steppify_bin(low_unc, True): \n", steppify_bin(low_unc, True)
	print "steppify_bin(low_unc, False): \n", steppify_bin(low_unc, False)
	print "low_unc: \n", low_unc
	print "alternative option: \n", np.array([low_unc, low_unc]).T.flatten()  #would give the same result as using steppify_bin(low_unc, False)
	

	#plot filled area (watch out with steppify... could use step?)
	ax.fill_between(binbounds.T.flatten(), steppify_bin(low_unc), steppify_bin(up_unc), edgecolor='grey', facecolor=fillcol, alpha=0.6)
	#ax.fill_between(x_axis, low_unc, up_unc, color=fillcol, alpha=0.3)


	#patch for current uncertainty (needed for legend)
	patches.append(matplotlib.patches.Rectangle((0,0), 0, 0, color=fillcol, label=uncname, alpha=0.6))
	ax.add_patch(patches[-1]) #add latest created patch to list of patches

	#set xlim in a way that there is no white-space on the left of the lowest bin / right side of the highest bin
	ax.set_xlim(binbounds[0,0], binbounds[1, -1])
	#print "bb[0,0]= %s and bb[1,-1]= %s"%(binbounds[0,0], binbounds[1,-1]) #just a check

	ax.set_xlabel(r'\mathrm{p_{T,avg} /GeV}')
	ax.set_ylabel(r'Zjet / Dijet', horizontalalginment='center', verticalalignment='top', y=0.5, rotation=90, labelpad=24)
	ax.set_title(r'%s'%bname)




#-----------------------------------------------------------------------------------------------#
#											  main-function										# 
#-----------------------------------------------------------------------------------------------#

def main():
	# start timer
	# measuring wall clock time - not necessary, just for fun.
	start_time = timeit.default_timer()

	# getting the parsed arguments
	args = args_parsing()

	# input filename
	rootfile = args['rootfile']
	#rootfile = TFile(args['rootfile'])
	print('\n')
	print("[ratios.py]: Taking input from %s" %rootfile)

	#prepare the lists (will later be arrays)
	bin_names = np.array(["yb0ys0", "yb0ys1", "yb0ys2", "yb1ys0", "yb1ys1", "yb2ys0"])

	#for the summary (all bins)
	gs = gridspec.GridSpec(3,3)
	fig_tot = plt.figure(figsize=(7,7), dpi=300)
	#ax1 = plt.subplot(gs[:-1,:])
	ax1_sum = plt.subplot()
	#patches = []

	#adjustments that are only done once for the summary plot
	#ax1_sum.axhline(y=1.0, xmin=0, xmax=1, linewidth=0.4, color=_ybysbin_col[bname], linestyle='dashed', dashes=(5,10))
	plt.yscale("log")
	plt.xscale("log")
	plt.grid(True, which="major", linestyle="dotted", color='k')

	ax1_sum.set_xlim(40, 1200)
	ax1_sum.set_ylim(0.0000001, 0.1)

 
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

		#plot the ratio for this bin
		plot_ratio(bname, ax1, ax1_sum, xax, ratio, err, low_bb, up_bb)


		plt.legend()
		ax1.legend(fontsize=14, numpoints=1, loc='upper left')
		fig.savefig("%s_ratio.png"%bname)


	#total plot --> summary
	#adjustments that are only done once for the summary plot
	ax1_sum.set_xlabel(r'$\mathrm{p_{T,avg}}$ /GeV', fontsize=12, horizontalalignment='right')
	ax1_sum.xaxis.set_label_coords(1.00, -0.06)
	#plt.xlabel(r'$\mathrm{p_{T,avg}}$ /GeV', horizontalalignment='right', loc='right')
	ax1_sum.set_ylabel(r'XS ratios: Zjet / Dijet', fontsize=12, horizontalalignment='center', verticalalignment='top', y=0.5, rotation=90, labelpad=24)
	ax1_sum.set_title(r'Ratio Summary', fontsize=14, loc='right', horizontalalignment='right')

	#set x-labels of summary plot manually --> more labels, scalar formatting
	locs, labels = plt.xticks()
	#plt.xticks((60, 100, 200, 400, 600, 1000), ('60', '100', '200', '400', '600', '1000'))
	print("xticks: locs, labels: ", locs, labels)
	locs, labels = plt.xticks()
	print("xticks: locs, labels: ", locs, labels)
	#ax1_sum.set_xlim(40,1200, auto=True)
	#set minor ticks from 1 to 1000 --> 2, 3, 4, 5, 6, 7, 8, 9, 20, 30, 40, 50, 60, 70, 80, 90, etc.
	ax1_sum.set_xticklabels(('', '', '', '', '', '', '', '', '', '', r'$40$','','60', '', '', '', '200', '', '400', '', '600', '', '', ''), minor=True)
	#set major ticks from 1 to 1000 --> 1, 10, 100, 1000
	ax1_sum.set_xticklabels(('1','10','100','1000'), minor=False)

	#ax1_sum.text(0.06, 0.92, "CMS", fontsize=26, weight='bold', ha='left', va='bottom', transform=ax1.transAxes)                        #text: CMS-logo
	#ax1_sum.text(0.06, 0.86, "private work", fontsize=20, fontstyle='italic', ha='left', va='bottom', transform=ax1.transAxes)          #text: 8TeV info
	ax1_sum.text(0.02, 1.01, "CMS", fontsize=17, weight='bold', ha='left', va='bottom', transform=ax1.transAxes)                        #text: CMS-logo
	ax1_sum.text(0.16, 1.01, "Preliminary", fontsize=13, fontstyle='italic', ha='left', va='bottom', transform=ax1.transAxes)          #text: 8TeV info





	ax1_sum.legend(fontsize=14, numpoints=1, loc='upper left')

	plt.legend()
	ax1_sum.legend(fontsize=14, numpoints=1, loc='upper left')
	fig_tot.savefig("ratio_summary.png")






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
