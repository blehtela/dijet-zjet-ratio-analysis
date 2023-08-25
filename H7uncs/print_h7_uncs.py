#!/usr/bin/env python2
#-*- coding:utf-8 -*-


#####################################################
#													#
#	Printing the relative uncertainties				#
#	on the Herwig7 samples, that were used			#
#	to determine the JES uncertainties.				#
#													#
#	print_h7_uncs.py								#
#	Created by B. Schillinger, 11.04.2020			#
#	Last modified: B.Schillinger, 11.04.2020		#
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
from matplotlib.legend_handler import HandlerBase	#for nice legend entries with two custom lines
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

_ybysbin_xind_up = {'yb0ys0':62, 'yb0ys1':54, 'yb0ys2':40, 'yb1ys0':50, 'yb1ys1':43, 'yb2ys0':37} #contraints coming from fnlo, applied on rootfiles (dijet)

_ybysbin_xind_up_zjetdata = {'yb0ys0':40, 'yb0ys1':30, 'yb0ys2':20, 'yb1ys0':37, 'yb1ys1':30, 'yb2ys0':26}	#constraints coming from rootfiles (zjet data) --> applied on rootfiles

_ybysbin_xind_up_zjetfnlo = {'yb0ys0':-22, 'yb0ys1':-24, 'yb0ys2':-20, 'yb1ys0':-13, 'yb1ys1':-13, 'yb2ys0':-11}	#contraints coming from rootfiles (zjet data) --> applied on fnlo [bins referring to fnlo binning]

#plotting only up to this limit
#only first guesses! check if adjustments desirable! (these have been introduced for dijet) --> yes, different for zjet
_ybysbin_xlim_up = {'yb0ys0':4000, 'yb0ys1':2600, 'yb0ys2':1200, 'yb1ys0':2100, 'yb1ys1':1400, 'yb2ys0':1000}
_ybysbin_xlim_up_zjet = {'yb0ys0':1200, 'yb0ys1':610, 'yb0ys2':220, 'yb1ys0':1000, 'yb1ys1':610, 'yb2ys0':410}


#-----------------------------------------------#
#			Class for nice legends				#
#-----------------------------------------------#
class AnyObjectHandler(HandlerBase):
	def create_artists(self, legend, orig_handle, x0, y0, width, height, fontsize, trans):
		l1 = plt.Line2D([x0,y0+width], [0.7*height,0.7*height], linestyle='solid', color=orig_handle[0])
		l2 = plt.Line2D([x0,y0+width], [0.3*height,0.3*height], linestyle=orig_handle[1], color=orig_handle[0])

		return [l1,l2]


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

	ax.set_xlabel(r'$\mathrm{p_{T,avg}}$ /GeV', fontsize=13, horizontalalignment='right')
	ax.xaxis.set_label_coords(1.00, -0.06)
	#plt.xlabel(r'$\mathrm{p_{T,avg}}$ /GeV', horizontalalignment='right', loc='right')
	ax.set_ylabel(r'XS Data, Theory and MC in fb/GeV', fontsize=12, horizontalalignment='center', verticalalignment='top', y=0.5, rotation=90, labelpad=24)



#function to plot the data/mc ratio for either dijet or zjet
#for individual ybys bins
#def plot_ratio(bname, ax, x_axis, ratio, ratio_err, m_color, m_alpha, process, labeling):
def plot_ratio(bname, ax, x_axis, ratio, ratio_err, low_bb, up_bb, m_color, m_alpha, labeling):
##low_bb and up_bb are used for drawing the "x-error" == binwidth

	x_min= low_bb[0]
	x_max = up_bb[-1]
	print("xmin: %s, xmax: %s" %(x_min, x_max))


	#for x-error-bars
	xerr_low = np.subtract(x_axis, low_bb)
	xerr_up = np.subtract(up_bb, x_axis)
	#zip(xerr_low, xerr_up) #-->only assigns [(low1, up1), (low2, up2), ...]
	#print "xerr_low: \n", xerr_low
	#print "xerr_up: \n", xerr_up

	xerrors = np.array([xerr_low, xerr_up])		#these just show the binwidth

	#plot the ratio as errorbar
	ax.errorbar(x_axis, ratio, xerr=xerrors, yerr=ratio_err, elinewidth=1.2, linewidth=1.0, marker=_ybysbin_marker[bname], ms=6, color=m_color, fillstyle='none', fmt='.', label=labeling)




#function to plot selected unc source(s) as band
#for individual ybys bins
#here: used for total JEC uncertainty
def plot_unc_filled(bname, ax, patches, x_axis, binbounds, low_unc, up_unc, linecol, fillcol, fillalpha, uncname, runperiod, generator):

	ax.set_yscale("linear")
	#ax_ratio.set_ylim(-1.5, 2.0)	#test
	#ax.set_ylim(0.0, 2.0)	#test
	#ax_ratio.get_legend().remove()

	#plot filled area
	ax.fill_between(binbounds.T.flatten(), steppify_bin(low_unc), steppify_bin(up_unc), edgecolor='grey', facecolor=fillcol, alpha=fillalpha)

	#patch for current uncertainty (needed for legend)
	#patches.append(matplotlib.patches.Rectangle((0,0), 0, 0, color=fillcol, label=uncname, alpha=0.6))
	patches.append(mpl.patches.Rectangle((0,0), 0, 0, color=fillcol, label=uncname, alpha=fillalpha))

	ax.add_patch(patches[-1])	#add latest created patch to list of patches

	#set xlim in a way that there is no white-space on the left of the lowest bin / right side of the highest bin
	##ax.set_xlim(binbounds[0,0], binbounds[1,-1])

	###########ax.set_xlim(30, 1200)	#testTESTtest

	#xstart, xend = ax_xs.get_xlim()	#get the xrange of the xs plot, set the same for the ratio plot
	#ax_ratio.set_xlim((xstart, xend))

#function for plotting the individual unc-sources
#plot the up and down ones seperately (despite them being symmetrised)
#def plot_indiv_uncs(bname, ax, patches, x_axis, binbounds, unc, linecol, linestyle, uncname, runperiod, generator):
def plot_indiv_uncs(bname, ax, patches, x_axis, binbounds, unc, linecol, linestyle, uncname, legbool=False):
	#check if legend entry required
	if(legbool):
		labeling=uncname
	else:
		labeling=None

	#need to use steppify_bin for correctly stepped plot (otherwise first step missing or other bugs)
	#ax.step(binbounds[1], unc, color=linecol, linestyle=linestyle, where='pre', label=labeling, alpha=0.9)
	ax.step(binbounds.T.flatten(), steppify_bin(unc), color=linecol, linestyle=linestyle, where='pre', label=labeling, alpha=0.9)



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
	parser.add_argument('dijet_rootfile', type=TFile, nargs='?',
						help='Required argument! File with data results that shall be plotted. (filename glob)')
	parser.add_argument('zjet_rootfile', type=TFile, nargs='?',
						help='Required argument! File with MC results that shall be plotted. (filename glob)')

	# Optional arguments
	parser.add_argument('--alphamarker', '-a', type=float, default=1.0, nargs='?', action='store',
						help='If set to True: draw marker in summary plot with given input alpha.')
	parser.add_argument('--verbose', '-v', default=False, action='store_true',
						help="Increase verbosity if chosen.")
	parser.add_argument('--slicing', '-s', default=False, action='store_true',
						help="Slice output range to available data and fnlo for each process individually.")



	if len(sys.argv)==1:
		parser.print_help(sys.stderr)
		sys.exit(1)

	# Parse arguments
	args = vars(parser.parse_args())

	# input filename
	dijet_rootfile = args['dijet_rootfile']
	zjet_rootfile = args['zjet_rootfile']
	m_alpha = args['alphamarker']
	verbose = args['verbose']
	slicing = args['slicing']

	print('\n')
	print("[print_h7_uncs.py]: Taking dijet input from %s" %dijet_rootfile)
	print("[print_h7_uncs.py]: Taking zjet input from %s" %zjet_rootfile)
	#print("[print_h7_uncs.py]: Taking theory input from %s" %fnlo_tables)

	ybys=None	#fix..


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
	#creating one ratio-plot per bin

	for i in range(bin_names.size): #should usually be 6
		bname = bin_names[i]
		if(verbose==True): print "[print_h7_uncs.py]: Now in bin ", bname

		#for the individual (ybys bin) ratio plot
		fig_ratio_bin = plt.figure(figsize=(9,7), dpi=300)
		ax_ratio_bin = plt.subplot(xscale="log", yscale="linear")
		patches_bin = []


		#reading the data and mc points in the current ybys bin from the rootfile
		#read (statistical?) errors on the histogram entries at the same time
		#furthermore: read x-axis and lower and upper bin bounds
		#dijet
		dijet, dijet_err, dijet_xax, dijet_low_bb, dijet_up_bb = read_rootfile(dijet_rootfile, "CrossSections_byBinwidth/%s/sumhisto_dijet_herwig7_%s"%(bname, bname))
		#zjet
		zjet, zjet_err, zjet_xax, zjet_low_bb, zjet_up_bb = read_rootfile(zjet_rootfile, "CrossSections_byBinwidth/%s/sumhisto_zjet_herwig7_%s"%(bname, bname))


		#-----------------------------------------------------------------------------------#
		#	Remove the low-pTavg bins, that are empty due to dijet trigger (ptavg < 56GeV)	#
		#	Additionally: restrictions from fnlo tables (see above in original version)		#
		#	--> numpy array slicing															#
		#-----------------------------------------------------------------------------------#
		#remove the lower and higher bins that are not present in fastnlo tables
		#get upper limit from dictionary _ybysbin_xind_up
		#data_xax = data_xax[7:62]	#bin 7 is included, bin 62 excluded
		if(slicing==True):
			dijet_xax= dijet_xax[11:_ybysbin_xind_up[bname]]
			dijet = dijet[11:_ybysbin_xind_up[bname]]
			dijet_err = dijet_err[11:_ybysbin_xind_up[bname]]
			dijet_low_bb = dijet_low_bb[11:_ybysbin_xind_up[bname]]
			dijet_up_bb = dijet_up_bb[11:_ybysbin_xind_up[bname]]

			#the new nentries:
			nentries = dijet.size

			#now for zjet
			#remove the lower bins that are not present in fastnlo tables
			#remove the higher bins that are not present in data rootfile
			#get upper limit from dictionary _ybysbin_xind_up_zjetdata
			#data_xax = data_xax[7:62]	#bin 7 is included, bin 62 excluded

			# NEW version: cut entries lower than 56 GeV to be equivalently displayed to dijet
			zjet_xax = zjet_xax[11:_ybysbin_xind_up_zjetdata[bname]]
			zjet = zjet[11:_ybysbin_xind_up_zjetdata[bname]]
			zjet_err = zjet_err[11:_ybysbin_xind_up_zjetdata[bname]]
			zjet_low_bb = zjet_low_bb[11:_ybysbin_xind_up_zjetdata[bname]]
			zjet_up_bb = zjet_up_bb[11:_ybysbin_xind_up_zjetdata[bname]]



		#-------------------------------------------------------------------#
		# 				Calculate the ratios stat.unc./XS					#
		#-------------------------------------------------------------------#

		#calculate the RELATIVE statistical uncertainties
		ratio_dijetdijet = np.divide(dijet, dijet) #this is 1.0 always.
		ratio_zjetzjet = np.divide(zjet, zjet)
		dijet_relerr = np.divide(dijet_err, dijet)
		zjet_relerr = np.divide(zjet_err, zjet)


		print "Relative statistical uncertainties in the Herwig7 Monte Carlo Sets."
		print "These H7 sets were used for calculating the JES uncs."
		print "Identify the ptavg-bins with too high stat. unc. --> remove these in the analysis."
		print "------------------------------------------------------------------------------------"
		print "------------------------------------------------------------------------------------"
		print "DIJET STAT UNCS in bin: %s"%bname
		print "------------------------------------------------------------------------------------"
		print "ptavg: low edge: ", dijet_low_bb
		print "XS in fb/GeV: ", dijet
		print "abs. stat. unc: ", dijet_err
		print "rel. stat. unc: ", dijet_relerr
		print "------------------------------------------------------------------------------------"
		print "------------------------------------------------------------------------------------"
		print "Z+JET STAT UNCS in bin: %s"%bname
		print "------------------------------------------------------------------------------------"
		print "ptavg: low edge: ", zjet_low_bb
		print "XS in fb/GeV: ", zjet
		print "abs. stat. unc: ", zjet_err
		print "rel. stat. unc: ", zjet_relerr
		print "------------------------------------------------------------------------------------"
		print "------------------------------------------------------------------------------------"
		print "\n"
	
		
	
		#-------------------------------------------------------#
		#		plotting the XS ratio of current ybys bin		#
		#-------------------------------------------------------#
		#--- Removed old code!!! can still be found in plot_xs_ratios.py ---#

		#make binbounds array (for flatten() later before using steppify())
		binbounds = np.array([dijet_low_bb, dijet_up_bb])
		binwidths = np.subtract(dijet_up_bb,dijet_low_bb)

		binbounds_zjet = np.array([zjet_low_bb, zjet_up_bb])
		binwidths_zjet = np.subtract(zjet_up_bb,zjet_low_bb)



		nentries = dijet_xax.size
		#make 20% line unc arrays:
		p20_up = np.full(nentries, 1+0.20)
		p20_dw = np.full(nentries, 1-0.20)
	
		moved_xax = np.add(dijet_xax,0.25*binwidths) 
		moved_xax_zjet = np.add(zjet_xax, 0.25*binwidths_zjet)
		#print "dijet_xax: ", dijet_xax
		#print "moved_xax: ", moved_xax

		plot_indiv_uncs(bname, ax_ratio_bin, patches_bin, dijet_xax, binbounds, p20_up, 'tomato', 'dashed', "20% rel.stat.unc")		#lumi
		plot_indiv_uncs(bname, ax_ratio_bin, patches_bin, dijet_xax, binbounds, p20_dw, 'tomato', 'dashed', "20% rel.stat.unc")	#lumi



		#remember that in matplotlib errorbar plots all the given errors (also downwards) have to be absolute (i.e. positive) values
		#therefore: take np.absolute()

		#plot the error bands (mc around 1) --> should rather use plot_indiv_uncs() here
		#band just for total unc!

		#plotting total uncertainty:
		#plot_unc_filled(bname, ax_ratio_bin, patches_bin, data_xax, binbounds, totunc_down, totunc_up, 'mediumseagreen', 'mediumseagreen', 0.6, "JES total", runperiod, generator) #JES unc total (new name)

		#lines -- TO DO --> adjust!! these are just copied from my doubleratio script
		#plot_indiv_uncs(bname, ax_ratio_bin, patches_bin, data_xax, binbounds, TotNoFlav_up, 'royalblue', 'solid', "jesTotalNoFlavor up (H7)")
		#plot_indiv_uncs(bname, ax_ratio_bin, patches_bin, data_xax, binbounds, TotNoFlav_down, 'royalblue', 'dashed', "jesTotalNoFlavor down (H7)")

		#new plotting:
		#ax_ratio_bin.errorbar(data_xax, ratio_datadata, xerr=None, yerr=np.absolute(ratio_data_err), elinewidth=1.2, linewidth=1.0, marker='.', ms=2, color=_ybysbin_col[bname], fillstyle='none', fmt='.', label="data (Run%s) stat. unc."%runperiod) #data stat. unc
		ax_ratio_bin.errorbar(dijet_xax, ratio_dijetdijet, xerr=None, yerr=np.absolute(dijet_relerr), elinewidth=1.2, linewidth=1.0, marker='.', ms=2, color='darkorange', fillstyle='none', fmt='.', label="Dijet relative statistical uncertainty (H7)") #data stat. unc --> in greenyellow, lawngreen, darkorange
		ax_ratio_bin.errorbar(moved_xax_zjet, ratio_zjetzjet, xerr=None, yerr=np.absolute(zjet_relerr), elinewidth=1.2, linewidth=1.0, marker='.', ms=2, color='lawngreen', fillstyle='none', fmt='.', label="Z+jet relative statistical uncertainty (H7)") #data stat. unc --> in greenyellow, lawngreen, darkorange


		#plot_ratio(bname, ax_ratio_bin, data_xax, ratio_mcdata, np.absolute(ratio_mc_err), binbounds[0,:], binbounds[1,:], _ybysbin_col[bname], m_alpha, "MC (%s) with MC stat. unc."%(genlabel))	#mc/data
		##plot_ratio(bname, ax_ratio_bin, data_xax, ratio_mcdata, np.absolute(ratio_mc_err), binbounds[0,:], binbounds[1,:], 'saddlebrown', m_alpha, "MC simulation (%s) with MC stat. unc."%(genlabel))	#mc/data --> in brown


		#save the figures: should be six in total
		#-------------------#
		# legend and save 	#
		#-------------------#
		#create extra legend for MC and theory
		#legend_elements = [Line2D([0],[0],color='saddlebrown', lw=4, label="MC %s in each ybys bin"%(generator)), Line2D([0],[0],color='gray', lw=4, label="LO Table %s in each ybys bin"%(pdfset)), Line2D([0],[0],color='black', lw=4, label="NLO Table %s in each ybys bin"%(pdfset))]
		
		#save plots for current ybys bin
		#xs plot:
		#plt.legend()
		'''
		ax_xs.legend(fontsize=10, numpoints=1, loc='lower left')
		first_legend = ax_xs.get_legend()
		ax_xs.add_artist(first_legend)
		second_legend = plt.legend(handles=legend_elements, loc='upper right')
		'''

		#add line at 1.0
		ax_ratio_bin.axhline(y=1.0, xmin=0, xmax=1, linewidth=0.6, color='grey', linestyle='dashed', dashes=(5,10))	#linewidth=0.4
		#ax_ratio_bin.grid(True, which="major", axis="y", color="k", linestyle="dotted", linewidth=0.8, alpha=0.2)
		ax_ratio_bin.grid(True, which="both", axis="y", color="k", linestyle="dotted", linewidth=0.8, alpha=0.2)
		ax_ratio_bin.grid(True, which="both", axis="x", color="k", linestyle="dotted", linewidth=0.8, alpha=0.2)

		#finalising (used to be done in plot_ratio() fct, but do it only once per bin, here)
		ax_ratio_bin.set_xlabel(r'$\mathrm{p_{T,avg}}$ /GeV', fontsize=13, horizontalalignment='right')
		ax_ratio_bin.xaxis.set_label_coords(1.00, -0.06)
		#ax_ratio_bin.set_ylabel(r'Ratios', fontsize=12, horizontalalignment='center', verticalalignment='top', y=0.5, rotation=90, labelpad=24)
		ax_ratio_bin.set_ylabel(r'Relative uncertainty', fontsize=13, horizontalalignment='center', verticalalignment='top', y=0.5, rotation=90, labelpad=24)


		ax_ratio_bin.set_title("Relative statistical uncertainty (Herwig 7 -- used for JES uncs, %s)"%bname, fontsize=12)
		ax_ratio_bin.legend(fontsize=10, numpoints=1, loc='upper left')


		#---------------- testing --------------#
		'''
		#NEW LEGEND HANDLING! --> make separate unc-legend
		ratiohandles, ratiolabels = ax_ratio_bin.get_legend_handles_labels()	#7 entries (totunc, scaleuncLO, scaleuncNLO, data, MC, LO, NLO)
		#order_entries = [3,4,5,6,0,1,2]
		order_entries = [3,4,5,6]
		ax_ratio_bin.legend([ratiohandles[idx] for idx in order_entries], [ratiolabels[idx] for idx in order_entries], fontsize=10, numpoints=1, loc='upper left')
		ratio_legend = ax_ratio_bin.get_legend()	#get the "old" legend containing all the ratios
		ax_ratio_bin.add_artist(ratio_legend)		#draw the old legend
		#NEW legend with title and different order, shorter labels:
		unc_legend = plt.legend([ratiohandles[0],("royalblue","dashed"),("tomato","dashed"),("khaki","dashed"),ratiohandles[1],ratiohandles[2]], [ratiolabels[0],"JES no flavor",'%s'%flavname, r"Luminosity ($\pm$ 2.5%)",ratiolabels[1],ratiolabels[2]], handler_map={tuple: AnyObjectHandler()}, loc='lower left', title='Uncertainties')	#make the unc legend
		'''
		#--------------end: testing-------------#

		xlocs = ax_xs.get_xticks(minor=True)
		xlabels = ax_xs.get_xticklabels(minor=True)
		#print "xlocs: "
		#print xlocs
		#print "xlabels: "
		#print xlabels

		#xlocs --> 48 minor label locations
		#xlabels --> 48 minor labels

		x_minticks = [50, 60, 70, 80, 90, 200, 300, 400, 500, 600, 700, 800, 900, 2000, 3000, 4000]
		x_minticklabels = ['','60', '', '', '', '200', '300', '', '500','', '700', '', '', '2000', '3000', '']
		x_majticks = [100, 1000]
		x_majticklabels = ['100', '1000']

		ax_ratio_bin.set_xticks(x_majticks, minor=False)	#set major ticks loc
		ax_ratio_bin.set_xticks(x_minticks, minor=True)		#set minor ticks loc
		ax_ratio_bin.set_xticklabels(x_majticklabels, minor=False, fontsize=12)	#set major ticks label
		ax_ratio_bin.set_xticklabels(x_minticklabels, minor=True, fontsize=12)	#set minor ticks label


		#ax_ratio_bin.set_yticks(np.arange(0.0,2.0, step=0.1))
		y_minticks = [0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.10, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 1.90]
		#y_minticklabels = ['', '0.20', '', '0.40', '', '0.60', '', '0.80', '', '', '1.20', '', '1.40', '', '1.60', '', '1.80', '']
		y_minticklabels = ['', '-0.80', '', '-0.60', '', '-0.40', '', '-0.20', '', '', '+0.20', '', '+0.40', '', '+0.60', '', '+0.80', '']
		y_majticks = [0.00, 1.00, 2.00]
		#y_majticklabels = ['0.00', '1.00', '2.00']
		y_majticklabels = ['-1.00', '0.00', '+1.00']
		ax_ratio_bin.set_yticks(y_majticks, minor=False)			#set major ticks loc
		ax_ratio_bin.set_yticks(y_minticks, minor=True)				#set minor ticks loc
		ax_ratio_bin.set_yticklabels(y_majticklabels, minor=False, fontsize=12)	#set major ticks label
		ax_ratio_bin.set_yticklabels(y_minticklabels, minor=True, fontsize=14)	#set minor titcks label

		#fontsize does not work for minor ticks...
		#ax_ratio_bin.tick_params(axis='y', which='minor', labelsize=14)

		for tick in ax_ratio_bin.yaxis.get_minor_ticks():
			tick.label.set_fontsize(10)

		#ax_ratio_bin.set_ylim(0.6, 1.4)	#test zoomed in
		ax_ratio_bin.set_ylim(0.0, 2.0)		#test

		#NEW: set xlim according to ybys bin --> those have been introduced for dijet
		'''
		if(process=="dijet"):
			ax_ratio_bin.set_xlim(50.0, _ybysbin_xlim_up[bname])
		else: #if zjet
			ax_ratio_bin.set_xlim(50.0, _ybysbin_xlim_up_zjet[bname])
		'''

		plt.tight_layout()
		#fig_xs.savefig("%s_xs_overview_theory_Run%s_%s.png"%(process, runperiod, generator))
		#fig_ratio_bin.savefig("%s_xs_ratios_%s.png"%(process, bname))

		if(slicing==True):
			ax_ratio_bin.set_xlim(50.0, _ybysbin_xlim_up[bname])
			fig_ratio_bin.savefig("rel_stat_uncs_MC_H7_%s_sliced.png"%(bname))
		else:
			fig_ratio_bin.savefig("rel_stat_uncs_MC_H7_%s.png"%(bname))


		#zoom in to 20% up and down
		ax_ratio_bin.set_ylim(0.8, 1.2)
		y_minticks = [0.82, 0.84, 0.86, 0.88, 0.92, 0.94, 0.96, 0.98, 1.02, 1.04, 1.06, 1.08, 1.12, 1.14, 1.16, 1.18]
		#y_minticklabels = ['0.82', '0.84', '0.86', '0.88', '0.92', '0.94', '0.96', '0.98', '1.02', '1.04', '1.06', '1.08', '1.12', '1.14', '1.16', '1.18']
		y_minticklabels = ['-0.18', '-0.16', '-0.14', '-0.12', '-0.08', '-0.06', '-0.04', '-0.02', '+0.02', '+0.04', '+0.06', '+0.08', '+0.12', '+0.14', '+0.16', '+0.18']
		y_majticks = [0.80, 0.90, 1.00, 1.10, 1.20]
		#y_majticklabels = ['0.80', '0.90', '1.00', '1.10', '1.20']
		y_majticklabels = ['-0.20', '-0.10', '0.00', '+0.10', '+0.20']
		ax_ratio_bin.set_yticks(y_majticks, minor=False)			#set major ticks loc
		ax_ratio_bin.set_yticks(y_minticks, minor=True)				#set minor ticks loc
		ax_ratio_bin.set_yticklabels(y_majticklabels, minor=False, fontsize=12)	#set major ticks label
		ax_ratio_bin.set_yticklabels(y_minticklabels, minor=True, fontsize=14)	#set minor titcks label

		#fontsize does not work for minor ticks...
		#ax_ratio_bin.tick_params(axis='y', which='minor', labelsize=14)

		for tick in ax_ratio_bin.yaxis.get_minor_ticks():
			tick.label.set_fontsize(10)

		plt.tight_layout()
		if(slicing==True):
			ax_ratio_bin.set_xlim(50.0, _ybysbin_xlim_up[bname])
			fig_ratio_bin.savefig("rel_stat_uncs_MC_H7_%s_zoomedTo20_sliced.png"%(bname))
		else:
			fig_ratio_bin.savefig("rel_stat_uncs_MC_H7_%s_zoomedTo20.png"%(bname))


	# stop timer
	stop_time = timeit.default_timer()
	timediff = stop_time-start_time
	print("[print_h7_uncs.py]: Elapsed time: %s sec = %s min" %(timediff, round(timediff/60., 2)))



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








