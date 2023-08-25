#!/usr/bin/env python2
#-*- coding:utf-8 -*-


#####################################################
#													#
#	XS ratios overview plots :						#
#	coparing data, theory, MC						#
#	individually for either dijet or zjet;			#
#	one plot per ybys bin.							#
#	Taking input from the ROOT file					#
#	and from fastNLO tables							#
#													#
#	This is the updated version:					#
#	> now data as denominator						#
#	> some changes in displaying the uncertainties	#
#													#
#													#
#	Original version can be found in:				#
#	plot_xs_ratios_old_DataOverTheory.py			#
#													#
#	Created by B. Schillinger, 25.02.2020			#
#	Last modified: B.Schillinger, 15.03.2020		#
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


	#ax.grid(True, which="major", axis="y", color="k", linestyle="dotted", linewidth=0.8, alpha=0.2)
	#ax.grid(True, which="both", axis="x", color="k", linestyle="dotted", linewidth=0.8, alpha=0.2)


	#plt.yscale("log")
	#plt.xscale("log")
	#ax.grid(True, which="both", axis="y", color="k", linestyle="dotted", linewidth=1, alpha=0.8)
	#plt.grid(True, which="major", linestyle="dotted", color='k', alpha=0.7)

	#ax.set_xlim(40, 1200)
	#ax.set_ylim(y_min+0.0000001, y_max+0.1)

	'''
	ax.set_xlabel(r'$\mathrm{p_{T,avg}}$ /GeV', fontsize=12, horizontalalignment='right')
	ax.xaxis.set_label_coords(1.00, -0.06)
	ax.set_ylabel(r'Ratios', fontsize=12, horizontalalignment='center', verticalalignment='top', y=0.5, rotation=90, labelpad=24)

	#ax.set_title(r'%s in %s (%s)'%(process,bname, _ybysbin_label[bname]))
	ax.legend(fontsize=10, numpoints=1, loc='upper left')

	plt.tick_params(axis='x', which='minor')				#?? not here but in main code
	ax.xaxis.set_minor_formatter(FormatStrFormatter("%.d"))	# same???
	'''	



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
	ax.set_xlim(30, 1200)	#test

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
	parser.add_argument('data_rootfile', type=TFile, nargs='?',
						help='Required argument! File with data results that shall be plotted. (filename glob)')
	parser.add_argument('mc_rootfile', type=TFile, nargs='?',
						help='Required argument! File with MC results that shall be plotted. (filename glob)')
	#parser.add_argument('fnlo_table', nargs='?',
	#					help='Required argument! FastNLO table to be evaluated.')
	parser.add_argument('fnlo_tables', type=str, nargs=6,	#nargs='+'
						help='Required argument! FastNLO tables for the six ybys bins to be evaluated.')
	parser.add_argument('mc_uncfile', type=TFile, nargs='?',
						help='Required argument! TFile with JEC uncertainties.')

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
	mc_uncfile = args['mc_uncfile']		#uncertainties for this analysis (dijet or zjet)
	process = args['process']			#specifies whether looking at dijet or zjet
	ybys = args['ybysbin']				#name of ybys bin
	m_alpha = args['alphamarker']
	pdfset = args['pdfset']
	runperiod = args['runperiod']	#optional argument, but good for labeling --> to be improved
	generator = args['generator']	#optional argument, pythia or herwig
	verbose = args['verbose']

	print args['asymmetric']

	#scale variation choice for uncertainty:
	if args['asymmetric']:
		scale_var_type = fastnlo.kAsymmetricSixPoint
		variation_type = 'scale uncertainty (6P)'
	else:
		scale_var_type = fastnlo.kSymmetricTwoPoint
		variation_type = 'scale uncertainty (2P)'


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
	#creating one ratio-plot per bin

	for i in range(bin_names.size): #should usually be 6
		bname = bin_names[i]
		if(verbose==True): print "[plot_xs_overview.py]: Now in bin ", bname

		#for the individual (ybys bin) ratio plot
		fig_ratio_bin = plt.figure(figsize=(9,7), dpi=300)
		ax_ratio_bin = plt.subplot(xscale="log", yscale="linear")
		patches_bin = []


		#reading the data and mc points in the current ybys bin from the rootfile
		#read (statistical?) errors on the histogram entries at the same time
		#furthermore: read x-axis and lower and upper bin bounds
		data, data_err, data_xax, data_low_bb, data_up_bb = read_rootfile(data_rootfile, "CrossSections_byBinwidth/%s/sumhisto_%s_Run%s_%s"%(bname, process, runperiod, bname))
		mc, mc_err, mc_xax, mc_low_bb, mc_up_bb = read_rootfile(mc_rootfile, "CrossSections_byBinwidth/%s/sumhisto_%s_%s_%s"%(bname, process, generator, bname))

		#fnloxs = read_fnlotable(fnlo_table)
		if(verbose==True): print "[plot_xs_overview.py]: Reading fnlo table: ", fnlo_tables[i]
		fnloxs_LO, scale_unc_LO, pdf_unc_LO = read_fnlotable(fnlo_tables[i], pdfset, 0, scale_var_type)
		fnloxs_NLO, scale_unc_NLO, pdf_unc_NLO = read_fnlotable(fnlo_tables[i], pdfset, 1, scale_var_type)

		#read uncertainties for dijet (obtained via herwig 7 MC set)
		#only have to read the first argument as the rest should be the same as for data
		unc_TotalNoFlavor = read_rootfile(mc_uncfile, "UncertaintyRatios/%s/UpDownDiv_TotalNoFlavor_%s"%(bname, bname))[0]
		#always read both flavor uncs, but for dijet plot QCD, for zjet plot ZJet
		unc_FlavorQCD = read_rootfile(mc_uncfile, "UncertaintyRatios/%s/UpDownDiv_FlavorQCD_%s"%(bname, bname))[0]
		unc_FlavorZJet = read_rootfile(mc_uncfile, "UncertaintyRatios/%s/UpDownDiv_FlavorZJet_%s"%(bname, bname))[0]

		#Calculate the uncertainties
		unc_size = unc_TotalNoFlavor.size	#should be the same as the entries of data
		ones_arr = np.full(unc_size, 1.0)
		TotNoFlav_up = 0.5*np.add(-1*ones_arr, unc_TotalNoFlavor)
		TotNoFlav_down = -1.0*TotNoFlav_up

		if(process=="dijet"):
			Flav_up = 0.5*np.add(-1*ones_arr, unc_FlavorQCD)
			flavname = "jesQCDFlavor"
		else:	#if zjet
			Flav_up = 0.5*np.add(-1*ones_arr, unc_FlavorZJet)
			flavname = "jesZJetFlavor"

		Flav_down = -1.0*Flav_up


		#Calculate total uncertainty
		delta_tot_up = np.sqrt(pow(TotNoFlav_up,2)+pow(Flav_up,2))			#quadratic sum
		delta_tot_down = np.sqrt(pow(TotNoFlav_down,2)+pow(Flav_down,2))

		print "##########################"
		print " The uncertainties in %s:"%bname
		print " TotNoFlav_up: \n", TotNoFlav_up
		print " TotNoFlav_down: \n", TotNoFlav_down
		print " Flav_up:		\n", Flav_up
		print " Flav_down:		\n", Flav_down
		print " delta_tot_up: \n", delta_tot_up
		print " delta_tot_down: \n", delta_tot_down
		print "---------------------------------"

		#add 1.0 for later drawing them around the 1.0 line in the ratio
		TotNoFlav_up = np.add(ones_arr, TotNoFlav_up)
		TotNoFlav_down = np.add(ones_arr, TotNoFlav_down)
		Flav_up = np.add(ones_arr, Flav_up)
		Flav_down = np.add(ones_arr, Flav_down)
		totunc_up = np.add(ones_arr, delta_tot_up)
		totunc_down = np.add(ones_arr, -1.0*delta_tot_down)
		

		#as long as we do not calculate pdf and scale uncs
		#TO DO: include ScaleUnc and PDFUnc calculation in read_fnlotable() function!!!
		nentries = fnloxs_LO.size

		#fnloLO_err = np.full(nentries,0)
		#fnloNLO_err = np.full(nentries,0)


		#check the absolute unc by giving output
		if(verbose==True):
			print "Scale unc NLO down: \n", scale_unc_NLO[0]
			print "Scale unc NLO up: \n", scale_unc_NLO[1]

		'''
		#calculate error using scale unc and pdf unc (quadratic sum)
		#fnlo_relerrNLO = np.sqrt(pow(rel_scale_unc_NLO, 2), pow(rel_pdf_unc_NLO, 2))
		#fnlo_errNLO = np.multiply(fnloxs_NLO, fnlo_relerrNLO)
		fnlo_errLO_down = np.sqrt(pow(scale_unc_LO[0,:],2)+pow(pdf_unc_LO[0,:],2))		#LO down
		fnlo_errLO_up = np.sqrt(pow(scale_unc_LO[1,:],2)+pow(pdf_unc_LO[1,:],2))		#LO up
		fnlo_errNLO_down = np.sqrt(pow(scale_unc_NLO[0,:],2)+pow(pdf_unc_NLO[0,:],2))	#NLO down
		fnlo_errNLO_up = np.sqrt(pow(scale_unc_NLO[1,:],2)+pow(pdf_unc_NLO[1,:],2))		#NLO up


		#fnloLO_err = zip(fnlo_errLO_down, fnlo_errLO_up)	#would make tuples
		#fnloNLO_err = zip(fnlo_errNLO_down, fnlo_errNLO_up)

		fnloLO_err = np.zeros([2, nentries])	#will contain absolute uncertainties
		fnloLO_err[0] = fnlo_errLO_down
		fnloLO_err[1] = fnlo_errLO_up

		fnloNLO_err = np.zeros([2, nentries])	#will contain absolute uncertainties
		fnloNLO_err[0] = fnlo_errNLO_down
		fnloNLO_err[1] = fnlo_errNLO_up
		'''		
		

		#convert from pb to fb (multiply by 1000) --> now already handled in read_fnlotable()
		#fnloxs_LO = 1000*fnloxs_LO
		#fnloxs_NLO = 1000*fnloxs_NLO


		#-----------------------------------------------------------------------------------#
		#	Remove the low-pTavg bins, that are empty due to dijet trigger (ptavg < 56GeV)	#
		#	Additionally: restrictions from fnlo tables (see above in original version)		#
		#	--> numpy array slicing															#
		#-----------------------------------------------------------------------------------#
		if(process=="dijet"):
			#remove the lower and higher bins that are not present in fastnlo tables
			#get upper limit from dictionary _ybysbin_xind_up
			#data_xax = data_xax[7:62]	#bin 7 is included, bin 62 excluded
			data_xax = data_xax[11:_ybysbin_xind_up[bname]]
			data = data[11:_ybysbin_xind_up[bname]]
			data_err = data_err[11:_ybysbin_xind_up[bname]]
			data_low_bb = data_low_bb[11:_ybysbin_xind_up[bname]]
			data_up_bb = data_up_bb[11:_ybysbin_xind_up[bname]]
			mc_err = mc_err[11:_ybysbin_xind_up[bname]]
			mc = mc[11:_ybysbin_xind_up[bname]]

			#do the same for the uncertainties!
			TotNoFlav_up = TotNoFlav_up[11:_ybysbin_xind_up[bname]]
			TotNoFlav_down = TotNoFlav_down[11:_ybysbin_xind_up[bname]]
			Flav_up = Flav_up[11:_ybysbin_xind_up[bname]]
			Flav_down = Flav_down[11:_ybysbin_xind_up[bname]]
			totunc_up = totunc_up[11:_ybysbin_xind_up[bname]]
			totunc_down = totunc_down[11:_ybysbin_xind_up[bname]]

			#now the fnlo stuff with different binning --> cut away the four lowest bins (< 56 GeV)
			fnloxs_LO = fnloxs_LO[4:]
			fnloxs_NLO = fnloxs_NLO[4:]
			#fnloLO_err = fnloLO_err[:,4:]		#take everything on first axis = up and down, and entries from index 4 onwards (>56GeV)
			#fnloNLO_err = fnloNLO_err[:,4:]	

			scale_unc_LO = scale_unc_LO[:,4:]
			scale_unc_NLO = scale_unc_NLO[:,4:]
			pdf_unc_LO = pdf_unc_LO[:,4:]
			pdf_unc_NLO = pdf_unc_NLO[:,4:]


			#cross sections and statistical errors
			#data_xax = data_xax[4:]		#keep everything above 56 GeV
			#systematic errrors (JEC)

			#the new nentries:
			nentries = data.size
		else:	#if zjet
			#remove the lower bins that are not present in fastnlo tables
			#remove the higher bins that are not present in data rootfile
			#get upper limit from dictionary _ybysbin_xind_up_zjetdata
			#data_xax = data_xax[7:62]	#bin 7 is included, bin 62 excluded
			data_xax = data_xax[7:_ybysbin_xind_up_zjetdata[bname]]
			data = data[7:_ybysbin_xind_up_zjetdata[bname]]
			data_err = data_err[7:_ybysbin_xind_up_zjetdata[bname]]
			data_low_bb = data_low_bb[7:_ybysbin_xind_up_zjetdata[bname]]
			data_up_bb = data_up_bb[7:_ybysbin_xind_up_zjetdata[bname]]
			mc_err = mc_err[7:_ybysbin_xind_up_zjetdata[bname]]
			mc = mc[7:_ybysbin_xind_up_zjetdata[bname]]

			#do the same for the uncertainties!
			TotNoFlav_up = TotNoFlav_up[7:_ybysbin_xind_up_zjetdata[bname]]
			TotNoFlav_down = TotNoFlav_down[7:_ybysbin_xind_up_zjetdata[bname]]
			Flav_up = Flav_up[7:_ybysbin_xind_up_zjetdata[bname]]
			Flav_down = Flav_down[7:_ybysbin_xind_up_zjetdata[bname]]
			totunc_up = totunc_up[7:_ybysbin_xind_up_zjetdata[bname]]
			totunc_down = totunc_down[7:_ybysbin_xind_up_zjetdata[bname]]

			#now the fnlo stuff with different binning --> cut away the higher bins that are no present in data
			fnloxs_LO = fnloxs_LO[:_ybysbin_xind_up_zjetfnlo[bname]]
			fnloxs_NLO = fnloxs_NLO[:_ybysbin_xind_up_zjetfnlo[bname]]
			#fnloLO_err = fnloLO_err[:,:_ybysbin_xind_up_zjetfnlo[bname]]		#take everything on first axis = up and down, and entries from index 4 onwards (>56GeV)
			#fnloNLO_err = fnloNLO_err[:,:_ybysbin_xind_up_zjetfnlo[bname]]	
			scale_unc_LO = scale_unc_LO[:,:_ybysbin_xind_up_zjetfnlo[bname]]		#take everything on first axis = up and down, and entries from index 4 onwards (>56GeV)
			scale_unc_NLO = scale_unc_NLO[:,:_ybysbin_xind_up_zjetfnlo[bname]]	
			pdf_unc_LO = pdf_unc_LO[:,:_ybysbin_xind_up_zjetfnlo[bname]]
			pdf_unc_NLO = pdf_unc_NLO[:,:_ybysbin_xind_up_zjetfnlo[bname]]


			#cross sections and statistical errors
			#data_xax = data_xax[4:]		#keep everything above 56 GeV
			#systematic errrors (JEC)

			#the new nentries:
			nentries = data.size

		#-------------------------------------------------------------------#
		# 	Calculate the ratios theory/data and SCALE their errors	by data #
		#-------------------------------------------------------------------#
		# so far: data/MC, data/LO, data/NLO --> should add for dijet: three MC sets, zjet: two MC sets
		##ratio_datamc = np.divide(data, mc)
		##ratio_datafnloLO = np.divide(data, fnloxs_LO)
		##ratio_datafnloNLO = np.divide(data, fnloxs_NLO)
		ratio_mcdata = np.divide(mc, data)
		print "data: --> ", data
		print "mc: -->", mc
		print("ratio_mcdata: \n ", ratio_mcdata)
		ratio_fnloLOdata = np.divide(fnloxs_LO, data)
		ratio_fnloNLOdata = np.divide(fnloxs_NLO, data)
		ratio_datadata = np.divide(data, data) #--> must be an array of 1.0 (used to draw data stat. error around 1.0)

		#Scale the errors by data:
		#data error will be drawn around 1.0 line
		#data_err = data_err
		# MC Stat. UNC #
		ratio_mc_err = np.divide(mc_err, data)
		ratio_data_err = np.divide(data_err, data)
		
		# SCALE UNC #
		ratio_fnloLO_scaleunc = np.zeros([2, nentries])
		ratio_fnloLO_scaleunc[0] = np.divide(scale_unc_LO[0],data)		#scale scaleunc by data
		ratio_fnloLO_scaleunc[1] = np.divide(scale_unc_LO[1],data)

		ratio_fnloNLO_scaleunc = np.zeros([2, nentries])
		ratio_fnloNLO_scaleunc[0] = np.divide(scale_unc_NLO[0],data)	#scale scaleunc by data
		ratio_fnloNLO_scaleunc[1] = np.divide(scale_unc_NLO[1],data)

		#PDF UNC#
		ratio_fnloLO_pdfunc = np.zeros([2, nentries])
		ratio_fnloLO_pdfunc[0] = np.divide(pdf_unc_LO[0],data)			#scale pdfunc by data
		ratio_fnloLO_pdfunc[1] = np.divide(pdf_unc_LO[1],data)

		ratio_fnloNLO_pdfunc = np.zeros([2, nentries])
		ratio_fnloNLO_pdfunc[0] = np.divide(pdf_unc_NLO[0],data)		#scale pdfunc by data
		ratio_fnloNLO_pdfunc[1] = np.divide(pdf_unc_NLO[1],data)
	
		#create values for scaleunc-band around theory/data ratio
		scaleunc_band_LO = np.zeros([2,nentries])
		#scaleunc_band_LO[0] = np.subtract(ratio_fnloLOdata, ratio_fnloLO_scaleunc[0])	#lower line in scaleunc band
		scaleunc_band_LO[0] = np.add(ratio_fnloLOdata, ratio_fnloLO_scaleunc[0])		#lower line in scaleunc band --> scaleunc itself is already neg or pos
		scaleunc_band_LO[1] = np.add(ratio_fnloLOdata, ratio_fnloLO_scaleunc[1])		#upper line in scaleunc band
		scaleunc_band_NLO = np.zeros([2,nentries])
		scaleunc_band_NLO[0] = np.add(ratio_fnloNLOdata, ratio_fnloNLO_scaleunc[0])#lower line in scaleunc band
		scaleunc_band_NLO[1] = np.add(ratio_fnloNLOdata, ratio_fnloNLO_scaleunc[1])		#upper line in scaleunc band

		print "SCALE UNC NLO", scale_unc_NLO
		print "ScaleuncNLO/data", ratio_fnloNLO_scaleunc
		print "ScaleuncBandNLO", scaleunc_band_NLO
		print "ratio NLO pdfunc", ratio_fnloNLO_pdfunc
		print "ratio NLO", ratio_fnloNLOdata

		#make lumi unc arrays:
		lumi_up = np.full(nentries, 1+0.025)
		lumi_dw = np.full(nentries, 1-0.025)
	

		'''
		ratio_fnloLO_err = np.zeros([2, nentries])
		ratio_fnloLO_err[0] = np.multiply(ratio_datafnloLO, ratio_fnloLO_relerr_dw)
		ratio_fnloLO_err[1] = np.multiply(ratio_datafnloLO, ratio_fnloLO_relerr_up)

		ratio_fnloNLO_err = np.zeros([2, nentries])
		ratio_fnloNLO_err[0] = np.multiply(ratio_datafnloNLO, ratio_fnloNLO_relerr_dw)
		ratio_fnloNLO_err[1] = np.multiply(ratio_datafnloNLO, ratio_fnloNLO_relerr_up)
		'''


		'''
		#try simple error propagation for absolute error of ratio:
		# ratio_err = Ratio*Sqrt((data_err/data)^2 + (mc_err/mc)^2)
		# calculate relative errors:
		data_relerr = np.divide(data_err, data)	# relative error
		mc_relerr = np.divide(mc_err, mc)		# relative error
		fnlo_LO_relerr_down = np.divide(fnloLO_err[0], fnloxs_LO)		#rel err
		fnlo_LO_relerr_up = np.divide(fnloLO_err[1], fnloxs_LO)			#rel err
		fnlo_NLO_relerr_down = np.divide(fnloNLO_err[0], fnloxs_NLO)	#rel err
		fnlo_NLO_relerr_up = np.divide(fnloNLO_err[1], fnloxs_NLO)		#rel err


		#calculate ratio relative error:
		ratio_relerr = np.sqrt(pow(data_relerr,2)+pow(mc_relerr,2))	#numpy.sqrt() returns sqrt for array elementwise
		ratio_fnloLO_relerr_up = np.sqrt(pow(data_relerr,2)+pow(fnlo_LO_relerr_up,2))	#end up with abs values...
		ratio_fnloLO_relerr_dw = np.sqrt(pow(data_relerr,2)+pow(fnlo_LO_relerr_down,2))
		ratio_fnloNLO_relerr_up = np.sqrt(pow(data_relerr,2)+pow(fnlo_NLO_relerr_up,2))
		ratio_fnloNLO_relerr_dw = np.sqrt(pow(data_relerr,2)+pow(fnlo_NLO_relerr_down,2))


		#calulate ratio absolute error:
		ratio_err = np.multiply(ratio_datamc, ratio_relerr)	#elementwise numpy multiplication of arrays

		ratio_fnloLO_err = np.zeros([2, nentries])
		ratio_fnloLO_err[0] = np.multiply(ratio_datafnloLO, ratio_fnloLO_relerr_dw)
		ratio_fnloLO_err[1] = np.multiply(ratio_datafnloLO, ratio_fnloLO_relerr_up)

		ratio_fnloNLO_err = np.zeros([2, nentries])
		ratio_fnloNLO_err[0] = np.multiply(ratio_datafnloNLO, ratio_fnloNLO_relerr_dw)
		ratio_fnloNLO_err[1] = np.multiply(ratio_datafnloNLO, ratio_fnloNLO_relerr_up)
		'''


		#-------------------------------------------------------#
		#		plotting the XS ratio of current ybys bin		#
		#-------------------------------------------------------#
		'''
	#	plot_xs(bname,ax_xs, data_xax, data, data_err, _ybysbin_col[bname], m_alpha, process, "Data Run%s in %s"%(runperiod, _ybysbin_label[bname]), True) #data
		plot_xs(bname, ax_xs, data_xax, mc, mc_err, 'saddlebrown', 4, m_alpha, process, "MC %s in %s"%(generator, bname), False)#mc
		plot_xs(bname, ax_xs, data_xax, fnloxs_LO, fnloLO_err, 'gray', 4, m_alpha, process, "LO Table %s in %s"%(pdfset, bname), False)#fnlo LO
		plot_xs(bname, ax_xs, data_xax, fnloxs_NLO, fnloLO_err, 'black', 4, m_alpha, process, "NLO Table %s in %s"%(pdfset, bname), False)#fnlo LO
		plot_xs(bname,ax_xs, data_xax, data, data_err, _ybysbin_col[bname], 6, m_alpha, process, "Data Run%s in %s"%(runperiod, _ybysbin_label[bname]), True) #data

		#ax_xs.set_title("%s XS in %s (Data Run%s, MC %s, FixedOrder Theory with %s)"%(process, bname, runperiod, generator, pdfset))
		ax_xs.set_title("%s XS overview (Data Run%s, MC %s, FixedOrder Theory with %s)"%(process, runperiod, generator, pdfset), fontsize=11)

		'''

		#make binbounds array (for flatten() later before using steppify())
		binbounds = np.array([data_low_bb, data_up_bb])

		#remember that in matplotlib errorbar plots all the given errors (also downwards) have to be absolute (i.e. positive) values
		#therefore: take np.absolute()


		#plot the error bands (mc around 1) --> should rather use plot_indiv_uncs() here
		#band just for total unc!
	#plotting total uncertainty:
	#	plot_unc_filled(bname, ax_doubleratio_bin, patches_bin, doubleratio_xax, binbounds, tot_down, tot_up, 'mediumseagreen', 'mediumseagreen', "zjet/dijet total unc (H7)", runperiod, generator)
		plot_unc_filled(bname, ax_ratio_bin, patches_bin, data_xax, binbounds, totunc_down, totunc_up, 'mediumseagreen', 'mediumseagreen', 0.6, "dijet total JES uncertainty (H7)", runperiod, generator)
		plot_unc_filled(bname, ax_ratio_bin, patches_bin, data_xax, binbounds, scaleunc_band_LO[0], scaleunc_band_LO[1], 'yellow', 'yellow', 0.1, "LO %s"%variation_type,runperiod,generator)

		plot_unc_filled(bname, ax_ratio_bin, patches_bin, data_xax, binbounds, scaleunc_band_NLO[0], scaleunc_band_NLO[1], 'tan', 'tan', 0.6, "NLO %s"%variation_type,runperiod,generator)



		#plot_unc_filled(bname, ax_ratio_bin, patches_bin, data_xax, binbounds, TotNoFlav_down, TotNoFlav_up, 'royalblue', 'royalblue', "jesTotalNoFlavor up (H7)", runperiod, generator)

		#lines -- TO DO --> adjust!! these are just copied from my doubleratio script
		plot_indiv_uncs(bname, ax_ratio_bin, patches_bin, data_xax, binbounds, TotNoFlav_up, 'royalblue', 'solid', "jesTotalNoFlavor up (H7)")
		plot_indiv_uncs(bname, ax_ratio_bin, patches_bin, data_xax, binbounds, TotNoFlav_down, 'royalblue', 'dashed', "jesTotalNoFlavor down (H7)")
		plot_indiv_uncs(bname, ax_ratio_bin, patches_bin, data_xax, binbounds, Flav_up, 'tomato', 'solid', "%s unc up (H7)"%flavname)		#new version
		plot_indiv_uncs(bname, ax_ratio_bin, patches_bin, data_xax, binbounds, Flav_down, 'tomato', 'dashed', "%s unc down (H7)"%flavname)	#new version
		#plot lumi:
		plot_indiv_uncs(bname, ax_ratio_bin, patches_bin, data_xax, binbounds, lumi_up, 'khaki', 'solid', "Luminosity (2.5 %)")		#lumi
		plot_indiv_uncs(bname, ax_ratio_bin, patches_bin, data_xax, binbounds, lumi_dw, 'khaki', 'dashed', "Luminosity (2.5 %)")	#lumi

		#should also plot Zjet?




		#plot the ratios (with statistical errors, and scale/pdf unc for theory)
		##plot_ratio(bname, ax_ratio_bin, data_xax, ratio_datamc, ratio_err, binbounds[0,:], binbounds[1,:], _ybysbin_col[bname], m_alpha, "data (Run%s) / MC (%s)"%(runperiod, generator))	#data/mc
		##plot_ratio(bname, ax_ratio_bin, data_xax, ratio_datafnloLO, ratio_fnloLO_err, binbounds[0,:], binbounds[1,:], 'gray', m_alpha, "data (Run%s) / LO theo. (%s)"%(runperiod, pdfset)) #data/LO
		##plot_ratio(bname, ax_ratio_bin, data_xax, ratio_datafnloNLO, ratio_fnloNLO_err, binbounds[0,:], binbounds[1,:], 'black', m_alpha, "data (Run%s) / NLO theo. (%s)"%(runperiod, pdfset)) #data/NLO

		#plot_ratio(bname, ax_ratio_bin, data_xax, ratio_datadata, np.absolute(ratio_data_err), binbounds[0,:], binbounds[1,:], _ybysbin_col[bname], m_alpha, "data (Run%s) stat. unc."%(runperiod))	#data/data  --> data stat. errors
		##ax_ratio_bin.errorbar(data_xax, ratio_datadata, xerr=np.array([np.subtract(data_xax,binbounds[0,:]),np.subtract(binbounds[1,:],data_xax)]), yerr=np.absolute(ratio_data_err), elinewidth=1.2, linewidth=1.0, marker='.', ms=2, color=_ybysbin_col[bname], fillstyle='none', fmt='.', label="data stat unc")
		ax_ratio_bin.errorbar(data_xax, ratio_datadata, xerr=None, yerr=np.absolute(ratio_data_err), elinewidth=1.2, linewidth=1.0, marker='.', ms=2, color=_ybysbin_col[bname], fillstyle='none', fmt='.', label="data (Run%s) stat. unc."%runperiod)


		plot_ratio(bname, ax_ratio_bin, data_xax, ratio_mcdata, np.absolute(ratio_mc_err), binbounds[0,:], binbounds[1,:], _ybysbin_col[bname], m_alpha, "MC (%s) / data (Run%s) with mc stat. unc."%(generator, runperiod))	#mc/data
		plot_ratio(bname, ax_ratio_bin, data_xax, ratio_fnloLOdata, np.absolute(ratio_fnloLO_pdfunc), binbounds[0,:], binbounds[1,:], 'gray', m_alpha, "LO theo. (%s) / data (Run%s) with PDF unc."%(pdfset, runperiod)) #LO/data
		plot_ratio(bname, ax_ratio_bin, data_xax, ratio_fnloNLOdata, np.absolute(ratio_fnloNLO_pdfunc), binbounds[0,:], binbounds[1,:], 'black', m_alpha, "NLO theo. (%s) / data (Run%s) with PDF unc."%(pdfset, runperiod)) #NLO/data





		#save the figures: should be six in total
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

		#add line at 1.0
		ax_ratio_bin.axhline(y=1.0, xmin=0, xmax=1, linewidth=0.6, color='grey', linestyle='dashed', dashes=(5,10))	#linewidth=0.4
		#ax_ratio_bin.grid(True, which="major", axis="y", color="k", linestyle="dotted", linewidth=0.8, alpha=0.2)
		ax_ratio_bin.grid(True, which="both", axis="y", color="k", linestyle="dotted", linewidth=0.8, alpha=0.2)
		ax_ratio_bin.grid(True, which="both", axis="x", color="k", linestyle="dotted", linewidth=0.8, alpha=0.2)

		#finalising (used to be done in plot_ratio() fct, but do it only once per bin, here)
		ax_ratio_bin.set_xlabel(r'$\mathrm{p_{T,avg}}$ /GeV', fontsize=12, horizontalalignment='right')
		ax_ratio_bin.xaxis.set_label_coords(1.00, -0.06)
		ax_ratio_bin.set_ylabel(r'Ratios', fontsize=12, horizontalalignment='center', verticalalignment='top', y=0.5, rotation=90, labelpad=24)

		#ax.set_title(r'%s in %s (%s)'%(process,bname, _ybysbin_label[bname]))
		#ax_ratio_bin.set_title("%s XS ratios in %s (Data Run%s, MC %s, FixedOrder Theory with %s)"%(process, _ybysbin_label[bname], runperiod, generator, pdfset), fontsize=11)
		ax_ratio_bin.set_title("%s XS ratios in %s (%s)"%(process, bname, _ybysbin_label[bname]), fontsize=12)
		ax_ratio_bin.legend(fontsize=10, numpoints=1, loc='upper left')

		#reordering of the legend entries
		ratiohandles, ratiolabels = ax_ratio_bin.get_legend_handles_labels()	#7 entries (totunc, scaleuncLO, scaleuncNLO, data, MC, LO, NLO)
		order_entries = [3,4,5,6,0,1,2]
		ax_ratio_bin.legend([ratiohandles[idx] for idx in order_entries], [ratiolabels[idx] for idx in order_entries], fontsize=10, numpoints=1, loc='upper left')
		ratio_legend = ax_ratio_bin.get_legend()	#get the "old" legend containing all the ratios
		ax_ratio_bin.add_artist(ratio_legend)		#draw the old legend
		unc_legend = plt.legend([("royalblue","dashed"),("tomato","dashed"),("khaki","dashed")], ["jesTotalNoFlavor unc. up, down (H7)",'%s unc. up, down (H7)'%flavname, r"Luminosity uncertainty ($\pm$ 2.5%)"], handler_map={tuple: AnyObjectHandler()}, loc='lower left')	#make the unc legend

		#plt.tick_params(axis='x', which='minor')				#?? not here but in main code
		#ax_ratio_bin.xaxis.set_minor_formatter(FormatStrFormatter("%.d"))	# same???
	

		xlocs = ax_xs.get_xticks(minor=True)
		xlabels = ax_xs.get_xticklabels(minor=True)
		print "xlocs: "
		print xlocs
		print "xlabels: "
		print xlabels

		#xlocs --> 48 minor label locations
		#xlabels --> 48 minor labels

		if(process=="dijet"):
			#x_minticks = [20, 30, 40, 50, 60, 70, 80, 90, 200, 300, 400, 500, 600, 700, 800, 900, 2000, 3000, 4000]
			#x_minticklabels = ['','30','','','60', '', '', '', '200', '', '', '500','', '700', '', '', '2000', '3000']
			x_minticks = [40, 50, 60, 70, 80, 90, 200, 300, 400, 500, 600, 700, 800, 900, 2000, 3000, 4000]
			x_minticklabels = ['','','60', '', '', '', '200', '', '', '500','', '700', '', '', '2000', '3000', '']
			x_majticks = [100, 1000]
			x_majticklabels = ['100', '1000']
			ax_ratio_bin.set_xticks(x_majticks, minor=False)	#set major ticks loc
			ax_ratio_bin.set_xticks(x_minticks, minor=True)		#set minor ticks loc
			ax_ratio_bin.set_xticklabels(x_majticklabels, minor=False)	#set major ticks label
			ax_ratio_bin.set_xticklabels(x_minticklabels, minor=True)	#set minor titcks label
		else:
			x_minticks = [40, 50, 60, 70, 80, 90, 200, 300, 400, 500, 600, 700, 800, 900]
			x_minticklabels = ['','','60', '', '', '', '200', '', '', '500','', '700', '', '']
			x_majticks = [100, 1000]
			x_majticklabels = ['100', '1000']
			ax_ratio_bin.set_xticks(x_majticks, minor=False)	#set major ticks loc
			ax_ratio_bin.set_xticks(x_minticks, minor=True)		#set minor ticks loc
			ax_ratio_bin.set_xticklabels(x_majticklabels, minor=False)	#set major ticks label
			ax_ratio_bin.set_xticklabels(x_minticklabels, minor=True)	#set minor titcks label


		#ax_ratio_bin.set_yticks(np.arange(0.0,2.0, step=0.1))
		y_minticks = [0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.10, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 1.90]
		y_minticklabels = ['', '0.20', '', '0.40', '', '0.60', '', '0.80', '', '', '1.20', '', '1.40', '', '1.60', '', '1.80', '']
		y_majticks = [0.00, 1.00, 2.00]
		y_majticklabels = ['0.00', '1.00', '2.00']
		ax_ratio_bin.set_yticks(y_majticks, minor=False)			#set major ticks loc
		ax_ratio_bin.set_yticks(y_minticks, minor=True)				#set minor ticks loc
		ax_ratio_bin.set_yticklabels(y_majticklabels, minor=False)	#set major ticks label
		ax_ratio_bin.set_yticklabels(y_minticklabels, minor=True)	#set minor titcks label

		for tick in ax_ratio_bin.yaxis.get_minor_ticks():
			tick.label.set_fontsize(8)

		#ax_ratio_bin.set_ylim(0.6, 1.4)	#test zoomed in
		ax_ratio_bin.set_ylim(0.0, 2.0)		#test


		plt.tight_layout()
		#fig_xs.savefig("%s_xs_overview_theory_Run%s_%s.png"%(process, runperiod, generator))
		#fig_ratio_bin.savefig("%s_xs_ratios_%s.png"%(process, bname))
		fig_ratio_bin.savefig("%s_xs_ratios_%s_Run%s_%s_%s.png"%(process, bname, runperiod, generator, pdfset))


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








