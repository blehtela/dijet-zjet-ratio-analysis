#!/usr/bin/env python2
#-*- coding:utf-8 -*-


#####################################################
#													#
#	Plotting the double ratio 						#
#	(Dijetdata/Zjetdata)/(DijetMC/ZjetMC)			#
#	Taking input from the ratio ROOT file			#
#													#
#	Created by B. Schillinger, 22.01.2020			#
#	Last modified: B.Schillinger, 26.02.2020		#
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

	print("[plots_DoubleRatio_Uncs.py]: Reading %s from %s" %(objname, rootfile))
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
def plot_doubleratio(bname, ax, ax_doubleratio_sum, x_axis, doubleratio, doubleratio_err, low_bb, up_bb, m_alpha, runperiod, generator):
##low_bb and up_bb are used for drawing the "x-error" == binwidth

	x_min= low_bb[0]
	x_max = up_bb[-1]
	print("xmin: %s, xmax: %s" %(x_min, x_max))
	'''
	y_min = min(ratio)
	y_max = max(ratio)
	print("ymin: %s, ymax: %s" %(y_min, y_max))
	'''

	#for x-error-bars
	xerr_low = np.subtract(x_axis, low_bb)
	xerr_up = np.subtract(up_bb, x_axis)
	#zip(xerr_low, xerr_up) #-->only assigns [(low1, up1), (low2, up2), ...]
	#print "xerr_low: \n", xerr_low
	#print "xerr_up: \n", xerr_up

	xerrors = np.array([xerr_low, xerr_up])


	#plot the double ratio data/MC (first zjet/dijet) errorbar
	#originally markersize was set to 10.
	ax.axhline(y=1.0, xmin=0, xmax=1, linewidth=0.4, color='grey', linestyle='dashed', dashes=(5,10))
	ax.errorbar(x_axis, doubleratio, xerr=xerrors, yerr=doubleratio_err, elinewidth=1.2, linewidth=1.0, marker=_ybysbin_marker[bname], ms=6, color=_ybysbin_col[bname], fillstyle='none', fmt='.', label=_ybysbin_label[bname])
	
	ax.grid(True, which="major", axis="y", color="k", linestyle="dotted", linewidth=0.8, alpha=0.2)
	ax.grid(True, which="both", axis="x", color="k", linestyle="dotted", linewidth=0.8, alpha=0.2)

	#does the following change anything at all if not applied to axes?
	#plt.yscale("log")
	#plt.xscale("log")

	#plt.grid(True, which="both", linestyle="dotted", color='k')


	#ax.set_xlim(40, 1200)
	#ax.set_ylim(y_min+0.0000001, y_max+0.1)

	ax.set_xlabel(r'$\mathrm{p_{T,avg}}$ /GeV', fontsize=12, horizontalalignment='right')
	ax.xaxis.set_label_coords(1.00, -0.06)
	#plt.xlabel(r'$\mathrm{p_{T,avg}}$ /GeV', horizontalalignment='right', loc='right')
	ax.set_ylabel(r'Double Ratio XS Data/MC', fontsize=12, horizontalalignment='center', verticalalignment='top', y=0.5, rotation=90, labelpad=24)

	if((runperiod!=None)and(generator!=None)):
		ax.set_title(r'Double Ratio data/MC (dijet/zjet) in %s (Run %s, %s)'%(bname,runperiod,generator))	#(first dijet/zjet)
	else:
		ax.set_title(r'Double Ratio in %s'%(bname))

	ax.legend(fontsize=14, numpoints=1, loc='upper left')

	'''
	if(runperiod!=None):
		ax.text(0.02, 1.01, "2018 data Run %s"%runperiod, fontsize=13, fontstyle='italic', ha='left', va='bottom', transform=ax.transAxes)          #text: Run period
	else:
		ax.text(0.02, 1.01, "2018 monte carlo", fontsize=13, fontstyle='italic', ha='left', va='bottom', transform=ax.transAxes)          			#text: monte carlo
	'''

	#add it to the summary plot
	ax_doubleratio_sum.errorbar(x_axis, doubleratio, xerr=0.0, yerr=doubleratio_err, elinewidth=1.2, linewidth=1.0, marker=_ybysbin_marker[bname], ms=6, color=_ybysbin_col[bname], fillstyle='none', fmt='.', label=_ybysbin_label[bname], alpha=m_alpha)

	plt.tick_params(axis='x', which='minor')
	ax.xaxis.set_minor_formatter(FormatStrFormatter("%.d"))
	#ax.xaxis.set_minor_formatter(FormatStrFormatter("%.1f"))



#function to plot selected unc source(s) as band
#for individual ybys bins
#here: used for total JEC uncertainty
def plot_unc_filled(bname, ax, patches, x_axis, binbounds, low_unc, up_unc, linecol, fillcol, uncname, runperiod, generator):

	ax.set_yscale("linear")
	#ax_ratio.set_ylim(-1.5, 2.0)	#test
	ax.set_ylim(0.0, 2.0)	#test
	#ax_ratio.get_legend().remove()

	#plot filled area
	ax.fill_between(binbounds.T.flatten(), steppify_bin(low_unc), steppify_bin(up_unc), edgecolor='grey', facecolor=fillcol, alpha=0.6)

	#patch for current uncertainty (needed for legend)
	#patches.append(matplotlib.patches.Rectangle((0,0), 0, 0, color=fillcol, label=uncname, alpha=0.6))
	patches.append(mpl.patches.Rectangle((0,0), 0, 0, color=fillcol, label=uncname, alpha=0.6))

	ax.add_patch(patches[-1])	#add latest created patch to list of patches

	#set xlim in a way that there is no white-space on the left of the lowest bin / right side of the highest bin
	##ax.set_xlim(binbounds[0,0], binbounds[1,-1])
	ax.set_xlim(30, 1200)	#test

	#xstart, xend = ax_xs.get_xlim()	#get the xrange of the xs plot, set the same for the ratio plot
	#ax_ratio.set_xlim((xstart, xend))

#function for plotting the individual unc-sources
#plot the up and down ones seperately (despite them being symmetrised)
def plot_indiv_uncs(bname, ax, patches, x_axis, binbounds, unc, linecol, linestyle, uncname, runperiod, generator):

	ax.step(binbounds[1], unc, color=linecol, linestyle=linestyle, where='pre', label=uncname, alpha=0.9)




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
	parser.add_argument('doubleratios_rootfile', type=TFile, nargs='?',
						help='Required argument! File with double ratio Data/MC (Dijet/Zjet for both) that shall be plotted. (filename glob)')
	parser.add_argument('uncsratios_rootfile', type=TFile, nargs='?',
						help='Required argument! File with uncertainty raitos that shall be plotted. (filename glob)')
	#parser.add_argument('process', type=str, nargs='?', choices=["zjet", "dijet"],
	#					help="Required argument! Specifies which analysis is considered: zjet for Z+jet analysis; dijet for Dijet analysis.") 

	# Optional arguments
	parser.add_argument('--alphamarker', '-a', type=float, default=1.0, nargs='?', action='store',
						help='If set to True: draw marker in summary plot with given input alpha.')

	parser.add_argument('--runperiod', '-r', type=str, nargs='?', default=None, action='store',
						help='Specifies which run period has been used: A, B, C, D')
	parser.add_argument('--generator', '-g', type=str, nargs='?', default=None, action='store',
						help='Specifies which generator has been used: herwig7, pythia8')

	parser.add_argument('--zjet-unc', dest='zjet_unc', type=TFile, nargs='?', default=None,
						help='TFile containing the uncertainty sources (symmetrised as up/down) for zjet.')
	parser.add_argument('--dijet-unc', dest='dijet_unc', type=TFile, nargs='?', default=None,
						help='TFile containing the uncertainty sources (symmetrised as up/down) for dijet.')
	parser.add_argument('--extrauncs', '-x', default=False, action='store_true',								#type=bool,
						help='If chosen, draw grey band for assumed uncorrelated uncertainties for Z+jet and Dijet.')


	#parser.add_argument('--alphamarker', '-a', default=False, action='store_true'
	#					help='If set to True: draw marker in summary plot with alpha=0.9.')

	if len(sys.argv)==1:
		parser.print_help(sys.stderr)
		sys.exit(1)

	# Parse arguments
	args = vars(parser.parse_args())

	# input filename
	ratio_rootfile = args['doubleratios_rootfile']
	uncs_rootfile = args['uncsratios_rootfile']
	#process = args['process']		#specifies whether looking at dijet or zjet
	m_alpha = args['alphamarker']
	runperiod = args['runperiod']	#optional argument, but good for labeling --> to be improved
	generator = args['generator']	#optional argument, pythia or herwig

	#optional
	zjetfile = args['zjet_unc']
	dijetfile = args['dijet_unc']

	#rootfile = TFile(args['rootfile'])
	print('\n')
	print("[plots_DoubleRatio_Uncs.py]: Taking DoubleRatio input from %s" %ratio_rootfile)
	print("[plots_DoubleRatio_Uncs.py]: Taking Uncertainty-Ratios input from %s" %uncs_rootfile)


	#prepare the lists (will later be arrays)
	bin_names = np.array(["yb0ys0", "yb0ys1", "yb0ys2", "yb1ys0", "yb1ys1", "yb2ys0"])


	#-------------------------------#
	#	preparing the summary plots	#
	#-------------------------------#
	#for the summary (all bins) of the doubleratio plots
	gs_doubleratio_sum = gridspec.GridSpec(3,3)
	fig_doubleratio_sum = plt.figure(figsize=(9,7), dpi=300)
	#ax1 = plt.subplot(gs[:-1,:])
	ax_doubleratio_sum = plt.subplot(xscale="log",yscale="linear")
	patches = []	#used for unc band



	#adjustments that are only done once for the summary plot
	#ax1_sum.axhline(y=1.0, xmin=0, xmax=1, linewidth=0.4, color=_ybysbin_col[bname], linestyle='dashed', dashes=(5,10))
	#plt.yscale("log")
	#plt.xscale("log")
	#plt.grid(True, which="major", linestyle="dotted", color='k')

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
	#loop through bins and create one plot for each ybys bin, containing the double ratio
	for j in range(0, 6):
		#current bin name:
		bname = bin_names[j]

		#for the individual (ybys bin) doubleratio plot
		fig_doubleratio_bin = plt.figure(figsize=(9,7), dpi=300)
		ax_doubleratio_bin = plt.subplot(xscale="log", yscale="linear")
		patches_bin = []

 
		#reading the ratio points in the current ybys bin from the rootfile
		#read (statistical?) errors on the histogram entries at the same time
		#read (from 2nd root file) JEC uncertainties
		#furthermore: read x-axis and lower and upper bin bounds
		doubleratio, doubleratio_err, doubleratio_xax, doubleratio_low_bb, doubleratio_up_bb = read_rootfile(ratio_rootfile, "Standard/%s/ratio_%s"%(bname, bname))
		#unc_TotalNoFlavor, _, unc_xax, unc_TotNoFlav_low_bb, unc_TotNoFlav_up_bb = read_rootfile(uncs_rootfile, "Standard/%s/mc_histo_%s"%(bname, bname))
		#unc_zjetZJetFlav_dijetQCDFlav, _, unc_xax, unc_Flav_low_bb, unc_Flav_up_bb = read_rootfile(uncs_rootfile, "")

		#read only the unc-ratio value (bin bounds etc. should be identical to the double ratio above)
		unc_TotalNoFlavor = read_rootfile(uncs_rootfile, "Standard/%s/ratio_TotalNoFlavor_%s"%(bname, bname))[0]
		#unc_zjetZJetFlav_dijetQCDFlav = read_rootfile(uncs_rootfile, "Standard/%s/zjetFlavZJet_dijetFlavQCD_%s"%(bname, bname))[0] #old version  with zjet/dijet
		unc_dijetQCDFlav_zjetZJetFlav = read_rootfile(uncs_rootfile, "Standard/%s/dijetFlavQCD_zjetFlavZJet_%s"%(bname, bname))[0]	#new version


		print("double ratios: \n %s"%doubleratio)


		#-----------------------#
		#		plotting		#
		#-----------------------#
		#creating the plots for the current ybys bin

		#create array with unc_up, and with unc_down (symmetrised errors)
		print("unc_TotalNoFlavor: ", unc_TotalNoFlavor)
		unc_size = unc_TotalNoFlavor.size	#size is an attribute of the array, not callable (i.e., no size() )
		ones_arr = np.full(unc_size, 1.0)

		
		#TOTAL NO FLAVOR
		delta_TotNoFlav = 0.5*np.absolute(np.add(ones_arr, (-1*unc_TotalNoFlavor)))	#delta_unc = 0.5*abs(unc-1)
		#preserve the sign: (do not use absolute value)
		delta_TotNoFlav_up = 0.5*(np.add(-1*ones_arr, unc_TotalNoFlavor))			#DeltaUnc upwards --> multiply it by -1 for downwards
		delta_TotNoFlav_down = -1.0*delta_TotNoFlav_up

		TotNoFlav_up = np.add(ones_arr, delta_TotNoFlav_up)		#1+DeltaUnc
		TotNoFlav_down = np.add(ones_arr, delta_TotNoFlav_down)	#1-DeltaUnc


		#FLAVOR
		delta_Flav = 0.5*np.absolute(np.add(ones_arr, (-1*unc_dijetQCDFlav_zjetZJetFlav)))
		delta_Flav_up = 0.5*(np.add(-1*ones_arr, unc_dijetQCDFlav_zjetZJetFlav))	#DeltaUnc upwards --> multiply it by -1 for downwards
		delta_Flav_down = -1.0*delta_Flav_up

		Flav_up = np.add(ones_arr, delta_Flav_up)				#1+DeltaUnc
		Flav_down = np.add(ones_arr, delta_Flav_down)			#1-DeltaUnc
	

		#calculate total uncertainty
		delta_tot = np.add(delta_TotNoFlav, delta_Flav) #--> MUSS QUADRATISCH ADDIERT WERDEN!!!! (siehe unten bei den extra uncs! da ist es richtig)

		tot_up = np.add(ones_arr, delta_tot)
		tot_down = np.add(ones_arr, -1.0*delta_tot)
		

		#make binbounds array (for flatten() later before using steppify())
		binbounds = np.array([doubleratio_low_bb, doubleratio_up_bb])



		#if extrauncs==true --> plot first the assumption of "uncorrelated dijet and zjet uncs" (grey band)
		if(args['extrauncs']==True):
			#read the files for Z+jet
			zjet_unc_TotalNoFlavor = read_rootfile(zjetfile, "UncertaintyRatios/%s/UpDownDiv_TotalNoFlavor_%s"%(bname, bname))[0]	
			zjet_unc_FlavorZJet = read_rootfile(zjetfile, "UncertaintyRatios/%s/UpDownDiv_FlavorZJet_%s"%(bname, bname))[0]	

			#Zjet uncs
			delta_zjet_TotNoFlav_down = 0.5*(np.add(-1*ones_arr, zjet_unc_TotalNoFlavor))
			delta_zjet_TotNoFlav_up	= -1.0*delta_zjet_TotNoFlav_down
			#delta_zjet_FlavZJet = 0.5*np.absolute(np.add(ones_arr, (-1*zjet_unc_FlavorZJet)))	#absolute value
			delta_zjet_FlavZJet_down = 0.5*(np.add(-1*ones_arr, zjet_unc_FlavorZJet))
			delta_zjet_FlavZJet_up	= -1.0*delta_zjet_FlavZJet_down
	
			#plus 1.0 to draw it "around MC" (which is denominator)
			zjet_TotNoFlav_up = np.add(ones_arr, delta_zjet_TotNoFlav_up)		#1+DeltaUnc
			zjet_TotNoFlav_down = np.add(ones_arr, delta_zjet_TotNoFlav_down)	#1-DeltaUnc
			zjet_FlavZJet_up = np.add(ones_arr, delta_zjet_FlavZJet_up)
			zjet_FlavZJet_down = np.add(ones_arr, delta_zjet_FlavZJet_down)

			#---------------
			#read the files for Dijet
			dijet_unc_TotalNoFlavor = read_rootfile(dijetfile, "UncertaintyRatios/%s/UpDownDiv_TotalNoFlavor_%s"%(bname, bname))[0]
			dijet_unc_FlavorQCD = read_rootfile(dijetfile, "UncertaintyRatios/%s/UpDownDiv_FlavorQCD_%s"%(bname, bname))[0]

			#Dijet uncs
			delta_dijet_TotNoFlav_down = 0.5*(np.add(-1*ones_arr, dijet_unc_TotalNoFlavor))
			delta_dijet_TotNoFlav_up = -1.0*delta_dijet_TotNoFlav_down
			delta_dijet_FlavQCD_down = 0.5*(np.add(-1*ones_arr, dijet_unc_FlavorQCD))
			delta_dijet_FlavQCD_up	= -1.0*delta_dijet_FlavQCD_down
	
			#plus 1.0 to draw it "around MC" (which is denominator)
			dijet_TotNoFlav_up = np.add(ones_arr, delta_dijet_TotNoFlav_up)		#1+DeltaUnc
			dijet_TotNoFlav_down = np.add(ones_arr, delta_dijet_TotNoFlav_down)	#1-DeltaUnc
			dijet_FlavQCD_up = np.add(ones_arr, delta_dijet_FlavQCD_up)
			dijet_FlavQCD_down = np.add(ones_arr, delta_dijet_FlavQCD_down)

			#-------------------
			# total uncorr unc
			#-------------------

			tot_zjet_up = np.add(delta_zjet_TotNoFlav_up, delta_zjet_FlavZJet_up)		#total zjet uncertainty up
			tot_zjet_down = np.add(delta_zjet_TotNoFlav_down, delta_zjet_FlavZJet_down) #total zjet uncertainty down

			tot_dijet_up = np.add(delta_dijet_TotNoFlav_up, delta_dijet_FlavQCD_up)			#total dijet uncertainty up
			tot_dijet_down = np.add(delta_dijet_TotNoFlav_down, delta_dijet_FlavQCD_down)	#total dijet uncertainty down

			delta_tot_uncorr_up = np.sqrt(np.add(tot_zjet_up*tot_zjet_up, tot_dijet_up*tot_dijet_up))		#total uncorr unc up = Sqrt[(zjetup^2)+(dijetup^2)]
			delta_tot_uncorr_down = np.sqrt(np.add(tot_zjet_down*tot_zjet_down, tot_dijet_down*tot_dijet_down))	#total uncorr unc down (should be the same as symmetrised...)

			#plotted around 1
			tot_uncorr_up = np.add(ones_arr, delta_tot_uncorr_up)				#1.0+tot_uncorr
			tot_uncorr_down = np.add(ones_arr, -1.0*delta_tot_uncorr_down)			#1.0-tot_uncorr

			'''
			delta_tot_zjet = np.add(delta_zjet_TotNoFlav, delta_zjet_FlavZJet)
			delta_tot_dijet = np.add(delta_dijet_TotNoFlav, delta_dijet_FlavQCD)
			delta_uncorr_tot = np.sqrt(
			'''

			#------------------------------------
			# plot the grey uncorr unc band first (below everything)
			#------------------------------------
			plot_unc_filled(bname, ax_doubleratio_bin, patches_bin, doubleratio_xax, binbounds, tot_uncorr_down, tot_uncorr_up, 'powderblue', 'powderblue', "Quadratic Sum (zjet, dijet)", runperiod, generator)



		#plotting total uncertainty:
		plot_unc_filled(bname, ax_doubleratio_bin, patches_bin, doubleratio_xax, binbounds, tot_down, tot_up, 'mediumseagreen', 'mediumseagreen', "zjet/dijet total unc (H7)", runperiod, generator)


		#plotting TotalNoFlavor uncertainty band
		#plot_unc_filled(bname, ax_doubleratio_bin, patches_bin, doubleratio_xax, binbounds, TotNoFlav_down, TotNoFlav_up, 'royalblue', 'royalblue', "zjet/dijet jesTotalNoFlavor (H7)", runperiod, generator)
		plot_indiv_uncs(bname, ax_doubleratio_bin, patches_bin, doubleratio_xax, binbounds, TotNoFlav_up, 'royalblue', 'solid', "jesTotalNoFlavor up (H7)", runperiod, generator)
		plot_indiv_uncs(bname, ax_doubleratio_bin, patches_bin, doubleratio_xax, binbounds, TotNoFlav_down, 'royalblue', 'dashed', "jesTotalNoFlavor down (H7)", runperiod, generator)


		#plotting zjetZJetFlavor_dijetQCDFlavor uncertainty band
		#plot_unc_filled(bname, ax_doubleratio_bin, patches_bin, doubleratio_xax, binbounds, Flav_down, Flav_up, 'tomato', 'tomato', "zjetZJetFlav unc / dijetQCDFlav unc (H7)", runperiod, generator)
		##plot_indiv_uncs(bname, ax_doubleratio_bin, patches_bin, doubleratio_xax, binbounds, Flav_up, 'tomato', 'solid', "zjetZJetFlav unc / dijetQCDFlav unc up (H7)", runperiod, generator) 		#old version
		##plot_indiv_uncs(bname, ax_doubleratio_bin, patches_bin, doubleratio_xax, binbounds, Flav_down, 'tomato', 'dashed', "zjetZJetFlav unc / dijetQCDFlav unc down (H7)", runperiod, generator)	#old version
		plot_indiv_uncs(bname, ax_doubleratio_bin, patches_bin, doubleratio_xax, binbounds, Flav_up, 'tomato', 'solid', "dijetQCDFlav unc / zjetZJetFlav unc up (H7)", runperiod, generator) 		#new version
		plot_indiv_uncs(bname, ax_doubleratio_bin, patches_bin, doubleratio_xax, binbounds, Flav_down, 'tomato', 'dashed', "dijetQCDFlav unc / zjetZJetFlav unc down (H7)", runperiod, generator)	#new version


	

		#plotting the doubleratio
		plot_doubleratio(bname, ax_doubleratio_bin, ax_doubleratio_sum, doubleratio_xax, doubleratio, doubleratio_err, doubleratio_low_bb, doubleratio_up_bb, m_alpha, runperiod, generator)



		'''
		#try simple error propagation for absolute error of ratio:
		# ratio_err = Ratio*Sqrt((data_err/data)^2 + (mc_err/mc)^2)
		# calculate relative errors:
		data_relerr = np.divide(data_err, data)	# relative error
		mc_relerr = np.divide(mc_err, mc)		# relative error

		#calculate ratio relative error:
		ratio_relerr = np.sqrt(pow(data_relerr,2)+pow(mc_relerr,2))	#numpy.sqrt() returns sqrt for array elementwise

		#calulate ratio absolute error:
		ratio_err = np.multiply(cur_ratio, ratio_relerr)	#elementwise numpy multiplication of arrays
		'''

		#print some axis info for testing
		#print(ax_doubleratio_bin.xaxis.get_minor_ticks()[-2].label1)
		#ax_doubleratio_bin.xaxis.get_minor_ticks()[-2].label1.set_visible(False)	#test to remove a label


		#-------------------#
		# legend and save 	#
		#-------------------#
		#save plots for current ybys bin
		#doubleratio plot:
		#plt.legend()
		ax_doubleratio_bin.legend(fontsize=10, numpoints=1, loc='lower left')
		if((runperiod!=None)and(generator!=None)):
			fig_doubleratio_bin.savefig("doubleratio_DijetZjet_%s_data_Run%s_mc_%s.png"%(bname,runperiod,generator))
		else:
			fig_doubleratio_bin.savefig("doubleratio_%s_data_mc.png"%(bname))

	
	#---------------------------------------#
	#	finalising the summary plots		#
	#---------------------------------------#

	'''
	#RATIO SUMMARY
	ax_ratio_sum.set_ylabel(r'Ratio Data/MC', fontsize=12, horizontalalignment='center', verticalalignment='top', y=0.5, rotation=90, labelpad=24)
	ystart, yend = ax_ratio_sum.get_ylim()
	ystart = np.floor(ystart)		#round to largest integer lower than original ystart
	yend = np.ceil(yend)			#round to lowest integer larger than original yend
	ax_ratio_sum.yaxis.set_ticks(np.arange(ystart, yend, 1.0))
	'''	
	

	# DOUBLERATIO SECTION SUMMARY
	#adjustments that are only done once for the summary plot
	ax_doubleratio_sum.set_xlim(30, 1200)
	ax_doubleratio_sum.set_xlabel(r'$\mathrm{p_{T,avg}}$ /GeV', fontsize=12, horizontalalignment='right')
	ax_doubleratio_sum.xaxis.set_label_coords(1.00, -0.06)
	#plt.xlabel(r'$\mathrm{p_{T,avg}}$ /GeV', horizontalalignment='right', loc='right')
	ax_doubleratio_sum.set_ylabel(r'XS ratio Data/MC (first dijet/zjet)', fontsize=12, horizontalalignment='center', verticalalignment='top', y=0.5, rotation=90, labelpad=24)
	if((runperiod!=None)and(generator!=None)):
		ax_doubleratio_sum.set_title(r'Double Ratio Summary (Run %s, %s)'%(runperiod,generator), fontsize=14, loc='right', horizontalalignment='right')
	else:
		ax_doubleratio_sum.set_title(r'Double Ratio Summary', fontsize=14, loc='right', horizontalalignment='right')

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
	ax_doubleratio_sum.legend(fontsize=14, numpoints=1, loc='lower left')

	plt.tight_layout()		#especially for the combi plot

	if((generator!=None)and(runperiod!=None)):
		fig_doubleratio_sum.savefig("doubleratio_DijetZjet_overview_mc_%s_data_Run%s.png"%(generator,runperiod))
	else:
		fig_doubleratio_sum.savefig("doubleratio_overview_mc_data.png")


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
