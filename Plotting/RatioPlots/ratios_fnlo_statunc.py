#!/usr/bin/env python2
#-*- coding:utf-8 -*-

#############################################################
#														   	#
# Plotting the Dijet / Z+jet ratios							#
# Taking input from fastNLO tables							#
# calculating the ratio!									#
#															#
#	This is an updated version of ratios_root.py			#															
#	It contains the following changes:						#
#	> Only display the ratio-range available in data		#
#	> Changed axis ranges accordingly						#
#	> Plus some smaller adjustments							#
#															#
#	> added statistical uncertainty							#
#	> removed  PDF and scale unc							#
#															#
#	--> Removed unused and commented-out lines				#
#	--> Only create a summary plot							#
#													   		#
# Created by B. Schillinger, 26.02.2020				 		#
# Last modified by B. Schillinger, 12.04.2020		   		#
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

import fastnlo
import ROOT
from ROOT import TFile, TH1D, TH1F
from ROOT import gROOT

#dictionary --> ybys bin to color
_ybysbin_col = {'yb0ys0':'forestgreen', 'yb0ys1':'mediumblue', 'yb0ys2':'orange', 'yb1ys0':'firebrick', 'yb1ys1':'deepskyblue', 'yb2ys0':'mediumpurple'}
_ybysbin_marker = {'yb0ys0':"o", 'yb0ys1':"^", 'yb0ys2':"s", 'yb1ys0':"d", 'yb1ys1':"P", 'yb2ys0':"v"}
_ybysbin_label = {'yb0ys0':r'$0 \leq y_b < 1$	$0 \leq y^{\ast} < 1$', 'yb0ys1':r'$0 \leq y_b < 1$	$1 \leq y^{\ast} < 2$', 'yb0ys2':r'$0 \leq y_b < 1$	$2 \leq y^{\ast} < 2.4$', 'yb1ys0':r'$1 \leq y_b < 2$	$0 \leq y^{\ast} < 1$', 'yb1ys1':r'$1 \leq y_b < 2$	$1 \leq y^{\ast} < 2$', 'yb2ys0':r'$2 \leq y_b < 2.4$	$0 \leq y^{\ast} < 1$'}

_ybysbin_xind_up_zjetdata = {'yb0ys0':40, 'yb0ys1':30, 'yb0ys2':20, 'yb1ys0':37, 'yb1ys1':30, 'yb2ys0':26}	#constraints coming from rootfiles (zjet data) --> applied on rootfiles


_ybysbin_xind_up_zjetfnlo = {'yb0ys0':-22, 'yb0ys1':-24, 'yb0ys2':-20, 'yb1ys0':-13, 'yb1ys1':-13, 'yb2ys0':-11}	#contraints coming from rootfiles (zjet data) --> applied on fnlo [bins referring to fnlo binning]


#----------------------------------------------------------------------------------------------------------------#
#												  function definitions											 # 
#----------------------------------------------------------------------------------------------------------------#


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
	

	#get the x-axis
	binbounds = np.array(fnlo.GetObsBinsBounds(0))
	x_axis = (binbounds.T[0]+binbounds.T[1])/2.0



	#return xs_array, entries_err, x_axis	#x_axis should be the same as for the root hists
	return xs_array, abs_scale_unc, abs_pdf_unc, binbounds


#function to read XS from fastnlo table
def read_fnlotable_simple(fnlotable, pdfset, order):
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

	#get the x-axis
	binbounds = np.array(fnlo.GetObsBinsBounds(0))
	x_axis = (binbounds.T[0]+binbounds.T[1])/2.0


	#return xs_array, entries_err, x_axis	#x_axis should be the same as for the root hists
	return xs_array, binbounds




def plot_ratio(bname, ax, ax_sum, x_axis, ratio, binerr, low_bb, up_bb, m_alpha, pdfset): ##low_bb and up_bb are used for drawing the "x-error" == binwidth

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

	print "in plot_ratio() function:"
	print "x_axis: ", x_axis
	print "low_bb: ", low_bb
	print "up_bb: ", up_bb

	print "TEST:  xerrors: ", xerrors

	ax.axhline(y=1.0, xmin=0, xmax=1, linewidth=0.4, color=_ybysbin_col[bname], linestyle='dashed', dashes=(5,10))
	ax.errorbar(x_axis, ratio, xerr=xerrors, yerr=binerr, elinewidth=1.2, linewidth=0.0, marker=_ybysbin_marker[bname], ms=12, color=_ybysbin_col[bname], fillstyle='none', fmt='.', label=_ybysbin_label[bname])	#label=bname
	plt.yscale("log")
	plt.xscale("log")

	plt.grid(True, which="major", linestyle="dotted", color='k')

	ax.set_xlim(40, 1200)
	#ax.set_ylim(y_min+0.0000001, y_max+0.1)
	ax.set_ylim(1000, 1000000)

	ax.set_xlabel(r'$\mathrm{p_{T,avg}}$ /GeV', fontsize=13, horizontalalignment='right')
	ax.xaxis.set_label_coords(1.00, -0.06)
	#plt.xlabel(r'$\mathrm{p_{T,avg}}$ /GeV', horizontalalignment='right', loc='right')
	ax.set_ylabel(r'XS ratios: Dijet / Zjet', fontsize=13, horizontalalignment='center', verticalalignment='top', y=0.5, rotation=90, labelpad=24)
	ax.set_title(r'%s'%bname)

	ax.legend(fontsize=13, numpoints=1, loc='upper right')
	ax.text(0.02, 1.01, "FixedOrder Theory Dijet/Z+jet (%s)"%pdfset, fontsize=13, fontstyle='italic', ha='left', va='bottom', transform=ax.transAxes)          #text: pdf


	#add it to the total plot
	ax_sum.errorbar(x_axis, ratio, xerr=xerrors, yerr=binerr, elinewidth=1.2, linewidth=0.0, marker=_ybysbin_marker[bname], ms=8, color=_ybysbin_col[bname], fillstyle='none', fmt='.', label=_ybysbin_label[bname], alpha=m_alpha)




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
	parser.add_argument('fnlo_tables_dijet', type=str, nargs=6,	#nargs='+'
						help='Required argument! FastNLO tables for DIJET for the six ybys bins to be evaluated.')
	parser.add_argument('fnlo_tables_zjet', type=str, nargs=6,	#nargs='+'
						help='Required argument! FastNLO tables for Z+JET for the six ybys bins to be evaluated.')

	parser.add_argument('datfiles_dijet', type=str, nargs=6,	#nargs='+'
						help='Required argument! Six .dat-files for DIJET for the six ybys bins, to evaluate the stat.uncs..')
	parser.add_argument('datfiles_zjet', type=str, nargs=6,	#nargs='+'
						help='Required argument! Six .dat-files for Z+JET for the six ybys bins, to evaluate the stat.uncs..')



	# Optional arguments
	parser.add_argument('--alphamarker', '-a', type=float, default=1.0, nargs='?', action='store',
						help='If set to True: draw marker in summary plot with given input alpha.')
	parser.add_argument('--pdfset', '-p', default='CT14nlo', nargs='?',
						help='PDF set for fnlo table evaluation. Default is CT14nlo.')

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
	fnlo_tables_dijet = args['fnlo_tables_dijet']	#this is a list
	fnlo_tables_zjet = args['fnlo_tables_zjet']		#this is a list
	datfiles_dijet = args['datfiles_dijet']			#this is a list
	datfiles_zjet = args['datfiles_zjet']			#this is a list

	ybys = args['ybysbin']				#name of ybys bin
	m_alpha = args['alphamarker']
	pdfset = args['pdfset']
	verbose = args['verbose']
	
	#rootfile = TFile(args['rootfile'])
	print('\n')
	print("[ratios_fnlo.py]: Taking dijet theory input from %s" %fnlo_tables_dijet)
	print("[ratios_fnlo.py]: Taking zjet theory input from %s" %fnlo_tables_zjet)


	#prepare the lists (will later be arrays)
	if(ybys==None):
		bin_names = np.array(["yb0ys0", "yb0ys1", "yb0ys2", "yb1ys0", "yb1ys1", "yb2ys0"])
	else:
		bin_names = np.array([ybys])



	#-------------------------------------------------------#
	#	preparing the plot for overview of all 6 ybys bin	#
	#-------------------------------------------------------#
	#for the summary (all bins) of the XS plots
	#Next-to-leading order
	gs_ratio_NLO = gridspec.GridSpec(3,3)
	fig_ratio_NLO = plt.figure(figsize=(8,7), dpi=300)
	ax_ratio_NLO = plt.subplot(xscale="log", yscale="log")

	#looping through the ybys bins
	#creating one ratio-plot per bin
	for i in range(bin_names.size): #should usually be 6
		bname = bin_names[i]
		if(verbose==True): print "[ratios_fnlo.py]: Now in bin ", bname

		#for the individual (ybys bin) ratio plot
		#fig_ratio_bin_LO = plt.figure(figsize=(9,7), dpi=300)
		#ax_ratio_bin_LO = plt.subplot(xscale="log", yscale="linear")
		fig_ratio_bin_NLO = plt.figure(figsize=(9,7), dpi=300)
		ax_ratio_bin_NLO = plt.subplot(xscale="log", yscale="linear")
		#patches_bin = []

		#fnloxs = read_fnlotable(fnlo_table)
		if(verbose==True): print "[ratios_fnlo.py]: Reading fnlo table: ", fnlo_tables_dijet[i]
		#read the tables
		dijet_fnloxs_NLO, dijet_binbounds_NLO = read_fnlotable_simple(fnlo_tables_dijet[i], pdfset, 1)
		zjet_fnloxs_NLO, zjet_binbounds_NLO = read_fnlotable_simple(fnlo_tables_zjet[i], pdfset, 1)

		#reat the datfiles
		dijet_dat_statunc = np.loadtxt(datfiles_dijet[i], usecols=4) 	#absolute scale uncertainty
		dijet_dat_xs = np.loadtxt(datfiles_dijet[i], usecols=3) 		#XS from datfile

		zjet_dat_statunc = np.loadtxt(datfiles_zjet[i], usecols=4)		#absolute scale uncertainty
		zjet_dat_xs = np.loadtxt(datfiles_zjet[i], usecols=3)			#XS from datfile

		#relative stat. unc. for individual analyses
		dijet_rel_statunc = np.divide(dijet_dat_statunc, dijet_dat_xs)
		zjet_rel_statunc = np.divide(zjet_dat_statunc, zjet_dat_xs)


		nentries = dijet_fnloxs_NLO.size

		

		#-------------------------------------------#
		# 	Calculate the ratios and their errors	#
		#-------------------------------------------#
		# dijet/zjet in LO and in NLO	
		ratio_NLO = np.divide(dijet_fnloxs_NLO,zjet_fnloxs_NLO)

		#calculate ratio relative error:
		ratio_NLO_relerr = np.sqrt(pow(dijet_rel_statunc,2)+pow(zjet_rel_statunc,2))	#end up with abs() values
		sliced_relerr = ratio_NLO_relerr[4:_ybysbin_xind_up_zjetfnlo[bname]]

		print "Relative statistical uncertainty on the ratio in %s: "%bname
		print ratio_NLO_relerr
		print "Sliced rel. stat. error: "
		print sliced_relerr

		#calculate ratio absolute error: (is symmetric)
		ratio_NLO_err = np.zeros([2, nentries])
		ratio_NLO_err[0] = np.multiply(ratio_NLO, ratio_NLO_relerr)
		ratio_NLO_err[1] = np.multiply(ratio_NLO, ratio_NLO_relerr)


		#x_axis = (binbounds.T[0]+binbounds.T[1])/2.0
		x_axis = (dijet_binbounds_NLO.T[0]+dijet_binbounds_NLO.T[1])/2.0

		#------------------------ SLICING ------------------------------#
		#NEW: Do some slicing to bring it to data-range! (50 to 1200 GeV)
		#restrictions coming from zjet-data, applied on fnlo
	
		print "dijet_binbounds_NLO before slicing", dijet_binbounds_NLO
	
		x_axis = x_axis[4:_ybysbin_xind_up_zjetfnlo[bname]]
		dijet_bb_NLO_trans = dijet_binbounds_NLO.T[:,4:_ybysbin_xind_up_zjetfnlo[bname]] 	#work-around
		dijet_binbounds_NLO = dijet_bb_NLO_trans.T											#get back to old bb-var

		print "dijet_binbounds_NLO after slicing", dijet_binbounds_NLO

		ratio_NLO = ratio_NLO[4:_ybysbin_xind_up_zjetfnlo[bname]]
		ratio_NLO_err = ratio_NLO_err[:,4:_ybysbin_xind_up_zjetfnlo[bname]]
		#--------------------END: SLICING-------------------------------#

		#plotting NLO
		plot_ratio(bname, ax_ratio_bin_NLO, ax_ratio_NLO, x_axis, ratio_NLO, ratio_NLO_err, dijet_binbounds_NLO.T[0], dijet_binbounds_NLO.T[1], m_alpha, pdfset)

	#---------------------------#
	#	legend, title and save	#
	#---------------------------#
	#ax_ratio.grid(True, which="both", axis="y", color="k", linestyle="dotted", linewidth=0.8, alpha=0.2)
	#Changed the ticks! --> for old version see ratios_fnlo.py
	x_minticks = [50, 60, 70, 80, 90, 200, 300, 400, 500, 600, 700, 800, 900, 2000]
	x_minticklabels = ['','60', '', '', '', '200', '', '', '500','', '700', '', '', '2000']
	x_majticks = [100, 1000]
	x_majticklabels = ['100', '1000']
	

	#---------------------------#
	#	Next-to-Leading Order	#
	#---------------------------#
	ax_ratio_NLO.grid(True, which="minor", axis="y", color="k", linestyle="dotted", linewidth=0.8, alpha=0.2)
	ax_ratio_NLO.grid(True, which="major", axis="y", color="k", linestyle="dotted", linewidth=0.8, alpha=0.6)
	ax_ratio_NLO.grid(True, which="both", axis="x", color="k", linestyle="dotted", linewidth=0.8, alpha=0.2)
	ax_ratio_NLO.set_xlabel(r'$\mathrm{p_{T,avg}}$ /GeV', fontsize=13, horizontalalignment='right')
	ax_ratio_NLO.xaxis.set_label_coords(1.00, -0.06)
	ax_ratio_NLO.set_ylabel(r'XS ratios: Dijet / Zjet', fontsize=13, horizontalalignment='center', verticalalignment='top', y=0.5, rotation=90, labelpad=24)

	ax_ratio_NLO.set_xticks(x_majticks, minor=False)	#set major ticks loc
	ax_ratio_NLO.set_xticks(x_minticks, minor=True)		#set minor ticks loc
	ax_ratio_NLO.set_xticklabels(x_majticklabels, minor=False, fontsize=12)	#set major ticks label
	ax_ratio_NLO.set_xticklabels(x_minticklabels, minor=True, fontsize=12)	#set minor titcks label

	#ytick labels
	for tick in ax_ratio_NLO.yaxis.get_major_ticks():
		tick.label.set_fontsize(12)

	#test for suitable y-axis range
	ax_ratio_NLO.set_ylim(1e3, 4*1e6)
	
	ax_ratio_NLO.legend(fontsize=12, numpoints=1, loc='upper right')
	ax_ratio_NLO.set_title("Cross Section Ratios: Dijet over Z+jet \n Fixed-Order Theory Calculations in NLO (%s)"%pdfset, fontsize=14)
	plt.tight_layout()
	fig_ratio_NLO.savefig("ratio_summary_fnlo_NLO_statuncs.png")

	#test with even smaller x-range
	ax_ratio_NLO.set_xlim(50,1200)
	plt.tight_layout()
	fig_ratio_NLO.savefig("ratio_summary_fnlo_NLO_smallxrange_50to1200_statuncs.png")


	# stop timer
	stop_time = timeit.default_timer()
	timediff = stop_time-start_time
	print("[ratios_fnlo.py]: Elapsed time: %s sec = %s min" %(timediff, round(timediff/60., 2)))



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








