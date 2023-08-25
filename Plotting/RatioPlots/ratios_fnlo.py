#!/usr/bin/env python2
#-*- coding:utf-8 -*-

#############################################################
#														   	#
# Plotting the Dijet / Z+jet ratios							#
# Taking input from fastNLO tables							#
# calculating the ratio!									#
#													   		#
# Created by B. Schillinger, 26.02.2020				 		#
# Last modified by B. Schillinger, 26.02.2020		   		#
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

	ax.axhline(y=1.0, xmin=0, xmax=1, linewidth=0.4, color=_ybysbin_col[bname], linestyle='dashed', dashes=(5,10))
	ax.errorbar(x_axis, ratio, xerr=xerrors, yerr=binerr, elinewidth=1.2, linewidth=0.0, marker=_ybysbin_marker[bname], ms=12, color=_ybysbin_col[bname], fillstyle='none', fmt='.', label=_ybysbin_label[bname])	#label=bname
	plt.yscale("log")
	plt.xscale("log")

	plt.grid(True, which="major", linestyle="dotted", color='k')

	ax.set_xlim(40, 1200)
	#ax.set_ylim(y_min+0.0000001, y_max+0.1)
	ax.set_ylim(1000, 1000000)

	ax.set_xlabel(r'$\mathrm{p_{T,avg}}$ /GeV', fontsize=12, horizontalalignment='right')
	ax.xaxis.set_label_coords(1.00, -0.06)
	#plt.xlabel(r'$\mathrm{p_{T,avg}}$ /GeV', horizontalalignment='right', loc='right')
	ax.set_ylabel(r'XS ratios: Dijet / Zjet', fontsize=12, horizontalalignment='center', verticalalignment='top', y=0.5, rotation=90, labelpad=24)
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

	# Optional arguments
	parser.add_argument('--alphamarker', '-a', type=float, default=1.0, nargs='?', action='store',
						help='If set to True: draw marker in summary plot with given input alpha.')
	parser.add_argument('--pdfset', '-p', default='CT14nlo', nargs='?',
						help='PDF set for fnlo table evaluation. Default is CT14nlo.')
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
	fnlo_tables_dijet = args['fnlo_tables_dijet']	#this is a list
	fnlo_tables_zjet = args['fnlo_tables_zjet']		#this is a list
	ybys = args['ybysbin']				#name of ybys bin
	m_alpha = args['alphamarker']
	pdfset = args['pdfset']
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
	#Leading order
	gs_ratio_LO = gridspec.GridSpec(3,3)
	fig_ratio_LO = plt.figure(figsize=(8,7), dpi=300)
	#ax1 = plt.subplot(gs[:-1,:])
	ax_ratio_LO = plt.subplot(xscale="log", yscale="log")
	#patches = []

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
		fig_ratio_bin_LO = plt.figure(figsize=(9,7), dpi=300)
		ax_ratio_bin_LO = plt.subplot(xscale="log", yscale="linear")
		fig_ratio_bin_NLO = plt.figure(figsize=(9,7), dpi=300)
		ax_ratio_bin_NLO = plt.subplot(xscale="log", yscale="linear")
		#patches_bin = []

		#fnloxs = read_fnlotable(fnlo_table)
		if(verbose==True): print "[ratios_fnlo.py]: Reading fnlo table: ", fnlo_tables_dijet[i]
		dijet_fnloxs_LO, dijet_scale_unc_LO, dijet_pdf_unc_LO, dijet_binbounds_LO = read_fnlotable(fnlo_tables_dijet[i], pdfset, 0, scale_var_type)
		dijet_fnloxs_NLO, dijet_scale_unc_NLO, dijet_pdf_unc_NLO, dijet_binbounds_NLO = read_fnlotable(fnlo_tables_dijet[i], pdfset, 1, scale_var_type)
		zjet_fnloxs_LO, zjet_scale_unc_LO, zjet_pdf_unc_LO, zjet_binbounds_LO = read_fnlotable(fnlo_tables_zjet[i], pdfset, 0, scale_var_type)
		zjet_fnloxs_NLO, zjet_scale_unc_NLO, zjet_pdf_unc_NLO, zjet_binbounds_NLO = read_fnlotable(fnlo_tables_zjet[i], pdfset, 1, scale_var_type)

		nentries = dijet_fnloxs_NLO.size

		dijet_fnloLO_err = np.full(nentries,0)
		dijet_fnloNLO_err = np.full(nentries,0)
		zjet_fnloLO_err = np.full(nentries,0)
		zjet_fnloNLO_err = np.full(nentries,0)

		#check the absolute unc by giving output
		if(verbose==True):
			print "Scale unc dijet NLO down: \n", dijet_scale_unc_NLO[0]
			print "Scale unc dijet NLO up: \n", dijet_scale_unc_NLO[1]

		#calculate error using scale unc and pdf unc (quadratic sum)
		#fnlo_relerrNLO = np.sqrt(pow(rel_scale_unc_NLO, 2), pow(rel_pdf_unc_NLO, 2))
		#fnlo_errNLO = np.multiply(fnloxs_NLO, fnlo_relerrNLO)
		dijet_fnlo_errLO_down = np.sqrt(pow(dijet_scale_unc_LO[0,:],2)+pow(dijet_pdf_unc_LO[0,:],2))	#LO down
		dijet_fnlo_errLO_up = np.sqrt(pow(dijet_scale_unc_LO[1,:],2)+pow(dijet_pdf_unc_LO[1,:],2))		#LO up
		dijet_fnlo_errNLO_down = np.sqrt(pow(dijet_scale_unc_NLO[0,:],2)+pow(dijet_pdf_unc_NLO[0,:],2))	#NLO down
		dijet_fnlo_errNLO_up = np.sqrt(pow(dijet_scale_unc_NLO[1,:],2)+pow(dijet_pdf_unc_NLO[1,:],2))	#NLO up
		zjet_fnlo_errLO_down = np.sqrt(pow(zjet_scale_unc_LO[0,:],2)+pow(zjet_pdf_unc_LO[0,:],2))		#LO down
		zjet_fnlo_errLO_up = np.sqrt(pow(zjet_scale_unc_LO[1,:],2)+pow(zjet_pdf_unc_LO[1,:],2))			#LO up
		zjet_fnlo_errNLO_down = np.sqrt(pow(zjet_scale_unc_NLO[0,:],2)+pow(zjet_pdf_unc_NLO[0,:],2))	#NLO down
		zjet_fnlo_errNLO_up = np.sqrt(pow(zjet_scale_unc_NLO[1,:],2)+pow(zjet_pdf_unc_NLO[1,:],2))		#NLO up


		dijet_fnloLO_err = np.zeros([2, nentries])	#will contain absolute uncertainties
		dijet_fnloLO_err[0] = dijet_fnlo_errLO_down
		dijet_fnloLO_err[1] = dijet_fnlo_errLO_up
		dijet_fnloNLO_err = np.zeros([2, nentries])	#will contain absolute uncertainties
		dijet_fnloNLO_err[0] = dijet_fnlo_errNLO_down
		dijet_fnloNLO_err[1] = dijet_fnlo_errNLO_up
		
		zjet_fnloLO_err = np.zeros([2, nentries])	#will contain absolute uncertainties
		zjet_fnloLO_err[0] = zjet_fnlo_errLO_down
		zjet_fnloLO_err[1] = zjet_fnlo_errLO_up
		zjet_fnloNLO_err = np.zeros([2, nentries])	#will contain absolute uncertainties
		zjet_fnloNLO_err[0] = zjet_fnlo_errNLO_down
		zjet_fnloNLO_err[1] = zjet_fnlo_errNLO_up
		


		#-------------------------------------------#
		# 	Calculate the ratios and their errors	#
		#-------------------------------------------#
		# dijet/zjet in LO and in NLO	
		ratio_LO = np.divide(dijet_fnloxs_LO, zjet_fnloxs_LO)
		ratio_NLO = np.divide(dijet_fnloxs_NLO,zjet_fnloxs_NLO)

		#try simple error propagation for absolute error of ratio:
		dijet_relerr_LO_dw = np.divide(dijet_fnloLO_err[0], dijet_fnloxs_LO)
		dijet_relerr_LO_up = np.divide(dijet_fnloLO_err[1], dijet_fnloxs_LO)
		dijet_relerr_NLO_dw = np.divide(dijet_fnloNLO_err[0], dijet_fnloxs_NLO)
		dijet_relerr_NLO_up = np.divide(dijet_fnloNLO_err[1], dijet_fnloxs_NLO)

		zjet_relerr_LO_dw = np.divide(zjet_fnloLO_err[0], zjet_fnloxs_LO)
		zjet_relerr_LO_up = np.divide(zjet_fnloLO_err[1], zjet_fnloxs_LO)
		zjet_relerr_NLO_dw = np.divide(zjet_fnloNLO_err[0], zjet_fnloxs_NLO)
		zjet_relerr_NLO_up = np.divide(zjet_fnloNLO_err[1], zjet_fnloxs_NLO)

		#calculate ratio relative error:
		ratio_LO_relerr_dw = np.sqrt(pow(dijet_relerr_LO_dw,2)+pow(zjet_relerr_LO_dw,2))	#end up with abs values
		ratio_LO_relerr_up = np.sqrt(pow(dijet_relerr_LO_up,2)+pow(zjet_relerr_LO_up,2))	#end up with abs values
		ratio_NLO_relerr_dw = np.sqrt(pow(dijet_relerr_NLO_dw,2)+pow(zjet_relerr_NLO_dw,2))	#end up with abs values
		ratio_NLO_relerr_up = np.sqrt(pow(dijet_relerr_NLO_up,2)+pow(zjet_relerr_NLO_up,2))	#end up with abs values

		#calculate ratio absolute error:
		ratio_LO_err = np.zeros([2, nentries])
		ratio_LO_err[0] = np.multiply(ratio_LO, ratio_LO_relerr_dw)
		ratio_LO_err[1] = np.multiply(ratio_LO, ratio_LO_relerr_up)

		ratio_NLO_err = np.zeros([2, nentries])
		ratio_NLO_err[0] = np.multiply(ratio_NLO, ratio_NLO_relerr_dw)
		ratio_NLO_err[1] = np.multiply(ratio_NLO, ratio_NLO_relerr_up)


		#x_axis = (binbounds.T[0]+binbounds.T[1])/2.0
		x_axis = (dijet_binbounds_NLO.T[0]+dijet_binbounds_NLO.T[1])/2.0

		#plotting LO
		plot_ratio(bname, ax_ratio_bin_LO, ax_ratio_LO, x_axis, ratio_LO, ratio_LO_err, dijet_binbounds_NLO.T[0], dijet_binbounds_NLO.T[1], m_alpha, pdfset)
		#plotting NLO
		plot_ratio(bname, ax_ratio_bin_NLO, ax_ratio_NLO, x_axis, ratio_NLO, ratio_NLO_err, dijet_binbounds_NLO.T[0], dijet_binbounds_NLO.T[1], m_alpha, pdfset)

	#---------------------------#
	#	legend, title and save	#
	#---------------------------#
	#ax_ratio.grid(True, which="both", axis="y", color="k", linestyle="dotted", linewidth=0.8, alpha=0.2)
	x_minticks = [30, 40, 50, 60, 70, 80, 90, 200, 300, 400, 500, 600, 700, 800, 900, 2000, 3000, 4000]
	x_minticklabels = ['','','','60', '', '', '', '200', '', '', '500','', '700', '', '', '2000', '3000', '']
	x_majticks = [100, 1000]
	x_majticklabels = ['100', '1000']
	
	#-------------------#
	#	Leading Order	#
	#-------------------#
	ax_ratio_LO.grid(True, which="minor", axis="y", color="k", linestyle="dotted", linewidth=0.8, alpha=0.2)
	ax_ratio_LO.grid(True, which="major", axis="y", color="k", linestyle="dotted", linewidth=0.8, alpha=0.6)
	ax_ratio_LO.grid(True, which="both", axis="x", color="k", linestyle="dotted", linewidth=0.8, alpha=0.2)
	ax_ratio_LO.set_xlabel(r'$\mathrm{p_{T,avg}}$ /GeV', fontsize=12, horizontalalignment='right')
	ax_ratio_LO.xaxis.set_label_coords(1.00, -0.06)
	ax_ratio_LO.set_ylabel(r'XS ratios: Dijet / Zjet', fontsize=12, horizontalalignment='center', verticalalignment='top', y=0.5, rotation=90, labelpad=24)

	ax_ratio_LO.set_xticks(x_majticks, minor=False)	#set major ticks loc
	ax_ratio_LO.set_xticks(x_minticks, minor=True)		#set minor ticks loc
	ax_ratio_LO.set_xticklabels(x_majticklabels, minor=False)	#set major ticks label
	ax_ratio_LO.set_xticklabels(x_minticklabels, minor=True)	#set minor titcks label
	
	ax_ratio_LO.legend(fontsize=12, numpoints=1, loc='upper right')
	ax_ratio_LO.set_title("Cross Section Ratios: Dijet over Z+jet \n Fixed Order Theory Calculations in LO (%s)"%pdfset, fontsize=12)
	plt.tight_layout()
	fig_ratio_LO.savefig("ratio_summary_fnlo_LO.png")

	#now zoom in:
	ax_ratio_LO.set_ylim(1e3, 1e7)	#some points won't be displayed
	ax_ratio_LO.text(0.02, 0.04, "zoomed in y-axis", fontsize=12, fontstyle='italic', ha='left', va='bottom', transform=ax_ratio_LO.transAxes)
	plt.tight_layout()
	fig_ratio_LO.savefig("ratio_summary_fnlo_LO_zoomed.png")
	


	#---------------------------#
	#	Next-to-Leading Order	#
	#---------------------------#
	ax_ratio_NLO.grid(True, which="minor", axis="y", color="k", linestyle="dotted", linewidth=0.8, alpha=0.2)
	ax_ratio_NLO.grid(True, which="major", axis="y", color="k", linestyle="dotted", linewidth=0.8, alpha=0.6)
	ax_ratio_NLO.grid(True, which="both", axis="x", color="k", linestyle="dotted", linewidth=0.8, alpha=0.2)
	ax_ratio_NLO.set_xlabel(r'$\mathrm{p_{T,avg}}$ /GeV', fontsize=12, horizontalalignment='right')
	ax_ratio_NLO.xaxis.set_label_coords(1.00, -0.06)
	ax_ratio_NLO.set_ylabel(r'XS ratios: Dijet / Zjet', fontsize=12, horizontalalignment='center', verticalalignment='top', y=0.5, rotation=90, labelpad=24)

	ax_ratio_NLO.set_xticks(x_majticks, minor=False)	#set major ticks loc
	ax_ratio_NLO.set_xticks(x_minticks, minor=True)		#set minor ticks loc
	ax_ratio_NLO.set_xticklabels(x_majticklabels, minor=False)	#set major ticks label
	ax_ratio_NLO.set_xticklabels(x_minticklabels, minor=True)	#set minor titcks label
	
	ax_ratio_NLO.legend(fontsize=12, numpoints=1, loc='upper right')
	ax_ratio_NLO.set_title("Cross Section Ratios: Dijet over Z+jet \n Fixed Order Theory Calculations in NLO (%s)"%pdfset, fontsize=12)
	plt.tight_layout()
	fig_ratio_NLO.savefig("ratio_summary_fnlo_NLO.png")

	#now zoom in:
	ax_ratio_NLO.set_ylim(1e3, 1e7)	#some points won't be displayed
	ax_ratio_NLO.text(0.02, 0.04, "zoomed in y-axis", fontsize=12, fontstyle='italic', ha='left', va='bottom', transform=ax_ratio_NLO.transAxes)
	plt.tight_layout()
	fig_ratio_NLO.savefig("ratio_summary_fnlo_NLO_zoomed.png")
	

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








