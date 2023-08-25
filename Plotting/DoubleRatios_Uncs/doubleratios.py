#!/usr/bin/env python2
#-*- coding:utf-8 -*-


#####################################################
#													#
#	Plotting the double ratio 						#
#	(DijetMC/ZjetMC)/(Dijetdata/Zjetdata)			#
#	Taking input from the ratio ROOT file			#
#													#
#	New version, where MC/data and NLO/data,		#
#	i.e. theory/data								#
#	instead of data/theory.							#
#	See "plots_DoubleRatio_Uncs.py" for old			#
#	version.										#
#													#
#	Created by B. Schillinger, 22.01.2020			#
#	Last modified: B.Schillinger, 15.03.2020		#
#	--> changes in pdfunc calculation!				#
#													#
#####################################################

# old script with wrong uncertainty handling (and commented-out the calculation of a total unc with pdf and scale) in "doubleratios_wrongpdfunc.py")

import argparse
import glob, os, sys
import string
import timeit
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

#position of the big and small legend in the ybys bins when using pythia8
_legpos_pythia = {'yb0ys0':['lower left','upper left'], 'yb0ys1':['lower left','upper right'], 'yb0ys2':['lower right','upper right'], 'yb1ys0':['lower left','upper left'], 'yb1ys1':['lower left','upper right'], 'yb2ys0':['lower left','upper left']}

#legend positions when using madgraph+P8
_legpos_madgraph = {'yb0ys0':['upper left','lower left'], 'yb0ys1':['lower left','upper right'], 'yb0ys2':['upper right','lower right'], 'yb1ys0':['lower left','upper left'], 'yb1ys1':['upper left','lower left'], 'yb2ys0':['upper left','lower left']}

# legend positions when using herwig7
_legpos_herwig = {'yb0ys0':['upper left','lower left'], 'yb0ys1':['upper left','lower left'], 'yb0ys2':['upper right','lower right'], 'yb1ys0':['upper left','lower left'], 'yb1ys1':['upper right','lower left'], 'yb2ys0':['upper left','lower left']}

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
	


#need to calculate the ratios for all the members 0 to 56 of CT14nlo
#----> no need to give out the pdf unc the other way.
#need to calculate the ratio dijet/zjet for EACH of the pdf members
#function for calculating the pdfunc correctly -- partly copied from read_fnlotable() (see below)
def calc_pdfuncs(fnlotable_dijet, fnlotable_zjet, pdfset, order, scale_var_type):
	#calculate the central cross sections and central ratio.
	fnlocentral_dijet = fastnlo.fastNLOLHAPDF(fnlotable_dijet, pdfset, 0)
	fnlocentral_zjet = fastnlo.fastNLOLHAPDF(fnlotable_zjet, pdfset, 0)

	fnlocentral_dijet.CalcCrossSection()
	dijetxs_central = 1000*np.array(fnlocentral_dijet.GetCrossSection())	#1000 to get correct unit (convert from pb to fb)

	fnlocentral_zjet.CalcCrossSection()
	zjetxs_central = 1000*np.array(fnlocentral_zjet.GetCrossSection())		#1000 to get correct unit (convert from pb to fb)

	ratio_central = np.divide(dijetxs_central, zjetxs_central)
	ratio_mem_list = []		#list that will contain the XS ratio for each of the pdf members from 0 to 56
	memdevs_list = []		#list that will contain the difference (deviation) of the dijet/zjet ratio calculated with member N to the one with member 0

	#go through all the 57 members (mem=0 --> central, mem=1-56 eigenvector sets 90%)
	for m in range(0,57):	#0 is just a cross check here
		print "<calc_pdfuncs()>: current member (total 0 to 56): ", m
		#create fnlo object
		fnlomem_dijet = fastnlo.fastNLOLHAPDF(fnlotable_dijet, pdfset, m)		#takes pdfmember m 
		fnlomem_zjet = fastnlo.fastNLOLHAPDF(fnlotable_zjet, pdfset, m)
		if(order==0):	#only LO desired
			fnlomem_dijet.SetContributionON(fastnlo.kFixedOrder, 0, True)
			fnlomem_dijet.SetContributionON(fastnlo.kFixedOrder, 1, False)
			fnlomem_zjet.SetContributionON(fastnlo.kFixedOrder, 0, True)
			fnlomem_zjet.SetContributionON(fastnlo.kFixedOrder, 1, False)
		elif(order==1):	#switch LO and NLO on
			fnlomem_dijet.SetContributionON(fastnlo.kFixedOrder, 0, True)
			fnlomem_dijet.SetContributionON(fastnlo.kFixedOrder, 1, True)
			fnlomem_zjet.SetContributionON(fastnlo.kFixedOrder, 0, True)
			fnlomem_zjet.SetContributionON(fastnlo.kFixedOrder, 1, True)
			
		#Calculate current XS (i.e. for this pdf member) for dijet and zjet
		fnlomem_dijet.CalcCrossSection()
		dijetxs_mem = 1000*np.array(fnlomem_dijet.GetCrossSection())	#convert from pb to fb (multiply by 1000)
		fnlomem_zjet.CalcCrossSection()
		zjetxs_mem = 1000*np.array(fnlomem_zjet.GetCrossSection())		#convert from pb to fb (multiply by 1000)

		#calculate current XS ratio dijet/zjet
		ratio_mem_cur = np.divide(dijetxs_mem, zjetxs_mem)
		ratio_mem_list.append(ratio_mem_cur)

		#calculate the relative deviation of this member-ratio to the central ratio
		#memratio_cur = np.divide(ratio_mem_cur, ratio_central)

		#append it to the list containing each member-ratio to the central ratio
		#memberratios_list.append(memratio_cur)

		#Calculate the deviation from the central member ratio!! (absolute difference)
		mem_dev = np.subtract(ratio_mem_cur, ratio_central)	#should be binwise (in ptavg) subtraction by numpy
		#append it to the list containing all the deviations
		memdevs_list.append(mem_dev)


	#make list an array:
	# ratio_mem_arr = [[Ratios pt arr mem0],[ratios pt arr mem1], ..., [ratios pt arr mem56]]
	ratio_mem_arr = np.array(ratio_mem_list)

	#make also this list an array
	#memberratios = np.array(memberratios_list)
	memdevs = np.array(memdevs_list)

	#separate the ratios in "up" and "down" uncertainty contributions (could do this in principle already in first loop)
	#---> for this: loop through members AND then through ptavg bins!
	#---> decide for each ptavg bin individually if it is a down- or up-variation caused by the new member

	#create the quadratic sums per ptavg bin!!
	nbins_ptavg = np.size(memdevs,1)
	updev_qsum = np.zeros([1,nbins_ptavg])	#one-dimensional array containing the quadratic sums of the deviations in each ptavg bin
	dwdev_qsum = np.zeros([1,nbins_ptavg])

	print "memdevs: \n", memdevs
	print "memdevs.size: ", memdevs.size
	#print "np.size(memdevs

	#loop through pdf-members
	for m in range(0,57):
		#loop through ptavg bins in current member m (would be 64 for the ones from the root-files, but smaller for these fnlo tables)
		#get size of axis 1 (ptavg bins per member) -- axis 0 corresponds to the amount of members
		for b in range(0,nbins_ptavg):
			if(memdevs[m,b]>=0.0):					#upwards variation --> put it to the "up-uncs" list
				#upuncs_list.append(memdevs[m])		#this is a list of arrays containing the absolute differences of ratios m - ratios 0
				updev_qsum[0,b] += pow(memdevs[m,b],2)	#add in quadrature (absolute ratio-deviation in current bin squared)
			elif(memdevs[m,b]<0.0):					#downwards variation --> put it to the "down-uncs" list
				#dwuncs_list.append(memdevs[m])	
				dwdev_qsum[0,b] += pow(memdevs[m,b],2)

	#upuncs = np.array(upuncs_list)
	#dwuncs = np.array(dwuncs_list)

	#could do the following already within the for-loop above and save some variables...
	#add up the differences quadratically, seperately for up and down --> #absolute total deviation
	#for m in range(0,57)
	#	updev_qsum += pow(upuncs[m,:],2)
	#	dwdev_qsum += pow(dwuncs[m,:],2)

	#take the sqrt -- should be done element-wise using numpy
	updev_sqrt = np.sqrt(updev_qsum)
	dwdev_sqrt = np.sqrt(dwdev_qsum)

	#get the relative pdfunc
	pdfunc_up = np.divide(updev_sqrt, ratio_central)
	pdfunc_dw = np.divide(dwdev_sqrt, ratio_central)

	#Give out the uncertainties!
	print "Tables: ", fnlotable_dijet, fnlotable_zjet
	print "---------------------------------------------------------------------"
	print "Absolute uncertainties: "
	print "updev_sqrt: \n", updev_sqrt
	print "dwdev_sqrt: \n", dwdev_sqrt
	print "Dijet/Zjet Ratio relative pdf uncertainties (dijet and zjet correlated):"
	print "pdfunc_up: \n", pdfunc_up
	print "pdfunc_dw: \n", pdfunc_dw

	return ratio_central, updev_sqrt, dwdev_sqrt, pdfunc_up, pdfunc_dw	#output 0, 1, 2, 3, 4



#--- NEW NEW NEW ----#
#This one works correctly AND gives out NLO and LO at once
def calc_pdfuncs_LO_NLO(fnlotable_dijet, fnlotable_zjet, pdfset, scale_var_type):
	#calculate the central cross sections and central ratio.
	ratio_mem_list = []		#list that will contain the XS ratio for each of the pdf members from 0 to 56
	memdevs_list = []		#list that will contain the difference (deviation) of the dijet/zjet ratio calculated with member N to the one with member 0

	#Create the output variables:
	abs_pdfunc_up_LO = -99999999.9
	abs_pdfunc_dw_LO = -99999999.9
	rel_pdfunc_up_LO = -99999999.9
	rel_pdfunc_dw_LO = -99999999.9
	abs_pdfunc_up_NLO = -99999999.9
	abs_pdfunc_dw_NLO = -99999999.9
	rel_pdfunc_up_NLO = -99999999.9
	rel_pdfunc_dw_NLO = -99999999.9

	orderlist = [0,1]

	#first do LO, then do NLO
	for order in orderlist:
		#go through all the 57 members (mem=0 --> central, mem=1-56 eigenvector sets 90%)
		for m in range(0,57):	#0 is just a cross check here
			print "<calc_pdfuncs()>: current member (total 0 to 56): ", m
			#create fnlo object
			fnlomem_dijet = fastnlo.fastNLOLHAPDF(fnlotable_dijet, pdfset, m)		#takes pdfmember m 
			fnlomem_zjet = fastnlo.fastNLOLHAPDF(fnlotable_zjet, pdfset, m)
			if(order==0):	#only LO desired
				fnlomem_dijet.SetContributionON(fastnlo.kFixedOrder, 0, True)
				fnlomem_dijet.SetContributionON(fastnlo.kFixedOrder, 1, False)
				fnlomem_zjet.SetContributionON(fastnlo.kFixedOrder, 0, True)
				fnlomem_zjet.SetContributionON(fastnlo.kFixedOrder, 1, False)
			elif(order==1):	#switch LO and NLO on
				fnlomem_dijet.SetContributionON(fastnlo.kFixedOrder, 0, True)
				fnlomem_dijet.SetContributionON(fastnlo.kFixedOrder, 1, True)
				fnlomem_zjet.SetContributionON(fastnlo.kFixedOrder, 0, True)
				fnlomem_zjet.SetContributionON(fastnlo.kFixedOrder, 1, True)
				
			#Calculate current XS (i.e. for this pdf member) for dijet and zjet
			fnlomem_dijet.CalcCrossSection()
			dijetxs_mem = 1000*np.array(fnlomem_dijet.GetCrossSection())	#convert from pb to fb (multiply by 1000)
			fnlomem_zjet.CalcCrossSection()
			zjetxs_mem = 1000*np.array(fnlomem_zjet.GetCrossSection())		#convert from pb to fb (multiply by 1000)

			#calculate current XS ratio dijet/zjet
			ratio_mem_cur = np.divide(dijetxs_mem, zjetxs_mem)
			ratio_mem_list.append(ratio_mem_cur)

			#Calculate the deviation from the central member ratio!! (absolute difference)
			#mem_dev = np.subtract(ratio_mem_cur, ratio_central)	#should be binwise (in ptavg) subtraction by numpy --- CAUTION: this always takes NLO central... ERROR
			mem_dev = np.subtract(ratio_mem_cur, ratio_mem_list[0])	#should be binwise (in ptavg) subtraction by numpy

			#append it to the list containing all the deviations
			memdevs_list.append(mem_dev)


		#make list an array:
		# ratio_mem_arr = [[Ratios pt arr mem0],[ratios pt arr mem1], ..., [ratios pt arr mem56]]
		ratio_mem_arr = np.array(ratio_mem_list)

		#make also this list an array
		#memberratios = np.array(memberratios_list)
		memdevs = np.array(memdevs_list)

		#separate the ratios in "up" and "down" uncertainty contributions (could do this in principle already in first loop)
		#---> for this: loop through members AND then through ptavg bins!
		#---> decide for each ptavg bin individually if it is a down- or up-variation caused by the new member

		#create the quadratic sums per ptavg bin!!
		nbins_ptavg = np.size(memdevs,1)
		updev_qsum = np.zeros([1,nbins_ptavg])	#one-dimensional array containing the quadratic sums of the deviations in each ptavg bin
		dwdev_qsum = np.zeros([1,nbins_ptavg])

		print "memdevs: \n", memdevs
		print "memdevs.size: ", memdevs.size
		#print "np.size(memdevs

		#loop through pdf-members
		for m in range(0,57):
			#loop through ptavg bins in current member m (would be 64 for the ones from the root-files, but smaller for these fnlo tables)
			#get size of axis 1 (ptavg bins per member) -- axis 0 corresponds to the amount of members
			for b in range(0,nbins_ptavg):
				if(memdevs[m,b]>=0.0):					#upwards variation --> put it to the "up-uncs" list
					#upuncs_list.append(memdevs[m])		#this is a list of arrays containing the absolute differences of ratios m - ratios 0
					updev_qsum[0,b] += pow(memdevs[m,b],2)	#add in quadrature (absolute ratio-deviation in current bin squared)
				elif(memdevs[m,b]<0.0):					#downwards variation --> put it to the "down-uncs" list
					#dwuncs_list.append(memdevs[m])	
					dwdev_qsum[0,b] += pow(memdevs[m,b],2)

		#upuncs = np.array(upuncs_list)
		#dwuncs = np.array(dwuncs_list)

		#could do the following already within the for-loop above and save some variables...
		#add up the differences quadratically, seperately for up and down --> #absolute total deviation
		#for m in range(0,57)
		#	updev_qsum += pow(upuncs[m,:],2)
		#	dwdev_qsum += pow(dwuncs[m,:],2)

		#take the sqrt -- should be done element-wise using numpy
		updev_sqrt = np.sqrt(updev_qsum)
		dwdev_sqrt = np.sqrt(dwdev_qsum)

		#get the relative pdfunc --> divide by central ratio (member 0)
		pdfunc_up = np.divide(updev_sqrt, ratio_mem_arr[0])
		pdfunc_dw = np.divide(dwdev_sqrt, ratio_mem_arr[0])

		#Give out the uncertainties!
		print "Tables: ", fnlotable_dijet, fnlotable_zjet
		print "Current Order (0=LO, 1=NLO): ", order
		print "---------------------------------------------------------------------"
		print "Absolute uncertainties: "
		print "updev_sqrt: \n", updev_sqrt
		print "dwdev_sqrt: \n", dwdev_sqrt
		print "Dijet/Zjet Ratio relative pdf uncertainties (dijet and zjet correlated):"
		print "pdfunc_up: \n", pdfunc_up
		print "pdfunc_dw: \n", pdfunc_dw
		
		if(order==0):
			abs_pdfunc_up_LO = updev_sqrt
			abs_pdfunc_dw_LO = dwdev_sqrt
			rel_pdfunc_up_LO = pdfunc_up
			rel_pdfunc_dw_LO = pdfunc_dw
		elif(order==1):
			abs_pdfunc_up_NLO = updev_sqrt
			abs_pdfunc_dw_NLO = dwdev_sqrt
			rel_pdfunc_up_NLO = pdfunc_up
			rel_pdfunc_dw_NLO = pdfunc_dw

	return  abs_pdfunc_dw_LO, abs_pdfunc_up_LO, abs_pdfunc_dw_NLO, abs_pdfunc_up_NLO, rel_pdfunc_dw_LO, rel_pdfunc_up_LO, rel_pdfunc_dw_NLO, rel_pdfunc_up_NLO	#output 0, 1, 2, 3, 4, 5, 6, 7



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

	#NEW: giving out the relative scale- and pdf uncertainties
	# switch the indices, so that 0 is the downwards and 1 the upwards unc
	rel_scaleunc = np.empty([2, len(xs_array)])
	rel_scaleunc[0] = rel_scale_unc[2]
	rel_scaleunc[1] = rel_scale_unc[1]
	
	rel_pdfunc = np.empty([2, len(xs_array)])
	rel_pdfunc[0] = rel_pdf_unc[2]
	rel_pdfunc[1] = rel_pdf_unc[1]
	

	#get the x-axis
	binbounds = np.array(fnlo.GetObsBinsBounds(0))
	x_axis = (binbounds.T[0]+binbounds.T[1])/2.0



	#return xs_array, entries_err, x_axis	#x_axis should be the same as for the root hists
	return xs_array, abs_scale_unc, abs_pdf_unc, rel_scaleunc, rel_pdfunc, binbounds



#plotting of the double ratio
#adjusted to how "plot_ratio()" is used in plot_xs_ratios.py
def plot_doubleratio(bname, ax_doubleratio, x_axis, doubleratio, doubleratio_err, low_bb, up_bb, m_color, m_alpha, labeling):
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
	xerrors = np.array([xerr_low, xerr_up])


	#plot the double ratio theory/data (first zjet/dijet) errorbar
	ax_doubleratio.errorbar(x_axis, doubleratio, xerr=xerrors, yerr=doubleratio_err, elinewidth=1.2, linewidth=1.0, marker=_ybysbin_marker[bname], ms=6, color=m_color, fillstyle='none', fmt='.', label=labeling)
	



# ---- old version of this function ---- #
#function to plot data and mc for either dijet or zjet into same plot
#for individual ybys bins
#def plot_xs(bname, ax, ax_xs_sum, x_axis, data, data_err, data_low_bb, data_up_bb, mc, mc_err, m_alpha, runperiod, generator):
def plot_doubleratio_old(bname, ax, ax_doubleratio_sum, x_axis, doubleratio, doubleratio_err, low_bb, up_bb, m_alpha, runperiod, generator):
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
	ax.set_ylabel(r'Double Ratio XS theory/data', fontsize=12, horizontalalignment='center', verticalalignment='top', y=0.5, rotation=90, labelpad=24)

	if((runperiod!=None)and(generator!=None)):
		ax.set_title(r'Double Ratio theory/data (dijet/zjet) in %s (Run %s, %s)'%(bname,runperiod,generator))	#(first dijet/zjet)
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
	#parser.add_argument('doubleratios_rootfile', type=TFile, nargs='?',
	#					help='Required argument! File with double ratio Data/MC (Dijet/Zjet for both) that shall be plotted. (filename glob)')
	parser.add_argument('data_rootfile', type=TFile, nargs='?',
						help='Required argument! File with data results dijet/zjet that shall be plotted. (filename glob)')
	parser.add_argument('mc_rootfile', type=TFile, nargs='?',
						help='Required argument! File with MC results dijet/zjet that shall be plotted. (filename glob)')
	parser.add_argument('fnlo_tables_dijet', type=str, nargs=6,	#nargs='+'
						help='Required argument! FastNLO dijet tables for the six ybys bins to be evaluated.')
	parser.add_argument('fnlo_tables_zjet', type=str, nargs=6,	#nargs='+'
						help='Required argument! FastNLO zjet tables for the six ybys bins to be evaluated.')
	
	parser.add_argument('uncsratios_rootfile', type=TFile, nargs='?',
						help='Required argument! File with uncertainty ratios dijet/zjet that shall be plotted. (filename glob)')



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

	parser.add_argument('--zjet-unc', dest='zjet_unc', type=TFile, nargs='?', default=None,
						help='TFile containing the uncertainty sources (symmetrised as up/down) for zjet.')
	parser.add_argument('--dijet-unc', dest='dijet_unc', type=TFile, nargs='?', default=None,
						help='TFile containing the uncertainty sources (symmetrised as up/down) for dijet.')
	parser.add_argument('--extrauncs', '-x', default=False, action='store_true',								#type=bool,
						help='If chosen, draw grey band for assumed uncorrelated uncertainties for Z+jet and Dijet.')

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
	#ratio_rootfile = args['doubleratios_rootfile']
	mc_rootfile = args['mc_rootfile']
	data_rootfile = args['data_rootfile']
	fnlo_tables_dijet = args['fnlo_tables_dijet']	#this is a list
	fnlo_tables_zjet = args['fnlo_tables_zjet']	#this is a list
	uncs_rootfile = args['uncsratios_rootfile']

	#optional -- for individual uncertainties
	zjetfile = args['zjet_unc']
	dijetfile = args['dijet_unc']

	ybys = args['ybysbin']				#name of ybys bin
	m_alpha = args['alphamarker']
	pdfset = args['pdfset']
	runperiod = args['runperiod']	#optional argument, but good for labeling --> to be improved
	generator = args['generator']	#optional argument, pythia or herwig
	verbose = args['verbose']

	#scale variation choice for uncertainty:
	if args['asymmetric']:
		scale_var_type = fastnlo.kAsymmetricSixPoint
		variation_type = 'scale uncertainty (6P)'
	else:
		scale_var_type = fastnlo.kSymmetricTwoPoint
		variation_type = 'scale uncertainty (2P)'



	#rootfile = TFile(args['rootfile'])
	print('\n')
	#print("[doubleratios.py]: Taking DoubleRatio input from %s" %ratio_rootfile)
	print("[doubleratios.py]: Taking Uncertainty-Ratios input from %s" %uncs_rootfile)


	#prepare the lists (will later be arrays)
	if(ybys==None):
		bin_names = np.array(["yb0ys0", "yb0ys1", "yb0ys2", "yb1ys0", "yb1ys1", "yb2ys0"])
	else:
		bin_names = np.array([ybys])



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
	#for j in range(0, 6):
	for i in range(bin_names.size):
		#current bin name:
		bname = bin_names[i]

		#for the individual (ybys bin) doubleratio plot
		fig_doubleratio_bin = plt.figure(figsize=(9,7), dpi=300)
		ax_doubleratio_bin = plt.subplot(xscale="log", yscale="linear")
		patches_bin = []

		#---------------------------------------#
		# reading the single ratios Dijet/Z+jet	#
		#	and creating the doubleratios		#
		#---------------------------------------#
		data, data_err, data_xax, data_low_bb, data_up_bb = read_rootfile(data_rootfile, "Standard/%s/ratio_%s"%(bname, bname), verbose)
		mc, mc_err, mc_xax, mc_low_bb, mc_up_bb = read_rootfile(mc_rootfile, "Standard/%s/ratio_%s"%(bname, bname), verbose)

		#read only the unc-ratio value (bin bounds etc. should be identical to the double ratio above)
		unc_TotalNoFlavor = read_rootfile(uncs_rootfile, "Standard/%s/ratio_TotalNoFlavor_%s"%(bname, bname), verbose)[0]
		#unc_zjetZJetFlav_dijetQCDFlav = read_rootfile(uncs_rootfile, "Standard/%s/zjetFlavZJet_dijetFlavQCD_%s"%(bname, bname))[0] #old version  with zjet/dijet
		unc_dijetQCDFlav_zjetZJetFlav = read_rootfile(uncs_rootfile, "Standard/%s/dijetFlavQCD_zjetFlavZJet_%s"%(bname, bname), verbose)[0]	#new version

		#print("double ratios: \n %s"%doubleratio)

		#Calculate the new doubleratio


		#-----------------------#
		#		errors			#	--> check this with plot_xs_ratios.py (if shorter version possible)
		#-----------------------#
		#create array with unc_up, and with unc_down (symmetrised errors)
		print("unc_TotalNoFlavor: ", unc_TotalNoFlavor)
		unc_size = unc_TotalNoFlavor.size	#size is an attribute of the array, not callable (i.e., no size() )
		ones_arr = np.full(unc_size, 1.0)

		#TOTAL NO FLAVOR
		delta_TotNoFlav = 0.5*np.absolute(np.add(ones_arr, (-1*unc_TotalNoFlavor)))	#delta_unc = 0.5*abs(unc-1)

		#preserve the sign: (do not use absolute value)
		delta_TotNoFlav_up = 0.5*(np.add(-1*ones_arr, unc_TotalNoFlavor))	#DeltaUnc upwards --> multiply it by -1 for downwards
		delta_TotNoFlav_down = -1.0*delta_TotNoFlav_up

		TotNoFlav_up = np.add(ones_arr, delta_TotNoFlav_up)		#1+DeltaUnc
		TotNoFlav_down = np.add(ones_arr, delta_TotNoFlav_down)	#1-DeltaUnc


		#FLAVOR
		delta_Flav = 0.5*np.absolute(np.add(ones_arr, (-1*unc_dijetQCDFlav_zjetZJetFlav)))	#not used
		delta_Flav_up = 0.5*(np.add(-1*ones_arr, unc_dijetQCDFlav_zjetZJetFlav))	#DeltaUnc upwards --> multiply it by -1 for downwards
		delta_Flav_down = -1.0*delta_Flav_up

		Flav_up = np.add(ones_arr, delta_Flav_up)				#1+DeltaUnc
		Flav_down = np.add(ones_arr, delta_Flav_down)			#1-DeltaUnc
	

		#calculate total uncertainty
		#delta_tot = np.add(delta_TotNoFlav, delta_Flav) #--> MUSS QUADRATISCH ADDIERT WERDEN!!!! (siehe unten bei den extra uncs! da ist es richtig)

		#tot_up = np.add(ones_arr, delta_tot)
		#tot_down = np.add(ones_arr, -1.0*delta_tot)
		delta_tot = np.sqrt(pow(delta_TotNoFlav,2)+pow(delta_Flav,2))	#quadratic sum
		totunc_up = np.add(ones_arr, delta_tot)
		totunc_down = np.add(ones_arr, -1.0*delta_tot)
		

		#-----------------------------------------------#
		#	get results on ratio from theory (fnlo)		#
		#-----------------------------------------------#
		if(verbose==True): print "[doubleratios.py]: Reading fnlo table: ", fnlo_tables_dijet[i]
		dijet_fnloxs_LO, dijet_scale_unc_LO, dijet_pdf_unc_LO, dijet_rel_scaleunc_LO, dijet_rel_pdfunc_LO, dijet_binbounds_LO = read_fnlotable(fnlo_tables_dijet[i], pdfset, 0, scale_var_type)
		dijet_fnloxs_NLO, dijet_scale_unc_NLO, dijet_pdf_unc_NLO, dijet_rel_scaleunc_NLO, dijet_rel_pdfunc_NLO, dijet_binbounds_NLO = read_fnlotable(fnlo_tables_dijet[i], pdfset, 1, scale_var_type)
		zjet_fnloxs_LO, zjet_scale_unc_LO, zjet_pdf_unc_LO, zjet_rel_scaleunc_LO, zjet_rel_pdfunc_LO, zjet_binbounds_LO = read_fnlotable(fnlo_tables_zjet[i], pdfset, 0, scale_var_type)
		zjet_fnloxs_NLO, zjet_scale_unc_NLO, zjet_pdf_unc_NLO, zjet_rel_scaleunc_NLO, zjet_rel_pdfunc_NLO, zjet_binbounds_NLO = read_fnlotable(fnlo_tables_zjet[i], pdfset, 1, scale_var_type)

		nentries = dijet_fnloxs_NLO.size
		oldNentries = dijet_fnloxs_NLO.size

		#-----------------------------------------------------------------------------------#
		#	Remove the low-pTavg bins, that are empty due to dijet trigger (ptavg < 56GeV)	#
		#	Additionally: restrictions from fnlo tables (see above in original version)		#
		#	--> numpy array slicing															#
		#-----------------------------------------------------------------------------------#
		# as everything will be compared to the ratio in data --> cut away everythig below 56 GeV
		# this contstraint comes from the dijet analysis, where the lowest trigger turn-on has this effect
		# additionally --> cut away the higher bins that are empty due to Z+jet analysis (in fnlo)
		data_xax = data_xax[11:_ybysbin_xind_up_zjetdata[bname]]
		data = data[11:_ybysbin_xind_up_zjetdata[bname]]
		data_err = data_err[11:_ybysbin_xind_up_zjetdata[bname]]
		data_low_bb = data_low_bb[11:_ybysbin_xind_up_zjetdata[bname]]
		data_up_bb = data_up_bb[11:_ybysbin_xind_up_zjetdata[bname]]
		mc_err = mc_err[11:_ybysbin_xind_up_zjetdata[bname]]
		mc = mc[11:_ybysbin_xind_up_zjetdata[bname]]
		ones_arr = ones_arr[11:_ybysbin_xind_up_zjetdata[bname]]

		#do the same for the uncertainties!
		TotNoFlav_up = TotNoFlav_up[11:_ybysbin_xind_up_zjetdata[bname]]
		TotNoFlav_down = TotNoFlav_down[11:_ybysbin_xind_up_zjetdata[bname]]
		Flav_up = Flav_up[11:_ybysbin_xind_up_zjetdata[bname]]
		Flav_down = Flav_down[11:_ybysbin_xind_up_zjetdata[bname]]
		totunc_up = totunc_up[11:_ybysbin_xind_up_zjetdata[bname]]
		totunc_down = totunc_down[11:_ybysbin_xind_up_zjetdata[bname]]

		#now the fnlo stuff with different binning --> cut away the four lowest bins (< 56 GeV)
		#example
		#fnloLO_err = fnloLO_err[:,4:]	#take everything on first axis = up and down, and entries from index 4 onwards (>56GeV)

		dijet_fnloxs_LO = dijet_fnloxs_LO[4:_ybysbin_xind_up_zjetfnlo[bname]]
		dijet_fnloxs_NLO = dijet_fnloxs_NLO[4:_ybysbin_xind_up_zjetfnlo[bname]]
		zjet_fnloxs_LO = zjet_fnloxs_LO[4:_ybysbin_xind_up_zjetfnlo[bname]]
		zjet_fnloxs_NLO = zjet_fnloxs_NLO[4:_ybysbin_xind_up_zjetfnlo[bname]]

		dijet_rel_scaleunc_LO = dijet_rel_scaleunc_LO[:,4:_ybysbin_xind_up_zjetfnlo[bname]]
		dijet_rel_scaleunc_NLO = dijet_rel_scaleunc_NLO[:,4:_ybysbin_xind_up_zjetfnlo[bname]]
		dijet_rel_pdfunc_LO = dijet_rel_pdfunc_LO[:,4:_ybysbin_xind_up_zjetfnlo[bname]]
		dijet_rel_pdfunc_NLO = dijet_rel_pdfunc_NLO[:,4:_ybysbin_xind_up_zjetfnlo[bname]]
		zjet_rel_scaleunc_LO = zjet_rel_scaleunc_LO[:,4:_ybysbin_xind_up_zjetfnlo[bname]]
		zjet_rel_scaleunc_NLO = zjet_rel_scaleunc_NLO[:,4:_ybysbin_xind_up_zjetfnlo[bname]]
		zjet_rel_pdfunc_LO = zjet_rel_pdfunc_LO[:,4:_ybysbin_xind_up_zjetfnlo[bname]]
		zjet_rel_pdfunc_NLO = zjet_rel_pdfunc_NLO[:,4:_ybysbin_xind_up_zjetfnlo[bname]]

		nentries = dijet_fnloxs_NLO.size

		#make binbounds array (for flatten() later before using steppify())
		#binbounds = np.array([doubleratio_low_bb, doubleratio_up_bb])
		binbounds = np.array([data_low_bb, data_up_bb])


		#-----------------------------------------------------------#
		# 	THEORY: (fastNLO) Calculate the ratios and their errors	#
		#-----------------------------------------------------------#
		#dijet/zjet errors (= errors on the fnlo-ratio)
		# RATIOS IN THEORY
		fnloLO = np.divide(dijet_fnloxs_LO, zjet_fnloxs_LO)
		fnloNLO = np.divide(dijet_fnloxs_NLO, zjet_fnloxs_NLO)

		#SCALE UNC
		scaleunc_abs_LO = np.empty([2, len(dijet_fnloxs_LO)])
		scaleunc_abs_NLO = np.empty([2, len(dijet_fnloxs_NLO)])
		scaleunc_abs_LO[0] = np.multiply(np.sqrt(pow(dijet_rel_scaleunc_LO[0,:],2)+pow(zjet_rel_scaleunc_LO[0,:],2)), fnloLO)	#LO scaleunc down
		scaleunc_abs_LO[1] = np.multiply(np.sqrt(pow(dijet_rel_scaleunc_LO[1,:],2)+pow(zjet_rel_scaleunc_LO[1,:],2)), fnloLO)	#LO scaleunc up 
		scaleunc_abs_NLO[0] = np.multiply(np.sqrt(pow(dijet_rel_scaleunc_NLO[0,:],2)+pow(zjet_rel_scaleunc_NLO[0,:],2)), fnloNLO)	#NLO scaleunc down
		scaleunc_abs_NLO[1] = np.multiply(np.sqrt(pow(dijet_rel_scaleunc_NLO[1,:],2)+pow(zjet_rel_scaleunc_NLO[1,:],2)), fnloNLO)	#NLO scaleunc up

		#add band with (fully) correlated pdf and scale unc


		#PDF UNC -- fully UNcorrelated case (but this is not true. therefore look at correlation, see below, using calc_pdfuncs())
		'''
		pdfunc_abs_LO = np.empty([2, len(dijet_fnloxs_LO)])
		pdfunc_abs_NLO = np.empty([2, len(dijet_fnloxs_NLO)])
		pdfunc_abs_LO[0] = np.multiply(np.sqrt(pow(dijet_rel_pdfunc_LO[0,:],2)+pow(zjet_rel_pdfunc_LO[0,:],2)), fnloLO)	#LO scaleunc down
		pdfunc_abs_LO[1] = np.multiply(np.sqrt(pow(dijet_rel_pdfunc_LO[1,:],2)+pow(zjet_rel_pdfunc_LO[1,:],2)), fnloLO)	#LO scaleunc up
		pdfunc_abs_NLO[0] = np.multiply(np.sqrt(pow(dijet_rel_pdfunc_NLO[0,:],2)+pow(zjet_rel_pdfunc_NLO[0,:],2)), fnloNLO)	#LO scaleunc down
		pdfunc_abs_NLO[1] = np.multiply(np.sqrt(pow(dijet_rel_pdfunc_NLO[1,:],2)+pow(zjet_rel_pdfunc_NLO[1,:],2)), fnloNLO)	#LO scaleunc up
		'''


		#Newest version: return of calc_pdfuncs_LO_NLO() :
		#return  abs_pdfunc_dw_LO, abs_pdfunc_up_LO, abs_pdfunc_dw_NLO, abs_pdfunc_up_NLO, rel_pdfunc_dw_LO, rel_pdfunc_up_LO, rel_pdfunc_dw_NLO, rel_pdfunc_up_NLO	#output 0, 1, 2, 3, 4, 5, 6, 7

		#NEW: read out correlated pdfunc!
		#return of calc_pdfuncs():  ratio_central, updev_sqrt, dwdev_sqrt, pdfunc_up, pdfunc_dw	#output 0, 1, 2, 3, 4
		pdfunc_abs_LO = np.empty([2, oldNentries])	#still have the old length, will be sliced after readout
		pdfunc_abs_NLO = np.empty([2, oldNentries])
		#pdfunc_abs_LO[1], pdfunc_abs_LO[0] = calc_pdfuncs(fnlo_tables_dijet[i], fnlo_tables_zjet[i], pdfset, 0, scale_var_type)[1:3]	##unc LO upwards, unc LO downwards
		#pdfunc_abs_LO[1] = calc_pdfuncs(fnlo_tables_dijet[i], fnlo_tables_zjet[i], pdfset, 0, scale_var_type)[1]	#unc LO upwards
		#pdfunc_abs_NLO[1] = calc_pdfuncs(fnlo_tables_dijet[i], fnlo_tables_zjet[i], pdfset, 1, scale_var_type)[1]	#unc NLO upwards
		#pdfunc_abs_NLO[1], pdfunc_abs_NLO[0] = calc_pdfuncs(fnlo_tables_dijet[i], fnlo_tables_zjet[i], pdfset, 1, scale_var_type)[1:3]	#unc NLO downwards
		pdfunc_abs_LO[0], pdfunc_abs_LO[1], pdfunc_abs_NLO[0], pdfunc_abs_NLO[1] = calc_pdfuncs_LO_NLO(fnlo_tables_dijet[i], fnlo_tables_zjet[i], pdfset, scale_var_type)[0:4] #NEW!!! only one call per ybys :)

		print "These are the calculated abs. pdfuncs: "
		print "pdfunc_abs_LO: \n", pdfunc_abs_LO
		print "pdfunc_abs_NLO: \n", pdfunc_abs_NLO


		#slice the new unc arrays correctly
		pdfunc_abs_LO = pdfunc_abs_LO[:,4:_ybysbin_xind_up_zjetfnlo[bname]]
		pdfunc_abs_NLO = pdfunc_abs_NLO[:,4:_ybysbin_xind_up_zjetfnlo[bname]]



		#testing:
		print "scaleunc_abs_LO:", scaleunc_abs_LO


		#----------- end of fnlo calculations of dijet/zjet ratio ---------------#

		#---------------------------------------------------------------------------#
		# 	Calculate the doubleratios theory/data and SCALE their errors	by data #
		#---------------------------------------------------------------------------#
		# so far: data/MC, data/LO, data/NLO --> should add for dijet: three MC sets, zjet: two MC sets
		##ratio_datamc = np.divide(data, mc)
		##ratio_datafnloLO = np.divide(data, fnloxs_LO)
		##ratio_datafnloNLO = np.divide(data, fnloxs_NLO)
		doubleratio_mcdata = np.divide(mc, data)
		print "data: --> ", data
		print "mc: -->", mc
		print("doubleratio_mcdata: \n ", doubleratio_mcdata)
		doubleratio_fnloLOdata = np.divide(fnloLO, data)		# adjusted to doubleratio naming
		doubleratio_fnloNLOdata = np.divide(fnloNLO, data)	# --- " ---
		doubleratio_datadata = np.divide(data, data) #--> must be an array of 1.0 (used to draw data stat. error around 1.0)

		#Scale the errors by data:
		#data error will be drawn around 1.0 line
		#data_err = data_err

		# MC and data STATISTIC UNC #
		doubleratio_mc_err = np.divide(mc_err, data)
		doubleratio_data_err = np.divide(data_err, data)
		
		# SCALE UNC #
		doubleratio_fnloLO_scaleunc = np.zeros([2, nentries])
		doubleratio_fnloLO_scaleunc[0] = np.divide(scaleunc_abs_LO[0],data)		#scale scaleunc by data
		doubleratio_fnloLO_scaleunc[1] = np.divide(scaleunc_abs_LO[1],data)

		doubleratio_fnloNLO_scaleunc = np.zeros([2, nentries])
		doubleratio_fnloNLO_scaleunc[0] = np.divide(scaleunc_abs_NLO[0],data)	#scale scaleunc by data
		doubleratio_fnloNLO_scaleunc[1] = np.divide(scaleunc_abs_NLO[1],data)

		#PDF UNC#
		doubleratio_fnloLO_pdfunc = np.zeros([2, nentries])
		doubleratio_fnloLO_pdfunc[0] = np.divide(pdfunc_abs_LO[0],data)			#scale pdfunc by data
		doubleratio_fnloLO_pdfunc[1] = np.divide(pdfunc_abs_LO[1],data)

		doubleratio_fnloNLO_pdfunc = np.zeros([2, nentries])
		doubleratio_fnloNLO_pdfunc[0] = np.divide(pdfunc_abs_NLO[0],data)		#scale pdfunc by data
		doubleratio_fnloNLO_pdfunc[1] = np.divide(pdfunc_abs_NLO[1],data)
	
		#create values for scaleunc-band around theory/data ratio
		scaleunc_band_LO = np.zeros([2,nentries])
		#scaleunc_band_LO[0] = np.subtract(ratio_fnloLOdata, ratio_fnloLO_scaleunc[0])	#lower line in scaleunc band
		scaleunc_band_LO[0] = np.add(doubleratio_fnloLOdata, -1.0*doubleratio_fnloLO_scaleunc[0])		#lower line in scaleunc band --> scaleunc itself is in this case always positive (when calculated from dijet and zjet)
		scaleunc_band_LO[1] = np.add(doubleratio_fnloLOdata, doubleratio_fnloLO_scaleunc[1])		#upper line in scaleunc band
		scaleunc_band_NLO = np.zeros([2,nentries])
		scaleunc_band_NLO[0] = np.add(doubleratio_fnloNLOdata, -1.0*doubleratio_fnloNLO_scaleunc[0])	#lower line in scaleunc band
		scaleunc_band_NLO[1] = np.add(doubleratio_fnloNLOdata, doubleratio_fnloNLO_scaleunc[1])		#upper line in scaleunc band

		print "SCALE UNC NLO", scaleunc_abs_NLO
		print "doubleratio ScaleuncNLO/data", doubleratio_fnloNLO_scaleunc
		print "ScaleuncBandNLO", scaleunc_band_NLO
		print "doubleratio NLO pdfunc/data", doubleratio_fnloNLO_pdfunc
		print "doubleratio NLO", doubleratio_fnloNLOdata

		#make lumi unc arrays:
		lumi_up = np.full(nentries, 1+0.025)
		lumi_dw = np.full(nentries, 1-0.025)
	

		#if extrauncs==true --> plot first the assumption of "uncorrelated dijet and zjet uncs" (grey band)
		if(args['extrauncs']==True):
			#read the files for Z+jet
			zjet_unc_TotalNoFlavor = read_rootfile(zjetfile, "UncertaintyRatios/%s/UpDownDiv_TotalNoFlavor_%s"%(bname, bname), verbose)[0]	
			zjet_unc_FlavorZJet = read_rootfile(zjetfile, "UncertaintyRatios/%s/UpDownDiv_FlavorZJet_%s"%(bname, bname), verbose)[0]	

			#slice immediately:
			zjet_unc_TotalNoFlavor = zjet_unc_TotalNoFlavor[11:_ybysbin_xind_up_zjetdata[bname]]
			zjet_unc_FlavorZJet = zjet_unc_FlavorZJet[11:_ybysbin_xind_up_zjetdata[bname]]

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
			dijet_unc_TotalNoFlavor = read_rootfile(dijetfile, "UncertaintyRatios/%s/UpDownDiv_TotalNoFlavor_%s"%(bname, bname), verbose)[0]
			dijet_unc_FlavorQCD = read_rootfile(dijetfile, "UncertaintyRatios/%s/UpDownDiv_FlavorQCD_%s"%(bname, bname), verbose)[0]

			#slice immediately:
			dijet_unc_TotalNoFlavor = dijet_unc_TotalNoFlavor[11:_ybysbin_xind_up_zjetdata[bname]]
			dijet_unc_FlavorQCD = dijet_unc_FlavorQCD[11:_ybysbin_xind_up_zjetdata[bname]]

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
			#in an earlier version of this script --> had bug: did not take quadratic sum. now it is fixed!
			tot_zjet_up = np.sqrt(np.add(delta_zjet_TotNoFlav_up*delta_zjet_TotNoFlav_up, delta_zjet_FlavZJet_up*delta_zjet_FlavZJet_up))				#total zjet uncertainty up ----> QUADRATIC SUM !!!
			tot_zjet_down = np.sqrt(np.add(delta_zjet_TotNoFlav_down*delta_zjet_TotNoFlav_down, delta_zjet_FlavZJet_down*delta_zjet_FlavZJet_down)) 	#total zjet uncertainty down ---> QUADRATIC SUM !!!

			tot_dijet_up = np.sqrt(np.add(delta_dijet_TotNoFlav_up*delta_dijet_TotNoFlav_up, delta_dijet_FlavQCD_up*delta_dijet_FlavQCD_up))			#total dijet uncertainty up
			tot_dijet_down = np.sqrt(np.add(delta_dijet_TotNoFlav_down*delta_dijet_TotNoFlav_down, delta_dijet_FlavQCD_down*delta_dijet_FlavQCD_down))	#total dijet uncertainty down


			delta_tot_uncorr_up = np.sqrt(np.add(tot_zjet_up*tot_zjet_up, tot_dijet_up*tot_dijet_up))			#total uncorr unc up = Sqrt[(zjetup^2)+(dijetup^2)]
			delta_tot_uncorr_down = np.sqrt(np.add(tot_zjet_down*tot_zjet_down, tot_dijet_down*tot_dijet_down))	#total uncorr unc down (should be the same as symmetrised...)

			#plotted around 1
			tot_uncorr_up = np.add(ones_arr, delta_tot_uncorr_up)				#1.0+tot_uncorr
			tot_uncorr_down = np.add(ones_arr, -1.0*delta_tot_uncorr_down)			#1.0-tot_uncorr


			#------------------------------------
			# plot the grey uncorr unc band first (below everything)
			#------------------------------------
			#plot_unc_filled(bname, ax_doubleratio_bin, patches_bin, doubleratio_xax, binbounds, tot_uncorr_down, tot_uncorr_up, 'powderblue', 'powderblue', "Quadratic Sum (zjet, dijet)", runperiod, generator)
			plot_unc_filled(bname, ax_doubleratio_bin, patches_bin, data_xax, binbounds, tot_uncorr_down, tot_uncorr_up, 'powderblue', 'powderblue', 0.6, "Quadratic Sum (zjet, dijet)", runperiod, generator)


		#-----------------------#
		#		plotting		#
		#-----------------------#
		#creating the plots for the current ybys bin
		#plotting total uncertainty:
		plot_unc_filled(bname, ax_doubleratio_bin, patches_bin, data_xax, binbounds, totunc_down, totunc_up, 'mediumseagreen', 'mediumseagreen', 0.6, "dijet/zjet total JES uncertainty (H7)", runperiod, generator)
		
		#draw scaleunc-bands
		plot_unc_filled(bname, ax_doubleratio_bin, patches_bin, data_xax, binbounds, scaleunc_band_LO[0], scaleunc_band_LO[1], 'yellow', 'yellow', 0.1, "LO %s"%variation_type,runperiod,generator)#scaleunc LO
		plot_unc_filled(bname, ax_doubleratio_bin, patches_bin, data_xax, binbounds, scaleunc_band_NLO[0], scaleunc_band_NLO[1], 'tan', 'tan', 0.6, "NLO %s"%variation_type,runperiod,generator)#scaleunc NLO



		#plotting TotalNoFlavor uncertainty band
		#plot_unc_filled(bname, ax_doubleratio_bin, patches_bin, doubleratio_xax, binbounds, TotNoFlav_down, TotNoFlav_up, 'royalblue', 'royalblue', "zjet/dijet jesTotalNoFlavor (H7)", runperiod, generator)
		plot_indiv_uncs(bname, ax_doubleratio_bin, patches_bin, data_xax, binbounds, TotNoFlav_up, 'royalblue', 'solid', "jesTotalNoFlavor up (H7)")#legbool=False (draw legend individually manually)
		plot_indiv_uncs(bname, ax_doubleratio_bin, patches_bin, data_xax, binbounds, TotNoFlav_down, 'royalblue', 'dashed', "jesTotalNoFlavor down (H7)")


		#plotting zjetZJetFlavor_dijetQCDFlavor uncertainty band
		#plot_unc_filled(bname, ax_doubleratio_bin, patches_bin, doubleratio_xax, binbounds, Flav_down, Flav_up, 'tomato', 'tomato', "zjetZJetFlav unc / dijetQCDFlav unc (H7)", runperiod, generator)
		##plot_indiv_uncs(bname, ax_doubleratio_bin, patches_bin, doubleratio_xax, binbounds, Flav_up, 'tomato', 'solid', "zjetZJetFlav unc / dijetQCDFlav unc up (H7)", runperiod, generator) 		#old version
		##plot_indiv_uncs(bname, ax_doubleratio_bin, patches_bin, doubleratio_xax, binbounds, Flav_down, 'tomato', 'dashed', "zjetZJetFlav unc / dijetQCDFlav unc down (H7)", runperiod, generator)	#old version
		plot_indiv_uncs(bname, ax_doubleratio_bin, patches_bin, data_xax, binbounds, Flav_up, 'tomato', 'solid', "dijetQCDFlav unc / zjetZJetFlav unc up (H7)") 		#new version
		plot_indiv_uncs(bname, ax_doubleratio_bin, patches_bin, data_xax, binbounds, Flav_down, 'tomato', 'dashed', "dijetQCDFlav unc / zjetZJetFlav unc down (H7)")	#new version

		#plot lumi:
		plot_indiv_uncs(bname, ax_doubleratio_bin, patches_bin, data_xax, binbounds, lumi_up, 'khaki', 'solid', "Luminosity (2.5 %)")		#lumi
		plot_indiv_uncs(bname, ax_doubleratio_bin, patches_bin, data_xax, binbounds, lumi_dw, 'khaki', 'dashed', "Luminosity (2.5 %)")	#lumi


		#plotting the doubleratio
		#-----------------------
		#plot_doubleratio(bname, ax_doubleratio_bin, ax_doubleratio_sum, data_xax, doubleratio, doubleratio_err, doubleratio_low_bb, doubleratio_up_bb, m_alpha, runperiod, generator)
		#plot the statistical error of the data around 1.0
		ax_doubleratio_bin.errorbar(data_xax, doubleratio_datadata, xerr=None, yerr=np.absolute(doubleratio_data_err), elinewidth=1.2, linewidth=1.0, marker='.', ms=2, color=_ybysbin_col[bname], fillstyle='none', fmt='.', label="data ratio (Run%s) stat. unc."%runperiod)#data stat. err around 1.0

		plot_doubleratio(bname, ax_doubleratio_bin, data_xax, doubleratio_mcdata, np.absolute(doubleratio_mc_err), binbounds[0,:], binbounds[1,:], _ybysbin_col[bname], m_alpha, "MC (%s) / data (Run%s) with mc stat. unc."%(generator, runperiod))#mc/data
		plot_doubleratio(bname, ax_doubleratio_bin, data_xax, doubleratio_fnloLOdata, np.absolute(doubleratio_fnloLO_pdfunc), binbounds[0,:], binbounds[1,:], 'gray', m_alpha, "LO theo. (%s) / data (Run%s) with PDF unc."%(pdfset, runperiod))#LO/data
		plot_doubleratio(bname, ax_doubleratio_bin, data_xax, doubleratio_fnloNLOdata, np.absolute(doubleratio_fnloNLO_pdfunc), binbounds[0,:], binbounds[1,:], 'black', m_alpha, "NLO theo. (%s) / data (Run%s) with PDF unc."%(pdfset, runperiod))#NLO/data

		#print some axis info for testing
		#print(ax_doubleratio_bin.xaxis.get_minor_ticks()[-2].label1)
		#ax_doubleratio_bin.xaxis.get_minor_ticks()[-2].label1.set_visible(False)	#test to remove a label


		#-------------------#
		# legend and save 	#
		#-------------------#
		#add line at 1.0
		ax_doubleratio_bin.axhline(y=1.0, xmin=0, xmax=1, linewidth=0.6, color='grey', linestyle='dashed', dashes=(5,10))	#linewidth=0.4
		#ax_ratio_bin.grid(True, which="major", axis="y", color="k", linestyle="dotted", linewidth=0.8, alpha=0.2)
		ax_doubleratio_bin.grid(True, which="both", axis="y", color="k", linestyle="dotted", linewidth=0.8, alpha=0.2)
		ax_doubleratio_bin.grid(True, which="both", axis="x", color="k", linestyle="dotted", linewidth=0.8, alpha=0.2)

		#finalising (used to be done in plot_ratio() fct, but do it only once per bin, here)
		ax_doubleratio_bin.set_xlabel(r'$\mathrm{p_{T,avg}}$ /GeV', fontsize=12, horizontalalignment='right')
		ax_doubleratio_bin.xaxis.set_label_coords(1.00, -0.06)
		ax_doubleratio_bin.set_ylabel(r'Double Ratios', fontsize=12, horizontalalignment='center', verticalalignment='top', y=0.5, rotation=90, labelpad=24)

		#ax.set_title(r'%s in %s (%s)'%(process,bname, _ybysbin_label[bname]))
		#ax_ratio_bin.set_title("%s XS ratios in %s (Data Run%s, MC %s, FixedOrder Theory with %s)"%(process, _ybysbin_label[bname], runperiod, generator, pdfset), fontsize=11)
		#ax_doubleratio_bin.set_title("%s XS ratios in %s (%s)"%(process, bname, _ybysbin_label[bname]), fontsize=12)
		ax_doubleratio_bin.set_title("XS double ratios in %s (%s)"%(bname, _ybysbin_label[bname]), fontsize=12)
		ax_doubleratio_bin.legend(fontsize=10, numpoints=1, loc='upper left')

		#legend position depending on ybys bin and mc generator in use
		if(generator=="pythia8"):
			legpos = _legpos_pythia[bname]	#list containing [big legend position, small legend position]
		elif(generator=="herwig7"):
			legpos = _legpos_herwig[bname]	#list
		elif(generator=="HT_madgraph_pythia8"):
			legpos = _legpos_madgraph[bname]

		#reordering of the legend entries
		#8 entries: uncorr.unc, totunc, scaleunc LO, scaleunc NLO, data stat.err, mc/data, LO/data, NLO/data
		ratiohandles, ratiolabels = ax_doubleratio_bin.get_legend_handles_labels()	#7 entries (totunc, scaleuncLO, scaleuncNLO, data, MC, LO, NLO)
		print "ratiohandles: ", ratiohandles
		print "ratiolabels: ", ratiolabels
		#order_entries = [3,4,5,6,0,1,2]
		order_entries = [4,5,6,7,0,1,2,3]
		#ax_doubleratio_bin.legend([ratiohandles[idx] for idx in order_entries], [ratiolabels[idx] for idx in order_entries], fontsize=10, numpoints=1, loc='upper left')
		ax_doubleratio_bin.legend([ratiohandles[idx] for idx in order_entries], [ratiolabels[idx] for idx in order_entries], fontsize=10, numpoints=1, loc=legpos[0])	#legpos depending on generator and bin

		ratio_legend = ax_doubleratio_bin.get_legend()	#get the "old" legend containing all the ratios
		ax_doubleratio_bin.add_artist(ratio_legend)		#draw the old legend
		#unc_legend = plt.legend([("royalblue","dashed"),("tomato","dashed"),("khaki","dashed")], ["jesTotalNoFlavor unc. up, down (H7)",'%s unc. up, down (H7)'%flavname, r"Luminosity uncertainty ($\pm$ 2.5%)"], handler_map={tuple: AnyObjectHandler()}, loc='lower left')	#make the unc legend
		#unc_legend = plt.legend([("royalblue","dashed"),("tomato","dashed"),("khaki","dashed")], ["jesTotalNoFlavor unc. up, down (H7)",'jesFlavor (dijetQCD/zjetZJet) unc. up, down (H7)', r"Luminosity uncertainty ($\pm$ 2.5%)"], handler_map={tuple: AnyObjectHandler()}, loc='lower left')	#make the unc legend
		unc_legend = plt.legend([("royalblue","dashed"),("tomato","dashed"),("khaki","dashed")], ["jesTotalNoFlavor unc. up, down (H7)",'jesFlavor (dijetQCD/zjetZJet) unc. up, down (H7)', r"Luminosity uncertainty ($\pm$ 2.5%)"], handler_map={tuple: AnyObjectHandler()}, loc=legpos[1])	#make the unc legend (depends on generator and bin)


		#plt.tick_params(axis='x', which='minor')				#?? not here but in main code
		#ax_ratio_bin.xaxis.set_minor_formatter(FormatStrFormatter("%.d"))	# same???
	

		xlocs = ax_doubleratio_bin.get_xticks(minor=True)
		xlabels = ax_doubleratio_bin.get_xticklabels(minor=True)
		print "xlocs: "
		print xlocs
		print "xlabels: "
		print xlabels

		#xlocs --> 48 minor label locations
		#xlabels --> 48 minor labels

		#use zjet binning and TO DO: cut lower part
		x_minticks = [40, 50, 60, 70, 80, 90, 200, 300, 400, 500, 600, 700, 800, 900]
		x_minticklabels = ['','','60', '', '', '', '200', '', '', '500','', '700', '', '']
		x_majticks = [100, 1000]
		x_majticklabels = ['100', '1000']
		ax_doubleratio_bin.set_xticks(x_majticks, minor=False)				#set major ticks loc
		ax_doubleratio_bin.set_xticks(x_minticks, minor=True)				#set minor ticks loc
		ax_doubleratio_bin.set_xticklabels(x_majticklabels, minor=False)	#set major ticks label
		ax_doubleratio_bin.set_xticklabels(x_minticklabels, minor=True)		#set minor titcks label

		#ax_ratio_bin.set_yticks(np.arange(0.0,2.0, step=0.1))
		y_minticks = [0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.10, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 1.90]
		y_minticklabels = ['', '0.20', '', '0.40', '', '0.60', '', '0.80', '', '', '1.20', '', '1.40', '', '1.60', '', '1.80', '']
		y_majticks = [0.00, 1.00, 2.00]
		y_majticklabels = ['0.00', '1.00', '2.00']
		ax_doubleratio_bin.set_yticks(y_majticks, minor=False)				#set major ticks loc
		ax_doubleratio_bin.set_yticks(y_minticks, minor=True)				#set minor ticks loc
		ax_doubleratio_bin.set_yticklabels(y_majticklabels, minor=False)	#set major ticks label
		ax_doubleratio_bin.set_yticklabels(y_minticklabels, minor=True)		#set minor titcks label

		for tick in ax_doubleratio_bin.yaxis.get_minor_ticks():
			tick.label.set_fontsize(8)

		#ax_ratio_bin.set_ylim(0.6, 1.4)	#test zoomed in
		ax_doubleratio_bin.set_ylim(0.0, 2.0)		#test


		plt.tight_layout()
		#fig_xs.savefig("%s_xs_overview_theory_Run%s_%s.png"%(process, runperiod, generator))
		#fig_ratio_bin.savefig("%s_xs_ratios_%s.png"%(process, bname))
		#fig_ratio_bin.savefig("%s_xs_ratios_%s_Run%s_%s_%s.png"%(process, bname, runperiod, generator, pdfset))

		#fig_doubleratio_bin.savefig("doubleratios_%s_Run%s_%s_%s.png"%(bname, runperiod, generator, pdfset))
		fig_doubleratio_bin.savefig("doubleratios_%s_Run%s_%s_%s_custom-legpos.png"%(bname, runperiod, generator, pdfset))

	
	#---------------------------------------#
	#	finalising the summary plots		#
	#---------------------------------------#
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
