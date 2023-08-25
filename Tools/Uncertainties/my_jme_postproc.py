#!/usr/bin/env python
#################################################
#						#
# Calculating JEC etc. using nanoAODtools	#
# See git repository for module authors		#
#						#
# This helper script:				# 
#	B.Schillinger, 08.01.2020		#
#						#
#################################################

#!!!
#NOTE: this script is under construction --> could be more analoguous to nano_postproc.py, but left out the arguments that were not used too often anyway, as running on skimmed data sets here
#!!!

import os, sys
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
from importlib import import_module
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor

#recommended module for jme uncertainties (see NanoAOD workbook on CMS twiki)
from PhysicsTools.NanoAODTools.postprocessing.modules.jme.jetmetHelperRun2 import *

if __name__ == "__main__":
	from optparse import OptionParser
	parser = OptionParser(usage="%prog [options] outputDir inputFiles")
	parser.add_option("-s", "--postfix",dest="postfix", type="string", default=None, help="Postfix which will be appended to the file name (default: _Friend for friends, _Skim for skims)")
	parser.add_option("-b", "--branch-selection",  dest="branchsel", type="string", default=None, help="Branch selection")
	parser.add_option("--bi", "--branch-selection-input",  dest="branchsel_in", type="string", default=None, help="Branch selection input")
	parser.add_option("--bo", "--branch-selection-output",  dest="branchsel_out", type="string", default=None, help="Branch selection output")
	parser.add_option("--friend",  dest="friend", action="store_true", default=False, help="Produce friend trees in output (current default is to produce full trees)")
	parser.add_option("--full",  dest="friend", action="store_false",  default=False, help="Produce full trees in output (this is the current default)")

	parser.add_option("-N", "--max-entries", dest="maxEntries", type="long",  default=None, help="Maximum number of entries to process from any single given input tree")
	parser.add_option("--first-entry", dest="firstEntry", type="long",  default=0, help="First entry to process in the three (to be used together with --max-entries)")

	#to decide which run period
	parser.add_option("-R", "--runperiod", dest="runperiod", type="string", default=None, help="Specify which runperiod is being looked at. -- for data.")
	parser.add_option("--mc", "--montecarlo", dest="montecarlo", action="store_true", default=False, help="Choose when looking at MC.")
	#parser.add_option("--data", dest="montecarlo", action="store_false", default=False, help="Choose if data ?")


	(options, args) = parser.parse_args()

	if len(args) < 2 :
	 parser.print_help()
	 sys.exit(1)
	outdir = args[0]; args = args[1:]


	#default settings:
	#createJMECorrector(isMC=True, dataYear=2016, runPeriod="B", jesUncert="Total", redojec=False, jetType = "AK4PFchs", noGroom=False, metBranchName="MET", applySmearing=True, isFastSim=False)


	#set desired jme corrections
	#---------------------------
	#check if these are correct!
	#do this on skimmed files (where JSON has already been applied)

	#dijet data:
	if((options.runperiod!=None)&&(options.montecarlo==False)):
		if(runperiod=="A"):
			jmeCorrections = createJMECorrector(False, "2018", "A", "Total", True, "AK4PFchs", False)	#RunA
		elif(runperiod=="B"):
			jmeCorrections = createJMECorrector(False, "2018", "B", "Total", True, "AK4PFchs", False)	#RunB
		elif(runperiod=="C"):
			jmeCorrections = createJMECorrector(False, "2018", "C", "Total", True, "AK4PFchs", False)	#RunC
		elif(runperiod=="D"):
			jmeCorrections = createJMECorrector(False, "2018", "D", "Total", True, "AK4PFchs", False)	#RunD
	#mc:
	elif(montecarlo==True):
		jmeCorrections = createJMECorrector(isMC=True, dataYear=2018, jesUncert="Total", redojec=True)	#or "All" for jesUncert?
	else:
		sys.exit("Something went wrong. Check runperiod or data/mc switch")



	#set the modules:
	modules = [jmeCorrections()]

	#check if any branch selection desired
	if options.branchsel!=None:
		options.branchsel_in = options.branchsel
		options.branchsel_out = options.branchsel

 
	#default setting postprocessor:
	#copied from source script:
	'''
	class PostProcessor :
	    def __init__(self,outputDir,inputFiles,cut=None,branchsel=None,modules=[],compression="LZMA:9",friend=False,postfix=None,
			 jsonInput=None,noOut=False,justcount=False,provenance=False,haddFileName=None,fwkJobReport=False,histFileName=None,histDirName=None, outputbranchsel=None,maxEntries=None,firstEntry=0,
	prefetch=False,longTermCache=False):
	'''


	#call the postprocessor			(What means provenance here? in jme example set to true..)
	#----------------------
	#calling only some keyword arguments --> the rest will go to default
	p=PostProcessor(outdir, args,
		branchsel = options.branchsel_in,
		modules = modules,
		friend = options.friend,
		postfix = options.postfix,
		provenance = True,
		maxEntries = options.maxEntries,
		firstEntry = options.firstEntry,
		outputbranchsel = options.branchsel_out)
	p.run()


