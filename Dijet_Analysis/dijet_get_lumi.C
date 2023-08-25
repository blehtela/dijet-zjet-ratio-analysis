#define DijetAnalysis_cxx
#define get_lumi_txt_cxx	//??
#include "DijetAnalysis.h"
#include <string>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "tdrstyle_mod14.C"
#include <TMath.h>
#include "Math/LorentzVector.h"
#include <iostream>
#include <fstream>

using namespace std;
using namespace ROOT::Math;


// script for creating a JSON file via reading "run" and "LuminosityBlock" from nanoAOD file.
// for later use with brilcalc



//void get_lumi_txt::Loop()
//void DijetAnalysis::Loop()
void DijetAnalysis::get_lumi_txt()
{

	setTDRStyle();
	fChain->SetBranchStatus("*",0);	//switch off all the branches
	fChain->SetBranchStatus("run",1);
	fChain->SetBranchStatus("luminosityBlock",1);

	TFile* json_output = new TFile("run_lumisec.txt", "RECREATE");
	json_output->cd();


	ofstream txt_file;
	txt_file.open("txt_file_test.txt");


	//create map for run number and lumisections
	map<int, set<pair<int, int>>> lumis;

	if (fChain == 0) return;
	Long64_t nentries = fChain->GetEntriesFast();
	cout << "nentries = " << nentries << endl << endl;			//information

	Long64_t nbytes = 0, nb = 0;
	///for (Long64_t jentry=0; jentry<nentries;jentry++) {}
	for (Long64_t jentry=0; jentry!=60; ++jentry){ //for testing with less entries/events

		Long64_t ientry = LoadTree(jentry);
		cout << "current entry: " << jentry << endl;
		if (jentry%10==0){cout << "------------------------------------------" << endl;}//for overview
		///if (jentry%100000==0){cout << "processed " << jentry << " of " << nentries << endl;}
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);   nbytes += nb;
		// if (Cut(ientry) < 0) continue;

		//set<int,int> cur_set = make_pair(luminosityBlock, jentry); //create pair of event index and corresponding luminosityBlock
		//lumis.insert(make_pair(run, cur_set));
		////lumis.insert(make_pair(run, make_pair(luminosityBlock, jentry)));
		int cur_lumisec = luminosityBlock;
		auto cur_pair = make_pair(cur_lumisec, jentry);
		int cur_run = run;
		/////lumis.insert(pair<int, set<int, int>>(cur_run, cur_pair));
		lumis[run].insert(make_pair(luminosityBlock, jentry));
		//cout << "lumis[run]: " << lumis[run].first << endl;


   }//event loop
	cout << "out of event loop" << endl;

	//now write contents of lumis (map) into .txt file:
	map<int, set<pair<int, int>>>::iterator it = lumis.begin();
	while(it != lumis.end()){
		//cout << "run number: " << it->first << " lumisec: " << it->second << endl;
		cout << "run" << endl;
		
		//within run-number-iteration, go through set of lumisections:
		//set<pair<int, int>>::iterator itlumi = (lumis->second).begin();
		set<pair<int, int>> curset = it->second;
		/////set<pair<int, int>>::iterator itlumi = it->second.begin();
		set<pair<int, int>>::iterator itlumi = curset.begin();

		txt_file << "--------------------" << endl;
		///set<pair<int, int>>::iterator itlumi = lumis.second.begin();
		//////while(itlumi != it->second.end()){
		while(itlumi != curset.end()){

		///while(itlumi != lumis.second.end()){
			//cout << "lumisection: " << itlumi->first << " event index: " << itlumi->second << endl;
			cout << "it " << Form("%d", itlumi) << endl;
			txt_file << "test: itlumi->first:: " << itlumi->first << "\n";
			break;
		}//while-loop for lumisection (inner loop)
		break;
	}//end writing-while-loop

	txt_file.close();
	return 0;
}

