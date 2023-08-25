#define DoubleMuon_RunA_chain_cxx
#include "DoubleMuon_RunA_chain.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

using namespace std;

// Created via MakeClass out of a TChain containing all the DoubleMuon files for RunA
// This script is mainly used for skimming
// Skimming: Loop through events: 	> Loop through all the events
//									>



void DoubleMuon_RunA_chain::Loop()
{
//   In a ROOT session, you can do:
//      root> .L DoubleMuon_RunA_chain.C
//      root> DoubleMuon_RunA_chain t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
   }
}

//function used for skimming the data set:
void DoubleMuon_RunA_chain::Skim(){
	//--------------------------//
	//	Preparing the Output	//
	//--------------------------//

	//create output file
	TFile *skimfile = new TFile("DoubleMuon_RunA_skimmed.root", "recreate");

	//forse 1st tree to be loaded
	fChain->LoadTree(0);

	//clone the original TTree
	TTree *newtree = fChain->GetTree()->CloneTree(0);

	Int_t nevents = (Int_t)fChain->GetEntries();
	Int_t nbytes = 0;
	Int_t nselected = 0;

	//all branches are switched on. (Addresses set in .h file)

	//some counters
	int count_nMuonFail = 0;		//count events with less than two muons (should not exist in this sample...)
	

	//create histogram containing meta-data
	const char* dataset_title = "/DoubleMuon/Run2018A-Nano1June2019-v1/NANOAOD";
	TH1D *metahist = new TH1D("metahist", Form("%s; status(#bf{0: failed} skim, #bf{1: passed} skim); Events", dataset_title), 2, -0.5, 1.5);

	metahist->GetXaxis()->SetNdivisions(2);
	metahist->SetLineWidth(2);
	metahist->SetLineColor(kGreen+2);

	//--------------//
	//	EVENT LOOP	//
	//--------------//
	//go through all the events (of all the files, in the TChain)
	//and do some sensible pre-selection
	///for(int ievent=0; ievent!=nevents; ++ievent){}
	for(int ievent=0; ievent!=100000; ++ievent){	//for testing
		nbytes += fChain->GetEntry(ievent);

		//for progress overview
		if(ievent%1000000==0){	//informs every millionths event
			cout << "]" << endl;
			cout << "Processed " << ievent << " of " << nevents << " events." << endl;
			cout << "[";
		}
		if(ievent%100000==0){cout << "#"; cout.flush();}	//progress bar grows after every 100000st event


		//throw away events with less than two muons
		if(nMuon<2){
			count_nMuonFail++;		//count that this event dropped out, because of too few muons
			metahist->Fill(0);		//increase the entries of failed events in the metadata histogram
			continue;				//do selection (by skipping such events)
		}//check if less than two muons.

		else{//if nMuon>1
			nselected++;			//count this event (it has passed, will be written to file)
			metahist->Fill(1);		//increase the entries of passed events in the metadata histogram
			newtree->Fill();		//fill newtree
		}//nMuon>1
	}//end of event loop


	//--------------//
	//	Finalising	//
	//--------------//
	//Style and Write, Output summary
	gStyle->SetTitleAlign(23);	//?? does this have any impact??

	//prepare an informative canvas containing metadata
	//will be stored in the .root file
	TCanvas *infocanv = new TCanvas("infocanv", "infocanv", 600, 400);
	TLatex *metainfo = new TLatex();
	metainfo->SetNDC();
	metainfo->SetTextSize(0.045);
	metainfo->DrawLatex(0.06, 0.92, Form("Dataset:  #scale[1.2]{#color[30]{/DoubleMuon/Run2018A-Nano1June2019-v1/NANOAOD}}"));
	metainfo->DrawLatex(0.06, 0.80, Form("Processed #scale[1.1]{#color[38]{%d}} events.", nevents));
	metainfo->DrawLatex(0.06, 0.74, Form("----------------------------------------------"));
	metainfo->DrawLatex(0.06, 0.68, Form("#scale[1.1]{#color[8]{%d}}  events selected in skim. (= passed)", nselected));
	metainfo->DrawLatex(0.06, 0.56, Form("#scale[1.1]{#color[46]{%d}}  events dropped out due to nMuon < 2.", count_nMuonFail));
	///metainfo->DrawLatex(0.06, 0.46, Form("#scale[1.1]{#color[46]{%d}}  events dropped out, because they did not pass any HLT jet trigger.", count_NoJetTrigger));
	///metainfo->DrawLatex(0.06, 0.40, Form("#scale[0.8]{Checked triggers: HLT_PFJet#color[42]{XXX} (with #color[42]{XXX} = #color[42]{40, 60, 80, 140, 200, 260, 320, 400, 450, 500, 550})}"));
	infocanv->Write();
	infocanv->Close();
	
	metahist->Write();
	newtree->Write(); 	//write the new "combined" tree, containing the skimmed files (as one file)
	skimfile->Close();	//close file

	//direct info: output for user
	cout << endl;
	cout << "Summary: " << endl;
	cout << "-----------------" << endl;
	cout << "Selected " << nselected << " events." << endl;
	cout << count_nMuonFail << " events dropped out due to nMuon<2." << endl;
	//cout << count_NoJetTrigger << " events dropped out, because they did not pass any jet trigger." << endl;
	cout << "Output .root file: " << skimfile->GetName() << endl;

}//end Skim()
