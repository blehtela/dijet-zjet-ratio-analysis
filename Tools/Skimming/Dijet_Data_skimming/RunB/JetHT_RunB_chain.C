#define JetHT_RunB_chain_cxx
#include "JetHT_RunB_chain.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

using namespace std;

// Created via MakeClass out of a TChain containing all the JetHT files for RunB
// This script is mainly used for skimming
// Skimming: Loop through events: 	> Check if nJet>2 --> if yes: keep, if no: discard
//									> Check if ANY of the jet-triggers fired --> if yes: keep, if no: discard



void JetHT_RunB_chain::Loop()
{
//   In a ROOT session, you can do:
//      root> .L JetHT_RunB_chain.C
//      root> JetHT_RunB_chain t
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


//function used for skimming the data set	(completely analogue to Run A skimming)
//only adjustment: Name of data set and output file etc.
void JetHT_RunB_chain::Skim(){

	//create output file
	TFile *skimfile = new TFile("JetHT_RunB_skimmed.root", "recreate");

	//force 1st tree to be loaded
	fChain->LoadTree(0);

	//clone the original TTree
	TTree *newtree = fChain->GetTree()->CloneTree(0);

	Int_t nevents = (Int_t)fChain->GetEntries();
	Int_t nbytes = 0;
	Int_t nselected = 0;

	//all branches are switched on. (Addresses set in .h file)

	//some counters
	int count_nJetFail = 0;			//count events with less than two jets
	int count_NoJetTrigger = 0;	//count events that dropped out because none of the HLT jet triggers has fired

	//create a histogram containing meta-data
	const char* dataset_title_1 = "/JetHT/Run2018B-Nano1June2019-v2/NANOAOD";
	const char* dataset_title_2 = "";
	const char* dataset_title_3 = "";
	//TH1D *metahist = new TH1D("metahist", Form("#splitline{%s}{#splitline{%s}{%s}}; status(#bf{0: failed} skim, #bf{1: passed} skim); Events", dataset_title_1, dataset_title_2, dataset_title_3), 2, -0.5, 1.5);
	TH1D *metahist = new TH1D("metahist", Form("%s; status(#bf{0: failed} skim, #bf{1: passed} skim); Events", dataset_title_1), 2, -0.5, 1.5);
	//--> should more bins be introduced? like "passed", "failed (nJet<2)", "failed (trigger)"

	metahist->GetXaxis()->SetNdivisions(2);
	metahist->SetLineWidth(2);
	metahist->SetLineColor(kGreen+2);


	//go through all the events (of all the files, in the TChain)
	//and do some sensible pre-selection
	for(int ievent=0; ievent!=nevents; ++ievent){
	///for(int ievent=0; ievent!=100000; ++ievent){}	//for testing
		nbytes += fChain->GetEntry(ievent);

		//for progress overview
		if(ievent%1000000==0){	//informs every millionths event
			cout << "]" << endl;
			cout << "Processed " << ievent << " of " << nevents << " events." << endl;
			cout << "[";
		}
		if(ievent%100000==0){cout << "#"; cout.flush();}	//progress bar grows after every 100000st event


		//throw away events with less than two jets.
		if(nJet<2){
			count_nJetFail++;	//count that this event dropped out, because of too few jets
			metahist->Fill(0);	//increse the entries of failed events in the metadata histogram
			continue;			//do selection (by skipping such events)
		}//check if less than 2 jets

		//check if current event passed any of the HLT jet triggers:
		bool passed_trg = (HLT_PFJet40 or HLT_PFJet60 or HLT_PFJet140 or HLT_PFJet200 or HLT_PFJet260 or HLT_PFJet320 or HLT_PFJet400 or HLT_PFJet450 or HLT_PFJet500 or HLT_PFJet550);

		if(nJet>1){
			if(passed_trg){
				nselected++;		//count this event (it has passed, will be written to file)
				metahist->Fill(1);	//increase the entries of passed events in the metadata histogram
				newtree->Fill();	//fill newtree
			}
			else{
				count_NoJetTrigger++;	//this event has not passed any trigger
				metahist->Fill(0);		//increase the entries of failed events in the metadata histogram
				continue;				//continue with next event in event-loop
			}//hasn't passed any trigger

			
		}//dijet case
	}//end of event loop


	//Style and Write, Output summary
	gStyle->SetTitleAlign(23);

	//prepare an informative canvas containing metadata
	//will be stored in the .root file
	TCanvas *infocanv = new TCanvas("infocanv", "infocanv", 600, 400);
	TLatex *metainfo = new TLatex();
	metainfo->SetNDC();
	metainfo->SetTextSize(0.045);
	metainfo->DrawLatex(0.06, 0.92, Form("Dataset:  #scale[1.2]{#color[30]{/JetHT/Run2018B-Nano1June2019-v2/NANOAOD}}"));
	metainfo->DrawLatex(0.06, 0.80, Form("Processed #scale[1.1]{#color[38]{%d}} events.", nevents));
	metainfo->DrawLatex(0.06, 0.74, Form("----------------------------------------------"));
	metainfo->DrawLatex(0.06, 0.68, Form("#scale[1.1]{#color[8]{%d}}  events selected in skim. (= passed)", nselected));
	metainfo->DrawLatex(0.06, 0.56, Form("#scale[1.1]{#color[46]{%d}}  events dropped out due to nJet < 2.", count_nJetFail));
	metainfo->DrawLatex(0.06, 0.46, Form("#scale[1.1]{#color[46]{%d}}  events dropped out, because they did not pass any HLT jet trigger.", count_NoJetTrigger));
	metainfo->DrawLatex(0.06, 0.40, Form("#scale[0.8]{Checked triggers: HLT_PFJet#color[42]{XXX} (with #color[42]{XXX} = #color[42]{40, 60, 80, 140, 200, 260, 320, 400, 450, 500, 550})}"));
	infocanv->Write();
	infocanv->Close();
	
	metahist->Write();
	newtree->Write(); 	//write the new "combined" tree, containing the skimmed files
	skimfile->Close();	//close file

	//direct info: output for user
	cout << endl;
	cout << "Summary: " << endl;
	cout << "-----------------" << endl;
	cout << "Selected " << nselected << " events." << endl;
	cout << count_nJetFail << " events dropped out due to nJet<2." << endl;
	cout << count_NoJetTrigger << " events dropped out, because they did not pass any jet trigger." << endl;
	cout << "Output .root file: " << skimfile->GetName() << endl;
}//end Skim()



