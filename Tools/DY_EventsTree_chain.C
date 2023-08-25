#define DY_EventsTree_chain_cxx
#include "DY_EventsTree_chain.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

using namespace std;
bool trigselect=false;	//if true: select via HLT_Mu23_Mu12 trigger. --> otherwise only use pT cuts --> no trigger!
//bool trigselect=true;	//just for this test

void DY_EventsTree_chain::Loop(){
//   In a ROOT session, you can do:
//      root> .L DY_EventsTree_chain.C
//      root> DY_EventsTree_chain t
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

	//Long64_t nentries = fChain->GetEntriesFast(); //original version, but yields to a large number..
	Long64_t nentries = fChain->GetEntries();

	std::cout << "Entries: " << nentries << std::endl;
	fChain->SetBranchStatus("*", 0);
	fChain->SetBranchStatus("Jet_pt", 1);
	fChain->SetBranchStatus("nMuon", 1);
	fChain->SetBranchStatus("Muon_pt",1);

	Long64_t nbytes = 0, nb = 0;
	for(Long64_t jentry=0; jentry<nentries;jentry++){
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);   nbytes += nb;
		// if (Cut(ientry) < 0) continue;

		if(jentry==200){
			std::cout << "Jet_pt of leading jet in event 200: " << Jet_pt[0] << std::endl;
		}

		//to test how muons are ordered in array
		if(jentry%1000000==0){
			if(nMuon>2){
				std::cout << Form("jentry: %lld --> Muon_pt[0]=%f, Muon_pt[1]=%f, Muon_pt[2]=%f", jentry, Muon_pt[0], Muon_pt[1], Muon_pt[2]) << endl;
			}
			else if(nMuon==2){
				std::cout << Form("jentry: %lld --> Muon_pt[0]=%f, Muon_pt[1]=%f", jentry, Muon_pt[0], Muon_pt[1]) << endl;
			}
			else{
				std::cout << Form("jentry: %lld --> Less than two muons.", jentry) << endl;
			}
		}
	}
}


void DY_EventsTree_chain::Skim_DY(){
	//make a TChain
	//use fChain coming from DY_EventsTree_chain.h file!
	//are they already all appended to this fChain due to the MakeClass call, that I made? -- yes.

	//create output file with combined skimmed MC-set
	TFile *skimfile = new TFile("DY_MC_skimmed.root", "recreate");
	//TFile *skimfile = new TFile("DY_MC_skimmed_TRIGGER_test.root", "recreate");
	//TFile *skimfile = new TFile("DY_MC_skimmed_DIMUON.root", "recreate");


	//force 1st tree to be loaded
	fChain->LoadTree(0);

	//clone the original TTree
	TTree *newtree = fChain->GetTree()->CloneTree(0);

	//just a check:
	//newtree->Print("p");

	Int_t nevents = (Int_t)fChain->GetEntries();
	Int_t nbytes = 0;
	Int_t nselected = 0;

	//set branch status (switch on branches that are used for the selection)
	//---> want ALL the branches as well in resulting .root file, so nothing switched off here!
	/*
	fChain->SetBranchStatus("*", 0);			//switch off all branches
	fChain->SetBranchStatus("nMuon", 1);		//switch on what is needed for skimming
	fChain->SetBranchStatus("Muon_mediumPromptId", 1);
	fChain->SetBranchStatus("HLT_Mu23_Mu12", 1);
	*/

	//some counters
	int dropout_nMuon = 0;
	int dropout_prompt_trigger = 0;

	//one histogram with meta-data
	//const char* dataset_title_2 = "/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18NanoAODv5-Nano1June2019_102X_upgrade2018_realistic_v19-v1/NANOAODSIM";

	const char* dataset_title_1 = "/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/";
	const char* dataset_title_2 = "RunIIAutumn18NanoAODv5-Nano1June2019_102X_upgrade2018_realistic_v19-v1/";
	const char* dataset_title_3 = "NANOAODSIM";
	TH1D *metahist = new TH1D("metahist", Form("#splitline{%s}{#splitline{%s}{%s}}; status(#bf{0: failed} due to #it{nMuon < 2}, #bf{1: passed}); Events", dataset_title_1, dataset_title_2, dataset_title_3), 2, -0.5, 1.5);
	metahist->GetXaxis()->SetNdivisions(2);
	metahist->SetLineWidth(2);
	metahist->SetLineColor(kGreen+2);

	//go through all the events (of all the files, in the TChain)
	//and do some sensible pre-selection
	for(int ievent=0; ievent!=nevents; ++ievent){
	///for(int ievent=0; ievent!=10000; ++ievent){}	//for testing
		nbytes += fChain->GetEntry(ievent);

		//for progress overview
		if(ievent%1000000==0){	//informs every millionths event
			cout << "]" << endl;
			cout << "Processed " << ievent << " of " << nevents << " events." << endl;
			cout << "[";
		}
		if(ievent%100000==0){cout << "#"; cout.flush();}	//progress bar grows after every 100000st event


		//throw away events with less than two muons!
		if(nMuon<2){
			dropout_nMuon++;	//count that this event dropped out, because of too few muons
			metahist->Fill(0);	//increase the entries of failed events in the metadata histogram
			continue; //do selection (by skipping such events)
		}//check if less than 2 muons


		if(nMuon>1){

			if(trigselect==true){ //usually false
				//check if this event has at least 2 muons that pass mediumPromptId
				Int_t n_muprompt = 0;
				for(int k=0; k!=nMuon; ++k){
					if(Muon_mediumPromptId[k]==true){
						n_muprompt += 1;		//count muon that passes requirement
					}
				}

				//check how many events with at least two muons drop out due to MediumPromptId selection:
				// ...or trigger
				bool muon_trigger = HLT_Mu23_Mu12;

				//throw away events that either have less than two prompt muons and/or do not pass trigger
				if((n_muprompt<2)||(muon_trigger==false)){
					dropout_prompt_trigger++;
					continue; //skip this event!
				}
			}	//trigger for skimming --> usually NOT USED when working on MC!
				//checking muprompt and mu-trigger

			/*
			//starting to throw away events where the muons have too low pt
			//assume that muons are sorted according to their pt (0-->highest pt, 1, 2,...,nMuon-1)
			for(int mu1_ind=0; mu1_ind!=(nMuon-1); ++mu1_ind){
				for(int mu2_ind=(mu1_ind+1); mu2_ind!=nMuon; ++mu2_ind){
					//find index of hardest and 2nd hardest muon
				}//inner loop through muons
			}//outer loop through muons

			//check the two hardest muons:
			if((Muon_pt[mu1]<23)||(Muon_pt[mu2]<12)){continue;}
			if((Muon_pt[0]<23)||(Muon_pt[1]<12)){continue;}	//working directly with indices 0 and 1
			*/
	

		}//selecting via pT of muon1 and muon2

		nselected++;
		metahist->Fill(1);	//increase the entries of passed events in the metadata histogram
		newtree->Fill();	//fill newtree
	}

	gStyle->SetTitleAlign(23);
	metahist->Write();
	newtree->Write(); 	//write the new "combined" tree, containing the skimmed files
	skimfile->Close();	//close file
	cout << endl;
	cout << "Summary: " << endl;
	cout << "-----------------" << endl;
	cout << "Selected " << nselected << " events." << endl;
	cout << dropout_nMuon << " events dropped out due to nMuon<2." << endl;
	if(trigselect==true){
		cout << dropout_prompt_trigger << " events dropped out due to #(Muon_mediumPropmtId)<2 and/or (HLT_Mu23_Mu12)==false." << endl << endl;
	}
	cout << "Output .root file: " << skimfile->GetName() << endl;

}//end Skim_DY()
