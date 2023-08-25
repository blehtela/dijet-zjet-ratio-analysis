#define Skim_DY_cxx
#include "ZjetAnalysis.h"	//contains classes declared via MakeClass()
//#include "ZjetClasses.h"	//contains Muon-, YbYs-, Histos- classes
#include <string>

using namespace std;


void Skim_DY(){
	//make a TChain
	TChain *chain = new TChain("ChainTTree");

	//append all the desired 53 files to the chain
	const char* datadir = "/Volumes/CMSdata/RunIIAutumn18NanoAODv5/DYJetsToLL_M-50_TuneCP5";
	chain->Add(Form("%s/05884C27-75AD-D340-B515-7017F9655675.root", datadir)); //#1
	chain->Add(Form("%s/0CA4B9C4-805D-C148-8281-D615F9DE8541.root", datadir)); //#2
	chain->Add(Form("%s/0F49C966-5F44-3D4F-AADF-F820A2EBF8A9.root", datadir)); //#3
	chain->Add(Form("%s/12C1D5AD-DFFB-F547-A634-17FE8AAB84B1.root", datadir)); //#4
	chain->Add(Form("%s/168D358A-B3B2-6849-9EF4-D2B6791A26AA.root", datadir)); //#5
	chain->Add(Form("%s/16B6B7CD-4310-A042-AB52-7DA8ADA22922.root", datadir)); //#6
	chain->Add(Form("%s/1A9BA6F1-F51D-F342-BB5D-F0F3B17ED70E.root", datadir)); //#7
	chain->Add(Form("%s/1C3AC8F7-987B-4D40-B002-767A2C65835B.root", datadir)); //#8
	chain->Add(Form("%s/26884FA0-B96A-1745-AA11-597C5168EF5E.root", datadir)); //#9
	chain->Add(Form("%s/274599AC-1636-3641-B09F-ECA42B8F63A4.root", datadir)); //#10

	chain->Add(Form("%s/2A9A7EDE-2249-2C44-AF6D-E44B83E8CBDF.root", datadir)); //#11
	chain->Add(Form("%s/3C0F69F9-2D31-6646-A1B0-FE021BE707C8.root", datadir)); //#12
	chain->Add(Form("%s/3C6B053B-C36F-BA4B-8A44-1BC2C291DC8D.root", datadir)); //#13
	chain->Add(Form("%s/425C6243-B617-5347-B015-F0908D757828.root", datadir)); //#14
	chain->Add(Form("%s/43744293-FC6D-D24E-8BB5-DEB0F975AAF5.root", datadir)); //#15
	chain->Add(Form("%s/4966A0E8-CC00-BA4D-AEFF-F37332AC1096.root", datadir)); //#16
	chain->Add(Form("%s/4989E9FC-CD00-D64C-8527-9AD2FF2546A8.root", datadir)); //#17
	chain->Add(Form("%s/49CF353A-3303-ED44-8A65-994191042C99.root", datadir)); //#18
	chain->Add(Form("%s/51F7E9BD-25F0-D941-BD9A-7F24FD6B1B04.root", datadir)); //#19
	chain->Add(Form("%s/572AD18F-BFA0-174E-9D28-0F75F965C147.root", datadir)); //#20

	chain->Add(Form("%s/5851642F-9C54-FC4E-A8B0-DEA3A10646E4.root", datadir)); //#21
	chain->Add(Form("%s/5AE16A97-FA54-0C45-8EDC-C1AA89D5B054.root", datadir)); //#22
	chain->Add(Form("%s/60C68FA7-5C29-2740-8242-E519DA6F5F10.root", datadir)); //#33
	chain->Add(Form("%s/6827CFDC-89DC-9147-98F8-B25D649B5857.root", datadir)); //#24
	chain->Add(Form("%s/698C3A61-1ADA-0446-A47D-D8192E2B5415.root", datadir)); //#25
	chain->Add(Form("%s/6BAFEC1D-F9BD-8C48-8A2C-D00229113414.root", datadir)); //#26
	chain->Add(Form("%s/761BC07D-C19C-C841-84E4-2766F1DCB60B.root", datadir)); //#27
	chain->Add(Form("%s/7AC63348-51B1-9649-B824-D04898C5BA5B.root", datadir)); //#28
	chain->Add(Form("%s/7BF6D85F-EE22-D840-A435-8CB817098E86.root", datadir)); //#29
	chain->Add(Form("%s/80D5E58B-2B5C-0440-9A1F-B6E2772FD7BD.root", datadir)); //#30

	chain->Add(Form("%s/81FFF806-71B3-CC44-AB43-714DBE4C9319.root", datadir)); //#31
	chain->Add(Form("%s/8E443B60-F9E9-5444-B37D-62902E68C0C3.root", datadir)); //#32
	chain->Add(Form("%s/8F3EEF08-F61E-4046-B140-B04B87602708.root", datadir)); //#33
	chain->Add(Form("%s/8FA629F5-385A-AD4A-BB6F-D0856E633712.root", datadir)); //#34
	chain->Add(Form("%s/932CE866-A30E-F34D-B0D5-4C4CEAA06CB8.root", datadir)); //#35
	chain->Add(Form("%s/948182F2-9993-C74D-B2EA-1D6E0098AD61.root", datadir)); //#36
	chain->Add(Form("%s/9F70ACE0-A9C2-494C-B0E5-42E7017ABF95.root", datadir)); //#37
	chain->Add(Form("%s/A1B3E169-6D65-E44E-B891-8F738CBB78AD.root", datadir)); //#38
	chain->Add(Form("%s/A5702444-A58D-364F-BF6C-EF28C9C52344.root", datadir)); //#39
	chain->Add(Form("%s/AB329578-42CC-4746-A15D-08E70CD2554E.root", datadir)); //#40

	chain->Add(Form("%s/AF265BB7-CF6C-8241-8DC2-F13BA8A9AD60.root", datadir)); //#41
	chain->Add(Form("%s/AF34E3F0-25B7-6644-B557-1428CF675FDC.root", datadir)); //#42
	chain->Add(Form("%s/B08CCF2F-A193-9640-AAE4-82F75CE2B436.root", datadir)); //#43
	chain->Add(Form("%s/B4D431A5-B57B-1D40-89AD-B6CB0C9F05E2.root", datadir)); //#44
	chain->Add(Form("%s/BB0AF882-527A-D641-91BE-4624222CCB17.root", datadir)); //#45
	chain->Add(Form("%s/C1450DF7-55A8-FE4F-B50B-9E844F9E103E.root", datadir)); //#46
	chain->Add(Form("%s/C1F8DBF7-4649-2F43-9E03-2161344D4663.root", datadir)); //#47
	chain->Add(Form("%s/C3DC24F1-8416-0240-AC54-D63A9760BF45.root", datadir)); //#48
	chain->Add(Form("%s/CFC96CE1-3B0B-3D45-8D85-C73EE0C8DAB3.root", datadir)); //#49
	chain->Add(Form("%s/DC937EFE-252C-814B-A4B2-743A368793D8.root", datadir)); //#50

	chain->Add(Form("%s/E3828699-7905-3142-A0A2-929E60406883.root", datadir)); //#51
	chain->Add(Form("%s/F690EE89-E028-C840-8C4B-3A3106316530.root", datadir)); //#52
	chain->Add(Form("%s/FC56B1DA-20B9-F14A-A2CF-2097B8095BEB.root", datadir)); //#53


	//set branch addresses
	chain->SetBranchAddress("Run", &myRun);
	chain->SetBranchAddress("Event", &myEvent);

	//set branch status (switch on branches that are used for the selection)
	//chain->SetBranchStatus();

	//create output file with combined skimmed MC-set
	TFile *skimfile = new TFile("DY_MC_skimmed.root", "recreate");

	//force 1st tree to be loaded
	chain->LoadTree(0);

	//clone the original TTree
	TTree *newtree = chain->GetTree()->CloneTree(0);

	//just a check:
	//newtree->Print("p");

	Int_t nevents = (Int_t)chain->GetEntries();
	Int_t nbytes = 0;
	Int_t nselected = 0;

	//set branch status (switch on branches that are used for the selection)
	chain->SetBranchStatus("*", 0);			//switch off all branches
	chain->SetBranchStatus("nMuon", 1);		//switch on what is needed for skimming
	chain->SetBranchStatus("Muon_mediumPromptId", 1);
	chain->SetBranchStatus("HLT_Mu23_Mu12", 1);

	//go through all the events (of all the files, in the TChain)
	//and do some sensible pre-selection
	for(int ievent=0; ievent!=nevents; ++ievent){
		nbytes += chain->GetEntry(ievent);
		if(nMuon<2){continue;}	//do selection (by skipping such events)
		if(nMuon>1){
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
			if((n_muprompt<2)||(muon_trigger==false)){continue;} //skip this event!
		}//checking muprompt and mu-trigger

		nselected++;
		newtree->Fill();	//fill newtree
	}

	newtree->Write(); 	//write the new "combined" tree, containing the skimmed files
	skimfile->Close();	//close file

}//end Skim_DY()
