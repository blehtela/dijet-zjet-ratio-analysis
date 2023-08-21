#define DijetAnalysis_cxx
#include "DijetAnalysis.h"
#include "FriendClass.h"		//for unc analysis
#include <string>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "tdrstyle_mod14.C"
#include <TMath.h>
#include "Math/LorentzVector.h"

//for JEC stuff
////#include "../../../jecsys/CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"		//should be absolute path?
////#include "../../../jecsys/CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

//#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h" 						//--> like this
//#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"


//for get_lumi_txt()
#include <iostream>
#include <fstream>

using namespace std;
using namespace ROOT::Math;

//TH1::SetDefaultSumw2(kTRUE);		//in order to ALWAYS store the bin errors

// Analysis of inclusive jet data
// --> Dijet Analysis
// Using pTavg, yboost, ystar

//verbose mode?
bool v = false;

//MONTE CARLO SWITCH (montecarlo==true) is already set in header file DijetAnalysis.h


// current choice for yboost / ystar bins:
// make folder for each of the six yboost-ystar-bins --> could of course be arranged differently

//--------------------------------------------------------------------------------------------------//
//											 CLASSES 												//
//--------------------------------------------------------------------------------------------------//

//class for booking and filling the various histograms (considering trigger decisions)
class Histos;

//class for handling yboost ystar binning
class YbYs {
	friend class Histos;
	private:
		string bname;			//bin name, e.g. "yb0ys0" for 0.0 <= yboost, ystar < 0.5, assigned in constructor
		//const char *histoname;	//not really needed if one leaves out the initial (triggerless) histogram
		float yb_low;			//lower yboost bin bound
		float yb_up;			//upper yboost bin bound
		float ys_low;			//lower ystar bin bound
		float ys_up;			//upper ystar bin bound

		TH1D* mc_histo;				//histogram for MC analysis
		TH1D* mc_histo_smearedJet;	//histogram for MC HERWIG 7 unc analysis

		TH1D* mc_genhisto;		//histogram for MC analysis --> generated spectrum
		TH2D* mc_genrec_noGenCuts;	//histogram as originally implemented
		TH2D* mc_gen_reco_ptav;		//histogram for MC analysis (for pt unfolding)
		TProfile* mc_profile;		//TProfile plot for checking reco and gen distributions
		TProfile* mc_profile_j1;	//TProfile for jet1 pt
		TProfile* mc_profile_j2;	//TProfile for jet2 pt
		TProfile* mc_profile_mergedj1j2;	//TProfile filled with jet1 pt and jet2 pt

		//stuff for herwig jme unc:
		TH1D* jt_ptraw_hist;
		TH1D* jt_ptnom_hist;
		TH1D* jt_corrJEC_hist;
		TH1D* jt_corrJER_hist;
		TH1D* jt_ptjerUp_hist;
		TH1D* jt_ptjerDw_hist;
		
		TH1D *sumhisto;			//summary histo for this ybys bin --> filled after trigger-analysis
		map<const char*, TH1D*> trglumi_hist_map;	//used in DijetAnalysis::triggerpaths() for accessing the ratio-histos --> not initialised per default --> must do this if needed!
		const char* ybys_dirname;	//because directory version did not work... --> pointer problems?
		//TDirectory *bdir;

	public:
		YbYs (float, float, float, float);					//constructor with lower and upper bin bound for yboost ystar
		YbYs (float, float, float, float, vector<string>);	//overloaded constructor for uncertainty histos
		//YbYs();								//new constructor --> used in normalise_mc() function, as only the folder structure is needed, not all the histos
		//YbYs (float, float);				//constructor with only LOWER bin bound for yb and ys
		pair <float, float> yb_get_bb() {return make_pair(yb_low,  yb_up);};		//get yboost bin bounds
		pair <float, float> ys_get_bb() {return make_pair(ys_low,  ys_up);};		//get yboost bin bounds
		string get_bname() {return bname;};		//to read bin name
		//TDirectory get_dir(){return bdir;};	//to get ybys directory
		const char* get_dirname(){return ybys_dirname;};	//getting name of directory
		map<string, TH1D*> unc_histos_up;			//map of unc histos (variation upwards)--> key=string uncname, mapped value=TH1D* histogram
		map<string, TH1D*> unc_histos_down;			//map of unc histos, downwards

		TH1D* get_mc_histo(){return mc_histo;}							//get the histogram with the MC
		TH1D* get_mc_histo_smearedJet(){return mc_histo_smearedJet;}	//only used in HERWIG 7 with smeared jet (Jet_pt_nom branch)

		TH1D* get_mc_genhisto(){return mc_genhisto;}		//get the histogram with the MC genjet spectrum
		TH2D* get_mc_genreco_noGenCuts(){return mc_genrec_noGenCuts;}	//get the histogram with gen-reco-ptav distribution, where cuts (pt, eta) have only be applied to RECO
		TH2D* get_mc_genreco(){return mc_gen_reco_ptav;}	//get the histogram with gen-reco-ptav distribution
		TProfile* get_mc_profile(){return mc_profile;}		//get the profile plot for ptavg
		TProfile* get_mc_profile_j1(){return mc_profile_j1;}	//get the profile plot for jet1
		TProfile* get_mc_profile_j2(){return mc_profile_j2;}	//get the profile plot for jet2
		TProfile* get_mc_profile_mergedj1j2(){return mc_profile_mergedj1j2;}	//get the profile plot for jet1-jet2-combi


		map<const char*, TH1D*> getmap(){return trglumi_hist_map;};	//return map where one can fill in ratio plots (used with DijetAnalysis::triggerpaths())
		void ybys_directories(TDirectory *standard_folder, vector<int> triggers, Histos curtrig, bool trgmode); //creates folder-hierarchy
		//void ybys_directories(TDirectory *standard_folder);	//creates folder-hierarchy in case of montecarlo!
		TDirectory* ybys_directories(TDirectory *standard_folder);	//creates folder-hierarchy in case of montecarlo! (version that returns TDirectory*)
		void ybys_fastdirs(TDirectory *standard_folder, vector<int> triggers, Histos curtrig); //for FastLoop()
		TDirectory* ybys_mkdir(TDirectory *standard_folder);		//only creates the ybys directory and returns pointer to it

		void ybys_uncdirs(TDirectory *unc_folder, vector<string> uncnames);	//used to create uncertainty histograms in Uncertainty folder
		void ybys_jmedirs(TDirectory *jme_folder);	//puts 6 jme histos in the additional jme uncs folder

		void trglumi_directories(TDirectory *standard_folder, vector<int> triggers, map<const char*, TH1D*> trglumi_hist_map); //used in function triggerpaths()
		//void fill_sumhisto(vector<Double_t> effvec, Histos curtrig);	//function for filling the summary histogram of this ybys bin
		TH1D* fill_sumhisto(map<int, Double_t> effmap, vector<int> triggers, Histos *curtrig, vector<Double_t> lumivec);//see above
		TH1D* fill_fastsumhisto(vector<int> triggers, Histos *curtrig, vector<Double_t> lumivec);	//used in DijetAnalysis::FastLoop() (simplified trigger handling)
};

//class for booking the histograms and creating trigger-directories
class Histos{
	friend class YbYs;
	private:
		//TH1D *jt_hist;
		//map<const int, TH1D*> jt_histos;		//map of trigger histograms --> key=int pt, mapped value=TH1D *histogram
		map<const int, bool*> triggers_bool;	//containing the trigger decisions	--> still has to be put in use!
		//vector<const char*> hjt_names;
		vector<const char*> trigger_names;		//containing names of the TBranches
	public:
		Histos ();				//default constructor
		Histos (vector<int>);	//Constructor with vector containing different triggers (pt(GeV) in int)
		Histos (vector<int> triggers, TFile *histo_file, const char* ybysname);	//overloaded constructor, used in triggerpaths() function for getting filled histos

		vector<const char*> hjt_names; //just for testing, usually private
		map<const int, TH1D*> jt_histos;		//map of trigger histograms --> key=int pt, mapped value=TH1D *histogram
		map<const int, TH1D*> trgobj_histos;	//map of "simulated" trg histos (combining previous HLT_PFJet trigger with current TrigObj_pt information)
												//used in order to show equivalence of both trigger fitting methods
		map<const int, TH1D*> n2_trgobj_histos;	//simulated "higher" triggers for the N/N-2 case --> this is the HLT_N  && TrigObj_pt[N+2] case, where trigger index = N+2
		//void trigger_histos(Etabins etabin, TH1D *jt_hist);	//function for creating and filling the histos associated with triggers --> currently not used
};



//--------------------------------------------------------------------------------------------------//
//											 CONSTRUCTORS											//
//--------------------------------------------------------------------------------------------------//
//Constructor for ybys bins without uncertainty histograms
YbYs::YbYs (float yboost_lo, float ystar_lo, float yboost_up, float ystar_up) {
	yb_low = yboost_lo;
	ys_low = ystar_lo;
	yb_up = yboost_up;	
	ys_up = ystar_up;

	bname = Form("yb%01.0fys%01.0f", yb_low, ys_low);
	if(v==true){cout << Form("Created yboost-ystar-bin with: %f <= yb < %f and %f <= ys < %f", yb_low, yb_up, ys_low, ys_up) << endl;}

	//customised pT bin bounds
	//the following are the usual:
	Double_t cust_bb[] = {10, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507, 548, 592, 638, 686, 737, 790, 846, 905, 967, 1032, 1101, 1172, 1248, 1327, 1410, 1497, 1588, 1684, 1784, 1890, 2000, 2116, 2238, 2366, 2500, 2640, 2787, 2941, 3103, 3273, 3450, 3832, 6076, 6389};

	/*
	//this are a little more bins in the central GeV region --> for triggerpaths()
	const int ncust_bb = 300;		//number of bins to be created
	Double_t minbb = 10;	//lowest bin bound
	Double_t maxbb = 6390;	//highest bin bound
	//Double_t *cust_bb = new Double_t[ncust_bb+1];	//array for bin bounds
	Double_t cust_bb[ncust_bb+1] = {};
	Double_t bmin_log = TMath::Log10(minbb);		//logarithm of lowest bb (--> later used as exponent)
	Double_t bmax_log = TMath::Log10(maxbb);		//logarithm of highest bb (---- "" ----)	
	Double_t dist_log = (bmax_log-bmin_log)/((Double_t)ncust_bb);	//distance of exponents (dist_log=0.3 --> 10^(1.0), 10^(1.3), 10^(1.6),...)
	//cout << "bmin_log:" << bmin_log << endl;
	for(int j=0;j!=ncust_bb+1;++j){
		Double_t cur_expo = bmin_log+(j*dist_log);	//next exponent x of 10
		cust_bb[j] = pow(10, cur_expo);				//bin bound --> 10^x
		cout << "current cust_bb[j] = " << cust_bb[j] << endl;	//test
	}
	//end of "a little more bins"
	*/

	int nbb = sizeof(cust_bb)/sizeof(cust_bb[0])-1;
	sumhisto = new TH1D(Form("sumhisto_%s", bname.c_str()), "; p_{T,avg} /GeV; XS in fb", nbb, cust_bb);	//summary histo


	//in case of monte carlo: prepare histo for XS / Events and one for gen-reco-ptavg distribution (for Unfolding purposes)
	if(montecarlo==true){
		mc_histo = new TH1D(Form("mc_histo_%s", bname.c_str()), "MC Dijet; p_{T,avg} /GeV; Events", nbb, cust_bb);
		mc_histo_smearedJet = new TH1D(Form("mc_histo_smearedJet_%s", bname.c_str()), "MC Dijet (smeared; used Jet_pt_nom); p_{T,avg} /GeV; Events", nbb, cust_bb); //only filled for herwig

		mc_genhisto = new TH1D(Form("mc_genhisto_%s", bname.c_str()), "MC Dijet GenJet spectrum; p^{gen}_{T,avg} /GeV; Events", nbb, cust_bb);
		mc_genrec_noGenCuts = new TH2D(Form("mc_genrec_noGenCuts_%s", bname.c_str()), "MC Dijet (no cuts on gen jets); p^{gen}_{T,avg} /GeV; p^{reco}_{T,avg} /GeV; Events", nbb, cust_bb, nbb, cust_bb);
		mc_gen_reco_ptav = new TH2D(Form("mc_gen_reco_ptavg_%s", bname.c_str()), "MC Dijet p^{reco}_{T,avg} vs. p^{gen}_{T,avg}; p^{gen}_{T,avg} /GeV; p^{reco}_{T,avg} /GeV; Events", nbb, cust_bb, nbb, cust_bb);
		mc_profile =  new TProfile(Form("mc_profile_%s", bname.c_str()), "TProfile of p_{T,avg}; p^{gen}_{T,avg}; #frac{p^{reco}_{T,avg}}{p^{gen}_{T,avg}}", nbb, cust_bb);	//only x bins given?
		mc_profile_j1 =  new TProfile(Form("mc_profile_j1_%s", bname.c_str()), "TProfile of p_{T,j1}; p^{gen}_{T,j1}; #frac{p^{reco}_{T,j1}}{p^{gen}_{T,j1}}", nbb, cust_bb);
		mc_profile_j2 =  new TProfile(Form("mc_profile_j2_%s", bname.c_str()), "TProfile of p_{T,j2}; p^{gen}_{T,j2}; #frac{p^{reco}_{T,j2}}{p^{gen}_{T,j2}}", nbb, cust_bb);
		mc_profile_mergedj1j2 =  new TProfile(Form("mc_profile_mergedj1j2_%s", bname.c_str()), "Merged TProfile of p_{T,j1} and p_{T,j2}; p^{gen}_{T,avg}; #frac{p^{reco}_{T,jets}}{p^{gen}_{T,avg}}", nbb, cust_bb);


		//setting the colors
		mc_profile->SetMarkerColor(kOrange-3);
		mc_profile_j1->SetMarkerColor(kGreen-2);
		mc_profile_j2->SetMarkerColor(kAzure+7);
		mc_profile_mergedj1j2->SetMarkerColor(kRed-6);

	}

	ybys_dirname = bname.c_str();
};


//Overloaded constructor for Yboost Ystar bin to create also unc histograms
YbYs::YbYs (float yboost_lo, float ystar_lo, float yboost_up, float ystar_up, vector<string> uncnames) {
	yb_low = yboost_lo;
	ys_low = ystar_lo;
	yb_up = yboost_up;	
	ys_up = ystar_up;

	bname = Form("yb%01.0fys%01.0f", yb_low, ys_low);
	if(v==true){cout << Form("Created yboost-ystar-bin with: %f <= yb < %f and %f <= ys < %f", yb_low, yb_up, ys_low, ys_up) << endl;}

	//customised pT bin bounds
	//the following are the usual:
	Double_t cust_bb[] = {10, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507, 548, 592, 638, 686, 737, 790, 846, 905, 967, 1032, 1101, 1172, 1248, 1327, 1410, 1497, 1588, 1684, 1784, 1890, 2000, 2116, 2238, 2366, 2500, 2640, 2787, 2941, 3103, 3273, 3450, 3832, 6076, 6389};

	
	int nbb = sizeof(cust_bb)/sizeof(cust_bb[0])-1;
	sumhisto = new TH1D(Form("sumhisto_%s", bname.c_str()), "; XS in fb; p_{T,avg} /GeV", nbb, cust_bb);	//summary histo

	//uncertainty histos (won't be normalised in loop! this has to be done in triggerpaths() function together with the data!!)
	for(int k=0; k!=uncnames.size(); ++k){
		TH1D *cur_unchist_up = new TH1D(Form("unc_%s_Up_%s", uncnames.at(k).c_str(), bname.c_str()), Form("%s Up; p_{T,avg} /GeV; Events", uncnames.at(k).c_str()), nbb, cust_bb);	//unc histo for current unc (k), directory will be set later
		TH1D *cur_unchist_down = new TH1D(Form("unc_%s_Down_%s", uncnames.at(k).c_str(), bname.c_str()), Form("%s Down; p_{T,avg} /GeV; Events", uncnames.at(k).c_str()), nbb, cust_bb);	//unc histo for current unc (k), directory will be set later

		unc_histos_up.insert(make_pair(uncnames.at(k), cur_unchist_up));
		unc_histos_down.insert(make_pair(uncnames.at(k), cur_unchist_down));
	}

	//additional unc histos (jme_folder)
	if((montecarlo==true)&&(strcmp(generator, "herwig")==0)){
		jt_ptraw_hist = new TH1D(Form("ptavg_Jet_pt_raw_%s", bname.c_str()), "Jet_pt_raw;", nbb, cust_bb);
		/*
		jt_ptnom_hist;
		jt_corrJEC_hist;
		jt_corrJER_hist;
		jt_ptjerUp_hist;
		jt_ptjerDw_hist;
		*/
	}
			

	
	//in case of monte carlo: prepare histo for XS / Events and one for gen-reco-ptavg distribution (for Unfolding purposes)
	if(montecarlo==true){
		mc_histo = new TH1D(Form("mc_histo_%s", bname.c_str()), "MC Dijet; p_{T,avg} /GeV; Events", nbb, cust_bb);
		mc_genhisto = new TH1D(Form("mc_genhisto_%s", bname.c_str()), "MC Dijet GenJet spectrum; p^{gen}_{T,avg} /GeV; Events", nbb, cust_bb);
		mc_genrec_noGenCuts = new TH2D(Form("mc_genrec_noGenCuts_%s", bname.c_str()), "MC Dijet (no cuts on gen jets); p^{gen}_{T,avg} /GeV; p^{reco}_{T,avg} /GeV; Events", nbb, cust_bb, nbb, cust_bb);
		mc_gen_reco_ptav = new TH2D(Form("mc_gen_reco_ptavg_%s", bname.c_str()), "MC Dijet p^{reco}_{T,avg} vs. p^{gen}_{T,avg}; p^{gen}_{T,avg} /GeV; p^{reco}_{T,avg} /GeV; Events", nbb, cust_bb, nbb, cust_bb);
		mc_profile =  new TProfile(Form("mc_profile_%s", bname.c_str()), "TProfile of p_{T,avg}; p^{gen}_{T,avg}; #frac{p^{reco}_{T,avg}}{p^{gen}_{T,avg}}", nbb, cust_bb);	//only x bins given?
		mc_profile_j1 =  new TProfile(Form("mc_profile_j1_%s", bname.c_str()), "TProfile of p_{T,j1}; p^{gen}_{T,j1}; #frac{p^{reco}_{T,j1}}{p^{gen}_{T,j1}}", nbb, cust_bb);
		mc_profile_j2 =  new TProfile(Form("mc_profile_j2_%s", bname.c_str()), "TProfile of p_{T,j2}; p^{gen}_{T,j2}; #frac{p^{reco}_{T,j2}}{p^{gen}_{T,j2}}", nbb, cust_bb);
		mc_profile_mergedj1j2 =  new TProfile(Form("mc_profile_mergedj1j2_%s", bname.c_str()), "Merged TProfile of p_{T,j1} and p_{T,j2}; p^{gen}_{T,avg}; #frac{p^{reco}_{T,jets}}{p^{gen}_{T,avg}}", nbb, cust_bb);


		mc_histo_smearedJet = new TH1D(Form("mc_histo_smearedJet_%s", bname.c_str()), "MC Dijet (smeared, used Jet_pt_nom); p_{T,avg} /GeV; Events", nbb, cust_bb);	//test: always write it, empty if not Herwig


		//setting the colors
		mc_profile->SetMarkerColor(kOrange-3);
		mc_profile_j1->SetMarkerColor(kGreen-2);
		mc_profile_j2->SetMarkerColor(kAzure+7);
		mc_profile_mergedj1j2->SetMarkerColor(kRed-6);

	}

	ybys_dirname = bname.c_str();
};


//Simple constructor for YbYs, only sets directories
//YbYs::YbYs(){};

//Constructor for Histos (should rather be done once and as member of Etabins?)
//TH1::AddDirectory(false); //just testing this
Histos::Histos (vector<int> triggers) { 		//should the constructor delete all existing histograms that will be recreated here? How to avoid potential memory leak? 
												//(histograms different, but histogram names identical for different eta-objects/eta-bins)
	//customised pT bin bounds
	//the following are the usual:
	//Double_t cust_bb[] = {10, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507, 548, 592, 638, 686, 737, 790, 846, 905, 967, 1032, 1101, 1172, 1248, 1327, 1410, 1497, 1588, 1684, 1784, 1890, 2000, 2116, 2238, 2366, 2500, 2640, 2787, 2941, 3103, 3273, 3450, 3832, 6076, 6389};

	///*
	//this are a little more bins in the central GeV region --> for triggerpaths()
	const int ncust_bb = 300;		//number of bins to be created
	Double_t minbb = 10;	//lowest bin bound
	Double_t maxbb = 6390;	//highest bin bound
	//Double_t *cust_bb = new Double_t[ncust_bb+1];	//array for bin bounds
	Double_t cust_bb[ncust_bb+1] = {};
	Double_t bmin_log = TMath::Log10(minbb);		//logarithm of lowest bb (--> later used as exponent)
	Double_t bmax_log = TMath::Log10(maxbb);		//logarithm of highest bb (---- "" ----)	
	Double_t dist_log = (bmax_log-bmin_log)/((Double_t)ncust_bb);	//distance of exponents (dist_log=0.3 --> 10^(1.0), 10^(1.3), 10^(1.6),...)
	//cout << "bmin_log:" << bmin_log;
	for(int j=0;j!=ncust_bb+1;++j){
		Double_t cur_expo = bmin_log+(j*dist_log);	//next exponent x of 10
		cust_bb[j] = pow(10, cur_expo);				//bin bound --> 10^x
		//cout << "current cust_bb[j] = " << cust_bb[j] << endl;	//test
	}
	//end of "a little more bins"
	//*/

	////cout << "sizeof(cust_bb): " << sizeof(cust_bb) << endl;
	////cout << "sizeof(cust_bb[0]): " << sizeof(cust_bb[0]) << endl;
	int nbb = sizeof(cust_bb)/sizeof(cust_bb[0])-1;
	
	/*
	for(int ind=0; ind!=nbb+1; ++ind){
		cout << "ind: " << ind << " binbound: " << cust_bb[ind] << endl;	
	}
	*/
	

	//hjt_names.clear();		//empty the perhaps prefilled (via previous calls) vector (but maintain the pointers --> how?)
	for(int j=0; j!=triggers.size(); ++j){ //or start from 1, because first one is "0" --> without trigger (here not taken into account... only matters in SetBranchStatus)
		//if(j==0){cout << "triggers.at(0)= " << triggers.at(j) << endl;}
		trigger_names.push_back(Form("HLT_PFJet%i", triggers.at(j)));
		hjt_names.push_back(Form("hpt_%i", triggers.at(j)));
		//cout << "hjt_names.back(): " << hjt_names.back() << endl;
		//histvec.at(j) = new TH1D(hjt_names.at(j), "jet_pt;p_{T};Events", 100, 10, 6000);
		//TH1D *curhist =  new TH1D(hjt_names.back(), "jet_pt;p_{T};Events", 100, 10, 6000);

		TH1D *curhist = new TH1D(hjt_names.back(), "jet_pt; p_{T,avg};Events", nbb, cust_bb);		//x-axis title etc. have to be adjusted to dijet --> do this in plotting script!
		TH1D *curhist_trgobj = new TH1D(Form("%s_trgobj", hjt_names.back()), "HLT Trigger N-1 && TrigObj_pt N;p_{T,avg};Events", nbb, cust_bb);
		TH1D *curhist_n2trgobj = new TH1D(Form("%s_n2trgobj", hjt_names.back()), "HLT Trigger N-2 && TrigObj_pt N;p_{T,avg};Events", nbb, cust_bb);

		//test with very small bins
		//TH1D *curhist = new TH1D(hjt_names.back(), "jet_pt;p_{T};Events", 63900, 10, 6400);	//book histogram with very fine binning

		jt_histos.insert(make_pair(triggers.at(j), curhist));
		trgobj_histos.insert(make_pair(triggers.at(j), curhist_trgobj));
		n2_trgobj_histos.insert(make_pair(triggers.at(j), curhist_n2trgobj));
		//gDirectory->ls(); //testing
	} 
	if(v==true){for(int j=0; j!=triggers.size(); ++j){cout << "hjt_names.at(j) = " << hjt_names.at(j) << " for j==" << j << endl;
											cout << "hjt_names.at(j) = " << addressof(hjt_names.at(j)) << " for j==" << j << endl;}}

};

//overloaded histos constructor, used in the triggerstudies function
//histograms in the jt_histos map etc. are filled according to given file
Histos::Histos (vector<int> triggers, TFile *histo_file, const char* ybysname) { 	

	for(int j=0; j!=triggers.size(); ++j){ //or start from 1, because first one is "0" --> without trigger (here not taken into account... only matters in SetBranchStatus)
		//if(j==0){cout << "triggers.at(0)= " << triggers.at(j) << endl;}
		trigger_names.push_back(Form("HLT_PFJet%i", triggers.at(j)));
		hjt_names.push_back(Form("hpt_%i", triggers.at(j)));

		//clone the already filled histograms from the input file histo_file
		TH1D *curhist;
		TH1D *curhist_trgobj;
		TH1D *curhist_n2trgobj;
		histo_file->GetObject(Form("Standard/%s/jt%d/%s", ybysname, triggers.at(j), hjt_names.back()), curhist);
		histo_file->GetObject(Form("Standard/%s/jt%d/%s_trgobj", ybysname, triggers.at(j), hjt_names.back()), curhist_trgobj);
		histo_file->GetObject(Form("Standard/%s/jt%d/%s_n2trgobj", ybysname, triggers.at(j), hjt_names.back()), curhist_n2trgobj);

		jt_histos.insert(make_pair(triggers.at(j), curhist));
		trgobj_histos.insert(make_pair(triggers.at(j), curhist_trgobj));
		n2_trgobj_histos.insert(make_pair(triggers.at(j), curhist_n2trgobj));
		//gDirectory->ls(); //testing
	} 
};

//default Histos constructor, used in DijetAnalysis::FastLoop() with set triggers and set turn-on points
Histos::Histos(){
	vector<int> triggers = {0, 15, 25, 40, 60, 80, 140, 200, 260, 320, 400, 450, 500, 550};		//for each of these there will be a folder in each ybys bin

	//customised pT bin bounds
	Double_t cust_bb[] = {10, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507, 548, 592, 638, 686, 737, 790, 846, 905, 967, 1032, 1101, 1172, 1248, 1327, 1410, 1497, 1588, 1684, 1784, 1890, 2000, 2116, 2238, 2366, 2500, 2640, 2787, 2941, 3103, 3273, 3450, 3832, 6076, 6389};
	int nbb = sizeof(cust_bb)/sizeof(cust_bb[0])-1;

	//hjt_names.clear();		//empty the perhaps prefilled (via previous calls) vector (but maintain the pointers --> how?)
	for(int j=0; j!=triggers.size(); ++j){ //or start from 1, because first one is "0" --> without trigger (here not taken into account... only matters in SetBranchStatus)
		//if(j==0){cout << "triggers.at(0)= " << triggers.at(j) << endl;}
		trigger_names.push_back(Form("HLT_PFJet%i", triggers.at(j)));
		hjt_names.push_back(Form("hpt_%i", triggers.at(j)));

		TH1D *curhist = new TH1D(hjt_names.back(), "jet_pt; p_{T,avg};Events", nbb, cust_bb);		//x-axis title etc. have to be adjusted to dijet --> do this in plotting script!
		jt_histos.insert(make_pair(triggers.at(j), curhist));
	} 
};




//--------------------------------------------------------------------------------------------------//
//											 MEMBER FUNCTIONS										//
//--------------------------------------------------------------------------------------------------//

//definition of the member functions of the class YbYs (used for the ybys binning)
//these functions don't need to get passed any arguments as they use the data members of their objects

//directory function --> for handling the directory structure and saving it to the "intermediate result output file"
//Histos curtrig are the "current(=belonging to this ybysbin) trigger histograms"
//trgmode --> if true: write the histograms (as they are copied from input file), otherwise just set the directory
void YbYs::ybys_directories(TDirectory *standard_folder, vector<int> triggers, Histos curtrig, bool trgmode=false){	//should also take the histograms (one per trigger) for this ybys-bin
	//create and handle directory structure of output file
	const char *dirname;  
	dirname = Form("%s", bname.c_str());
	standard_folder->mkdir(dirname);
	TDirectory *Dirname = standard_folder->GetDirectory(dirname);
	Dirname->cd();
	//assign sumhisto to this ybys dir
	sumhisto->SetDirectory(Dirname);

	//add trigger-subdirectories
	//curtrig = trigger histograms which belong to current ybys-bin
	for (int j=0; j!=triggers.size(); ++j){
		const char *tdir; //name of trigger directory
		tdir = Form("jt%i", triggers.at(j));
		Dirname->mkdir(tdir);
		TDirectory *tdirname = Dirname->GetDirectory(tdir);
		tdirname->cd();
		curtrig.jt_histos[triggers.at(j)]->SetDirectory(tdirname); 			//move corresponding histogram to correct trigger-directory in current eta-directory
		curtrig.trgobj_histos[triggers.at(j)]->SetDirectory(tdirname);		//move histo to correct trigger dir in current ybys dir
		curtrig.n2_trgobj_histos[triggers.at(j)]->SetDirectory(tdirname);	// ---- "" ----
		if(trgmode){
			curtrig.jt_histos[triggers.at(j)]->Write();
			curtrig.trgobj_histos[triggers.at(j)]->Write();
			curtrig.n2_trgobj_histos[triggers.at(j)]->Write();
		}
	}
	
	//change back to standard_folder
	standard_folder->cd();
};

//overloaded function --> in case of MC, do not need trigger-folders:
//directly create TH1D here in correct folder
///void YbYs::ybys_directories(TDirectory *standard_folder){
TDirectory* YbYs::ybys_directories(TDirectory *standard_folder){
	//create and handle directory structure of output file
	const char *dirname;  
	dirname = Form("%s", bname.c_str());
	standard_folder->mkdir(dirname);
	TDirectory *ybysDir = standard_folder->GetDirectory(dirname);
	ybysDir->cd();
	
	//created histogram for this ybys bin in constructor, now put it to correct dir
	mc_histo->SetDirectory(ybysDir);
	mc_histo_smearedJet->SetDirectory(ybysDir);
	mc_genhisto->SetDirectory(ybysDir);
	mc_genrec_noGenCuts->SetDirectory(ybysDir);
	mc_gen_reco_ptav->SetDirectory(ybysDir);
	mc_profile->SetDirectory(ybysDir);
	mc_profile_j1->SetDirectory(ybysDir);
	mc_profile_j2->SetDirectory(ybysDir);
	mc_profile_mergedj1j2->SetDirectory(ybysDir);
	standard_folder->cd();

	//as also in data analysis --> set "extra directory" for additional plots (Zmass etc.)
	//const char *extradir_name = Form("additional_plots_%s", bname.c_str());
	//ybysDir->mkdir(extradir_name);
	//TDirectory *extraDir = ybysDir->GetDirectory(extradir_name);
	//extraDir->cd();

	//change back to standard_folder:
	standard_folder->cd();
	
	return ybysDir;
};

//implementation for FastLoop() function
void YbYs::ybys_fastdirs(TDirectory *standard_folder, vector<int> triggers, Histos curtrig){	//should also take the histograms (one per trigger) for this ybys-bin
	//create and handle directory structure of output file
	const char *dirname;  
	dirname = Form("%s", bname.c_str());
	standard_folder->mkdir(dirname);
	TDirectory *Dirname = standard_folder->GetDirectory(dirname);
	Dirname->cd();
	//assign sumhisto to this ybys dir
	sumhisto->SetDirectory(Dirname);

	//add trigger-subdirectories
	//curtrig = trigger histograms which belong to current ybys-bin
	for (int j=0; j!=triggers.size(); ++j){
		const char *tdir; //name of trigger directory
		tdir = Form("jt%i", triggers.at(j));
		Dirname->mkdir(tdir);
		TDirectory *tdirname = Dirname->GetDirectory(tdir);
		tdirname->cd();
		curtrig.jt_histos[triggers.at(j)]->SetDirectory(tdirname); 			//move corresponding histogram to correct trigger-directory in current eta-directory
	}
	
	//change back to standard_folder
	standard_folder->cd();
};



//this function ONLY creates the ybys directory and returns it
//no directories set for any histograms
TDirectory* YbYs::ybys_mkdir(TDirectory *standard_folder){
	//create and handle directory structure of output file
	const char *dirname;  
	dirname = Form("%s", bname.c_str());
	standard_folder->mkdir(dirname);
	TDirectory *ybysDir = standard_folder->GetDirectory(dirname);
	ybysDir->cd();
	standard_folder->cd();
	return ybysDir;
};


//for unc histos
void YbYs::ybys_uncdirs(TDirectory *unc_folder, vector<string> uncnames){
	//create the unc histos in the uncertainty_folder
	const char *dirname;
	dirname = Form("%s", bname.c_str());
	unc_folder->mkdir(dirname);
	TDirectory *ybysDir = unc_folder->GetDirectory(dirname);
	ybysDir->cd();

	//put the unc histos to this directory --> should these be separated into two folders Uncs_Up and Uncs_Down? (so far in the same)
	for(int k=0; k!=uncnames.size(); ++k){
		unc_histos_up[uncnames.at(k)]->SetDirectory(ybysDir);
		unc_histos_down[uncnames.at(k)]->SetDirectory(ybysDir);
	}
};


void YbYs::ybys_jmedirs(TDirectory *jme_folder){
	const char *dirname;
	dirname = Form("%s", bname.c_str());
	jme_folder->mkdir(dirname);
	TDirectory *ybysDir_jme = jme_folder->GetDirectory(dirname);
	ybysDir_jme->cd();

	//put histograms created in constructor to their correct directory
	jt_ptraw_hist->SetDirectory(ybysDir_jme);
	jt_ptnom_hist->SetDirectory(ybysDir_jme);
	jt_corrJEC_hist->SetDirectory(ybysDir_jme);
	jt_corrJER_hist->SetDirectory(ybysDir_jme);
	jt_ptjerUp_hist->SetDirectory(ybysDir_jme);
	jt_ptjerDw_hist->SetDirectory(ybysDir_jme);

	jme_folder->cd();
};








//currently not used? --> could make things easier... check it out.
void YbYs::trglumi_directories(TDirectory *standard_folder, vector<int> triggers, map<const char*, TH1D*> trglumi_hist_map){
	//create and handle directory structure of output file
	const char *dirname;  
	dirname = Form("%s", bname.c_str());
	standard_folder->mkdir(dirname);
	TDirectory *Dirname = standard_folder->GetDirectory(dirname);
	Dirname->cd();
	//bdir = Dirname;

	for(int j=1; j!=triggers.size(); ++j){ //or start from 1, because first one is "0" --> without trigger
		int curtrig = triggers.at(j); 		//current trigger
		int prevtrig = triggers.at(j-1); 	//previous trigger
		TH1D *curhist =  new TH1D(Form("jet%i_jet%i", curtrig, prevtrig), Form("trigger comparison;p_{T,avg};#frac{jet%i / lumi%i}{jet%i / lumi%i}", curtrig, curtrig, prevtrig, prevtrig), 62900, 10, 6300);
		curhist->SetDirectory(Dirname);
		trglumi_hist_map.insert(make_pair(Form("jet%i_jet%i", curtrig, prevtrig), curhist));
	} 
};






//###### ###### ###### ######  UNDER CONSTRUCTION ###### ###### ###### ######//
// ------------ this is a draft ---------------- CHECK AGAIN! --------------//
//function for filling the summary histogram of this ybys bin
//make a map instead of a vector? --> effmap[trig], mappe eff99Xval auf trigger-GeV
///void YbYs::fill_sumhisto(vector<Double_t> effvec, Histos curtrig){
//uses results of trigger-fitting method 1  (combi fits of N/(N-2) and N/(N-1))
TH1D* YbYs::fill_sumhisto(map<int, Double_t> effmap, vector<int> triggers, Histos *curtrig, vector<Double_t> lumivec){	//really want to work with pointers? rather Histos?
	//cout << "In fill_sumhisto() function." << endl;

	//get old histograms, members of "curtrig" (one per trigger)
	///int nbins = sumhisto->GetNbinsX();
	int nbins_sumhist = sumhisto->GetNbinsX();								//number of bins in summary histogram
	int nbins_trghist = curtrig->jt_histos[triggers.at(0)]->GetNbinsX();	//number of bins in trigger histos (here retrieved from 1st one, 0GeV)

	int donebins = 0;	//bins of the old histogram that have already been taken care of
	for(int bin_ind=nbins_sumhist; bin_ind!=0; --bin_ind){	//loop through bins of new histogram (sumhisto)
	///for(int bin_ind=nbins; bin_ind!=0; --bin_ind){
		///cout << "in for-loop through bins." << endl;
		///double curpt = sumhisto->GetXaxis()->GetBinCenter(bin_ind);
		double curpt = sumhisto->GetXaxis()->GetBinLowEdge(bin_ind);

		//now loop through effmap, if this trigger doesn't fit, go through next lowest.
		//zähle trgnr von oben runter
		//for(int trgnr=triggers.size()-1; trgnr!=4; --trgnr){	//for now just go down to trg 80 GeV. (ind=5)
		for(int trgnr=triggers.size()-1; trgnr!=2; --trgnr){	//go down to trg 40 GeV. (ind=3) --> newly implemented (extrapolation)
			///cout << "trgnr= " << trgnr << " and triggers.at(trgnr): " << triggers.at(trgnr) << " GeV." << endl;
			if(bname=="yb0ys0"){cout << Form("curpt= %f and effmap[triggers.at(trgnr)]= %f", curpt, effmap[triggers.at(trgnr)]) << endl;}
			if(curpt>effmap[triggers.at(trgnr)]){
				//cout << "Bincont: " << curtrig->jt_histos[triggers.at(trgnr)]->GetBinContent(bin_ind) << endl;
				//sumhisto->SetBinContent(bin_ind, curtrig->jt_histos[triggers.at(trgnr)]->GetBinContent(bin_ind));

				//fill all the old histogram bins above this threshold into the corresponding new bigger sumhisto bin
				//loop through the bins of the old histogram, break before loop would run into "next lower new bin"
				double new_bincont = 0;		//will contain all the old bin contents that contribute to the new (wider) bin
				double new_binerr_quad = 0;	//will contain the quadratic sum before taking the sqrt of the statistical errors of all old bins contributing to the new bin
				double new_binerr = 0;		//will contain the new statistic error for the current bin
				for(int oldbin_ind=(nbins_trghist-donebins); oldbin_ind!=0; --oldbin_ind){
					double old_ptlow = curtrig->jt_histos[triggers.at(trgnr)]->GetBinLowEdge(oldbin_ind);
					if(old_ptlow >= curpt){	//calculating new bin content
						new_bincont += (curtrig->jt_histos[triggers.at(trgnr)]->GetBinContent(oldbin_ind));
						//sum up the bin errors in quadrature (for each of the new bins)
						new_binerr_quad += pow((curtrig->jt_histos[triggers.at(trgnr)]->GetBinError(oldbin_ind)), 2);
						donebins +=1;	//one more old bin used. count.
					}
					else{break;}	//break, because this new bin is done now.
				}
				///sumhisto->SetBinContent(bin_ind, (curtrig->jt_histos[triggers.at(trgnr)]->GetBinContent(bin_ind))/lumivec.at(trgnr-1));
				sumhisto->SetBinContent(bin_ind, new_bincont/lumivec.at(trgnr-1));	//set new bin content, normalised by luminosity

				//take the square root of the individual bin errors added in quadrature to get the error for the current new bin
				new_binerr = TMath::Sqrt(new_binerr_quad);
				//divide sqrt(sum_squared) by corresponding luminosity in the current bin and set the bin error for the current bin
				sumhisto->SetBinError(bin_ind, new_binerr/lumivec.at(trgnr-1));	//set new bin content, normalised by luminosity

				break; //have found correct trigger.
			}
			else{
				//if(	//if we have not checked all of them yet
				continue;	//try range of next trigger
				//else if(){} //if curpt lies lower than lowest trigger eff threshold --> fill it below?? keep it unchanged??
			}
		}
	}
	cout << "sumhisto->GetEntries: " << sumhisto->GetEntries() << endl;
	sumhisto->SetLineWidth(3);
	sumhisto->SetLineColor(kAzure);
	sumhisto->SetMarkerColor(kAzure);
	sumhisto->Draw();
	//gROOT->cd(ybys_dirname);
	sumhisto->Write();
	return sumhisto;	//otherwise problem with writing within loop... write the returned histo.
	//getting center (or lower bin bound?) of current bin
	//go through all the bins (starting from highest)
	//if
		//sumhisto->SetBinContent(curtrig.jt_histos[triggers.at(j)]->GetBinContent());
		//sumhisto->Fill();
	//curtrig.jt_histos[triggers.at(j)]->Fill(cur_ptavg);	//not ready yet... this is testing!
};



//fill histos in DijetAnalysis::FastLoop()
TH1D* YbYs::fill_fastsumhisto(vector<int> triggers, Histos *curtrig, vector<Double_t> lumivec){	//really want to work with pointers? rather Histos?

	//fixed effmap
	//map containing the int-value of the trigger as the key and the 99% eff. turn-on point of the trigger as the value
	//the turn-on points have been obtained using the emulation method (N1 and N2 combi) for Run D 2018
	map<int, Double_t> effmap = {{40, 51.6}, {60, 77.0}, {80, 104.4}, {140, 176.1}, {200, 242.5}, {260, 310.8}, {320, 383.8}, {400, 467.8}, {450, 521.5}, {500, 568.7}, {550, 621.6}};

	//get old histograms, members of "curtrig" (one per trigger)
	///int nbins = sumhisto->GetNbinsX();
	int nbins_sumhist = sumhisto->GetNbinsX();								//number of bins in summary histogram
	int nbins_trghist = curtrig->jt_histos[triggers.at(0)]->GetNbinsX();	//number of bins in trigger histos (here retrieved from 1st one, 0GeV)

	//check if they are the same... should be!
	if(nbins_sumhist==nbins_trghist){
		cout << "Summary histogram and trigger histograms have same amount of entries! Go on! :)" << endl;
	}
	else{
		cout << Form("Entries of sumhisto (%d) and trghistos (%d) differ! Aborting!", nbins_sumhist, nbins_trghist) << endl;
		exit(EXIT_FAILURE);
	}
	
	for(int bin_ind=nbins_sumhist; bin_ind!=0; --bin_ind){	//loop through bins of new histogram (sumhisto) == at the same time through old histos
		///cout << "in for-loop through bins." << endl;
		double curpt = sumhisto->GetXaxis()->GetBinLowEdge(bin_ind);

		//now loop through effmap, if this trigger doesn't fit, go through next lowest.
		//zähle trgnr von oben runter
		//for(int trgnr=triggers.size()-1; trgnr!=4; --trgnr){	//for now just go down to trg 80 GeV. (ind=5)
		for(int trgnr=triggers.size()-1; trgnr!=2; --trgnr){	//go down to trg 40 GeV. (ind=3) --> newly implemented (extrapolation)
			///cout << "trgnr= " << trgnr << " and triggers.at(trgnr): " << triggers.at(trgnr) << " GeV." << endl;
			if(bname=="yb0ys0"){cout << Form("curpt= %f and effmap[triggers.at(trgnr)]= %f", curpt, effmap[triggers.at(trgnr)]) << endl;}
			if(curpt>effmap[triggers.at(trgnr)]){//check if current low edge of the summary histo is still above the current trigger's 99% eff point
				double new_bincont = curtrig->jt_histos[triggers.at(trgnr)]->GetBinContent(bin_ind);	//get the bin content of the trigger histo for this bin
				sumhisto->SetBinContent(bin_ind, new_bincont/lumivec.at(trgnr-1));						//set new bin content, normalised by luminosity
				double new_binerr = curtrig->jt_histos[triggers.at(trgnr)]->GetBinError(bin_ind);		//get the bin error of the trghisto for this bin
				sumhisto->SetBinError(bin_ind, new_binerr/lumivec.at(trgnr-1));							//set new bin error (=old bin error in trghiso), normalised by luminosity
				break; //have found correct trigger. Move to next lower bin (to be filled) in sumhisto.
			}
			else{//if the current trigger is not efficient yet at current pt
				//if(	//if we have not checked all of them yet
				continue;	//try range of next trigger
				//else if(){} //if curpt lies lower than lowest trigger eff threshold --> fill it below?? keep it unchanged??
			}
		}
	}
	cout << "sumhisto->GetEntries: " << sumhisto->GetEntries() << endl;
	sumhisto->SetLineWidth(3);
	sumhisto->SetLineColor(kAzure);
	sumhisto->SetMarkerColor(kAzure);
	sumhisto->Draw();
	//gROOT->cd(ybys_dirname);
	sumhisto->Write();
	return sumhisto;	//otherwise problem with writing within loop... write the returned histo.
};



//--------------------------------------------------------------------------------------------------//
//											 OTHER FUNCTIONS										//
//--------------------------------------------------------------------------------------------------//

//overloaded drawstyle function for either TH1D or TGraphErrors
//function to customize how histo/graph is drawn in rootfile
void SetDrawStyleWrite(TH1D *object, const char *drawcmd="", const char *title="", const char *xtitle="p_{T,avg}", const char *ytitle="", Color_t linecolor=kGray+2, Int_t linewidth=2, Style_t markerstyle=kFullCircle, Color_t markercolor=kGray+2, Int_t markersize=1, bool morelog=false, Double_t xTitleOffset=1.2, bool writing=false){
	object->Draw(drawcmd);
	object->SetTitle(title);
	object->GetXaxis()->SetTitle(xtitle);
	object->GetYaxis()->SetTitle(ytitle);
	object->GetXaxis()->SetTitleSize(0.05);		//current default for all the histos
	object->GetYaxis()->SetTitleSize(0.05);		//current default for all the histos
	object->GetXaxis()->SetLabelSize(0.04);		//current default for all the histos
	object->GetYaxis()->SetLabelSize(0.04);		//current default for all the histos

	object->GetXaxis()->SetTitleOffset(xTitleOffset);

	object->SetLineColor(linecolor);
	object->SetLineWidth(linewidth);

	object->SetMarkerStyle(markerstyle);
	object->SetMarkerColor(markercolor);
	object->SetMarkerSize(markersize);

	if(morelog){
		object->GetXaxis()->SetMoreLogLabels();
		object->GetXaxis()->SetNoExponent();
	}

	if(writing){
		object->Write();
	}
}

void SetDrawStyleWrite(TGraphErrors *object, const char *drawcmd="", const char *title="", const char *xtitle="p_{T,avg}", const char *ytitle="", Color_t linecolor=kGray+2, Int_t linewidth=2, Style_t markerstyle=kFullCircle, Color_t markercolor=kGray+2, Int_t markersize=2, bool morelog=false, Double_t xTitleOffset=1.2){
	object->Draw(drawcmd);
	object->SetTitle(title);
	object->GetXaxis()->SetTitle(xtitle);
	object->GetYaxis()->SetTitle(ytitle);
	object->GetXaxis()->SetTitleSize(0.05);		//current default for all the histos
	object->GetYaxis()->SetTitleSize(0.05);		//current default for all the histos
	object->GetXaxis()->SetTitleOffset(xTitleOffset);

	object->SetLineColor(linecolor);
	object->SetLineWidth(linewidth);

	object->SetMarkerStyle(markerstyle);
	object->SetMarkerColor(markercolor);
	object->SetMarkerSize(markersize);

	if(morelog){
		object->GetXaxis()->SetMoreLogLabels();
		object->GetXaxis()->SetNoExponent();
	}

	object->Write();
}

//function to initialize fitting function
//make sure (outside function), that fctname == name of the fit in code
TF1* InitFit(const char* fctname, const char* function, Double_t range_low, Double_t range_up, Color_t linecolor, vector<Double_t> &init_params, vector<const char*> &par_names){
	TF1 *fit = new TF1(fctname, function, range_low, range_up);
	fit->SetLineColor(linecolor);
	fit->SetLineWidth(2);

	//set initial parameters:
	for(int i=0; i!=init_params.size(); ++i){
		fit->SetParameter(i, init_params.at(i));
		fit->SetParName(i, par_names.at(i));
	}
	
	return fit;
}

//function to initialize 2 parameter fitting function
TF1* InitFit2P(const char* fctname, const char* function, Double_t range_low, Double_t range_up, Color_t linecolor, Double_t init_p0, Double_t init_p1, const char* par0_name="p0", const char* par1_name="p1"){
	TF1 *fit = new TF1(fctname, function, range_low, range_up);
	fit->SetLineColor(linecolor);
	fit->SetLineWidth(2);

	//set initial parameters:
	fit->SetParameters(init_p0, init_p1);
	fit->SetParNames(par0_name, par1_name);

	/*
	for(int i=0; i!=init_params.size(); ++i){
		fit->SetParameter(i, init_params.at(i));
		fit->SetParName(i, par_names.at(i));
	}
	*/
	
	return fit;
}



//function to do fit and store resulting parameters
void DoFit(TH1D *data, TF1 *fitfct, Double_t xlow, Double_t xup, vector<Double_t> &pnomvec, vector<Double_t> &pwidvec){
	data->Fit(fitfct);
}

void DoFit(TGraphErrors *data){}	//implementation for multigraph-fitting



//function to store fit-parameters and their errors for later access
//always call this twice per fit --> once for pnom, once for pwid
void StoreParam(Int_t nom, Double_t cur_p, Double_t cur_err, vector<Double_t> &par_x, vector<Double_t> &par_y, vector<Double_t> &par_err){
	par_x.push_back(nom);
	par_y.push_back(cur_p/nom);
	par_err.push_back(cur_err/nom);
}

//can be used for 2-parameter fit (pnom, pwid) directly
void Store2Params(Int_t nom, TF1* fit, vector<Double_t> &pnom_x, vector<Double_t> &pnom_y, vector<Double_t> &pnom_err, vector<Double_t> &pwid_x, vector<Double_t> &pwid_y, vector<Double_t> &pwid_err){
	//get the parameters from the fit:
	//could be held more general with [0] and [1]
	Int_t pnom_ind = fit->GetParNumber("pnom");
	Int_t pwid_ind = fit->GetParNumber("pwid");
	///Double_t plateau = fit->GetParameter("plateau");
	Double_t p_nom = fit->GetParameter("pnom");
	Double_t p_wid = fit->GetParameter("pwid");

	//append to x-axis and y-axis vectors for the parameter-stability-check
	pnom_x.push_back(nom);			//nominal value of trigger
	pnom_y.push_back(p_nom/nom); 	//pnom/nom
	pnom_err.push_back((fit->GetParError(pnom_ind))/nom);

	pwid_x.push_back(nom); 			//could also plot against sqrt(nom), but then use other y-values
	pwid_y.push_back(p_wid/nom);	//pwid/nom
	pwid_err.push_back((fit->GetParError(pwid_ind))/nom);
}


//function to create and cut multigraph
//only returns the cut multigraph.
//TMultigraph MakeMultigraphCut(TH1D* N1hist, TH1D* N2hist, Double_t lowxN1, Double_t upxN2){
TMultiGraph* MakeMultigraphCut(TH1D* N1hist, TH1D* N2hist, TF1* N1fit, TF1* N2fit, int curtrig){
	//happens still within trigger(fitting-) loop
	//need int curtrig only for setting the correct title of the mg
	//create a multigraph to use N/N-2 for the lower part of the turn-on curve and N/N-1 result for the part near the plateau
	//mg contains both histograms, "converted" to TGraphErrors graphs
	TMultiGraph *mg = new TMultiGraph();

	int nbinsN1 = N1hist->GetNbinsX(); // == nbinsN2 --> number of bins is set at initial creation and equals for these histos
	/////int nbinsN1 = N1hist->GetEntries();
	//--> problem: this includes bins with 0 bin content. Remove those. (-->how and when?)
	Double_t n1hist_x[nbinsN1];
	Double_t n1hist_y[nbinsN1];
	Double_t n1hist_err[nbinsN1];
	Double_t n2hist_x[nbinsN1];
	Double_t n2hist_y[nbinsN1];
	Double_t n2hist_err[nbinsN1];

	for(int bin=1; bin!=nbinsN1+1; ++bin){
		n1hist_x[bin]=N1hist->GetXaxis()->GetBinCenter(bin);
		n1hist_y[bin]=N1hist->GetBinContent(bin);
		n1hist_err[bin]=N1hist->GetBinError(bin);

		n2hist_x[bin]=N2hist->GetXaxis()->GetBinCenter(bin); //--> should actually be the same as n2hist_x values
		n2hist_y[bin]=N2hist->GetBinContent(bin);
		n2hist_err[bin]=N2hist->GetBinError(bin);
	}
	TGraphErrors *N1hist_graph = new TGraphErrors(nbinsN1, &n1hist_x[0], &n1hist_y[0], 0, &n1hist_err[0]);
	TGraphErrors *N2hist_graph = new TGraphErrors(nbinsN1, &n2hist_x[0], &n2hist_y[0], 0, &n2hist_err[0]); //same amount of bins as n1hist

	N1hist_graph->SetLineColor(kGreen+2);
	N1hist_graph->SetMarkerStyle(kOpenTriangleUp);
	N1hist_graph->SetMarkerColor(kGreen+2);
	N1hist_graph->SetMarkerSize(1.2);
	N2hist_graph->SetLineColor(kBlue+2);
	N2hist_graph->SetMarkerStyle(kOpenCircle);
	N2hist_graph->SetMarkerColor(kBlue+2);
	N2hist_graph->SetMarkerSize(1.2);

	mg->Add(N1hist_graph);	//Add N/N-1 histo
	mg->Add(N2hist_graph);	//Add N/N-2 histo

	mg->SetTitle(Form("N/N-1 and N/N-2 trigger-ratios; p_{T,avg}; Ratio of trg%d to (N-1) and (N-2)", curtrig)); //triggers.at(ind) instead of curtrig
	mg->Draw("APZ");

	//create new multigraph, containing only the points wanted for the fit
	//remove all the points in N1hist_graph, where x-value is lower than N1fit->GetX(0.4) 
	//remove all the points in N2hist_graph, where x-value is higher than N2fit->GetX(0.8)
	Double_t lowxN1 = N1fit->GetX(0.4);
	Double_t upxN2	= N2fit->GetX(0.8);

	if(v==true){ // testing --> should be checked (not only if v==true)... do the correct points get removed? are too many lost? --> adjusting cuts?
		cout << "remove all N1-points lower than x: lowxN1 ===>>> " << lowxN1 << endl;
		cout << "remove all N2-points higher than x: upxN2 ===>>> " << upxN2 << endl;
	}
	//--> this procedure leads to problems when previous fit already had problems!!! --> needs to be improved and changed!

//just testing for now
//if(ybysnr==0 && ind==13){}
	for(int bin_i=0; bin_i!=N1hist_graph->GetN(); ++bin_i){
		//cout << "X: " << N1hist_graph->GetX()[bin_i] << endl;
	}
	//cout << "Now remove points." << endl;
	for(int bin_i=0; bin_i!=(N1hist_graph->GetN()); ++bin_i){
	//for(int bin_i=0; bin_i!=nbinsN1;){ //try without set increment
		//--bin_i;
		//cout << "point index = " << bin_i << "|| X-value: " << N1hist_graph->GetX()[bin_i] << endl;
		//if((N1hist_graph->GetX()[bin_i])<lowxN1){
		while((N1hist_graph->GetX()[bin_i])<lowxN1){	//--> THIS DOES NOT WORK WHEN EARLY FLUCTUATIONS APPEAR	 --> better use some "if" condition.
			//cout << Form("%f is lower than %f", N1hist_graph->GetX()[bin_i], lowxN1) << endl;
			//cout << Form("X-value of bin_i=%d before removal of point: X= %f", bin_i, N1hist_graph->GetX()[bin_i]) << endl;
			N1hist_graph->RemovePoint(bin_i);
			//cout << Form("X-value of bin_i=%d AFTER removal of point: X= %f", bin_i, N1hist_graph->GetX()[bin_i]) << endl;

			//need to set current index lower / check this index again, 
			//because all the points above bin_i moved to a new index (ind_new = ind_old-1)
			//--bin_i;
		}
		//else{bin_i++;}	//only increment if points are not "moving downwards" indexwise anymore
	}
	//due to new indices, and the fact that N1hist_graph contains now less points than before -->
	//need an extra loop for N2hist_graph!
	//cout << "N2hist_graph->GetN() --> " << N2hist_graph->GetN() << endl;
	for(int bin_i=0; bin_i<(N2hist_graph->GetN()); ++bin_i){ //otherwise end up in infinite loop... (so use "<" instead of "!=")
		//cout << Form("In bin_i= %d, X-value of N2hist_graph: %f // reminder of upxN2: %f", bin_i, (N2hist_graph->GetX()[bin_i]), upxN2) << endl;

		if((N2hist_graph->GetX()[bin_i])>upxN2){
		//while((N2hist_graph->GetX()[bin_i])>upxN2){
			//cout << Form("X-value of bin_i=%d before removal of point: X= %f", bin_i, N2hist_graph->GetX()[bin_i]) << endl;

			N2hist_graph->RemovePoint(bin_i);
			//cout << "cut from above. ##### " << endl;

			//cout << Form("X-value of bin_i=%d AFTER removal of point: X= %f", bin_i, N2hist_graph->GetX()[bin_i]) << endl;
			--bin_i;
		}
	}


	//CHECK what points are still in graphs:
	int leftpoints = N1hist_graph->GetN();
	//cout << "Still in N1hist_graph. (Cut lower part):" << endl;
	//cout << "----------------------------------------" << endl;
	//cout << leftpoints << " points left out of " << nbinsN1 << endl;
	for(int bin_i=0; bin_i!=N1hist_graph->GetN(); ++bin_i){
		//cout << "X: " << N1hist_graph->GetX()[bin_i] << endl;
	}
//}

	TMultiGraph *mg_cut = new TMultiGraph();
	//////TCanvas *mgcut_canv = new TCanvas("mgcut_canv", "mgcut_canv", 600, 400);
	mg_cut->Add(N1hist_graph);		//add cut graph to new multigraph
	mg_cut->Add(N2hist_graph);		//add cut graph to new multigraph

	mg_cut->SetTitle(Form("Selected Ranges of: N/N-1 and N/N-2 trigger-ratios; p_{T,avg}; Ratio of trg%d to (N-1) and (N-2)", curtrig)); //triggers.at(ind) instead of curtrig (when not separated to function like here)
	mg_cut->Draw("APZ");

	gStyle->SetOptStat(1);

	return mg_cut;	//return the ready-cut multigraph
}


//function to calculate uncertainty on resulting trigger-99%-threshold
//drawing same turn-on fct again, but with p0 shifted to p0+err_p0
//--> taking difference between original curve 99% point and new curve 99% point as error on that threshold
//do this once for each trigger --> get error
pair<Int_t, Double_t> FitUncCalc(Int_t trigger,TF1* turnon, Double_t pnom, Double_t pwid, Double_t pwid_err){
	TF1* new_turnon = new TF1(Form("new_turnon_%d", trigger), Form("0.5*(1+TMath::Erf((x-%f)/(%f)))", pnom, pwid+pwid_err), 20, 800);	//new turn-on curve with pwid shifted by its error
	Double_t orig99 = turnon->GetX(0.99);		//getting the old 99% trigger threshold
	Double_t new99 = new_turnon->GetX(0.99);	//getting the new 99% trigger threshold

	Double_t diff99 = TMath::Abs(orig99-new99);	//difference in threshold for triggers.at(ind)

	return make_pair(trigger, diff99);
}


//function to draw the summary of the fit parameters (pnom, pwid) --> not needed, is handled via SetDrawStyleWrite
//function to get the graph for the fit parameters
//usually called twice: once for pnom, once for pwid
//NOT NEEDED -- SIMPLY USE ADDRESS OF VECTOR TO CALL TGRAPHERRORS CONSTRUCTOR
//TGraphErrors ParamGraph(vector<Double_t> &par_x, vector<Double_t> &par_y, vector<Double_t> &par_err, ){}
	//create (draw) graph for checking the parameters p_nom and p_wid
	//make arrays to use TGraph constructor
	///Int_t ntrig = par_x.size();
	///Double_t xpar[ntrig];
	///Double_t ypar[ntrig];
	///for(int j=0; j!=ntrig; ++j){
		///xpar[j]=par_x.at(j);
		///ypar[j]=par_y.at(j);
	///}


//--------------------------------------------------------------------------------------------------//
//											 ANALYSIS LOOP											//
//--------------------------------------------------------------------------------------------------//



void DijetAnalysis::Loop(){	
	TH1::SetDefaultSumw2(kTRUE);		//in order to ALWAYS store the bin errors


	//insert analysis loop as it used to be in triggerpaths
	setTDRStyle();
	std::string name_h = "hpt";
	fChain->SetBranchStatus("*",0); 			//exclude / switch off all the branches
	fChain->SetBranchStatus("nJet",1); 			//include this branch in analysis	
	fChain->SetBranchStatus("Jet_pt",1); 		//include this branch as well
	fChain->SetBranchStatus("Jet_eta",1); 		//has to be included if eta-bins shall be checked...
	fChain->SetBranchStatus("Jet_phi",1);		//for TLorentzVector
	fChain->SetBranchStatus("Jet_mass",1);		//for TLorentzVector

	//for noise reduction (MET checks)
	fChain->SetBranchStatus("ChsMET_pt",1);		//MET in event
	fChain->SetBranchStatus("ChsMET_sumEt",1);	//Et sum in event
	fChain->SetBranchStatus("MET_pt",1);		//MET in event
	fChain->SetBranchStatus("MET_sumEt",1);	//Et sum in event


	if(montecarlo){
		cout << "This is a MC analysis. Switching on MC-specific branches." << endl;
		fChain->SetBranchStatus("Generator_weight",1);
		fChain->SetBranchStatus("nGenJet",1);
		fChain->SetBranchStatus("GenJet_pt",1);
		fChain->SetBranchStatus("GenJet_eta",1);
		fChain->SetBranchStatus("GenJet_phi",1);

		fChain->SetBranchStatus("Generator_binvar",1);	//for pthat checks
		fChain->SetBranchStatus("LHE_HT",1);			//for HT samples
	}
	else{//if data
		cout << Form("This is an analysis of dijet data. Using data from Run period: %s", runperiod)  << endl;	//output does not work...
	}

	//for additional trigger studies
	//Int_t j1=-1;	//triggerobject index of leading jet (arbitrary initial value)
	//Int_t j2=-1;	//triggerobject index of second leading jet	(arbitrary initival value)	--> is this okay?
	fChain->SetBranchStatus("nTrigObj",1);
	fChain->SetBranchStatus("TrigObj_id",1);
	fChain->SetBranchStatus("TrigObj_pt",1);


	//TFile* interim_output = new TFile("trg_lumi_output.root", "RECREATE"); //for storing the histograms
	//either set lumi here or not at all... currently lumi is only used in triggerpaths() anyway for normalisation of the final histogram

	///TFile* friendfile = TFile::Open("~/Documents/JetAnalysis/master_thesis_git/master_thesis/Tools/Uncertainties/mc_QCD_Pt-15to7000_Flat_0BD9CF76-F784-5C4C-896E-4A7D97B04661_all_jmeMCtest_singlefile100000_Friend.root");
	///TTree *friendtree = (TTree*)friendfile->Get("Friends");

	TFile* interim_output;
	if(montecarlo==false){//if data
		//check which runperiod
		//prepare output file (= after analysis BEFORE trigger analysis or creating plots)
		if(strcmp(runperiod,"A")==0){
			cout << "Run A" << endl;
			interim_output = new TFile("dijet_interim_output_data_RunA.root", "RECREATE");
		}//endif runA

		else if(strcmp(runperiod,"B")==0){
			cout << "Run B" << endl;
			interim_output = new TFile("dijet_interim_output_data_RunB.root", "RECREATE");
		}//endif runB

		else if(strcmp(runperiod,"C")==0){
			cout << "Run C" << endl;
			interim_output = new TFile("dijet_interim_output_data_RunC.root", "RECREATE");
		}//endif runC

		else if(strcmp(runperiod,"D")==0){
			cout << "Run D" << endl;
			interim_output = new TFile("dijet_interim_output_data_RunD.root", "RECREATE");
		}//endif runD

		else{//if no runperiod specified
			interim_output = new TFile("dijet_interim_output_data.root", "RECREATE"); //output file after analysis BEFORE trigger analysis or creating plots
		}
	}
	else{ //if MC
		if(strcmp(generator, "herwig")==0){
			interim_output = new TFile("dijet_interim_output_mc_herwig7.root", "RECREATE"); //output file if herwig7

			//fChain->AddFriend("Friends", friendfile);
			//fChain->Friends->SetBranchStatus("*", 0);

			//FriendClass(fFriends);
			//FriendClass::Init(fFriends);

			fFriends->SetBranchStatus("nJet", 0);
			//fFriends->SetBranchStatus("Jet_pt_raw", 1);
			//fFriends->SetBranchStatus("Jet_pt_nom", 1);
			fFriends->SetBranchStatus("Jet_pt_uncUp", 1);
			fFriends->SetBranchStatus("Jet_pt_uncDw", 1);
			fFriends->SetBranchStatus("Jet_pt*", 1);
		}
		else if(strcmp(generator, "pythia")==0){
			interim_output = new TFile("dijet_interim_output_mc_pythia8.root", "RECREATE"); //output file if pythia8
		}
		else if(strcmp(generator, "HTpythia")==0){
			cout << "This is a HT pythia8 sample." << endl;
			interim_output = new TFile("dijet_interim_output_mc_HT_madgraph_pythia8.root", "RECREATE"); //output file if pythia8 with MADGRAPH in HT bins
		}
		else{//if no generator specified
			interim_output = new TFile("dijet_interim_output_mc.root", "RECREATE"); //output file after analysis BEFORE trigger analysis or creating plots
		}
	}
	interim_output->cd();
	interim_output->mkdir("Standard"); //create directory in output file

	//create subfolder "Standard" in outputfile
	TDirectory *standard_folder = interim_output->GetDirectory("Standard"); 
	gROOT->GetListOfBrowsables()->Add(standard_folder, "Standard");
	//standard_folder->cd();


	//#### Preparing uncertainty analysis ###//
	//create subfolder "Uncertainties" in outputfile
	interim_output->cd();
	interim_output->mkdir("Uncertainties");
	TDirectory *unc_folder = interim_output->GetDirectory("Uncertainties");
	gROOT->GetListOfBrowsables()->Add(unc_folder, "Uncertainties");

	//create these always, but only use them in case of mc herwig
	/*
	TH1D* jt_ptraw_hist;
	TH1D* jt_ptnom_hist;
	TH1D* jt_corrJEC_hist;
	TH1D* jt_corrJER_hist;
	TH1D* jt_ptjerUp_hist;
	TH1D* jt_ptjerDw_hist;
	if(strcmp(generator, "herwig")==0){	//additional folder if herwig
		interim_output->cd();
		interim_output->mkdir("jmeHelperRun2");
		TDirectory *jme_folder = interim_output->GetDirectory("jmeHelperRun2");
		gROOT->GetListOfBrowsables()->Add(jme_folder, "jmeHelperRun2");

		//go later through ybys bins and call ybys_jmedirs() function to set directories
	}
	*/

	standard_folder->cd();	//change back to Standard directory

	//const char *a = algo.c_str();

	//create vector containing string names of uncertainty sources
	//these names have been extracted from Autumn18_V19_MC_UncertaintySources_AK4PF.txt
	//but should be the same names also for chs etc.
	vector<string> uncnames;
	uncnames.push_back("AbsoluteStat");					//ind:  0
	uncnames.push_back("AbsoluteScale");				//ind:  1
	uncnames.push_back("AbsoluteSample");				//ind:  2
	uncnames.push_back("AbsoluteFlavMap");				//ind:  3
	uncnames.push_back("AbsoluteMPFBias");				//ind:  4
	uncnames.push_back("Fragmentation");				//ind:  5
	uncnames.push_back("SinglePionECAL");				//ind:  6
	uncnames.push_back("SinglePionHCAL");				//ind:  7
	uncnames.push_back("FlavorQCD");					//ind:  8
	uncnames.push_back("TimePtEta");					//ind:  9
	uncnames.push_back("RelativeJEREC1");				//ind: 10
	uncnames.push_back("RelativeJEREC2");				//ind: 11
	uncnames.push_back("RelativeJERHF");				//ind: 12
	uncnames.push_back("RelativePtBB");					//ind: 13
	uncnames.push_back("RelativePtEC1");				//ind: 14
	uncnames.push_back("RelativePtEC2");				//ind: 15
	uncnames.push_back("RelativePtHF");
	uncnames.push_back("RelativeBal");
	uncnames.push_back("RelativeSample");
	uncnames.push_back("RelativeFSR");
	uncnames.push_back("RelativeStatFSR");				//ind: 20
	uncnames.push_back("RelativeStatEC");
	uncnames.push_back("RelativeStatHF");
	uncnames.push_back("PileUpDataMC");
	uncnames.push_back("PileUpPtRef");
	uncnames.push_back("PileUpPtBB");					//ind: 25
	uncnames.push_back("PileUpPtEC1");
	uncnames.push_back("PileUpPtEC2");
	uncnames.push_back("PileUpPtHF");
	uncnames.push_back("PileUpMuZero");
	uncnames.push_back("PileUpEnvelope");				//ind: 30
	uncnames.push_back("SubTotalPileUp");
	uncnames.push_back("SubTotalRelative");
	uncnames.push_back("SubTotalPt");
	uncnames.push_back("SubTotalScale");
	uncnames.push_back("SubTotalAbsolute");				//ind: 35
	uncnames.push_back("SubTotalMC");
	uncnames.push_back("Total");
	uncnames.push_back("TotalNoFlavor");
	uncnames.push_back("TotalNoTime");
	uncnames.push_back("TotalNoFlavorNoTime");			//ind: 40
	uncnames.push_back("FlavorZJet");
	uncnames.push_back("FlavorPhotonJet");
	uncnames.push_back("FlavorPureGluon");
	uncnames.push_back("FlavorPureQuark");
	uncnames.push_back("FlavorPureCharm");				//ind: 45
	uncnames.push_back("FlavorPureBottom");
	uncnames.push_back("TimeRunA");
	uncnames.push_back("TimeRunB");
	uncnames.push_back("TimeRunC");
	uncnames.push_back("TimeRunD");						//ind: 50
	uncnames.push_back("CorrelationGroupMPFInSitu");
	uncnames.push_back("CorrelationGroupIntercalibration");
	uncnames.push_back("CorrelationGroupbJES");
	uncnames.push_back("CorrelationGroupFlavor");
	uncnames.push_back("CorrelationGroupUncorrelated");	//ind: 55

	//make a map (56 branches up?) -->  does not make sense to create this map for each event. 
	/*
	map<string, TBranch> friend_branches_map;
	friend_branches_map.insert(make_pair(uncnames.at(0), Friends->Jet_pt_jesAbsoluteStatUp));
	friend_branches_map.insert(make_pair(uncnames.at(1), Friends->Jet_pt_jesAbsoluteScaleUp));
	friend_branches_map.insert(make_pair(uncnames.at(2), Friends->Jet_pt_jesAbsoluteSampleUp));
	friend_branches_map.insert(make_pair(uncnames.at(3), Friends->Jet_pt_jesAbsoluteFlavMapUp));
	friend_branches_map.insert(make_pair(uncnames.at(4), Friends->Jet_pt_jesAbsoluteMPFBiasUp));
	friend_branches_map.insert(make_pair(uncnames.at(5), Friends->Jet_pt_jesFragmentationUp));
	friend_branches_map.insert(make_pair(uncnames.at(6), Friends->Jet_pt_jesSinglePionECALUp));
	friend_branches_map.insert(make_pair(uncnames.at(7), Friends->Jet_pt_jesSinglePionHCALUp));
	friend_branches_map.insert(make_pair(uncnames.at(8), Friends->Jet_pt_jesFlavorQCDUp));
	friend_branches_map.insert(make_pair(uncnames.at(9), Friends->Jet_pt_jesTimePtEtaUp));
	friend_branches_map.insert(make_pair(uncnames.at(10), Friends->Jet_pt_jesRelativeJEREC1Up));
	friend_branches_map.insert(make_pair(uncnames.at(11), Friends->Jet_pt_jesRelativeJEREC2Up));
	friend_branches_map.insert(make_pair(uncnames.at(12), Friends->Jet_pt_jesRelativeJERHFUp));
	friend_branches_map.insert(make_pair(uncnames.at(13), Friends->Jet_pt_jesRelativePtBBUp));
	friend_branches_map.insert(make_pair(uncnames.at(14), Friends->Jet_pt_jesRelativePtEC1Up));
	friend_branches_map.insert(make_pair(uncnames.at(15), Friends->Jet_pt_jesRelativePtEC2Up));
	friend_branches_map.insert(make_pair(uncnames.at(16), Friends->Jet_pt_jesRelativePtHFUp));
	friend_branches_map.insert(make_pair(uncnames.at(17), Friends->Jet_pt_jesRelativeBalUp));
	friend_branches_map.insert(make_pair(uncnames.at(18), Friends->Jet_pt_jesRelativeSampleUp));
	friend_branches_map.insert(make_pair(uncnames.at(19), Friends->Jet_pt_jesRelativeFSRUp));
	friend_branches_map.insert(make_pair(uncnames.at(20), Friends->Jet_pt_jesRelativeStatFSRUp));
	friend_branches_map.insert(make_pair(uncnames.at(21), Friends->Jet_pt_jesRelativeStatECUp));
	friend_branches_map.insert(make_pair(uncnames.at(22), Friends->Jet_pt_jesRelativeStatHFUp));
	friend_branches_map.insert(make_pair(uncnames.at(23), Friends->Jet_pt_jesPileUpDataMCUp));
	friend_branches_map.insert(make_pair(uncnames.at(24), Friends->Jet_pt_jesPileUpPtRefUp));
	friend_branches_map.insert(make_pair(uncnames.at(25), Friends->Jet_pt_jesPileUpPtBBUp));
	friend_branches_map.insert(make_pair(uncnames.at(26), Friends->Jet_pt_jesPileUpPtEC1Up));
	friend_branches_map.insert(make_pair(uncnames.at(27), Friends->Jet_pt_jesPileUpPtEC2Up));
	friend_branches_map.insert(make_pair(uncnames.at(28), Friends->Jet_pt_jesPileUpPtHFUp));
	friend_branches_map.insert(make_pair(uncnames.at(29), Friends->Jet_pt_jesPileUpMuZeroUp));
	friend_branches_map.insert(make_pair(uncnames.at(30), Friends->Jet_pt_jesPileUpEnvelopeUp));
	friend_branches_map.insert(make_pair(uncnames.at(31), Friends->Jet_pt_jesSubTotalPileUpUp));
	friend_branches_map.insert(make_pair(uncnames.at(32), Friends->Jet_pt_jesSubTotalRelativeUp));
	friend_branches_map.insert(make_pair(uncnames.at(33), Friends->Jet_pt_jesSubTotalPtUp));
	friend_branches_map.insert(make_pair(uncnames.at(34), Friends->Jet_pt_jesSubTotalScaleUp));
	friend_branches_map.insert(make_pair(uncnames.at(35), Friends->Jet_pt_jesSubTotalAbsoluteUp));
	friend_branches_map.insert(make_pair(uncnames.at(36), Friends->Jet_pt_jesSubTotalMCUp));
	friend_branches_map.insert(make_pair(uncnames.at(37), Friends->Jet_pt_jesTotalUp));
	friend_branches_map.insert(make_pair(uncnames.at(38), Friends->Jet_pt_jesTotalNoFlavorUp));
	friend_branches_map.insert(make_pair(uncnames.at(39), Friends->Jet_pt_jesTotalNoTimeUp));
	friend_branches_map.insert(make_pair(uncnames.at(40), Friends->Jet_pt_jesTotalNoFlavorNoTime));
	friend_branches_map.insert(make_pair(uncnames.at(41), Friends->Jet_pt_jesFlavorZJetUp));
	friend_branches_map.insert(make_pair(uncnames.at(42), Friends->Jet_pt_jesFlavorPhotonJetUp));
	friend_branches_map.insert(make_pair(uncnames.at(43), Friends->Jet_pt_jesFlavorPureGluonUp));
	friend_branches_map.insert(make_pair(uncnames.at(44), Friends->Jet_pt_jesFlavorPureQuarkUp));
	friend_branches_map.insert(make_pair(uncnames.at(45), Friends->Jet_pt_jesFlavorPureCharmUp));
	friend_branches_map.insert(make_pair(uncnames.at(46), Friends->Jet_pt_jesFlavorPureBottomUp));
	friend_branches_map.insert(make_pair(uncnames.at(47), Friends->Jet_pt_jesTimeRunAUp));
	friend_branches_map.insert(make_pair(uncnames.at(48), Friends->Jet_pt_jesTimeRunBUp));
	friend_branches_map.insert(make_pair(uncnames.at(49), Friends->Jet_pt_jesTimeRunCUp));
	friend_branches_map.insert(make_pair(uncnames.at(50), Friends->Jet_pt_jesTimeRunDUp));
	friend_branches_map.insert(make_pair(uncnames.at(51), Friends->Jet_pt_jesCorrelationGroupMPFInSituUp));
	friend_branches_map.insert(make_pair(uncnames.at(52), Friends->Jet_pt_jesCorrelationGroupIntercalibrationUp));
	friend_branches_map.insert(make_pair(uncnames.at(53), Friends->Jet_pt_jesCorrelationGroupbJESUp));
	friend_branches_map.insert(make_pair(uncnames.at(54), Friends->Jet_pt_jesCorrelationGroupFlavorUp));
	friend_branches_map.insert(make_pair(uncnames.at(55), Friends->Jet_pt_jesCorrelationGroupUncorrelatedUp));
	*/

	/*
	//create uncertainty sources vector (see jecsys package)
	vector<JetCorrectionUncertainty*> uncs_vec(uncnames.size());
	for(unsigned int k=0; k!=uncnames.size(); ++k){
		//file where the uncs are stored
		string sfilename = Form("/Documents/JetAnalysis/unc_JEC_files/Autumn18_V19_MC/Autumn18_V19_MC_UncertaintySources_AK4PF.txt");
		//or with choice of algorithm
		//string sfilename = Form("/Documents/JetAnalysis/unc_JEC_files/Autumn18_V19_MC/Autumn18_V19_MC_UncertaintySources_AK4PF%s.txt", algo);

		const char* sfname = sfilename.c_str();		//.txt file name now as a char
		const char* src	= uncnames.at(k).c_str();		//uncertainty source name as char, snames.[k].c_str()

		cout << "sfname: " << sfname << endl;
		cout << "src: " << src << endl;

		JetCorrectorParameters *p = new JetCorrectorParameters(sfname, src);
		JetCorrectionUncertainty *unc = new JetCorrectionUncertainty(*p);
		uncs_vec[k] = unc;	//set value of the uncertainty vector for current unc source
	}
	*/

	//### end of preparing unc analysis ###//



	Int_t genbins = 106;
	//Double_t genmin = fChain->GetMinimum("Generator_weight");
	//Double_t genmax = fChain->GetMaximum("Generator_weight");
	Double_t genmin = -1.06;
	Double_t genmax = 1.06;

	TH1D* mc_genweight = new TH1D("mc_genweight", "Distribution before any selection; Generator_weight; entries", genbins, genmin, genmax);	//for normalising later - histo with generator weights
	mc_genweight->SetDirectory(standard_folder);	//this histo is only created once for the whole data set


	//could shorten this vector (starting from 40), but need to adjust trigger indexing everywhere...
	vector<int> triggers = {0, 15, 25, 40, 60, 80, 140, 200, 260, 320, 400, 450, 500, 550};

	//lumivec should be taken from an extra file!!!!
	///vector<Double_t> lumivec = {0, 0, 0.000414643, 0.001737711, 0.006541972, 0.097552471, 0.418082018, 0.868855667, 3.475422669, 6.950845338, 13.901690676, 111.213525407, 111.213525407};	//now for testing --> should be more automatised / generalised [units are /pb]

	//new lumi (from own json files, but after skimming, used brilcalc) units are /fb
	//vector<Double_t> lumivec = {};

	if(montecarlo==false){
		//switching on selected TBranches corresponding to chosen triggers:
		for (int j=1; j!=triggers.size(); ++j){ //start from 1, because "0" means "without trigger"
			char *trigname = Form("HLT_PFJet%i", triggers.at(j));
			fChain->SetBranchStatus(trigname, 1);
		}
	}

	YbYs *ybys_pointer;	//should maybe handle this differently
	//call constructor YbYs(float yb_lo, float ys_lo, float yb_up, float ys_up)
	//ybys_pointer = new YbYs[6]{{0.0,0.0,1.0,1.0}, {0.0, 1.0, 1.0, 2.0}, {0.0, 2.0, 1.0, 3.0}, {1.0, 0.0, 2.0, 1.0}, {1.0, 1.0, 2.0, 2.0}, {2.0, 0.0, 3.0, 1.0}};	
	//new outermost bin bounds --> go only up to 2.4 for both yb, ys
	//ybys_pointer = new YbYs[6]{{0.0,0.0,1.0,1.0}, {0.0, 1.0, 1.0, 2.0}, {0.0, 2.0, 1.0, 2.4}, {1.0, 0.0, 2.0, 1.0}, {1.0, 1.0, 2.0, 2.0}, {2.0, 0.0, 2.4, 1.0}}; //new, without unc-directory

	//new and contains unc_folder stuff
	ybys_pointer = new YbYs[6]{{0.0,0.0,1.0,1.0, uncnames}, {0.0, 1.0, 1.0, 2.0, uncnames}, {0.0, 2.0, 1.0, 2.4, uncnames}, {1.0, 0.0, 2.0, 1.0, uncnames}, {1.0, 1.0, 2.0, 2.0, uncnames}, {2.0, 0.0, 2.4, 1.0, uncnames}}; 

	cout << "this is a test." << endl;


	vector <Histos*> histosvec;
	for (int ybysnr=0; ybysnr!=6; ++ybysnr){
		YbYs cur_ybysbin = ybys_pointer[ybysnr];
		if(montecarlo==false){//only if data==true --> montecarlo==false
			histosvec.push_back(new Histos(triggers));
			Histos curtrig = *histosvec.back();
			YbYs cur_ybysbin = ybys_pointer[ybysnr];
			cur_ybysbin.ybys_directories(standard_folder, triggers, curtrig);	//creates folder structure
			cur_ybysbin.ybys_uncdirs(unc_folder, uncnames);						//setting the unc histos to the correct folder
		}
		else{//if MC
			cur_ybysbin.ybys_directories(standard_folder);
			cur_ybysbin.ybys_uncdirs(unc_folder, uncnames);						//same procedure for mc
		}//if MC
	}

	//counters
	int count_nodijet = 0; 		//count number of events that are no dijet events, i.e. nJet < 2
	int count_NoBin_ybys = 0; 	//count number of events that have no ybys bin, where they could be assigned to (i.e. yb or ys > 3.0)
	int count_BadRapidity = 0;	//count number of events where the jets do not fulfill the rapidity requirements
	int count_selected = 0;		//count number of selected events: implemented for MC
	int count_noGen1 = 0;		//no matching gen jet for leading jet
	int count_noGen2 = 0;		//leading jet has matching gen jet, but second leading jet does not...
	int count_jetmatch = 0;		//count events with two matched reco and gen leading jets

	int count_pthat_fail_100_400 = 0;	//new test



	//the actual Loop()
	if (fChain == 0) return;
	Long64_t nentries = fChain->GetEntriesFast();
	cout << "nentries = " << nentries << endl << endl;			//information

	Long64_t nbytes = 0, nb = 0;
	//starting event loop
	for (Long64_t jentry=0; jentry<nentries;jentry++) {
	////for (Long64_t jentry=0; jentry!=100000; ++jentry){} //for testing with less entries/events

		Long64_t ientry = LoadTree(jentry);
		if(v==true){cout << "current entry: " << jentry << endl;}
		///if (jentry%10==0){cout << "------------------------------------------" << endl;}//for overview
		if (jentry%100000==0){cout << "processed " << jentry << " of " << nentries << endl;}
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);   nbytes += nb;
		
		map<const int, Bool_t> tdecisions;

		if(montecarlo==false){
			tdecisions.insert(make_pair(triggers.at(0), true)); //-->no trigger
			tdecisions.insert(make_pair(triggers.at(1), HLT_PFJet15)); //could this be handled in a more general way? (see attempt and problem above)
			tdecisions.insert(make_pair(triggers.at(2), HLT_PFJet25));
			tdecisions.insert(make_pair(triggers.at(3), HLT_PFJet40));
			tdecisions.insert(make_pair(triggers.at(4), HLT_PFJet60));
			tdecisions.insert(make_pair(triggers.at(5), HLT_PFJet80));
			tdecisions.insert(make_pair(triggers.at(6), HLT_PFJet140));
			tdecisions.insert(make_pair(triggers.at(7), HLT_PFJet200));
			tdecisions.insert(make_pair(triggers.at(8), HLT_PFJet260));
			tdecisions.insert(make_pair(triggers.at(9), HLT_PFJet320));
			tdecisions.insert(make_pair(triggers.at(10), HLT_PFJet400));
			tdecisions.insert(make_pair(triggers.at(11), HLT_PFJet450));
			tdecisions.insert(make_pair(triggers.at(12), HLT_PFJet500));
			tdecisions.insert(make_pair(triggers.at(13), HLT_PFJet550));
		}
		else{
			mc_genweight->Fill(Generator_weight);		//fill this histo with generator weight of current event, before any cuts!
		}

		Int_t j1=-1;	//triggerobject index of leading jet (arbitrary initial value)
		Int_t j2=-1;	//triggerobject index of second leading jet	(arbitrary initival value)	--> is this okay?


		//check if Event (jentry) is a dijet event. --> otherwise do nothing
		if(nJet>1){	//required for dijet ! --> only look at two leading jets

			//add MET check for reducing noise
			//bool good_met = ChsMET_pt<(0.3*ChsMET_sumEt);		 	//Chs variable not existent in MC sets
			bool good_met = MET_pt<(0.3*MET_sumEt);			//exists in both data and mc

			//check pthat in case of MC (set it to true for data)
			bool good_pthat;	//change this name to good_binvar in the future! (as when looking at HT samples, this is not pthat, but HT...)
			//if(montecarlo&&(Generator_binvar<100)){}
			if(montecarlo){//do this always, for every event
				//if(Jet_pt[0]<1.5*Generator_binvar){good_pthat=true;}
				if(Jet_pt[0]<2.0*Generator_binvar){good_pthat=true;}		//looser than 1.5 criterion
				else if(strcmp(generator, "HTpythia")==0){
					if(Jet_pt[0]<2.0*(LHE_HT/2.0)){good_pthat=true;}	//HT sample does have "Generator_binvar=-1" for each entry...
					else{good_pthat=false;}
				}
				else{good_pthat=false;}
			}
			else{//if data
				good_pthat=true;
			}


			Double_t j1_pt	= Jet_pt[0];
			Double_t j2_pt	= Jet_pt[1];
			Double_t j1_eta	= Jet_eta[0];
			Double_t j2_eta	= Jet_eta[1];
			Double_t j1_phi	= Jet_phi[0];
			Double_t j2_phi	= Jet_phi[1];
			Double_t j1_mass = Jet_mass[0];
			Double_t j2_mass = Jet_mass[1];
			//use Lorentzvector:
			PtEtaPhiMVector jet1_vec(j1_pt, j1_eta, j1_phi, j1_mass);
			PtEtaPhiMVector jet2_vec(j2_pt, j2_eta, j2_phi, j2_mass);

			//calculate current ptavg, yb, ys
			Double_t j1_y = jet1_vec.Rapidity();
			Double_t j2_y = jet2_vec.Rapidity();
			Double_t cur_ptavg = 0.5*(Jet_pt[0]+Jet_pt[1]);
			Double_t cur_yboost = 0.5*TMath::Abs(j1_y+j2_y);
			Double_t cur_ystar = 0.5*TMath::Abs(j1_y-j2_y);

			if((good_pthat==false)&&(cur_ptavg<400)&&(cur_ptavg>100)){count_pthat_fail_100_400++;}	//TEST
			//if(good_pthat==false){count_pthat_fail_100_400++;}	//TEST


			//additional variables that will only be used for HERWIG 7 (same eta, phi and mass as jet1 and jet2, but pt is smeared)
			PtEtaPhiMVector jet1_smeared(Jet_pt_nom[0], j1_eta, j1_phi, j1_mass);
			PtEtaPhiMVector jet2_smeared(Jet_pt_nom[1], j2_eta, j2_phi, j2_mass);


			//if MC also gen-ptavg and other gen properties
			bool matched=false;		//if matching gen jet for both leading jets found
			//the three observables in gen
			Double_t gen_ptavg;
			Double_t gen_yboost;
			Double_t gen_ystar;
			//more helper variables in gen (rapidity) and Lorentzvectors
			Double_t gen1_y;
			Double_t gen2_y;
			PtEtaPhiMVector gen1_vec;
			PtEtaPhiMVector gen2_vec;
			//other bools on gen
			bool good_genrapidity;
			bool good_genpt;

			//add offline cut on Jet-Rapidity --> both jets must have Abs(y)<5.0 and at least one shall have Abs(y)<2.5
			Double_t absYj1 = TMath::Abs(j1_y);
			Double_t absYj2 = TMath::Abs(j2_y);
			//bool good_rapidity = (absYj1<5)&&(absYj2<5)&&(absYj1<2.5 || absYj2<2.5);	//old cuts
			bool good_rapidity = (absYj1<2.4)&&(absYj2<2.4);	//NEW cuts!

			//add offline cut on Jet-pT --> both jets must fulfill pT > 30 GeV.
			//bool good_pt = (j1_pt>30)&&(j2_pt>30);
			bool good_pt = (j1_pt>20)&&(j2_pt>20);		//set down to 20 pt; check strict criterion later

			//for additional trigger studies
			//loop over all the trigger objects in current event
			//find the two hardest trigger objects in this event and store their trigger object index (as j1, j2)
			//cout << "nTrigObj: " << nTrigObj << endl;
			//use pt
			Double_t pt1=0;
			Double_t pt2=0;
			if(montecarlo==false){
				for(int obj_ind=0; obj_ind!=nTrigObj; ++obj_ind){
					//cout << "in loop. obj_ind=" << obj_ind << endl;
					//check that the ID corresponds to Jet (TrigObj_id == 1)
					if(TrigObj_id[obj_ind]==1){
						//if(j1==-1){j1=obj_ind;}	//setting initial j1
						//how shall j2 be set initially?	--> arbitrarily ... will change anyway in case it's wrongly assigned.
						Double_t cur_objpt = TrigObj_pt[obj_ind];	//pt of current trigger object
						//if(cur_objpt>TrigObj_pt[j1]){}
						//if(j1==-1 || cur_objpt>TrigObj_pt[j1]){}
						if(cur_objpt>pt1){
							j2 = j1;
							j1 = obj_ind;	//current trigger object harder than previously found hardest one
							pt2 = pt1;
							pt1 = cur_objpt;
							
							//cout << "changed j1 to current obj_ind." << endl;
						}
						//else if((cur_objpt<TrigObj_pt[j1])&&(cur_objpt>TrigObj_pt[j2])){}
						//else if(cur_objpt>TrigObj_pt[j2]){}
						else if(cur_objpt>pt2){
							pt2 = cur_objpt;
							j2 = obj_ind;
							//cout << "changed j2 to current obj_ind." << endl;
						}
						//else{
						//	continue;	//current trigger object had smaller pt than the ones investigated before
						//}
					}
				
					else{	//this was not a jet (does not have correct ID)
						//cout << "not a jet." << endl;
						//continue;
					}
				}
			}//if MC==false --> IF DATA
			else{ //if MC==true --> do jet matching
				Int_t gen1_ind=-1;		//index of leading genjet
				Int_t gen2_ind=-1;		//index of second genjet
				//go through all the genjets and find one that matches with leading jet
				for(int gen_ind=0;gen_ind!=nGenJet;++gen_ind){
					//Double_t genjet_eta = GenJet_eta[gen_ind];
					//Double_t genjet_phi = GenJet_phi[gen_ind];
					Double_t del_Eta_sqr1 = TMath::Power(TMath::Abs(j1_eta-GenJet_eta[gen_ind]), 2);	//squared
					Double_t del_Phi_sqr1 = TMath::Power(TMath::Abs(j1_phi-GenJet_phi[gen_ind]), 2);

					Double_t delR_j1 = TMath::Sqrt(del_Eta_sqr1+del_Phi_sqr1);	//check current genjet with jet1
					if(delR_j1<0.2){
						gen1_ind = gen_ind;	
						//now search matching genjet to jet2
						for(int g=0; g!=nGenJet; ++g){
							Double_t del2_Eta = TMath::Power(TMath::Abs(j2_eta-GenJet_eta[g]), 2);
							Double_t del2_Phi = TMath::Power(TMath::Abs(j2_phi-GenJet_phi[g]), 2);

							Double_t delR_j2 = TMath::Sqrt(del2_Eta+del2_Phi);	//check current genjet with jet2
							if(delR_j2<0.2){
								gen2_ind = g;
								if(gen1_ind!=gen2_ind){	//MATCHED!
									matched = true;
									gen1_vec = PtEtaPhiMVector(GenJet_pt[gen1_ind], GenJet_eta[gen1_ind], GenJet_phi[gen1_ind], GenJet_mass[gen1_ind]);
									gen2_vec = PtEtaPhiMVector(GenJet_pt[gen2_ind], GenJet_eta[gen2_ind], GenJet_phi[gen2_ind], GenJet_mass[gen2_ind]);
									gen1_y = gen1_vec.Rapidity();
									gen2_y = gen2_vec.Rapidity();

									gen_ptavg = 0.5*(GenJet_pt[gen1_ind]+GenJet_pt[gen2_ind]);
									gen_yboost = 0.5*TMath::Abs(gen1_y+gen2_y);
									gen_ystar = 0.5*TMath::Abs(gen1_y-gen2_y);
									good_genrapidity = (TMath::Abs(gen1_y)<2.4)&&(TMath::Abs(gen2_y)<2.4);
									good_genpt = (gen1_vec.pt()>30)&&(gen2_vec.pt()>30);

									count_jetmatch++;	//matched both leading recojets to genjets!
								}
								break; //found matching genjet
							}//end if(delR_j2<0.2)
							else if(g==nGenJet-1){ //did not find matching genjet for second jet
								count_noGen2++;//second leading jet has no matching gen jet...
							}
						}//end of second genjet loop
					}//end if(delR_j1<0.2)
					else if(gen_ind==nGenJet-1){
						count_noGen1++;//leading jet has no matching gen jet...
					}
				}//end of first genjet loop
			}//end if MC==true


			//check if rapidity requirement is fulfilled
			//check at the same time whether pT requirement is fulfilled
			//generator weight
			if(good_rapidity&&good_pt&&good_met&&good_pthat){	//now two new extra checks: good_met and good_pthat
			//if(good_rapidity&&good_pt){}

				//set the generator weight (has to be done already here, because of different procedure for HT samples
				Double_t gen_weight;
				if(strcmp(generator, "HTpythia")==0){//HT sample
					Double_t xsHT;	//this values are taken from the jetphys settings.h_template file and are (as it seems) in pb!
					Double_t nevHT; //sum of event weights --> does not have to be an integer!
					if((50<=LHE_HT)&&(LHE_HT<100)){
						//nevHT = 1319894;			//no. of generated entries
						nevHT = 1319844.8;			//sum of event weights
						//xsHT = 23700000.0;
						///xsHT = 195476190.0;		//calculated with parabol fit on 2016 over 2018 XS from jetphys settings.h_template file...
						xsHT = 19380000;	//commented out in jetphys

						gen_weight = xsHT/nevHT;
					}
					else if((100<=LHE_HT)&&(LHE_HT<200)){
						//nevHT = 1474667;			//no. of generated entries
						nevHT = 1474261.5;	//sum of event weights
						///xsHT = 23700000.0;
						xsHT = 19380000;	//commented out in jetphys

						gen_weight = xsHT/nevHT;
					}
					else if((200<=LHE_HT)&&(LHE_HT<300)){
						//nevHT = 2055720;			//no. of generated entries
						nevHT = 2054219.0;			//sum of event weights
						///xsHT = 1547000.0;
						xsHT = 1559000;	//commented out in jetphys
						gen_weight = xsHT/nevHT;
					}
					else if((300<=LHE_HT)&&(LHE_HT<500)){
						//nevHT = 1251844;			//no. of generated entries
						nevHT = 1250271.9;			//sum of event weights
						///xsHT = 322600.0;
						xsHT = 311900;	//commented out in jetphys
						gen_weight = xsHT/nevHT;
					}
					else if((500<=LHE_HT)&&(LHE_HT<700)){
						//nevHT = 1048024;			//no. of generated entries
						nevHT = 1045993.6;			//sum of event weights
						///xsHT = 29980.0;
						xsHT = 29070;	//commented out in jetphys
						gen_weight = xsHT/nevHT;
					}
					else if((700<=LHE_HT)&&(LHE_HT<1000)){
						//nevHT = 1225358;			//no. of generated entries
						nevHT = 1222024.0;			//sum of event weights
						///xsHT = 6334.0;
						xsHT = 5962;	//commented out in jetphys
						gen_weight = xsHT/nevHT;
					}
					else if((1000<=LHE_HT)&&(LHE_HT<1500)){
						//nevHT = 1465418;			//no. of generated entries
						nevHT = 1459472.4;			//sum of event weights
						///xsHT = 1088.0;
						xsHT = 1005;	//commented out in jetphys
						gen_weight = xsHT/nevHT;
					}
					else if((1500<=LHE_HT)&&(LHE_HT<2000)){
						//nevHT = 1031687;			//no. of generated entries
						nevHT = 1024810.8;			//sum of event weights
						///xsHT = 99.11;
						xsHT = 101.8;	//commented out in jetphys
						gen_weight = xsHT/nevHT;
					}
					else if(2000<=LHE_HT){
						//nevHT = 186052;			//no. of generated entries
						nevHT = 183974.07;			//sum of event weights
						///xsHT = 20.23;
						xsHT = 20.54;	//commented out in jetphys
						gen_weight = xsHT/nevHT;
					}
					else{
						cout << "Does this not have LHE_HT ???" << endl;
					}
				}//end: if ht sample
				else{
					gen_weight = Generator_weight;		//MC event weight to be given to fillhistos() function
				}
				//Double_t gen_weight = Generator_weight;		//MC event weight to be given to fillhistos() function
				//check ybys bin --> put this into function!!
				for (int ybysnr=0; ybysnr!=6; ++ybysnr){
					YbYs cur_ybysbin = ybys_pointer[ybysnr];
					float yb_lo = cur_ybysbin.yb_get_bb().first;
					float yb_up = cur_ybysbin.yb_get_bb().second;
					float ys_lo = cur_ybysbin.ys_get_bb().first;
					float ys_up = cur_ybysbin.ys_get_bb().second;
					if(v==true){cout << "yb_lo: " << yb_lo << " ## yb_up: " << yb_up << " ## ys_lo: " << ys_lo << " ## ys_up: " << ys_up << endl;}
					if((yb_lo<=cur_yboost && cur_yboost<yb_up) and (ys_lo<=cur_ystar && cur_ystar<ys_up)){	//-->better to implement function InInterval(value, low, up) that checks if value in interval... return(lo<=x & x<up)
						if(v==true){cout << "Event is in ybys bin: " << cur_ybysbin.get_bname() << endl;}

						if(montecarlo==false){//only for data:
							//check the pT > 30GeV cut:
							if((j1_pt>30)&&(j2_pt>30)){
								//check trigger
								for(int trignr=0; trignr!=triggers.size(); ++trignr){
									if(tdecisions[triggers.at(trignr)]){	//==true
										//fill into histo of current trigger
										histosvec.at(ybysnr)->jt_histos[triggers.at(trignr)]->Fill(cur_ptavg);	//fill ptavg trigger histogram

										//if(trignr>2){}	//this would count them several times for each trigger path the event fired...
										if(trignr==0){//count only once.
											count_selected++; //count this event (in case trigger 40 or higher)
										}

										//check if at least one of the leading jets' pt is higher than next trigger-threshold
										if((trignr!=triggers.size()-1) && (j1 != -1) && (TrigObj_pt[j1] >= triggers.at(trignr+1))){
										//if((TrigObj_pt[j1] >= triggers.at(trignr+1)) || (TrigObj_pt[j2] >= triggers.at(trignr+1))){}
											//fill "simulated" histo of next higher trigger by applying additional condition
											histosvec.at(ybysnr)->trgobj_histos[triggers.at(trignr+1)]->Fill(cur_ptavg);
										}
										//check if at least one of the leading jets' pt is higher than the next next trigger-threshold (N+2)
										if((trignr<triggers.size()-2) && (j1 != -1) && (TrigObj_pt[j1] >= triggers.at(trignr+2))){
											histosvec.at(ybysnr)->n2_trgobj_histos[triggers.at(trignr+2)]->Fill(cur_ptavg);
										}
										continue;
									}//trigger true?
									else{//trigger false
										continue;
									}//trigger false
								}//check triggers
							}//endif both pt>30	 (is else --> continue needed?)
						}//end if DATA (if mc==false)
						else if(montecarlo==true){
							if((j1_pt>30)&&(j2_pt>30)){
								cur_ybysbin.get_mc_histo()->Fill(cur_ptavg, gen_weight);	//apply generator weight (as this is MC) (might contain unmatched ones)
							}

							if(strcmp(generator, "herwig")==0){
								if((jet1_smeared.pt()>30)&&(jet2_smeared.pt()>30)){
									Double_t smeared_ptavg = 0.5*(jet1_smeared.pt()+jet2_smeared.pt());
									cur_ybysbin.get_mc_histo_smearedJet()->Fill(smeared_ptavg, gen_weight);
								}
								//need to "loop through the branches" of the Friends TTree -- don't have to loop. just take the values. (from the map)
								//The branch usually has the same name as the Unc-histo, just added an "Up" or "Down"
								for(unsigned int k=0; k!=uncnames.size(); ++k){
									/*
									if(jentry%10000==0){
										cout << "k==" << k << endl;
										cout << "unc = " << uncnames.at(k) << endl;
										//cout << "uncs_vec.at(k) = " << uncs_vec.at(k) << endl;
									}
									*/

									//Get the corresponding branches
									//these already contain the varied Jet_pt value
									const char* unc_up_branchname = Form("Jet_pt_jes%sUp", uncnames.at(k).c_str());
									const char* unc_down_branchname = Form("Jet_pt_jes%sDown", uncnames.at(k).c_str());


									//Get the variations
									///double newpt1 = friend_branches_map[uncnames.at(k)][0];
									///double newpt2 = friend_branches_map[uncnames.at(k)][1];

									//double newpt1 = fFriends->Jet_pt_uncUp[k][0];
									//double newpt2 = fFriends->Jet_pt_uncUp[k][1];
									double newpt1_up = Jet_pt_uncUp[k][0];
									double newpt2_up = Jet_pt_uncUp[k][1];
									double newpt1_down = Jet_pt_uncDw[k][0];
									double newpt2_down = Jet_pt_uncDw[k][1];


									//how to get new yb ys values?? 
									//new ptavg
									double newptavg_up = 0.5*(newpt1_up+newpt2_up);
									double newptavg_down = 0.5*(newpt1_down+newpt2_down);

									/*
									if(jentry%10000==0){		//testing
										cout << "old pt1: " << jet1_vec.pt() << endl;	
										cout << "new pt1 up: " << newpt1_up << endl;
										cout << "new pt1 down: " << newpt1_down << endl;
									}
									*/

									//fill corresponding histogram
									//is this done within the ORIGINAL trg bin? yeah, right?
									//so added as innermost loop when selecting?
									//For this certain triggerbin, for this certain ybysbin, do the VARIATION?
									if((newpt1_up>30)&&(newpt2_up>30)){	//check if conditions still fulfilled
										cur_ybysbin.unc_histos_up[uncnames.at(k)]->Fill(newptavg_up, Generator_weight);
									}
									if((newpt1_down>30)&&(newpt2_down>30)){
										cur_ybysbin.unc_histos_down[uncnames.at(k)]->Fill(newptavg_down, Generator_weight);
									}
								}//unc sources loop

								//now fill as well Jet_pt_raw, Jet_pt_nom from uncfile (additional information)

							}//if HERWIG


							if(matched){//if both leading jets have gen-match
								///cur_ybysbin.get_mc_genreco_noWeights_noGenCuts()->Fill(gen_ptavg, cur_ptavg);	//fill without any cuts on gen jets

								//fill AND APPLY THE WEIGHTS!
								cur_ybysbin.get_mc_genreco_noGenCuts()->Fill(gen_ptavg, cur_ptavg, gen_weight);	//fill without any cuts on gen jets, but applied WEIGHTS

								//fill the TProfile
								//cur_ybysbin.get_mc_profile()->Fill(gen_ptavg, cur_ptavg/gen_ptavg, gen_weight);	//what does the weight do here?
							

								//now: same selection cuts on gen jets as on the reco jets:
								if(good_genrapidity&&good_genpt){//if gen jets have good rapidity (eta<2.4) and good pt (pt>30GeV)
									if((yb_lo<=gen_yboost && gen_yboost<yb_up) and (ys_lo<=gen_ystar && gen_ystar<ys_up)){
										///cur_ybysbin.get_mc_genreco_noWeights()->Fill(gen_ptavg, cur_ptavg);	//selection on gen and reco, but no weights
	
										cur_ybysbin.get_mc_genreco()->Fill(gen_ptavg, cur_ptavg, gen_weight);	//also applying weights

										if(TMath::Abs((cur_ptavg/gen_ptavg)-1)<0.5){
											cur_ybysbin.get_mc_profile()->Fill(gen_ptavg, cur_ptavg/gen_ptavg, gen_weight);	//SAME SELECTION ON GEN AS ON RECO
											//should the following be done somewhere else instead? (i.e. other if condition?)
											cur_ybysbin.get_mc_profile_mergedj1j2()->Fill(gen_ptavg, (jet1_vec.pt())/(gen_ptavg), gen_weight); //filling for jet1
											cur_ybysbin.get_mc_profile_mergedj1j2()->Fill(gen_ptavg, (jet2_vec.pt())/(gen_ptavg), gen_weight); //filling for jet2
										}
										if(TMath::Abs((jet1_vec.pt())/(gen1_vec.pt())-1)<0.5){
											cur_ybysbin.get_mc_profile_j1()->Fill(gen1_vec.pt(), (jet1_vec.pt())/(gen1_vec.pt()), gen_weight);	//leading gen and reco jets
										}
										if(TMath::Abs((jet2_vec.pt())/(gen2_vec.pt())-1)<0.5){
											cur_ybysbin.get_mc_profile_j2()->Fill(gen2_vec.pt(), (jet2_vec.pt())/(gen2_vec.pt()), gen_weight);	//second leading jets
										}
									}//end: if gen ybys in current ybys bin
								}//end: if good gen-rapidity && good gen-pt
					
							}//endif matched
							count_selected++;	//selected event
						}//end: if MC
						break; //break out of ybys-bin-iteration, because event has now been assigned to its corresponding ybys bin (don't need to check further ybys bins, it is unique)
					}//if dijet-event in ybys bin
					else{	//if no ybys bin was suitable
						if(ybysnr==5){
							if(v==true){cout << "last ybys bin ---> no suitable ybys bin found" << endl;}
							count_NoBin_ybys++;	//no suitable ybys bin found
							break;
						}
						continue;
					}//no suitable ybys bin
				}//loop through ybys-bins
			}//end: if good_rapidity && good_pt
			else{
				count_BadRapidity++;	//this dijet-event does not fulfill rapidity requirements
			}
		}//if (nJet>1)
		else{//if nJet<2
			count_nodijet++;	//this is not a dijet event
		}

		
		if(montecarlo){
			//#################################################################//
			//Loop for the GenJet spectrum -- identical to reco procedure
			//check if Event (jentry) is a dijet event. --> otherwise do nothing
			if(nGenJet>1){	//required for dijet ! --> only look at two leading jets
				//use Lorentzvector:
				PtEtaPhiMVector genjet1_vec(GenJet_pt[0], GenJet_eta[0], GenJet_phi[0], GenJet_mass[0]);
				PtEtaPhiMVector genjet2_vec(GenJet_pt[1], GenJet_eta[1], GenJet_phi[1], GenJet_mass[1]);

				//calculate current ptavg, yb, ys
				Double_t genj1_y = genjet1_vec.Rapidity();
				Double_t genj2_y = genjet2_vec.Rapidity();
				Double_t gen_ptavg = 0.5*(GenJet_pt[0]+GenJet_pt[1]);
				Double_t gen_yboost = 0.5*TMath::Abs(genj1_y+genj2_y);
				Double_t gen_ystar = 0.5*TMath::Abs(genj1_y-genj2_y);

				//add offline cut on GenJet-Rapidity --> both jets must have Abs(y)<5.0 and at least one shall have Abs(y)<2.5
				Double_t absYgenj1 = TMath::Abs(genj1_y);
				Double_t absYgenj2 = TMath::Abs(genj2_y);
				//bool good_rapidity = (absYj1<5)&&(absYj2<5)&&(absYj1<2.5 || absYj2<2.5);	//old cuts
				bool good_rapidity_gen = (absYgenj1<2.4)&&(absYgenj2<2.4);	//NEW cuts!
				//add offline cut on GenJet-pT --> both jets must fulfill pT > 30 GeV.
				bool good_pt_gen = (GenJet_pt[0]>30)&&(GenJet_pt[1]>30);

				//check if rapidity requirement is fulfilled
				//check at the same time whether pT requirement is fulfilled
				//generator weight
				if(good_rapidity_gen&&good_pt_gen){
					//check ybys bin --> put this into function!!
					for (int ybysnr=0; ybysnr!=6; ++ybysnr){
						YbYs cur_ybysbin = ybys_pointer[ybysnr];
						float yb_lo = cur_ybysbin.yb_get_bb().first;
						float yb_up = cur_ybysbin.yb_get_bb().second;
						float ys_lo = cur_ybysbin.ys_get_bb().first;
						float ys_up = cur_ybysbin.ys_get_bb().second;
						if((yb_lo<=gen_yboost && gen_yboost<yb_up) and (ys_lo<=gen_ystar && gen_ystar<ys_up)){
							cur_ybysbin.get_mc_genhisto()->Fill(gen_ptavg, Generator_weight);	//MC-->generator weight
							break; //break out of ybys-bin-iteration, because event has now been assigned to its corresponding ybys bin (don't need to check further ybys bins, it is unique)
						}//if dijet-event in ybys bin
						else{	//if no ybys bin was suitable
							if(ybysnr==5){
								break;
							}
							continue;
						}//no suitable ybys bin
					}//loop through ybys-bins
				}//end: if good_rapidity && good_pt
				else{
				}
			}//if (nGenJet>1)
			else{//if nGenJet<2
			}
			//#################################################################//
		}//end: if montecarlo (for gen spectrum)


	}//end of event-loop




	//finalising the response matrix histogram, in case of MC	(for unfolding)
	if(montecarlo){
		//normalise the filled mc_gen_reco_ptav histogram to unity	(BUT: store both versions!)
		//loop through ybys bins
		for (int ybysnr=0; ybysnr!=6; ++ybysnr){
			YbYs cur_ybysbin = ybys_pointer[ybysnr];

			//clone the original response matrix to create a normalised one
			TH2D* normresp_genreco = (TH2D*)cur_ybysbin.get_mc_genreco()->Clone("normalised_response");		//clone
			normresp_genreco->SetDirectory(cur_ybysbin.get_mc_genreco()->GetDirectory());					//set directory
			normresp_genreco->GetZaxis()->SetTitle("fraction");												//new Z-axis title (changed from Events to fraction)

			//same for the response matrix without genjet-cuts
			TH2D* normresp_genreco_noGenCuts = (TH2D*)cur_ybysbin.get_mc_genreco_noGenCuts()->Clone("normalised_response_noGenCuts");
			normresp_genreco_noGenCuts->SetDirectory(cur_ybysbin.get_mc_genreco_noGenCuts()->GetDirectory());
			normresp_genreco_noGenCuts->GetZaxis()->SetTitle("fraction");


			//genbins are on the x-axis, recobins on the y-axis (these are the same for both matrices)
			Int_t nGenBins = cur_ybysbin.get_mc_genreco()->GetNbinsX();
			Int_t nRecoBins = cur_ybysbin.get_mc_genreco()->GetNbinsY();

			//loop through all the genbins (ROOT-histograms start at binindex 1 !)
			for(int genind=1; genind!=nGenBins+1; ++genind){
				Double_t genBinCont = 0;			//count the bin content in current genbin --> use TH2D::Integral()
				Double_t genBinCont_noGenCuts = 0;	//same for the histo without gen cuts

				//TH2D::Integral() returns sum of bin contents in that range
				//loop through all the recobins and count
				for(int recoind=1; recoind!=nRecoBins+1; ++recoind){
					//genBinCont+=cur_ybysbin.get_mc_genreco()->GetBinContent(genind, recoind);
					if(recoind==1){	//calculate the integral only in first iteration (afterwards it changes as bincontents change!)
						genBinCont = cur_ybysbin.get_mc_genreco()->Integral(genind, genind, 1, nRecoBins);	//calculate integral
						genBinCont_noGenCuts = cur_ybysbin.get_mc_genreco_noGenCuts()->Integral(genind, genind, 1, nRecoBins);	//calculate integral for "noGenCuts" version of matrix
						//if(isnan(genBinCont)==true){break;}
						///if(genBinCont==0){break;} //<-- not used here anymore, as now we check TWO histos at the same time, cannot break, as 2nd histo might have content
						if((genBinCont==0)&&(genBinCont_noGenCuts==0)){break;} // break if genbin of both histos is empty (--> unnecessary to loop through all the remaining reco bins)
					}
					if(genBinCont!=0){//only do this if this genbin is not empty
						Double_t cur_binc = cur_ybysbin.get_mc_genreco()->GetBinContent(genind, recoind);		//get current bin content
						//cur_ybysbin.get_mc_genreco()->SetBinContent(genind, recoind, cur_binc/genBinCont);	//normalise bin content
						normresp_genreco->SetBinContent(genind, recoind, cur_binc/genBinCont);					//normalise bin content
					}
					if(genBinCont_noGenCuts!=0){//only do this if this genbin is not empty
						Double_t cur_binc = cur_ybysbin.get_mc_genreco_noGenCuts()->GetBinContent(genind, recoind);		//get current bin content
						normresp_genreco_noGenCuts->SetBinContent(genind, recoind, cur_binc/genBinCont_noGenCuts);		//normalise bin content
					}
				}
			}
		}//end ybys loop
	}//end: if montecarlo (for mc_gen_reco_ptav histogram normalisation)



	interim_output->Write();

	//Quick summary:
	cout << "Summary after filling the histograms: " << endl;
	cout << "--------------------------------------" << endl;
	cout << "Number of events: " << nentries << endl;
	cout << "Dropped out as nJet<2: " << count_nodijet << endl;
	cout << "Dropped out due to j1, j2 bad rapidity or bad pt: " << count_BadRapidity << endl;
	cout << "Dropped out as no ybys bin suitable: " << count_NoBin_ybys << endl;
	cout << "-----------------------------------------------" << endl;
	cout << "Selected Events: " << count_selected << endl << endl;

	cout << "Events where leading jet has no matching gen jet: " << count_noGen1 << endl;
	cout << "Events where leading jet matched, but no genjet for 2nd jet: " << count_noGen2 << endl << endl;
	cout << "Events where genjet has been matched for both leading recojets: " << count_jetmatch << endl << endl; 

	//adjust this to mc / data case:
	cout << "Saved results to file: " << interim_output->GetName() << endl;
	cout << "Check: count_pthat_fail_100_400 = " << count_pthat_fail_100_400 << endl;

}



//--------------------------------------------------------------------------------------------------//
//											 CREATE JSON FILE										//
//--------------------------------------------------------------------------------------------------//
//new attempt to create function which helps with finding the lumisections and runnumbers for each event
//originally planned to implement this in seperate file --> see get_lumi_txt.C
void DijetAnalysis::get_lumi_txt(){
// function for creating a JSON file via reading "run" and "LuminosityBlock" from nanoAOD file.
// for later use with brilcalc
	TH1::SetDefaultSumw2(kTRUE);		//in order to ALWAYS store the bin errors

	cout << "Creating JSON file for later usage with brilcalc (--> lumi etc)." << endl;

	setTDRStyle();
	fChain->SetBranchStatus("*",0);	//switch off all the branches
	fChain->SetBranchStatus("run",1);
	fChain->SetBranchStatus("luminosityBlock",1);

	//first tests - relict
	/////TFile* json_output = new TFile("run_lumisec.txt", "RECREATE");
	/////json_output->cd();


	ofstream txt_file;
	txt_file.open("txt_file_test.txt");

	//check runperiod and create json file (this function is only called for data anyway --> no if(data) needed)
	//check which runperiod
	ofstream json_file;
	if((strcmp(runperiod,"A")==0)or(strcmp(runperiod,"B")==0)or(strcmp(runperiod,"C")==0)or(strcmp(runperiod,"D")==0)){
		cout << "Run " << runperiod  << endl;
		json_file.open(Form("json_file_Run%s.txt", runperiod));
	}
	else{//if no runperiod specified
		json_file.open("json_file.txt");
	}

	//ofstream json_file;
	//json_file.open("json_file.txt");

	//for testing
	ofstream test_json;
	test_json.open("test_json.txt");


	//create map for run number and lumisections
	//key = run numbers, value = set of pairs (containing lumisections in this run# and corresponding event#s)
	///map<int, set<pair<int, int>>> lumis;

	//map<int, set<pair<int, vector<int>>>> lumis;

	//here without jentry... just mapping of run# and lumisec#s:
	map<int, set<int>> lumis;


	if (fChain == 0) return;
	Long64_t nentries = fChain->GetEntriesFast();
	cout << "nentries = " << nentries << endl << endl;			//information

	Long64_t nbytes = 0, nb = 0;
	for (Long64_t jentry=0; jentry<nentries;jentry++) {
	///for (Long64_t jentry=0; jentry!=160; ++jentry){} //for testing with less entries/events

		Long64_t ientry = LoadTree(jentry);
		if(v==true){cout << "current entry: " << jentry << endl;}
		///if (jentry%10==0){cout << "------------------------------------------" << endl;}//for overview
		if (jentry%100000==0){cout << "processed " << jentry << " of " << nentries << endl;}
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);   nbytes += nb;
		// if (Cut(ientry) < 0) continue;

		//set<int,int> cur_set = make_pair(luminosityBlock, jentry); //create pair of event index and corresponding luminosityBlock
		//lumis.insert(make_pair(run, cur_set));
		////lumis.insert(make_pair(run, make_pair(luminosityBlock, jentry)));
		int cur_lumisec = luminosityBlock;
		auto cur_pair = make_pair(cur_lumisec, jentry);
		int cur_run = run;
		if(v==true){
			cout << "luminosityBlock: " << luminosityBlock << endl;
			cout << "run: " << run << endl;
		}
		/////lumis.insert(pair<int, set<int, int>>(cur_run, cur_pair));

		//insert lumisec# in set which corresponds to key "run":
		///lumis[run].insert(make_pair(luminosityBlock, jentry));

		//without jentry:
		lumis[run].insert(luminosityBlock);
		//cout << "lumis[run]: " << lumis[run].first << endl;
		


   	}//event loop
	//cout << "out of event loop" << endl;

	//now write contents of lumis (map) into .txt file:
	if(v==true){
		cout << "size of created map is: " << lumis.size() << endl;
		cout << "size of inner set is: " << lumis[run].size() << endl;
	}
	///////map<int, set<pair<int, int>>>::const_iterator it = lumis.begin();
	map<int, set<int>>::const_iterator it = lumis.begin();
	//while(it != lumis.end()){}
	json_file << "{";
	for(it; it!=lumis.end(); ++it){
		//cout << "run number: " << it->first << " lumisec: " << it->second << endl;
		//cout << "run" << endl;
		
		//within run-number-iteration, go through set of lumisections:
		//set<pair<int, int>>::iterator itlumi = (lumis->second).begin();
		//set<pair<int, int>> curset = it->second;
		/////set<pair<int, int>>::iterator itlumi = it->second.begin();
		//set<pair<int, int>>::const_iterator itlumi = curset.begin();

		txt_file << " --------------------------" << endl;
		txt_file << Form(" |  this is run#: %d  |", it->first) << endl;
		txt_file << " --------------------------" << endl;

		json_file << Form("\"%d\": [[", it->first);
		test_json << Form("\"%d\": [[", it->first); //testing

		///set<pair<int, int>>::iterator itlumi = lumis.second.begin();
		//////while(itlumi != it->second.end()){}
		//while(itlumi != curset.end()){}
		///////for(itlumi; itlumi!=curset.end(); ++itlumi){}
		set<int>::const_iterator it_lumis = it->second.begin();
		int cur_lumisec = (*it_lumis)-1;
		int old_lumisec = *it_lumis;	//initial value (first lower bound of lumisection range)
		json_file << *it_lumis; //initial value
		for(it_lumis; it_lumis!=(it->second).end(); ++it_lumis){

		///while(itlumi != lumis.second.end()){}
			//cout << "lumisection: " << itlumi->first << " event index: " << itlumi->second << endl;
			/////cout << "it " << Form("%d", itlumi) << endl;
			/////txt_file << "test: itlumi->first:: " << itlumi->first << "\n";
			
			//cur_lumisec is current lumisection --> to check if this is lumi range or single section
			//cout << "current *it_lumis: " << *it_lumis << endl;
			if((cur_lumisec+1)!=*it_lumis){
				//cout << "cur_lumisec+1: " << cur_lumisec+1 << " *it_lumis: " << *it_lumis << endl;
				//if new lumisection starts and this is not a continuous range --> close previous range, open new one
				json_file << Form(", %d], [%d", old_lumisec,  *it_lumis);
				//cur_lumisec = *it_lumis;
			}
			else{
				//range of lumisections continues, don't close brackets yet
				//continue;
			}

			cur_lumisec = *it_lumis;
			old_lumisec = *it_lumis;
			txt_file << Form("(counter: \?) lumisection#: %d", *it_lumis) << endl;
			test_json << Form("[%d,%d],", *it_lumis, *it_lumis); //test!
			//break;
		}//while-loop for lumisection (inner loop) //--> is now a FOR LOOP
		txt_file << endl;

		//last comma and space only if this is not overall last entry
		//reverse iterator (via rbegin()) to the reverse beginning of the container == its last element
		map<int, set<int>>::reverse_iterator rit = lumis.rbegin();
		if(it->first!=rit->first){
			json_file << Form(", %d]], ", cur_lumisec); //test! writing the lastbin bound of last lumisection in current run-number
		}
		else{	//no comma when we reach the overall end
			json_file << Form(", %d]]", cur_lumisec);
		}
		//break;
	}//end writing-while-loop // --> is now a FOR LOOP

	json_file << "}";
	test_json << "}"; //test
	txt_file.close();
	return 0;

	cout << "Saved lumisection results to JSON file: run_lumisec.txt" << endl;
}//end function get_lumi_txt()



//###################################################################################################################
//																													#
//							Current Main Function --> handles original histograms + trigger fits					#
//							--> under construction --> move histogram filling in Loop function!						#
//																													#
//###################################################################################################################

void DijetAnalysis::triggerpaths(){
	TH1::SetDefaultSumw2(kTRUE);		//in order to ALWAYS store the bin errors
	setTDRStyle();
	std::string name_h = "hpt";

	//for now moved to Loop() function
	/*
	fChain->SetBranchStatus("*",0); //exclude / switch off all the branches
	fChain->SetBranchStatus("nJet",1); //include this branch in analysis	
	fChain->SetBranchStatus("Jet_pt",1); //include this branch as well
	fChain->SetBranchStatus("Jet_eta",1); //has to be included if eta-bins shall be checked...
	fChain->SetBranchStatus("Jet_phi",1);	//for TLorentzVector
	fChain->SetBranchStatus("Jet_mass",1);	//for TLorentzVector

	fChain->SetBranchStatus("Generator_weight",1);

	//for additional trigger studies
	//Int_t j1=-1;	//triggerobject index of leading jet (arbitrary initial value)
	//Int_t j2=-1;	//triggerobject index of second leading jet	(arbitrary initival value)	--> is this okay?
	fChain->SetBranchStatus("nTrigObj",1);
	fChain->SetBranchStatus("TrigObj_id",1);
	fChain->SetBranchStatus("TrigObj_pt",1);

	*/


	//AT THE SAME TIME SET LUMI-VECTOR ACCORDING TO RUNPERIOD!
	//new lumi (from own json files, but after skimming, used brilcalc) units are /fb
	vector<Double_t> lumivec;

	//open input file (histo_file)
	const char* inputfile_name;
	//if(montecarlo==false){ inputfile_name = "dijet_interim_output_data.root"; }
	if(montecarlo==false){//if data
		//check which runperiod
		if(strcmp(runperiod,"A")==0){
			cout << "Run A" << endl;
			inputfile_name = "dijet_interim_output_data_RunA.root";
			lumivec = {0.000000084, 0.000000084, 0.000065758, 0.000088949, 0.002453039, 0.010118712, 0.044313793, 0.115597030, 0.397308191, 0.910868208, 1.780016341, 13.977332725, 13.977332725};
		}//endif runA
		else if(strcmp(runperiod,"B")==0){
			cout << "Run B" << endl;
			inputfile_name = "dijet_interim_output_data_RunB.root";
			lumivec = {0, 0, 0.000026808, 0.000110281, 0.000415176, 0.005846909, 0.025058183, 0.055140569, 0.220562275, 0.441124550, 0.882249099, 7.057992794, 7.057992794};
		}//endif runB
		else if(strcmp(runperiod,"C")==0){
			cout << "Run C" << endl;
			inputfile_name = "dijet_interim_output_data_RunC.root";
			lumivec = {0.000000063, 0.000000063, 0.000027608, 0.000107731, 0.000405575, 0.005647855, 0.024204982, 0.054976951, 0.216547076, 0.431973909, 0.862827576, 6.894778911, 6.894778911};
		}//endif runC
		else if(strcmp(runperiod,"D")==0){
			cout << "Run D" << endl;
			inputfile_name = "dijet_interim_output_data_RunD.root";
			lumivec = {0.000000041, 0.000000041, 0.000119557, 0.000495470, 0.001865299, 0.026760713, 0.114682242, 0.248467482, 0.991655052, 1.982571812, 3.964405332, 31.710074612, 31.710074612};
		}//endif runD
		else{//if no runperiod specified
			inputfile_name = "dijet_interim_output_data.root";
		}
	}
	else{ inputfile_name = "dijet_interim_output_mc.root"; } //if MC

	TFile *inputfile = TFile::Open(inputfile_name);

	//TFile* trglumi_output = new TFile("trg_lumi_output.root", "RECREATE"); //for storing the histograms
	TFile* trglumi_output;
	if(montecarlo==false){//if data
		//trglumi_output = new TFile("trg_lumi_output_data.root", "RECREATE"); //output file after analysis BEFORE creating plots
		trglumi_output = new TFile(Form("trg_lumi_output_data_Run%s.root", runperiod), "RECREATE"); //TEST!!
	}
	else{ //if MC
		trglumi_output = new TFile("trg_lumi_output_mc.root", "RECREATE"); //output file after analysis BEFORE creating plots
	}
	trglumi_output->cd();
	trglumi_output->mkdir("Standard"); //create directory in output file

	//create subfolder "Standard" in outputfile
	TDirectory *standard_folder = trglumi_output->GetDirectory("Standard"); 
	gROOT->GetListOfBrowsables()->Add(standard_folder, "Standard");
	standard_folder->cd();



	vector<int> triggers = {0, 15, 25, 40, 60, 80, 140, 200, 260, 320, 400, 450, 500, 550};
	///vector<Double_t> lumivec = {0, 0, 0.000414643, 0.001737711, 0.006541972, 0.097552471, 0.418082018, 0.868855667, 3.475422669, 6.950845338, 13.901690676, 111.213525407, 111.213525407};	//now for testing --> should be more automatised / generalised [units are /pb]




	/*
	if(montecarlo==false){
		//switching on selected TBranches corresponding to chosen triggers:
		for (int j=1; j!=triggers.size(); ++j){ //start from 1, because "0" means "without trigger"
			char *trigname = Form("HLT_PFJet%i", triggers.at(j));
			fChain->SetBranchStatus(trigname, 1);
		}
	}
	*/

	YbYs *ybys_pointer;	//should maybe handle this differently
	//call constructor YbYs(float yb_lo, float ys_lo, float yb_up, float ys_up)
	///ybys_pointer = new YbYs[6]{{0.0,0.0,1.0,1.0}, {0.0, 1.0, 1.0, 2.0}, {0.0, 2.0, 1.0, 3.0}, {1.0, 0.0, 2.0, 1.0}, {1.0, 1.0, 2.0, 2.0}, {2.0, 0.0, 3.0, 1.0}};	

	//new outermost bin bounds --> go only up to 2.4 for both yb, ys
	ybys_pointer = new YbYs[6]{{0.0,0.0,1.0,1.0}, {0.0, 1.0, 1.0, 2.0}, {0.0, 2.0, 1.0, 2.4}, {1.0, 0.0, 2.0, 1.0}, {1.0, 1.0, 2.0, 2.0}, {2.0, 0.0, 2.4, 1.0}}; //new!!!





	vector <Histos*> histosvec;
	for (int ybysnr=0; ybysnr!=6; ++ybysnr){
		YbYs cur_ybysbin = ybys_pointer[ybysnr];
		if(montecarlo==false){//only if data==true --> montecarlo==false
			histosvec.push_back(new Histos(triggers, inputfile, cur_ybysbin.get_bname().c_str()));
			Histos curtrig = *histosvec.back();
			YbYs cur_ybysbin = ybys_pointer[ybysnr];
			cur_ybysbin.ybys_directories(standard_folder, triggers, curtrig, true);	//creates folder structure	//"true" stands for "trgmode=true" --> writes histograms
		}
		else{//if MC
			cur_ybysbin.ybys_directories(standard_folder);
		}//if MC
	}

	//counters
	int count_nodijet = 0; 		//count number of events that are no dijet events, i.e. nJet < 2
	int count_NoBin_ybys = 0; 	//count number of events that have no ybys bin, where they could be assigned to (i.e. yb or ys > 3.0)
	int count_BadRapidity = 0;	//count number of events where the jets do not fulfill the rapidity requirements
	int count_selected = 0;		//count number of selected events: implemented for MC


	/* // HANDLED LATER WHEN NEEDED!!!!! 	
	//map<const char*, TH1D*> trglumi_hist_map;		//not needed here, as trglumi_hist_map is now a member of class YbYs
	for (int ybysnr=0; ybysnr!=6; ++ybysnr){
		YbYs cur_ybysbin = ybys_pointer[ybysnr];
		cur_ybysbin.trglumi_directories(standard_folder, triggers, cur_ybysbin.getmap());	//creates folder structure
		//cur_ybysbin.trglumi_directories(standard_folder, triggers, trglumi_hist_map);	//creates folder structure
	}
	*/

	//should write this in an extra general function, as this is used in almost every DijetAnalysis:: function() here.
	/*
	if (fChain == 0) return;
	Long64_t nentries = fChain->GetEntriesFast();
	cout << "nentries = " << nentries << endl << endl;			//information

	Long64_t nbytes = 0, nb = 0;
	//starting event loop
	for (Long64_t jentry=0; jentry<nentries;jentry++) {
	///for (Long64_t jentry=0; jentry!=160; ++jentry){} //for testing with less entries/events

		Long64_t ientry = LoadTree(jentry);
		if(v==true){cout << "current entry: " << jentry << endl;}
		///if (jentry%10==0){cout << "------------------------------------------" << endl;}//for overview
		if (jentry%100000==0){cout << "processed " << jentry << " of " << nentries << endl;}
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);   nbytes += nb;
		
		map<const int, Bool_t> tdecisions;

		if(montecarlo==false){
			tdecisions.insert(make_pair(triggers.at(0), true)); //-->no trigger
			tdecisions.insert(make_pair(triggers.at(1), HLT_PFJet15)); //could this be handled in a more general way? (see attempt and problem above)
			tdecisions.insert(make_pair(triggers.at(2), HLT_PFJet25));
			tdecisions.insert(make_pair(triggers.at(3), HLT_PFJet40));
			tdecisions.insert(make_pair(triggers.at(4), HLT_PFJet60));
			tdecisions.insert(make_pair(triggers.at(5), HLT_PFJet80));
			tdecisions.insert(make_pair(triggers.at(6), HLT_PFJet140));
			tdecisions.insert(make_pair(triggers.at(7), HLT_PFJet200));
			tdecisions.insert(make_pair(triggers.at(8), HLT_PFJet260));
			tdecisions.insert(make_pair(triggers.at(9), HLT_PFJet320));
			tdecisions.insert(make_pair(triggers.at(10), HLT_PFJet400));
			tdecisions.insert(make_pair(triggers.at(11), HLT_PFJet450));
			tdecisions.insert(make_pair(triggers.at(12), HLT_PFJet500));
			tdecisions.insert(make_pair(triggers.at(13), HLT_PFJet550));
		}

		Int_t j1=-1;	//triggerobject index of leading jet (arbitrary initial value)
		Int_t j2=-1;	//triggerobject index of second leading jet	(arbitrary initival value)	--> is this okay?


		//check if Event (jentry) is a dijet event. --> otherwise do nothing
		if(nJet>1){	//required for dijet ! --> only look at two leading jets
			Double_t j1_pt	= Jet_pt[0];
			Double_t j2_pt	= Jet_pt[1];
			Double_t j1_eta	= Jet_eta[0];
			Double_t j2_eta	= Jet_eta[1];
			Double_t j1_phi	= Jet_phi[0];
			Double_t j2_phi	= Jet_phi[1];
			Double_t j1_mass = Jet_mass[0];
			Double_t j2_mass = Jet_mass[1];
			//use Lorentzvector:
			PtEtaPhiMVector jet1_vec(j1_pt, j1_eta, j1_phi, j1_mass);
			PtEtaPhiMVector jet2_vec(j2_pt, j2_eta, j2_phi, j2_mass);

			//calculate current ptavg, yb, ys
			Double_t j1_y = jet1_vec.Rapidity();
			Double_t j2_y = jet2_vec.Rapidity();
			Double_t cur_ptavg = 0.5*(Jet_pt[0]+Jet_pt[1]);
			Double_t cur_yboost = 0.5*TMath::Abs(j1_y+j2_y);
			Double_t cur_ystar = 0.5*TMath::Abs(j1_y-j2_y);

			//add offline cut on Jet-Rapidity --> both jets must have Abs(y)<5.0 and at least one shall have Abs(y)<2.5
			Double_t absYj1 = TMath::Abs(j1_y);
			Double_t absYj2 = TMath::Abs(j2_y);
			bool good_rapidity = (absYj1<5)&&(absYj2<5)&&(absYj1<2.5 || absYj2<2.5);

			//add offline cut on Jet-pT --> both jets must fulfill pT > 30 GeV.
			bool good_pt = (j1_pt>30)&&(j2_pt>30);

			//for additional trigger studies
			//loop over all the trigger objects in current event
			//find the two hardest trigger objects in this event and store their trigger object index (as j1, j2)
			//cout << "nTrigObj: " << nTrigObj << endl;
			//use pt
			Double_t pt1=0;
			Double_t pt2=0;
			if(montecarlo==false){
				for(int obj_ind=0; obj_ind!=nTrigObj; ++obj_ind){
					//cout << "in loop. obj_ind=" << obj_ind << endl;
					//check that the ID corresponds to Jet (TrigObj_id == 1)
					if(TrigObj_id[obj_ind]==1){
						//if(j1==-1){j1=obj_ind;}	//setting initial j1
						//how shall j2 be set initially?	--> arbitrarily ... will change anyway in case it's wrongly assigned.
						Double_t cur_objpt = TrigObj_pt[obj_ind];	//pt of current trigger object
						//if(cur_objpt>TrigObj_pt[j1]){
						//if(j1==-1 || cur_objpt>TrigObj_pt[j1]){
						if(cur_objpt>pt1){
							j2 = j1;
							j1 = obj_ind;	//current trigger object harder than previously found hardest one
							pt2 = pt1;
							pt1 = cur_objpt;
							
							//cout << "changed j1 to current obj_ind." << endl;
						}
						//else if((cur_objpt<TrigObj_pt[j1])&&(cur_objpt>TrigObj_pt[j2])){
						//else if(cur_objpt>TrigObj_pt[j2]){
						else if(cur_objpt>pt2){
							pt2 = cur_objpt;
							j2 = obj_ind;
							//cout << "changed j2 to current obj_ind." << endl;
						}
						//else{
						//	continue;	//current trigger object had smaller pt than the ones investigated before
						//}
					}
				
					else{	//this was not a jet (does not have correct ID)
						//cout << "not a jet." << endl;
						//continue;
					}
				}
			}//if MC==false --> IF DATA

			//check if rapidity requirement is fulfilled
			//check at the same time whether pT requirement is fulfilled
			//generator weight
			if(good_rapidity&&good_pt){
				Double_t gen_weight = Generator_weight;		//MC event weight to be given to fillhistos() function
				//check ybys bin --> put this into function!!
				for (int ybysnr=0; ybysnr!=6; ++ybysnr){
					YbYs cur_ybysbin = ybys_pointer[ybysnr];
					float yb_lo = cur_ybysbin.yb_get_bb().first;
					float yb_up = cur_ybysbin.yb_get_bb().second;
					float ys_lo = cur_ybysbin.ys_get_bb().first;
					float ys_up = cur_ybysbin.ys_get_bb().second;
					if(v==true){cout << "yb_lo: " << yb_lo << " ## yb_up: " << yb_up << " ## ys_lo: " << ys_lo << " ## ys_up: " << ys_up << endl;}
					if((yb_lo<=cur_yboost && cur_yboost<yb_up) and (ys_lo<=cur_ystar && cur_ystar<ys_up)){	//-->better to implement function InInterval(value, low, up) that checks if value in interval... return(lo<=x & x<up)
						if(v==true){cout << "Event is in ybys bin: " << cur_ybysbin.get_bname() << endl;}

						if(montecarlo==false){//only for data:
							//check trigger
							for(int trignr=0; trignr!=triggers.size(); ++trignr){
								if(tdecisions[triggers.at(trignr)]){	//==true
									//fill into histo of current trigger
									histosvec.at(ybysnr)->jt_histos[triggers.at(trignr)]->Fill(cur_ptavg);	//fill ptavg trigger histogram

									//check if at least one of the leading jets' pt is higher than next trigger-threshold
									if((trignr!=triggers.size()-1) && (j1 != -1) && (TrigObj_pt[j1] >= triggers.at(trignr+1))){
									//if((TrigObj_pt[j1] >= triggers.at(trignr+1)) || (TrigObj_pt[j2] >= triggers.at(trignr+1))){
										//fill "simulated" histo of next higher trigger by applying additional condition
										histosvec.at(ybysnr)->trgobj_histos[triggers.at(trignr+1)]->Fill(cur_ptavg);
									}
									//check if at least one of the leading jets' pt is higher than the next next trigger-threshold (N+2)
									if((trignr<triggers.size()-2) && (j1 != -1) && (TrigObj_pt[j1] >= triggers.at(trignr+2))){
										histosvec.at(ybysnr)->n2_trgobj_histos[triggers.at(trignr+2)]->Fill(cur_ptavg);
									}
									continue;
								}//trigger true?
								else{//trigger false
									continue;
								}//trigger false
							}//check triggers
						}//end if DATA (if mc==false)
						else if(montecarlo==true){
							cur_ybysbin.get_mc_histo()->Fill(cur_ptavg, gen_weight);//apply generator weight (as this is MC)
							count_selected++;	//selected event
						}
						break; //break out of ybys-bin-iteration, because event has now been assigned to its corresponding ybys bin (don't need to check further ybys bins, it is unique)
					}//if dijet-event in ybys bin
					else{	//if no ybys bin was suitable
						if(ybysnr==5){
							if(v==true){cout << "last ybys bin ---> no suitable ybys bin found" << endl;}
							count_NoBin_ybys++;	//no suitable ybys bin found
							break;
						}
						continue;
					}//no suitable ybys bin
				}//loop through ybys-bins
			}//end: if good_rapidity && good_pt
			else{
				count_BadRapidity++;	//this dijet-event does not fulfill rapidity requirements
			}
		}//if (nJet>1)
		else{//if nJet<2
			count_nodijet++;	//this is not a dijet event
		}
	}//end of event-loop
	trglumi_output->Write();

	//Quick summary:
	cout << "Summary after filling the histograms: " << endl;
	cout << "--------------------------------------" << endl;
	cout << "Number of events: " << nentries << endl;
	cout << "Dropped out as nJet<2: " << count_nodijet << endl;
	cout << "Dropped out due to j1, j2 bad rapidity: " << count_BadRapidity << endl;
	cout << "Dropped out as no ybys bin suitable: " << count_NoBin_ybys << endl;
	cout << "-----------------------------------------------" << endl;
	cout << "Selected Events: " << count_selected << endl;

	*/

	//FOLLOWING TRIGGER FITS ONLY FOR DATA --> if montecarlo==false
	
	if(montecarlo==false){
		//define some vectors with values used in the fit
		//vector<int> triggers = {0, 15, 25, 40, 60, 80, 140, 200, 260, 320, 400, 450, 500, 550};
		//first trigger considered is 25 GeV (divided by 15GeV trigger) [due to luminosity = 0, trg25/trg15 and trg40/trg25 don't help]
		//vector<int> xlow_vec = {10, 26, 30, 40, 60, 100, 140, 200, 300, 360, 400, 460};
		//vector<int> xup_vec = {60, 100, 180, 210, 220, 360, 400, 550, 900, 920, 960, 1200};

		//now only start when lower trigger in ratio 90% efficient (guess)
		vector<int> xlow_vec = {10, 26, 30, 40, 60, 100, 140, 200, 300, 360, 400, 460};
		vector<int> xup_vec = {60, 100, 180, 210, 220, 360, 400, 550, 900, 920, 960, 1200};


		//for test where only yb0ys0 trigger thresholds are used (for all the ybys bins)
		//this is an additional map, should instead improve triggerfits and use the original "trg_effmap" (individual for each ybys bin) for the sumhisto
		map<int, Double_t> trg_effmap_yb0ys0;	//mapping trigger threshold (in pTavg) onto nominal trigger value (int)

		//Go again through different ybys bins and check the trigger ratios
		for (int ybysnr=0; ybysnr!=6; ++ybysnr){
			YbYs cur_ybysbin = ybys_pointer[ybysnr];
			float yb_lo = cur_ybysbin.yb_get_bb().first;
			float yb_up = cur_ybysbin.yb_get_bb().second;
			float ys_lo = cur_ybysbin.ys_get_bb().first;
			float ys_up = cur_ybysbin.ys_get_bb().second;

			
			//handle directories in output file --> write ratio plots in corresponding ybys folder
			//create new directory for the trigger fits
			const char *dirname = Form("%s", cur_ybysbin.get_bname().c_str());
			const char *ratio_dirname = Form("trigger_ratios");
			standard_folder->cd(dirname);
			//standard_folder->mkdir(dirname);
			TDirectory *ybysDir = standard_folder->GetDirectory(dirname);
			//ybysDir->cd();
			ybysDir->mkdir(ratio_dirname);
			TDirectory *ratioDir = ybysDir->GetDirectory(ratio_dirname);



			//map<const char*, TH1D*> cur_ratiohist_map = cur_ybysbin.getmap();	//remove this... too complicated and unnecessary


			//later: create graphs that are used to check which parameter is more stable in the fit
			//here: create the vectors that will be used to store x- and y- coordinates of the graphs
			//pnom = parameter which is approx. nominal value of trigger threshold
			//pwid = parameter which is approx. width of errorfunction (turn-on-curve) (=sqrt(pnom))
			vector<Double_t> curx_pnom;
			vector<Double_t> cury_pnom;
			vector<Double_t> curx_pwid;
			vector<Double_t> cury_pwid;

			vector<Double_t> pnom_err;
			vector<Double_t> pwid_err;

			//same procedure as well for multigraph (mg_cut)
			vector<Double_t> mg_pnom_X;
			vector<Double_t> mg_pnom_Y;
			vector<Double_t> mg_pwid_X;
			vector<Double_t> mg_pwid_Y;

			vector<Double_t> mg_pnom_err;
			vector<Double_t> mg_pwid_err;

			//and once more for the 1param fit in the end
			vector<Double_t> pnom_1p_X;
			vector<Double_t> pnom_1p_Y;
			vector<Double_t> pwid_1p_X;
			vector<Double_t> pwid_1p_Y;

			vector<Double_t> pnom_1p_err;
			vector<Double_t> pwid_1p_err;


			//now divide histograms
			//start from index "2", because histogram with index "0" is jt0 (--> no trigger) and histogram with index "1" will be the first denominator in the ratios (-->divide Nth histo by (N-1)th histo)
		
			//start at trigger-index "4" which is jt60 (ratio to ind=3 --> jt40), because jt15 (ind=1) and jt25 (ind=2) have lumi=0
			//luminosity values start from index "0" nevertheless
			//for(int ind=2; ind!=triggers.size(); ++ind){}
			for(int ind=4; ind!=triggers.size(); ++ind){
				TH1D* prehist = (TH1D*)histosvec.at(ybysnr)->jt_histos[triggers.at(ind-1)]->Clone(Form("yb%01.0fys%01.0f_trg%d", yb_lo, ys_lo, triggers.at(ind-1))); //previous histogram (N-1)
				TH1D* curhist = (TH1D*)histosvec.at(ybysnr)->jt_histos[triggers.at(ind)]->Clone(Form("yb%01.0fys%01.0f_trg%d", yb_lo, ys_lo, triggers.at(ind)));		//current histogram (N), that will be divided by hist(N-1)


				//newly added: also check N/N-2 --> ratio to prevprevtrigger
				TH1D* prevprev_ratio = (TH1D*)histosvec.at(ybysnr)->jt_histos[triggers.at(ind)]->Clone(Form("yb%01.0fys%01.0f_trg%d_to_trg%d", yb_lo, ys_lo, triggers.at(ind), triggers.at(ind-2))); //current histogram (N)
				TH1D* preprehist = (TH1D*)histosvec.at(ybysnr)->jt_histos[triggers.at(ind-2)]->Clone(Form("yb%01.0fys%01.0f_trg%d_asPrePreTrg", yb_lo, ys_lo, triggers.at(ind-2))); //preprevious histogram (N-2)


				/*
				TH1D* curhist = cur_ratiohist_map[Form("jet%d_jet%d", triggers.at(ind), triggers.at(ind-1))];	//curhist, already created earlier, will now be filled with ratio
				int curbins = prehist->GetNbinsX();
				for(int k=0; k!=curbins+1; ++k){
					Double_t currbincont = (histosvec.at(ybysnr)->jt_histos[triggers.at(ind)])->GetBinContent(k);
					curhist->SetBinContent(k, currbincont);
				}
				*/
		
				//setting the ratio plots' directory to ratioDir, which has been created for this current ybysbin (done one level higher in ybys loop)
				curhist->SetDirectory(ratioDir);
				prevprev_ratio->SetDirectory(ratioDir);
				ratioDir->cd();		//change to directory where the trigger-ratios and corresponding fits are stored

				//Taking care of N/N-1
				curhist->Sumw2();			//important for error handling  --> see TH1D::Divide() in ROOT documentation
				prehist->Sumw2();
				curhist->SetLineWidth(2);
				curhist->Scale(1./(lumivec.at(ind-1)));	//scale current histo by lumi
				prehist->Scale(1./(lumivec.at(ind-2)));	//scale previous histogram by lumi
				curhist->Divide(curhist, prehist, 1, 1, "B");				//divide current histo by previous one (both already scaled by lumi), binominal statistics

				//Taking care of N/N-2
				prevprev_ratio->Sumw2();
				preprehist->Sumw2();
				prevprev_ratio->Scale(1./(lumivec.at(ind-1)));	//scale CURRENT triggerhisto by lumi (luminidices are shifted by one compared to trigger, as triggervec includes jt0 --> none)
				preprehist->Scale(1./(lumivec.at(ind-3)));	//scale preprehist by corresponding lumi (first one is ind=4-3=1 --> lumi25)
				prevprev_ratio->Divide(prevprev_ratio, preprehist, 1, 1, "B");		//divide current histo by preprevious one (Binomial statistics)

				
				//adjusting plotting for N/N-1
				SetDrawStyleWrite(curhist, "", "", "p_{T,avg}", Form("#frac{trg%d / lumi%d}{trg%d / lumi%d}", triggers.at(ind), triggers.at(ind), triggers.at(ind-1), triggers.at(ind-1)));

				//finding suitable y-axis or using default?
				//does not work without canvas...
				//take care of such things in EXTRA plotting script
				/*
				TLegend* leg = tdrLeg(0.7, 0.2, 0.9, 0.5);
				leg->SetTextFont(42);
		
				TLatex *bininfo = new TLatex();
				bininfo->SetNDC();
				bininfo->SetTextSize(0.045);
				bininfo->DrawLatex(0.7, 0.82, Form("#font[42]{yb%01.0fys%01.0f}", yb_lo, ys_lo)); //font 12 looks nicer, but does not get saved correctly...

				//Add horizontal line at y=1.0 (canvas needs coordinate system to place line into!)
				//TLine *line = new TLine(divcanv->GetUxmin(), 1.0, divcanv->GetUxmax(), 1.0);
				TLine *line = new TLine(10, 1.0, 6300, 1.0);
				line->SetLineColorAlpha(kAzure+5, 0.36);
				line->SetLineStyle(9);
				line->Draw("SAME");
				//line->Draw();
				*/

				//adjusting plotting for N/N-2
				//just basic plotting, later combined, nicer plot... this is just testing
				//for now --> no canvas --> add it to the N/N-1 plot canvas later, or create new, shared one
				SetDrawStyleWrite(prevprev_ratio, "", "", "p_{T,avg}", Form("#frac{trg%d / lumi%d}{trg%d / lumi%d}", triggers.at(ind), triggers.at(ind), triggers.at(ind-2), triggers.at(ind-2)));


				//############################### FITTING ##############################//
				//Now: Fit erf() onto curve
				//set range, where to fit (so far just a guess..)
				Double_t rangelow = triggers.at(ind-1);
				Double_t rangeup = 1.5*triggers.at(ind);

				// Approx. the values of the fitting parameters: 
				// 	taken as initial fit parameter values
				//	[0]-->efficiency at saturation, here assumed to be "1.0"
				//	[1]-->position of turn-on (50% eff. // nominal value), 
				//	[2]-->(1./width)=Sqrt([1]) or 5% at high pT

				// 3-param-fit
				// try this as well! (so far only 2-param-fit done)
				TF1 *erf = new TF1("erf", "[0]*0.5*(1+TMath::Erf((x-[1])/[2]))", rangelow, rangeup);
				erf->SetParameter(0, 1.0); //--> set parameter [0] to 1.
				erf->SetParameters(1.0, triggers.at(ind), TMath::Sqrt(triggers.at(ind))); //set initial p1 and p2 (see above) [thinking of better initial values?]

				// iterate --> start fit of next fitting round at X, where previous Trigger 90% eff. (found in previous fitting round)

				// 2-param-fit (set [0] to 1.0, renaming of [1]-->[0] and [2]-->[1])
				TF1 *erf2p_test = InitFit2P("erf2p_test", "0.5*(1+TMath::Erf((x-[0])/[1]))", rangelow, rangeup, kRed, triggers.at(ind), TMath::Sqrt(triggers.at(ind)), "pnom", "pwid");

				//curhist->Fit("erf","R"); //-->"R" = take range from constructor
				//curhist->Fit("erf", "", "", xlow_vec.at(ind-2), xup_vec.at(ind-2));	//3param-fit
				////curhist->Fit("erf2p_test", "", "", xlow_vec.at(ind-2), xup_vec.at(ind-2));		//2param-fit



				//######### METHOD A --> fit from 30% of current trigger up to 99.2% of current trigger #################//
				//Double_t cur_xlow_new = xlow_vec.at(ind-2);
				//Double_t cur_xup_new = xup_vec.at(ind-2);
				Double_t cur_xlow_new = rangelow;
				Double_t cur_xup_new = rangeup;

				Double_t prevtrig90 = rangelow;
				for(int it_ind=0; it_ind!=5; ++it_ind){
					curhist->Fit("erf2p_test", "Q", "", cur_xlow_new, cur_xup_new); //"Q" for quiet --> minimal output
					//curhist->Fit("erf2p_test", "", "", prevtrig90, cur_xup_new);

					cur_xlow_new = erf2p_test->GetX(0.30);		//will start fit of next iteration at 30% effic. of trigger (is not the same as 90% eff. of prev. trg.)
					cur_xup_new = erf2p_test->GetX(0.992);		//will end fit of next iteration at 99.2% effic. of current trigger
					//prevtrig90 = erf2p_test->GetX(0.90); 		//DOES NOT work here, because we are WITHIN triggerloop.. need to fit outside of this loop
					
				}

				//Double_t prev90eff = erf2p_test->GetX(0.90); // --> use this in next round of trigger-loop as a strat value!

				//get the parameters from the fit:
				Store2Params(triggers.at(ind), erf2p_test, curx_pnom, cury_pnom, pnom_err, curx_pwid, cury_pwid, pwid_err);

				//write histogram plus its fit
				curhist->Write(Form("%s_Fit",curhist->GetName()));
			}//end trigger-loop


			//###################### METHOD B --> fit from 80% of previous trigger up to 99.92% of current trigger ##############//
			//##### additionally: N/N-2
			//new fitting loop for using information of previous trigger-fit for the next trigger
			//use very big range initially
			Double_t rangelow = 40.00;
			Double_t rangeup = 630.00;

			Double_t N2rangelow = 30.00;
			Double_t N2rangeup = 630.00;
			vector<Double_t> N2rangelow_vec={31.0, 48.0, 70.0};	//needs to start with X-position: 25GeV trigger @ 80% (guessed with fit)
			vector<Double_t> N1rangeup_vec;	//stores upper fitrange-limit obtained from N1-fit (=N1fit)

			TDirectory *curdirFit = standard_folder->GetDirectory(Form("%s/trigger_ratios", cur_ybysbin.get_bname().c_str()));
			standard_folder->cd();
			//start at ind=4 (==>> 60GeV-trigger (over 40GeV --> ind-1))
			for(int ind=4; ind!=triggers.size(); ++ind){
				//initial fitting range for current "trigger-pair"
				//use very big range initially
				//Double_t rangelow = 40.00;
				Double_t rangeup = 630.00;	//rangeup will be set according to current trigger, rangelow has been obtained from previous trigger
				Double_t N2rangeup = 630.00;

				TH1D* N1hist; 
				trglumi_output->GetObject(Form("Standard/%s/trigger_ratios/yb%01.0fys%01.0f_trg%d_Fit", cur_ybysbin.get_bname().c_str(), yb_lo, ys_lo, triggers.at(ind)), N1hist);

				//get the N/N-2 histo
				TH1D* N2hist;
				const char* prevprevname =Form("yb%01.0fys%01.0f_trg%d_to_trg%d", yb_lo, ys_lo, triggers.at(ind), triggers.at(ind-2));
				trglumi_output->GetObject(Form("Standard/%s/trigger_ratios/%s", cur_ybysbin.get_bname().c_str(), prevprevname), N2hist);

				// 2-param-fit for fitting of N/N-1 (set [0] to 1.0, renaming of [1]-->[0] and [2]-->[1])
				TF1 *N1fit = InitFit2P("N1fit", "0.5*(1+TMath::Erf((x-[0])/[1]))", rangelow, rangeup, kGreen+2, triggers.at(ind), TMath::Sqrt(triggers.at(ind)), "pnom", "pwid");
				
				//same 2-param-fit function for fitting of N/N-2
				TF1 *N2fit = InitFit2P("N2fit", "0.5*(1+TMath::Erf((x-[0])/[1]))", N2rangelow, N2rangeup, kBlue+2, triggers.at(ind), TMath::Sqrt(triggers.at(ind)), "pnom", "pwid");

				//fit iteratively
				for(int it_ind=0; it_ind!=5; ++it_ind){
					N1hist->Fit("N1fit", "Q", "", rangelow, rangeup); //"0" for not drawing immediately, "Q" for quiet mode ("R" for taking range as defined in InitFit2P())
					//change upper fitting range for next iteration
					rangeup = N1fit->GetX(0.9992); 	//fit until 99.2% eff. of current trigger (test with 99.92%)

					N2hist->Fit("N2fit", "Q", "", N2rangelow, N2rangeup);
					N2rangeup = N2fit->GetX(0.9992);
				}
				//rangelow = N1fit->GetX(0.90); 	//next iteration will be started at 90% eff. of previous trigger
				rangelow = N1fit->GetX(0.80); 	//next (trigger-)iteration will be started at 80% eff. of previous trigger
				if(ind>4){
					N2rangelow_vec.push_back(N2fit->GetX(0.80));	//will not be used in next, but next to next iteration //first iteration (ind=4) is trtrgg60/trg25, which is still nonsense
				}
				N1rangeup_vec.push_back(N1fit->GetX(0.9992));	//first one (ind=4) is 60GeV trigger at 99.92% eff.
				N2rangelow = N2rangelow_vec.at(ind-4);		//explain: ---


				TDirectory *curdirFit = standard_folder->GetDirectory(Form("%s/trigger_ratios", cur_ybysbin.get_bname().c_str()));
				N1hist->SetDirectory(curdirFit);	//why does this alone (without cd) have no effect on where the histo is written??
				curdirFit->cd();
				N1hist->Write(Form("%s_NewFit", N1hist->GetName()));
				N2hist->SetDirectory(curdirFit);
				//N2hist->Write("", TObject::kOverwrite);
				N2hist->Write(Form("%s_test", N2hist->GetName()));

				//still within trigger loop
				//create a multigraph to use N/N-2 for the lower part of the turn-on curve and N/N-1 result for the part near the plateau
				//mg contains both histograms, "converted" to TGraphErrors graphs and cut corresponding to their "good" regions
				TMultiGraph *mg_cut = MakeMultigraphCut(N1hist, N2hist, N1fit, N2fit, triggers.at(ind));

	/*
				//still within trigger loop
				//create a multigraph to use N/N-2 for the lower part of the turn-on curve and N/N-1 result for the part near the plateau
				//mg contains both histograms, "converted" to TGraphErrors graphs
				TMultiGraph *mg = new TMultiGraph();

				int nbinsN1 = N1hist->GetNbinsX(); // == nbinsN2 --> number of bins is set at initial creation and equals for these histos
				/////int nbinsN1 = N1hist->GetEntries();
				//--> problem: this includes bins with 0 bin content. Remove those. (-->how and when?)
				Double_t n1hist_x[nbinsN1];
				Double_t n1hist_y[nbinsN1];
				Double_t n1hist_err[nbinsN1];
				Double_t n2hist_x[nbinsN1];
				Double_t n2hist_y[nbinsN1];
				Double_t n2hist_err[nbinsN1];

				for(int bin=1; bin!=nbinsN1+1; ++bin){
					n1hist_x[bin]=N1hist->GetXaxis()->GetBinCenter(bin);
					n1hist_y[bin]=N1hist->GetBinContent(bin);
					n1hist_err[bin]=N1hist->GetBinError(bin);

					n2hist_x[bin]=N2hist->GetXaxis()->GetBinCenter(bin); //--> should actually be the same as n2hist_x values
					n2hist_y[bin]=N2hist->GetBinContent(bin);
					n2hist_err[bin]=N2hist->GetBinError(bin);
				}
				TGraphErrors *N1hist_graph = new TGraphErrors(nbinsN1, &n1hist_x[0], &n1hist_y[0], 0, &n1hist_err[0]);
				TGraphErrors *N2hist_graph = new TGraphErrors(nbinsN1, &n2hist_x[0], &n2hist_y[0], 0, &n2hist_err[0]); //same amount of bins as n1hist

				N1hist_graph->SetLineColor(kGreen+2);
				N1hist_graph->SetMarkerStyle(kOpenTriangleUp);
				N1hist_graph->SetMarkerColor(kGreen+2);
				N1hist_graph->SetMarkerSize(1.2);
				N2hist_graph->SetLineColor(kBlue+2);
				N2hist_graph->SetMarkerStyle(kOpenCircle);
				N2hist_graph->SetMarkerColor(kBlue+2);
				N2hist_graph->SetMarkerSize(1.2);

				mg->Add(N1hist_graph);	//Add N/N-1 histo
				mg->Add(N2hist_graph);	//Add N/N-2 histo

				mg->SetTitle(Form("N/N-1 and N/N-2 trigger-ratios; p_{T,avg}; Ratio of trg%d to (N-1) and (N-2)", triggers.at(ind)));
				mg->Draw("APZ");

				//create new multigraph, containing only the points wanted for the fit
				//remove all the points in N1hist_graph, where x-value is lower than N1fit->GetX(0.4) 
				//remove all the points in N2hist_graph, where x-value is higher than N2fit->GetX(0.8)
				Double_t lowxN1 = N1fit->GetX(0.4);
				Double_t upxN2	= N2fit->GetX(0.8);

				if(v==true){ // testing --> should be checked (not only if v==true)... do the correct points get removed? are too many lost? --> adjusting cuts?
					cout << "remove all N1-points lower than x: lowxN1 ===>>> " << lowxN1 << endl;
					cout << "remove all N2-points higher than x: upxN2 ===>>> " << upxN2 << endl;
				}
				//--> this procedure leads to problems when previous fit already had problems!!! --> needs to be improved and changed!

			//just testing for now
			//if(ybysnr==0 && ind==13){}
				for(int bin_i=0; bin_i!=N1hist_graph->GetN(); ++bin_i){
					//cout << "X: " << N1hist_graph->GetX()[bin_i] << endl;
				}
				//cout << "Now remove points." << endl;
				for(int bin_i=0; bin_i!=(N1hist_graph->GetN()); ++bin_i){
				//for(int bin_i=0; bin_i!=nbinsN1;){ //try without set increment
					//--bin_i;
					//cout << "point index = " << bin_i << "|| X-value: " << N1hist_graph->GetX()[bin_i] << endl;
					//if((N1hist_graph->GetX()[bin_i])<lowxN1){
					while((N1hist_graph->GetX()[bin_i])<lowxN1){	//--> THIS DOES NOT WORK WHEN EARLY FLUCTUATIONS APPEAR	 --> better use some "if" condition.
						//cout << Form("%f is lower than %f", N1hist_graph->GetX()[bin_i], lowxN1) << endl;
						//cout << Form("X-value of bin_i=%d before removal of point: X= %f", bin_i, N1hist_graph->GetX()[bin_i]) << endl;
						N1hist_graph->RemovePoint(bin_i);
						//cout << Form("X-value of bin_i=%d AFTER removal of point: X= %f", bin_i, N1hist_graph->GetX()[bin_i]) << endl;

						//need to set current index lower / check this index again, 
						//because all the points above bin_i moved to a new index (ind_new = ind_old-1)
						//--bin_i;
					}
					//else{bin_i++;}	//only increment if points are not "moving downwards" indexwise anymore
				}
				//due to new indices, and the fact that N1hist_graph contains now less points than before -->
				//need an extra loop for N2hist_graph!
				//cout << "N2hist_graph->GetN() --> " << N2hist_graph->GetN() << endl;
				for(int bin_i=0; bin_i<(N2hist_graph->GetN()); ++bin_i){ //otherwise end up in infinite loop... (so use "<" instead of "!=")
					//cout << Form("In bin_i= %d, X-value of N2hist_graph: %f // reminder of upxN2: %f", bin_i, (N2hist_graph->GetX()[bin_i]), upxN2) << endl;

					if((N2hist_graph->GetX()[bin_i])>upxN2){
					//while((N2hist_graph->GetX()[bin_i])>upxN2){
						//cout << Form("X-value of bin_i=%d before removal of point: X= %f", bin_i, N2hist_graph->GetX()[bin_i]) << endl;

						N2hist_graph->RemovePoint(bin_i);
						//cout << "cut from above. ##### " << endl;

						//cout << Form("X-value of bin_i=%d AFTER removal of point: X= %f", bin_i, N2hist_graph->GetX()[bin_i]) << endl;
						--bin_i;
					}
				}


				//CHECK what points are still in graphs:
				int leftpoints = N1hist_graph->GetN();
				//cout << "Still in N1hist_graph. (Cut lower part):" << endl;
				//cout << "----------------------------------------" << endl;
				//cout << leftpoints << " points left out of " << nbinsN1 << endl;
				for(int bin_i=0; bin_i!=N1hist_graph->GetN(); ++bin_i){
					//cout << "X: " << N1hist_graph->GetX()[bin_i] << endl;
				}
			//}

				TMultiGraph *mg_cut = new TMultiGraph();
				//////TCanvas *mgcut_canv = new TCanvas("mgcut_canv", "mgcut_canv", 600, 400);
				mg_cut->Add(N1hist_graph);		//add cut graph to new multigraph
				mg_cut->Add(N2hist_graph);		//add cut graph to new multigraph

				mg_cut->SetTitle(Form("Selected Ranges of: N/N-1 and N/N-2 trigger-ratios; p_{T,avg}; Ratio of trg%d to (N-1) and (N-2)", triggers.at(ind)));
				mg_cut->Draw("APZ");
	*/

				//check ranges again... maybe rename them, so that it is clearer what ranges are used.
				//following function is used to fit the multigraph made out of N/N-1 and N/N-2 histograms
				TF1 *combifit = InitFit2P("combifit", "0.5*(1+TMath::Erf((x-[0])/[1]))", N2rangelow, rangeup, kRed+2, triggers.at(ind), TMath::Sqrt(triggers.at(ind)), "pnom", "pwid");

				if(v==true){
					cout << Form("N2rangelow: %f // rangeup: %f", N2rangelow, rangeup) << endl;
					cout << Form("N2fit->GetX(0.3): %f // N1fit->GetX(0.992): %f", N2fit->GetX(0.3), N1fit->GetX(0.992)) << endl;
				}			
				//mg_cut->Fit("combifit", "Q", "", N2rangelow_vec.at(ind-4), N2rangeup); //reuse previously calculated fittin granges (see N2fit procedure)
				
				//fit from 80% eff. of trg(N-2) up to 99.92% eff. of trg(N-1)
				//this does only make sense for ind>4 --> starting at ind=5 --> 80GeV trigger to 60 (ind=4) and 40 (ind=3)
				//if(ind>4){
					mg_cut->Fit("combifit", "Q", "", N2rangelow_vec.at(ind-4), N1rangeup_vec.at(ind-4)); //reuse previously calculated fitting ranges (see N2fit procedure)
				//}
				/*
				else{
					//use other upper range for Fit, but this won't be sensible anyway, because trg25 is not existent
					mg_cut->Fit("combifit", "Q", "", N2rangelow_vec.at(ind-4), N2rangeup); 
				}
				*/

				//get the parameters from the combifit:
				///Double_t combifit_pnom = combifit->GetParameter(0);
				///Double_t combifit_pwid= combifit->GetParameter(1);

				//append to x-axis and y-axis vectors for the parameter-stability-check

				//for plot of pnom/nom versus nom
				//only append AFTER trg60 (because in N2 method trg60/trg25 is not usefule as trg25 does not exist)
				if(triggers.at(ind)!=60){ //== ind>4
					Store2Params(triggers.at(ind), combifit, mg_pnom_X, mg_pnom_Y, mg_pnom_err, mg_pwid_X, mg_pwid_Y, mg_pwid_err);
				}

				mg_cut->Write(Form("multigraph_CUT_%s", N1hist->GetName()), TObject::kOverwrite);
				///mg_cut->Write(Form("multigraph_CUT_Fit_%s", N1hist->GetName()));	/TESTATESTATESTA!
				combifit->Draw("SAME");

			}//end of fitting loop	
			//trglumi_output->Write();

			if(v==true){
				//CHECK what is in the N2rangelow_vec.
				cout << "-----------------------------" << endl;
				cout << "CHECK N2rangelow_vec and N1rangeup_vec: " << endl;
				cout << "N2rangelow_vec" << endl;
				for(int vecin=0; vecin!=N2rangelow_vec.size(); ++vecin){
					cout << Form("#vecind: %d // #entry: %f ", vecin, N2rangelow_vec.at(vecin)) << endl;
				}
				cout << "N1rangeup_vec" << endl;
				for(int vecin=0; vecin!=N1rangeup_vec.size(); ++vecin){
					cout << Form("#vecind: %d // #entry: %f ", vecin, N1rangeup_vec.at(vecin)) << endl;
				}
			}

			// CHECKING THE PARAMETERS OF 2PARAM FIT!
			
			//##############	method A: N/N-1		################//
			//draw graph for checking the parameters p_nom and p_sqrt
			//make arrays to use TGraph constructor

			//**//TCanvas *pnom_canv = new TCanvas("pnom_canv", "pnom_canv", 600, 400);

			//TGraphErrors *cur_graph_nom = new TGraph(&curx_pnom, &cury_pnom);
			Int_t ntrig = curx_pnom.size();
			Double_t xpnom[ntrig];
			Double_t ypnom[ntrig];
			Double_t xpwid[ntrig];
			Double_t ypwid[ntrig];
			for(int j=0; j!=ntrig; ++j){
				xpnom[j]=curx_pnom.at(j);
				xpwid[j]=curx_pwid.at(j);
				ypnom[j]=cury_pnom.at(j);
				ypwid[j]=cury_pwid.at(j);
			}

			//TGraph *cur_graph_nom = new TGraph(ntrig, xpnom, ypnom);
			TGraphErrors *cur_graph_nom = new TGraphErrors(ntrig, &curx_pnom[0], &cury_pnom[0], 0, &pnom_err[0]);
			cur_graph_nom->SetName(Form("%s_pnom", (cur_ybysbin.get_bname()).c_str()));

			//specific style choices in drawing function for better comparison with combifit (choose different style and color)
			SetDrawStyleWrite(cur_graph_nom, "AP", "p_nom/nom (N1 method --> N/N-1)", "nom [GeV]", "p_nom / nom", kAzure-2, 2, kOpenSquare, kAzure-2, 2);	//plotting of pnom and Write()

			TGraphErrors *cur_graph_wid = new TGraphErrors(ntrig, &curx_pwid[0], &cury_pwid[0], 0, &pwid_err[0]);
			cur_graph_wid->SetName(Form("%s_pwid", (cur_ybysbin.get_bname()).c_str()));

			SetDrawStyleWrite(cur_graph_wid, "AP", "p_wid/nom (N1 method --> N/N-1)", "nom [GeV]", "p_wid / nom", kAzure-2, 2, kOpenSquare, kAzure-2, 2);	//plotting of pwid and Write()



			//##############	method C: N/N-1 and N/N-2 combined	##########//
			// same as above --> graph for checking the fit paramters pnom and pwid
			// create arrays that will be used in TGraph constructor
			//**//TCanvas *combi_pnom_canv = new TCanvas("combi_pnom_canv", "combi_pnom_canv", 600, 400);

			//Int_t ntrig = curx_pnom.size(); //-->stays the same throughout the whole analysis here
			//remove trg60, because it cannot be used (N2--> trg60/trg25 and trg25 is not existent)


			TGraphErrors *graph_combi_nom = new TGraphErrors(ntrig-1, &mg_pnom_X[0], &mg_pnom_Y[0], 0, &mg_pnom_err[0]);	//start at ind=0, but leave out trg60 already earlier
			graph_combi_nom->SetName(Form("combi_%s_pnom", (cur_ybysbin.get_bname()).c_str()));
			SetDrawStyleWrite(graph_combi_nom, "AP", "p_nom/nom (combifit N2 and N1)", "nom [GeV]", "p_nom / nom", kRed-2, 2, kOpenCircle, kRed-2, 2);


			TGraphErrors *graph_combi_wid = new TGraphErrors(ntrig-1, &mg_pwid_X[0], &mg_pwid_Y[0], 0, &mg_pwid_err[0]);	//start at ind=0, leave out trg60
			graph_combi_wid->SetName(Form("combi_%s_pwid", (cur_ybysbin.get_bname()).c_str()));
			SetDrawStyleWrite(graph_combi_wid, "AP", "p_wid/nom (combifit N2 and N1)", "nom [GeV]", "p_wid / nom", kRed-2, 2, kOpenCircle, kRed-2, 2);

				
			//1-param-fit
			//fitting pnom/nom curve in order to reduce 2param fit to 1param fit --> obtain pnom, which will be fixed in upcoming fit
			TF1 *pnomFit = InitFit2P("pnomFit", "[0]+[1]*x", 70, 560, kCyan-2, 1.04, -0.00014, "intercept", "slope");

			cout << "############## FITTING PNOM ###############" << endl;
			graph_combi_nom->Fit("pnomFit", "RQ", ""); //Q for quiet, R for prespecified range
			graph_combi_nom->Write("", TObject::kOverwrite);
			gStyle->SetOptFit(1);

			//now fitting pwid/nom --> later used to get pwid values for 40GeV and 60GeV in order to reproduce those turn-on curves as well
			TF1 *pwidFit = InitFit2P("pwidFit", "[0]+[1]*x", 70, 460, kCyan-2, 0.18, -0.000084, "intercept", "slope");

			cout << "############## FITTING PWID ###############" << endl;
			graph_combi_wid->Fit("pwidFit", "RQ", ""); //Q for quiet, R for prespecified range
			graph_combi_wid->Write("", TObject::kOverwrite);
			gStyle->SetOptFit(1);
		

			/*	
			//summary of the turnon-curves --> prepare canvas and legend
			TCanvas *sum_turnon = new TCanvas("sum_turnon", "sum_turnon", 600, 400); //summary canvas for all the turnon-curves of current ybysbin

			//make a legend for the turnon-curves
			TLegend *turnon_leg = tdrLeg(0.76, 0.2, 0.96, 0.5);

			//turnon_leg->SetEntrySeparation(0.5); //does not work as desired..
			//gStyle->SetLegendTextSize(0.002);
			turnon_leg->SetTextFont(42);
			turnon_leg->SetTextAlign(31);
			*/


			//save the obtained 0.99 efficiency pTavg-threshold to the following vector
			vector<Double_t> trg_eff99_vec;
			map<int, Double_t> trg_effmap;	//mapping trigger threshold (in pTavg) onto nominal trigger value (int)
			

			//redo combifit, but this time fix the pnom to the value obtained via pnomFit:
			///for(int ind=5; ind!=triggers.size(); ++ind){}	//start at ind=5 (80GeV trigger --> 80/40 and 80/60)
			for(int ind=3; ind!=triggers.size(); ++ind){	//start at ind=3 (40GeV trigger --> get fixed parameters from previous fit - 40 and 60 from extrapolation)
				//get already created combi multigraph (by name)
				//fit again, this time fix parameter pnom
				Double_t nom = triggers.at(ind);
				Double_t pnom_fix = (pnomFit->Eval(nom))*nom;	//--> should be possible also for 40GeV and 60GeV --> extrapolation to lower pT
				//cout << "TEST: Eval at 40 and 60: " << (pnomFit->Eval(40))*40 << " -- " << (pnomFit->Eval(60))*60 << endl;

				if(v==true){cout << "current pnom from fit: " << pnom_fix << endl;}
				//TF1* fit_1param = new TF1("fit_1param", "[0]*0.5*(1+TMath::Erf((x-[1])/[2]))");
				TF1 *fit_1param = new TF1("fit_1param", Form("0.5*(1+TMath::Erf((x-%f)/[0]))", pnom_fix));
				fit_1param->SetParameter(0,TMath::Sqrt(triggers.at(ind)));
				fit_1param->SetParNames("pwid");


				//fix p0 to 1.0 (=plateau), fix p1 to previous fitresult (=pnom_fix)
				///fit_1param->FixParameter(0, 1.0);
				///fit_1param->FixParameter(1, pnom_fix);

				//only for trigger 80GeV and upwards
				TMultiGraph* combi_graph;
				if(triggers.at(ind)>60){
					trglumi_output->GetObject(Form("Standard/%s/trigger_ratios/multigraph_CUT_yb%01.0fys%01.0f_trg%d", cur_ybysbin.get_bname().c_str(), yb_lo, ys_lo, triggers.at(ind)), combi_graph);
					//TMultiGraph* combi_graph = (TMultiGraph*)trglumi_output->Get(Form("Standard/%s/trigger_ratios/multigraph_CUT_yb%01.0fys%01.0f_trg%d", cur_ybysbin.get_bname().c_str(), yb_lo, ys_lo, triggers.at(ind)))->Clone();
			
					combi_graph->SetName(Form("combi_%s_trg%d", cur_ybysbin.get_bname().c_str(), triggers.at(ind)));
					combi_graph->Draw("AP");

					fit_1param->SetLineColor(kOrange+7);
					fit_1param->SetLineWidth(3);
					gStyle->SetOptFit(1);	//trying if this helps with the display in TBrowser()
					combi_graph->Fit("fit_1param", "Q", "", N2rangelow_vec.at(ind-4), N1rangeup_vec.at(ind-4));
					gStyle->SetOptFit(1);

					combi_graph->Write(Form("1param_combi_%s_trg%d", cur_ybysbin.get_bname().c_str(), triggers.at(ind)));

					fit_1param->Draw("SAME");
				}

				//##### UNDER CONSTRUCTION 12.8.2019	####
				/*
	if(triggers.at(ind)!=60){ //== ind>4
				//collecting the values for plotting the parameters later (where from 40GeV and 60GeV trigger turn-on curves will be reconstructed)
				Double_t 1paramfit_pwid = fit_1param->GetParameter(0);
				pnom_1p_X.push_back(triggers.at(ind));	//nominal value of trigger
				pnom_1p_Y.push_back(pnom_fix/(triggers.at(ind))); //pnom/nom --> here pnom is not obtained from fit_1param, but is input to that fit

				//for plot of pwid/nom versus nom (--> is the more relevant plot here)
				pwid_1p_X.push_back(triggers.at(ind)); //could also plot against sqrt(nom), but then use other y-values
				//cury_pwid.push_back((pow(p_sqrt,2))/(triggers.at(ind))); //(p_sqrt^2)/nom
				pwid_1p_Y.push_back(1paramfit_pwid/(triggers.at(ind))); //(p_sqrt^2)/nom //TRY HERE WITHOUT SQUARING P2

				//y-errors for pnom and pwid plot (stored in vectors)
				///pnom_1p_err.push_back((combifit->GetParError(0))/triggers.at(ind));	--> HOW IS THIS ERROR CALCULATED? from pnom_fit
				pwid_1p_err.push_back((fit_1param->GetParError(0))/triggers.at(ind));
				}
				*/

				//get threshold for trigger efficiency > 0.99 for current trigger
				//write it to vector containing all those "bin bounds"

				//getting pwid for turn-on curve...
				Double_t pwid_fix;	//width used for drawing turn-on curve (Erf)
				if(triggers.at(ind)<80){	//for 40GeV and 60GeV trigger
					pwid_fix = (pwidFit->Eval(nom))*nom;	//obtained from fit of pwid (graph_combi_wid)
				}
				else{
					pwid_fix = fit_1param->GetParameter(0);
				}
				
				//draw function which takes third parameter from fit and fixed 1st and 2nd parameter as in fit_1param
				//TF1 *cur_turnon = new TF1(Form("cur_turnon_%d", triggers.at(ind)), Form("0.5*(1+TMath::Erf((x-%f)/(%f)))", pnom_fix, fit_1param->GetParameter(0)), 40, 800); //range just fixed for testing
				TF1 *cur_turnon = new TF1(Form("cur_turnon_%d", triggers.at(ind)), Form("0.5*(1+TMath::Erf((x-%f)/(%f)))", pnom_fix, pwid_fix), 20, 800); //range just fixed for testing

				//TF1 *fit_1p_res = combi_graph->GetFunction("fit_1param");
				//Double_t cur_X99eff = fit_1p_res->GetX(0.99);
				Double_t cur_X99eff = cur_turnon->GetX(0.99);
				trg_eff99_vec.push_back(cur_X99eff);	//append this to efficiency-Xpoint-vector
				trg_effmap.insert(make_pair(triggers.at(ind),cur_X99eff));

				if(ybysnr==0){trg_effmap_yb0ys0.insert(make_pair(triggers.at(ind),cur_X99eff));}	//test! trg_effmap_yb0ys0 for all ybys bins

				if(v==true){
					cout << endl;
					cout << Form("looking at trigger: %d GeV", triggers.at(ind)) << endl;
					cout << "Eff. > 0.99 at " << cur_X99eff << endl;
					cout << "----------------------------------------" << endl;
				}

				vector<Color_t> cols = {kYellow+1, kGray+1, kGreen+2, kBlue-2, kOrange-3, kRed+2, kAzure+1, kViolet+6, kSpring+5, kViolet-5, kCyan-3};
				//cur_turnon->SetTitle(Form("Trigger Turn-On Curve: %d GeV;p_{T,avg} /GeV;Trigger-Ratio to (N-1) and (N-2)", triggers.at(ind)));
				cur_turnon->SetTitle(";p_{T,avg} /GeV;Trigger-Ratio N to (N-1) and (N-2)");
				cur_turnon->SetLineWidth(3);
				///cout << "color vec size:" << cols.size() << endl;
				cur_turnon->SetLineColor(cols.at(ind-3));	//because first trigger is 40GeV (at ind=3)

				cur_turnon->Write();	//write to root-file

				/*
				sum_turnon->cd();	//change to summary canvas
				cur_turnon->DrawCopy();
				sum_turnon->Modified();
				sum_turnon->Update();
				*/
				///TDirectory *curdirFit = standard_folder->GetDirectory(Form("%s/trigger_ratios", cur_ybysbin.get_bname().c_str()));

				// add an entry to legend for the turnon-curves
				//TLegend *turnon_leg = tdrLeg(0.8, 0.1, 0.98, 0.4);
				///cout << "Separation: " << turnon_leg->GetEntrySeparation() << endl;

				///////////turnon_leg->AddEntry(cur_turnon->GetName(), Form("%d GeV", triggers.at(ind)), "l");

			} //end triggerloop (started at ind=5 --> 80GeV trigger (80/60 and 80/40))
			//write summary canvas to root-file
			/*
			sum_turnon->cd();
			turnon_leg->Draw();
			sum_turnon->Write(Form("sum_turnon_%s", cur_ybysbin.get_bname().c_str()));
			sum_turnon->Close();
			*/

			//use the results from the 1paramfit and the drawn Turn-On Curves!
			//fill new histogram with ptavg values by using the new "99%-eff. thresholds" and the previously filled histograms for each trigger
			//create a new TH1D histogram for current ybys bin --> shall this be a class object? --> YES.
			//TH1D *sumhisto
			trglumi_output->cd(cur_ybysbin.get_dirname());
			//below: how it should be done, but currently testing with same trg_effmap for all the ybys bins, see further below
			////TH1D* cur_sumhist = cur_ybysbin.fill_sumhisto(trg_effmap, triggers, histosvec.at(ybysnr), lumivec);	//also give it lumis for normalisation
			////cur_sumhist->Write();
		
			//test: take ALWAYS trg_effmap of central ybys bin (used until trigger fits have been improved)
			//just an intermediate solution
			TH1D* cur_sumhist = cur_ybysbin.fill_sumhisto(trg_effmap_yb0ys0, triggers, histosvec.at(ybysnr), lumivec);	//also give it lumis for normalisation
			cur_sumhist->Write();
				
		}//end ybys-loop
		cout << "Saved results to file: trg_lumi_output.root" << endl;
		//trglumi_output->Write();	//writes all the currently open objects (so doubles some)



	//---------------------------------------------------------------------------------------------------------------------------------------------------//
	//|																																					|//
	//|														ALTERNATIVE	FITTING METHOD																	|//
	//|											Take (HLT_n + TrigObj_ptCUT_n+1) versus HLT_n															|//
	//|																																					|//
	//---------------------------------------------------------------------------------------------------------------------------------------------------//

		//################### NEW: CHECK ALTERNATIVE TRIGGER FIT METHOD - under construction! ####################//
		//alternative trigger-fitting method (TO-DO: write general functions to carry out the fitting more modular)

		//for summary of the alternative method concerning the central bin
		map<int, Double_t> trg_effmap_yb0ys0_alt;		//mapping trigger threshold (in pTavg) onto nominal trigger value (int)
		map<int, Double_t> trg_effmap_yb0ys0_3p_alt;	//trigger thresholds obtained after one 3-parameter fit of the emulation method (N/N-1)

		//go once more through ybys bins
		for (int ybysnr=0; ybysnr!=6; ++ybysnr){
			YbYs cur_ybysbin = ybys_pointer[ybysnr];
			float yb_lo = cur_ybysbin.yb_get_bb().first;
			float yb_up = cur_ybysbin.yb_get_bb().second;
			float ys_lo = cur_ybysbin.ys_get_bb().first;
			float ys_up = cur_ybysbin.ys_get_bb().second;

			
			//handle directories in output file --> write alternative ratio plots in corresponding ybys folder
			//create new directory for the alternative trigger fits
			const char *dirname = Form("%s", cur_ybysbin.get_bname().c_str());
			const char *alternative_dirname = Form("alternative_trigger_fits");
			standard_folder->cd(dirname);
			//standard_folder->mkdir(dirname);
			TDirectory *ybysDir = standard_folder->GetDirectory(dirname);
			ybysDir->cd();
			ybysDir->mkdir(alternative_dirname);
			TDirectory *altDir = ybysDir->GetDirectory(alternative_dirname);


			//later: create graphs that are used to check which parameter is more stable in the fit
			//here: create the vectors that will be used to store x- and y- coordinates of the graphs
			//pnom = parameter which is approx. nominal value of trigger threshold
			//pwid = parameter which is approx. width of errorfunction (turn-on-curve) (=sqrt(pnom))
			// for initial 3-param fit
			vector<Double_t> pnom_3p_x;
			vector<Double_t> pnom_3p_y;
			vector<Double_t> pwid_3p_x;
			vector<Double_t> pwid_3p_y;

			vector<Double_t> pnom_3p_err;
			vector<Double_t> pwid_3p_err;

			// again for 2-param fit
			vector<Double_t> pnom_2p_x;
			vector<Double_t> pnom_2p_y;
			vector<Double_t> pwid_2p_x;
			vector<Double_t> pwid_2p_y;

			vector<Double_t> pnom_2p_err;
			vector<Double_t> pwid_2p_err;


			//same procedure as well for multigraph (mg_cut)
			vector<Double_t> mg_pnom_X;
			vector<Double_t> mg_pnom_Y;
			vector<Double_t> mg_pwid_X;
			vector<Double_t> mg_pwid_Y;

			vector<Double_t> mg_pnom_err;
			vector<Double_t> mg_pwid_err;

			//and once more for the 1param fit in the end
			vector<Double_t> pnom_1p_X;
			vector<Double_t> pnom_1p_Y;
			vector<Double_t> pwid_1p_X;
			vector<Double_t> pwid_1p_Y;

			vector<Double_t> pnom_1p_err;
			vector<Double_t> pwid_1p_err;


			//now divide histograms
			//leave out the lumi-scaling so far

			//always divide "simulated trigger" N by "HLT_PFJet trigger" (N-1)
			//start at trigger-index "4" which is jt60 (ratio to ind=3 --> jt40), because jt15 (ind=1) and jt25 (ind=2) have lumi=0
			//luminosity values start from index "0" nevertheless
			//for(int ind=2; ind!=triggers.size(); ++ind){}
			for(int ind=4; ind!=triggers.size(); ++ind){
				TH1D* prehist = (TH1D*)histosvec.at(ybysnr)->jt_histos[triggers.at(ind-1)]->Clone(Form("yb%01.0fys%01.0f_trg%d", yb_lo, ys_lo, triggers.at(ind-1))); //previous histogram (N-1)
				TH1D* curhist = (TH1D*)histosvec.at(ybysnr)->trgobj_histos[triggers.at(ind)]->Clone(Form("yb%01.0fys%01.0f_simutrg%d", yb_lo, ys_lo, triggers.at(ind)));		//current ("simulated") histogram (N), that will be divided by hist(N-1)
				//for later fits: also prepare the N/N-2 ratio in this directory
				//newly added: also check N/N-2 --> ratio to prevprevtrigger
				TH1D* prevprev_ratio = (TH1D*)histosvec.at(ybysnr)->n2_trgobj_histos[triggers.at(ind)]->Clone(Form("yb%01.0fys%01.0f_simutrg%d_to_trg%d", yb_lo, ys_lo, triggers.at(ind), triggers.at(ind-2))); //current histogram (N)
				//simutrg_N2ratio
				TH1D* preprehist = (TH1D*)histosvec.at(ybysnr)->jt_histos[triggers.at(ind-2)]->Clone(Form("yb%01.0fys%01.0f_trg%d_asPrePreTrg", yb_lo, ys_lo, triggers.at(ind-2))); //preprevious histogram (N-2)


				//setting the ratio plots' directory to altDir, which has been created for the alternative-method-results of this current ybysbin
				curhist->SetDirectory(altDir);
				prevprev_ratio->SetDirectory(altDir);
				altDir->cd();		//change to directory where the alternative trigger fits are stored


				//Taking care of N/N-1
				curhist->Sumw2();	//important for error handling  --> see TH1D::Divide() in ROOT documentation
				prehist->Sumw2();
				curhist->Divide(curhist, prehist, 1, 1, "B");	//divide current histo by previous one (no lumi scaling needed)
				//TH1D* simutrg_N1ratio = curhist->Divide(curhist, prehist, 1, 1, "B");		//divide current histo by previous one (both already scaled by lumi)
				//simutrg_N1ratio->SetDirectory(altDir);

				//adjusting plotting for N/N-1
				SetDrawStyleWrite(curhist, "", "", "p_{T,avg}", Form("#frac{trg%d && TrigObj_pt(trg%d)}{trg%d}", triggers.at(ind-1), triggers.at(ind), triggers.at(ind-1)), kGray+2, 2, kOpenCircle, kGray+2, true);

				//Taking care of N/N-2
				prevprev_ratio->Sumw2();
				preprehist->Sumw2();
				prevprev_ratio->Divide(prevprev_ratio, preprehist, 1, 1, "B"); //divide current histo by preprevious one (Binomial statistics) (no lumi scaling needed)

				//adjusting plotting for N/N-2
				SetDrawStyleWrite(prevprev_ratio, "", "", "p_{T,avg}", Form("#frac{trg%d && TrigObj_pt(trg%d)}{trg%d}", triggers.at(ind-2), triggers.at(ind), triggers.at(ind-2)), kGray+2, 2, kOpenCircle, kGray+2, true);

	//-----
	//-------

				//#### first fits of this method --> for now just very simple --> three parameters and taking very rough fitting range (guessed!!!) ####//
				//Now: Fit erf() onto curve
				//set range, where to fit (so far just a guess..)
				Double_t rangelow = triggers.at(ind-1);
				Double_t rangeup = 1.5*triggers.at(ind);

				// Approx. the values of the fitting parameters: 
				// 	taken as initial fit parameter values
				//	[0]-->efficiency at saturation, here assumed to be "1.0"
				//	[1]-->position of turn-on (50% eff. // nominal value), 
				//	[2]-->(1./width)=Sqrt([1]) or 5% at high pT

				// 3-param-fit
				TF1 *erf = new TF1("erf", "[0]*0.5*(1+TMath::Erf((x-[1])/[2]))", rangelow, rangeup);
				//erf->SetParameter(0, 1.0); //--> set parameter [0] to 1.
				erf->SetParameters(1.0, triggers.at(ind), TMath::Sqrt(triggers.at(ind))); //set initial p1 and p2 (see above)
				erf->SetParLimits(0,0,1.0);
				erf->SetParNames("plateau", "pnom", "pwid");
				
				erf->SetLineWidth(2);
				curhist->Fit("erf", "BQ", "", rangelow, rangeup); //"Q" for quiet --> minimal output //"B" for "bounds" --> using limit for parameter

				Double_t X99eff_3p_alt = erf->GetX(0.99);
				if(ybysnr==0){trg_effmap_yb0ys0_3p_alt.insert(make_pair(triggers.at(ind),X99eff_3p_alt));}	//test! trg_effmap_yb0ys0_3p_alt for all ybys bins


				//get the parameters from the fit:
				Store2Params(triggers.at(ind), erf, pnom_3p_x, pnom_3p_y, pnom_3p_err, pwid_3p_x, pwid_3p_y, pwid_3p_err);

				curhist->Write(Form("%s_3pFit_N1",curhist->GetName()));
				gStyle->SetOptFit(1);

				//for now (will access this later, need to write it first):
				prevprev_ratio->Write(Form("%s_N2", prevprev_ratio->GetName()));

			}//end trigger loop

			// CHECKING THE PARAMETERS OF 3PARAM FIT of alternative fitting version!
			
			//##############	method A: N/N-1		################//
			//draw graph for checking the parameters p_nom and p_wid
			// naming: 'alt' == alternative

			Int_t ntrig_3p = pnom_3p_x.size();
			TGraphErrors *cur_graph_nom_3p = new TGraphErrors(ntrig_3p, &pnom_3p_x[0], &pnom_3p_y[0], 0, &pnom_3p_err[0]);
			cur_graph_nom_3p->SetName(Form("%s_pnom_alt_3p", (cur_ybysbin.get_bname()).c_str()));
			SetDrawStyleWrite(cur_graph_nom_3p, "AP", "p_nom/nom (N1 method (alternative, 3param) --> N/N-1)", "nom [GeV]", "p_nom / nom", kAzure-2, 2, kOpenSquare, kAzure-2, 2);	//plotting of pnom and Write()

			TGraphErrors *cur_graph_wid_3p = new TGraphErrors(ntrig_3p, &pwid_3p_x[0], &pwid_3p_y[0], 0, &pwid_3p_err[0]);
			cur_graph_wid_3p->SetName(Form("%s_pwid_alt_3p", (cur_ybysbin.get_bname()).c_str()));
			SetDrawStyleWrite(cur_graph_wid_3p, "AP", "p_wid/nom (N1 method (alternative, 3param) --> N/N-1)", "nom [GeV]", "p_wid / nom", kAzure-2, 2, kOpenSquare, kAzure-2, 2);	//plotting of pwid and Write()

			//cur_graph_nom_3p->Write("", TObject::kOverwrite);
			//cur_graph_wid_3p->Write("", TObject::kOverwrite);


	//------%%----

			//------------------------------------------------------------------------------------------------------------------------------//
			//				N/N-1 and N/N-2 fits (seperately) for the emulation-method (=alternative method)								//
			//------------------------------------------------------------------------------------------------------------------------------//
			//--> completely analogous to Method one (HLT ratio)

			//###################### METHOD B --> fit from 80% of previous trigger up to 99.92% of current trigger ##############//
			//##### additionally: N/N-2
			//new fitting loop for using information of previous trigger-fit for the next trigger
			//use very big range initially
			Double_t rangelow = 40.00;
			Double_t rangeup = 630.00;

			Double_t N2rangelow = 30.00;
			Double_t N2rangeup = 630.00;
			vector<Double_t> N2rangelow_vec={31.0, 48.0, 70.0};	//needs to start with X-position: 25GeV trigger @ 80% (guessed with fit)
			vector<Double_t> N1rangeup_vec;	//stores upper fitrange-limit obtained from N1-fit (=N1fit)

			TDirectory *curdirFit = standard_folder->GetDirectory(Form("%s/alternative_trigger_fits", cur_ybysbin.get_bname().c_str()));
			standard_folder->cd();
			//start at ind=4 (==>> 60GeV-trigger (over 40GeV --> ind-1))
			for(int ind=4; ind!=triggers.size(); ++ind){
				//initial fitting range for current "trigger-pair"
				//use very big range initially
				//Double_t rangelow = 40.00;
				Double_t rangeup = 630.00;	//rangeup will be set according to current trigger, rangelow has been obtained from previous trigger
				Double_t N2rangeup = 630.00;

				//get the N/N-1 histo
				TH1D* N1hist; 
				trglumi_output->GetObject(Form("Standard/%s/alternative_trigger_fits/yb%01.0fys%01.0f_simutrg%d_3pFit_N1", cur_ybysbin.get_bname().c_str(), yb_lo, ys_lo, triggers.at(ind)), N1hist);

				//get the N/N-2 histo
				TH1D* N2hist;
				const char* prevprevname =Form("yb%01.0fys%01.0f_simutrg%d_to_trg%d_N2", yb_lo, ys_lo, triggers.at(ind), triggers.at(ind-2));
				trglumi_output->GetObject(Form("Standard/%s/alternative_trigger_fits/%s", cur_ybysbin.get_bname().c_str(), prevprevname), N2hist);

				// 2-param-fit for fitting of N/N-1 (set [0] to 1.0, renaming of [1]-->[0] and [2]-->[1])
				TF1 *N1fit = InitFit2P("N1fit", "0.5*(1+TMath::Erf((x-[0])/[1]))", rangelow, rangeup, kGreen+2, triggers.at(ind), TMath::Sqrt(triggers.at(ind)), "pnom", "pwid");
				
				//same 2-param-fit function for fitting of N/N-2
				TF1 *N2fit = InitFit2P("N2fit", "0.5*(1+TMath::Erf((x-[0])/[1]))", N2rangelow, N2rangeup, kBlue+2, triggers.at(ind), TMath::Sqrt(triggers.at(ind)), "pnom", "pwid");

				//fit iteratively
				for(int it_ind=0; it_ind!=5; ++it_ind){
					N1hist->Fit("N1fit", "Q", "", rangelow, rangeup); //"0" for not drawing immediately, "Q" for quiet mode ("R" for taking range as defined in InitFit2P())
					//cout << "fitted N1fit to N1hist." << endl;
					//change upper fitting range for next iteration
					rangeup = N1fit->GetX(0.9992); 	//fit until 99.2% eff. of current trigger (test with 99.92%)

					N2hist->Fit("N2fit", "Q", "", N2rangelow, N2rangeup);
					N2rangeup = N2fit->GetX(0.9992);
				}
				//rangelow = N1fit->GetX(0.90); 	//next iteration will be started at 90% eff. of previous trigger
				rangelow = N1fit->GetX(0.80); 	//next (trigger-)iteration will be started at 80% eff. of previous trigger
				if(ind>4){
					N2rangelow_vec.push_back(N2fit->GetX(0.80));	//will not be used in next, but next to next iteration //first iteration (ind=4) is trtrgg60/trg25, which is still nonsense
				}
				N1rangeup_vec.push_back(N1fit->GetX(0.9992));	//first one (ind=4) is 60GeV trigger at 99.92% eff.
				N2rangelow = N2rangelow_vec.at(ind-4);		//explain: ---


				TDirectory *curdirFit = standard_folder->GetDirectory(Form("%s/alternative_trigger_fits", cur_ybysbin.get_bname().c_str()));
				N1hist->SetDirectory(curdirFit);	//why does this alone (without cd) have no effect on where the histo is written??
				curdirFit->cd();
				N1hist->Write(Form("%s_N1fit", N1hist->GetName()), TObject::kOverwrite);
				N2hist->SetDirectory(curdirFit);
				//N2hist->Write("", TObject::kOverwrite);
				N2hist->Write(Form("%s_N2fit", N2hist->GetName()), TObject::kOverwrite);


				Store2Params(triggers.at(ind), N1fit, pnom_2p_x, pnom_2p_y, pnom_2p_err, pwid_2p_x, pwid_2p_y, pwid_2p_err);	//not really necessary...


			//------------------------------------------------------------------------------------------------------------------------------//
			//					Multigraph-Fit (mg = combine N1&N2) for the emulation-method (=alternative method)							//
			//------------------------------------------------------------------------------------------------------------------------------//
				//--> completely analogous to Method one (HLT ratio)

				//still within trigger loop
				//create a multigraph to use N/N-2 for the lower part of the turn-on curve and N/N-1 result for the part near the plateau
				//mg contains both histograms, "converted" to TGraphErrors graphs and cut corresponding to their "good" regions
				TMultiGraph *mg_cut = MakeMultigraphCut(N1hist, N2hist, N1fit, N2fit, triggers.at(ind));

				//check ranges again... maybe rename them, so that it is clearer what ranges are used.
				//following function is used to fit the multigraph made out of N/N-1 and N/N-2 histograms
				TF1 *combifit = InitFit2P("combifit", "0.5*(1+TMath::Erf((x-[0])/[1]))", N2rangelow, rangeup, kRed+2, triggers.at(ind), TMath::Sqrt(triggers.at(ind)), "pnom", "pwid");

				//fit from 80% eff. of trg(N-2) up to 99.92% eff. of trg(N-1)
				//this does only make sense for ind>4 --> starting at ind=5 --> 80GeV trigger to 60 (ind=4) and 40 (ind=3)
				//if(ind>4){
				mg_cut->Fit("combifit", "Q", "", N2rangelow_vec.at(ind-4), N1rangeup_vec.at(ind-4)); //reuse previously calculated fitting ranges (see N2fit procedure)
				//}

				//get the parameters from the combifit:
				//append to x-axis and y-axis vectors for the parameter-stability-check
				//for plot of pnom/nom versus nom
				//only append AFTER trg60 (because in N2 method trg60/trg25 is not usefule as trg25 does not exist)
				if(triggers.at(ind)!=60){ //== ind>4
						Store2Params(triggers.at(ind), combifit, mg_pnom_X, mg_pnom_Y, mg_pnom_err, mg_pwid_X, mg_pwid_Y, mg_pwid_err);
				}

				mg_cut->Write(Form("multigraph_CUT_%s", N1hist->GetName()), TObject::kOverwrite);
				combifit->Draw("SAME");


			}//end of (trigger) fitting loop


			// CHECKING THE PARAMETERS OF 2PARAM FIT! --> Preparing the final 1-param-fit
			
			//##############	method A: N/N-1		################//
			//draw graph for checking the parameters p_nom and p_sqrt
			Int_t ntrig = pnom_2p_x.size();
			TGraphErrors *cur_graph_nom = new TGraphErrors(ntrig, &pnom_2p_x[0], &pnom_2p_y[0], 0, &pnom_2p_err[0]);
			cur_graph_nom->SetName(Form("%s_2p_pnom_N1", (cur_ybysbin.get_bname()).c_str()));
			//specific style choices in drawing function for better comparison with combifit (choose different style and color)
			SetDrawStyleWrite(cur_graph_nom, "AP", "p_nom/nom (N1 method --> N/N-1, 2param)", "nom [GeV]", "p_nom / nom", kAzure-2, 2, kOpenCross, kAzure-2, 2);	//plotting of pnom and Write()

			TGraphErrors *cur_graph_wid = new TGraphErrors(ntrig, &pwid_2p_x[0], &pwid_2p_y[0], 0, &pwid_2p_err[0]);
			cur_graph_wid->SetName(Form("%s_2p_pwid_N1", (cur_ybysbin.get_bname()).c_str()));
			SetDrawStyleWrite(cur_graph_wid, "AP", "p_wid/nom (N1 method --> N/N-1, 2param)", "nom [GeV]", "p_wid / nom", kAzure-2, 2, kOpenCross, kAzure-2, 2);	//plotting of pwid and Write()


			//############# method B: N/N-2 --> no separate validation, as not used.
			//could do the same as above for method A: N/N-1 done.

			//##############	method C: N/N-1 and N/N-2 combined	##########//
			// same as above --> graph for checking the fit paramters pnom and pwid

			//Int_t ntrig = curx_pnom.size(); //-->stays the same throughout the whole analysis here
			//remove trg60, because it cannot be used (N2--> trg60/trg25 and trg25 is not existent)

			TGraphErrors *graph_combi_nom = new TGraphErrors(ntrig-1, &mg_pnom_X[0], &mg_pnom_Y[0], 0, &mg_pnom_err[0]);	//start at ind=0, but leave out trg60 already earlier
			graph_combi_nom->SetName(Form("combi_%s_pnom", (cur_ybysbin.get_bname()).c_str()));
			SetDrawStyleWrite(graph_combi_nom, "AP", "p_nom/nom (combifit N2 and N1)", "nom [GeV]", "p_nom / nom", kRed-2, 2, kOpenDiamond, kRed-2, 2);

			TGraphErrors *graph_combi_wid = new TGraphErrors(ntrig-1, &mg_pwid_X[0], &mg_pwid_Y[0], 0, &mg_pwid_err[0]);	//start at ind=0, leave out trg60
			graph_combi_wid->SetName(Form("combi_%s_pwid", (cur_ybysbin.get_bname()).c_str()));
			SetDrawStyleWrite(graph_combi_wid, "AP", "p_wid/nom (combifit N2 and N1)", "nom [GeV]", "p_wid / nom", kRed-2, 2, kOpenDiamond, kRed-2, 2);

				
			//1-param-fit
			//fitting pnom/nom curve in order to reduce 2param fit to 1param fit --> obtain pnom, which will be fixed in upcoming fit
			TF1 *pnomFit = InitFit2P("pnomFit", "[0]+[1]*x", 70, 560, kCyan-2, 1.04, -0.00014, "intercept", "slope");
			graph_combi_nom->Fit("pnomFit", "RQ", ""); //Q for quiet, R for prespecified range
			graph_combi_nom->Write("", TObject::kOverwrite);
			gStyle->SetOptFit(1);

			//now fitting pwid/nom --> later used to get pwid values for 40GeV and 60GeV in order to reproduce those turn-on curves as well
			TF1 *pwidFit = InitFit2P("pwidFit", "[0]+[1]*x", 70, 460, kCyan-2, 0.18, -0.000084, "intercept", "slope");
			cout << "############## FITTING PWID ###############" << endl;
			graph_combi_wid->Fit("pwidFit", "RQ", ""); //Q for quiet, R for prespecified range
			graph_combi_wid->Write("", TObject::kOverwrite);
			gStyle->SetOptFit(1);
		
			
			//save the (soon) obtained 0.99 efficiency pTavg-threshold to the following vector
			vector<Double_t> trg_eff99_vec;
			map<int, Double_t> trg_effmap_alt;	//mapping trigger threshold (in pTavg) onto nominal trigger value (int)
			

			//redo combifit, but this time fix the pnom to the value obtained via pnomFit:
			///for(int ind=5; ind!=triggers.size(); ++ind){}	//start at ind=5 (80GeV trigger --> 80/40 and 80/60)
			for(int ind=3; ind!=triggers.size(); ++ind){	//start at ind=3 (40GeV trigger --> get fixed parameters from previous fit - 40 and 60 from extrapolation)
				//get already created combi multigraph (by name)
				//fit again, this time fix parameter pnom
				Double_t nom = triggers.at(ind);
				Double_t pnom_fix = (pnomFit->Eval(nom))*nom;	//--> should be possible also for 40GeV and 60GeV --> extrapolation to lower pT
				//cout << "TEST: Eval at 40 and 60: " << (pnomFit->Eval(40))*40 << " -- " << (pnomFit->Eval(60))*60 << endl;

				if(v==true){cout << "current pnom from fit: " << pnom_fix << endl;}
				//TF1* fit_1param = new TF1("fit_1param", "[0]*0.5*(1+TMath::Erf((x-[1])/[2]))");
				TF1 *fit_1param = new TF1("fit_1param", Form("0.5*(1+TMath::Erf((x-%f)/[0]))", pnom_fix));
				fit_1param->SetParameter(0,TMath::Sqrt(triggers.at(ind)));
				fit_1param->SetParNames("pwid");


				//fix p0 to 1.0 (=plateau), fix p1 to previous fitresult (=pnom_fix)
				///fit_1param->FixParameter(0, 1.0);
				///fit_1param->FixParameter(1, pnom_fix);

				//only for trigger 80GeV and upwards
				TMultiGraph* combi_graph;
				if(triggers.at(ind)>60){
					trglumi_output->GetObject(Form("Standard/%s/alternative_trigger_fits/multigraph_CUT_yb%01.0fys%01.0f_simutrg%d", cur_ybysbin.get_bname().c_str(), yb_lo, ys_lo, triggers.at(ind)), combi_graph);
					//TMultiGraph* combi_graph = (TMultiGraph*)trglumi_output->Get(Form("Standard/%s/trigger_ratios/multigraph_CUT_yb%01.0fys%01.0f_trg%d", cur_ybysbin.get_bname().c_str(), yb_lo, ys_lo, triggers.at(ind)))->Clone();
			
					combi_graph->SetName(Form("combi_%s_simutrg%d", cur_ybysbin.get_bname().c_str(), triggers.at(ind)));
					combi_graph->Draw("AP");

					fit_1param->SetLineColor(kOrange+7);
					fit_1param->SetLineWidth(3);
					gStyle->SetOptFit(1);	//trying if this helps with the display in TBrowser()
					combi_graph->Fit("fit_1param", "Q", "", N2rangelow_vec.at(ind-4), N1rangeup_vec.at(ind-4));
					gStyle->SetOptFit(1);

					combi_graph->Write(Form("1param_combi_%s_simutrg%d", cur_ybysbin.get_bname().c_str(), triggers.at(ind)));
					fit_1param->Draw("SAME");
				}


				//get threshold for trigger efficiency > 0.99 for current trigger
				//write it to vector containing all those "bin bounds"

				//getting pwid for turn-on curve...
				Double_t pwid_fix;	//width used for drawing turn-on curve (Erf)
				if(triggers.at(ind)<80){	//for 40GeV and 60GeV trigger
					pwid_fix = (pwidFit->Eval(nom))*nom;	//obtained from fit of pwid (graph_combi_wid)
				}
				else{
					pwid_fix = fit_1param->GetParameter(0);
				}
				
				//draw function which takes third parameter from fit and fixed 1st and 2nd parameter as in fit_1param
				//TF1 *cur_turnon = new TF1(Form("cur_turnon_%d", triggers.at(ind)), Form("0.5*(1+TMath::Erf((x-%f)/(%f)))", pnom_fix, fit_1param->GetParameter(0)), 40, 800); //range just fixed for testing
				TF1 *cur_turnon = new TF1(Form("cur_turnon_%d", triggers.at(ind)), Form("0.5*(1+TMath::Erf((x-%f)/(%f)))", pnom_fix, pwid_fix), 20, 800); //range just fixed for testing

				//TF1 *fit_1p_res = combi_graph->GetFunction("fit_1param");
				//Double_t cur_X99eff = fit_1p_res->GetX(0.99);
				Double_t cur_X99eff = cur_turnon->GetX(0.99);
				trg_eff99_vec.push_back(cur_X99eff);	//append this to efficiency-Xpoint-vector
				trg_effmap_alt.insert(make_pair(triggers.at(ind),cur_X99eff));

				if(ybysnr==0){trg_effmap_yb0ys0_alt.insert(make_pair(triggers.at(ind),cur_X99eff));}	//test! trg_effmap_yb0ys0_alt for all ybys bins

				if(v==true){
					cout << endl;
					cout << Form("looking at trigger: %d GeV", triggers.at(ind)) << endl;
					cout << "Eff. > 0.99 at " << cur_X99eff << endl;
					cout << "----------------------------------------" << endl;
				}

				vector<Color_t> cols = {kYellow+1, kGray+1, kGreen+2, kBlue-2, kOrange-3, kRed+2, kAzure+1, kViolet+6, kSpring+5, kViolet-5, kCyan-3};
				//cur_turnon->SetTitle(Form("Trigger Turn-On Curve: %d GeV;p_{T,avg} /GeV;Trigger-Ratio to (N-1) and (N-2)", triggers.at(ind)));
				cur_turnon->SetTitle(";p_{T,avg} /GeV;Trigger-Ratio N(emulated) to (N-1) and (N-2)");
				cur_turnon->SetLineWidth(3);
				cur_turnon->SetLineColor(cols.at(ind-3));	//because first trigger is 40GeV (at ind=3)

				cur_turnon->Write();	//write to root-file


			} //end triggerloop (started at ind=5 --> 80GeV trigger (80/60 and 80/40))
			//write summary canvas to root-file

			//use the results from the 1paramfit and the drawn Turn-On Curves!
			//fill new histogram with ptavg values by using the new "99%-eff. thresholds" and the previously filled histograms for each trigger
			//create a new TH1D histogram for current ybys bin --> shall this be a class object? --> YES.
			//TH1D *sumhisto
			trglumi_output->cd(cur_ybysbin.get_dirname());
			//below: how it should be done, but currently testing with same trg_effmap for all the ybys bins, see further below
			////TH1D* cur_sumhist = cur_ybysbin.fill_sumhisto(trg_effmap, triggers, histosvec.at(ybysnr), lumivec);	//also give it lumis for normalisation
			////cur_sumhist->Write();
		
			//test: take ALWAYS trg_effmap of central ybys bin (used until trigger fits have been improved)
			//just an intermediate solution

			//so far not done with this method:
			TH1D* cur_sumhist_trgEmulMeth = cur_ybysbin.fill_sumhisto(trg_effmap_yb0ys0_alt, triggers, histosvec.at(ybysnr), lumivec);	//also give it lumis for normalisation
			cur_sumhist_trgEmulMeth->Write();
		
				
		cout << "Saved results to file: trg_lumi_output.root" << endl;
		//trglumi_output->Write();	//writes all the currently open objects (so doubles some)

		//------%%----

		}//end ybys loop


		//summarising:
		cout << endl << endl;
		cout << "Comparing the RATIO METHOD with the EMULATION METHOD (in central bin yb0y0)" << endl;
		cout << "----------------------------------------------------------------------------" << endl;
		cout << left << setw(10) << "Trigger: |"
		<< setw(20) << " 99% threshold RATIO |" 
		<< setw(20) << " 99% threshold 3P fit (emul.) |"
		<< setw(20) << " 99% threshold EMULATION (same proc. as for ratio)" << endl;
		for(int trgnr=3; trgnr!=triggers.size(); ++trgnr){
			cout << left << setw(20) << Form(" % 3.d GeV |", triggers.at(trgnr))
			<< setw(30) << Form("% 3.1f GeV |", trg_effmap_yb0ys0[triggers.at(trgnr)])
			<< setw(20) << Form("% 3.1f GeV |", trg_effmap_yb0ys0_3p_alt[triggers.at(trgnr)])
			<< setw(20) << Form("% 3.1f GeV", trg_effmap_yb0ys0_alt[triggers.at(trgnr)]) << endl;
		}
		cout << "-----------------------------------------------------------------------------" << endl;

		//just for now to check the trigger results again the final bin bounds:
		vector<int> finalbins = {10, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507, 548, 592, 638, 686, 737, 790, 846, 905, 967, 1032, 1101, 1172, 1248, 1327, 1410, 1497, 1588, 1684, 1784, 1890, 2000, 2116, 2238, 2366, 2500, 2640, 2787, 2941, 3103, 3273, 3450, 3832, 6076, 6389};

		cout << endl;
		cout << "Final binning (after trigger studies with finer binning!):" << endl;
		cout << "[";
		int nbb = finalbins.size();
		for(int i=0; i!=nbb-1; ++i){
			cout << Form("%d, ", finalbins.at(i));
		}
		cout << Form("%d]", finalbins.at(nbb-1)) << endl;



		//TESTING FOR OUTOUT OF TRIGGERTHRESHOLDS AS TXT
		ofstream trg_turnon_file;
		if((strcmp(runperiod,"A")==0)or(strcmp(runperiod,"B")==0)or(strcmp(runperiod,"C")==0)or(strcmp(runperiod,"D")==0)){
			cout << "Run " << runperiod  << endl;
			trg_turnon_file.open(Form("trg_turnons_Run%s.txt", runperiod));
		}
		else{//if no runperiod specified
			trg_turnon_file.open("trg_turnons.txt");
		}

		//now write results of trigger studies into .txt file (ratio-method and emulation-method):
		trg_turnon_file << "These Results refer to Run: " << runperiod << endl;
		trg_turnon_file << "Comparing the RATIO METHOD with the EMULATION METHOD (in central bin yb0ys0)" << endl;
		trg_turnon_file << "----------------------------------------------------------------------------" << endl;
		trg_turnon_file << left << setw(10) << "Trigger: |"
		<< setw(20) << " 99% threshold RATIO |" 
		<< setw(20) << " 99% threshold 3P fit (emul.) |"
		<< setw(20) << " 99% threshold EMULATION (same proc. as for ratio)" << endl;
		for(int trgnr=3; trgnr!=triggers.size(); ++trgnr){
			trg_turnon_file << left << setw(20) << Form(" % 3.d GeV |", triggers.at(trgnr))
			<< setw(30) << Form("% 3.1f GeV |", trg_effmap_yb0ys0[triggers.at(trgnr)])
			<< setw(20) << Form("% 3.1f GeV |", trg_effmap_yb0ys0_3p_alt[triggers.at(trgnr)])
			<< setw(20) << Form("% 3.1f GeV", trg_effmap_yb0ys0_alt[triggers.at(trgnr)]) << endl;
		}
		trg_turnon_file << "-----------------------------------------------------------------------------" << endl << endl;

		//announce the final binning also in this info-file (same as above in output):
		trg_turnon_file << "Final binning (after trigger studies with finer binning!):" << endl;
		trg_turnon_file << "[";
		for(int i=0; i!=nbb-1; ++i){
			trg_turnon_file << Form("%d, ", finalbins.at(i));
		}
		trg_turnon_file << Form("%d]", finalbins.at(nbb-1)) << endl;
		trg_turnon_file.close();

	}//end IF DATA ( == if(montecarlo==false))

}//end function triggerpaths()


// triggerpaths() function without ybys binning
// --> removed as not used.



//------------------------------------------//
//		function for unfolding				//
//------------------------------------------//
//takes MC Loop() output file dijet_interim_output_mc.root
//takes data triggerpaths() output file trg_lumi_output_data_Run[ABCD].root
//does unfolding using the response matrix in the mc-input file
void DijetAnalysis::unfold(){
	TH1::SetDefaultSumw2(kTRUE);		//in order to ALWAYS store the bin errors
	setTDRStyle();

	//open input file (histo_file after trigger-studies) and MC-file
	//usually this should be a data-file, but for cross-checks of unfolding procedure can as well be a MC file
	//here don't catch the case of "no runperiod given"... simply assume there is one.
	const char* inputfile_name;
	const char* mcfile_name;
	if(montecarlo==false){//if data
		//check which runperiod
		cout << "Run " << runperiod << endl;
		inputfile_name = Form("trg_lumi_output_data_Run%s.root", runperiod);
	}
	else{ 
		cout << "MC file." << endl;
		if(strcmp(generator, "herwig")==0){
			inputfile_name = "dijet_interim_output_mc_herwig7.root"; //if MC
		}
		else if(strcmp(generator, "pythia")==0){
			inputfile_name = "dijet_interim_output_mc_pythia8.root"; //if MC
		}
		else if(strcmp(generator, "HTpythia")==0){
			inputfile_name = "dijet_interim_output_mc_HT_madgraph_pythia8.root"; //if MC
		}
		else{
			inputfile_name = "dijet_interim_output_mc.root"; //if MC
		}
	}

	//MC input file --> where the response matrix for the unfolding is stored
	if(strcmp(generator, "herwig")==0){
		mcfile_name = "dijet_interim_output_mc_herwig7.root"; //if MC
	}
	else if(strcmp(generator, "pythia")==0){
		mcfile_name = "dijet_interim_output_mc_pythia8.root"; //if MC
	}
	else if(strcmp(generator, "HTpythia")==0){
		mcfile_name = "dijet_interim_output_mc_HT_madgraph_pythia8.root"; //if MC
	}

	else{
		mcfile_name = "dijet_interim_output_mc.root"; //if MC
	}

	//open input file
	TFile *inputfile = TFile::Open(inputfile_name);	//contains input
	TFile *mcfile = TFile::Open(mcfile_name);		//contains response matrix

	//Create output file
	TFile* unfolding_output;
	if(montecarlo==false){//if data
		unfolding_output = new TFile(Form("unfolding_output_data_Run%s.root", runperiod), "RECREATE");
	}
	else{ //if MC
		if(strcmp(generator, "herwig")==0){
			unfolding_output = new TFile("unfolding_output_mc_herwig7.root", "RECREATE");
		}
		else if(strcmp(generator, "pythia")==0){
			unfolding_output = new TFile("unfolding_output_mc_pythia8.root", "RECREATE");
		}
		else if(strcmp(generator, "HTpythia")==0){
			unfolding_output = new TFile("unfolding_output_mc_HT_madgraph_pythia8.root", "RECREATE");
		}
		else{//if no generator specified
			unfolding_output = new TFile("unfolding_output_mc.root", "RECREATE");
		}
	}


	//create subfolder "Standard" in outputfile
	unfolding_output->cd();
	unfolding_output->mkdir("Standard"); //create directory in output file
	TDirectory *standard_folder = unfolding_output->GetDirectory("Standard"); 
	gROOT->GetListOfBrowsables()->Add(standard_folder, "Standard");
	standard_folder->cd();

	//getting the response matrix
	//so far just testing with central bin --> change to loop through ybys bins
	///inputfile->cd("Standard/yb0ys0")
	TH2D* resp;
	mcfile->GetObject("Standard/yb0ys0/mc_gen_reco_ptavg_yb0ys0", resp);

	//Create TUnfold instance
	TUnfold *unf = new TUnfold(resp, TUnfold::kHistMapOutputHoriz);

	//getting input histogram
	TH1D* inp;
	inputfile->GetObject("Standard/yb0ys0/alternative_trigger_fits/sumhisto_yb0ys0", inp);

	//setting the input
	unf->SetInput(inp);

	//do the unfolding
	unf->DoUnfold(0);

	//get output
	TH1D* hist_data_result;	//--> clone from input data?
	TH2D* hist_data_rho;	//??
	//histogram with unfolded data:
	unf->GetOutput(hist_data_result);
	//TH1D* hist_data_result = unf->GetOutput("Unfolded");

	//histogram (2D) of correlatin coefficient
	unf->GetRhoIJ(hist_data_rho);


}//end function unfold()




//function for normalising the MC histograms
void DijetAnalysis::normalise_mc(){
	TH1::SetDefaultSumw2(kTRUE);		//in order to ALWAYS store the bin errors
	setTDRStyle();

	fChain->SetBranchStatus("Generator_weight", 1);	//generator weight
	Long64_t nentries = fChain->GetEntriesFast();	//get no. of entries

	//vectors containing the style for each ybys bin (as later used in the RatioAnalysis.C script)
    vector<Color_t> cols = {kYellow+1, kGray+2, kGreen+2, kBlue-2, kOrange-3, kRed+2, kAzure+1, kViolet+6, kSpring+5, kViolet-5, kCyan-3};
    vector<Style_t> markers = {kOpenCircle, kOpenTriangleUp, kOpenSquare, kOpenDiamond, kOpenCross, kOpenTriangleDown};


	//open input MC file and create output file
	const char* inputfile_name;
	const char* mc_inputhist;				//name of the histogram to be normalised
	TFile* normalised_mc_output;			//dijet depends on generator

	//prepare normalisation factors
	Int_t nevents_tot;		//total number of generated MC events
	Double_t xs_mc_sample;	//cross section of the mc sample used (in pb)

	cout << endl;
	//Actions according to generator that has been used
	//-------------------------------------------------
	//choose input file containing the histograms
	//create output file with proper name
	//set the normalisation factor
	if(strcmp(generator, "herwig")==0){
		inputfile_name = "dijet_interim_output_mc_herwig7.root";
		mc_inputhist = "mc_histo_smearedJet";					//if herwig --> take smeared Jet

		normalised_mc_output = new TFile("dijet_normalised_output_mc_herwig7.root", "RECREATE");
		//nevents_tot = 19464000;		//now assigned via mc_genweight histo
		xs_mc_sample = 1370000000.0;	//currently taken from the pythia flat sample, as not found for herwig sample

		cout << Form("[DijetAnalysis::normalise_mc()]: Looking at a Herwig7 sample.") << endl;
		cout << "[DijetAnalysis::normalise_mc()]: Cross Section of MC sample (used from pythia8 flat sample): " << xs_mc_sample << " pb"<< endl;
	}
	else if(strcmp(generator, "pythia")==0){
		inputfile_name = "dijet_interim_output_mc_pythia8.root";
		mc_inputhist = "mc_histo";

		normalised_mc_output = new TFile("dijet_normalised_output_mc_pythia8.root", "RECREATE");
		//nevents_tot = 19708000;		//now assigned via mc_genweight histo
		xs_mc_sample = 1370000000.0; 	//in pb? (von XSDB)

		cout << "[DijetAnalysis::normalise_mc()]: Looking at a Pythia8 sample." << endl;
		cout << "[DijetAnalysis::normalise_mc()]: Cross Section of this MC sample: " << xs_mc_sample << " pb" << endl;
	}
	else if(strcmp(generator, "HTpythia")==0){
		inputfile_name = "dijet_interim_output_mc_HT_madgraph_pythia8.root";
		mc_inputhist = "mc_histo";

		normalised_mc_output = new TFile("dijet_normalised_output_mc_HT_madgraph_pythia8.root", "RECREATE");
		xs_mc_sample = 1370000000.0; 	//in pb? (von XSDB, taken from the other pythia8 sample)

		cout << "[DijetAnalysis::normalise_mc()]: Looking at a Pythia8 madgraph sample with HT binning." << endl;
		cout << "[DijetAnalysis::normalise_mc()]: Cross Section of this MC sample: " << xs_mc_sample << " pb" << endl;
	}
	else{
		inputfile_name = "dijet_interim_output_mc.root";
		normalised_mc_output = new TFile("dijet_normalised_output_mc.root", "RECREATE");
		cout << "ERROR: No generator specified in DijetAnalysis.h !!! -- Aborting." << endl;
		exit(EXIT_FAILURE);
	}

	cout << endl;

	//open input file
	TFile *inputfile = TFile::Open(inputfile_name);	//contains input = MC histogram (not yet normalised)

	//create subfolder "Standard" in outputfile
	normalised_mc_output->cd();					//change to output file
	normalised_mc_output->mkdir("Standard"); 	//create directory in output file
	TDirectory *standard_folder = normalised_mc_output->GetDirectory("Standard"); 
	gROOT->GetListOfBrowsables()->Add(standard_folder, "Standard");
	standard_folder->cd();


	//get the genweight histogram
	TH1D* mc_genweight;
	inputfile->GetObject(Form("Standard/mc_genweight"), mc_genweight);		//histogram with genweights for original mc set (before cuts and selection)
	cout << "[DijetAnalysis::normalise_mc()]: Total number of generated events: " << mc_genweight->GetEntries() << endl;


	//create the "ybXysX" folders in outputfile
	YbYs *ybys_pointer;	//contains way more than needed here
	//new outermost bin bounds --> go only up to 2.4 for both yb, ys
	ybys_pointer = new YbYs[6]{{0.0,0.0,1.0,1.0}, {0.0, 1.0, 1.0, 2.0}, {0.0, 2.0, 1.0, 2.4}, {1.0, 0.0, 2.0, 1.0}, {1.0, 1.0, 2.0, 2.0}, {2.0, 0.0, 2.4, 1.0}}; //new!!!
	for (int ybysnr=0; ybysnr!=6; ++ybysnr){
		YbYs cur_ybysbin = ybys_pointer[ybysnr];
        const char* ybys_bname = cur_ybysbin.get_bname().c_str();				//get current ybys bin name
		TDirectory* ybys_Dir = cur_ybysbin.ybys_directories(standard_folder);	//create ybys-directory
		///cur_ybysbin.ybys_directories(standard_folder);	//create ybys-directory


		//Get the all the MC histos for dijet
		//will then additionally add correctly normalised ones
		TH1D* cur_recohisto;	//histo containing reco distribution before normalising
		TH1D* cur_genhisto;		//histo containing gen distribution (same selection as on reco, using GenJet instead of Jet Branches)
		TH2D* cur_response_noGenCuts;	//response matrix (reco, gen in ptavg) without any gen cuts applied
		TH2D* cur_norm_resp_noGenCuts;	//response matrix (no gen cuts) normalised to number of generated jets
		TH2D* cur_response;				//response matrix, including same cuts on reco and gen
		TH2D* cur_norm_resp;			//normalised response matrix

		//get the objects:
		////inputfile->GetObject(Form("Standard/%s/mc_histo_%s", ybys_bname, ybys_bname), cur_recohisto);						//reco histo
		inputfile->GetObject(Form("Standard/%s/%s_%s", ybys_bname, mc_inputhist, ybys_bname), cur_recohisto);						//reco histo (mc_histo_%s for pythia, mc_histo_smearedJet_%s for herwig)

		inputfile->GetObject(Form("Standard/%s/mc_genhisto_%s", ybys_bname, ybys_bname), cur_genhisto);						//gen histo
		inputfile->GetObject(Form("Standard/%s/mc_genrec_noGenCuts_%s", ybys_bname, ybys_bname), cur_response_noGenCuts);	//resp matrix w/o gen cuts
		inputfile->GetObject(Form("Standard/%s/normalised_response_noGenCuts", ybys_bname), cur_norm_resp_noGenCuts);		//normalised resp matrix w/o gen cuts
		inputfile->GetObject(Form("Standard/%s/mc_gen_reco_ptavg_%s", ybys_bname, ybys_bname), cur_response);				//response matrix with gen and reco cuts
		inputfile->GetObject(Form("Standard/%s/normalised_response", ybys_bname), cur_norm_resp);							//response matrix with gen and reco cuts

		//append the ybys bin name to the normalised response matrices name:
		cur_norm_resp_noGenCuts->SetName(Form("%s_%s", cur_norm_resp_noGenCuts->GetName(), ybys_bname));
		cur_norm_resp->SetName(Form("%s_%s", cur_norm_resp->GetName(), ybys_bname));

		//Clone the histograms that shall be normalised
		TH1D* norm_recohisto = (TH1D*)cur_recohisto->Clone();
		TH1D* norm_genhisto = (TH1D*)cur_genhisto->Clone();

		//change name of the unnormalised histograms:
		cur_recohisto->SetName(Form("%s_NotNormalised", cur_recohisto->GetName()));
		cur_genhisto->SetName(Form("%s_NotNormalised", cur_genhisto->GetName()));

		//Normalise the histograms (factor is already set according to generator):
		//NEW: get the sum of weights directly from the inputfile!
		Double_t weights_mean = mc_genweight->GetMean();
		nevents_tot	= mc_genweight->GetEntries();
		Double_t sum_event_weights = nevents_tot*weights_mean;

		//reco histo
		///norm_recohisto->Scale(1./nevents_tot);		//divide by total number of generated events
		norm_recohisto->Scale(1./sum_event_weights);
		norm_recohisto->Scale(xs_mc_sample);		//multiply by cross section of this mc sample
		norm_recohisto->Scale(1000);				//to get end value in fb (--> would pb be more suitable?)
		//norm_recohisto->Scale(58.83);				//this is the total luminosity of 2018 data in /fb ...not needed here.

		//same for gen histo
		///norm_genhisto->Scale(1./nevents_tot);
		norm_genhisto->Scale(1./sum_event_weights);
		norm_genhisto->Scale(xs_mc_sample);
		norm_genhisto->Scale(1000);
		//norm_genhisto->Scale(58.83);

		cout << Form("[DijetAnalysis::normalise_mc()]: Normalised histograms in %s bin", ybys_bname) << endl;
		cout << "[DijetAnalysis::normalise_mc()]: Storing to file..." << endl;
		//store all the histograms to the outputfile
		ybys_Dir->cd();

		SetDrawStyleWrite(norm_recohisto, "", norm_recohisto->GetTitle(), "p_{T,avg} /GeV", "XS in fb", cols.at(ybysnr+2), 2, markers.at(ybysnr), cols.at(ybysnr+2), 2, true, 0.9, true);
		norm_recohisto->GetYaxis()->SetTitleOffset(0.9);
		//norm_recohisto->Write();
		norm_genhisto->Write();
		cur_recohisto->Write();
		cur_genhisto->Write();
		cur_response_noGenCuts->Write();
		cur_norm_resp_noGenCuts->Write();
		cur_response->Write();
		cur_norm_resp->Write();


	 	//vectors containing the style for each ybys bin
    	//vector<Color_t> cols = {kYellow+1, kGray+2, kGreen+2, kBlue-2, kOrange-3, kRed+2, kAzure+1, kViolet+6, kSpring+5, kViolet-5, kCyan-3};
    	//vector<Style_t> markers = {kOpenCircle, kOpenTriangleUp, kOpenSquare, kOpenDiamond, kOpenCross, kOpenTriangleDown};

		///void SetDrawStyleWrite(TH1D *object, const char *drawcmd="", const char *title="", const char *xtitle="p_{T,avg}", const char *ytitle="", Color_t linecolor=kGray+2, Int_t linewidth=2, Style_t markerstyle=kFullCircle, Color_t markercolor=kGray+2, Int_t markersize=1, bool morelog=false, Double_t xTitleOffset=1.2, bool writing=false){}
	
	}

	cout << endl;
	cout << "[DijetAnalysis::normalise_mc()]: Results can be found in file: " << normalised_mc_output->GetName() << endl;



}//end function normalise_mc()


//function to divide uncUp by uncDown for each uncertainty source (taking the input file coming from the Loop function)
void DijetAnalysis::uncs_UpDownDivide_mc(){
	TH1::SetDefaultSumw2(kTRUE);		//in order to ALWAYS store the bin errors

	
	//recreate the uncnames vector --> could also do that just once as a global variable?
	vector<string> uncnames;
	uncnames.push_back("AbsoluteStat");
	uncnames.push_back("AbsoluteScale");
	uncnames.push_back("AbsoluteSample");
	uncnames.push_back("AbsoluteFlavMap");
	uncnames.push_back("AbsoluteMPFBias");
	uncnames.push_back("Fragmentation");
	uncnames.push_back("SinglePionECAL");
	uncnames.push_back("SinglePionHCAL");
	uncnames.push_back("FlavorQCD");
	uncnames.push_back("TimePtEta");
	uncnames.push_back("RelativeJEREC1");
	uncnames.push_back("RelativeJEREC2");
	uncnames.push_back("RelativeJERHF");
	uncnames.push_back("RelativePtBB");
	uncnames.push_back("RelativePtEC1");
	uncnames.push_back("RelativePtEC2");
	uncnames.push_back("RelativePtHF");
	uncnames.push_back("RelativeBal");
	uncnames.push_back("RelativeSample");
	uncnames.push_back("RelativeFSR");
	uncnames.push_back("RelativeStatFSR");
	uncnames.push_back("RelativeStatEC");
	uncnames.push_back("RelativeStatHF");
	uncnames.push_back("PileUpDataMC");
	uncnames.push_back("PileUpPtRef");
	uncnames.push_back("PileUpPtBB");
	uncnames.push_back("PileUpPtEC1");
	uncnames.push_back("PileUpPtEC2");
	uncnames.push_back("PileUpPtHF");
	uncnames.push_back("PileUpMuZero");
	uncnames.push_back("PileUpEnvelope");
	uncnames.push_back("SubTotalPileUp");
	uncnames.push_back("SubTotalRelative");
	uncnames.push_back("SubTotalPt");
	uncnames.push_back("SubTotalScale");
	uncnames.push_back("SubTotalAbsolute");
	uncnames.push_back("SubTotalMC");
	uncnames.push_back("Total");
	uncnames.push_back("TotalNoFlavor");
	uncnames.push_back("TotalNoTime");
	uncnames.push_back("TotalNoFlavorNoTime");
	uncnames.push_back("FlavorZJet");
	uncnames.push_back("FlavorPhotonJet");
	uncnames.push_back("FlavorPureGluon");
	uncnames.push_back("FlavorPureQuark");
	uncnames.push_back("FlavorPureCharm");
	uncnames.push_back("FlavorPureBottom");
	uncnames.push_back("TimeRunA");
	uncnames.push_back("TimeRunB");
	uncnames.push_back("TimeRunC");
	uncnames.push_back("TimeRunD");
	uncnames.push_back("CorrelationGroupMPFInSitu");
	uncnames.push_back("CorrelationGroupIntercalibration");
	uncnames.push_back("CorrelationGroupbJES");
	uncnames.push_back("CorrelationGroupFlavor");
	uncnames.push_back("CorrelationGroupUncorrelated");


	//take care of input and output file
	const char* inputfile_name = "dijet_interim_output_mc_herwig7.root";
	TFile *inputfile = TFile::Open(inputfile_name);												//open inputfile
	TFile *divided_uncs_mc_output = new TFile("dijet_divided_uncsources_mc_herwig7.root", "RECREATE");	//create output file

	cout << "[DijetAnalysis::uncs_UpDownDivide_mc()]: Getting the Up Down uncsources from inputfile: " << inputfile->GetName() << endl;
	cout << endl;

	//create subfolder "Standard" in outputfile
	divided_uncs_mc_output->cd();					//change to output file
	divided_uncs_mc_output->mkdir("UncertaintyRatios"); 		//create directory in output file
	TDirectory *standard_folder = divided_uncs_mc_output->GetDirectory("UncertaintyRatios"); 
	gROOT->GetListOfBrowsables()->Add(standard_folder, "UncertaintyRatios");
	standard_folder->cd();

	
	//create the "ybXysX" folders in outputfile
	YbYs *ybys_pointer;	//contains way more than needed here
	//new outermost bin bounds --> go only up to 2.4 for both yb, ys
	ybys_pointer = new YbYs[6]{{0.0,0.0,1.0,1.0}, {0.0, 1.0, 1.0, 2.0}, {0.0, 2.0, 1.0, 2.4}, {1.0, 0.0, 2.0, 1.0}, {1.0, 1.0, 2.0, 2.0}, {2.0, 0.0, 2.4, 1.0}}; //new!!!
	for (int ybysnr=0; ybysnr!=6; ++ybysnr){
		YbYs cur_ybysbin = ybys_pointer[ybysnr];
        const char* ybys_bname = cur_ybysbin.get_bname().c_str();				//get current ybys bin name
		TDirectory* ybys_Dir = cur_ybysbin.ybys_directories(standard_folder);	//create ybys-directory
		///cur_ybysbin.ybys_directories(standard_folder);	//create ybys-directory


		//loop through the uncertainty sources for creating the unc ratios
		for(unsigned int k=0; k!=uncnames.size(); ++k){
			cout << Form("[DijetAnalysis::uncs_UpDownDivide_mc()]: Looking at unc source: %s", uncnames.at(k).c_str()) << endl;

			TH1D* cur_uncUphist;	//histo containing the Up variation of current unc source
			TH1D* cur_uncDownhist;	//histo containing the Down variation of current unc source
			//--> currently not copied to new file --> should this be imlemented?

			//get the objects:
			inputfile->GetObject(Form("Uncertainties/%s/unc_%s_Up_%s", ybys_bname, uncnames.at(k).c_str(), ybys_bname), cur_uncUphist);		//up
			inputfile->GetObject(Form("Uncertainties/%s/unc_%s_Down_%s", ybys_bname, uncnames.at(k).c_str(), ybys_bname), cur_uncDownhist);	//down

			//do not forget to store the uncs Sumw2()
			cur_uncUphist->Sumw2();
			cur_uncDownhist->Sumw2();

			//Clone the histograms (unc up) that will be divided by the other (unc down)
			TH1D* cur_unc_UpDownDiv = (TH1D*)cur_uncUphist->Clone(Form("UpDownDiv_%s_%s", uncnames.at(k).c_str(), ybys_bname));
			cur_unc_UpDownDiv->Sumw2();
			

			//Do the division:
			cur_unc_UpDownDiv->Divide(cur_uncDownhist);
		
			cout << "[DijetAnalysis::uncs_UpDownDivide_mc()]: Storing to file..." << endl;
			//store all the histograms to the outputfile
			ybys_Dir->cd();
			//SetDrawStyleWrite(norm_recohisto, "", norm_recohisto->GetTitle(), "p_{T,avg} /GeV", "XS in fb", cols.at(ybysnr+2), 2, markers.at(ybysnr), cols.at(ybysnr+2), 2, true, 0.9, true);
			//norm_recohisto->GetYaxis()->SetTitleOffset(0.9);
			//norm_recohisto->Write();
			//norm_genhisto->Write();
			cur_unc_UpDownDiv->GetYaxis()->SetTitleOffset(0.9);
			cur_unc_UpDownDiv->SetTitle(Form("%s Up/Down", uncnames.at(k).c_str()));
			cur_unc_UpDownDiv->GetYaxis()->SetTitle("Ratio");
			cur_unc_UpDownDiv->Write();

		}//end of unc sources for-loop
		cout << Form("[DijetAnalysis::uncs_UpDownDivide_mc()]: Divided unc histograms in %s bin", ybys_bname) << endl;
		cout << "---------------------------------------------------------------------------------" << endl;
		cout << endl;

	}//end of ybys loop

	//is the following necessary? --No.
	//divided_uncs_mc_output->Write();

	cout << "[DijetAnalysis::uncs_UpDownDivide_mc()]: Results can be found in file: " << divided_uncs_mc_output->GetName() << endl;

}//end of function uncs_UpDownDivide_mc()


//Function to divide each bin (in XS histograms) by binwidth, in order to end up with fb/GeV as XS unit
void DijetAnalysis::DivideBinwidth(){
	TH1::SetDefaultSumw2(kTRUE);		//in order to ALWAYS store the bin errors

	//styleparams
    vector<Color_t> cols = {kYellow+1, kGray+2, kGreen+2, kBlue-2, kOrange-3, kRed+2, kAzure+1, kViolet+6, kSpring+5, kViolet-5, kCyan-3};
    vector<Style_t> markers = {kOpenCircle, kOpenTriangleUp, kOpenSquare, kOpenDiamond, kOpenCross, kOpenTriangleDown};


	//choose the file to be used
	const char* histname;
	const char* mc_inputhist;		//only used if mc (because of different naming in pythia(not smeared) and herwig(smeared))
	const char* inputfile_name;
	TFile* binwidth_output;
	if(montecarlo==false){//if data
		//check which runperiod
		//prepare output file (= after analysis BEFORE trigger analysis or creating plots)
		if(strcmp(runperiod,"A")==0){
			cout << "Run A" << endl;
			///inputfile_name = "trg_lumi_output_data_RunA.root";			//RunA -- these files have a binning problem!
			inputfile_name = "dijet_fastloop_output_data_RunA.root";		//RunA	 -- created with FastLoop()
			binwidth_output = new TFile("dijet_xs_BinwidthDivided_data_RunA.root", "RECREATE");
			histname = "sumhisto_dijet_RunA";
		}//endif runA

		else if(strcmp(runperiod,"B")==0){
			cout << "Run B" << endl;
			//inputfile_name = "trg_lumi_output_data_RunB.root";			//RunB
			inputfile_name = "dijet_fastloop_output_data_RunB.root";		//RunB	 -- created with FastLoop()
			binwidth_output = new TFile("dijet_xs_BinwidthDivided_data_RunB.root", "RECREATE");
			histname = "sumhisto_dijet_RunB";
		}//endif runB

		else if(strcmp(runperiod,"C")==0){
			cout << "Run C" << endl;
			///inputfile_name = "trg_lumi_output_data_RunC.root";			//RunC
			inputfile_name = "dijet_fastloop_output_data_RunC.root";		//RunC	 -- created with FastLoop()
			binwidth_output = new TFile("dijet_xs_BinwidthDivided_data_RunC.root", "RECREATE");
			histname = "sumhisto_dijet_RunC";
		}//endif runC

		else if(strcmp(runperiod,"D")==0){
			cout << "Run D" << endl;
			///inputfile_name = "trg_lumi_output_data_RunD.root";			//RunD
			inputfile_name = "dijet_fastloop_output_data_RunD.root";		//RunD	 -- created with FastLoop()
			binwidth_output = new TFile("dijet_xs_BinwidthDivided_data_RunD.root", "RECREATE");
			histname = "sumhisto_dijet_RunD";
		}//endif runD

		else{//if no runperiod specified
			cout << "ERROR: If this is data set -- specify run period." << endl;
			exit(EXIT_FAILURE);
		}
	}
	else{ //if MC
		if(strcmp(generator, "herwig")==0){
			inputfile_name = "dijet_normalised_output_mc_herwig7.root";	//herwig7
			binwidth_output = new TFile("dijet_xs_BinwidthDivided_mc_herwig7.root", "RECREATE"); //output file if herwig7
			mc_inputhist = "mc_histo_smearedJet";
			histname = "sumhisto_dijet_herwig7";
		}
		else if(strcmp(generator, "pythia")==0){
			inputfile_name = "dijet_normalised_output_mc_pythia8.root";	//pythia8
			binwidth_output = new TFile("dijet_xs_BinwidthDivided_mc_pythia8.root", "RECREATE"); //output file if pythia8
			mc_inputhist = "mc_histo";
			histname = "sumhisto_dijet_pythia8";
		}
		else if(strcmp(generator, "HTpythia")==0){
			//inputfile_name = "dijet_normalised_output_mc_HT_madgraph_pythia8.root";	//pythia8 --> this normalisation is already done in Loop() for each HT bin
			inputfile_name = "dijet_interim_output_mc_HT_madgraph_pythia8.root";	//--> use directly the Loop() output
			binwidth_output = new TFile("dijet_xs_BinwidthDivided_mc_HT_madgraph_pythia8.root", "RECREATE"); //output file if pythia8 with madgraph in HT bins
			mc_inputhist = "mc_histo";
			histname = "sumhisto_dijet_HT_madgraph_pythia8";
		}
		else{//if no generator specified
			cout << "ERROR: If this is mc set -- specify generator." << endl;
			exit(EXIT_FAILURE);
		}
	}

	TFile *inputfile = TFile::Open(inputfile_name);												//open inputfile

	cout << "[DijetAnalysis::DivideBinwidth()]: Getting the XS histograms from inputfile: " << inputfile->GetName() << endl;
	cout << endl;

	//create subfolder "Standard" in outputfile
	binwidth_output->cd();					//change to output file
	binwidth_output->mkdir("CrossSections_byBinwidth"); 		//create directory in output file
	TDirectory *standard_folder = binwidth_output->GetDirectory("CrossSections_byBinwidth"); 
	gROOT->GetListOfBrowsables()->Add(standard_folder, "CrossSections_byBinwidth");
	standard_folder->cd();

	
	//create the "ybXysX" folders in outputfile
	YbYs *ybys_pointer;	//contains way more than needed here
	//new outermost bin bounds --> go only up to 2.4 for both yb, ys
	ybys_pointer = new YbYs[6]{{0.0,0.0,1.0,1.0}, {0.0, 1.0, 1.0, 2.0}, {0.0, 2.0, 1.0, 2.4}, {1.0, 0.0, 2.0, 1.0}, {1.0, 1.0, 2.0, 2.0}, {2.0, 0.0, 2.4, 1.0}}; //new!!!
	for (int ybysnr=0; ybysnr!=6; ++ybysnr){
		YbYs cur_ybysbin = ybys_pointer[ybysnr];
        const char* ybys_bname = cur_ybysbin.get_bname().c_str();				//get current ybys bin name
		TDirectory* ybys_Dir = cur_ybysbin.ybys_mkdir(standard_folder);			//create ybys-directory
		///cur_ybysbin.ybys_directories(standard_folder);	//create ybys-directory

		TH1D* cur_xshisto;		//histo containing the cross section for current ybys bin

		//get the object:
		if(montecarlo==false){
			//inputfile->GetObject(Form("Standard/%s/alternative_trigger_fits/sumhisto_%s", ybys_bname, ybys_bname), cur_xshisto);		//cross section histo to be normalised (this is the trg_lumi_output histo, which has a binning problem!
			inputfile->GetObject(Form("Standard/%s/sumhisto_%s", ybys_bname, ybys_bname), cur_xshisto);		//cross section histo to be normalised
		}
		else{
			inputfile->GetObject(Form("Standard/%s/%s_%s", ybys_bname, mc_inputhist, ybys_bname), cur_xshisto);							//cross section histo to be normalised

			//for HT binned histogram --> scale by factor 1000 as it is so far in pb and want it in fb
			if(strcmp(generator, "HTpythia")==0){cur_xshisto->Scale(1000);}	//this hopefully just affects the new histo as the old file won't be written in this function? (check!)
		}


		//Clone the histogram
		TH1D* cur_xshisto_bindiv = (TH1D*)cur_xshisto->Clone(Form("%s_%s", histname, ybys_bname));
		cur_xshisto_bindiv->Sumw2();
			
		//loop through the bins of the cloned histogram and divide each bincontent by the corresponding binwidth
		int nbins = cur_xshisto_bindiv->GetNbinsX();
		for(int j=1; j!=nbins+1; ++j){		//start indexing from 1 because of ROOT (0 is underflow bin)
			cur_xshisto_bindiv->SetBinContent(j, 0);	//new histo: set content to zero
			cur_xshisto_bindiv->SetBinError(j, 0);		//new histo: set binerror to zero
			Double_t curcont = cur_xshisto->GetBinContent(j);				//original histo: get bin content
			Double_t curwidth = cur_xshisto->GetXaxis()->GetBinWidth(j);	//original histo: get binwidth
			Double_t curerr = cur_xshisto->GetBinError(j);					//original histo: get bin error

			//calculate the new bincontent and -error
			//fill new content and error to bin in new histogram
			cur_xshisto_bindiv->SetBinContent(j, curcont/curwidth);
			cur_xshisto_bindiv->SetBinError(j, curerr/curwidth);


		}//end of loop through histogram (ptavg) bins

		cout << "[DijetAnalysis::DivideBinwidth()]: Set new bin contents, divided by binwidth." << endl;
		cout << "[DijetAnalysis::DivideBinwidth()]: Done with ybysbin: " << ybys_bname << endl;

		//store the histogram to the outputfile
		ybys_Dir->cd();
		cur_xshisto_bindiv->GetYaxis()->SetTitleOffset(0.9);
		cur_xshisto_bindiv->GetYaxis()->SetTitle("XS in fb/GeV");
		cur_xshisto_bindiv->SetLineColor(cols.at(ybysnr+2));
		cur_xshisto_bindiv->SetMarkerColor(cols.at(ybysnr+2));
		cur_xshisto_bindiv->SetMarkerStyle(markers.at(ybysnr));
		if(montecarlo==false){
			cur_xshisto_bindiv->SetTitle(Form("Dijet XS Run %s", runperiod));
		}
		else{
			cur_xshisto_bindiv->SetTitle(Form("Dijet XS %s", generator));
		}

		cur_xshisto_bindiv->Write();

		cout << Form("[DijetAnalysis::DivideBinwidth()]: Writing to file.") << endl;
		cout << "---------------------------------------------------------------------------------" << endl;
		cout << endl;
	}//end of ybys loop

	cout << "[DijetAnalysis::DivideBinwidth()]: Results can be found in file: " << binwidth_output->GetName() << endl;

}//end of function DivideBinwidth()


//function for speeding up the trigger-analysis stuff by taking fixed turn-on point values
//currently only implemented for data, as there are no triggers for MC
void DijetAnalysis::FastLoop(){
	TH1::SetDefaultSumw2(kTRUE);		//in order to ALWAYS store the bin errors

	//styleparams
    vector<Color_t> cols = {kYellow+1, kGray+2, kGreen+2, kBlue-2, kOrange-3, kRed+2, kAzure+1, kViolet+6, kSpring+5, kViolet-5, kCyan-3};
    vector<Style_t> markers = {kOpenCircle, kOpenTriangleUp, kOpenSquare, kOpenDiamond, kOpenCross, kOpenTriangleDown};

	//the following is partly copied from the original Loop() function
	//insert analysis loop as it used to be in triggerpaths
	setTDRStyle();
	std::string name_h = "hpt";
	fChain->SetBranchStatus("*",0); 			//exclude / switch off all the branches
	fChain->SetBranchStatus("nJet",1); 			//include this branch in analysis	
	fChain->SetBranchStatus("Jet_pt",1); 		//include this branch as well
	fChain->SetBranchStatus("Jet_eta",1); 		//has to be included if eta-bins shall be checked...
	fChain->SetBranchStatus("Jet_phi",1);		//for TLorentzVector
	fChain->SetBranchStatus("Jet_mass",1);		//for TLorentzVector

	//for noise reduction (MET checks)
	fChain->SetBranchStatus("ChsMET_pt",1);		//MET in event
	fChain->SetBranchStatus("ChsMET_sumEt",1);	//Et sum in event
	fChain->SetBranchStatus("MET_pt",1);		//MET in event
	fChain->SetBranchStatus("MET_sumEt",1);	//Et sum in event

	cout << Form("This is an analysis of dijet data. Using data from Run period: %s", runperiod)  << endl;	//output does not work...
	cout << "Calling DijetAnalysis::FastLoop() for quick dijet data anlysis (fixed trigger turn-ons)." << endl;

	//AT THE SAME TIME SET LUMI-VECTOR ACCORDING TO RUNPERIOD!
	//new lumi (from own json files, but after skimming, used brilcalc) units are /fb
	vector<Double_t> lumivec;
	TFile* interim_output;
	if(montecarlo==false){//if data --- always has to be data for this function to run
		//check which runperiod
		//prepare output file (= after analysis BEFORE trigger analysis or creating plots)
		if(strcmp(runperiod,"A")==0){
			cout << "Run A" << endl;
			interim_output = new TFile("dijet_fastloop_output_data_RunA.root", "RECREATE");
			lumivec = {0.000000084, 0.000000084, 0.000065758, 0.000088949, 0.002453039, 0.010118712, 0.044313793, 0.115597030, 0.397308191, 0.910868208, 1.780016341, 13.977332725, 13.977332725};
		}//endif runA
		else if(strcmp(runperiod,"B")==0){
			cout << "Run B" << endl;
			interim_output = new TFile("dijet_fastloop_output_data_RunB.root", "RECREATE");
			lumivec = {0, 0, 0.000026808, 0.000110281, 0.000415176, 0.005846909, 0.025058183, 0.055140569, 0.220562275, 0.441124550, 0.882249099, 7.057992794, 7.057992794};
		}//endif runB
		else if(strcmp(runperiod,"C")==0){
			cout << "Run C" << endl;
			interim_output = new TFile("dijet_fastloop_output_data_RunC.root", "RECREATE");
			lumivec = {0.000000063, 0.000000063, 0.000027608, 0.000107731, 0.000405575, 0.005647855, 0.024204982, 0.054976951, 0.216547076, 0.431973909, 0.862827576, 6.894778911, 6.894778911};
		}//endif runC
		else if(strcmp(runperiod,"D")==0){
			cout << "Run D" << endl;
			interim_output = new TFile("dijet_fastloop_output_data_RunD.root", "RECREATE");
			lumivec = {0.000000041, 0.000000041, 0.000119557, 0.000495470, 0.001865299, 0.026760713, 0.114682242, 0.248467482, 0.991655052, 1.982571812, 3.964405332, 31.710074612, 31.710074612};
		}//endif runD
		else{//if no runperiod specified
			cout << "Error: Need to specify runperiod in DijetAnalysis.h for this function to run. Abort!" << endl;
			exit(EXIT_FAILURE);
		}
	}
	else{ //if MC
		cout << "Error: This script only runs on dijet DATA! (runperiod A, B, C or D). Abort!" << endl;
		exit(EXIT_FAILURE);
	}
	interim_output->cd();
	interim_output->mkdir("Standard"); //create directory in output file

	//create subfolder "Standard" in outputfile
	TDirectory *standard_folder = interim_output->GetDirectory("Standard"); 
	gROOT->GetListOfBrowsables()->Add(standard_folder, "Standard");
	standard_folder->cd();	//change back to Standard directory

	//could shorten this vector (starting from 40), but need to adjust trigger indexing everywhere...
	vector<int> triggers = {0, 15, 25, 40, 60, 80, 140, 200, 260, 320, 400, 450, 500, 550};

	if(montecarlo==false){
		//switching on selected TBranches corresponding to chosen triggers:
		for (int j=1; j!=triggers.size(); ++j){ //start from 1, because "0" means "without trigger"
			char *trigname = Form("HLT_PFJet%i", triggers.at(j));
			fChain->SetBranchStatus(trigname, 1);
		}
	}

	YbYs *ybys_pointer;	//should maybe handle this differently
	//call constructor YbYs(float yb_lo, float ys_lo, float yb_up, float ys_up)
	//new without unc_folder stuff
	ybys_pointer = new YbYs[6]{{0.0,0.0,1.0,1.0}, {0.0, 1.0, 1.0, 2.0}, {0.0, 2.0, 1.0, 2.4}, {1.0, 0.0, 2.0, 1.0}, {1.0, 1.0, 2.0, 2.0}, {2.0, 0.0, 2.4, 1.0}}; 

	vector <Histos*> histosvec;
	for (int ybysnr=0; ybysnr!=6; ++ybysnr){
		YbYs cur_ybysbin = ybys_pointer[ybysnr];
		histosvec.push_back(new Histos());
		Histos curtrig = *histosvec.back();
		cur_ybysbin.ybys_fastdirs(standard_folder, triggers, curtrig);	//creates folder structure
	}

	//counters
	int count_nodijet = 0; 		//count number of events that are no dijet events, i.e. nJet < 2
	int count_NoBin_ybys = 0; 	//count number of events that have no ybys bin, where they could be assigned to (i.e. yb or ys > 3.0)
	int count_BadRapidity = 0;	//count number of events where the jets do not fulfill the rapidity requirements
	int count_selected = 0;		//count number of selected events: implemented for MC


	//the actual Loop()
	if (fChain == 0) return;
	Long64_t nentries = fChain->GetEntriesFast();
	cout << "nentries = " << nentries << endl << endl;			//information

	Long64_t nbytes = 0, nb = 0;
	//starting event loop
	for (Long64_t jentry=0; jentry<nentries;jentry++) {
	////for (Long64_t jentry=0; jentry!=100000; ++jentry){} //for testing with less entries/events

		Long64_t ientry = LoadTree(jentry);
		if(v==true){cout << "current entry: " << jentry << endl;}
		///if (jentry%10==0){cout << "------------------------------------------" << endl;}//for overview
		if (jentry%100000==0){cout << "processed " << jentry << " of " << nentries << endl;}
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);   nbytes += nb;
		
		map<const int, Bool_t> tdecisions;

		tdecisions.insert(make_pair(triggers.at(0), true)); //-->no trigger
		tdecisions.insert(make_pair(triggers.at(1), HLT_PFJet15)); //could this be handled in a more general way? (see attempt and problem above)
		tdecisions.insert(make_pair(triggers.at(2), HLT_PFJet25));
		tdecisions.insert(make_pair(triggers.at(3), HLT_PFJet40));
		tdecisions.insert(make_pair(triggers.at(4), HLT_PFJet60));
		tdecisions.insert(make_pair(triggers.at(5), HLT_PFJet80));
		tdecisions.insert(make_pair(triggers.at(6), HLT_PFJet140));
		tdecisions.insert(make_pair(triggers.at(7), HLT_PFJet200));
		tdecisions.insert(make_pair(triggers.at(8), HLT_PFJet260));
		tdecisions.insert(make_pair(triggers.at(9), HLT_PFJet320));
		tdecisions.insert(make_pair(triggers.at(10), HLT_PFJet400));
		tdecisions.insert(make_pair(triggers.at(11), HLT_PFJet450));
		tdecisions.insert(make_pair(triggers.at(12), HLT_PFJet500));
		tdecisions.insert(make_pair(triggers.at(13), HLT_PFJet550));


		//check if Event (jentry) is a dijet event. --> otherwise do nothing
		if(nJet>1){	//required for dijet ! --> only look at two leading jets

			//add MET check for reducing noise
			//bool good_met = ChsMET_pt<(0.3*ChsMET_sumEt);		 	//Chs variable not existent in MC sets
			bool good_met = MET_pt<(0.3*MET_sumEt);			//exists in both data and mc

			//check pthat in case of MC (set it to true for data)
			bool good_pthat;
			if(montecarlo){//do this always, for every event
				if(Jet_pt[0]<2.0*Generator_binvar){good_pthat=true;}	//looser than 1.5 criterion
				else{good_pthat=false;}
			}
			else{//if data
				good_pthat=true;
			}


			Double_t j1_pt	= Jet_pt[0];
			Double_t j2_pt	= Jet_pt[1];
			Double_t j1_eta	= Jet_eta[0];
			Double_t j2_eta	= Jet_eta[1];
			Double_t j1_phi	= Jet_phi[0];
			Double_t j2_phi	= Jet_phi[1];
			Double_t j1_mass = Jet_mass[0];
			Double_t j2_mass = Jet_mass[1];
			//use Lorentzvector:
			PtEtaPhiMVector jet1_vec(j1_pt, j1_eta, j1_phi, j1_mass);
			PtEtaPhiMVector jet2_vec(j2_pt, j2_eta, j2_phi, j2_mass);

			//calculate current ptavg, yb, ys
			Double_t j1_y = jet1_vec.Rapidity();
			Double_t j2_y = jet2_vec.Rapidity();
			Double_t cur_ptavg = 0.5*(Jet_pt[0]+Jet_pt[1]);
			Double_t cur_yboost = 0.5*TMath::Abs(j1_y+j2_y);
			Double_t cur_ystar = 0.5*TMath::Abs(j1_y-j2_y);

			//add offline cut on Jet-Rapidity --> both jets must have Abs(y)<5.0 and at least one shall have Abs(y)<2.5
			Double_t absYj1 = TMath::Abs(j1_y);
			Double_t absYj2 = TMath::Abs(j2_y);
			//bool good_rapidity = (absYj1<5)&&(absYj2<5)&&(absYj1<2.5 || absYj2<2.5);	//old cuts
			bool good_rapidity = (absYj1<2.4)&&(absYj2<2.4);	//NEW cuts!

			//add offline cut on Jet-pT --> both jets must fulfill pT > 30 GeV.
			//bool good_pt = (j1_pt>30)&&(j2_pt>30);
			bool good_pt = (j1_pt>20)&&(j2_pt>20);		//set down to 20 pt; check strict criterion later

			//check if rapidity requirement is fulfilled
			//check at the same time whether pT requirement is fulfilled
			//generator weight
			if(good_rapidity&&good_pt&&good_met&&good_pthat){	//now two new extra checks: good_met and good_pthat
			//if(good_rapidity&&good_pt){}
				//check ybys bin --> put this into function!!
				for (int ybysnr=0; ybysnr!=6; ++ybysnr){
					YbYs cur_ybysbin = ybys_pointer[ybysnr];
					float yb_lo = cur_ybysbin.yb_get_bb().first;
					float yb_up = cur_ybysbin.yb_get_bb().second;
					float ys_lo = cur_ybysbin.ys_get_bb().first;
					float ys_up = cur_ybysbin.ys_get_bb().second;
					if((yb_lo<=cur_yboost && cur_yboost<yb_up) and (ys_lo<=cur_ystar && cur_ystar<ys_up)){	//-->better to implement function InInterval(value, low, up) that checks if value in interval... return(lo<=x & x<up)
						if(montecarlo==false){//only for data:
							//check the pT > 30GeV cut:
							if((j1_pt>30)&&(j2_pt>30)){
								//check trigger
								for(int trignr=0; trignr!=triggers.size(); ++trignr){
									if(tdecisions[triggers.at(trignr)]){	//==true
										//fill into histo of current trigger
										histosvec.at(ybysnr)->jt_histos[triggers.at(trignr)]->Fill(cur_ptavg);	//fill ptavg trigger histogram
										if(trignr==0){//count only once.
											count_selected++; //count this event (in case trigger 40 or higher)
										}
										continue;
									}//trigger true?
									else{//trigger false
										continue;
									}//trigger false
								}//check triggers
							}//endif both pt>30	 (is else --> continue needed?)
						}//end if DATA (if mc==false)
						else if(montecarlo==true){
							cout << "ERROR: this is apparantly MC set ! Should not even get into this place of loop!!! Abort!!!" << endl;
							exit(EXIT_FAILURE);
						}//end: if MC
						break; //break out of ybys-bin-iteration, because event has now been assigned to its corresponding ybys bin (don't need to check further ybys bins, it is unique)
					}//if dijet-event in ybys bin
					else{	//if no ybys bin was suitable
						if(ybysnr==5){
							if(v==true){cout << "last ybys bin ---> no suitable ybys bin found" << endl;}
							count_NoBin_ybys++;	//no suitable ybys bin found
							break;
						}
						continue;
					}//no suitable ybys bin
				}//loop through ybys-bins
			}//end: if good_rapidity && good_pt
			else{
				count_BadRapidity++;	//this dijet-event does not fulfill rapidity requirements
			}
		}//if (nJet>1)
		else{//if nJet<2
			count_nodijet++;	//this is not a dijet event
		}

	}//end of event-loop



	//now fill the summary histogram for each ybys bin, based on the given trigger-turnon map (taken from yb0ys0 results)
	//loop again through ybys bins
	for (int ybysnr=0; ybysnr!=6; ++ybysnr){
		YbYs cur_ybysbin = ybys_pointer[ybysnr];
		TH1D* cur_sumhist = cur_ybysbin.fill_fastsumhisto(triggers, histosvec.at(ybysnr), lumivec);	//also give it lumis for normalisation
		cur_sumhist->Write();
	}//end ybys-loop
	

	interim_output->Write();

	//Quick summary:
	cout << "Summary after filling the histograms: " << endl;
	cout << "--------------------------------------" << endl;
	cout << "Number of events: " << nentries << endl;
	cout << "Dropped out as nJet<2: " << count_nodijet << endl;
	cout << "Dropped out due to j1, j2 bad rapidity or bad pt: " << count_BadRapidity << endl;
	cout << "Dropped out as no ybys bin suitable: " << count_NoBin_ybys << endl;
	cout << "-----------------------------------------------" << endl;
	cout << "Selected Events: " << count_selected << endl << endl;


	//adjust this to mc / data case:
	cout << "Saved results to file: " << interim_output->GetName() << endl;

}//end of function FastLoop()


/*

//////// trying out code for uncertainty analysis ////////////
// -- 	Later to be implemented in Loop() function 		--	//
//for each event: need to know pt and eta of jets 1 and 2
//in case of ptavg --> cannot do this on the histograms... need to do it in event loop,
//because one has to vary the pt and eta of BOTH the leading jets, which are not accessible anymore after Loop()

//const char *a = algo.c_str();

//create vector containing string names of uncertainty sources
//these names have been extracted from Autumn18_V19_MC_UncertaintySources_AK4PF.txt
//but should be the same names also for chs etc.
vector<string> snames;
snames.push_back("AbsoluteStat");
snames.push_back("AbsoluteScale");
snames.push_back("AbsoluteSample");
snames.push_back("AbsoluteFlavMap");
snames.push_back("AbsoluteMPFBias");
snames.push_back("Fragmentation");
snames.push_back("SinglePionECAL");
snames.push_back("SinglePionHCAL");
snames.push_back("FlavorQCD");
snames.push_back("TimePtEta");
snames.push_back("RelativeJEREC1");
snames.push_back("RelativeJEREC2");
snames.push_back("RelativeJERHF");
snames.push_back("RelativePtBB");
snames.push_back("RelativePtEC1");
snames.push_back("RelativePtEC2");
snames.push_back("RelativePtHF");
snames.push_back("RelativeBal");
snames.push_back("RelativeSample");
snames.push_back("RelativeFSR");
snames.push_back("RelativeStatFSR");
snames.push_back("RelativeStatEC");
snames.push_back("RelativeStatHF");
snames.push_back("PileUpDataMC");
snames.push_back("PileUpPtRef");
snames.push_back("PileUpPtBB");
snames.push_back("PileUpPtEC1");
snames.push_back("PileUpPtEC2");
snames.push_back("PileUpPtHF");
snames.push_back("PileUpMuZero");
snames.push_back("PileUpEnvelope");
snames.push_back("SubTotalPileUp");
snames.push_back("SubTotalRelative");
snames.push_back("SubTotalPt");
snames.push_back("SubTotalScale");
snames.push_back("SubTotalAbsolute");
snames.push_back("SubTotalMC");
snames.push_back("Total");
snames.push_back("TotalNoFlavor");
snames.push_back("TotalNoTime");
snames.push_back("TotalNoFlavorNoTime");
snames.push_back("FlavorZJet");
snames.push_back("FlavorPhotonJet");
snames.push_back("FlavorPureGluon");
snames.push_back("FlavorPureQuark");
snames.push_back("FlavorPureCharm");
snames.push_back("FlavorPureBottom");
snames.push_back("TimeRunA");
snames.push_back("TimeRunB");
snames.push_back("TimeRunC");
snames.push_back("TimeRunD");
snames.push_back("CorrelationGroupMPFInSitu");
snames.push_back("CorrelationGroupIntercalibration");
snames.push_back("CorrelationGroupbJES");
snames.push_back("CorrelationGroupFlavor");
snames.push_back("CorrelationGroupUncorrelated");


//Create new folder (and histo? for this source) (in each of the ybys- and trigger-bins)
//loop through ybys bins
//loop through trigger folders
TDirectory *uncdir = curtrigdir->mkdir("uncertainty_histos");	//make uncertainty histo directory

//create uncertainty sources vector (see jecsys package)
vector<JetCorrectionUncertainty*> uncs(snames.size());
for(unsigned int k=0; k!=snames.size(); ++k){
	//file where the uncs are stored
	string sfilename = Form("/Documents/JetAnalysis/unc_JEC_files/Autumn18_V19_MC/Autumn18_V19_MC_UncertaintySources_AK4PF.txt");
	//or with choice of algorithm
	//string sfilename = Form("/Documents/JetAnalysis/unc_JEC_files/Autumn18_V19_MC/Autumn18_V19_MC_UncertaintySources_AK4PF%s.txt", algo);

	const char* sfname = sfilename.c_str();		//.txt file name now as a char
	const char* src	= snames.at(k).c_str();		//uncertainty source name as char, snames.[k].c_str()

	JetCorrectorParameters *p = new JetCorrectorParameters(sfname, src);
	JetCorrectionUncertainty *unc = new JetCorrectionUncertainty(*p);
	uncs[k] = unc;	//set value of the uncertainty vector for current unc source

	//Create new histo for this source (in each of the ybys- and trigger-bins)
}

//do the variation
//Loop through eta?! (cannot do this in my case..)

//BEFORE going into "filling loop" a second time, change the pt of jet1 and jet2
//--> this means adding an extra loop after the nJet>1 selection

//looping again through unc sources. (== variation by variation... create new histos)
for(unsigned int k=0; k!=snames.size(); ++k){
	JetCorrectionUncertainty *unc = uncs.at(k); //uncs[k]

	//check usage of the following functions:
	//calculate uncertainty for the current pt values of the jets
	unc->setJetPt(jet1.pt());
	unc->setJetEta(jet1.eta());
	double unc_j1 = unc->getUncertainty(true);

	unc->setJetPt(jet2.pt());
	unc->setJetEta(jet2.eta());
	double unc_j2 = unc->getUncertainty(true);

	//Calculate new ptavg with new pt values (need old pt values for that!!!)
	newpt1 = jet1.pt()*(1+unc_j1);	//?? like this ??
	newpt2 = ...

	//how to get new yb ys values?? 

	//new ptavg
	newptavg = 0.5*(newpt1+newpt2)


	//fill corresponding histogram
	//is this done within the ORIGINAL trg bin? yeah, right?
	//so added as innermost loop when selecting?
	//For this certain triggerbin, for this certain ybysbin, do the VARIATION?
	histosvec.....unc.at(k) ->Fill(newptavg)


}


*/

//test:
		/*
									//for UNC ANALYSIS -- this was the original attempt --> now: only do this in MC
									//looping again through unc sources. (== variation by variation... create new histos)
									for(unsigned int k=0; k!=uncnames.size(); ++k){
										///JetCorrectionUncertainty *unc = uncs_vec.at(k); //uncs[k]
										//JetCorrectionUncertainty *unc = uncs_vec[k]; assert(unc);	//like in drawSourceCorrelations.C of jecsys


										//calculate relative uncertainty for the current pt values of the jets
										unc->setJetPt(jet1_vec.pt());
										unc->setJetEta(jet1_vec.eta());
										double unc_j1 = unc->getUncertainty(true);

										unc->setJetPt(jet2_vec.pt());
										unc->setJetEta(jet2_vec.eta());
										double unc_j2 = unc->getUncertainty(true);

										//Calculate new ptavg with new pt values (need old pt values for that!!!)
										double newpt1 = jet1_vec.pt()*(1+unc_j1);
										double newpt2 = jet2_vec.pt()*(1+unc_j2);

										//how to get new yb ys values?? 
										//new ptavg
										double newptavg = 0.5*(newpt1+newpt2);

										//fill corresponding histogram
										//is this done within the ORIGINAL trg bin? yeah, right?
										//so added as innermost loop when selecting?
										//For this certain triggerbin, for this certain ybysbin, do the VARIATION?
										cur_ybysbin.unc_histos[uncnames.at(k)]->Fill(newptavg);
									}
									*/





//--------------------------------------------------------------------------------------------------//
//											 PLOTTING												//
//--------------------------------------------------------------------------------------------------//

// taken care of in separate plotting script
