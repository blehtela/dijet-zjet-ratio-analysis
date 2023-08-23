#define ZjetAnalysis_cxx
#include "ZjetAnalysis.h"
//#include "ZjetClasses.h" //contains Muon-, YbYs-, Histos- classes (TO DO)
#include <string>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "tdrstyle_mod14.C"
#include <TMath.h>
#include "Math/LorentzVector.h"

//for get_lumi_txt()
#include <iostream>
#include <fstream>

using namespace std;
using namespace ROOT::Math;

// Analysis of DiMuon data
// --> Z+jet analysis
// use 0.5*(p_{T,Z}+p_{T,jet1}), yboost, ystar

//verbose mode --> TO BE IMPLEMENTED PROPERLY
bool v = false;
//bool v = true;
bool muratio_mode = false;	//if pTmu1/pTmu2 over pTZ is evaluated or not.

// not necessary if MC-data-switch already set in .h file
///bool montecarlo = false;	//to switch between MC and data analysis
///bool montecarlo = true;

//flag for either ybys binning or alternative binning --> not used yet. for now just ybys
string mode = "ybys";
//const char mode = "eta";

//later used for keeping track of the jetindices of the muons, in case there are several Z-candidates in the event
//typedef map<PtEtaPhiMVector, pair<Int_t, Int_t>> JetIndexMap;
//key is index of Z-Boson candidate in Zcands_vec, value is jetindex of first and of second muon
typedef map<int, pair<Int_t, Int_t>> JetIndexMap;


//setTDRStyle();	//should be enough in Loop() actually...

//--------------------------------------------------------------------------------------------------//
//											 CLASSES 												//
//--------------------------------------------------------------------------------------------------//

//class for booking and filling the various histograms (considering trigger decisions)
class Histos;

//class for storing muon variables
//stored properties of Muon: charge, jetIdx, pT, eta, phi, mass
class Muon{
	private:
		Int_t mu_charge;	//charge of muon --> either +1 or -1
		Int_t mu_jetind;	//jet index of current muon (important for later jet selection)
		Double_t mu_pt;		//transversal momentum of muon
		Double_t mu_eta;	//eta of muon
		Double_t mu_phi;	//phi of muon
		Double_t mu_mass;	//mass of muon
		Int_t mu_muind;		//index of the muon in the original MuonArray
	public:
		//constructor
		Muon(Int_t charge, Int_t jet_index, Double_t pt, Double_t eta, Double_t phi, Double_t mass, Int_t muon_index) 
			: mu_charge(charge),
			mu_jetind(jet_index),
			mu_pt(pt),
			mu_eta(eta),
			mu_phi(phi),
			mu_mass(mass),
			mu_muind(muon_index){};
		//functions to access private members of class
		Int_t charge(){return mu_charge;};
		Int_t jet_index(){return mu_jetind;};
		Double_t pt(){return mu_pt;};
		Double_t eta(){return mu_eta;};
		Double_t phi(){return mu_phi;};
		Double_t mass(){return mu_mass;};
		Int_t mu_index(){return mu_muind;};
};

//class for handling yboost ystar binning
class YbYs {
	friend class Histos;
	private:
		string bname;			//bin name, e.g. "yb0ys0" for 0.0 <= yboost, ystar < 0.5, assigned in constructor
		float yb_low;			//lower yboost bin bound
		float yb_up;			//upper yboost bin bound
		float ys_low;			//lower ystar bin bound
		float ys_up;			//upper ystar bin bound

		//if(montecarlo==true){TH1D* mc_histo;};		//histogram for MC analysis
		TH1D* mc_histo;				//histogram for MC analysis
		TH1D* mc_histo_smearedJet;	//histogram for MC HERWIG 7 unc analysis

		//now added analoguously to dijet analysis:
		TH1D* mc_genhisto;		//histogram for MC analysis --> generated spectrum
		TH2D* mc_genrec_noGenCuts;	//histogram as originally implemented
		TH2D* mc_gen_reco_ptav;		//histogram for MC analysis (for pt unfolding)
		TProfile* mc_profile;		//TProfile plot for checking reco and gen distributions
		TProfile* mc_profile_j1;	//TProfile for jet1 pt
		TProfile* mc_profile_j2;	//TProfile for jet2 pt	-- not sensible in Z+jet ...
		TProfile* mc_profile_mergedj1j2;	//TProfile filled with jet1 pt and jet2 pt

		//stuff for herwig jme unc:	(NOT used yet)
		TH1D* jt_ptraw_hist;
		TH1D* jt_ptnom_hist;
		TH1D* jt_corrJEC_hist;
		TH1D* jt_corrJER_hist;
		TH1D* jt_ptjerUp_hist;
		TH1D* jt_ptjerDw_hist;
	

		//TH1D zmass_histo;					//histogram with Zmass distribution
		TH1D* zmass_histo_ptr;				//pointer to histogram containing Zmass distribution for THIS YbYsbin!
		TH1D* zpt_histo_ptr;				//pointer to histogram containing Z-pT distribution ---- " " ----
		TH2D* muratio_ptr;					//pointer to the histogram containing ZpT and mu1/mu2 pt --> lego plot
		//TH1D* zjet_sumhisto_ptr;			//pointer to summary histogram of ZjetAnalysis --> Number of events devided by luminosity!

		TH1D* mu1mu2_ptRatio_ptr;			//pointer to histo containing (mu1 pT)/(mu2 pT)
		vector<Double_t> zpt_vec{};			//vector containing pT values of Zbosons in this ybys bin --> why does it not work with pointers?
		//deque<Double_t> zpt_vec;
		vector<Double_t> *mu1pt_vec;			//vector containing pT values of muon1 used for Z-reconstr. in this ybys bin
		vector<Double_t> *mu2pt_vec;			//vector containing pT values of muon2 used for Z-reconstr. in this ybys bin
		vector<Double_t> mupt_ratio_vec;		//vector containing pT-ratio of muon1 and muon2 for the muons used for Z-reconstr. in this ybysbin
		//deque<Double_t> mupt_ratio_vec;
		TGraph *mu1mu2_ratio_graph;

		map<const char*, TH1D*> trglumi_hist_map;	//used in DijetAnalysis::triggerpaths() for accessing the ratio-histos --> not initialised per default --> must do this if needed!
		const char* ybys_dirname;	//is this ever used? if not --> kick it out

	public:
		YbYs (float, float, float, float);	//constructor with lower and upper bin bound for yboost ystar
		YbYs (float, float, float, float, vector<string>);	//overloaded constructor for uncertainty histos

		//YbYs (float, float);				//constructor with only LOWER bin bound for yb and ys
		//vector <Double_t> get_zpt_vec(){return zpt_vec;};	//return by reference? 
		//void pushpoint(Double_t zpt_value){this->zpt_vec.push_back(zpt_value);}	//for doing push_back on private vector zpt_vec
		void push_Zpt(Double_t zpt_value){zpt_vec.push_back(zpt_value);} //cout << "pushed back zpt." << endl;}	//for doing push_back on private vector zpt_vec
		void push_Muratio(Double_t muratio_value){mupt_ratio_vec.push_back(muratio_value);} 
		vector <Double_t>* get_mu1pt_vec(){return mu1pt_vec;};
		vector <Double_t>* get_mu2pt_vec(){return mu2pt_vec;};
		//vector <Double_t> get_mupt_ratio_vec(){return mupt_ratio_vec;};

		TH1D* get_zmass_histo(){return zmass_histo_ptr;};					//get histogram with Zmass distribution
		TH1D* get_zpt_histo(){return zpt_histo_ptr;};						//get histogram with Zpt distribution
		TH2D* get_muratio_histo(){return muratio_ptr;};						//get 2D-histogram with Zpt and mu1/mu2 pt
		//TH1D* get_sumhist(){return zjet_sumhisto_ptr;};						//return pointer to the summary histogram
		TH1D* get_mu1mu2_ptRatio(){return mu1mu2_ptRatio_ptr;};				//get pT-Ratio of mu1/mu2 (histo)
		pair <float, float> yb_get_bb() {return make_pair(yb_low,  yb_up);};		//get yboost bin bounds
		pair <float, float> ys_get_bb() {return make_pair(ys_low,  ys_up);};		//get yboost bin bounds
		string get_bname() {return bname;};	//to read bin name

		const char* get_dirname(){return ybys_dirname;};	//getting name of directory
		map<string, TH1D*> unc_histos_up;			//map of unc histos (variation upwards)--> key=string uncname, mapped value=TH1D* histogram
		map<string, TH1D*> unc_histos_down;			//map of unc histos, downwards

		TH1D* get_mc_histo(){return mc_histo;}							//get the histogram with the MC
		TH1D* get_mc_histo_smearedJet(){return mc_histo_smearedJet;}	//only used in HERWIG 7 with smeared jet (Jet_pt_nom branch)
		TH1D* get_mc_genhisto(){return mc_genhisto;}					//get the histogram with the MC genjet spectrum
		TH2D* get_mc_genreco_noGenCuts(){return mc_genrec_noGenCuts;}	//get the histogram with gen-reco-ptav distribution, where cuts (pt, eta) have only be applied to RECO
		TH2D* get_mc_genreco(){return mc_gen_reco_ptav;}	//get the histogram with gen-reco-ptav distribution
		TProfile* get_mc_profile(){return mc_profile;}		//get the profile plot for ptavg
		TProfile* get_mc_profile_j1(){return mc_profile_j1;}	//get the profile plot for jet1
		TProfile* get_mc_profile_j2(){return mc_profile_j2;}	//get the profile plot for jet2
		TProfile* get_mc_profile_mergedj1j2(){return mc_profile_mergedj1j2;}	//get the profile plot for jet1-jet2-combi



		//TH1D* zmass_histo_ptr;

		map<const char*, TH1D*> getmap(){return trglumi_hist_map;};	//return map where one can fill in ratio plots (used with DijetAnalysis::triggerpaths())
		void ybys_directories(TDirectory *standard_folder, vector<int> triggers, Histos curtrig); //creates folder-hierarchy
		//void ybys_directories(TDirectory *standard_folder);	//creates folder-hierarchy in case of montecarlo!
		TDirectory* ybys_directories(TDirectory *standard_folder);	//creates folder-hierarchy in case of montecarlo! (version that returns TDirectory*)
		TDirectory* ybys_mkdir(TDirectory *standard_folder);		//only creates the ybys directory and returns pointer to it

		void ybys_uncdirs(TDirectory *unc_folder, vector<string> uncnames);	//used to create uncertainty histograms in Uncertainty folder
		void ybys_jmedirs(TDirectory *jme_folder);	//puts 6 jme histos in the additional jme uncs folder


		void trglumi_directories(TDirectory *standard_folder, vector<int> triggers, map<const char*, TH1D*> trglumi_hist_map); //used in function triggerpaths()
		//void muratio_graph()const;	//creates TGraph mu1mu2_ratio_graph, using zpt_vec and mupt_ratio_vec
		TGraph* muratio_graph();
};

//class for booking the histograms and creating trigger-directories
class Histos{
	//friend class Etabins;
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
		vector<const char*> hjt_names; //just for testing, usually private
		map<const int, TH1D*> jt_histos;		//map of trigger histograms --> key=int pt, mapped value=TH1D *histogram
		//void trigger_histos(Etabins etabin, TH1D *jt_hist);	//function for creating and filling the histos associated with triggers --> currently not used
};



//--------------------------------------------------------------------------------------------------//
//											 CONSTRUCTORS											//
//--------------------------------------------------------------------------------------------------//


/// UNDER CONSTRUCTION !!!! --> compare it with DijetAnalysis.C and tidy up that script at the same time as building this! Careful! ///

//Constructor for Yboost Ystar bin
YbYs::YbYs (float yboost_lo, float ystar_lo, float yboost_up, float ystar_up) {
	yb_low = yboost_lo;
	ys_low = ystar_lo;
	yb_up = yboost_up;	
	ys_up = ystar_up;

	bname = Form("yb%01.0fys%01.0f", yb_low, ys_low);
	if(v==true){cout << Form("Created yboost-ystar-bin with: %f <= yb < %f and %f <= ys < %f", yb_low, yb_up, ys_low, ys_up) << endl;}

	//in case of monte carlo: prepare histo for XS / Events
	/*
	if(montecarlo==true){
		Double_t cust_bb[] = {10, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507, 548, 592, 638, 686, 737, 790, 846, 905, 967, 1032, 1101, 1172, 1248, 1327, 1410, 1497, 1588, 1684, 1784, 1890, 2000, 2116, 2238, 2366, 2500, 2640, 2787, 2941, 3103, 3273, 3450, 3832, 6076, 6389};
		int nbb = sizeof(cust_bb)/sizeof(cust_bb[0])-1;
		mc_histo = new TH1D(Form("mc_histo_%s", bname.c_str()), "MC Z+jet; p_{T,avg} /GeV; Events", nbb, cust_bb);
		if(strcmp(generator, "herwig")==0){
			mc_histo_smearedJet = new TH1D(Form("mc_histo_smearedJet_%s", bname.c_str()), "MC Z+jet (smeared; used Jet_pt_nom); p_{T,avg} /GeV; Events", nbb, cust_bb);
		}
	}
	*/

	Double_t cust_bb[] = {10, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507, 548, 592, 638, 686, 737, 790, 846, 905, 967, 1032, 1101, 1172, 1248, 1327, 1410, 1497, 1588, 1684, 1784, 1890, 2000, 2116, 2238, 2366, 2500, 2640, 2787, 2941, 3103, 3273, 3450, 3832, 6076, 6389};
	int nbb = sizeof(cust_bb)/sizeof(cust_bb[0])-1;
	

	//copied from other constructor...	(to avoid problems with class member functions)
	//in case of monte carlo: prepare histo for XS / Events and one for gen-reco-ptavg distribution (for Unfolding purposes)
	if(montecarlo==true){
		mc_histo = new TH1D(Form("mc_histo_%s", bname.c_str()), "MC Z+jet; p_{T,avg} /GeV; Events", nbb, cust_bb);
		mc_genhisto = new TH1D(Form("mc_genhisto_%s", bname.c_str()), "MC Z+jet GenJet spectrum; p^{gen}_{T,avg} /GeV; Events", nbb, cust_bb);
		mc_genrec_noGenCuts = new TH2D(Form("mc_genrec_noGenCuts_%s", bname.c_str()), "MC Z+jet (no cuts on gen jets); p^{gen}_{T,avg} /GeV; p^{reco}_{T,avg} /GeV; Events", nbb, cust_bb, nbb, cust_bb);
		mc_gen_reco_ptav = new TH2D(Form("mc_gen_reco_ptavg_%s", bname.c_str()), "MC Z+jet p^{reco}_{T,avg} vs. p^{gen}_{T,avg}; p^{gen}_{T,avg} /GeV; p^{reco}_{T,avg} /GeV; Events", nbb, cust_bb, nbb, cust_bb);
		mc_profile =  new TProfile(Form("mc_profile_%s", bname.c_str()), "TProfile of p_{T,avg}; p^{gen}_{T,avg}; #frac{p^{reco}_{T,avg}}{p^{gen}_{T,avg}}", nbb, cust_bb);	//only x bins given?
		mc_profile_j1 =  new TProfile(Form("mc_profile_j1_%s", bname.c_str()), "TProfile of p_{T,j1}; p^{gen}_{T,j1}; #frac{p^{reco}_{T,j1}}{p^{gen}_{T,j1}}", nbb, cust_bb);
		mc_profile_j2 =  new TProfile(Form("mc_profile_j2_%s", bname.c_str()), "TProfile of p_{T,j2}; p^{gen}_{T,j2}; #frac{p^{reco}_{T,j2}}{p^{gen}_{T,j2}}", nbb, cust_bb);
		mc_profile_mergedj1j2 =  new TProfile(Form("mc_profile_mergedj1j2_%s", bname.c_str()), "Merged TProfile of p_{T,j1} and p_{T,j2}; p^{gen}_{T,avg}; #frac{p^{reco}_{T,jets}}{p^{gen}_{T,avg}}", nbb, cust_bb);

		/*
		if(strcmp(generator, "herwig")==0){
			mc_histo_smearedJet = new TH1D(Form("mc_histo_smearedJet_%s", bname.c_str()), "MC Z+jet (smeared; used Jet_pt_nom); p_{T,avg} /GeV; Events", nbb, cust_bb);
		}
		*/
		mc_histo_smearedJet = new TH1D(Form("mc_histo_smearedJet_%s", bname.c_str()), "MC Z+jet (smeared, used Jet_pt_nom); p_{T,avg} /GeV; Events", nbb, cust_bb);	//test: always write it, empty if not Herwig

	

		//setting the colors
		mc_profile->SetMarkerColor(kOrange-3);
		mc_profile_j1->SetMarkerColor(kGreen-2);
		mc_profile_j2->SetMarkerColor(kAzure+7);
		mc_profile_mergedj1j2->SetMarkerColor(kRed-6);

	}

	ybys_dirname = bname.c_str();



	zmass_histo_ptr = new TH1D(Form("zmass_histo_%s", bname.c_str()), "Zmass (cut: 76 GeV < m_{Z} < 106 GeV); m_{Z} /GeV; Events", 120, 60, 120);
	zmass_histo_ptr->SetLineWidth(2);
	zmass_histo_ptr->SetLineColor(kAzure-2);

	//pointer to histogram containing Zpt distribution for current ybys bin
	zpt_histo_ptr = new TH1D(Form("zpt_histo_%s", bname.c_str()), "Z-Boson p_{T} distribution; p_{T,Z} /GeV; Events", 6000, 0, 6000);
	zpt_histo_ptr->SetLineWidth(2);
	zpt_histo_ptr->SetLineColor(kGreen+2);

	//pointer to histogram containing Zpt and mu1/mu2 pt distribution for current ybys bin
	Double_t ratio_bb[] = {1, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 8.0, 9.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0};
	int nbb_ratio = sizeof(ratio_bb)/sizeof(ratio_bb[0])-1;

	muratio_ptr = new TH2D(Form("zpt_muratio_%s", bname.c_str()), "p_{T,Z} and  p_{T}-Ratio of #mu_{1}/#mu_{2}; p_{T,Z} /GeV; p_{T,#mu_{1}}/p_{T,#mu_{2}}; Events", 500, 0, 1000, nbb_ratio, ratio_bb);
	muratio_ptr->SetLineWidth(1);
	muratio_ptr->SetFillColor(kGreen+2);

	muratio_ptr->GetXaxis()->SetTitleSize(0.04);
	muratio_ptr->GetYaxis()->SetTitleSize(0.04);
	muratio_ptr->GetZaxis()->SetTitleSize(0.04);

	muratio_ptr->GetXaxis()->SetTitleOffset(1.6);
	muratio_ptr->GetYaxis()->SetTitleOffset(1.4);
	muratio_ptr->GetZaxis()->SetTitleOffset(1.0);

	muratio_ptr->GetXaxis()->SetLabelSize(0.03);
	muratio_ptr->GetYaxis()->SetLabelSize(0.03);
	muratio_ptr->GetZaxis()->SetLabelSize(0.03);

	muratio_ptr->Draw("LEGO2 0"); // --> quite nice layout
	//c1->SetLogz();
	

	//pointer to histogram containing pT ratio of mu1/mu2 -->should be handled differently... Ratio on y-axis, pTZ on x-axis?
	mu1mu2_ptRatio_ptr = new TH1D(Form("mu1mu2_ptRatio_histo_%s", bname.c_str()), "p_{T}-Ratio of mu1/mu2; p_{T,mu1}/p_{T,mu2}; Events", 6000, 0, 6000);
	mu1mu2_ptRatio_ptr->SetLineWidth(2);
	mu1mu2_ptRatio_ptr->SetLineColor(kAzure);

	//reserving space for the vectors
	///zpt_vec.reserve(180000);
	///mupt_ratio_vec.reserve(180000);

	//TCanvas *zmass_canv = tdrCanvas("zmass_canv", zmass_histo_ptr, 2, 11, kSquare);
	TLine *line_locut = new TLine(76.0, 0, 76.0, 6000);
	TLine *line_upcut = new TLine(106.0, 0, 106.0, 6000);
	//line_locut->Draw("SAME");
	//line_upcut->Draw("SAME");
	line_locut->Draw();
	line_upcut->Draw();
};


// Constructor in use:
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
	//sumhisto = new TH1D(Form("sumhisto_%s", bname.c_str()), "; XS in fb; p_{T,avg} /GeV", nbb, cust_bb);	//summary histo
	zmass_histo_ptr = new TH1D(Form("zmass_histo_%s", bname.c_str()), "Zmass (cut: 76 GeV < m_{Z} < 106 GeV); m_{Z} /GeV; Events", 120, 60, 120);
	zmass_histo_ptr->SetLineWidth(2);
	zmass_histo_ptr->SetLineColor(kAzure-2);

	//pointer to histogram containing Zpt distribution for current ybys bin
	zpt_histo_ptr = new TH1D(Form("zpt_histo_%s", bname.c_str()), "Z-Boson p_{T} distribution; p_{T,Z} /GeV; Events", 6000, 0, 6000);
	zpt_histo_ptr->SetLineWidth(2);
	zpt_histo_ptr->SetLineColor(kGreen+2);


	//remove the following at some point... won't switch it on anyway
	//pointer to histogram containing Zpt and mu1/mu2 pt distribution for current ybys bin
	Double_t ratio_bb[] = {1, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 8.0, 9.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0};
	int nbb_ratio = sizeof(ratio_bb)/sizeof(ratio_bb[0])-1;

	muratio_ptr = new TH2D(Form("zpt_muratio_%s", bname.c_str()), "p_{T,Z} and  p_{T}-Ratio of #mu_{1}/#mu_{2}; p_{T,Z} /GeV; p_{T,#mu_{1}}/p_{T,#mu_{2}}; Events", 500, 0, 1000, nbb_ratio, ratio_bb);
	muratio_ptr->SetLineWidth(1);
	muratio_ptr->SetFillColor(kGreen+2);

	muratio_ptr->GetXaxis()->SetTitleSize(0.04);
	muratio_ptr->GetYaxis()->SetTitleSize(0.04);
	muratio_ptr->GetZaxis()->SetTitleSize(0.04);

	muratio_ptr->GetXaxis()->SetTitleOffset(1.6);
	muratio_ptr->GetYaxis()->SetTitleOffset(1.4);
	muratio_ptr->GetZaxis()->SetTitleOffset(1.0);

	muratio_ptr->GetXaxis()->SetLabelSize(0.03);
	muratio_ptr->GetYaxis()->SetLabelSize(0.03);
	muratio_ptr->GetZaxis()->SetLabelSize(0.03);

	muratio_ptr->Draw("LEGO2 0"); // --> quite nice layout
	//c1->SetLogz();
	

	//pointer to histogram containing pT ratio of mu1/mu2 -->should be handled differently... Ratio on y-axis, pTZ on x-axis?
	mu1mu2_ptRatio_ptr = new TH1D(Form("mu1mu2_ptRatio_histo_%s", bname.c_str()), "p_{T}-Ratio of mu1/mu2; p_{T,mu1}/p_{T,mu2}; Events", 6000, 0, 6000);
	mu1mu2_ptRatio_ptr->SetLineWidth(2);
	mu1mu2_ptRatio_ptr->SetLineColor(kAzure);



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
		mc_histo = new TH1D(Form("mc_histo_%s", bname.c_str()), "MC Z+jet; p_{T,avg} /GeV; Events", nbb, cust_bb);
		mc_genhisto = new TH1D(Form("mc_genhisto_%s", bname.c_str()), "MC Z+jet GenJet spectrum; p^{gen}_{T,avg} /GeV; Events", nbb, cust_bb);
		mc_genrec_noGenCuts = new TH2D(Form("mc_genrec_noGenCuts_%s", bname.c_str()), "MC Z+jet (no cuts on gen jets); p^{gen}_{T,avg} /GeV; p^{reco}_{T,avg} /GeV; Events", nbb, cust_bb, nbb, cust_bb);
		mc_gen_reco_ptav = new TH2D(Form("mc_gen_reco_ptavg_%s", bname.c_str()), "MC Z+jet p^{reco}_{T,avg} vs. p^{gen}_{T,avg}; p^{gen}_{T,avg} /GeV; p^{reco}_{T,avg} /GeV; Events", nbb, cust_bb, nbb, cust_bb);
		mc_profile =  new TProfile(Form("mc_profile_%s", bname.c_str()), "TProfile of p_{T,avg}; p^{gen}_{T,avg}; #frac{p^{reco}_{T,avg}}{p^{gen}_{T,avg}}", nbb, cust_bb);	//only x bins given?
		mc_profile_j1 =  new TProfile(Form("mc_profile_j1_%s", bname.c_str()), "TProfile of p_{T,j1}; p^{gen}_{T,j1}; #frac{p^{reco}_{T,j1}}{p^{gen}_{T,j1}}", nbb, cust_bb);
		mc_profile_j2 =  new TProfile(Form("mc_profile_j2_%s", bname.c_str()), "TProfile of p_{T,j2}; p^{gen}_{T,j2}; #frac{p^{reco}_{T,j2}}{p^{gen}_{T,j2}}", nbb, cust_bb);
		mc_profile_mergedj1j2 =  new TProfile(Form("mc_profile_mergedj1j2_%s", bname.c_str()), "Merged TProfile of p_{T,j1} and p_{T,j2}; p^{gen}_{T,avg}; #frac{p^{reco}_{T,jets}}{p^{gen}_{T,avg}}", nbb, cust_bb);

		/*
		if(strcmp(generator, "herwig")==0){
			mc_histo_smearedJet = new TH1D(Form("mc_histo_smearedJet_%s", bname.c_str()), "MC Z+jet (smeared; used Jet_pt_nom); p_{T,avg} /GeV; Events", nbb, cust_bb);
		}
		*/
		mc_histo_smearedJet = new TH1D(Form("mc_histo_smearedJet_%s", bname.c_str()), "MC Z+jet (smeared, used Jet_pt_nom); p_{T,avg} /GeV; Events", nbb, cust_bb);	//test: always write it, empty if not Herwig

	

		//setting the colors
		mc_profile->SetMarkerColor(kOrange-3);
		mc_profile_j1->SetMarkerColor(kGreen-2);
		mc_profile_j2->SetMarkerColor(kAzure+7);
		mc_profile_mergedj1j2->SetMarkerColor(kRed-6);

	}

	ybys_dirname = bname.c_str();
};







//Constructor for Histos (should rather be done once and as member of Etabins/YbYs?)
//TH1::AddDirectory(false); //just testing this
Histos::Histos (vector<int> triggers) { 		//should the constructor delete all existing histograms that will be recreated here? How to avoid potential memory leak? 
												//(histograms different, but histogram names identical for different eta-objects/eta-bins)
	//customised pT bin bounds
	//the following are the usual:
	Double_t cust_bb[] = {10, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507, 548, 592, 638, 686, 737, 790, 846, 905, 967, 1032, 1101, 1172, 1248, 1327, 1410, 1497, 1588, 1684, 1784, 1890, 2000, 2116, 2238, 2366, 2500, 2640, 2787, 2941, 3103, 3273, 3450, 3832, 6076, 6389};

	/*
	//these are a little more bins in the central GeV region --> for triggerpaths()
	const int ncust_bb = 300;		//number of bins to be created
	Double_t minbb = 10;	//lowest bin bound
	Double_t maxbb = 6390;	//highest bin bound
	//Double_t *cust_bb = new Double_t[ncust_bb+1];	//array for bin bounds
	Double_t cust_bb[ncust_bb+1] = {};
	Double_t bmin_log = TMath::Log10(minbb);		//logarithm of lowest bb (--> later used as exponent)
	Double_t bmax_log = TMath::Log10(maxbb);		//logarithm of highest bb (---- "" ----)	
	Double_t dist_log = (bmax_log-bmin_log)/((Double_t)ncust_bb);	//distance of exponents (dist_log=0.3 --> 10^(1.0), 10^(1.3), 10^(1.6),...)
	cout << "bmin_log:" << bmin_log;
	for(int j=0;j!=ncust_bb+1;++j){
		Double_t cur_expo = bmin_log+(j*dist_log);	//next exponent x of 10
		cust_bb[j] = pow(10, cur_expo);				//bin bound --> 10^x
		//cout << "current cust_bb[j] = " << cust_bb[j] << endl;	//test
	}
	//end of "a little more bins"
	*/

	//cout << "sizeof(cust_bb): " << sizeof(cust_bb) << endl;
	//cout << "sizeof(cust_bb[0]): " << sizeof(cust_bb[0]) << endl;
	int nbb = sizeof(cust_bb)/sizeof(cust_bb[0])-1;
	

	//hjt_names.clear();		//empty the perhaps prefilled (via previous calls) vector (but maintain the pointers --> how?)
	for(int j=0; j!=triggers.size(); ++j){ //or start from 1, because first one is "0" --> without trigger (here not taken into account... only matters in SetBranchStatus)
		//if(j==0){cout << "triggers.at(0)= " << triggers.at(j) << endl;}
		trigger_names.push_back(Form("HLT_PFJet%i", triggers.at(j)));
		hjt_names.push_back(Form("hpt_%i", triggers.at(j)));
		TH1D *curhist = new TH1D(hjt_names.back(), "Z+jet; p_{T,avg};Events", nbb, cust_bb);		//x-axis title etc. have to be adjusted to dijet --> do this in plotting script!

		//test with very small bins
		//TH1D *curhist = new TH1D(hjt_names.back(), "jet_pt;p_{T};Events", 63900, 10, 6400);	//book histogram with very fine binning

		jt_histos.insert(make_pair(triggers.at(j), curhist));
		//gDirectory->ls(); //testing
	} 
};


//--------------------------------------------------------------------------------------------------//
//											 MEMBER FUNCTIONS										//
//--------------------------------------------------------------------------------------------------//

//definition of the member functions of the class YbYs
//these functions don't need to get passed any arguments as they use the data members of their objects
//there used to be a eta_bin_histohandler() function, which became obsolete --> see "ReprodOutput_classes_testing.C"

//the two functions (see DijetAnalysis.C) could be generalised to one general "make_directories" function handling ybys or eta
//new function --> for handling the directory structure and saving it to the "intermediate result output file"
void YbYs::ybys_directories(TDirectory *standard_folder, vector<int> triggers, Histos curtrig){	//should also take the histograms (one per trigger) for this ybys-bin
	//create and handle directory structure of output file
	const char *dirname;  
	dirname = Form("%s", bname.c_str());
	standard_folder->mkdir(dirname);
	TDirectory *Dirname = standard_folder->GetDirectory(dirname);
	Dirname->cd();

	//add trigger-subdirectories
	//curtrig = trigger histograms which belong to current ybys-bin
	for (int j=0; j!=triggers.size(); ++j){
		const char *tdir; //name of trigger directory
		tdir = Form("jt%i", triggers.at(j));
		Dirname->mkdir(tdir);
		TDirectory *tdirname = Dirname->GetDirectory(tdir);
		tdirname->cd();
		curtrig.jt_histos[triggers.at(j)]->SetDirectory(tdirname); //move corresponding histogram to correct trigger-directory in current eta-directory
	}

	//make directory for additional plots (Zmass, etc...)
	Double_t cust_bb[] = {10, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507, 548, 592, 638, 686, 737, 790, 846, 905, 967, 1032, 1101, 1172, 1248, 1327, 1410, 1497, 1588, 1684, 1784, 1890, 2000, 2116, 2238, 2366, 2500, 2640, 2787, 2941, 3103, 3273, 3450, 3832, 6076, 6389};
	int nbb = sizeof(cust_bb)/sizeof(cust_bb[0])-1;

	const char *extradir_name = Form("additional_plots_%s", bname.c_str());
	Dirname->mkdir(extradir_name);
	TDirectory *extraDir = Dirname->GetDirectory(extradir_name);
	extraDir->cd();
	//zmass_histo_ptr = new TH1D("zmass_histo", "Zmass; m_{Z} /GeV; Events", nbb, cust_bb);

	zmass_histo_ptr->SetDirectory(extraDir);
	zpt_histo_ptr->SetDirectory(extraDir);
	muratio_ptr->SetDirectory(extraDir);
	mu1mu2_ptRatio_ptr->SetDirectory(extraDir);

	//mu1mu2_ratio_graph->SetDirectory(extraDir); 
	//extraDir->Append(mu1mu2_ratio_graph); //TGraph has no "SetDirectory()" function (to remove it: dir->Remove(TGraph))
	//zmass_histo_ptr->Write();
	//cout << "zmass histo.." << endl;
	//change back to standard_folder:
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
	const char *extradir_name = Form("additional_plots_%s", bname.c_str());
	ybysDir->mkdir(extradir_name);
	TDirectory *extraDir = ybysDir->GetDirectory(extradir_name);
	extraDir->cd();

	zmass_histo_ptr->SetDirectory(extraDir);
	zpt_histo_ptr->SetDirectory(extraDir);
	muratio_ptr->SetDirectory(extraDir);			// leave this out in this version
	mu1mu2_ptRatio_ptr->SetDirectory(extraDir);	// leave this out in this version


	//change back to standard_folder:
	standard_folder->cd();

	return ybysDir;
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






void YbYs::trglumi_directories(TDirectory *standard_folder, vector<int> triggers, map<const char*, TH1D*> trglumi_hist_map){
	//create and handle directory structure of output file
	const char *dirname;  
	dirname = Form("%s", bname.c_str());
	standard_folder->mkdir(dirname);
	TDirectory *Dirname = standard_folder->GetDirectory(dirname);
	Dirname->cd();

	for(int j=1; j!=triggers.size(); ++j){ //or start from 1, because first one is "0" --> without trigger
		int curtrig = triggers.at(j); 		//current trigger
		int prevtrig = triggers.at(j-1); 	//previous trigger
		TH1D *curhist =  new TH1D(Form("jet%i_jet%i", curtrig, prevtrig), Form("trigger comparison;p_{T,avg};#frac{jet%i / lumi%i}{jet%i / lumi%i}", curtrig, curtrig, prevtrig, prevtrig), 62900, 10, 6300);
		curhist->SetDirectory(Dirname);
		trglumi_hist_map.insert(make_pair(Form("jet%i_jet%i", curtrig, prevtrig), curhist));
	} 
}

//void YbYs::muratio_graph()const{
TGraph* YbYs::muratio_graph(){
	int npoints = zpt_vec.size();
	Double_t xvals[npoints], yvals[npoints];	//make arrays for x- and y-values of TGraph
	cout << endl;
	cout << Form("In %s found %d Z+jet events.", bname.c_str(), npoints) << endl;

	for(int k=0; k!=npoints; ++k){
		xvals[k] = zpt_vec.at(k);
		yvals[k] = mupt_ratio_vec.at(k);
	}
	cout << "npoints, sizeof(xvals), sizeof(yvals): " << npoints << " " << sizeof(xvals)/sizeof(xvals[0]) << " " << sizeof(yvals)/sizeof(yvals[0]) << endl;
	//TGraph *curgraph = new TGraph(npoints, xvals, yvals);
	mu1mu2_ratio_graph = new TGraph(npoints, &xvals[0], &yvals[0]);
	//mu1mu2_ratio_graph = new TGraph(npoints, xvals, yvals);

	//curgraph->SetDirectory()	//happens before already
	mu1mu2_ratio_graph->SetTitle(Form("Ratio p_{T,mu1}/p_{T,mu2} in %s", bname.c_str()));
	mu1mu2_ratio_graph->SetName(Form("muratio_graph_%s", bname.c_str()));
	mu1mu2_ratio_graph->GetXaxis()->SetTitle("p_{T,Z} /GeV");
	//mu1mu2_ratio_graph->GetYaxis()->SetTitle("#frac{p_{T,#mu_{1}}}{p_{T,#mu_{2}}}");
	mu1mu2_ratio_graph->GetYaxis()->SetTitle("p_{T,#mu_{1}}/p_{T,#mu_{2}}");
	mu1mu2_ratio_graph->GetYaxis()->SetTitleSize(0.04);
	mu1mu2_ratio_graph->GetYaxis()->SetLabelSize(0.04);
	mu1mu2_ratio_graph->GetXaxis()->SetTitleSize(0.04);
	mu1mu2_ratio_graph->GetXaxis()->SetLabelSize(0.04);

	mu1mu2_ratio_graph->SetLineColor(kOrange);
	mu1mu2_ratio_graph->SetMarkerStyle(kOpenTriangleUp);
	mu1mu2_ratio_graph->SetMarkerColor(kBlue-2);
	mu1mu2_ratio_graph->Draw("APC*");	//here already?

	return mu1mu2_ratio_graph;
}


//--------------------------------------------------------------------------------------------------//
//											 OTHER FUNCTIONS										//
//--------------------------------------------------------------------------------------------------//
//function to reconstruct Zboson from two muons mu1 and mu2
PtEtaPhiMVector Zreco(Muon* mu1, Muon* mu2){
	//first muon
	Double_t mu1_pt = mu1->pt();
	Double_t mu1_eta = mu1->eta();
	Double_t mu1_phi = mu1->phi();
	Double_t mu1_mass = mu1->mass();
	Int_t mu1_jetindex = mu1->jet_index();
	//second muon
	Double_t mu2_pt = mu2->pt();
	Double_t mu2_eta = mu2->eta();
	Double_t mu2_phi = mu2->phi();
	Double_t mu2_mass = mu2->mass();
	Int_t mu2_jetindex = mu2->jet_index();

	//4-vectors for both muons
	PtEtaPhiMVector muon1(mu1_pt, mu1_eta, mu1_phi, mu1_mass);
	PtEtaPhiMVector muon2(mu2_pt, mu2_eta, mu2_phi, mu2_mass);

	//reconstruct Zboson
	PtEtaPhiMVector Zboson = muon1 + muon2;

	return Zboson;	//return Zboson
}

//function to reconstruct the leading jet
//need jetindex of the two muons to make sure that those are not confused for jet
//	DOES NOT WORK YET --> TROUBLE WITH ARRAYS...
/*
PtEtaPhiMVector jetreco(int mu1_jetindex, int mu2_jetindex){
	//find the leading jet which Index != mu1_jetindex and != mu2_jetindex
	Double_t jet_pt;
	Double_t jet_eta;
	Double_t jet_phi;
	Double_t jet_mass;
	for(int k=0; k!=nJet; ++k){	//loop for obtaining properties of leading jet
		if(mu1_jetindex!=k && mu2_jetindex!=k){
			jet_pt = Jet_pt[k];	
			jet_eta = Jet_eta[k];
			jet_phi = Jet_phi[k];
			jet_mass = Jet_mass[k];
			break;	//escape inner for-loop as we found the leading jet
		}
		else{continue;}
	}//end loop for jet-properties

	//4-vector for leading jet
	PtEtaPhiMVector jet(jet_pt, jet_eta, jet_phi, jet_mass);
	return jet;
}
*/

//histo filling function:
///void fillhistos(PtEtaPhiMVector Zboson, PtEtaPhiMVector jet, vector<YbYs> ybys_pointer, vector<int> triggers, vector<Histos*> histosvec, 
void fillhistos(PtEtaPhiMVector Zboson, PtEtaPhiMVector jet, YbYs* ybys_pointer, vector<int> triggers, vector<Histos*> histosvec, 
				vector<bool> cur_trigdec, int jentry, int &count_zjet, int &count_NoBin_ybys, Double_t mupt_ratio){
	//calculate current ptavg, yb, ys
	//need leading jet (which is not a muon!) for this
	Double_t Z_y = Zboson.Rapidity();
	Double_t jet_y = jet.Rapidity();
	Double_t cur_ptavg = 0.5*(Zboson.Pt()+jet.Pt());
	Double_t cur_yboost = 0.5*TMath::Abs(Z_y+jet_y);
	Double_t cur_ystar = 0.5*TMath::Abs(Z_y-jet_y);

	Double_t Zmass = Zboson.M();	//mass of Z-Boson
	Double_t Zpt = Zboson.Pt();		//pT of Z-Boson

	if(v==true){
		//print info
		//cout << endl;
		cout << "#####################################################################" << endl;
		cout << "Variables in Event " << jentry << endl;
		cout << "---------------------------" << endl;
		cout << Form("Zmass = %f", Zboson.M()) << endl;
		cout << Form("ptavg= %f, yb= %f, ys= %f", cur_ptavg, cur_yboost, cur_ystar) << endl << endl;
	}

	//check ybys bin
	for (int ybysnr=0; ybysnr!=6; ++ybysnr){
		YbYs cur_ybysbin = ybys_pointer[ybysnr];
		//YbYs cur_ybysbin = ybys_pointer.at(ybysnr);

		float yb_lo = cur_ybysbin.yb_get_bb().first;
		float yb_up = cur_ybysbin.yb_get_bb().second;
		float ys_lo = cur_ybysbin.ys_get_bb().first;
		float ys_up = cur_ybysbin.ys_get_bb().second;
		if(v==true){
			cout << "yb_lo: " << yb_lo << " ## yb_up: " << yb_up << " ## ys_lo: " << ys_lo << " ## ys_up: " << ys_up << endl;
		}
		if((yb_lo<=cur_yboost && cur_yboost<yb_up) and (ys_lo<=cur_ystar && cur_ystar<ys_up)){	
		//-->better to implement function InInterval(value, low, up) that checks if value in interval... return(lo<=x & x<up)
			if(v==true){
				cout << "Event is in ybys bin: " << cur_ybysbin.get_bname() << endl;
			}
			//plot mass to distribution
			//cout << "zmass histo at: " << cur_ybysbin.get_zmass_histo() << endl;
			cur_ybysbin.get_zmass_histo()->Fill(Zmass); 	//works because working with a pointer here... --> cur_ybysbin itself is not the original object
			cur_ybysbin.get_zpt_histo()->Fill(Zpt);
			cur_ybysbin.get_muratio_histo()->Fill(Zpt, mupt_ratio);	//filling 2D histogram --> lego plot
			//cur_ybysbin.zmass_histo->Fill(Zmass);
			//cur_ybysbin.get_zpt_vec().push_back(Zpt);
			if(muratio_mode){
				cur_ybysbin.push_Zpt(Zpt);	//does not work, would just be applied on copy (cur_ybysbin)
				ybys_pointer[ybysnr].push_Zpt(Zpt);//need to work with the original object pointed to by ybys_pointer[ybysnr]
				ybys_pointer[ybysnr].push_Muratio(mupt_ratio);	//append ratio to mu1pt/mu2pt vector
			}

			//cur_ybysbin.pushpoint(2.3);	//test
			//cur_ybysbin.get_mupt_ratio_vec().push_back(mupt_ratio);	//append ratio to mu1pt/mu2pt vector

			//check trigger
			for(int trignr=0; trignr!=triggers.size(); ++trignr){ //triggers are checked per event, so correspond to leading jet (???)
				if(cur_trigdec.at(trignr)){	//==true
					histosvec.at(ybysnr)->jt_histos[triggers.at(trignr)]->Fill(cur_ptavg);	//fill ptavg trigger histogram
					if(trignr==0){
						count_zjet += 1;// count z+jet event that has been within our phase space of choice
						if(v==true){cout << "NOW ZJET! " << endl << endl;}
					}
					continue;
				}//trigger true?
				else{//trigger false
					continue;
				}//trigger false
			}//check triggers
			break; //break out of ybys-bin-iteration, because event has now been assigned to its corresponding ybys bin (don't need to check further ybys bins, it is unique)
		}//if zjet-event in ybys bin
		else{	//if no ybys bin was suitable
			if(ybysnr==5){
				if(v==true){cout << "last ybys bin ---> no suitable ybys bin found" << endl;}
				count_NoBin_ybys += 1;
				break;
			}
			continue;
		}//no suitable ybys bin
	}//loop through ybys-bins
}

//fillhistos() version for Monte Carlo (MC)
//could also include this in above version of fillhistos(), but need to set default values for trigger-stuff...
//last options are just called in case of herwig 7 -- uncertainty analysis stuff
void fillhistos(PtEtaPhiMVector Zboson, PtEtaPhiMVector jet, YbYs* ybys_pointer, int jentry, int &count_zjet, int &count_NoBin_ybys, Double_t mupt_ratio, Double_t gen_weight){

	//calculate current ptavg, yb, ys
	//need leading jet (which is not a muon!) for this
	Double_t Z_y = Zboson.Rapidity();
	Double_t jet_y = jet.Rapidity();
	Double_t cur_ptavg = 0.5*(Zboson.Pt()+jet.Pt());
	Double_t cur_yboost = 0.5*TMath::Abs(Z_y+jet_y);
	Double_t cur_ystar = 0.5*TMath::Abs(Z_y-jet_y);

	Double_t Zmass = Zboson.M();	//mass of Z-Boson
	Double_t Zpt = Zboson.Pt();		//pT of Z-Boson

	//check ybys bin
	for (int ybysnr=0; ybysnr!=6; ++ybysnr){
		//YbYs cur_ybysbin = ybys_pointer[ybysnr];
		//YbYs cur_ybysbin = ybys_pointer.at(ybysnr);

		YbYs cur_ybysbin = ybys_pointer[ybysnr];	//this is copying the object pointed to by ybys_pointer[ybysnr] (including corresp. pointers)
		///YbYs* cur_ybysbin = &ybys_pointer[ybysnr];	//this is now just renaming the pointer!


		float yb_lo = cur_ybysbin.yb_get_bb().first;
		float yb_up = cur_ybysbin.yb_get_bb().second;
		float ys_lo = cur_ybysbin.ys_get_bb().first;
		float ys_up = cur_ybysbin.ys_get_bb().second;
		//if(jentry%100000==0){
		//	cout << "yb_lo: " << yb_lo << " ## yb_up: " << yb_up << " ## ys_lo: " << ys_lo << " ## ys_up: " << ys_up << endl;
		//}


		if((yb_lo<=cur_yboost && cur_yboost<yb_up) and (ys_lo<=cur_ystar && cur_ystar<ys_up)){	
		//-->better to implement function InInterval(value, low, up) that checks if value in interval... return(lo<=x & x<up)
			//plot mass to distribution
			cur_ybysbin.get_zmass_histo()->Fill(Zmass); 	//works because working with a pointer here... --> cur_ybysbin itself is not the original object
			cur_ybysbin.get_zpt_histo()->Fill(Zpt);
			cur_ybysbin.get_muratio_histo()->Fill(Zpt, mupt_ratio);	//filling 2D histogram --> lego plot
			if(muratio_mode){
				cur_ybysbin.push_Zpt(Zpt);	//does not work, would just be applied on copy (cur_ybysbin)
				ybys_pointer[ybysnr].push_Zpt(Zpt);//need to work with the original object pointed to by ybys_pointer[ybysnr]
				ybys_pointer[ybysnr].push_Muratio(mupt_ratio);	//append ratio to mu1pt/mu2pt vector
			}

			//do not check any triggers (as there are none in MC)
			//fill current ybys histogram with MC
			cur_ybysbin.get_mc_histo()->Fill(cur_ptavg, gen_weight);	//apply generator weight (as this is MC)
			//count_zjet += 1;												//count z+jet event that has been within our phase space of choice
			break; //break out of ybys-bin-iteration, because event has now been assigned to its corresponding ybys bin (don't need to check further ybys bins, it is unique)
		}//if zjet-event in ybys bin
		else{	//if no ybys bin was suitable
			if(ybysnr==5){
				if(v==true){cout << "last ybys bin ---> no suitable ybys bin found" << endl;}
				count_NoBin_ybys += 1;
				break;
			}
			continue;
		}//no suitable ybys bin
	}//loop through ybys-bins
}

//overloaded a second time -- for herwig and uncertainties
void fillhistos(PtEtaPhiMVector Zboson, PtEtaPhiMVector jet, YbYs* ybys_pointer, int jentry, int &count_zjet, int &count_NoBin_ybys, Double_t mupt_ratio, Double_t gen_weight, vector<string> uncnames, Float_t Jet_pt_uncUp[56][47], Float_t Jet_pt_uncDw[56][47], int jet_jetindex){

	//calculate current ptavg, yb, ys
	//need leading jet (which is not a muon!) for this
	Double_t Z_y = Zboson.Rapidity();
	Double_t jet_y = jet.Rapidity();
	Double_t cur_ptavg = 0.5*(Zboson.Pt()+jet.Pt());
	Double_t cur_yboost = 0.5*TMath::Abs(Z_y+jet_y);
	Double_t cur_ystar = 0.5*TMath::Abs(Z_y-jet_y);

	Double_t Zmass = Zboson.M();	//mass of Z-Boson
	Double_t Zpt = Zboson.Pt();		//pT of Z-Boson

	//check ybys bin
	for (int ybysnr=0; ybysnr!=6; ++ybysnr){
		//YbYs cur_ybysbin = ybys_pointer[ybysnr];
		//YbYs cur_ybysbin = ybys_pointer.at(ybysnr);

		YbYs cur_ybysbin = ybys_pointer[ybysnr];	//this is copying the object pointed to by ybys_pointer[ybysnr] (including corresp. pointers)
		///YbYs* cur_ybysbin = &ybys_pointer[ybysnr];	//this is now just renaming the pointer!


		float yb_lo = cur_ybysbin.yb_get_bb().first;
		float yb_up = cur_ybysbin.yb_get_bb().second;
		float ys_lo = cur_ybysbin.ys_get_bb().first;
		float ys_up = cur_ybysbin.ys_get_bb().second;
		//if(jentry%100000==0){
		//	cout << "yb_lo: " << yb_lo << " ## yb_up: " << yb_up << " ## ys_lo: " << ys_lo << " ## ys_up: " << ys_up << endl;
		//}


		if((yb_lo<=cur_yboost && cur_yboost<yb_up) and (ys_lo<=cur_ystar && cur_ystar<ys_up)){	
		//-->better to implement function InInterval(value, low, up) that checks if value in interval... return(lo<=x & x<up)
			//plot mass to distribution
			cur_ybysbin.get_zmass_histo()->Fill(Zmass); 	//works because working with a pointer here... --> cur_ybysbin itself is not the original object
			cur_ybysbin.get_zpt_histo()->Fill(Zpt);
			cur_ybysbin.get_muratio_histo()->Fill(Zpt, mupt_ratio);	//filling 2D histogram --> lego plot
			if(muratio_mode){
				cur_ybysbin.push_Zpt(Zpt);	//does not work, would just be applied on copy (cur_ybysbin)
				ybys_pointer[ybysnr].push_Zpt(Zpt);//need to work with the original object pointed to by ybys_pointer[ybysnr]
				ybys_pointer[ybysnr].push_Muratio(mupt_ratio);	//append ratio to mu1pt/mu2pt vector
			}

			//do not check any triggers (as there are none in MC)
			//fill current ybys histogram with MC
			if(jet.pt()>30){	//this if-cond is only needed in herwig with the smeared jet (otherwise always checked before calling fillhistos())
				//cur_ybysbin.get_mc_histo()->Fill(cur_ptavg, gen_weight);	//apply generator weight (as this is MC)
				cur_ybysbin.get_mc_histo_smearedJet()->Fill(cur_ptavg, gen_weight);	//apply generator weight (as this is MC)
				count_zjet += 1;												//count z+jet event that has been within our phase space of choice
			}

			//now do uncertainty analysis if Herwig7 MC set
			//This event has passed event selection --> do uncertainty analysis (if herwig)
			if(strcmp(generator, "herwig")==0){	//kind of unnecessary, as this function is only called with herwig...
				//need to "loop through the branches" of the Friends TTree -- don't have to loop. just take the values. (from the map)
				//The branch usually has the same name as the Unc-histo, just added an "Up" or "Down"
				for(unsigned int k=0; k!=uncnames.size(); ++k){
					/*
					if(jentry%10000==0){
						cout << "k==" << k << endl;
						cout << "unc = " << uncnames.at(k) << endl;
					}
					*/

					//Get the corresponding branches
					//these already contain the varied Jet_pt value
					const char* unc_up_branchname = Form("Jet_pt_jes%sUp", uncnames.at(k).c_str());
					const char* unc_down_branchname = Form("Jet_pt_jes%sDown", uncnames.at(k).c_str());

					//Get the variations -- only vary leading jet (Z stays unchanged)
					//varied leading jet pt
					double newpt1_up = Jet_pt_uncUp[k][jet_jetindex];
					double newpt1_down = Jet_pt_uncDw[k][jet_jetindex];

					//unchanged Zboson pt:
					double Zpt = Zboson.Pt();

					//new ptavg
					double newptavg_up = 0.5*(newpt1_up+Zpt);
					double newptavg_down = 0.5*(newpt1_down+Zpt);

					/*
					if(jentry%10000==0){		//testing
						cout << "old pt1: " << jet.pt() << endl;	
						cout << "new pt1 up: " << newpt1_up << endl;
						cout << "new pt1 down: " << newpt1_down << endl;
					}
					*/

					//check if the leading jet still fulfills pT > 30 GeV
					if(newpt1_up>30){
						//fill corresponding histogram
						//fill the histograms (could also write function instead)
						cur_ybysbin.unc_histos_up[uncnames.at(k)]->Fill(newptavg_up, gen_weight);
					}
					if(newpt1_down>30){
						cur_ybysbin.unc_histos_down[uncnames.at(k)]->Fill(newptavg_down, gen_weight);
					}


				}//unc sources loop

				//now fill as well Jet_pt_raw, Jet_pt_nom from uncfile (additional information)

			}//if mc HERWIG (end of unc-loop)

			break; //break out of ybys-bin-iteration, because event has now been assigned to its corresponding ybys bin (don't need to check further ybys bins, it is unique)
		}//if zjet-event in ybys bin
		else{	//if no ybys bin was suitable
			if(ybysnr==5){
				if(v==true){cout << "last ybys bin ---> no suitable ybys bin found" << endl;}
				count_NoBin_ybys += 1;
				break;
			}
			continue;
		}//no suitable ybys bin
	}//loop through ybys-bins
}






//copied from my dijet analysis script: (there also a TGraphErrors version)
//should be used more extensively here as well
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





//--------------------------------------------------------------------------------------------------//
//											 ANALYSIS LOOP											//
//--------------------------------------------------------------------------------------------------//

/*
void ZjetAnalysis::Loop(){
{
//   In a ROOT session, you can do:
//      root> .L ZjetAnalysis.C
//      root> ZjetAnalysis t
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
}
*/

void ZjetAnalysis::Loop(){
	TH1::SetDefaultSumw2(kTRUE);      //in order to ALWAYS store the bin errors
	setTDRStyle();
	std::string name_h = "hpt";
	fChain->SetBranchStatus("*",0); //exclude / switch off all the branches
	fChain->SetBranchStatus("nJet",1); //include this branch in analysis	
	fChain->SetBranchStatus("Jet_pt",1); //include this branch as well
	fChain->SetBranchStatus("Jet_eta",1); //has to be included if eta-bins shall be checked...
	fChain->SetBranchStatus("Jet_phi",1);	//for TLorentzVector
	fChain->SetBranchStatus("Jet_mass",1);	//for TLorentzVector

	fChain->SetBranchStatus("nMuon",1);
	fChain->SetBranchStatus("Muon_mediumPromptId",1);
	fChain->SetBranchStatus("Muon_charge",1);
	fChain->SetBranchStatus("Muon_jetIdx",1);
	fChain->SetBranchStatus("Muon_pt",1);
	fChain->SetBranchStatus("Muon_eta",1);
	fChain->SetBranchStatus("Muon_phi",1);
	fChain->SetBranchStatus("Muon_mass",1);
	fChain->SetBranchStatus("HLT_Mu18_Mu9");
	fChain->SetBranchStatus("HLT_Mu23_Mu12",1);
	fChain->SetBranchStatus("Muon_pfIsoId",1);

	//for noise reduction (MET checks)
	fChain->SetBranchStatus("ChsMET_pt",1);		//MET in event
	fChain->SetBranchStatus("ChsMET_sumEt",1);	//Et sum in event
	fChain->SetBranchStatus("MET_pt",1);		//MET in event
	fChain->SetBranchStatus("MET_sumEt",1);	//Et sum in event



	if(montecarlo){
		cout << "This is a MC analysis. Switching on MC-specific branches." << endl;
		fChain->SetBranchStatus("Generator_weight",1);

		//the following is taken from the DijetAnalysis.C, not used yet in this script
		//fChain->SetBranchStatus("nGenJet",1);
		//fChain->SetBranchStatus("GenJet_pt",1);
		//fChain->SetBranchStatus("GenJet_eta",1);
		//fChain->SetBranchStatus("GenJet_phi",1);

		//fChain->SetBranchStatus("Generator_binvar",1);	//for pthat checks -- only in dijet

	}
	else{//if data
		cout << Form("This is an analysis of Z+jet data. Using data from Run period: %s", runperiod) << endl;
	}


	//-------------- building folder structure --------------//
	// output file where folder structure will be placed in
	// also luminosity of HLT_Mu23_Mu12 is assigned here (according to runperiod) //
	double lumi_mutrg;		//lumi of HLT_Mu23_Mu12 trigger [in 1/fb]
	TFile* interim_output;
	if(montecarlo==false){//if data
		//check which runperiod
		//prepare output file (= after analysis BEFORE creating plots)
		if(strcmp(runperiod,"A")==0){
			cout << "Run A" << endl;
			interim_output = new TFile("zjet_interim_output_data_RunA.root", "RECREATE");
			lumi_mutrg = 0.640070681;
		}//endif runA

		else if(strcmp(runperiod,"B")==0){
			cout << "Run B" << endl;
			interim_output = new TFile("zjet_interim_output_data_RunB.root", "RECREATE");
			lumi_mutrg = 0.352869651;
		}//endif runB

		else if(strcmp(runperiod,"C")==0){
			cout << "Run C" << endl;
			interim_output = new TFile("zjet_interim_output_data_RunC.root", "RECREATE");
			lumi_mutrg = 0.345769690;
		}//endif runC

		else if(strcmp(runperiod,"D")==0){
			cout << "Run D" << endl;
			interim_output = new TFile("zjet_interim_output_data_RunD.root", "RECREATE");
			lumi_mutrg = 1.587709963;
		}//endif runD

		else{//if no runperiod specified
		interim_output = new TFile("zjet_interim_output_data.root", "RECREATE"); //output file after analysis BEFORE creating plots
		}
	}
	else{ //if MC (for Z+jet only used pythia so far, now added herwig7)
		if(strcmp(generator, "herwig")==0){
			interim_output = new TFile("zjet_interim_output_mc_herwig7.root", "RECREATE"); //output file if herwig7

			fFriends->SetBranchStatus("nJet", 0);
			//fFriends->SetBranchStatus("Jet_pt_raw", 1);
			//fFriends->SetBranchStatus("Jet_pt_nom", 1);
			fFriends->SetBranchStatus("Jet_pt_uncUp", 1);
			fFriends->SetBranchStatus("Jet_pt_uncDw", 1);
			fFriends->SetBranchStatus("Jet_pt*", 1);
		}
		else if(strcmp(generator, "pythia")==0){
			interim_output = new TFile("zjet_interim_output_mc_pythia8.root", "RECREATE"); //output file after analysis BEFORE creating plots
		}
		else{//if no generator specified -- should not happen ??
			interim_output = new TFile("zjet_interim_output_mc.root", "RECREATE"); //output file after analysis BEFORE trigger analysis or creating plots
		}
	}
	interim_output->cd();		//change to this outputfile
	TDirectory *curdir = gDirectory; //get current directory
	interim_output->mkdir("Standard"); //create directory in output file

	//create subfolder "Standard" in outputfile
	TDirectory *standard_folder = interim_output->GetDirectory("Standard"); 
	gROOT->GetListOfBrowsables()->Add(standard_folder, "Standard");
	standard_folder->cd();

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



	//---------------- First preparation for UncSources --------------//
	//create vector containing string names of uncertainty sources
	//these names have been extracted from Autumn18_V19_MC_UncertaintySources_AK4PF.txt
	//but should be the same names also for chs etc.
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





	//later used in normalise_mc() function --> generator weights
	//Int_t genbins = 100;
	//Double_t genmin = fChain->GetMinimum("Generator_weight");
	//Double_t genmax = fChain->GetMaximum("Generator_weight");
	Double_t genmin = -1.06;	//maximum x in mc_genweight histo
	Double_t genmax = 1.06;		//minimum x in mc_genweight histo
	Int_t	genbins = 106;		//number of bins in mc_genweight histo

	TH1D* mc_genweight = new TH1D("mc_genweight", "Distribution before any selection; Generator_weight; entries", genbins, genmin, genmax);	//for normalising later - histo with generator weights
	mc_genweight->SetDirectory(standard_folder);	//this histo is only created once for the whole data set



	//vector containing the different trigger-pt-values (same for each eta-bin)
	//if this was etabin-specific, it would be better to make it an Etabin-class member, but it is general here.
	//vector<int> triggers = {0, 40, 60, 80, 140, 200, 260, 320, 400, 450, 500};
	vector<int> triggers = {0, 15, 25, 40, 60, 80, 140, 200, 260, 320, 400, 450, 500, 550};

	//switching on selected TBranches corresponding to chosen triggers:
	for (int j=1; j!=triggers.size(); ++j){ //start from 1, because "0" means "without trigger"
		char *trigname = Form("HLT_PFJet%i", triggers.at(j));
		fChain->SetBranchStatus(trigname, 1);
	}	

	// counters
	int count_jetless = 0; //count number of events without a jet (nJet==0)
	int count_higheta = 0; //count number of events where leading jet has eta > 4.7
	int count_less2mu = 0; //count number of events that cannt have a Z, i.e. nMuon < 2
	int count_NoBin_ybys = 0; //count number of events that have no ybys bin, where they could be assigned to (i.e. yb or ys > 3.0)
	int count_NoMediumPromptId = 0; //count how many muons do not pass Muon_mediumPromptId selection
	int count_dropout_MediumPromptId = 0; //count how many events (!) with nMuon>1 drop out due to Muon_mediumPromptId selection
	
	int count_zjet = 0;		//count how many z+jet events actually "survive" selection

	int count_No12Z = 0;	//count events that passed Muon_mediumPromptId selection, but Z cannot be reconstructed with muon 1 and 2, due to charge
	int count_wrongmass = 0; //events that passed all previous cuts, but have wrong Z-mass

	//testing!
	int n_0mu = 0;	//count events without muon
	int n_1mu = 0;	//count events with 1 muon
	int n_2mu = 0;	//count events with 2 muons
	int n_3mu = 0;	//count events with 3 muons
	int n_4mu = 0;	//count events with 4 muons
	int n_5mu = 0;	//count events with 5 or more muons

	int prompt2 = 0;	//count events with 2 muons passing mediumPromptId
	int prompt3 = 0;	//count events with 3 muons passing mediumPromptId
	int prompt4 = 0;	//count events with 4 muons passing mediumPromptId
	int prompt5 = 0;	//count events with 5 or more  muons passing mediumPromptId

	int n_success_3rd = 0;	//count events where Z boson was successfully reconstructed when using 3rd muon (mu1+mu3 or mu2+mu3) in event

	//muon-trigger dropout
	int count_dropout_MuonTrigger = 0;	//count Z+jet events that drop out due to the chosen muon_trigger
	int count_fail_prompt_OR_trg = 0;	//count Z+jet events that have two muons, but fail because of Muon_mediumPromptId or HLT_Mu23_Mu12 (or both!)

	int count_noZcands = 0;		//count number of events where no good Zcandidate has been found
	int count_manyZcands = 0;	//count number of events where more than one good Zcand has been found

	//jet dropout
	//int count_softjet = 0;		//count number of events where leading jet is too soft (< 30 GeV)
	//int count_highYjet = 0;		//count number of events where leading jet has too high rapidity (abs(y)>=3.0)
	int count_badJet = 0;			//counts events that dropped out due to bad leading jet (combines and replaces the two unused counters above)

	//Zcand dropout
	int count_badZrapidity = 0;		//count events that drop out due to bad Zrapidity
	int count_badZpt = 0;			//count events that drop out due to pTZ <= 30 GeV

	//new
	int count_noGen1 = 0;		//no matching gen jet for leading jet
	int count_jetmatch = 0;		//count events with matched reco and gen leading jet
	int count_20pt30 = 0;		//if everything was fine, but the jetpt is between 20 and 30 GeV (late drop-out)



//########### NEED TO STILL ADJUST THE FOLLOWING TO Z+JET ANALYSIS ! --> EVENT SELECTION ETC. ###################//
	//go through all the events and check muon conditions and others

	YbYs *ybys_pointer;	//should maybe handle this differently
	//call constructor YbYs(float yb_lo, float ys_lo, float yb_up, float ys_up)
	//updated max y=2.4 and added unc_folder stuff
	ybys_pointer = new YbYs[6]{{0.0,0.0,1.0,1.0, uncnames}, {0.0, 1.0, 1.0, 2.0, uncnames}, {0.0, 2.0, 1.0, 2.4, uncnames}, {1.0, 0.0, 2.0, 1.0, uncnames}, {1.0, 1.0, 2.0, 2.0, uncnames}, {2.0, 0.0, 2.4, 1.0, uncnames}};	

	//TEST
	/*
	vector<YbYs> ybys_pointer_pre;
	ybys_pointer_pre.push_back(new YbYs(0.0,0.0,1.0,1.0));	// --> to get vector of pointers
	ybys_pointer_pre.push_back(YbYs(0.0, 1.0, 1.0, 2.0));
	ybys_pointer_pre.push_back(YbYs(0.0, 2.0, 1.0, 3.0));
	const vector<YbYs> ybys_pointer = ybys_pointer_pre;
	*/

	vector <Histos*> histosvec;	//has to be declared anyway...
	for (int ybysnr=0; ybysnr!=6; ++ybysnr){
		YbYs cur_ybysbin = ybys_pointer[ybysnr];
		if(montecarlo==false){//only if data==true --> montecarlo==false
			histosvec.push_back(new Histos(triggers));
			Histos curtrig = *histosvec.back();
			//YbYs cur_ybysbin = ybys_pointer[ybysnr];
			///YbYs cur_ybysbin = ybys_pointer.at(ybysnr);
			cur_ybysbin.ybys_directories(standard_folder, triggers, curtrig);	//creates folder structure
			cur_ybysbin.ybys_uncdirs(unc_folder, uncnames);						//setting the unc histos to the correct folder
		}
		else{//if MC
			cur_ybysbin.ybys_directories(standard_folder);
			cur_ybysbin.ybys_uncdirs(unc_folder, uncnames);						//same procedure for mc
		}//if MC
	}

	if (fChain == 0) return;

	Long64_t nentries;
	if(montecarlo){nentries = fChain->GetEntries();}	//otherwise results in wrong number of entries (as TChain is used in this MC)
	else{//if data and just one file
	nentries = fChain->GetEntriesFast();
	}

	//Long64_t nentries = fChain->GetEntriesFast();	//"normal" version, when looking at single file instead of TChain
	if(v==true){cout << "nentries = " << nentries << endl;}		//information

	//analysis loop
	cout << endl;
	Long64_t nbytes = 0, nb = 0;
	for (Long64_t jentry=0; jentry<nentries;++jentry){
	///for (Long64_t jentry=0; jentry<10000; ++jentry){} //for testing with less entries/events
		Long64_t ientry = LoadTree(jentry);
		///if (jentry%10==0){cout << "------------------------------------------" << endl;}//for overview
		if (jentry%100000==0){cout << "]" << endl << "processed " << jentry << " of " << nentries << endl;} //just for testing -> always. can delete this afterwards
		if(jentry%100000==0){cout << "["; cout.flush();}
		if(jentry%10000==0){cout << "#"; cout.flush();}
		if(v==true){
			if (jentry%100000==0){cout << "processed " << jentry << " of " << nentries << endl;}
		}
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);	nbytes += nb;
		// if (Cut(ientry) < 0) continue;


		vector<Bool_t> cur_trigdec; //has to be declared anyway (because later used in function call... but only if data..)
		if(montecarlo==false){
			//### Check TRIGGER DECISIONS ### ---> REMOVE THE SINGLE-JET TRIGGERS! --> not needed here!!! --> also in folder structure of output: REMOVE IT!
			//take care of trigger decisions in current event (triggers are per-event-variables --> how to handle this for lower-pt jets?)
			//for now: do it by hand... (could, of course, be much more elegant. this is a testing-relict.)
			map<const int, Bool_t> tdecisions;

			// use different triggering in DiMuon-Analysis! (HOW???)
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


			//vector<Bool_t> cur_trigdec; //vector of triggerdecisions for this event (index jentry)! (can be avoided using tdecisions directly!)
			for(int trignr=0; trignr!=triggers.size(); ++trignr){
				cur_trigdec.push_back(tdecisions[triggers.at(trignr)]); //unnecessary... instead: read directly from tdecisions
			}
		}//if montecarlo==false --> only if data
		else{
			mc_genweight->Fill(Generator_weight);		//fill this histo with generator weight of current event, before any cuts!
		}




		if(v==true){cout << Form("nJet in Event %lld is: %d || nMuon is: %d", jentry, nJet, nMuon) << endl;}

		//just a test
		if(nMuon==0){n_0mu += 1;}
		if(nMuon==1){n_1mu += 1;}
		if(nMuon==2){n_2mu += 1;}
		if(nMuon==3){n_3mu += 1;}
		if(nMuon==4){n_4mu += 1;}
		if(nMuon>4){n_5mu += 1;}
	


		//also need to check if event contains at least one jet APART FROM MUONS!
		//in nanoAOD muons are also counted in the jet-variables... --> need to check jetId!
		//check if event (jentry) contains at least two muons!
		if(nMuon>1){	//required for Z-boson! --> only look at events containing a Z --> if(nMuon>1 && nJet>2){}

			//add MET check for reducing noise
			bool good_met = MET_pt<(0.3*MET_sumEt);			//exists in both data and mc

			//loop over muons in the event and count the ones that pass Muon_mediumPromptId = 1
			//store mu-index in list, in order to know which of the muons meet the requirements
			Int_t n_muprompt  = 0;
			vector<Int_t> idx_GoodMuons;
			for(int k=0; k!=nMuon; ++k){
				if(Muon_mediumPromptId[k]==true){
					n_muprompt += 1;		//count muon that passes requirement
					idx_GoodMuons.push_back(k);	//store index (within this current event) of good muon
				}
				else{
					count_NoMediumPromptId += 1;	//count muon that does not pass the MediumPromptId
				}
			}

			//check how many events with at least two muons drop out due to MediumPromptId selection:
			if(n_muprompt<2){
				count_dropout_MediumPromptId += 1;
				if(v==true){cout << Form("Event %lld failed Muon_mediumPromptId selection.", jentry) << endl;}
			}
			
			/*
			if(n_muprompt<2){	//less than two muons (i.e. 0 or 1) passed mediumPromptId selection
				if(Muon_mediumPromptId[0]==false || Muon_mediumPromptId[1]==false){
					count_dropout_MediumPromptId += 1;
					if(v==true){cout << Form("Event %lld failed Muon_mediumPromptId selection.", jentry) << endl;}
				}
			}
			else{	//two or more muons passed mediumPromptId selection
				if(Muon_mediumPromptId[idx_GoodMuons.at(0)]==false || Muon_mediumPromptId[idx_GoodMuons.at(1)]==false){//that does not make sense..
					count_dropout_MediumPromptId += 1;
					if(v==true){cout << Form("Event %lld failed Muon_mediumPromptId selection.", jentry) << endl;}
				}
			}
			*/

			//if(n_muprompt>2){cout << Form("n_muprompt in event %lld is: %d", jentry, n_muprompt) << endl;}
			if(n_muprompt==2){prompt2 += 1;}
			if(n_muprompt==3){prompt3 += 1;}
			if(n_muprompt==4){prompt4 += 1;}
			if(n_muprompt>4){prompt5 += 1;}

			//check all the muons, different amount for MC and data (therefore different uppermost index)
			int n_mu;	//number of muons to be checked
			//int n_mu1;	//uppermost index of outer loop (1st muon)
			//int n_mu2;	//uppermost index of inner loop (2nd muon)
			if(montecarlo){ //if MC
				n_mu = nMuon;
				//n_mu1 = nMuon-1;
				//n_mu2 = nMuon;
			}
			else{	//if data
				n_mu = n_muprompt;
				//n_mu1 = n_muprompt-1;
				//n_mu2 = n_muprompt;
			}


			//go on with processing events where at least two muons have passed the mediumPromptId
			//check the MET condition (good_met) in data and MC
			//AT THE SAME TIME --> only go on if triggercondition is true!
			bool muon_trigger = HLT_Mu23_Mu12;
			//if((n_muprompt>1)&&(muon_trigger)){
			//if MC --> do not check muprompt or mutrg --> just go on.
			if((montecarlo&&good_met)||((n_muprompt>1)&&muon_trigger&&good_met)){
				//Muon *muons_ptr;	//declare it and assign length later? (array of several pointers to the different muons)
				Muon *muons_ptr[n_mu]; //create object of class Muon with properties of the current muon
				
				//muons_ptr contains all muons of current event which passed mediumPromptId selection
				for(int m=0; m!=n_mu; ++m){
					//properties of Muon: charge, jetIdx, pT, eta, phi, mass
					int mu_index;
					if(montecarlo==false){ //if data
						mu_index = idx_GoodMuons.at(m);		//index of muon in current event (to keep it more readable..)
					}
					else{mu_index = m;}	//if MC
					muons_ptr[m] = new Muon{Muon_charge[mu_index], Muon_jetIdx[mu_index], Muon_pt[mu_index], Muon_eta[mu_index], Muon_phi[mu_index], Muon_mass[mu_index], mu_index};
				}
				

				/*
				//if(n_muprompt==2){}
				//############################## TESTING !!! JUST A DRAFT!!! #################################################//
				if(n_muprompt==2){
					Muon* muon1 = muons_ptr[0];
					Muon* muon2 = muons_ptr[1];
					//check if leading muons have opposite sign of charge
					//int Zcharge = muons_ptr[0]->charge()+muons_ptr[1]->charge();	//might not be the optimal procedure... ?
					//and if their invariant mass is about Zmass (see below)
					int Zcharge = muon1->charge()+muon2->charge();
					
					if(Zcharge==0){
						Int_t mu1_jetindex = muons_ptr[0]->jet_index();
						Int_t mu2_jetindex = muons_ptr[1]->jet_index();


					}//end if Zcharge==0

				}
				//else{} //possibility of taking third muon into account, in case first two do not yield "good Z"
				//###########################################################################################################//
				*/

				//test of mu1/mu2 in terms of pT	--> CANNOT DO THIS BEFORE EVENT HAS BEEN SELECTED!
				for (int ybysnr=0; ybysnr!=6; ++ybysnr){
					YbYs cur_ybysbin = ybys_pointer[ybysnr];
					///YbYs cur_ybysbin = ybys_pointer.at(ybysnr);
					cur_ybysbin.get_mu1mu2_ptRatio()->Fill(Muon_pt[0]/Muon_pt[1]);
				}

				//#######################################################
				// ### --- From here onwards: change to new kind of loop where all the muon-combinations of this event are checked for possible Z-candidates --- ###
				//#######################################################
				
				//JetIndexMap jetIdmap;	//jet index map: (Zboson, (mu1_jetid, mu2_jetid))
				//JetIndexMap Zcands_map;	//jet index map: (Zboson, (mu1_jetid, mu2_jetid))

				vector<PtEtaPhiMVector> Zcands_vec;
				JetIndexMap jetId_map;	//jet index map: (Zboson_index, (mu1_jetid, mu2_jetid))

				//n_mu is different for MC and data
				//loop over all the muons in the event --> check all possible combinations (mu_j, mu_k) for "good Z candidates" (double loop)
				//for(int mu1_ind=0; mu1_ind!=(n_muprompt-1); ++mu1_ind){
					//for(int mu2_ind=(mu1_ind+1); mu2_ind!=n_muprompt; ++mu2_ind){
				for(int mu1_ind=0; mu1_ind!=(n_mu-1); ++mu1_ind){
					for(int mu2_ind=(mu1_ind+1); mu2_ind!=n_mu; ++mu2_ind){

						Muon *muon1 = muons_ptr[mu1_ind];
						Muon *muon2 = muons_ptr[mu2_ind];

						//need to anyway keep track of jetindex...
						Int_t mu1_jetindex = muon1->jet_index();
						Int_t mu2_jetindex = muon2->jet_index();

						PtEtaPhiMVector Zcand = Zreco(muon1, muon2);
						int Zcand_charge = muon1->charge()+muon2->charge();	//use this variable for better readability
						//for isolation --> need to use the corresponding indices!!!!!!!!!!
						bool good_isolation = (Muon_pfIsoId[muon1->mu_index()]>1)&&(Muon_pfIsoId[muon2->mu_index()]>1);
						//bool good_isolation = (Muon_pfIsoId[0]>1)&&(Muon_pfIsoId[1]>1);	//check muon isolation --> 2=loose (for both, goes from 0 to 6 --> 'very very loose' to 'very very tight'
						bool good_eta =	(TMath::Abs(muon1->eta())<2.4)&&(TMath::Abs(muon2->eta())<2.4);	//check that both muons have abs(Eta)<2.4

						//pf-isolation only for data or also for MC?
						//test
						/*
						if(((TMath::Abs(Zcand.M()-91.2))<20)&&(Zcand_charge==0)){
							cout << "Zcand.M(): " << Zcand.M() << endl;
							cout << "mu1iso: " << (unsigned int) Muon_pfIsoId[muon1->mu_index()] << " index0: " << (unsigned int) Muon_pfIsoId[0] << endl;
							cout << "mu2iso: " << (unsigned int) Muon_pfIsoId[muon2->mu_index()] << " index1: " << (unsigned int) Muon_pfIsoId[1] << endl;
							cout << "mu1eta: " << muon1->eta() << endl;
							cout << "Muon Index mu1: " << muon1->mu_index() << endl;
						}
						*/

						//check if charge is correct and muon pT-selection and Zmass-window AND CHECK ISOLATION! (shall additional counters be created? --> not for now.)
						//bool isgoodZ = (Zcand_charge==0)&&(min(muon1->pt(), muon2->pt())>13)&&(max(muon1->pt(), muon2->pt())>24)&&((TMath::Abs(Zcand.M()-91.2))<20)&&good_isolation&&good_eta;

						//introduce cut on Z-pt !!! --> analoguous to DijetAnalysis --> pT_Z > 30 GeV
						bool isgoodZ = (Zcand.pt()>30)&&(Zcand_charge==0)&&(min(muon1->pt(), muon2->pt())>13)&&(max(muon1->pt(), muon2->pt())>24)&&((TMath::Abs(Zcand.M()-91.2))<20)&&good_isolation&&good_eta;
						//just counter:
						if(Zcand.pt()<=30){count_badZpt++;}

						if(isgoodZ){//only reconstruct Z when it is a good candidate?
							///PtEtaPhiMVector Zcand = Zreco(muons_ptr[mu1_ind], muons_ptr[mu2_ind]);
							////goodZs.insert(Zcand);
							//need to anyway keep track of jetindex...
							Int_t mu1_jetindex = muons_ptr[mu1_ind]->jet_index();
							Int_t mu2_jetindex = muons_ptr[mu2_ind]->jet_index();

							Zcands_vec.push_back(Zcand);
							int Z_idx = Zcands_vec.size()-1;
							jetId_map.insert(make_pair(Z_idx, make_pair(mu1_jetindex, mu2_jetindex)));
							//Zcands_map.insert(make_pair(Zcand, make_pair(mu1_jetindex, mu2_jetindex)));	//first attempt..
						}
						else{continue;}
					}//inner loop through prompt muons
				}//outer loop through prompt muons


				//should there be a separate vector storing the Zcands?
				//go on with the analysis of the Z-candidates --> check for Z+jet events
				int nZcands = Zcands_vec.size();	//number of good Zcands
				if(nZcands>0){	//have found at least one Z-candidate
					//for now: only use events where EXACTLY one good Z candidate has been found
					if(nZcands==1){
						//take reconstructed Z candidate as Zboson of this event
						int Z_idx = 0;
						PtEtaPhiMVector Zboson = Zcands_vec.at(Z_idx);

						//check if Zcand has good rapidity
						bool good_Zrapidity = TMath::Abs(Zboson.Rapidity())<2.4;

						if(good_Zrapidity){
							//jet index for each muon
							Int_t mu1_jetindex = jetId_map[Z_idx].first;
							Int_t mu2_jetindex = jetId_map[Z_idx].second;

							//now --> find the leading jet which Index != mu1_jetindex and != mu2_jetindex
							// for this implement and use jetreco(int mu1_jetindex, int mu2_jetindex) (returns PtEtaPhiMVector jet)
							Double_t jet_pt;
							Double_t jet_eta;
							Double_t jet_phi;
							Double_t jet_mass;
							int jet_jetindex;	//for storing leading jet's index
							for(int k=0; k!=nJet; ++k){	//loop for obtaining properties of leading jet
								if(mu1_jetindex!=k && mu2_jetindex!=k){
									jet_pt = Jet_pt[k];	
									jet_eta = Jet_eta[k];
									jet_phi = Jet_phi[k];
									jet_mass = Jet_mass[k];
									jet_jetindex = k;				//store jet-index of leading jet!
									break;	//escape inner for-loop as we found the leading jet
								}
								else{continue;}
							}//end loop for jet-properties

							//4-vector for leading jet
							PtEtaPhiMVector jet(jet_pt, jet_eta, jet_phi, jet_mass);

							//4-vector for leading jet ---> does not work yet with that function... need to somehow pass the Float_t arrays (Branches)
							//PtEtaPhiMVector jet = jetreco(mu1_jetindex, mu2_jetindex, Jet_pt, Jet_eta, Jet_phi, Jet_mass);

							//now selection cuts on the leading jet! pT and rapidity	(how to handle triggers?)
							//bool goodJet = (jet.pt()>30)&&(TMath::Abs(jet.Rapidity())<2.4);	//adjusted rapidity cut
							bool goodJet = (jet.pt()>20)&&(TMath::Abs(jet.Rapidity())<2.4);	//adjusted rapidity cut //(21.01.)--> set down pt cut, will be re-checked before filling though


							if(goodJet){
								//should filling of zpt_vec etc.  become part of fillhistos() ?? --> YES. (ybys_pointer is anyway used there)
								// make ratio mu1pt/mu2pt --> will be passed on to fillhistos()
								// ratio mu1pt/mu2pt:
								Double_t mupt_ratio = (muons_ptr[0]->pt())/(muons_ptr[1]->pt());
								// or simply give mu1pt and mu2pt as Double_t to fillhistos()
								
								
								//if the Z+jet event passes the muon_trigger is already checked before entering the current scope
								//histo filling function:
								//void fillhistos(PtEtaPhiMVector Zboson, PtEtaPhiMVector jet, YbYs* ybys_pointer, vector<int> triggers, vector<Histos*> histosvec, 
													//vector<bool> cur_trigdec, int jentry, int &count_zjet, int &count_NoBin_ybys, Double_t mupt_ratio)

								if(montecarlo==false){//if DATA
									if(jet.pt()>30){	//stricter jet.pt() cut --> check if jet fulfills also the desired requirement of minimum_pT 30 GeV
										fillhistos(Zboson, jet, ybys_pointer, triggers, histosvec, cur_trigdec, jentry, count_zjet, count_NoBin_ybys, mupt_ratio);
									}
									else{
										count_20pt30++;
									}
								}
								else if(montecarlo==true){	//if MONTE CARLO //if MC call fillhistos without trigger-options
									//void fillhistos(PtEtaPhiMVector Zboson, PtEtaPhiMVector jet, YbYs* ybys_pointer, int jentry, int &count_zjet, int &count_NoBin_ybys, Double_t mupt_ratio, Double_t gen_weight)
									Double_t gen_weight = Generator_weight;		//MC event weight to be given to fillhistos() function
									//fillhistos(Zboson, jet, ybys_pointer, jentry, count_zjet, count_NoBin_ybys, mupt_ratio, gen_weight);

									//might contain unmatched leading jets! (matching only checked later)
									if(strcmp(generator, "herwig")==0){
										//call modified fill-histos --> call it with smeared jet (if pt is > 30 will be checked within fillhistos() here)
										PtEtaPhiMVector smearedJet(Jet_pt_nom[jet_jetindex], jet_eta, jet_phi, jet_mass);
										fillhistos(Zboson, smearedJet, ybys_pointer, jentry, count_zjet, count_NoBin_ybys, mupt_ratio, gen_weight, uncnames, Jet_pt_uncUp, Jet_pt_uncDw, jet_jetindex);
										if(jet.pt()>30){//also fill the unsmeared jet (no unc-calculation)
											fillhistos(Zboson, jet, ybys_pointer, jentry, count_zjet, count_NoBin_ybys, mupt_ratio, gen_weight);
										}
										else{
											count_20pt30++;
										}
									}
									else{	//normal call if pythia
										if(jet.pt()>30){	//check if desired condition is met
											fillhistos(Zboson, jet, ybys_pointer, jentry, count_zjet, count_NoBin_ybys, mupt_ratio, gen_weight);
										}
										else{
											count_20pt30++;
										}
									}


									/*
									//This event has passed event selection --> do uncertainty analysis (if herwig)
									if(strcmp(generator, "herwig")==0){
										//need to "loop through the branches" of the Friends TTree -- don't have to loop. just take the values. (from the map)
										//The branch usually has the same name as the Unc-histo, just added an "Up" or "Down"
										for(unsigned int k=0; k!=uncnames.size(); ++k){
											if(jentry%10000==0){
												cout << "k==" << k << endl;
												cout << "unc = " << uncnames.at(k) << endl;
											}

											//Get the corresponding branches
											//these already contain the varied Jet_pt value
											const char* unc_up_branchname = Form("Jet_pt_jes%sUp", uncnames.at(k).c_str());
											const char* unc_down_branchname = Form("Jet_pt_jes%sDown", uncnames.at(k).c_str());

											//Get the variations -- only vary leading jet (Z stays unchanged)
											//varied leading jet pt
											double newpt1_up = Jet_pt_uncUp[k][0];
											double newpt1_down = Jet_pt_uncDw[k][0];

											//unchanged Zboson pt:
											double Zpt = Zboson.Pt();

											//new ptavg
											double newptavg_up = 0.5*(newpt1_up+Zpt);
											double newptavg_down = 0.5*(newpt1_down+Zpt);

											if(jentry%10000==0){		//testing
												cout << "old pt1: " << jet.pt() << endl;	
												cout << "new pt1 up: " << newpt1_up << endl;
												cout << "new pt1 down: " << newpt1_down << endl;
											}

											//fill corresponding histogram
											//fill the histograms (could also write function instead)
											cur_ybysbin.unc_histos_up[uncnames.at(k)]->Fill(newptavg_up, Generator_weight);
											cur_ybysbin.unc_histos_down[uncnames.at(k)]->Fill(newptavg_down, Generator_weight);

										}//unc sources loop

										//now fill as well Jet_pt_raw, Jet_pt_nom from uncfile (additional information)

									}//if mc HERWIG (end of unc-loop)
									*/


//-------------the following is still experimental -------------------//
									// JET MATCHING
									//check if leading jet has a matched gen jet (if yes, fill TH2D afterwards)
									bool matched = false;
									//the three observables in gen
									Double_t gen_ptavg;
									Double_t gen_yboost;
									Double_t gen_ystar;
									//more helper variables in gen (rapidity) and Lorentzvectors
									Double_t genjet_y;				//in dijet called: gen1_y
									PtEtaPhiMVector genjet_vec;
									//other bools on gen
									bool good_genrapidity;
									bool good_genpt;
									Int_t genjet1_ind=-1;		//index of leading genjet

									//go through all the genjets and find one that matches with leading jet
									for(int gen_ind=0;gen_ind!=nGenJet;++gen_ind){
										//Double_t genjet_eta = GenJet_eta[gen_ind];
										//Double_t genjet_phi = GenJet_phi[gen_ind];
										Double_t del_Eta_sqr = TMath::Power(TMath::Abs(jet.Eta()-GenJet_eta[gen_ind]), 2);	//squared
										Double_t del_Phi_sqr = TMath::Power(TMath::Abs(jet.Phi()-GenJet_phi[gen_ind]), 2);

										Double_t delR_j = TMath::Sqrt(del_Eta_sqr+del_Phi_sqr);	//check current genjet with jet1
										if(delR_j<0.2){ //MATCHED!
											genjet1_ind = gen_ind;	
											matched = true;
											genjet_vec = PtEtaPhiMVector(GenJet_pt[genjet1_ind], GenJet_eta[genjet1_ind], GenJet_phi[genjet1_ind], GenJet_mass[genjet1_ind]);
											genjet_y = genjet_vec.Rapidity();

											gen_ptavg = 0.5*(GenJet_pt[genjet1_ind]+Zboson.Pt());
											gen_yboost = 0.5*TMath::Abs(genjet_y+Zboson.Rapidity());
											gen_ystar = 0.5*TMath::Abs(genjet_y-Zboson.Rapidity());
											good_genrapidity = (TMath::Abs(genjet_y)<2.4)&&(good_Zrapidity);
											good_genpt = (genjet_vec.pt()>30);	//here already clear that Zpt>30GeV

											count_jetmatch++;	//matched leading recojet to a genjet!
										}//end if(delR_j1<0.2)
										else if(gen_ind==nGenJet-1){
											count_noGen1++;//leading jet has no matching gen jet...
										}
									}//end of genjet loop


									// to do (under construction, taken from dijet)
									if(matched){//if leading jet has a gen-match
										Double_t cur_ptavg = 0.5*(Zboson.Pt()+jet.Pt());


										//loop through ybys bins
										for (int ybysnr=0; ybysnr!=6; ++ybysnr){
											YbYs cur_ybysbin = ybys_pointer[ybysnr];
											float yb_lo = cur_ybysbin.yb_get_bb().first;
											float yb_up = cur_ybysbin.yb_get_bb().second;
											float ys_lo = cur_ybysbin.ys_get_bb().first;
											float ys_up = cur_ybysbin.ys_get_bb().second;
						
											///cur_ybysbin.get_mc_genreco_noWeights_noGenCuts()->Fill(gen_ptavg, cur_ptavg);	//fill without any cuts on gen jets

											//fill AND APPLY THE WEIGHTS!
											cur_ybysbin.get_mc_genreco_noGenCuts()->Fill(gen_ptavg, cur_ptavg, gen_weight);	//fill without any cuts on gen jets, but applied WEIGHTS


										//matched gen-spectrum? (leading genjet and Zboson)
			/*
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
			*/


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
														cur_ybysbin.get_mc_profile_mergedj1j2()->Fill(gen_ptavg, (jet.pt())/(gen_ptavg), gen_weight); //filling for leading jet
														cur_ybysbin.get_mc_profile_mergedj1j2()->Fill(gen_ptavg, (Zboson.pt())/(gen_ptavg), gen_weight); //filling for Z boson
													}
													if(TMath::Abs((jet.pt())/(genjet_vec.pt())-1)<0.5){
														cur_ybysbin.get_mc_profile_j1()->Fill(genjet_vec.pt(), (jet.pt())/(genjet_vec.pt()), gen_weight);	//leading gen and reco jets
													}
												}//end: if gen ybys in current ybys bin
											}//end: if good gen-rapidity && good gen-pt
										}//end of ybys-loop for the genjet
								
									}//endif matched
	
//---------------------------end of experimental matching-part --------------------------//

								}//end if MC
							}//if good jet
							else{
								count_badJet += 1;	//leading jet does not fulfill requirements (=>need pT>30GeV and abs(y)<3.0)
							}
						}//if good Z-rapidity
						else{
							count_badZrapidity +=1;	//Z rapidity too high
						}
					}//if nZcands==1
					else{//more than one good Z candidate
						count_manyZcands += 1;
					}
				}//end if(Zcands_vec.size()>0)
				else{//have not found any possible (good) Z candidate in this event
					count_noZcands +=1;
				}
				

				/*
				//OLD CODE FOLLOWS:
				//check if leading muons have opposite sign of charge
				int Zcharge = muons_ptr[0]->charge()+muons_ptr[1]->charge();	//might not be the optimal procedure... ?
				//and if their invariant mass is about Zmass (see below)

				if(Zcharge==0){
					// find "leading muons" --> their invariant mass should be in some window around the Z-mass (CHECK THIS!!)
					//need to anyway keep track of jetindex...
					Int_t mu1_jetindex = muons_ptr[0]->jet_index();
					Int_t mu2_jetindex = muons_ptr[1]->jet_index();
					//reconstruct Zboson
					PtEtaPhiMVector Zboson = Zreco(muons_ptr[0], muons_ptr[1]);
					//check if reconstructed Zboson lies within mass window m_Z +- 15 GeV (for m_Z = 91 GeV)
					// 76GeV < Zmass < 106GeV
					Double_t Zmass = Zboson.M();
					if(76.0<Zmass && Zmass<106.0){
						//now --> find the leading jet which Index != mu1_jetindex and != mu2_jetindex
						// for this implement and use jetreco(int mu1_jetindex, int mu2_jetindex) (returns PtEtaPhiMVector jet)
						Double_t jet_pt;
						Double_t jet_eta;
						Double_t jet_phi;
						Double_t jet_mass;
						for(int k=0; k!=nJet; ++k){	//loop for obtaining properties of leading jet
							if(mu1_jetindex!=k && mu2_jetindex!=k){
								jet_pt = Jet_pt[k];	
								jet_eta = Jet_eta[k];
								jet_phi = Jet_phi[k];
								jet_mass = Jet_mass[k];
								break;	//escape inner for-loop as we found the leading jet
							}
							else{continue;}
						}//end loop for jet-properties

						//4-vector for leading jet
						PtEtaPhiMVector jet(jet_pt, jet_eta, jet_phi, jet_mass);

						//4-vector for leading jet ---> does not work yet with that function... need to somehow pass the Float_t arrays (Branches)
						//PtEtaPhiMVector jet = jetreco(mu1_jetindex, mu2_jetindex, Jet_pt, Jet_eta, Jet_phi, Jet_mass);

						//should filling of zpt_vec etc.  become part of fillhistos() ?? --> YES. (ybys_pointer is anyway used there)
						// make ratio mu1pt/mu2pt --> will be passed on to fillhistos()
						// ratio mu1pt/mu2pt:
						Double_t mupt_ratio = (muons_ptr[0]->pt())/(muons_ptr[1]->pt());
						// or simply give mu1pt and mu2pt as Double_t to fillhistos()
						
						
						//if the Z+jet event passes the muon_trigger is already checked before entering the current scope
						//histo filling function:
						//void fillhistos(PtEtaPhiMVector Zboson, PtEtaPhiMVector jet, YbYs* ybysbins_ptr, vector<int> triggers, vector<Histos*> histosvec)
						fillhistos(Zboson, jet, ybys_pointer, triggers, histosvec, cur_trigdec, jentry, count_zjet, count_NoBin_ybys, mupt_ratio);
					}//if Zmass okay
					else{//Zmass wrong
						if(v==true){cout << "Event yields to Z-boson with wrong mass." << endl;}
						count_wrongmass += 1;
					}//Zmass wrong
				}//if Zcharge==0
				else if(Zcharge!=0){	//if the event has 3 suitable muons, check if the third one can be matched
					if(v==true){cout << "Check possible third muon.... --> TO BE IMPLEMENTED!" << endl;}	//IMPLEMENT THIS CHECK!
					count_No12Z += 1;	//count: Z not reconstructed with first and second muon
					if(n_muprompt>2 && muons_ptr[0]->charge()+muons_ptr[2]->charge()==0){ //check first with third muon 
					//as 1st and 2nd muon have same charge, this works autmatically for 2nd and 3rd if it works for 1st and 3rd...

						// FOLLOWING PROCEDURE equals what has been done before --> could write a function, but have to pass arrays "Jet_pt" etc.

						//reconstruct Zboson with (1st and 3rd) and (2nd and 3rd) --> then: check mass.
						//reconstruct Zboson from 1st and 3rd muon
						PtEtaPhiMVector Zboson_13 = Zreco(muons_ptr[0], muons_ptr[2]);
						Double_t Zmass_13 = Zboson_13.M();

						//reconstruct Zboson from 2nd and 3rd muon
						PtEtaPhiMVector Zboson_23 = Zreco(muons_ptr[1], muons_ptr[2]);
						Double_t Zmass_23 = Zboson_23.M();

						//which mass lies closer to Zmass?
						Double_t Mdiff_13 = TMath::Abs(Zmass_13-91.2);
						Double_t Mdiff_23 = TMath::Abs(Zmass_23-91.2);

						PtEtaPhiMVector Zboson;	//prepare Zboson
						Int_t mu1_jetindex;		//declare jet index of muon1 of Z
						Int_t mu2_jetindex;		//declare jet index of muon2 of Z
						Double_t mupt_ratio;
						if(Mdiff_13<Mdiff_23){
							//need to anyway keep track of jetindex
							mu1_jetindex = muons_ptr[0]->jet_index();
							mu2_jetindex = muons_ptr[2]->jet_index();
							Zboson = Zboson_13;
							mupt_ratio = (muons_ptr[0]->pt())/(muons_ptr[2]->pt());	//will only be used if in Zmass window -> in fillhistos()
							if(Mdiff_13<15.0){ if(v==true){cout << Form("Event: %lld --> Combination (mu1, mu3) successful.", jentry) << endl;}}
						}
						else{
							mu1_jetindex = muons_ptr[1]->jet_index();
							mu2_jetindex = muons_ptr[2]->jet_index();
							Zboson = Zboson_23;
							mupt_ratio = (muons_ptr[1]->pt())/(muons_ptr[2]->pt());	//will only be used if in Zmass window -> in fillhistos()
							if(Mdiff_23<15.0){ if(v==true){cout << Form("Event: %lld --> Combination (mu2, mu3) successful.", jentry) << endl;}}
						}

						//check if reconstructed Zboson lies within mass window m_Z +- 15 GeV (for m_Z = 91 GeV)
						// 76GeV < Zmass < 106GeV
						Double_t Zmass = Zboson.M();

						if(76.0<Zmass && Zmass<106.0){
							if(v==true){ //a little verbose
								cout << Form("Mass of Zboson:  %f", Zmass) << endl;
								cout << Form("Mdiff_13: %f || Mdiff_23: %f", Mdiff_13, Mdiff_23) << endl << endl;
							}
							n_success_3rd += 1;

							//now --> find the leading jet which Index != mu1_jetindex and != mu2_jetindex
							// for this implement and use jetreco(int mu1_jetindex, int mu2_jetindex) (returns PtEtaPhiMVector jet)
							Double_t jet_pt;
							Double_t jet_eta;
							Double_t jet_phi;
							Double_t jet_mass;
							for(int k=0; k!=nJet; ++k){	//loop for obtaining properties of leading jet
								if(mu1_jetindex!=k && mu2_jetindex!=k){
									jet_pt = Jet_pt[k];	
									jet_eta = Jet_eta[k];
									jet_phi = Jet_phi[k];
									jet_mass = Jet_mass[k];
									break;	//escape inner for-loop as we found the leading jet
								}
								else{continue;}
							}//end loop for jet-properties

							//4-vector for leading jet
							PtEtaPhiMVector jet(jet_pt, jet_eta, jet_phi, jet_mass);

							//4-vector for leading jet ---> does not work yet with that function... need to somehow pass the Float_t arrays (Branches)
							//PtEtaPhiMVector jet = jetreco(mu1_jetindex, mu2_jetindex, Jet_pt, Jet_eta, Jet_phi, Jet_mass);

							if(HLT_Mu18_Mu9){	//event passes muon pT trigger
								//histo filling function:
								//void fillhistos(PtEtaPhiMVector Zboson, PtEtaPhiMVector jet, YbYs* ybysbins_ptr, vector<int> triggers, vector<Histos*> histosvec)
								fillhistos(Zboson,jet, ybys_pointer, triggers, histosvec, cur_trigdec, jentry, count_zjet, count_NoBin_ybys, mupt_ratio);
							}
							else{	//event does not pass muon pT trigger
								count_dropout_MuonTrigger += 1;
							}
						}//end if correct Zmass
						else{//Zmass wrong
							//cout << "combination of 1st and 3rd muon leads to wrong mass." << endl;
							count_wrongmass += 1;
						}//Zmass wrong
					}//if 1st and 3rd muon lead to correct charge 

					//the following is always true if the first "if" here is true...
					if(n_muprompt>2 && muons_ptr[1]->charge()+muons_ptr[2]->charge()==0){ //check second with third muon
					}

				}//if third muon has to be used (because Zcharge!=0)
				*/

			}//if n_muprompt>1 AND muon_trigger==true
			//else if((muon_trigger==false)||(n_muprompt<2)){
			else if((montecarlo==false)&&((muon_trigger==false)||(n_muprompt<2))){
				count_fail_prompt_OR_trg +=1;	//event dropped out because of Muon_mediumPromptId or HLT_Mu23_Mu12 (or both!)
				if(muon_trigger==false){
					count_dropout_MuonTrigger +=1;
				}
			}
		}//if (nMuon>1)
		else{
			count_less2mu += 1;	//count events with less than two muons
			//count_jetless += 1;	//does not contain jet
		}//no Z, because less than 2 muons
	}//loop through events

	if(montecarlo==false){
		//### SUMMARY ###
		//normalise histogram --> divide filled histograms by luminosity
		//double lumi_mutrg = 13.004544791;	//obtained with brilcalc for HLT_Mu23_Mu12 trigger (in /pb)
		for (int ybysnr=0; ybysnr!=6; ++ybysnr){
			YbYs cur_ybysbin = ybys_pointer[ybysnr];
			const char* cur_bname = cur_ybysbin.get_bname().c_str();
			TDirectory *extraDir = standard_folder->GetDirectory(Form("%s/additional_plots_%s", cur_bname, cur_bname));			//TO BE CHANGED --> put it in normal ybysdir ..not in "additional"
			extraDir->cd();
			TH1D* sumhisto = (TH1D*)(histosvec.at(ybysnr)->jt_histos[triggers.at(0)])->Clone(Form("sumhisto_%s", cur_bname));
			sumhisto->GetXaxis()->SetTitle("p_{T,avg} /GeV");
			sumhisto->GetYaxis()->SetTitle("N/lumi");	//actually: sigma^3/(yb ys ptav) in pb/GeV
			sumhisto->Scale(1./lumi_mutrg);
			standard_folder->cd();
		}
	}



	//summary: mu1/mu2 over Zpt --> could write this into a function. (now just testing)
	//write it as member-function!
	if(muratio_mode){
		for (int ybysnr=0; ybysnr!=6; ++ybysnr){
			YbYs cur_ybysbin = ybys_pointer[ybysnr];	//this is copying the object pointed to by ybys_pointer[ybysnr] (including corresp. pointers)
			///YbYs* cur_ybysbin = &ybys_pointer[ybysnr];	//this is now just renaming the pointer!

			///YbYs cur_ybysbin = ybys_pointer.at(ybysnr);
			//cur_ybysbin->pushpoint(23.2);
			//cur_ybysbin->pushpoint(26.2);

			//ybys_pointer[ybysnr].pushpoint(5.3);
			//ybys_pointer[ybysnr].pushpoint(9.3);

			//TGraph* curgraph = cur_ybysbin.muratio_graph();
			TGraph* curgraph = ybys_pointer[ybysnr].muratio_graph();

			//cur_ybysbin->muratio_graph();
			//if(ybysnr>0){ybys_pointer[ybysnr-1].muratio_graph(); cout << "test." << endl;}
			curgraph->Write();
		}
	}

	/*
	for (int ybysnr=0; ybysnr!=6; ++ybysnr){
		YbYs cur_ybysbin = ybys_pointer[ybysnr];
		int npoints = cur_ybysbin.get_zpt_vec()->size();
		Double_t xvals[npoints], yvals[npoints];	//make arrays for x- and y-values of TGraph
		cout << Form("In %s found %d Z+jet events.", cur_ybysbin.get_bname().c_str(), npoints) << endl;
		for(int k=0; k!=npoints; ++k){
			xvals[k] = cur_ybysbin.get_zpt_vec()->at(k);
			yvals[k] = cur_ybysbin.get_mupt_ratio_vec()->at(k);

		}
		//TGraph *curgraph = new TGraph(npoints, xvals, yvals);
		TGraph *curgraph = cur_ybysbin.get_graph();
		curgraph->SetDirectory()
		curgraph->Draw("AC*");
	}
	*/



			//------ RESPONSE MATRIX STUFF COPIED FROM DIJET CODE --> implement! ------//
	/*
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
	*/






	interim_output->Write();

	//cout << endl << "Amount of events with (nJet>0) and (Abs(Jet_eta[0])>4.7): " << count_higheta << endl;
	//cout << "Amount of events with (nJet==0): " << count_jetless << endl << endl;
	//cout << "Amount of events that are no dijet events (nJet<2): " << count_nodijet << endl;
	cout << endl;
	cout << "		Analysis SUMMARY:			" << endl;
	cout << "------------------------------------------------------------------" << endl;
	cout << "Amount of events with less than 2 muons: " << count_less2mu << endl << endl;
	//cout << "From all the events with at least 2 muons, this many muons have not passed the MediumPromptId selection: " << count_NoMediumPromptId << endl;
	cout << "Number of Events with nMuon>1 that dropped out due to MediumPromptId selection: " << count_dropout_MediumPromptId << endl;
	//cout << "(Q_mu1+Q_mu2 != 0) --> Z cannot be reconstructed from 1st and 2nd muon in this many cases: " <<  count_No12Z << endl;
	//cout << "Events that remain after previous cuts, but lead to Zmass outside selected range: " << count_wrongmass << endl;
	///cout << "Events that dropped out due to muon pT trigger: " << count_dropout_MuonTrigger << endl << endl;
	cout << "Events that dropped out due to Muon_mediumPromptId or HLT_Mu23_Mu12 (or both): " << count_fail_prompt_OR_trg << endl;
	cout << "Amount of Z+jet events that are not within chosen  ybys bins: " << count_NoBin_ybys << endl;
	cout << "Events, where rapidity of Zboson was too high (|y|>=2.4): " << count_badZrapidity << endl;
	cout << "(Events, where pT of Zboson was too low (pT<=30GeV): " << count_badZpt << " --> does not have to be only dropout reason.)" <<  endl;
	cout << "Events, where leading jet did not fulfill criteria: " << count_badJet << endl;
	cout << "Events that dropped out late, passing everything, but 20GeV < pTjet < 30GeV : " << count_20pt30 << endl << endl;

	cout << "--> Amount of Z+jet events that have passed selection criteria: " << count_zjet << endl << endl;;

	//just for testing purposes
	cout << "Events without any muons: " << n_0mu << endl;
	cout << "Events containing one muon: " << n_1mu << endl;
	cout << "Events containing two muons: " << n_2mu << endl;
	cout << "Events containing three muons: " << n_3mu << endl;
	cout << "Events containing four muons: " << n_4mu << endl;
	cout << "Events containing five or more muons: " << n_5mu << endl;
	cout << "--------------------------------------------------" << endl;
	cout << "events with 2 prompt muons: " << prompt2 << endl;
	cout << "events with 3 prompt muons: " << prompt3 << endl;
	cout << "events with 4 prompt muons: " << prompt4 << endl;
	cout << "events with 5 or more prompt muons: " << prompt5 << endl << endl;

	//3rd muon tests
	//cout << "Amount of events where 3rd muon brought success (with 1st or 2nd): " << n_success_3rd << endl << endl;
	cout << "Amount of events (nMuon>1, n_muprompt>1) that passed muon_trigger, but where no good Z candidate has been found: " << count_noZcands << endl;
	cout << "Amount of events, where more than one good Z candidate has been found: " << count_manyZcands << endl << endl;


} //end of Loop()





//--------------------------------------------------------------------------------------------------//
//											 CREATE JSON FILE										//
//--------------------------------------------------------------------------------------------------//
//new attempt to create function which helps with finding the lumisections and runnumbers for each event
//originally planned to implement this in seperate file --> see get_lumi_txt.C
//copied from my DijetAnalysis.C
//here not used so far.
void ZjetAnalysis::get_lumi_txt(){
// function for creating a JSON file via reading "run" and "LuminosityBlock" from nanoAOD file.
// for later use with brilcalc
	TH1::SetDefaultSumw2(kTRUE);      //in order to ALWAYS store the bin errors
	setTDRStyle();
	fChain->SetBranchStatus("*",0);	//switch off all the branches
	fChain->SetBranchStatus("run",1);
	fChain->SetBranchStatus("luminosityBlock",1);

	TFile* json_output = new TFile("run_lumisec.txt", "RECREATE");
	json_output->cd();


	ofstream txt_file;
	txt_file.open("txt_file_test.txt");

	//check runperiod and create json file (this function is only called for data anyway --> no if(data) needed)
	//check which runperiod
	ofstream json_file;
	if((strcmp(runperiod,"A")==0)or(strcmp(runperiod,"B")==0)or(strcmp(runperiod,"C")==0)or(strcmp(runperiod,"D")==0)){
		cout << "Run " << runperiod  << endl;
		json_file.open(Form("json_file_zjet_Run%s.txt", runperiod));
	}
	else{//if no runperiod specified
		json_file.open("json_file_zjet.txt");
	}


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




//function for normalising the MC histograms
//taken from Dijet
void ZjetAnalysis::normalise_mc(){
	TH1::SetDefaultSumw2(kTRUE);      //in order to ALWAYS store the bin errors
	setTDRStyle();

	fChain->SetBranchStatus("Generator_weight", 1);	//generator weight
	Long64_t nentries = fChain->GetEntriesFast();	//get no. of entries

	//vectors containing the style for each ybys bin (as later used in the RatioAnalysis.C script)
    vector<Color_t> cols = {kYellow+1, kGray+2, kGreen+2, kBlue-2, kOrange-3, kRed+2, kAzure+1, kViolet+6, kSpring+5, kViolet-5, kCyan-3};
    vector<Style_t> markers = {kOpenCircle, kOpenTriangleUp, kOpenSquare, kOpenDiamond, kOpenCross, kOpenTriangleDown};


	//open input MC file and create output file
	const char* inputfile_name;
	const char* mc_inputhist;				//name of the histogram to be normalised
	TFile* normalised_mc_output;			//in dijet depends on generator, here look only at pythia

	//prepare normalisation factors
	Int_t nevents_tot;		//total number of generated MC events
	Double_t xs_mc_sample;	//cross section of the mc sample used (in pb)

	cout << endl;

	//Actions according to generator that has been used (Z+jet: added herwig7 sample, as used for uncs, but XS on XSDB currently only available for pythia8)
	//-------------------------------------------------
	//choose input file containing the histograms
	//create output file with proper name
	//set the normalisation factor

	if(strcmp(generator, "herwig")==0){
		inputfile_name = "zjet_interim_output_mc_herwig7.root";
		mc_inputhist = "mc_histo_smearedJet";					//if herwig --> take smeared Jet

		normalised_mc_output = new TFile("zjet_normalised_output_mc_herwig7.root", "RECREATE");
		//nevents_tot = 19464000;		//now assigned via mc_genweight histo
		xs_mc_sample = 6077.22;			// in pb, currently taken from the pythia flat sample, as not found for herwig sample

		cout << Form("[ZjetAnalysis::normalise_mc()]: Looking at a Herwig7 sample.") << endl;
		cout << "[ZjetAnalysis::normalise_mc()]: Cross Section of MC sample (used from pythia8 flat sample): " << xs_mc_sample << " pb"<< endl;
	}
	else if(strcmp(generator, "pythia")==0){
		inputfile_name = "zjet_interim_output_mc_pythia8.root";
		mc_inputhist = "mc_histo";

		normalised_mc_output = new TFile("zjet_normalised_output_mc_pythia8.root", "RECREATE");
		//nevents_tot = 19708000;		//now assigned via mc_genweight histo
		xs_mc_sample = 6077.22; 		//in pb (von XSDB)

		cout << "[ZjetAnalysis::normalise_mc()]: Looking at a Pythia8 sample." << endl;
		cout << "[ZjetAnalysis::normalise_mc()]: Cross Section of this MC sample: " << xs_mc_sample << " pb" << endl;
	}
	else{
		inputfile_name = "zjet_interim_output_mc.root";
		normalised_mc_output = new TFile("zjet_normalised_output_mc.root", "RECREATE");
		cout << "ERROR: No generator specified in ZjetAnalysis.h !!! -- Aborting." << endl;
		exit(EXIT_FAILURE);
	}


	/*
	//simpler version in Z+jet (as only pythia in use), above: see how it is done otherwise
	inputfile_name = "zjet_interim_output_mc_pythia8.root";
	normalised_mc_output = new TFile("zjet_normalised_output_mc_pythia8.root", "RECREATE");
	//nevents_tot = 19708000;		//now assigned via mc_genweight histo
	xs_mc_sample = 6077.22;		//in pb?? (NNLO von XSDB --> check all, there are also LO) -- this one is the updated one in twiki
	//xs_mc_sample = 5343.00;			//in pb? (eine LO von XSDB, is probably the older one?)
	//xs_mc_sample = 1370000000.0; 	//in pb? (von XSDB)
	cout << "[ZjetAnalysis::normalise_mc()]: Cross Section of this MC sample: " << xs_mc_sample << " pb" << endl;
	*/


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
	cout << "[ZjetAnalysis::normalise_mc()]: Total number of generated events: " << mc_genweight->GetEntries() << endl;


	//create the "ybXysX" folders in outputfile
	YbYs *ybys_pointer;	//contains way more than needed here
	//new outermost bin bounds --> go only up to 2.4 for both yb, ys
	ybys_pointer = new YbYs[6]{{0.0,0.0,1.0,1.0}, {0.0, 1.0, 1.0, 2.0}, {0.0, 2.0, 1.0, 2.4}, {1.0, 0.0, 2.0, 1.0}, {1.0, 1.0, 2.0, 2.0}, {2.0, 0.0, 2.4, 1.0}}; //new!!!
	for (int ybysnr=0; ybysnr!=6; ++ybysnr){
		YbYs cur_ybysbin = ybys_pointer[ybysnr];
        const char* ybys_bname = cur_ybysbin.get_bname().c_str();				//get current ybys bin name
		TDirectory* ybys_Dir = cur_ybysbin.ybys_directories(standard_folder);	//create ybys-directory
		///cur_ybysbin.ybys_directories(standard_folder);	//create ybys-directory

		const char *extradir_name = Form("additional_plots_%s", ybys_bname);
		ybys_Dir->mkdir(extradir_name);
		TDirectory *extraDir = ybys_Dir->GetDirectory(extradir_name);
		extraDir->cd();

		//Get the all the MC histos for zjet (this differs from the ones in dijet)
		//will then additionally add correctly normalised ones
		TH1D* cur_recohisto;	//histo containing reco distribution before normalising
		//TH1D* cur_genhisto;		//histo containing gen distribution (same selection as on reco, using GenJet instead of Jet Branches)
		//TH2D* cur_response_noGenCuts;	//response matrix (reco, gen in ptavg) without any gen cuts applied
		//TH2D* cur_norm_resp_noGenCuts;	//response matrix (no gen cuts) normalised to number of generated jets
		//TH2D* cur_response;				//response matrix, including same cuts on reco and gen
		//TH2D* cur_norm_resp;			//normalised response matrix

		//get the objects:
		//inputfile->GetObject(Form("Standard/%s/mc_histo_%s", ybys_bname, ybys_bname), cur_recohisto);						//reco histo
		inputfile->GetObject(Form("Standard/%s/%s_%s", ybys_bname, mc_inputhist, ybys_bname), cur_recohisto);				//reco histo (for H7: smeared Jet)

		
		/*	//these are so far only present in dijet --> but should be added here as well for unfolding
		inputfile->GetObject(Form("Standard/%s/mc_genhisto_%s", ybys_bname, ybys_bname), cur_genhisto);						//gen histo
		inputfile->GetObject(Form("Standard/%s/mc_genrec_noGenCuts_%s", ybys_bname, ybys_bname), cur_response_noGenCuts);	//resp matrix w/o gen cuts
		inputfile->GetObject(Form("Standard/%s/normalised_response_noGenCuts", ybys_bname), cur_norm_resp_noGenCuts);		//normalised resp matrix w/o gen cuts
		inputfile->GetObject(Form("Standard/%s/mc_gen_reco_ptavg_%s", ybys_bname, ybys_bname), cur_response);				//response matrix with gen and reco cuts
		inputfile->GetObject(Form("Standard/%s/normalised_response", ybys_bname), cur_norm_resp);							//response matrix with gen and reco cuts
		*/

		/*
		//append the ybys bin name to the normalised response matrices name:
		cur_norm_resp_noGenCuts->SetName(Form("%s_%s", cur_norm_resp_noGenCuts->GetName(), ybys_bname));
		cur_norm_resp->SetName(Form("%s_%s", cur_norm_resp->GetName(), ybys_bname));
		*/

		//Clone the histograms that shall be normalised
		TH1D* norm_recohisto = (TH1D*)cur_recohisto->Clone();
		///TH1D* norm_genhisto = (TH1D*)cur_genhisto->Clone();

		//change name of the unnormalised histograms:
		cur_recohisto->SetName(Form("%s_NotNormalised", cur_recohisto->GetName()));
		///cur_genhisto->SetName(Form("%s_NotNormalised", cur_genhisto->GetName()));

		//Normalise the histograms (factor is already set according to generator):
		//NEW: get the sum of weights directly from the inputfile!
		Double_t weights_mean = mc_genweight->GetMean();
		nevents_tot	= mc_genweight->GetEntries();
		Double_t sum_event_weights = nevents_tot*weights_mean;
		cout << "sum_event_weights: " << sum_event_weights << endl;	//test
		///Double_t sum_event_weights = 100108280;			//just for now as I had a mistake in mc_genweight

		//reco histo
		///norm_recohisto->Scale(1./nevents_tot);		//divide by total number of generated events
		norm_recohisto->Scale(1./sum_event_weights);
		norm_recohisto->Scale(xs_mc_sample);		//multiply by cross section of this mc sample
		norm_recohisto->Scale(1000);				//to get end value in fb (--> would pb be more suitable?)
		//norm_recohisto->Scale(58.83);				//this is the total luminosity of 2018 data in /fb ...not needed here.

		/*
		//same for gen histo
		///norm_genhisto->Scale(1./nevents_tot);
		norm_genhisto->Scale(1./sum_event_weights);
		norm_genhisto->Scale(xs_mc_sample);
		norm_genhisto->Scale(1000);
		//norm_genhisto->Scale(58.83);
		*/

		cout << Form("[ZjetAnalysis::normalise_mc()]: Normalised histograms in %s bin", ybys_bname) << endl;
		cout << "[ZjetAnalysis::normalise_mc()]: Storing to file..." << endl;
		//store all the histograms to the outputfile
		ybys_Dir->cd();

		SetDrawStyleWrite(norm_recohisto, "", norm_recohisto->GetTitle(), "p_{T,avg} /GeV", "XS in fb", cols.at(ybysnr+2), 2, markers.at(ybysnr), cols.at(ybysnr+2), 2, true, 0.9, true);
		norm_recohisto->GetYaxis()->SetTitleOffset(0.9);
		//norm_recohisto->Write();
		/*
		norm_genhisto->Write();
		cur_recohisto->Write();
		cur_genhisto->Write();
		cur_response_noGenCuts->Write();
		cur_norm_resp_noGenCuts->Write();
		cur_response->Write();
		cur_norm_resp->Write();
		*/


	 	//vectors containing the style for each ybys bin
    	//vector<Color_t> cols = {kYellow+1, kGray+2, kGreen+2, kBlue-2, kOrange-3, kRed+2, kAzure+1, kViolet+6, kSpring+5, kViolet-5, kCyan-3};
    	//vector<Style_t> markers = {kOpenCircle, kOpenTriangleUp, kOpenSquare, kOpenDiamond, kOpenCross, kOpenTriangleDown};

		///void SetDrawStyleWrite(TH1D *object, const char *drawcmd="", const char *title="", const char *xtitle="p_{T,avg}", const char *ytitle="", Color_t linecolor=kGray+2, Int_t linewidth=2, Style_t markerstyle=kFullCircle, Color_t markercolor=kGray+2, Int_t markersize=1, bool morelog=false, Double_t xTitleOffset=1.2, bool writing=false){
	
	}

	cout << endl;
	cout << "[ZjetAnalysis::normalise_mc()]: Results can be found in file: " << normalised_mc_output->GetName() << endl;

}//end function normalise_mc()



//function to divide uncUp by uncDown for each uncertainty source (taking the input file coming from the Loop function)
void ZjetAnalysis::uncs_UpDownDivide_mc(){
	TH1::SetDefaultSumw2(kTRUE);      //in order to ALWAYS store the bin errors/
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
	const char* inputfile_name = "zjet_interim_output_mc_herwig7.root";
	TFile *inputfile = TFile::Open(inputfile_name);												//open inputfile
	TFile *divided_uncs_mc_output = new TFile("zjet_divided_uncsources_mc_herwig7.root", "RECREATE");	//create output file

	cout << "[ZjetAnalysis::uncs_UpDownDivide_mc()]: Getting the Up Down uncsources from inputfile: " << inputfile->GetName() << endl;
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
	//ybys_pointer = new YbYs[6]{{0.0,0.0,1.0,1.0}, {0.0, 1.0, 1.0, 2.0}, {0.0, 2.0, 1.0, 2.4}, {1.0, 0.0, 2.0, 1.0}, {1.0, 1.0, 2.0, 2.0}, {2.0, 0.0, 2.4, 1.0}}; //new!!!
	ybys_pointer = new YbYs[6]{{0.0,0.0,1.0,1.0, uncnames}, {0.0, 1.0, 1.0, 2.0, uncnames}, {0.0, 2.0, 1.0, 2.4, uncnames}, {1.0, 0.0, 2.0, 1.0, uncnames}, {1.0, 1.0, 2.0, 2.0, uncnames}, {2.0, 0.0, 2.4, 1.0, uncnames}};	

	for (int ybysnr=0; ybysnr!=6; ++ybysnr){
		YbYs cur_ybysbin = ybys_pointer[ybysnr];
        const char* ybys_bname = cur_ybysbin.get_bname().c_str();				//get current ybys bin name
		TDirectory* ybys_Dir = cur_ybysbin.ybys_directories(standard_folder);	//create ybys-directory
		///cur_ybysbin.ybys_directories(standard_folder);	//create ybys-directory


		//loop through the uncertainty sources for creating the unc ratios
		for(unsigned int k=0; k!=uncnames.size(); ++k){
			cout << Form("[ZjetAnalysis::uncs_UpDownDivide_mc()]: Looking at unc source: %s", uncnames.at(k).c_str()) << endl;

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
		
			cout << "[ZjetAnalysis::uncs_UpDownDivide_mc()]: Storing to file..." << endl;
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
		cout << Form("[ZjetAnalysis::uncs_UpDownDivide_mc()]: Divided unc histograms in %s bin", ybys_bname) << endl;
		cout << "---------------------------------------------------------------------------------" << endl;
		cout << endl;

	}//end of ybys loop

	//is the following necessary? --No.
	//divided_uncs_mc_output->Write();

	cout << "[ZjetAnalysis::uncs_UpDownDivide_mc()]: Results can be found in file: " << divided_uncs_mc_output->GetName() << endl;

}//end of function uncs_UpDownDivide_mc()


//Function to divide each bin (in XS histograms) by binwidth, in order to end up with fb/GeV as XS unit
void ZjetAnalysis::DivideBinwidth(){
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
			inputfile_name = "zjet_interim_output_data_RunA.root";		//RunA
			binwidth_output = new TFile("zjet_xs_BinwidthDivided_data_RunA.root", "RECREATE");
			histname = "sumhisto_zjet_RunA";
		}//endif runA

		else if(strcmp(runperiod,"B")==0){
			cout << "Run B" << endl;
			inputfile_name = "zjet_interim_output_data_RunB.root";		//RunB
			binwidth_output = new TFile("zjet_xs_BinwidthDivided_data_RunB.root", "RECREATE");
			histname = "sumhisto_zjet_RunB";
		}//endif runB

		else if(strcmp(runperiod,"C")==0){
			cout << "Run C" << endl;
			inputfile_name = "zjet_interim_output_data_RunC.root";		//RunC
			binwidth_output = new TFile("zjet_xs_BinwidthDivided_data_RunC.root", "RECREATE");
			histname = "sumhisto_zjet_RunC";
		}//endif runC

		else if(strcmp(runperiod,"D")==0){
			cout << "Run D" << endl;
			inputfile_name = "zjet_interim_output_data_RunD.root";		//RunD
			binwidth_output = new TFile("zjet_xs_BinwidthDivided_data_RunD.root", "RECREATE");
			histname = "sumhisto_zjet_RunD";
		}//endif runD

		else{//if no runperiod specified
			cout << "ERROR: If this is data set -- specify run period." << endl;
			exit(EXIT_FAILURE);
		}
	}
	else{ //if MC
		if(strcmp(generator, "herwig")==0){
			inputfile_name = "zjet_normalised_output_mc_herwig7.root";	//herwig7
			binwidth_output = new TFile("zjet_xs_BinwidthDivided_mc_herwig7.root", "RECREATE"); //output file if herwig7
			mc_inputhist = "mc_histo_smearedJet";
			histname = "sumhisto_zjet_herwig7";
		}
		else if(strcmp(generator, "pythia")==0){
			inputfile_name = "zjet_normalised_output_mc_pythia8.root";	//pythia8
			binwidth_output = new TFile("zjet_xs_BinwidthDivided_mc_pythia8.root", "RECREATE"); //output file if pythia8
			mc_inputhist = "mc_histo";
			histname = "sumhisto_zjet_pythia8";
		}
		else{//if no generator specified
			cout << "ERROR: If this is mc set -- specify generator." << endl;
			exit(EXIT_FAILURE);
		}
	}

	TFile *inputfile = TFile::Open(inputfile_name);												//open inputfile

	cout << "[ZjetAnalysis::DivideBinwidth()]: Getting the XS histograms from inputfile: " << inputfile->GetName() << endl;
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
			inputfile->GetObject(Form("Standard/%s/additional_plots_%s/sumhisto_%s", ybys_bname, ybys_bname, ybys_bname), cur_xshisto);		//cross section histo to be normalised
		}
		else{
			inputfile->GetObject(Form("Standard/%s/%s_%s", ybys_bname, mc_inputhist, ybys_bname), cur_xshisto);									//cross section histo to be normalised
		}

		//if herwig, additionally the smearedJet histogram --> no. smearedJet histo is the default for herwig.
		//if(strcmp(generator, "herwig")==0){}
	
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

		cout << "[ZjetAnalysis::DivideBinwidth()]: Set new bin contents, divided by binwidth." << endl;
		cout << "[ZjetAnalysis::DivideBinwidth()]: Done with ybysbin: " << ybys_bname << endl;

		//store the histogram to the outputfile
		ybys_Dir->cd();
		cur_xshisto_bindiv->GetYaxis()->SetTitleOffset(0.9);
		cur_xshisto_bindiv->GetYaxis()->SetTitle("XS in fb/GeV");
		cur_xshisto_bindiv->SetLineColor(cols.at(ybysnr+2));
		cur_xshisto_bindiv->SetMarkerColor(cols.at(ybysnr+2));
		cur_xshisto_bindiv->SetMarkerStyle(markers.at(ybysnr));
		if(montecarlo==false){
			cur_xshisto_bindiv->SetTitle(Form("Z+jet XS Run %s", runperiod));
		}
		else{
			cur_xshisto_bindiv->SetTitle(Form("Z+jet XS %s", generator));
		}

		cur_xshisto_bindiv->Write();

		cout << Form("[ZjetAnalysis::DivideBinwidth()]: Writing to file.") << endl;
		cout << "---------------------------------------------------------------------------------" << endl;
		cout << endl;
	}//end of ybys loop

	cout << "[ZjetAnalysis::DivideBinwidth()]: Results can be found in file: " << binwidth_output->GetName() << endl;

}//end of function DivideBinwidth()



