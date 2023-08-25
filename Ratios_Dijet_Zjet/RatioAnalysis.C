//#######################################################################################################################//
//#																														#//
//#	This script takes the previously filled histograms form the Dijet-Analysis and the Z+jet-Analysis.					#//
//#	The intermediate results of the single analyses have been stored in ROOT files.										#//
//#	Here those ROOT files are processed further, carrying out the Ratio analysis.										#//
//#	Results are again written into a ROOT file, which can later be passed on to plotting scripts (yet to be created).	#//
//#																														#//
//#######################################################################################################################//

#define RatioAnalysis_cxx
#include "RatioAnalysis.h"
#include <string>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "tdrstyle_mod14.C" //use this for canvases etc.
#include <TMath.h>

using namespace std;

//verbose mode
bool v = true;

//runperiod -- has to be set manually
//use this variable also to indicate that it's MC

//const char* runperiod="A";
//const char* runperiod="B";
//const char* runperiod="C";
const char* runperiod="D";
//const char* runperiod="mc";

//const char* generator="pythia";
//const char* generator="herwig";
//const char* generator="HTpythia";
const char* generator="data";


//--------------------------------------------------------------------------------------------------//
//											 CLASSES 												//
//--------------------------------------------------------------------------------------------------//

//class for handling yboost ystar binning -- folder structure and ybys properties for "histogram-infos"
class YbYs {
	friend class Ratio;		//Ratio class is a friend --> it will use the ybys information
	private:
		string bname;			//bin name, e.g. "yb0ys0" for 0.0 <= yboost, ystar < 0.5, assigned in constructor
		float yb_low;			//lower yboost bin bound
		float yb_up;			//upper yboost bin bound
		float ys_low;			//lower ystar bin bound
		float ys_up;			//upper ystar bin bound

		//TH1D* ratiohisto_ptr;	//pointer to summary histogram

	public:
		YbYs (float, float, float, float);	//constructor with lower and upper bin bound for yboost ystar
		//TH1D* get_ratiohist(){return ratiohisto_ptr;};	//return pointer to the summary (ratio) histogram
		pair <float, float> yb_get_bb() {return make_pair(yb_low,  yb_up);};	//get yboost bin bounds
		pair <float, float> ys_get_bb() {return make_pair(ys_low,  ys_up);};	//get yboost bin bounds
		string get_bname() {return bname;};										//to read bin name

		//void ybys_directories(TDirectory *standard_folder); //creates folder-hierarchy
		TDirectory* ybys_directories(TDirectory *standard_folder); //creates folder-hierarchy, returns pointer to created ybys directory

};

//class for handling the resulting ratio histograms
//taking info on the ybys bin for complementing the histo with it.
class Ratio {
	//friend class YbYs;
	private:
		string histname;			//name of the histogram, will be using the ybys bname
		const char* nomhist_name;	//name of the nominator histogram (just for info purposes)
		TH1D* ratiohisto_ptr;		//pointer to summary histogram

	public:
		Ratio (YbYs, TDirectory*, TH1D*);			//constructor (will also set directory of ratiohisto) --> takes nominator (1st histo)
		void make_ratio(TH1D* denhist);				//function to compute ratio of two histograms --> takes denominator (2nd histo)
		void set_style(YbYs curbin, Color_t lc, Width_t lw, Color_t mc, Style_t ms, string ytitle);		//for setting histo style parameters (linecolor, width, etc.)
		void set_ratiotitle(string title);			//for setting the ratio title manually
		void set_histname(string histname);			//for changing the name of the ratiohistogram
		string get_histname(){return histname;};	//to read histname

};


//--------------------------------------------------------------------------------------------------//
//											 CONSTRUCTORS											//
//--------------------------------------------------------------------------------------------------//

//Constructor for Yboost Ystar bin
YbYs::YbYs (float yboost_lo, float ystar_lo, float yboost_up, float ystar_up) {
	yb_low = yboost_lo;
	ys_low = ystar_lo;
	yb_up = yboost_up;	
	ys_up = ystar_up;

	bname = Form("yb%01.0fys%01.0f", yb_low, ys_low);
	if(v==true){
		cout << Form("Created yboost-ystar-bin with: %01.2f <= yb < %01.2f and %01.2f <= ys < %01.2f", yb_low, yb_up, ys_low, ys_up) << endl;
	}
};

//Constructor for Ratio (ratio result, histogram)
Ratio::Ratio(YbYs curbin, TDirectory* ybysDir, TH1D* nomhist){
	histname = Form("ratio_%s", curbin.bname.c_str());
	ratiohisto_ptr = (TH1D*)nomhist->Clone(histname.c_str());	//clone numerator histo
	ratiohisto_ptr->SetDirectory(ybysDir);	//set ybys directory
	nomhist_name = nomhist->GetName();
};


//--------------------------------------------------------------------------------------------------//
//											 MEMBER FUNCTIONS										//
//--------------------------------------------------------------------------------------------------//

//function for handling ybys directories
//void YbYs::ybys_directories(TDirectory *standard_folder){
TDirectory* YbYs::ybys_directories(TDirectory *standard_folder){
	//create and handle directory structure of output file
	const char *dirname;  
	dirname = Form("%s", bname.c_str());
	standard_folder->mkdir(dirname);
	TDirectory *ybysDir = standard_folder->GetDirectory(dirname);
	ybysDir->cd();
	standard_folder->cd();	//change back to standard_folder

	return ybysDir;
};

//function for obtaining the ratio of nominator histogram (see constructor of Ratio) and denominator histo
void Ratio::make_ratio(TH1D* denhist){
	ratiohisto_ptr->Divide(denhist);
	//ratiohisto_ptr->Write();	//write to output file	-- not necessary, write later.

	if(v==true){
		cout << Form("Divided %s by %s.", nomhist_name, denhist->GetName()) << endl;
		cout << "Resulting histogram: " << histname.c_str() << endl;
	}
};

void Ratio::set_style(YbYs curbin, Color_t lc, Width_t lw, Color_t mc, Style_t ms, string ytitle){
	//general style parameters
	const char* title = curbin.bname.c_str();
	const char* ybname = Form("%01.0f #leg y_{b} < %01.0f", curbin.yb_low, curbin.yb_up);
	const char* ysname = Form("%f #leg y* < %f", curbin.ys_low, curbin.ys_up);

	ratiohisto_ptr->SetLineWidth(lw);
	ratiohisto_ptr->SetLineColor(lc);
	ratiohisto_ptr->SetMarkerStyle(ms);
	ratiohisto_ptr->SetMarkerColor(mc);
	//ratiohisto_ptr->SetMarkerSize();
	ratiohisto_ptr->SetTitle(title);
	ratiohisto_ptr->GetYaxis()->SetTitle(ytitle.c_str());

	//to be adjusted:
	ratiohisto_ptr->GetYaxis()->SetTitleSize(0.04);
	ratiohisto_ptr->GetYaxis()->SetTitleOffset(1.1);

	//write ybys bin info as TLatex text
	//should there be a canvas?
	/*
	TLatex* bininfo = new TLatex();
	bininfo->SetNDC();
	bininfo->SetTextSize(0.045);
	bininfo->DrawLatex(0.20, 0.86, ybname);
	bininfo->DrawLatex(0.20, 0.80, ysname);
	*/
	//--> handle such things in extra plotting script in python!

};

//function to set the ratio-title manually 
//(instead of just having the ybys binname there)
void Ratio::set_ratiotitle(string title){
	ratiohisto_ptr->SetTitle(title.c_str());
};

void Ratio::set_histname(string histname){
	ratiohisto_ptr->SetName(histname.c_str());	
};


//--------------------------------------------------------------------------------------------------//
//											 OTHER FUNCTIONS										//
//--------------------------------------------------------------------------------------------------//


//--------------------------------------------------------------------------------------------------//
//											 ANALYSIS/MAIN CODE										//
//--------------------------------------------------------------------------------------------------//

void RatioAnalysis::Analyse(){
	TH1::SetDefaultSumw2(kTRUE);		//in order to ALWAYS store the bin errors

	setTDRStyle();
	extraText="private work";


	// Create outputfile and choose inputfiles considering Run period:
	// Create file for saving the ratio histograms
	TFile* ratiofile;

	// Open prepared output files from dijet and zjet analysis
	TFile* dijetfile;
	TFile* zjetfile;
	const char* mc_inputhist;				//name of the histogram to be used if MC (smearedJet for H7) 


	//check which runperiod
	//prepare output file (= after analysis BEFORE creating nice plots)
	//the input files here contain histograms that are not yet normalised by binwidth.
	//this is no problem, because that factor would drop out in the ratio anyway.
	if(strcmp(runperiod,"A")==0){
		cout << "Run A" << endl;
		//output file
		ratiofile = new TFile("ratio_output_data_RunA.root", "RECREATE");	//for saving the ratios
		
		//input files
		///dijetfile = TFile::Open("input_files_data/dijet_interim_output_data_RunA.root", "READ");	//Dijet
		///dijetfile = TFile::Open("input_files_data/trg_lumi_output_data_RunA.root", "READ");	//Dijet
		dijetfile = TFile::Open("input_files_data/dijet_fastloop_output_data_RunA.root", "READ");	//Dijet
		zjetfile = TFile::Open("input_files_data/zjet_interim_output_data_RunA.root", "READ");		//Z+jet
	}//endif runA

	else if(strcmp(runperiod,"B")==0){
		cout << "Run B" << endl;
		//output file
		ratiofile = new TFile("ratio_output_data_RunB.root", "RECREATE");	//for saving the ratios
		
		//input files
		//dijetfile = TFile::Open("input_files_data/dijet_interim_output_data_RunB.root", "READ");	//Dijet
		//dijetfile = TFile::Open("input_files_data/trg_lumi_output_data_RunB.root", "READ");	//Dijet
		dijetfile = TFile::Open("input_files_data/dijet_fastloop_output_data_RunB.root", "READ");	//Dijet
		zjetfile = TFile::Open("input_files_data/zjet_interim_output_data_RunB.root", "READ");		//Z+jet
	}//endif runB

	else if(strcmp(runperiod,"C")==0){
		cout << "Run C" << endl;
		//output file
		ratiofile = new TFile("ratio_output_data_RunC.root", "RECREATE");	//for saving the ratios
		
		//input files
		//dijetfile = TFile::Open("input_files_data/dijet_interim_output_data_RunC.root", "READ");	//Dijet
		//dijetfile = TFile::Open("input_files_data/trg_lumi_output_data_RunC.root", "READ");	//Dijet
		dijetfile = TFile::Open("input_files_data/dijet_fastloop_output_data_RunC.root", "READ");	//Dijet
		zjetfile = TFile::Open("input_files_data/zjet_interim_output_data_RunC.root", "READ");		//Z+jet
	}//endif runC

	else if(strcmp(runperiod,"D")==0){
		cout << "Run D" << endl;
		//output file
		ratiofile = new TFile("ratio_output_data_RunD.root", "RECREATE");	//for saving the ratios
		
		//input files
		//dijetfile = TFile::Open("input_files_data/dijet_interim_output_data_RunD.root", "READ");	//Dijet
		//dijetfile = TFile::Open("input_files_data/trg_lumi_output_data_RunD.root", "READ");	//Dijet
		dijetfile = TFile::Open("input_files_data/dijet_fastloop_output_data_RunD.root", "READ");	//Dijet
		zjetfile = TFile::Open("input_files_data/zjet_interim_output_data_RunD.root", "READ");		//Z+jet
	}//endif runD

	else if(strcmp(runperiod,"mc")==0){	//IF IT IS MONTECARLO
		cout << "This is a MC ratio." << endl;
		if(strcmp(generator,"herwig")==0){
			cout << "Herwig 7" << endl;
			//output file
			ratiofile = new TFile("ratio_output_MC_herwig7.root", "RECREATE");	//for saving the ratios
			
			//input files
			dijetfile = TFile::Open("input_files_mc/dijet_normalised_output_mc_herwig7.root", "READ");		//Dijet
			zjetfile = TFile::Open("input_files_mc/zjet_normalised_output_mc_herwig7.root", "READ");		//Z+jet
			mc_inputhist = "mc_histo_smearedJet";					//if herwig --> take smeared Jet
		}
		else if(strcmp(generator,"pythia")==0){
			cout << "Pythia 8" << endl;
			//output file
			ratiofile = new TFile("ratio_output_MC_pythia8.root", "RECREATE");	//for saving the ratios
			
			//input files
			dijetfile = TFile::Open("input_files_mc/dijet_normalised_output_mc_pythia8.root", "READ");		//Dijet
			zjetfile = TFile::Open("input_files_mc/zjet_normalised_output_mc_pythia8.root", "READ");		//Z+jet
			mc_inputhist = "mc_histo";

		}
		else if(strcmp(generator, "HTpythia")==0){
			cout << "Pythia 8 (now madgraph also for dijet!)" << endl;
			//output file
			ratiofile = new TFile("ratio_output_MC_HT_madgraph_pythia8.root", "RECREATE");	//for saving the ratios

			//input files (zjet same as for "pure pythia", as here madgraph included already)
			dijetfile = TFile::Open("input_files_mc/dijet_interim_output_mc_HT_madgraph_pythia8.root", "READ"); //dijet
			zjetfile = TFile::Open("input_files_mc/zjet_normalised_output_mc_pythia8.root", "READ");			//Z+jet
			mc_inputhist = "mc_histo";
		}
		else{
			//output file
			ratiofile = new TFile("ratio_output_MC.root", "RECREATE");	//for saving the ratios
			
			//input files
			dijetfile = TFile::Open("input_files_mc/dijet_interim_output_mc.root", "READ");		//Dijet
			zjetfile = TFile::Open("input_files_mc/zjet_interim_output_mc.root", "READ");		//Z+jet
		}
	}//endif MC

	else{//if no runperiod specified, not even mc
		ratiofile = new TFile("ratio_output.root", "RECREATE");			//ratio output
		dijetfile = TFile::Open("dijet_interim_output.root", "READ");	//Dijet
		zjetfile = TFile::Open("zjet_interim_output.root", "READ");		//Z+jet
	}


	if(v==true){
		if (dijetfile->IsOpen()){
			printf("Opened dijetfile.\n");
		}
		if (zjetfile->IsOpen()){
			printf("Opened zjetfile.\n");
		}
	}

	//exit(0); //for testing

	//get standard directory of each of the input files --> not necessary
	//TDirectory *stddir = (TDirectory*)dijetfile->GetObject("Standard");
	///TDirectory *dijet_stddir;
	///dijetfile->GetObject("Standard", dijet_stddir);


	//-------------- building folder structure of output file --------------//
	// output file where folder structure will be placed in
	ratiofile->cd(); //change to output file
	TDirectory *curdir = gDirectory; //get current directory
	ratiofile->mkdir("Standard"); //create directory in output file

	//create subfolder "Standard" in outputfile
	TDirectory *standard_folder = ratiofile->GetDirectory("Standard"); 
	gROOT->GetListOfBrowsables()->Add(standard_folder, "Standard");

	standard_folder->cd();

	//make ybys bins
	YbYs *ybys_pointer;
	//call constructor YbYs(float yb_lo, float ys_lo, float yb_up, float ys_up)
	//ybys_pointer = new YbYs[6]{{0.0,0.0,1.0,1.0}, {0.0, 1.0, 1.0, 2.0}, {0.0, 2.0, 1.0, 3.0}, {1.0, 0.0, 2.0, 1.0}, {1.0, 1.0, 2.0, 2.0}, {2.0, 0.0, 3.0, 1.0}};
	//adjusted outermost bin bounds --> now go only until yb=2.4 or ys=2.4
	ybys_pointer = new YbYs[6]{{0.0,0.0,1.0,1.0}, {0.0, 1.0, 1.0, 2.0}, {0.0, 2.0, 1.0, 2.4}, {1.0, 0.0, 2.0, 1.0}, {1.0, 1.0, 2.0, 2.0}, {2.0, 0.0, 2.4, 1.0}};

	//vectors containing the style for each ybys bin
	vector<Color_t> cols = {kYellow+1, kGray+2, kGreen+2, kBlue-2, kOrange-3, kRed+2, kAzure+1, kViolet+6, kSpring+5, kViolet-5, kCyan-3};
	vector<Style_t> markers = {kOpenCircle, kOpenTriangleUp, kOpenSquare, kOpenDiamond, kOpenCross, kOpenTriangleDown};

	//go through ybys bins, make ratios
	for (int ybysnr=0; ybysnr!=6; ++ybysnr){
		YbYs cur_ybysbin = ybys_pointer[ybysnr];	//cur_ybysbin is not a pointer, but the ybysbin itself!!
		const char* ybys_bname = cur_ybysbin.get_bname().c_str();

		//ybys directories -- folder structure	(actually so far only creates the ybys bin folder)
		TDirectory* ybys_Dir = cur_ybysbin.ybys_directories(standard_folder);

		cout << "ybys name: " << ybys_bname << endl;
		//RATIOS!
		//get summary histograms from dijet and zjet analyses
		//using "GetObject(<objname>, <newobj>);" or "Get(<objname>)->Clone(<newobjname>);" ?
		//should organise that better in the dijet-output (put result in different folder)
	
		////dijetfile->ls();
		////cout << "Dijetfile obj: " << dijetfile->Get(Form("Standard/%s/trigger_ratios/sum_histo_%s", ybys_bname, ybys_bname)) << endl;
		////TH1D* dijethisto = (TH1D*)dijetfile->Get(Form("Standard/%s/trigger_ratios/sum_histo_%s", ybys_bname, ybys_bname))->Clone(Form("dijet_%s", ybys_bname));
		////TH1D* zjethisto = (TH1D*)zjetfile->Get(Form("Standard/%s/additional_plots_%s/sum_histo_%s", ybys_bname, ybys_bname, ybys_bname))->Clone(Form("zjet_%s", ybys_bname));

		TH1D* dijethisto;
		TH1D* zjethisto;
		if(strcmp(runperiod,"mc")==0){	//IF IT IS MONTECARLO
			//dijetfile->GetObject(Form("Standard/%s/mc_histo_%s", ybys_bname, ybys_bname), dijethisto);
			dijetfile->GetObject(Form("Standard/%s/%s_%s", ybys_bname, mc_inputhist, ybys_bname), dijethisto); //for H7: take smearedJet histo
			//zjetfile->GetObject(Form("Standard/%s/mc_histo_%s", ybys_bname, ybys_bname), zjethisto);
			zjetfile->GetObject(Form("Standard/%s/%s_%s", ybys_bname, mc_inputhist, ybys_bname), zjethisto); //for H7: take smearedJet histo
		}
		else{ //if it is data, no matter which runperiod
			//dijetfile->GetObject(Form("Standard/%s/trigger_ratios/sumhisto_%s", ybys_bname, ybys_bname), dijethisto);
			dijetfile->GetObject(Form("Standard/%s/sumhisto_%s", ybys_bname, ybys_bname), dijethisto);
			zjetfile->GetObject(Form("Standard/%s/additional_plots_%s/sumhisto_%s", ybys_bname, ybys_bname, ybys_bname), zjethisto);
		}
	
		dijethisto->SetName(Form("dijet_%s", ybys_bname));
		zjethisto->SetName(Form("zjet_%s", ybys_bname));

		//the following should be implemented in a safer way, i.e. happen long before this script...
		if((strcmp(generator, "HTpythia")==0)&&(strcmp(runperiod,"mc")==0)){
			dijethisto->Scale(1000);	//as it is still in pb if its not taken from the binwidth-noramlised file
		}
			



		/*
		// TESTING without clone..
		dijetfile->cd(Form("Standard/%s/trigger_ratios", ybys_bname));
		//gDirectory->ls();
		TH1D* dijethisto = (TH1D*)gDirectory->Get(Form("sumhisto_%s", ybys_bname));

		zjetfile->cd(Form("Standard/%s/additional_plots_%s", ybys_bname, ybys_bname));
		//gDirectory->ls();
		TH1D* zjethisto = (TH1D*)gDirectory->Get(Form("sumhisto_%s", ybys_bname));
		*/

		ybys_Dir->cd();
		dijethisto->Write();
		zjethisto->Write();
		


		//Doing this before Divide() [in case Errors shall be used later?]
		dijethisto->Sumw2();
		zjethisto->Sumw2();

		//TH1D* dijethisto = (TH1D*)dijetfile->Get(Form("Standard/%s/trigger_ratios/sum_histo_%s", ybys_bname, ybys_bname));
		//TH1D* zjethisto = (TH1D*)zjetfile->Get(Form("Standard/%s/additional_plots_%s/sum_histo_%s", ybys_bname, ybys_bname, ybys_bname));


		//cout << "zjethisto " << zjethisto << endl;


		
		//create Ratio object for current ybys bin
		Ratio cur_ratio(cur_ybysbin, ybys_Dir, dijethisto);	//in this case the dijet histogram is the numerator
		Ratio cur_ratio_inverted(cur_ybysbin, ybys_Dir, zjethisto); //here z+jet as numerator
		cur_ratio_inverted.set_histname(Form("ratio_ZjetOverDijet_%s", ybys_bname));


		//make ratio
		string ytitle = "Dijet / Z+jet";			//name of the y-axis
		string ytitle_inv = "Z+jet / Dijet";		//name of the y-axis of inverted histo
		cur_ratio.make_ratio(zjethisto);			//divide dijethisto by zjethisto
		cur_ratio_inverted.make_ratio(dijethisto);	//divide zjethisto by dijethisto

		//setting the histogram style parameters (very basic so far -> is no canvas)
		cur_ratio.set_style(cur_ybysbin, cols.at(ybysnr+2), 2, cols.at(ybysnr+2), markers.at(ybysnr), ytitle);
		cur_ratio_inverted.set_style(cur_ybysbin, cols.at(ybysnr+2), 2, cols.at(ybysnr+2), markers.at(ybysnr), ytitle_inv);
	}

	ratiofile->Write();

}


//function using two previous outputs of RatioAnalysis() function (one for data, one for mc)
//creating a ratio of ratios
//runperiod and generator variables need to be set
void RatioAnalysis::RatioOfRatios(){
	TH1::SetDefaultSumw2(kTRUE);		//in order to ALWAYS store the bin errors

	setTDRStyle();

	// Create outputfile and choose inputfiles considering Run period and generator
	// Create file for saving the ratio of ratios histograms
	TFile* ror_outputfile;	//ratio of ratios output file
	ror_outputfile = new TFile(Form("ratio_output_RatioDataMC_Run%s_%s.root", runperiod, generator), "RECREATE");

	// Open prepared ratio files from ratio analysis
	TFile* mc_inputfile;
	TFile* data_inputfile;

	//check which runperiod
	//prepare output file (= after analysis BEFORE creating nice plots)
	if(strcmp(runperiod,"A")==0){
		cout << "Run A" << endl;
		//data input ratio file
		data_inputfile = TFile::Open("ratio_output_data_RunA.root", "READ");	//data ratios dijet/zjet
		data_inputfile->cd("Standard");
		gDirectory->ls();
	}//endif runA
	else if(strcmp(runperiod,"B")==0){
		cout << "Run B" << endl;
		//data input ratio file
		data_inputfile = TFile::Open("ratio_output_data_RunB.root", "READ");		//data ratios dijet/zjet
	}//endif runB
	else if(strcmp(runperiod,"C")==0){
		cout << "Run C" << endl;
		//data input ratio file
		data_inputfile = TFile::Open("ratio_output_data_RunC.root", "READ");		//data ratios dijet/zjet
	}//endif runC
	else if(strcmp(runperiod,"D")==0){
		cout << "Run D" << endl;
		//data input ratio file
		data_inputfile = TFile::Open("ratio_output_data_RunD.root", "READ");		//data ratios dijet/zjet
	}//endif runD

	else if(strcmp(runperiod,"mc")==0){	//IF IT IS MONTECARLO
		cout << "runperiod = " << runperiod << endl;
		cout << "But this function needs both: a valid runperiod and valid generator choice." << endl;
		cout << "Aborting." << endl;
		exit(EXIT_FAILURE);
	}

	//check with generator
	if(strcmp(generator,"herwig")==0){
		cout << "Herwig 7" << endl;
		//mc input ratio file
		mc_inputfile = TFile::Open("ratio_output_MC_herwig7.root", "READ");		//mc ratios dijet/zjet
	}
	else if(strcmp(generator,"pythia")==0){
		cout << "Pythia 8" << endl;
		//mc input ratio file
		mc_inputfile = TFile::Open("ratio_output_MC_pythia8.root", "READ");			//mc ratios dijet/zjet
	}
	else if(strcmp(generator, "HTpythia")==0){
		cout << "Pythia 8 (now madgraph also for dijet!)" << endl;
		//mc input ratio file
		mc_inputfile = TFile::Open("ratio_output_MC_HT_madgraph_pythia8.root", "READ");	//mc ratios
	}

	else{
		cout << "generator = " << generator << endl;
		cout << "But this function needs both: a valid runperiod and valid generator choice." << endl;
		cout << "Aborting." << endl;
		exit(EXIT_FAILURE);
	}

	
	if(v==true){
		if (data_inputfile->IsOpen()){
			printf("Opened data input file.\n");
		}
		if (mc_inputfile->IsOpen()){
			printf("Opened mc input file.\n");
		}
	}


	//-------------- building folder structure of output file --------------//
	// output file where folder structure will be placed in
	ror_outputfile->cd(); 				//change to output file
	TDirectory *curdir = gDirectory; 	//get current directory
	ror_outputfile->mkdir("Standard"); 	//create directory in output file

	//create subfolder "Standard" in outputfile
	TDirectory *standard_folder = ror_outputfile->GetDirectory("Standard"); 
	gROOT->GetListOfBrowsables()->Add(standard_folder, "Standard");
	standard_folder->cd();

	//make ybys bins
	YbYs *ybys_pointer;
	//call constructor YbYs(float yb_lo, float ys_lo, float yb_up, float ys_up)
	//ybys_pointer = new YbYs[6]{{0.0,0.0,1.0,1.0}, {0.0, 1.0, 1.0, 2.0}, {0.0, 2.0, 1.0, 3.0}, {1.0, 0.0, 2.0, 1.0}, {1.0, 1.0, 2.0, 2.0}, {2.0, 0.0, 3.0, 1.0}};
	//adjusted outermost bin bounds --> now go only until yb=2.4 or ys=2.4
	ybys_pointer = new YbYs[6]{{0.0,0.0,1.0,1.0}, {0.0, 1.0, 1.0, 2.0}, {0.0, 2.0, 1.0, 2.4}, {1.0, 0.0, 2.0, 1.0}, {1.0, 1.0, 2.0, 2.0}, {2.0, 0.0, 2.4, 1.0}};

	//vectors containing the style for each ybys bin
	vector<Color_t> cols = {kYellow+1, kGray+2, kGreen+2, kBlue-2, kOrange-3, kRed+2, kAzure+1, kViolet+6, kSpring+5, kViolet-5, kCyan-3};
	vector<Style_t> markers = {kOpenCircle, kOpenTriangleUp, kOpenSquare, kOpenDiamond, kOpenCross, kOpenTriangleDown};

	//go through ybys bins, make ratios
	for (int ybysnr=0; ybysnr!=6; ++ybysnr){
		YbYs cur_ybysbin = ybys_pointer[ybysnr];	//cur_ybysbin is not a pointer, but the ybysbin itself!!
		const char* ybys_bname = cur_ybysbin.get_bname().c_str();

		//ybys directories -- folder structure	(actually so far only creates the ybys bin folder)
		TDirectory* ybys_Dir = cur_ybysbin.ybys_directories(standard_folder);

		cout << "ybys name: " << ybys_bname << endl;
		//RATIOS!
		//get summary histograms from dijet and zjet analyses
		//using "GetObject(<objname>, <newobj>);" or "Get(<objname>)->Clone(<newobjname>);" ?
		//should organise that better in the dijet-output (put result in different folder)
	
		////dijetfile->ls();
		////cout << "Dijetfile obj: " << dijetfile->Get(Form("Standard/%s/trigger_ratios/sum_histo_%s", ybys_bname, ybys_bname)) << endl;
		////TH1D* dijethisto = (TH1D*)dijetfile->Get(Form("Standard/%s/trigger_ratios/sum_histo_%s", ybys_bname, ybys_bname))->Clone(Form("dijet_%s", ybys_bname));
		////TH1D* zjethisto = (TH1D*)zjetfile->Get(Form("Standard/%s/additional_plots_%s/sum_histo_%s", ybys_bname, ybys_bname, ybys_bname))->Clone(Form("zjet_%s", ybys_bname));

		//these are just copied
		TH1D* data_dijethisto;
		TH1D* data_zjethisto;
		TH1D* mc_dijethisto;
		TH1D* mc_zjethisto;

		//these will be used	(using Ratio class!)
		TH1D* data_ratio;
		TH1D* mc_ratio;

		//get all the histograms
		data_inputfile->GetObject(Form("Standard/%s/dijet_%s", ybys_bname, ybys_bname), data_dijethisto);
		data_inputfile->GetObject(Form("Standard/%s/zjet_%s", ybys_bname, ybys_bname), data_zjethisto);
		data_inputfile->GetObject(Form("Standard/%s/ratio_%s", ybys_bname, ybys_bname), data_ratio);
		mc_inputfile->GetObject(Form("Standard/%s/dijet_%s", ybys_bname, ybys_bname), mc_dijethisto);
		mc_inputfile->GetObject(Form("Standard/%s/zjet_%s", ybys_bname, ybys_bname), mc_zjethisto);
		mc_inputfile->GetObject(Form("Standard/%s/ratio_%s", ybys_bname, ybys_bname), mc_ratio);

		//setting all the names (apart from the final ratio of ratios --> will be set in Ratio class)
		data_dijethisto->SetName(Form("data_dijet_%s", ybys_bname));
		data_zjethisto->SetName(Form("data_zjet_%s", ybys_bname));
		mc_dijethisto->SetName(Form("mc_dijet_%s", ybys_bname));
		mc_zjethisto->SetName(Form("mc_zjet_%s", ybys_bname));
		data_ratio->SetName(Form("data_ratio_%s", ybys_bname));
		mc_ratio->SetName(Form("mc_ratio_%s", ybys_bname));



		//write the ones that are just copied: (in data / mc folder?)
		ybys_Dir->cd();
		data_dijethisto->Write();
		data_zjethisto->Write();
		data_ratio->Write();
		mc_dijethisto->Write();
		mc_zjethisto->Write();
		mc_ratio->Write();


		//Doing this before Divide() [in case Errors shall be used later?]
		data_ratio->Sumw2();
		mc_ratio->Sumw2();


		//create Ratio object for current ybys bin
		Ratio cur_ratio(cur_ybysbin, ybys_Dir, data_ratio);	//in this case the data histogram is the numerator

		//make ratio
		string ytitle = "Data/MC (both first Dijet/Z+jet)";	//name of the y-axis
		cur_ratio.make_ratio(mc_ratio);	//divide data ratio histo by mc ratio histo

		//setting the histogram style parameters (very basic so far -> is no canvas)
		cur_ratio.set_style(cur_ybysbin, cols.at(ybysnr+2), 2, cols.at(ybysnr+2), markers.at(ybysnr), ytitle);
	}

	ror_outputfile->Write();
}// end RatioOfRatios() function


//function to calculate the ratio of the uncertainty ratios (zjet_uncRatio)/(dijet_uncRatio)
//only takes MC Herwig 7 (for zjet and dijet)
void RatioAnalysis::RatioOfUncs(){
	TH1::SetDefaultSumw2(kTRUE);		//in order to ALWAYS store the bin errors

	setTDRStyle();

	// Create outputfile and choose inputfiles --> MC herwig 7 (zjet, dijet)
	// Create file for saving the ratio of uncratios histograms
	TFile* unc_doubleratios_outputfile;	//ratio of ratios output file
	unc_doubleratios_outputfile = new TFile(Form("ratio_output_UncDoubleRatios_mc_herwig7.root"), "RECREATE");

	// Open prepared uncertainty-ratio (DivUpDown) files from each analysis
	TFile* zjet_inputfile;
	TFile* dijet_inputfile;

	//prepare output file (= after analysis BEFORE creating nice plots)
	//zjet input unc-ratio file (= zjet mc herwig 7 uncertainty ratio file)
	zjet_inputfile = TFile::Open("input_files_uncs/zjet_divided_uncsources_mc_herwig7.root", "READ");	//mc uncertainty sources ratios zjet
	zjet_inputfile->cd("UncertaintyRatios");

	//dijet input unc-ratio file (= dijet mc herwig 7 uncertainty ratio file)
	dijet_inputfile = TFile::Open("input_files_uncs/dijet_divided_uncsources_mc_herwig7.root", "READ");	//mc uncertainty sources ratios dijet
	dijet_inputfile->cd("UncertaintyRatios");

	/*
	if(strcmp(runperiod,"mc")==0){	//IF IT IS MONTECARLO
		cout << "runperiod = " << runperiod << endl;
		cout << "But this function needs both: a valid runperiod and valid generator choice." << endl;
		cout << "Aborting." << endl;
		exit(EXIT_FAILURE);
	}
	*/
	
	if(v==true){
		if (zjet_inputfile->IsOpen()){
			printf("Opened zjet input file.\n");
		}
		if (dijet_inputfile->IsOpen()){
			printf("Opened dijet input file.\n");
		}
	}


	//-------------- building folder structure of output file --------------//
	// output file where folder structure will be placed in
	unc_doubleratios_outputfile->cd(); 					//change to output file
	TDirectory *curdir = gDirectory; 					//get current directory

	//create subfolder "Standard" in outputfile
	unc_doubleratios_outputfile->mkdir("Standard"); 	//create directory in output file
	TDirectory *standard_folder = unc_doubleratios_outputfile->GetDirectory("Standard"); 
	gROOT->GetListOfBrowsables()->Add(standard_folder, "Standard");
	standard_folder->cd();

	//make ybys bins
	YbYs *ybys_pointer;
	//call constructor YbYs(float yb_lo, float ys_lo, float yb_up, float ys_up)
	//adjusted outermost bin bounds --> now go only until yb=2.4 or ys=2.4
	ybys_pointer = new YbYs[6]{{0.0,0.0,1.0,1.0}, {0.0, 1.0, 1.0, 2.0}, {0.0, 2.0, 1.0, 2.4}, {1.0, 0.0, 2.0, 1.0}, {1.0, 1.0, 2.0, 2.0}, {2.0, 0.0, 2.4, 1.0}};

	//vectors containing the style for each ybys bin
	vector<Color_t> cols = {kYellow+1, kGray+2, kGreen+2, kBlue-2, kOrange-3, kRed+2, kAzure+1, kViolet+6, kSpring+5, kViolet-5, kCyan-3};
	vector<Style_t> markers = {kOpenCircle, kOpenTriangleUp, kOpenSquare, kOpenDiamond, kOpenCross, kOpenTriangleDown};


	//Vector with the Uncertainty Sources names:
	vector<string> uncnames;
	uncnames.push_back("AbsoluteStat");
	uncnames.push_back("AbsoluteScale");
	uncnames.push_back("AbsoluteSample");
	uncnames.push_back("AbsoluteFlavMap");
	uncnames.push_back("AbsoluteMPFBias");
	uncnames.push_back("Fragmentation");
	uncnames.push_back("SinglePionECAL");
	uncnames.push_back("SinglePionHCAL");
	//uncnames.push_back("FlavorQCD");			//--> has to be treated separately in this ratio
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
	//uncnames.push_back("FlavorZJet");			//--> has to be treated separately in this ratio
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


	//go through ybys bins, make ratios
	for (int ybysnr=0; ybysnr!=6; ++ybysnr){
		YbYs cur_ybysbin = ybys_pointer[ybysnr];	//cur_ybysbin is not a pointer, but the ybysbin itself!!
		const char* ybys_bname = cur_ybysbin.get_bname().c_str();

		//ybys directories -- folder structure	(actually so far only creates the ybys bin folder)
		TDirectory* ybys_Dir = cur_ybysbin.ybys_directories(standard_folder);
		cout << "ybys name: " << ybys_bname << endl;

		//loop through all the uncertainty names --> for ZJET --> divide them by corresponding dijet uncratio
		//always the same unc sources, except for FlavorQCD (dijet) with FlavorZJet (zjet) and the other vice versa
		for(int k=0; k!=uncnames.size(); ++k){	//should have length 54, as FlavorQCD and FlavorZJet are treated separately
			//these are just copied
			TH1D* zjet_cur_uncratio_hist;
			TH1D* dijet_cur_uncratio_hist;

			//these will be used	(using Ratio class!)
			TH1D* data_ratio;
			TH1D* mc_ratio;

			//get all the histograms (they will later be cloned in constructor of "Ratio object"
			zjet_inputfile->GetObject(Form("UncertaintyRatios/%s/UpDownDiv_%s_%s", ybys_bname, uncnames.at(k).c_str(), ybys_bname), zjet_cur_uncratio_hist);
			dijet_inputfile->GetObject(Form("UncertaintyRatios/%s/UpDownDiv_%s_%s", ybys_bname, uncnames.at(k).c_str(), ybys_bname), dijet_cur_uncratio_hist);

			//setting all the names (apart from the final ratio of ratios --> will be set in Ratio class)
			zjet_cur_uncratio_hist->SetName(Form("ratio_%s_%s", uncnames.at(k).c_str(), ybys_bname));
			dijet_cur_uncratio_hist->SetName(Form("ratio_%s_%s", uncnames.at(k).c_str(), ybys_bname));


			//Doing this before Divide() [in case Errors shall be used later?]
			zjet_cur_uncratio_hist->Sumw2();
			dijet_cur_uncratio_hist->Sumw2();


			//create Ratio object for current ybys bin and current Uncertainty Source
			//given histogram is the nominator, so here always the Z+jet unc(ratio) histo
			///Ratio cur_ratio(cur_ybysbin, ybys_Dir, zjet_cur_uncratio_hist);	//in this case the zjet histogram is the nominator
			Ratio cur_ratio(cur_ybysbin, ybys_Dir, dijet_cur_uncratio_hist);	//in this case the dijet histogram is the numerator


			//make ratio
			//string ytitle = "Z+jet/Dijet (both first UncUp/UncDown)";	//name of the y-axis
			string ytitle = "Dijet/Z+jet (both first UncUp/UncDown)";	//name of the y-axis
			//string ratio_title = Form("%s", uncnames.at(k).c_str());	//title
			string ratio_title = uncnames.at(k);						//title
			string ratio_histname = Form("ratio_%s_%s", uncnames.at(k).c_str(), cur_ybysbin.get_bname().c_str());
			//cur_ratio.make_ratio(dijet_cur_uncratio_hist);				//divide zjet uncratio histo by mc dijet uncratio histo
			cur_ratio.make_ratio(zjet_cur_uncratio_hist);				//divide dijet uncratio histo by mc zjet uncratio histo


			//setting the histogram style parameters (very basic so far -> is no canvas)
			cur_ratio.set_style(cur_ybysbin, cols.at(ybysnr+2), 2, cols.at(ybysnr+2), markers.at(ybysnr), ytitle);
			cur_ratio.set_ratiotitle(ratio_title);
			cur_ratio.set_histname(ratio_histname);

		}//end of loop through unc sources

		//now add the two missing ratios: FlavorQCD(dijet)/FlavorZjet(zjet) and FloverZjet(dijet)/FlavorQCD(zjet)
		TH1D* zjet_flavorqcd;
		TH1D* zjet_flavorzjet;
		TH1D* dijet_flavorzjet;
		TH1D* dijet_flavorqcd;

		//get them
		zjet_inputfile->GetObject(Form("UncertaintyRatios/%s/UpDownDiv_FlavorQCD_%s", ybys_bname, ybys_bname), zjet_flavorqcd);
		zjet_inputfile->GetObject(Form("UncertaintyRatios/%s/UpDownDiv_FlavorZJet_%s", ybys_bname, ybys_bname), zjet_flavorzjet);
		dijet_inputfile->GetObject(Form("UncertaintyRatios/%s/UpDownDiv_FlavorZJet_%s", ybys_bname, ybys_bname), dijet_flavorzjet);
		dijet_inputfile->GetObject(Form("UncertaintyRatios/%s/UpDownDiv_FlavorQCD_%s", ybys_bname, ybys_bname), dijet_flavorqcd);

		//store error-bars
		zjet_flavorqcd->Sumw2();
		zjet_flavorzjet->Sumw2();
		dijet_flavorzjet->Sumw2();
		dijet_flavorqcd->Sumw2();

		/*
		Ratio zjetFlavQCD_dijetFlavZJet(cur_ybysbin, ybys_Dir, zjet_flavorqcd);		//create the ratio FlavorQCD (zjet) / FlavorZJet (dijet)
		Ratio zjetFlavZJet_dijetFlavQCD(cur_ybysbin, ybys_Dir, zjet_flavorzjet);	//create the ratio FlavorZJet (zjet) / FlavorQCD (dijet)
		*/

		Ratio dijetFlavZJet_zjetFlavQCD(cur_ybysbin, ybys_Dir, dijet_flavorzjet);		//create the ratio FlavorQCD (zjet) / FlavorZJet (dijet)
		Ratio dijetFlavQCD_zjetFlavZJet(cur_ybysbin, ybys_Dir, dijet_flavorqcd);	//create the ratio FlavorZJet (zjet) / FlavorQCD (dijet)
		

	
		//y-axis titles and overall titles
		//string ytitle = "Z+jet/Dijet (both first UncUp/UncDown)";	//name of the y-axis
		string ytitle = "Dijet/Z+jet (both first UncUp/UncDown)";	//name of the y-axis


		///string title_zjetQCD = "FlavorQCD (Zjet) / FlavorZJet (Dijet)";
		///string title_zjetZJet = "FlavorZJet (Zjet) / FlavorQCD (Dijet)";
		string title_dijetZJet = "FlavorZJet (Dijet) / FlavorQCD (Zjet)";
		string title_dijetQCD = "FlavorQCD (Dijet) / FlavorZJet (Zjet)";


		//make ratios
		///zjetFlavQCD_dijetFlavZJet.make_ratio(dijet_flavorzjet);
		///zjetFlavZJet_dijetFlavQCD.make_ratio(dijet_flavorqcd);
		dijetFlavZJet_zjetFlavQCD.make_ratio(zjet_flavorqcd);
		dijetFlavQCD_zjetFlavZJet.make_ratio(zjet_flavorzjet);

		//setting some simple style parameters
		dijetFlavZJet_zjetFlavQCD.set_style(cur_ybysbin, cols.at(ybysnr+2), 2, cols.at(ybysnr+2), markers.at(ybysnr), ytitle);
		dijetFlavQCD_zjetFlavZJet.set_style(cur_ybysbin, cols.at(ybysnr+2), 2, cols.at(ybysnr+2), markers.at(ybysnr), ytitle);
		dijetFlavZJet_zjetFlavQCD.set_ratiotitle(title_dijetZJet);
		dijetFlavQCD_zjetFlavZJet.set_ratiotitle(title_dijetQCD);

		/*
		zjetFlavQCD_dijetFlavZJet.set_style(cur_ybysbin, cols.at(ybysnr+2), 2, cols.at(ybysnr+2), markers.at(ybysnr), ytitle);
		zjetFlavZJet_dijetFlavQCD.set_style(cur_ybysbin, cols.at(ybysnr+2), 2, cols.at(ybysnr+2), markers.at(ybysnr), ytitle);
		zjetFlavQCD_dijetFlavZJet.set_ratiotitle(title_zjetQCD);
		zjetFlavZJet_dijetFlavQCD.set_ratiotitle(title_zjetZJet);
		*/

		//set histogram names (for writing)
		dijetFlavZJet_zjetFlavQCD.set_histname(Form("dijetFlavZJet_zjetFlavQCD_%s", cur_ybysbin.get_bname().c_str()));
		dijetFlavQCD_zjetFlavZJet.set_histname(Form("dijetFlavQCD_zjetFlavZJet_%s", cur_ybysbin.get_bname().c_str()));

		/*
		zjetFlavQCD_dijetFlavZJet.set_histname(Form("zjetFlavQCD_dijetFlavZJet_%s", cur_ybysbin.get_bname().c_str()));
		zjetFlavZJet_dijetFlavQCD.set_histname(Form("zjetFlavZJet_dijetFlavQCD_%s", cur_ybysbin.get_bname().c_str()));
		*/


	}//end of ybys loop

	unc_doubleratios_outputfile->Write();

	cout << endl;
	cout << "-----------------------------------------------------------------------" << endl;
	cout << Form("Results can be found in: %s", unc_doubleratios_outputfile->GetName()) << endl;
}// end RatioOfUncs() function




//---------------------------------------------------------------------------------------------------------//
//other procedure for opening files... why?
	/*
	TFile* dijetfile = (TFile*)gROOT->GetListOfFiles()->FindObject("dijet_interim_output.root");	//Dijet
	TFile* zjetfile = (TFile*)gROOT->GetListOfFiles()->FindObject("zjet_interim_output.root");	//Zjet

	if (!dijetfile || !dijetfile->IsOpen()){
		cout << "first if condition " << endl;
		dijetfile = TFile::Open("dijet_interim_output.root", "READ");
	
		cout << "after TFile::Open()" << endl;
		cout << "get name: " << dijetfile->GetName() << endl;
		dijetfile->IsOpen();
	}
	if (dijetfile->IsOpen()){
		printf("Opened dijetfile.\n");
	}
	*/


