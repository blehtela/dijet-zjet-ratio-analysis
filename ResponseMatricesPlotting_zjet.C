//#include "ReprodOutput.h"
#include <string>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "tdrstyle_mod14.C" //use this for canvases etc.
#include <TMath.h>

using namespace std;

// Script for plotting of histograms obtained from a .root file
//--------------------------------------------------------------

//new:
void ResponseMatricesPlotting_zjet(){
	setTDRStyle();
	extraText="private work";

	// Open output file
	TFile* inputfile = (TFile*)gROOT->GetListOfFiles()->FindObject("zjet_matrices.root"); //could be handled via input? (so that plotting can be used for any root-file)
	cout << "after inputfile " << endl;
	cout << "&inputfile " << &inputfile << endl;
	if (!inputfile || !inputfile->IsOpen()){
		//inputfile = new TFile("interim_output.root", "READ"); //does not make sense to create this file, want to read from it!
		cout << "first if condition " << endl;
		inputfile = TFile::Open("zjet_matrices.root", "READ");
		cout << "after TFile::Open()" << endl;
		cout << "get name: " << inputfile->GetName() << endl;
		inputfile->IsOpen();
	}
	if (inputfile->IsOpen()){
		printf("Opened file.\n");
	}

	//exit(0); //for testing


	//change to directory in input file, containing hist0 --> draw hist0 to canvas and save it in outputfile (same directory and histo name)
	//list of etabins:
	//vector<const char*> etabin_names;

	/*
	//get standard directory in input file
	//TDirectory *stddir = (TDirectory*)inputfile->GetObject("Standard");
	TDirectory *stddir;
	inputfile->GetObject("Standard", stddir);
	
	plotfile->cd(); 		//need to be in output file to start following function
	//const char* plotname = "plot";		//initial plotname --> will be updated within the following function
	string plotname = "plot";
	//int histcount = 0;						//to count histograms that are redrawn 
	copydir(plotfile, stddir, plotname, histcount);
	

	//create canvas for drawing the histograms (contained in file f)
	//TCanvas *canv0 = tdrCanvas();
	
	*/
	
	//create trigger-summary plot for each eta-bin (should be generalised, this is just test version)
	//vector<const char*> etabins = {"Eta_00_05", "Eta_05_10", "Eta_10_15", "Eta_15_20", "Eta_20_25", "Eta_25_30", "Eta_30_32",  "Eta_32_47"};
	vector<const char*> etabins = {"yb0ys0", "yb0ys1", "yb0ys2", "yb1ys0", "yb1ys1", "yb2ys0" };
	vector<Double_t> upperx = {1000, 700, 700, 600, 600, 600};

	//plot one TH2D for each bin
	for(int k=0; k!=etabins.size(); ++k){
		TDirectory* etadir; 
		//plotfile->GetObject(Form("Standard/%s", etabins.at(k)), etadir);
		inputfile->GetObject(Form("Standard/%s", etabins.at(k)), etadir);

		string plotname = Form("plot_zjet_%s", etabins.at(k));

		const char* etaname = etadir->GetName();
		string expr = "hpt";	//expression that has to be contained in histogram-name in order for histo to appear in summary plot
		vector<TH1D*> histosvec;
		//gethistos(etadir, expr, histosvec);	//fills histograms in vector
		
		TH2D* curhist;
		inputfile->GetObject(Form("Standard/%s/normalised_response", etabins.at(k)), curhist);


		//TCanvas *curcanv = tdrCanvas("curcanv", curhist, 2, 11, kSquare); 
		//TCanvas *curcanv = tdrCanvas("curcanv", curhist, 2, 0, kSquare); 
		TCanvas *curcanv = new TCanvas("curcanv", "curcanv", 600, 600); 

		curcanv->cd();
		curcanv->SetLogx();
		curcanv->SetLogy();
		curcanv->SetGridx();
		curcanv->SetGridy();

		//curcanv->SetLeftMargin(0.2);
		//curcanv->SetRightMargin(0.2);
		curcanv->SetMargin(0.16, 0.16, 0.16, 0.16); //left, right, bottom, top
	
		//gStyle->SetOptStat("nemr");	//draw statistics box with histogram name and number of entries
		//gStyle->SetStatStyle(0);	//transparent
		//gStyle->SetStatStyle(111);	//not transparent
		//curhist->SetStats(kTRUE);
		curhist->SetStats(kFALSE);
		curhist->Draw("COLZ");
	
		curhist->GetXaxis()->SetRangeUser(50,upperx.at(k));
		curhist->GetYaxis()->SetRangeUser(50,upperx.at(k));

		//curhist->Draw();
		//tdrDraw(curhist, "L, SAME, ][");
		//tdrDraw(curhist, "HIST, ][");
		//curhist->SetLineColor(linecol);
		//curhist->SetLineWidth(3);
		//curhist->SetMarkerSize(0.1);
		//curhist->GetXaxis()->SetMoreLogLabels();
		curhist->GetXaxis()->SetNoExponent();
		curhist->GetYaxis()->SetNoExponent();
		curhist->GetXaxis()->SetMoreLogLabels();
		curhist->GetYaxis()->SetMoreLogLabels();

		//curhist->GetXaxis()->SetTitle("p_{T} [GeV]");
		curhist->GetXaxis()->SetTitleSize(0.03);
		curhist->GetXaxis()->SetLabelSize(0.03);
		curhist->GetXaxis()->SetTitleOffset(1.2);
		//curhist->GetYaxis()->SetTitle("Events");
		curhist->GetYaxis()->SetTitleSize(0.03);
		curhist->GetYaxis()->SetLabelSize(0.03);
		curhist->GetYaxis()->SetTitleOffset(1.8);
		//curcanv->SetLeftMargin(0.22);
		curhist->GetZaxis()->SetTitleSize(0.03);
		curhist->GetZaxis()->SetLabelSize(0.03);
		curhist->GetZaxis()->SetTitleOffset(1.2);
	
		//CMS_logo->Draw("X")		//see tdrStyle file
		//curhist->tdrDraw(self, "");
		///curhist->SetFillColor(kRed);
		///curhist->SetFillStyle(3001);
		//curhist->SetFillStyle(0);
		//curhist->Write("", TObject::kOverwrite);	//overwrite existing histogram

		/*
		gPad->Update();
		TPaveStats *stat = (TPaveStats*)curhist->FindObject("stats");	//get statistics box
		stat->SetX1NDC(0.8);		//set position of statistics box
		stat->SetX2NDC(0.94);
		stat->SetY1NDC(0.8);
		stat->SetY2NDC(0.92);
		//curhist->SetStats(kFALSE);		//remove original statbox
		*/
		
		curcanv->Update();
		//curcanv->Print(Form("%s.pdf", plotname));
		curcanv->SaveAs(Form("%s.pdf", plotname.c_str()));
		//delete curcanv;
		curcanv->Close();
	
		//plotfile->cd();
		//trigger_summary(histosvec, etaname);
	}



	inputfile->Close();		//close file
	//plotfile->Close();
}


// FROM PREVIOUS VERSION OF REPRODOUTPUT.C
//void Etabins::eta_bin_plotting(TDirectory *standard_folder) { //due to eta_directories --> need to change to correct directory here!!!
//new:
/*
void Etabins::eta_bin_plotting(TFile interim_output, TDirectory *standard_folder){
	TFile output_file = TFile::Open("interim_output.root");
};
*/
	/*
	//create and handle directory structure of output file
	const char *dirname;  
	dirname = Form("Eta_%s", bname.c_str());
	standard_folder->mkdir(dirname);
	TDirectory *Dirname = standard_folder->GetDirectory(dirname);
	Dirname->cd();
	*/

	/* ///
	eta_canv->cd();
	eta_histo->SetLineColor(38);
	gStyle->SetOptStat("nemr");	//draw statistics box with histogram name and number of entries
	eta_canv->SetLogx();
	eta_canv->SetLogy();
	eta_histo->Draw();	//draw histogram
	eta_histo->Write("", TObject::kOverwrite); //with options to avoid duplicate histograms
	eta_canv->Update();
	eta_canv->Print(Form("%s", plotname.c_str()));
	//eta_canv->Print(Form("%s.pdf", canvname)); //somehow causes problems... undefined behaviour in naming (char problem?)
	*/ ///
//};


