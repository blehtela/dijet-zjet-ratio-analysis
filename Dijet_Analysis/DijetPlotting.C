#include "DijetAnalysis.h"
#include <string>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "tdrstyle_mod14.C" //use this for canvases etc.
#include <TMath.h>

using namespace std;

// Script for plotting of histograms obtained from a .root file
//--------------------------------------------------------------

//adjusted to dijet variables pTavg, yboost, ystar

//make counter global
int histcount = 0;
string mode = "";

//void redraw_histo(TH1D* curhist, Color_t linecol, const char* plotname); //declaration just for using it in copydir()
void redraw_histo(TH1D* curhist, Color_t linecol, string plotname); //declaration just for using it in copydir()
void redraw_histo_ybys(TH1D* curhist, Color_t linecol, string plotname, string yb_bname, string ys_bname);
bool checkname(string &histname, string &expr);

pair<string,string> return_ybys_names(string foldername){
	pair<string,string> ybys_name("0 #leq y_{b} < 1", "0 #leq y* < 1");
	string yb1="yb1", yb2="yb2", ys1="ys1", ys2="ys2";

	if(checkname(foldername, yb1)){
		ybys_name.first = "1 #leq y_{b} < 2";
	}
	else if(checkname(foldername, yb2)){
		ybys_name.first = "2 #leq y_{b} < 3";
	}

	if(checkname(foldername, ys1)){
		ybys_name.second = "1 #leq y* < 2";
	}
	else if(checkname(foldername, ys2)){
		ybys_name.second = "2 #leq y* < 3";
	}

	return ybys_name;
}


pair<string,string> ybys_name;
//function to copy existing directory from input file in outputfile --> add the (newly created) canvases with histos following the original directory structure
//similar to ROOT example: tutorials/io/copyFiles.C
//void copydir(TFile *plotfile, TDirectory *origdir, const char* plotname){
//in addition, this function takes a "summary vector" to store all trigger histos (and plot them later into one canvas)
void copydir(TFile *plotfile, TDirectory *origdir, string plotname, int &histcount){
	cout << "INITIAL HISTCOUNT: " << histcount << endl;
	cout << "ls on origdir: " << endl;
	cout << "---------------" << endl;
	//origdir->ls();		//show directory content
	//plotfile->cd();		//change to output file
	cout << "ls on plotfile: " << endl;
	cout << "----------------" << endl;
	//plotfile->ls();
	TDirectory *plotfiledir = gDirectory; //current directory in outputfile
	//TDirectory *outdir = plotfile->mkdir(origdir->GetName());		//create directory in outputfile which corresponds to dir in input file
	TDirectory *outdir = plotfiledir->mkdir(origdir->GetName());		//create directory in outputfile which corresponds to dir in input file
	outdir->cd();		//change to "output-directory" where current histogram will be saved in

	//check for ybys folder:
	string foldname = (string)outdir->GetName();
	string yb_str = "yb";
	if(checkname(foldname, yb_str)){		//check if current (output) directory is a "ybys" directory
		cout << "Changed mode to 'ybys'" << endl;
		ybys_name = return_ybys_names(foldname);
		mode="ybys";
	}

	//for naming the .pdf file (just append all the ROOT directory-names)
	plotname.append(Form("_%s", outdir->GetName()));
	cout << "plotname is now: >>  " << plotname << endl;

	//now iterate through subdirectories of origdir
	//example: origdir refers to certain Etabin --> subdirectories could refer to different triggers
	//or: origdir refers to Standard directory --> subdirectories could refer to different etabins
	TKey *key;
	TIter nextkey(origdir->GetListOfKeys());		//iterate through keys in original directory (from input file)
	while ((key=(TKey*)nextkey())) {
		const char *classname = key->GetClassName();	//what kind of object is it? --> check class (TH1D, TDirectory, ...)
		TClass *cla = gROOT->GetClass(classname);		//get corresponding class
		if (!cla) continue;
		if (cla->InheritsFrom(TDirectory::Class())){	//if current object is a directory, execute this
			origdir->cd(key->GetName());				//get name of this directory in input file
			//if origdir is called "Standard" or "Eta.." --> change into that directory
			TDirectory *subdir = gDirectory;			//create subdirectory, which will be created in plotfile (output file)
			outdir->cd();								//go to outputfile (main directory)
			cout << "histcount here: " << histcount << endl;
			copydir(plotfile, subdir, plotname, histcount);		//calls itself --> same procedure as before, this time for subdirectory
			outdir->cd();
		}
		//else if(cla->InheritsFrom(TTree::Class())){
		else{	//executed if current object is not a directory (but for instance a histogram)
			origdir->cd();
			TObject *obj = key->ReadObj();	//read object from input directory in inputfile
			outdir->cd();
			obj->Write();					//write object to output directory in outputfile

			//check if it is a histogram that shall be redrawn:
			const char *clname = key->GetClassName();
			TClass *objcla = gROOT->GetClass(clname);
			if(objcla->InheritsFrom(TH1D::Class())){
				cout << "--> Found a TH1D object. " << endl;
				((TH1D*)obj)->SetStats(kFALSE);		//remove stats box before copying (create a new one in redraw_histo())
				++histcount;
				//redraw_histo((TH1D*)obj, kBlue+2, key->GetName());
				//to stop early (only while testing)
				if(histcount>1){continue;}//exit(0);} //JUST FOR TESTING
				if(mode=="ybys"){
					redraw_histo_ybys((TH1D*)obj, kAzure+5, plotname, ybys_name.first, ybys_name.second);
				}
				else{
					redraw_histo((TH1D*)obj, kGreen+2, plotname);
				}
				cout << "Plotted histogram no. " << histcount << endl;
			}	
			cout << "histcount after increasing no. : " << histcount << endl;
			delete obj;
		}
		cout << "histcount after ELSE: " << histcount << endl;
		//return histcount;
	}
	outdir->SaveSelf(kTRUE);
	origdir->cd(); //necessary?
}


//function to create canvas for the existing histos, draw them on the canvas and adjust the plotting style
//currently not general enough --> adjusted to example with only one histogram per trigger and always "events over pT"
//see next function: redraw_histo_ybys()
//void redraw_histos(){
void redraw_histo(TH1D* curhist, Color_t linecol, string plotname){
	
	//set y errors to zero (otherwise they are plotted) [to do: leave option to draw them]
	int nbins = curhist->GetNbinsX();
	for(int k=0; k!=nbins; ++k){
		curhist->SetBinError(k, 0);
	}

	//TCanvas *curcanv = tdrCanvas("curcanv", curhist, 2, 11, kSquare); 
	TCanvas *curcanv = tdrCanvas("curcanv", curhist, 4, 0, kSquare); 
	curcanv->cd();
	curcanv->SetLogx();
	curcanv->SetLogy();
	gStyle->SetOptStat("nemr");	//draw statistics box with histogram name and number of entries
	//gStyle->SetStatStyle(0);	//transparent
	//gStyle->SetStatStyle(111);	//not transparent
	curhist->SetStats(kTRUE);


	//curhist->Draw();
	tdrDraw(curhist, "L, SAME, ][");
	//tdrDraw(curhist, "HIST, ][");
	curhist->SetLineColor(linecol);
	curhist->SetLineWidth(3);
	curhist->SetMarkerSize(0.1);
	//curhist->GetXaxis()->SetMoreLogLabels();
	curhist->GetXaxis()->SetNoExponent();
	curhist->GetXaxis()->SetTitle("p_{T} [GeV]");
	curhist->GetXaxis()->SetTitleSize(0.05);
	curhist->GetXaxis()->SetLabelSize(0.04);
	curhist->GetXaxis()->SetTitleOffset(1.2);
	curhist->GetYaxis()->SetTitle("Events");
	curhist->GetYaxis()->SetTitleSize(0.05);
	curhist->GetYaxis()->SetLabelSize(0.04);
	curhist->GetYaxis()->SetTitleOffset(1.3);
	//curcanv->SetLeftMargin(0.22);
	
	//CMS_logo->Draw("X")		//see tdrStyle file
	//curhist->tdrDraw(self, "");
	///curhist->SetFillColor(kRed);
	///curhist->SetFillStyle(3001);
	curhist->SetFillStyle(0);
	curhist->Write("", TObject::kOverwrite);	//overwrite existing histogram

	gPad->Update();
	TPaveStats *stat = (TPaveStats*)curhist->FindObject("stats");	//get statistics box
	stat->SetX1NDC(0.8);		//set position of statistics box
	stat->SetX2NDC(0.94);
	stat->SetY1NDC(0.8);
	stat->SetY2NDC(0.92);
	//curhist->SetStats(kFALSE);		//remove original statbox
	
	curcanv->Update();
	//curcanv->Print(Form("%s.pdf", plotname));
	curcanv->SaveAs(Form("%s.pdf", plotname.c_str()));
	//delete curcanv;
	curcanv->Close();
}

//for yboost ystar histos
void redraw_histo_ybys(TH1D* curhist, Color_t linecol, string plotname, string yb_bname, string ys_bname){
	//set y errors to zero (otherwise they are plotted) [to do: leave option to draw them]
	int nbins = curhist->GetNbinsX();
	for(int k=0; k!=nbins; ++k){
		curhist->SetBinError(k, 0);
	}

	//TCanvas *curcanv = tdrCanvas("curcanv", curhist, 2, 11, kSquare); 
	TCanvas *curcanv = tdrCanvas("curcanv", curhist, 4, 0, kSquare); 
	curcanv->cd();
	curcanv->SetLogx();
	curcanv->SetLogy();
	curcanv->SetGrid(1,1);
	gStyle->SetOptStat("nemr");	//draw statistics box with histogram name and number of entries
	//gStyle->SetStatStyle(0);	//transparent
	//gStyle->SetStatStyle(111);	//not transparent
	curhist->SetStats(kTRUE);


	//curhist->Draw();
	tdrDraw(curhist, "E, SAME, ][");
	//tdrDraw(curhist, "HIST, ][");
	curhist->SetLineColor(linecol);
	curhist->SetLineWidth(3);
	curhist->SetMarkerSize(0.1);
	gPad->SetGrid(1,1);
	//curhist->GetXaxis()->SetMoreLogLabels();
	curhist->GetXaxis()->SetNoExponent();
	curhist->GetXaxis()->SetTitle("p_{T,avg} [GeV]");
	curhist->GetXaxis()->SetTitleSize(0.05);
	curhist->GetXaxis()->SetLabelSize(0.04);
	curhist->GetXaxis()->SetTitleOffset(1.2);
	curhist->GetYaxis()->SetTitle("Events");
	curhist->GetYaxis()->SetTitleSize(0.05);
	curhist->GetYaxis()->SetLabelSize(0.04);
	curhist->GetYaxis()->SetTitleOffset(1.3);
	//curcanv->SetLeftMargin(0.22);
	
	//CMS_logo->Draw("X")		//see tdrStyle file
	//curhist->tdrDraw(self, "");
	///curhist->SetFillColor(kRed);
	///curhist->SetFillStyle(3001);
	curhist->SetFillStyle(0);
	curhist->Write("", TObject::kOverwrite);	//overwrite existing histogram

	gPad->Update();
	TPaveStats *stat = (TPaveStats*)curhist->FindObject("stats");	//get statistics box
	stat->SetX1NDC(0.8);		//set position of statistics box
	stat->SetX2NDC(0.94);
	stat->SetY1NDC(0.8);
	stat->SetY2NDC(0.92);
	//curhist->SetStats(kFALSE);		//remove original statbox

	TLatex* bininfo = new TLatex();                                                                                                                                                                                     
    bininfo->SetNDC();
    bininfo->SetTextSize(0.045);
	bininfo->SetTextFont(42);
    bininfo->DrawLatex(0.20, 0.86, yb_bname.c_str());
    bininfo->DrawLatex(0.20, 0.80, ys_bname.c_str());
	
	curcanv->Update();
	//curcanv->Print(Form("%s.pdf", plotname));
	curcanv->SaveAs(Form("%s.pdf", plotname.c_str()));
	//delete curcanv;
	curcanv->Close();
}

bool checkname(string &histname, string &expr){
	cout << Form("look for %s in %s. This is returned: --> ", expr.c_str(), histname.c_str()) << endl;
	cout << "--> " << boolalpha << (histname.find(expr)!=string::npos) << endl;
	return histname.find(expr)!=string::npos;
}

//vector<TH1D*> gethistos(TDirectory *curdir, string &expr, vector<TH1D*> &histvec){
void gethistos(TDirectory *curdir, string &expr, vector<TH1D*> &histvec){
	//find all the histograms in given directory curdir which include expression expr
	curdir->cd();
	cout << endl << "ls on curdir = " << curdir->GetName() << endl;
	cout << "--------------------------------" << endl;
	curdir->ls();

	//vector<TH1D*> histvec;
	TKey *key;
	TIter nextkey(curdir->GetListOfKeys());
	while((key=(TKey*)nextkey())){ 		//goes on until there are no keys left
		const char *classname = key->GetClassName();
		TClass *cla = gROOT->GetClass(classname);
		if(!cla) continue;
		if(cla->InheritsFrom(TH1D::Class())){
			TH1D* newhist = (TH1D*)key->ReadObj();
			newhist->SetStats(kFALSE);			//remove stats box
			string histname = (string)newhist->GetName();
			cout << "This is a histogram." << endl;
			if(checkname(histname, expr)){
				histvec.push_back((TH1D*)newhist->Clone());
			}
			else{
				cout << "Not a histogram that has been looked for." << endl;	
			}
			delete newhist;
		}
		else if(cla->InheritsFrom(TDirectory::Class())){
			curdir->cd(key->GetName());
			TDirectory *subdir = gDirectory;
			gethistos(subdir, expr, histvec);
		}
		else{continue;}
	}
	//return histvec;
}


//function for the trigger summary (turn-on-curves) ==>> TO BE IMPROVED: do not only save this as canvas --> also as histograms! (more options in TBrowser)
//go through outputfile, this function takes a certain (Eta-) directory as parameter and summarises all histograms from there
//void trigger_summary(vector<TH1D*> histvec, vector<string> histnames){
//void trigger_summary(TDirectory *curdir, vector<TH1D*> &histvec){
void trigger_summary(vector<TH1D*> &histvec, const char* title){		//use directory only for naming.. could instead
	//vector<Style_t> markerstyles = {kFullCircle, kFullSquare, kFullTriangleUp, kFullDiamond, kFullStar, kFullCross, kFullCrossX};
	vector<Style_t> markerstyles = {20, 21, 22, 33, 29, 34, 47, 23};
	//histvec.push_back(((TH1D*)obj));	//watch out when generalising --> this is adjusted to "only triggers"-case

	//customised binning
	Double_t cust_bb[] = {10, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507, 548, 592, 638, 686, 737, 790, 846, 905, 967, 1032, 1101, 1172, 1248, 1327, 1410, 1497, 1588, 1684, 1784, 1890, 2000, 2116, 2238, 2366, 2500, 2640, 2787, 2941, 3103, 3273, 3450, 3832, 6076, 6389};

	TH1D* dummyhist = new TH1D("dummyhist", ";p_{T};Events", 64, cust_bb);
	//TH1D* dummyhist = new TH1D("dummyhist", ";p_{T};Events", 5990, 10, 6000);
	TCanvas *sumcanv = tdrCanvas("sumcanv", dummyhist, 2, 0, kSquare);
	gStyle->SetPalette(kStarryNight); //how big are the steps in colorpalette? (from histo to histo)

	//dummyhist->SetAxisRange(0.1, 1E5, "Y");	//avoid such hard coded ranges...
	//loop to find suitable range for histogram: (should initial values be changed?)
	Double_t minc=0.1;
	Double_t maxc=100;		//don't choose this too high
	for(int ind=0; ind!=histvec.size(); ++ind){
		cout << "minc= " << minc << " and maxc= " << maxc << endl;
		//cout.precision(20);
		cout << "histogram: " << histvec.at(ind)->GetName() << endl;
		cout << histvec.at(ind)->GetBinContent(histvec.at(ind)->GetMinimumBin()) << endl;
		//Double_t tmp_min = histvec.at(ind)->GetMinimum();	//get lowest value in current histo (often zero, so not useful here)
		Double_t tmp_max = histvec.at(ind)->GetMaximum();	//get highest value in current histo

		int nbins = histvec.at(ind)->GetNbinsX();
		for(int bini=1; bini!=nbins+1; ++bini){	//get lowest (non-zero) value
			Double_t bincon = histvec.at(ind)->GetBinContent(bini);	//get bin content
			if(bincon!=0 and bincon<minc){		//check if bin content is lower than current minimum value
				minc = bincon;
			}
		}
		
		//if(tmp_min<minc){
		//	minc = tmp_min;
		//}

		if(tmp_max>maxc){	//check if highest bin content in current histogram is higher than current maximum value
			maxc = tmp_max;
		}
		//no "continue", as this loop is so short
	}
	//calculate log_10 of minimum bin content
	Double_t potmin = TMath::Log10(minc);
	Double_t potmax = TMath::Log10(maxc);
	cout << "After calculation: potmin= " << potmin << " and potmax= " << potmax << endl;
	potmin = floor(potmin);				//round power down
	potmax = ceil(potmax);					//round power up

	//if(potmin<0){potmin=floor(potmin);}
	//else{potmin=ceil(potmin);}
	cout << "POTMIN: " << potmin << endl;
	cout << "POTMAX: " << potmax << endl;
	Double_t ymin = pow(10, potmin);
	Double_t ymax = pow(10, potmax);
	dummyhist->SetAxisRange(ymin, ymax, "Y");	//set y-range according to histogram contents
	cout << "&&&&&&&&&&&&&&&&&&&&&" << endl;
	cout << "ymin= " << ymin << "  and ymax= " << ymax << endl;
	sumcanv->SetLogx();
	sumcanv->SetLogy();
	//sumcanv->SetGrid();
	//gPad->SetGrid();
	//dummyhist->GetXaxis()->SetMoreLogLabels();
	dummyhist->GetXaxis()->SetNoExponent();
	dummyhist->SetTitle(title);		//indicating for which bin this is the summary


	TLegend* leg = tdrLeg(0.7, 0.2, 0.9, 0.5);
	leg->SetTextFont(42);
	
	TLatex *bininfo = new TLatex();
	bininfo->SetNDC();
	bininfo->SetTextSize(0.045);
	bininfo->DrawLatex(0.7, 0.82, Form("#font[42]{%s}", title)); //font 12 looks nicer, but does not get saved correctly...

	sumcanv->cd();
	cout << "Size of histvec: " << histvec.size() << endl;
	for(int ind=0; ind!=histvec.size(); ++ind){
		//histvec.at(ind)->Draw("SAME PLC PMC");
		histvec.at(ind)->SetLineWidth(2);
		histvec.at(ind)->Draw("SAME PLC PMC HIST PC");
		histvec.at(ind)->Print();
		if(ind<8){
			histvec.at(ind)->SetMarkerStyle(markerstyles.at(ind));
		}
		else{
			histvec.at(ind)->SetMarkerStyle(markerstyles.at(ind-8));
		}
		histvec.at(ind)->SetMarkerSize(1);
		leg->AddEntry(histvec.at(ind), histvec.at(ind)->GetName(), "lp");
		cout << "Current histo's name: " << histvec.at(ind)->GetName() << endl;
		cout << "Bin center of bin_index=1: --> " << histvec.at(ind)->GetXaxis()->GetBinCenter(1) << endl;
		cout << "added new histo to summary!" << endl;
		sumcanv->Draw(); //test if it shows
		sumcanv->Modified();
		sumcanv->Update();
		//sumcanv->SaveAs(Form("test_index_%d.pdf", ind));
		//histvec.at(ind)->SaveAs(Form("%s.pdf", (histnames.at(ind)).c_str()));
	}

	//Add grid manually.
	TPad *grid = new TPad("grid", "", 0, 0, 1, 1);
	grid->Draw("SAME");
	grid->cd();
	grid->SetGrid();
	grid->SetFillStyle(4000);

	leg->Draw("SAME");
	leg->SetTextSize(0.025);
	cout << "legend: GetTextFont(): " << leg->GetTextFont() << endl;
	//leg->SetTextFont(132);
	leg->SetTextFont(12);
	sumcanv->cd();
	//sumcanv->SetGrid(1,1);			//switch on grid
	sumcanv->Modified();
	sumcanv->Update();
	sumcanv->Draw();
	sumcanv->SaveAs(Form("summary_trigger_%s.pdf", title));
	sumcanv->Modified();
	sumcanv->Update();
	sumcanv->Write(Form("summary_canvas_%s", title)); //for testing
	//sumcanv->Close();
}


//function for plotting trigger efficiency
//takes unscaled histograms from root outputfile, scales each by its corresponding integrated luminosity
//devides it by the previous (scaled) triggerhistogram
void trg_lumi(vector<TH1D*> &histvec, vector<Double_t> &lumivec, const char* bintitle){
	//marker styles, as there are many graphs in one single plot
	vector<Style_t> markerstyles = {20, 21, 22, 33, 29, 34, 47, 23};
	//customised binning
	Double_t cust_bb[] = {10, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507, 548, 592, 638, 686, 737, 790, 846, 905, 967, 1032, 1101, 1172, 1248, 1327, 1410, 1497, 1588, 1684, 1784, 1890, 2000, 2116, 2238, 2366, 2500, 2640, 2787, 2941, 3103, 3273, 3450, 3832, 6076, 6389};

	TH1D* dummyhist = new TH1D("dummyhist", ";p_{T,avg};#frac{trg_{N} / lumi_{N}}{trg_{N-1} / lumi_{N-1}}", 64, cust_bb);
	//TH1D* dummyhist = new TH1D("dummyhist", ";p_{T,avg};#frac{trg_{N} / lumi_{N}}{trg_{N-1} / lumi_{N-1}}", 5990, 10, 6000); 
	TCanvas *sumcanv = tdrCanvas("sumcanv", dummyhist, 2, 0, kSquare);

	gStyle->SetPalette(kStarryNight); //how big are the steps in colorpalette? (from histo to histo)

	
	TLegend* leg = tdrLeg(0.7, 0.2, 0.9, 0.5);
	//TLegend* leg = tdrLeg(1.01, 0.2, 1.21, 0.5);
	leg->SetTextFont(42);
	
	TLatex *bininfo = new TLatex();
	bininfo->SetNDC();
	bininfo->SetTextSize(0.045);
	bininfo->DrawLatex(0.7, 0.82, Form("#font[42]{%s}", bintitle)); //font 12 looks nicer, but does not get saved correctly...

	sumcanv->cd();
	cout << "Size of histvec: " << histvec.size() << endl;
	
	//check y-values minimum and maximum
	Double_t minc=0.1;
	Double_t maxc=1; 	//don't choose this too high, as it will be checked if other maximum is HIGHER than this
	//starting from "1", because need to divide by the (N-1)th ("0"th) histogram
	//starting from "2", because SEE ABOVE and first histogram is jt0 --> no trigger --> leave this out here
	for(int ind=2; ind!=histvec.size(); ++ind){
	///for(int ind=2; ind!=5; ++ind){	//test
		//histvec.at(ind)->Draw("SAME PLC PMC");
		cout << "CURRENT INDEX: " << ind << endl;
		//TH1D* prevhist = (TH1D*)histvec.at(ind-1)->Clone(Form("prevhist_%d", ind));	//previous histogram (N-1)
		//TH1D* divhist = (TH1D*)histvec.at(ind)->Clone(Form("divhist_%d", ind));	//current histogram (N), that will be divided by hist(N-1)
		TH1D* prevhist = (TH1D*)histvec.at(ind-1)->Clone(Form("prevhist_%d", ind));	//previous histogram (N-1)
		TH1D* divhist = (TH1D*)histvec.at(ind)->Clone(Form("divhist_%d", ind));	//current histogram (N), that will be divided by hist(N-1)
		divhist->SetLineWidth(2);
		prevhist->Scale(1./(lumivec.at(ind-2)));	//scale previous histogram by lumi
		divhist->Scale(1./lumivec.at(ind-1));			//scale current histo by lumi
		divhist->Divide(prevhist);					//divide current histogram by previous one (both already scaled by lumi)

		//for finding suitable y-axis range
		Double_t tmp_max = divhist->GetMaximum();	//get highest value in current ratiohisto
		int nbins = divhist->GetNbinsX();
		for(int bini=1; bini!=nbins+1; ++bini){
			Double_t bincon = divhist->GetBinContent(bini);
			if(bincon!=0 and bincon<minc){ 						//check if bin content is lower than current minimum value, is allowed to be zero, because y-axis is not log.
				minc = bincon;			
			}
		}
		if(tmp_max>maxc){							//check if highest bin content in current histogram is higher than current maximum y-value
			maxc = tmp_max;
		}
		
		divhist->Draw("SAME PLC PMC HIST PC");
		divhist->Print();
		if(ind<9){
			divhist->SetMarkerStyle(markerstyles.at(ind-1));
		}
		else{
			divhist->SetMarkerStyle(markerstyles.at(ind-9));
		}
		divhist->SetMarkerSize(1);
		leg->AddEntry(divhist, Form("%s / %s", divhist->GetName(), prevhist->GetName()), "lp");
		cout << "Current histo's name: " << divhist->GetName() << endl;
		cout << "Bin center of bin_index=1: --> " << divhist->GetXaxis()->GetBinCenter(1) << endl;
		cout << "added new histo to summary!" << endl;
		sumcanv->Draw(); //test if it shows
		sumcanv->Modified();
		sumcanv->Update();
		//sumcanv->SaveAs(Form("test_index_%d.pdf", ind));
		//histvec.at(ind)->SaveAs(Form("%s.pdf", (histnames.at(ind)).c_str()));
	}

	//calculate log_10 of minimum bin content
	Double_t potmin = TMath::Log10(minc);
	Double_t potmax = TMath::Log10(maxc);
	cout << "After calculation: potmin= " << potmin << " and potmax= " << potmax << endl;
	potmin = floor(potmin);				//round power down
	potmax = ceil(potmax);					//round power up

	//if(potmin<0){potmin=floor(potmin);}
	//else{potmin=ceil(potmin);}
	cout << "POTMIN: " << potmin << endl;
	cout << "POTMAX: " << potmax << endl;
	Double_t ymin = pow(10, potmin);
	Double_t ymax = pow(10, potmax);
	dummyhist->SetAxisRange(ymin, ymax, "Y");	//set y-range according to histogram contents
	

	//dummyhist->SetAxisRange(minc, maxc, "Y");	//set y-range according to histogram contents
	sumcanv->SetLogx();
	sumcanv->SetLogy();
	//dummyhist->GetXaxis()->SetMoreLogLabels();
	dummyhist->GetXaxis()->SetNoExponent();
	dummyhist->SetTitle(bintitle);		//indicating for which bin this is the summary


	//Add horizontal line at y=1.0
	TLine *line = new TLine(sumcanv->GetUxmin(), 1.0, sumcanv->GetUxmax(), 1.0);
	line->SetLineColorAlpha(kAzure+5, 0.36);
	line->SetLineStyle(9);
	line->Draw("SAME");


	//Add grid manually.
	TPad *grid = new TPad("grid", "", 0, 0, 1, 1);
	grid->Draw("SAME");
	grid->cd();
	grid->SetGrid();
	grid->SetFillStyle(4000);

	leg->Draw("SAME");
	leg->SetTextSize(0.025);
	cout << "legend: GetTextFont(): " << leg->GetTextFont() << endl;
	//leg->SetTextFont(132);
	leg->SetTextFont(12);
	sumcanv->cd();
	//sumcanv->SetGrid(1,1);			//switch on grid
	//dummyhist->GetYaxis()->SetTitle("Events");
	dummyhist->GetYaxis()->SetTitleSize(0.03);
	dummyhist->GetYaxis()->SetLabelSize(0.04);
	dummyhist->GetYaxis()->SetTitleOffset(2.6);
	dummyhist->GetXaxis()->SetTitleSize(0.03);
	dummyhist->GetXaxis()->SetLabelSize(0.04);
	sumcanv->SetLeftMargin(0.22);
	//sumcanv->SetRightMargin(0.16);
	sumcanv->Modified();
	sumcanv->Update();
	sumcanv->Draw();
	sumcanv->SaveAs(Form("summary_trglumi_%s.pdf", bintitle));
	sumcanv->Modified();
	sumcanv->Update();
	sumcanv->Write(Form("summary_trglumi_canvas_%s", bintitle)); //for testing
	//sumcanv->Close();


}//end function trg_lumi





//new:
//void HistoPlotting(){
void DijetPlotting(){
	setTDRStyle();
	extraText="private work";

	// Create file for saving the histograms
	TFile* plotfile = new TFile("plots_test_new.root", "RECREATE");

	// Open output file
	TFile* inputfile = (TFile*)gROOT->GetListOfFiles()->FindObject("interim_output_new.root"); //could be handled via input? (so that plotting can be used for any root-file)
	cout << "after inputfile " << endl;
	cout << "&inputfile " << &inputfile << endl;
	if (!inputfile || !inputfile->IsOpen()){
		//inputfile = new TFile("interim_output.root", "READ"); //does not make sense to create this file, want to read from it!
		cout << "first if condition " << endl;
		inputfile = TFile::Open("interim_output_new.root", "READ");
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
	


	//create trigger-summary plot for each eta-bin (should be generalised, this is just test version)
	vector<const char*> etabins = {"Eta_00_05", "Eta_05_10", "Eta_10_15", "Eta_15_20", "Eta_20_25", "Eta_25_30", "Eta_30_32",  "Eta_32_47"};

	//list of ybys bins
	vector<const char*> ybysbins = {"yb0ys0", "yb0ys1", "yb0ys2", "yb1ys0", "yb1ys1", "yb2ys0"};
	
	/*
	vector<TDirectory*> etadirs;
	vector<const char*> etadirnames;
	for(int k=0, k!=etabins.size(), ++k){
		TDirectory* etadir; 
		plotfile->GetObject(Form("Standard/%s", etabins.at(k)), etadir);
		etadirs.push_back(etadir);
		etadirnames.push_back(etadir->GetName());
	}
	*/

	if(mode!="ybys"){	//only executed for other modes than ybys
		for(int k=0; k!=etabins.size(); ++k){
			TDirectory* etadir; 
			plotfile->GetObject(Form("Standard/%s", etabins.at(k)), etadir);
			const char* etaname = etadir->GetName();
			string expr = "hpt";	//expression that has to be contained in histogram-name in order for histo to appear in summary plot
			vector<TH1D*> histosvec;
			gethistos(etadir, expr, histosvec);	//fills histograms in vector
		
			plotfile->cd();
			trigger_summary(histosvec, etaname);
		}
	}
	else if (mode=="ybys"){
		vector<Double_t> lumivec = {0, 0, 0.000414643, 0.001737711, 0.006541972, 0.097552471, 0.418082018, 0.868855667, 3.475422669, 6.950845338, 13.901690676, 111.213525407, 111.213525407};	//now for testing --> should be more automatised / generalised [units are /pb]
		for(int k=0; k!=ybysbins.size(); ++k){
			TDirectory* ybysdir;
			plotfile->GetObject(Form("Standard/%s", ybysbins.at(k)), ybysdir);
			const char* ybysname = ybysdir->GetName();
			string expr = "hpt";
			vector <TH1D*> histosvec;
			gethistos(ybysdir, expr, histosvec);

			plotfile->cd();
			trg_lumi(histosvec, lumivec, ybysname);
		}
	}


	/*
	string expr = "hpt";	//expression that has to be contained in histogram-name in order for histo to appear in summary plot
	//vector<TH1D*> histosvec = gethistos(etadir, expr);
	vector<TH1D*> histosvec; 
	gethistos(etadir, expr, histosvec);	//fills histograms in vector
	
	plotfile->cd();
	trigger_summary(histosvec, etaname);

	//TCanvas* sumcanv = trigger_summary(etadir);
	//sumcanv->SaveAs(Form("summary_trigger.pdf"));
	//sumcanv->Close();
	*/


	inputfile->Close();		//close file
	plotfile->Close();
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


