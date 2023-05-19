// nice dictionary tutorials
// https://github.com/eguiraud/root_dictionaries_tutorial
// https://root.cern/manual/integrate_root_into_my_cmake_project/

#include "TFile.h"
#include "TTree.h"
#include "TInterpreter.h"
#include <vector>
#include <iostream>
#include "TROOT.h"
#include "TApplication.h"
#include "TCanvas.h"
#include "TObjArray.h"

#include "TProfile.h"
#include "TF1.h"
#include "TH2.h"
#include "TLine.h"



using std::cout;
using std::endl;

void pulse_heightDRS(TString file);

int main(int argc, char **argv) {

  if (argc<2){
    cout << "No file name entered" << endl;
    return(1);
  }
  TString file(argv[1]);

  // This allows you to view ROOT-based graphics in your C++ program
  // If you don't want view graphics (eg just read/process/write data files), 
  // this can be ignored
  TApplication theApp("App", &argc, argv);

  pulse_heightDRS(file);
  
  
  // view graphics in ROOT if we are in an interactive session
  if (!gROOT->IsBatch()) {
    std::cout << "To exit, quit ROOT from the File menu of the plot (or use control-C)" << std::endl;
    theApp.SetIdleTimer(30,".q");  // set up a failsafe timer to end the program
    theApp.Run(true);
  }
  return 0;

}


void pulse_heightDRS(TString file){
  cout << "Processing: " << file << endl;

  auto tf=new TFile(file);
  auto tree=(TTree*)tf->Get("waves");

  // figure out the signal pulatity and which channel to plot from the filename
  double polarity=1.0;
  if (file.Contains("TIA")) polarity=-1.0;
  int ampid=2;
  if (file.Contains("C6")) ampid=6;
  if (file.Contains("C4")) ampid=4;
  if (file.Contains("C0")) ampid=0;
  
  std::vector<std::vector<float> > *chs=0;
  std::vector<float> *ttime=0;

  float min;
  Long_t event;
  tree->SetBranchAddress("chs", &chs);
  tree->SetBranchAddress("time",&ttime);
  tree->SetBranchAddress("min",&min);


  auto tcsum = new TCanvas("tcsum","summary");
  tcsum->Divide(2,2);
  tcsum->cd(1);
  tree->Draw("chs[8]:time");

  tcsum->cd(2);
  TString ampst = TString::Format("chs[%d]:time",ampid);
  if (polarity<0) tree->Draw("-1.0*"+ampst);
  else tree->Draw(ampst);

  // demonstrate alignment scheme and determine the mean pulse shape
  int TrigSet=820;  // fix nominal trigger position
  int LEN=1000;     // depends on data file, make this dynamic in the future
  auto htrig = new TH2F("htrig","Trigger aligned;sample no.;ADC",LEN,0,LEN,1300,1000,2300);
  auto hprofv= new TProfile("hprofv","Average waveform;sample no.;ADC",LEN,0,LEN);

  
  for (int i=0; i<tree->GetEntries(); ++i){
    tree->GetEntry(i);
    // find shift in peak location due to trigger
    int imin = int(min);
    int dt=imin-TrigSet;
    // cout << dt << endl;
    for (int n=0; n<LEN; ++n) {
      if ((n+dt)>=0 && (n+dt)<LEN){
	htrig->Fill(n,(*chs)[8][n+dt]);
	hprofv->Fill(n,polarity*((*chs)[ampid][n+dt]));
      }
    }
  }
  tcsum->cd(3);
  htrig->DrawCopy();

  
  tcsum->cd(4);


  // find peak and baseline
  int ipeak =  hprofv->GetMinimumBin();
  cout << "peak at " << ipeak << endl;

  // start end window for pulse integral
  int istart = ipeak-25;  // tuned by hand
  int istop = ipeak+25;

  hprofv->Fit("pol0","","",0,istart);
  hprofv->DrawCopy();
  double baseline = hprofv->GetFunction("pol0")->GetParameter(0);
  double ymin = hprofv->GetMinimum();
  double ymax = hprofv->GetMaximum();
  int meanpulse =  abs(baseline-hprofv->GetMinimum());

  
  auto l1= new TLine();
  l1->DrawLine(istart,ymin,istart,ymax);
  auto l2= new TLine (istop,ymin,istop,ymax);
  l2->Draw();


  auto tcsum2 = new TCanvas("tcsum2","summary2");
  tcsum2->Divide(1,2);
  
  auto phd = new TH1F("phd","PulseHeights;ADC;frequency",3*meanpulse+50,-50,3*meanpulse);
  auto pid = new TH1F("pid","PulseIntegral;A.U.;frequency",meanpulse*(istop-istart)/40,
		      -meanpulse*(istop-istart)/10,meanpulse*(istop-istart));

   
  for (int ievt=0; ievt<tree->GetEntries(); ++ievt){
    tree->GetEntry(ievt);
    int imin = int(min);
    int dt=imin-TrigSet;
    double height=-polarity*((*chs)[ampid][ipeak+dt]-polarity*baseline);
    phd->Fill(height);
    double sum=0;

    for (int n=istart; n<istop; ++n) {
      double V=-polarity*((*chs)[ampid][n+dt]-polarity*baseline);
      sum+=(V); 
    }
    pid->Fill(sum);
  }

  tcsum2->cd(1);
  phd->DrawCopy();
  tcsum2->cd(2);
  pid->DrawCopy();

  TString plotfile=file;
  int idx=plotfile.Last('/');
  if ( idx>-1 ) plotfile.Remove(0,idx+1);
  TString outfile=plotfile;
  outfile.ReplaceAll(".root","_out.root");
 
  // generate output file names
  // note PDFs are HUGE b/c of 2d scatter plots, so use png
  plotfile.ReplaceAll(".root","_1.png");
  tcsum->Print(plotfile);
  plotfile.ReplaceAll("_1.png","_2.png");
  tcsum2->Print(plotfile);

  auto tfout = new TFile(outfile,"RECREATE");

  phd->Write();
  pid->Write();
  hprofv->Write();
  tfout->Write();
  tfout->Close();
}  






