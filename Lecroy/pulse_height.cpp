// 400 samples left of peak

#include "TMath.h"

void pulse_height(TString file="C2--CV-40_54V-partslaser--5k--00000.root",int ch=3,
		  bool savePlot=false){

  const double dTbaseline = 30e-9;  // time to sample before peak (new amplifier)
  const double integrateT = 10e-9;  // integrate from tPeak-dTbaseline to tPeak+integrateT
  auto tf=new TFile(file);
  auto tree=(TTree*)tf->Get("waves");
  // output file
  TString fnout=file;
  fnout.ReplaceAll(".root",".pdf");
  int pos=fnout.Last('/');
  if ( pos>-1 ) fnout.Remove(0,pos+1);

  gStyle->SetOptStat(0);

  // find the size of the sample buffer
  auto br=tree->GetBranch("time");
  auto tl=br->GetLeaf("time");
  int LEN = tl->GetNdata();
  cout << "Processing buffers of length: " << LEN << endl;
  cout << "Number of buffers: " << tree->GetEntries() << endl;
  
  Double_t *time = new Double_t[LEN];
  Double_t *pool = new Double_t[LEN*4]; 
  Double_t **volts = new Double_t *[4];
  for (int i = 0; i < 4; i++) {
    volts[i] = &(pool[LEN*i]);
  }


  Double_t startx;
  Long_t event;
  tree->SetBranchAddress("time",time);
  tree->SetBranchAddress("volts",pool);
  tree->SetBranchAddress("startx",&startx);

  
  auto tcsum = new TCanvas("tcsum","summary");
  tcsum->Divide(2,2);
  tcsum->cd(1);
  TString m;
  m.Form("volts[%d]*1000:(time-startx)*1e9", ch);
  std::cout << m << std::endl; 
  tree->Draw(m,"event==0","",1);
  //tree->Draw("volts[3]*1000:(time-startx)*1e9","event==0","",1);
  double sampleTime=0;

  // determine the pulse shape
  auto hprof= new TProfile("hprof","Average waveform;sample no.;mV",LEN,0,LEN);
  for (int i=0; i<tree->GetEntries(); ++i){
    tree->GetEntry(i);
    for (int n=0; n<LEN; ++n) hprof->Fill(n,volts[ch][n]*1000);
    if (i==0) sampleTime = time[1]-time[0];
  }
  cout << "sampling time = " << sampleTime << endl;

  // find baseline and peak
  tcsum->cd(2);
  hprof->Fit("pol0","","",0,1000);
  double baseline = hprof->GetFunction("pol0")->GetParameter(0);
  double min = hprof->GetMinimum();
  double max = hprof->GetMaximum();
  double polarity=-1.0;
  int ipeak = hprof->GetMinimumBin();
  if (fabs(max)>fabs(min)) {
   ipeak = hprof->GetMaximumBin();
    polarity=1.0;
  }
  double ymax=fabs(volts[ch][ipeak]);
  
  // start end window for pulse integral
  //int istart = ipeak-90;  // tuned by hand
  //int istop = ipeak+100;

  // # bins before the peak  to do the baseline subtraction
  int iBLS = (int)(dTbaseline/sampleTime);
  // # bins after the peak to integrate
  int iInt = (int)(integrateT/sampleTime);

  int istart = ipeak-iBLS;  // tuned by hand
  int istop = ipeak+iInt;
  auto l1= new TLine(istart,min,istart,max);
  l1->SetLineStyle(2);
  l1->Draw();
  auto l2= new TLine (istop,min,istop,max);
  l2->SetLineStyle(2);
  l2->Draw();
  auto l3= new TLine (ipeak,min,ipeak,max);
  l3->SetLineStyle(2);
  l3->SetLineColor(kRed);
  l3->Draw();
    
  // find scale factor for the integral to normalize it to the pHd
  double iScale = hprof->GetMaximum() / hprof->Integral(istart,istop);



  TString plotfile=file;
  int loc=plotfile.Last('/');
  if ( loc>-1 ) plotfile.Remove(0,loc+1);
  TString outfile=plotfile;
  outfile.ReplaceAll(".root","_out.root");
  auto tfout = new TFile(outfile,"RECREATE");

  //auto phd = new TH1F("phd","PulseHeights;mV;frequency",160,-0.5,ymax*1000*2);
  auto phd = new TH1F("phd","PulseHeights;mV;frequency",256,-10.0,80);
  //auto pid = new TH1F("pid","PulseIntegral;A.U.;frequency",150,-200,4000);
  //auto pid = new TH1F("pid","PulseIntegral;A.U.;frequency",150,-200,30000);
  auto pid = (TH1F*)(phd->Clone("pid"));
  pid->SetTitle("PulseIntegral/%Delta t;A.U.;frequency");

  //double V0=0;
  //double Vmax=4;
  //int Vbins=500;
  //double Vwid=(Vmax-V0)/Vbins;
  //auto hscan = new TH1I("hscan","Threshold Scan;mV;frequency",Vbins,V0,Vmax);
   
  
  for (int ievt=0; ievt<tree->GetEntries(); ++ievt){
    tree->GetEntry(ievt);
    baseline = volts[ch][ipeak-iBLS]*1000;
    double height=polarity*(volts[ch][ipeak]*1000-baseline);
    phd->Fill(height);
    double sum=0;
    // full pulse integral histogram and threshold scan
    //    for (int n=0; n<LEN; ++n){
    //      double V=-(volts[ch][n]*1000-baseline);
    //      if (n>=istart && n<istop) sum+=(V);      
    //      if (V<V0 || V>=Vmax) continue;
    //      int ibin = (int)((V-V0)/Vwid) + 1;
    //      for (int ib=1; ib<=ibin; ++ib) hscan->AddBinContent(ib);
    //    }

    for (int n=istart; n<istop; ++n) {
      double V=polarity*(volts[ch][n]*1000-baseline);
      sum+=(V); 
    }
    pid->Fill(sum * iScale);   // new
  }

  tcsum->cd(3);
  phd->DrawCopy();
  tcsum->cd(4);
  pid->DrawCopy();

  if (savePlot) tcsum->SaveAs(fnout);
  
  delete[] time;
  delete[] pool; 

  tfout->Write();
  tfout->Close();
  
  return;
  
  // DCR analysis
  // loop over N buffers at a time to read a long sample window
  // Use TSpectrum to extract peaks
  int NBUF=50;
  auto buffer = new TH1F("buffer","buffer",LEN*NBUF,0,LEN*NBUF);
  int bbins=buffer->GetNbinsX();
  int offset=0;
  int ievt=0;
  while (ievt<NBUF){
    tree->GetEntry(ievt);
    for (int n=0; n<LEN; ++n){
       double V=-(volts[ch][n]*1000-baseline);
       buffer->SetBinContent(n+offset,V);
    }
    offset+=LEN;
    ievt++;
  }

  auto tcbuf=new TCanvas("rawbuffer","raw buffer");

  double *source = new double[bbins];
  for (int n=0; n<bbins; ++n) source[n]=buffer->GetBinContent(n+1);
  auto back=(TH1F*)buffer->Clone();
  auto TS=new TSpectrum(100*NBUF);

  TS->Background(source,bbins,150,TSpectrum::kBackDecreasingWindow,
		TSpectrum::kBackOrder2,kTRUE,
		TSpectrum::kBackSmoothing15,kFALSE);
  
   for (int n = 0; n < bbins; n++) back->SetBinContent(n + 1,source[n]+baseline);  //HACK
   back->SetLineColor(kRed);
   buffer->Draw();
   back->Draw("SAME L");
   delete[] source;
   auto tcbuf2=new TCanvas("Backgroundsubtracted","Background subtracted");
   auto buffer2 = (TH1F*)buffer->Clone("buffer2");
   buffer2->Add(back,-1.0);
   buffer2->Draw();

   int nrebin=20;
   buffer2->Rebin(nrebin);
   int npeaks = TS->Search(buffer2,2,"nobackground noMarkov",0.05);
   cout << npeaks << " peaks found" << endl;
   Double_t *xpeaks;
   xpeaks = TS->GetPositionX();
   Double_t *ypeaks;
   ypeaks = TS->GetPositionY();
   /// hacks ///
   double ymin = 6; 
   double t0=0;
   double dt=1;   // ns
   ///
   /// sort in time
   auto idx = new int[npeaks];
   TMath::Sort(npeaks,xpeaks,idx,false);

   auto hdt=new TH1I("hdt","hdt",100,0,3000);
   double last=xpeaks[idx[0]];
   double sum=0;
   int count=0;
   for (int n=1; n<npeaks; ++n){
     int i=idx[n];
     if (ypeaks[i]<ymin) continue;
     double dt=(xpeaks[i]-last)*nrebin*0.05;
     //cout << last << " " << xpeaks[i] << endl;
     hdt->Fill(dt);
     sum+=dt;
     count++;
     last=xpeaks[i];
   }
   
   auto tdt=new TCanvas("dt","dt");
   hdt->Fit("expo","","",300,3000);

   double meanDC=sum/count * nrebin * 0.05;
   cout << " Mean dt = " << sum/count * nrebin * 0.05 << " ns" << endl;
   cout << " DCR = " << 1/(meanDC*1e-9) << " Hz" << endl;
  
}

