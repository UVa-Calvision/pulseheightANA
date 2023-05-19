import ROOT as r
from math import sqrt
from math import log
from math import pi
    
class PeakFitter():
    def __init__(self,hPhD):
        self.hPhD=hPhD
        self.ts=r.TSpectrum()
        self.xyPeaks=[]

    def fft(self):
        self.hfft=self.hPhD.FFT(0,"RE")
        self.hfft.SetBinContent(1,0)  # suppress DC component
        self.hfft.SetBinContent(2,0)
        self.hfft.SetBinContent(self.hfft.GetNbinsX(),0)
        self.hfft.SetBinContent(self.hfft.GetNbinsX()-1,0)
        max=self.hfft.GetMaximumBin()
        # respect the Nyquist limit
        if max>self.hfft.GetNbinsX()/2: max = self.hfft.GetNbinsX()-max
        rangefft=self.hfft.GetXaxis().GetXmax()-self.hfft.GetXaxis().GetXmin()
        rangedata=self.hPhD.GetXaxis().GetXmax()-self.hPhD.GetXaxis().GetXmin()
        self.fftperiod = rangedata/max # estimate of peak to peak distance
        
    def FindPeaks(self):
        self.fft()
        # find the zero peak
        self.hPhD.Fit("gaus","L","",0-self.fftperiod/2,0+self.fftperiod/2)
        fcn=self.hPhD.GetFunction("gaus")
        self.peakWid=fcn.GetParameter(2)/self.hPhD.GetBinWidth(1) # peak width in units of bins
        # for overlapping peaks, setting the search width smaller
        # appears to help finding the peaks
        self.peakWid=self.peakWid*0.75
        print ("peakwid for TSpectrum Search:",self.peakWid)
        self.npeaks=self.ts.Search(self.hPhD,self.peakWid)
        
        print ("TSpectrum found",self.npeaks,"peaks")
        xvals=self.ts.GetPositionX()
        yvals=self.ts.GetPositionY()
        for i in range(self.npeaks):  # store peaks as a list of 2 element lists
            self.xyPeaks.append([xvals[i],yvals[i]])
        self.xyPeaks.sort()
        del self.xyPeaks[8:]   # limit analysis to first 8 peaks (including 0pe)

        # check to see if the 0PE peak was found
        fcn=self.hPhD.GetFunction("gaus")
        npe0mu=fcn.GetParameter(1)
        npe0sig=fcn.GetParameter(2)
        if abs(npe0mu-self.xyPeaks[0][0])>npe0sig:
            print ("Warning: Noise peak is missing, attempting to fix...")
            xaxis=self.hPhD.GetXaxis()
            npe0bin=xaxis.FindBin(npe0mu)
            npe0y=self.hPhD.GetBinContent(npe0bin)
            print (self.xyPeaks)
            self.xyPeaks.insert(0,[npe0mu,npe0y])
            print (self.xyPeaks)
            self.npeaks=self.npeaks+1
        return self.npeaks

    # warning must call FindPeaks first
    def FitPhD(self):
        #print("load status=",r.gSystem.Load("PEFitter_C.so"))
        r.gSystem.CompileMacro("PEFitter.C","k")
        self.npefcn=r.TF1("npefcn",r.fcn,0,10,15)
        self.npefcn.Print()
        
        parnames=["nfit","a0","mu0","sig0","enf","gain"] # other names appended below
        xmin=self.hPhD.GetXaxis().GetXmin()
        xmax=self.hPhD.GetXaxis().GetXmax()
        #print "fcn min max",xmin,xmax
        self.npefcn.SetRange(xmin,xmax)
        self.npefcn.SetNpx(self.hPhD.GetNbinsX())

        npePeaks=len(self.xyPeaks)-1  # remove 0pe peak from count

        fcn0=self.hPhD.GetFunction("gaus")
        mu0=fcn0.GetParameter(1)  # noise peak mean
        sig0=fcn0.GetParameter(2)
        ymax=self.hPhD.GetMaximum()
        gain=self.xyPeaks[2][0]-self.xyPeaks[1][0] # approx gain as dist btwn peak 2 and peak 1

        self.npefcn.FixParameter(0,npePeaks) # of peaks to fit 0=noise only, 1=noise+1pe, ....
        self.npefcn.SetParameter(1,self.xyPeaks[0][1]) # noise peak normalization "a0"
        self.npefcn.SetParLimits(1,0,ymax)
        self.npefcn.SetParameter(2,mu0)
        self.npefcn.SetParameter(3,sig0/gain) # use noise peak width as fraction of gain
        self.npefcn.SetParLimits(2,mu0-2*sig0,mu0+2*sig0)
        self.npefcn.SetParameter(4,sig0/gain) # enf -- starting guess, again as fraction of gain
        self.npefcn.SetParameter(5,gain) # gain approx, dist btwn peak1&0
        
        for i in range(npePeaks):
            self.npefcn.SetParameter(6+i,self.xyPeaks[1+i][1]) # heights of peaks 1...n
            self.npefcn.SetParLimits(6+i,0,ymax)
            parnames.append("a"+str(i+1))
            #print "a"+str(i+1),self.xyPeaks[1+i][1]
        for i in range(len(parnames)): self.npefcn.SetParName(i,parnames[i])
        for i in range(len(parnames),self.npefcn.GetNpar()): self.npefcn.FixParameter(i,0)  # unused parameters
        xend = self.xyPeaks[-1][0]+sig0*1  # end of fit range

        for i in range(len(parnames)):
            print(self.npefcn.GetParName(i),self.npefcn.GetParameter(i))

        c3=r.TCanvas("c3","c3")
        self.hPhD.Draw()
        self.npefcn.Draw("same")
        c3.Update()
        input("Press Enter to continue...")
        c4=r.TCanvas("c4","c4")
        self.npefcn.SetLineWidth(2)
        self.hPhD.Fit("npefcn","","",xmin,xend)
        #self.hPhD.Draw()
        c4.Update()
        input("Press Enter to continue...")
        #self.fcn0=r.TF1("fcn0","gaus",xmin,xmax)
        #self.fcn0.SetParameters(npefcn.GetParameter(1),npefcn.GetParameter(2),
        #                        npefcn.GetParameter(3)*npefcn.GetParameter(5))
        #self.fcn0.SetRange(xmin,xmax)
        #self.fcn0.SetLineStyle(2)
        #self.hPhD.Draw();
        #npefcn.Draw("same")
    

    def GetNoise(self):
        return self.npefcn.GetParameter("sig0")
    def GetENF(self):
        return self.npefcn.GetParameter("enf")
    def GetGain(self):
        return self.npefcn.GetParameter("gain")
    

