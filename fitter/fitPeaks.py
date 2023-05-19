import ROOT as r
from peakFitter import *
import sys, os


if len(sys.argv)<2:
    print("No input file given")
    sys.exit()

tf=r.TFile(sys.argv[1])
hLight=tf.Get("phd")
hLight.Rebin()
hLight.Draw()
input("Press Enter to continue...")


ana=PeakFitter(hLight.Clone())
ana.FindPeaks() # this is required before doing the fitting
ana.FitPhD() # do a nice fit to the peaks

print("Results in units of ADC counts")
print("Gain: ",ana.GetGain())
print("Noise: ",ana.GetNoise())
print("ENF: ",ana.GetENF())


