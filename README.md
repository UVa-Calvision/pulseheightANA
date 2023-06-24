# ANA-DRS
Analysis code for data taken with DRS module

```
# Structure of the DRS data TTree
******************************************************************************
*Tree    :waves     : CV Laser Studies DRS                                   *
*Entries :     4558 : Total =       187163819 bytes  File  Size =  108766047 *
*        :          : Tree compression factor =   1.72                       *
******************************************************************************
*Br    0 :time      : vector<float>                                          *
*Entries :     4558 : Total  Size=   18752485 bytes  File Size  =     936459 *
*Baskets :      188 : Basket Size=    1292800 bytes  Compression=  20.02     *
*............................................................................*
*Br    1 :chs       : vector<vector<float> >                                 *
*Entries :     4558 : Total  Size=  168377370 bytes  File Size  =  107805335 *
*Baskets :     1243 : Basket Size=   25600000 bytes  Compression=   1.56     *
*............................................................................*
*Br    2 :event     : event/i                                                *
*Entries :     4558 : Total  Size=      19009 bytes  File Size  =       6846 *
*Baskets :        4 : Basket Size=      51200 bytes  Compression=   2.71     *
*............................................................................*
*Br    3 :min       : min/f                                                  *
*Entries :     4558 : Total  Size=      14433 bytes  File Size  =       5135 *
*Baskets :        4 : Basket Size=      51200 bytes  Compression=   2.72     *
*............................................................................*
```
* time: sample time, roughly in units of 0.2ns
* chs[0...7]: the data for module 1 (1st 8 channels), each of these is a vector<float>[with 1000 samples]
* chs[8]: digitized trigger data, a vector<float>[with 1000 samples]
* event: event number
* min: this is the location of minimal y-value of the trigger pulse.  We use this to align the data



## build the DRSPulseAna executable
This code will open a data fine from the DRS and plot the pulse height spectrum.

mkdir build <br>
cd build <br>
cmake .. <br>
make <br>
run the code with an input file: <br>
./DRSPulseAna A1-C2-L1440-DG-5k.root <br>
This creates an output file: A1-C2-L1440-DG-5k_out.root

Note: Pulse test generates the following plots:
* A plot overlaying all the trigger pulses
* A plot overlaying all the signal channel buffers
* A profile histogram of all the trigger pulses after alignment
* A profile histogram of all the signal buffers after alignment

From the signal histogram we determine a baseline and the peak position. Then the pulse high distribution is shown in two histograms (after baseline subtraction):
* The pulse height distribution of samples at the signal peaking location
* The distribution of pulse integrals approximately from the start of the pulse to before the polarity change (remember we are AC coupled).  The integration range is shown in the 4th plot above.

  
## fitting the peaks
In the fitter directory: <br>
python fitPeaks.py ../build/A1-C2-L1440-DG-5k_out.root
  
