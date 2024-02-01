The code 'draw_IMF_solar_wind.py' in this documentary is used to plot the solar wind parameters, IMF components and geomagnetic indices displayed in Figure 1.

The solar wind parameters and IMF components along with the SYM-H index are restored in 'OMNI_HRO2_1MIN_2074319.csv' which is downloaded from https://cdaweb.gsfc.nasa.gov. 

The other geomagnetic indices are restored in DST.txt and Kp_ap.txt, the raw data is downloaded https://kp.gfz-potsdam.de/en/data and https://wdc.kugi.kyoto-u.ac.jp/dstae/index.html.

The code also calculates the inclination angle of the ICME flux rope in YOZ plane, using the Mean Variance Analysis (MVA) method. The source code of this method is written by *** Harry Wheeler *** and is from https://github.com/harrywheeleriv/minvarly/tree/master.
