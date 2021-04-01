# Gorkha-Nepal-Earthquake-Amplification-Signal-Processing
This repository contains ground motion statistics on the M7.8 2015 Gorkha, Nepal earthquake. The MATLAB code used to process the raw time series data is included.
The aim of this project was to quantify the amount of ground motion amplification that occurred at two seismic sites in Nepal.
# time series analysis
The earthquake data was first loaded, then had the linear trend and mean removed from each of the 3 station datasets.
A butterworth filter between 0.05Hz-20Hz was used to filter out the desired frequency content above the nyquist rate.
The data is plotted as a time series so that the arrival times (P,S waves, Love, Rayleigh waves) could be identified on the seismograms.

