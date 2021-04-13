# Gorkha NepalEarthquake Amplification and Signal Processing
This repository contains ground motion statistics on the M7.8 2015 Gorkha, Nepal earthquake. The MATLAB code used to process the raw time series data is included.
The aim of this project was to quantify the amount of ground motion amplification that occurred at two seismic sites in Nepal.

# Time series analysis
The earthquake data was first loaded, then had the linear trend and mean removed from each of the 3 station datasets.
A butterworth filter between 0.05Hz-20Hz was used to filter out the desired frequency content above the nyquist rate.
The data is plotted as a time series so that the arrival times (P,S waves, Love, Rayleigh waves) could be identified on the seismograms.

The peak ground motions (PGA,PGV,PGD) are computed to determine the peak ground motions at the seismic stations.
The data is Fourier transformed to transform the data from the time domain into the frequency domain.
This allows for the power spectrum to be computed, giving information on the frequencies with the greatest amplitudes.

# HVSR

Using the Fourier transformed data, the horizontal-vertical spectral ratio is computed. 
This HVSR gives an indication of the amplification of the surface waves caused by the topography and subsurface soils.
A high HVSR ratio means a large amount of amplification occured in the region, such as in a sedimentary basin.


