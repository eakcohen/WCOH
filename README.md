# WCOH
Wavelet Coherence

twcoh.m is the main piece of code and calculates the temporally smoothed wavelet coherence for a pair of time series for a given vector of scales at all possible times.
 
calc_cwt.m is called by twcoh.m to compute the cwts
 
temporal_smooth.m is called by twcoh.m to compute the smoothed auto and cross scalograms
 
calcarange.m gives the minimum and maximum allowed scales due to aliasing and wrap-around. This is called in twcoh by can also be used by the user to know the possible range of scale values before inputting the scale vector into twcoh
 
calcdof.m computes the degrees of freedom associated with the wavelet coherence estimator at a particular scale
