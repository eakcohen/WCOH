function [amin,amax] = calcarange(delt,kappa,N)
% Authors: Dr Edward A K Cohen and Prof Andrew T Walden, Department of
% Mathematics, Imperial College London
% Date: 18 Aug 2015
% 
% This function computes the minimum and maximum values of scale that
% the temporally smoothed wavelet coherence can be calculated at
% and follows the method as outlined in 
% E.A.K. Cohen and A.T. Walden, "A Statistical Study of Temporally Smoothed
% Wavelet Coherence", IEEE Transactions on Signal Processing, Vol. 58, No.
% 6, pp 2964 - 2973, 2010.
%
% based on Morlet wavelet with shape parameter d=1;
% Inputs:       delt        sampling frequency
%               kappa       smoothing parameter
%               N           length of time series (number of points)
%
% Outputs:      amin        minimum permissible scale
%               amax        maximum permissible scale

% f'(d) is frequency at which wavelet frequency response is 10e-6 * peak
% value
d = 1;
fd = 1.8;
% Nyquist freq
fn = 1./(2*delt);
% min scale to prevent wavelet aliasing
amin = fd/fn;
% max scale value to prevent wraparound
amax = delt*(N-5)/(6*d + 2*kappa);
