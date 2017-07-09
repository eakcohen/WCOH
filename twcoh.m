function [coherence] = twcoh(X,Y,a,delt,kappa)
% Authors: Dr Edward A K Cohen and Prof Andrew T Walden, Department of
% Mathematics, Imperial College London
% Date: 18 Aug 2015
% 
% This function computes the temporally smoothed wavelet coherence for a
% pair of time series and follows the method as outlined in 
% E.A.K. Cohen and A.T. Walden, "A Statistical Study of Temporally Smoothed
% Wavelet Coherence", IEEE Transactions on Signal Processing, Vol. 58, No.
% 6, pp 2964 - 2973, 2010.
% Morlet wavelet used is with shape parameter d=1
% 
% Inputs:   X           Time Series 1
%           Y           Time Series 2
%           a           vector of scales to evaluate wavelet coherence
%           delt        Sampling interval
%           kappa       Smoothing Parameter
%
% Outputs:  coherence   2D array of dimension length(a) by length(X)
%                       with (i,j)th entry corresponding to the wavelet coherence 
%                       at scale a(i) and time (j-1)*delt. A value of NaN
%                       indicates the coherence is undefined for that
%                       scale/time pair.

if nargin < 3
    error('Not enough inputs - two time series and scales required')
elseif nargin == 3
    %default sampling interval set
    delt = 1;
end
% set wavelet shape parameter d
d = 1;
% determine number of data points
N = length(X);
NY = length(Y);
% error message if different lengths
if N ~= NY
    error('Time series X and Y are different lengths. Abort')
end

% check min and max a values
[amin,amax] = calcarange(delt,kappa,N)
if min(a) < amin || max(a) > amax 
    error('at least one value of a falls outside minimum or maximum allowed values')
end

% perform continuous wavelet tranform on each signal
[cwt_X] = calc_cwt(X,delt,a,d);
[cwt_Y] = calc_cwt(Y,delt,a,d);

% compute temporally smoothed scalograms
[smoothed_scaloXX] = temporal_smooth(cwt_X,cwt_X,a,kappa,delt,d);
[smoothed_scaloYY] = temporal_smooth(cwt_Y,cwt_Y,a,kappa,delt,d);
[smoothed_scaloXY] = temporal_smooth(cwt_X,cwt_Y,a,kappa,delt,d);
        
% compute coherence
coherence = (abs(smoothed_scaloXY).^2)./(smoothed_scaloXX.*smoothed_scaloYY);

