function [smoothed_scalo] = temporal_smooth(cwt_x,cwt_y,a,kappa,delt,d)
% Authors: Dr Edward A K Cohen and Prof Andrew T Walden, Department of
% Mathematics, Imperial College London
% Date: 18 Aug 2015
% 
% This function computes the temporally smoothed scalogram (auto or cross)
% from two continuous wavelet transforms. 
% and follows the method as outlined in 
% E.A.K. Cohen and A.T. Walden, "A Statistical Study of Temporally Smoothed
% Wavelet Coherence", IEEE Transactions on Signal Processing, Vol. 58, No.
% 6, pp 2964 - 2973, 2010.
%
% Inputs:       cwt_x           continuous wavelet transform of first
%                               time series
%               cwt_y           continuous wavelet transform of second
%                               time series
%               a               scales on which CWT is evaluated
%               kappa           smoothing parameter
%               delt            sampling frequency
%               d               Morlet wavelet shape parameter
%
%       ***If autoscalogram then cwt_x and cwt_y are the same****
%
% Outputs:      smoothed_scalo       2D array of dimension cwt_x
%                                    with (i,j)th element being the smoothed 
%                                    scalogram evaluated at scale a(i) and time point 
%                                    (j-1)*delt. A value of NaN indicates the 
%                                    smoothed scalogram is undefined for that scale/time pair.

% time series length
N = size(cwt_x,2);
% normalise scales
a_0 = a./delt;

% M for all a   
M = ceil(kappa*abs(a_0));
% nu for all a
nu = ceil(3*abs(a_0)*d);

smoothed_scalo = NaN(size(cwt_x));

for jj = 1:length(a)

    % smallest b value that scalogram can be calculated for that a
    b_start = nu(jj) + M(jj);
    % largest b value that scalogram can be calculated for that a
    b_end = N - 1 - nu(jj) - M(jj);

    for b = b_start:b_end;
        % smoothing width
        tau = [b-M(jj):b+M(jj)];
        % smooth
        smoothed_scalo(jj,b) = mean(cwt_x(jj,tau).*conj(cwt_y(jj,tau)));
    end
end