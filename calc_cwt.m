function [cwt_x] = calc_cwt(x,delt,a,d)
% Authors: Dr Edward A K Cohen and Prof Andrew T Walden, Department of
% Mathematics, Imperial College London
% Date: 18 Aug 2015
% 
% This function computes the continuous wavelet transform using a Morlet
% wavelet on the time series x at scales given in vector a
% and follows the method as outlined in 
% E.A.K. Cohen and A.T. Walden, "A Statistical Study of Temporally Smoothed
% Wavelet Coherence", IEEE Transactions on Signal Processing, Vol. 58, No.
% 6, pp 2964 - 2973, 2010.
%
% Inputs:       x           time series data
%               delt        sampling frequency
%               a           scales on which CWT is to be evaluated
%               d           Morlet wavelet shape parameter
%
% Outputs:      cwt_x       2D array of dimension length(a) by length(x)
%                           with (i,j)th element being the CWT of x
%                           evaluated at scale a(i) and time point
%                           (j-1)*delt. A value of NaN indicates the cwt is
%                           undefined for that scale/time pair.

% calculate length of time series
N = length(x);

% Fourier transfrom
x_hat = (fft(x));

% Calculate nyquist frequency
Fn = 1/delt;

% calc frequency vector for fft data
l = [0:N-1];
f = l./N;

% calc scale coefficient a_{0}
a_0 = a./delt;

% How many scales
nscales = length(a);

% calc f_{0} vector with which to calculate wavelet transform frequency responses
% corresponding to scale coefficient a_{0}
f_pos = f(1:N/2);
f_neg = (f(N/2+1:end)-1);

f_0 = [f_pos f_neg];
f_0= repmat(a_0',1,N).*repmat(f_0,nscales,1);

% Initiate the wavelet transform matrices
cwt_x = zeros(nscales,N);


    % Calculate frequency response of wavelet at frequencies in f_{0}
    Psi = (pi^(1/4))*((2*d)^(1/2))*exp(-2*(d*pi*(f_0-1)).^2);
    
    % Calculate the products
    x_hat_Psi     = repmat(x_hat,length(a_0),1).*Psi;

    % For each scale
    for jj = 1:length(a_0);

        % calculate the cwts for each type of wavelet
        cwt_x(jj,:)       = (sqrt(abs(a_0(jj)))/(N*sqrt(delt))) * ifft(x_hat_Psi(jj,:));

    end