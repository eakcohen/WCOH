function DoF = calcdof(a,delt,kappa)
% Authors: Dr Edward A K Cohen and Prof Andrew T Walden, Department of
% Mathematics, Imperial College London
% Date: 18 Aug 2015
%
% This function computes the degrees of freedom in the wavelet coherence
% esimator and follows the method as outlined in 
% E.A.K. Cohen and A.T. Walden, "A Statistical Study of Temporally Smoothed
% Wavelet Coherence", IEEE Transactions on Signal Processing, Vol. 58, No.
% 6, pp 2964 - 2973, 2010.
%
% based on Morlet wavelet with shape parameter d=1;
%
% Inputs:   delt        sampling frequency
%           kappa       smoothing parameter
%           a           scale (single value only)    
%
% Outputs:  DoF         Degrees of freedom for wavelet coherence estimator

% use morlet wavelet with shape parameter d
d = 1;
% calc a0 
a0 = a./delt;

% compute smoothing width M (number of points each side used in smoothing)
M = ceil(kappa*a0);

% total number of cwt values used in smoothed scalo
NB = 2*M + 1;

% compute nu
nu = ceil(3*a0*d);

% number of time series points used in each cwt
NS = 2*nu + 1;

a0sigma = nu/3;

% set up time vector (unitless)
tg=[(-(NS-1)/2):((NS-1)/2)];

% compute taper g
g = exp((-1/2)*((tg)/(a0sigma)).^2);

% sum of squared values
C = sum(g.*g);
% normalised
g = g/sqrt(C);

% shift factor
q = 1;

G=zeros(NP,NB);
for j=1:NB
    G([(j-1).*q+1:(j-1).*q+NS],j)= transpose(g);
end
B=G./sqrt(NB);
% outer product matrix
OP=B*transpose(B);

% find largest NB eigenvalues and eigenvectors
options.disp=0; % suppresses output
[u,D,flag]=eigs(OP,NB,'la',options);
if flag ~= 0
    error(' failure of eigenvector routine eigs');
end
% the vectors have norm unity:
%
% now sort out polarity of eigenvectors

% eigenvalues
lambda=diag(D);
weights = lambda*NB;

%%% calculate degrees of freedom
DoF = (sum(weights)^2)/sum(weights.^2);
