function [ v ] = fun_compute_edof( w , M , N , overlap )
% Inputs:
%   w       - window (windowed timeseries of FFT length)
%   M       - FFT length 
%   N       - total number of points
%   overlap - overlap in % used to compute number of new points in each FFT
%
% Outputs:
%   v       - equivalent degrees of freedom in a PSD estimate
%
% Notes:
%   This uses the recommendations from the report "PSD Computations Using Welch’s Method" 
%   written by Otis M. Solomon, Jr., in 1991 (Sandia Report, SAND91-1533 UC-706).
%   It uses the concept of equivalent degrees of freedom introduced in the original paper by Peter D. Welch:
%   "The use of fast Fourier transform for the estimation of power spectra: A method based on time averaging
%   over short, modified periodograms",  IEEE Transactions on Audio and Electroacoustics, 15 (2), 1967. 
%
% September 10, 2020
% Kévin Martins - kevin.martins@u-bordeaux.fr

  % Initialisation
  s = 0.0;
  
  % Number of new points in each FFT
  S = fix(M*(100-overlap)/100);
  k = 1 + (N - M)/S;
  
  % Loop over lag of window correlation
  for i = 1:k-1
    if (i*S < M-1)
      s = s + (k-i)/k * fun_compute_rho( w , M , i , overlap );
    end
  end
  
  % Equivalent degrees of freedom
  v = 2.0*k / (1.0 + 2.0*s);
end

%% fun_compute_rho computes the window autocorrelation function 
function [ r ] = fun_compute_rho( w , M , k , overlap )
% Inputs:
%   w       - window
%   M       - FFT length 
%   k       - lag  of  window  correlation 
%   overlap - overlap in % used to compute number of new points in each FFT 
%
% Outputs:
%   r       - autocorrelation of window

  % Initialisation
  r = 0.0; 
  Pw = 0.0;
  
  % Number of new points in each FFT
  S = fix(M*(100-overlap)/100);
  
  % Calculation of window autocorrelation
  for i = 1:M
    Pw = Pw + w(i)*w(i)/M;
  end
  for i =1:M-k*S-1
    r = r + w(i)*w(i+k*S);
  end
  r = sqrt(r/(M*Pw));
end

