function [ krms ] = fun_compute_krms( h0 , f , P , B )
% Compute root-mean square wavenumker k following Herbers et al. (2000), Eq. 12
%
% Inputs:
%   h0 - local water depth
%   f  - frequency array 
%   P  - power spectrum [m^2], not a density
%   B  - power bispectrum [m^3], not a density
%
% Notes: f, P and B arrays are two-sided in regards to f; f should be centered around 0 Hz
%
% Outputs: 
%   krms - root-mean square wavenumker k [1/m], size of input f
%
% March 31, 2020
% Kévin Martins - kevin.martins@u-bordeaux.fr
  
  % Initialisation
  krms = NaN*f;
  nmid  = (length(f)-1)/2 + 1; % Middle frequency (f = 0)
  g     = 9.81;                % Gravity
  
  % Transforming P and B in densities
  df = abs(f(2)-f(1));
  P  = P/df;
  B  = B/df^2;
  
  for fi = 1:length(f)
    % Initialisation of the cumulative sum and wavenumbers
    sumtmp         = 0; 
    krms(fi)       = 2*pi*f(fi)/sqrt(g*h0); % Linear part
    
    % Loop over frequencies
    for ff = 1:length(f)
      ifr3 = nmid + (fi-nmid) - (ff-nmid);
      if (ifr3 >= 1) && (ifr3 <= length(f))
        sumtmp = sumtmp + df*real(B(ff,ifr3));
      end
    end
    
    % Non-linear terms
    Beta_fr = h0*(2*pi*f(fi))^2/(3*g);
    Beta_am = 3*sumtmp/(2*h0*P(fi));
    
    % Non-linear wavenumber
    krms(fi) = krms(fi)*sqrt( 1 + Beta_fr - Beta_am);
  end
end
  
