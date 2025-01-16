function [ krms ] = fun_compute_krms( h0 , f , P , B , varargin )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the root-mean square wavenumker krms following Herbers et al. (2000, Eq. 12), 
% based on the Boussinesq theory of Herbers and Burton (1997). In this work, the authors
% neglect fourth-order frequency terms, which can improve the linear predictions in deeper 
% water depths. This is an option here.
% 
% Inputs:
%   h0 - local water depth [m]
%   f  - frequency array [Hz]
%   P  - power spectrum [m^2], not a density
%   B  - power bispectrum [m^3], not a density
%   option - 'second' or 'fourth' order frequency dispersion term
%
% Outputs: 
%   krms - root-mean square wavenumker k [1/m], size of input f
%
% Notes: f, P and B arrays are two-sided with regards to f; f should be centered around 0 Hz
%
% Last update on January 15, 2025
% Kévin Martins - kevin.martins@cnrs.fr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % Checking inputs
  if nargin < 5
    order = 'second';
  elseif nargin == 5
    order = char(varargin{1});
    if and(~strcmp(order,'second'),~strcmp(order,'fourth')), error('Wrong input for order of frequency dispersion effects.'); end
  end  
  
  % Initialisation
  nmid = (length(f)-1)/2 + 1; % Middle frequency (f = 0)
  g    = 9.81;                % Gravity
  krms = 2*pi*f/sqrt(g*h0);   % Linear, shallow water wavenumber
  
  % Transforming P and B in densities
  df = abs(f(2)-f(1));
  P  = P/df;
  B  = B/df^2;
  
  % Looping over frequencies
  for fi = 1:length(f)
    % Initialisation of the cumulative sum and wavenumbers
    sumtmp   = 0;
    
    % Loop over frequencies
    for ff = 1:length(f)
      ifr3 = nmid + (fi-nmid) - (ff-nmid);
      if (ifr3 >= 1) && (ifr3 <= length(f))
        sumtmp = sumtmp + df*real(B(ff,ifr3));
      end
    end
    
    % Non-linear terms
    Beta_fr  = h0*(2*pi*f(fi))^2/(3*g);
    Beta_fr2 = h0^2*(2*pi*f(fi))^4/(6*g)^2;
    Beta_am  = 3*sumtmp/(2*h0*P(fi));
    
    % Non-linear wavenumber
    if strcmp(order,'second')
      krms(fi) = krms(fi)*sqrt( 1 + Beta_fr - Beta_am);
    else
      krms(fi) = krms(fi)*sqrt( 1 + Beta_fr + Beta_fr2 - Beta_am);
    end
  end

  return
end
  
