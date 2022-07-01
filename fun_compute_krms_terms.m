function [ gamma_fr , gamma_fr2 , gamma_am ] = fun_compute_krms_terms( f , P , B )
% Compute non-linear frequency and amplitude dispersion terms contributing to the 
% root-mean square wavenumker krms following Herbers et al. (2002).
% Compared to his Eq. 12, we retain frequency terms at order mu^2, in 
% order to improve the linear dispersive properties of the Boussinesq approximation.
%
% Inputs:
%   f  - frequency array [Hz]
%   P  - power spectrum [m^2], not a density
%   B  - power bispectrum [m^3], not a density
%
% Notes: f, P and B arrays are two-sided in regards to f; f should be centered around 0 Hz
%
% Outputs: 
%   gamma_fr - frequency dispersion term, corresponding to beta_fr/h, size of input f
%   gamma_fr2 - frequency dispersion term, corresponding to beta_fr2/h^2, size of input f
%   gamma_am - amplitude dispersion term, corresponding to h*beta_am, size of input f
%
% 1 July 2022
% Kévin Martins - kevin.martins@u-bordeaux.fr
  
  % Initialisation
  gamma_fr  = NaN*f;
  gamma_fr2 = NaN*f;
  gamma_am  = NaN*f;
  nmid      = (length(f)-1)/2 + 1; % Middle frequency (f = 0)
  g         = 9.806;               % Gravity
  
  % Transforming P and B in densities
  df = abs(f(2)-f(1));
  P  = P/df;
  B  = B/df^2;
  
  for fi = 1:length(f)
    % Initialisation of the cumulative sum and wavenumbers
    sumtmp = 0;
    
    % Loop over frequencies
    for ff = 1:length(f)
      ifr3 = nmid + (fi-nmid) - (ff-nmid);
      if (ifr3 >= 1) && (ifr3 <= length(f))
        sumtmp = sumtmp + df*real(B(ff,ifr3));
      end
    end
    
    % Non-linear terms
    gamma_fr(fi)  = (2*pi*f(fi))^2/(3*g);
    gamma_fr2(fi) = (2*pi*f(fi))^4/(6*g)^2;
    gamma_am(fi)  = 3*sumtmp/(2*P(fi));
  end
end
  
