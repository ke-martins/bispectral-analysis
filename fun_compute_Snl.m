function [ Snl ] = fun_compute_Snl( h0 , f , B )
% Compute Snl following Herbers and Burton (1997), see Eq. 14
%
% Inputs:
%   h0 - local water depth
%   f  - frequency array 
%   B  - power bispectrum [m^3], not a density
%
% Notes: f, P and B arrays are two-sided in regards to f; f should be centered around 0 Hz
%
% Outputs: 
%   Snl - source term for non-linear energy transfers between triads [m^2], size of input f
%
% March 31, 2020
% Kévin Martins - kevin.martins@u-bordeaux.fr
  
  % Initialisation
  Snl  = NaN*f;
  nmid = (length(f)-1)/2 + 1; % Middle frequency (f = 0)
  g    = 9.81;                % Gravity
  
  % Transforming B in densities
  df = abs(f(2)-f(1));
  B  = B/df^2;
  
  for fi = 1:length(f)
    % Initialisation of the cumulative sum and wavenumbers
    sumtmp = 0;
    
    % Loop over frequencies
    for ff = 1:length(f)
      ifr3 = nmid + (fi-nmid) - (ff-nmid);
      if (ifr3 >= 1) && (ifr3 <= length(f))
        sumtmp = sumtmp + W(f(fi),h0)*imag(B(ff,ifr3))*df;
      end
    end
    
    % Non-linear Source term
    Snl(fi) = sumtmp;
  end
end

function res = W(f,h)
  res = (3 * pi * f) / h;
end
  
