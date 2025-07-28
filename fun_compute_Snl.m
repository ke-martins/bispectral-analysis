function [ Snl ] = fun_compute_Snl( h0 , f , B )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the source term 'Snl' for the non-linear energy transfers between triads following Herbers and Burton (1997).
% As discussed in Martins et al. (2021, https://doi.org/10.1016/j.coastaleng.2021.103917), we take the conjugate of the 
% bispectrum B as computed in this toolbox with fun_compute_bispectrum.m in order to be consistent with their Eq. 14.
%
% Inputs:
%   h0 - local mean water depth [m]
%   f  - frequency array [Hz]
%   B  - power bispectrum [m^3], not a density
%
% Notes: f and B arrays are two-sided in regards to f; f should be centered around 0 Hz
%
% Outputs: 
%   Snl - source term for non-linear energy transfers between triads [m^2], size of input f
%
% Last update on July 28, 2025
% Kévin Martins - kevin.martins@cnrs.fr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % Initialisation
  Snl  = NaN*f;
  nmid = (length(f)-1)/2 + 1; % Middle frequency (f = 0)

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
        sumtmp = sumtmp + (3*pi*f(fi)/h0)*imag(conj(B(ff,ifr3)))*df;
      end
    end
    
    % Source term for non-linear energy transfers
    Snl(fi) = sumtmp;
  end
  return
end
