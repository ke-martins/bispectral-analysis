function [ gamma_fr , gamma_fr2 , gamma_am ] = fun_compute_krms_terms( f , P , B , varargin )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the non-linear frequency and amplitude dispersion terms contributing to the 
% root-mean square wavenumker krms following Herbers et al. (2000).
% Compared to their Eq. 12, we retain frequency terms at order mu^2, in 
% order to improve the linear dispersive properties in deeper water depth. This function was written for 
% Boussinesq-based depth inversion applications. See the following reference for more details:
%
%    Martins, K., Bonneton, P., de Viron, O., Turner, I. L., Harley, M. D., & Splinter, K. (2023). 
%    New Perspectives for Nonlinear Depth‐Inversion of the Nearshore Using Boussinesq Theory. 
%    Geophysical Research Letters, 50(2), e2022GL100498.
%
% Inputs:
%   f  - frequency array [Hz]
%   P  - power spectrum [m^2], not a density
%   B  - power bispectrum [m^3], not a density
%   fc - optional cutoff frequency [Hz], in case the raw timeseries was noisy
%
% Outputs: 
%   gamma_fr  - frequency dispersion term, corresponding to beta_fr/h, size of input f
%   gamma_fr2 - fourth order frequency dispersion term, corresponding to beta_fr2/h^2, size of input f
%   gamma_am  - amplitude dispersion term, corresponding to h*beta_am, size of input f
%
% Notes: f, P and B arrays are two-sided in regards to f; f should be centered around 0 Hz
%
% January 23, 2025
% Kévin Martins - kevin.martins@cnrs.fr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % Checking inputs
  if nargin < 4
    fc = max(f);
  elseif nargin == 4
    fc = varargin{1};
  end  
  
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
  
  for ifn = 1:length(f)
    % Initialisation of the cumulative sum and wavenumbers
    sumtmp = 0;
    
    % Loop over frequencies
    for ifp = 1:length(f)
      ifn_p = nmid + (ifn-nmid) - (ifp-nmid);
      if (ifn_p >= 1) && (ifn_p <= length(f))
        if and(abs(f(ifn_p)) > df,abs(f(ifp)) > df)
          sumtmp = sumtmp + df*real(B(ifp,ifn_p));
        end
      end
    end
    
    % Non-linear terms
    gamma_fr(ifn)  = (2*pi*f(ifn))^2/(3*g);
    gamma_fr2(ifn) = (2*pi*f(ifn))^4/(6*g)^2;
    gamma_am(ifn)  = 3*sumtmp/(2*P(ifn));
  end

  return
end
  
