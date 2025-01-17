function [ data ] = fun_compute_spectrum_mat( x , fs , overlap , wind )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Direct FFT-based estimation of surface elevation spectral densities [m^2/Hz].
% The function is written for data organised by blocks in a matrix, 
% which can be handy to accomodate for gappy series (e.g., lidar-derived).
% The first dimension of this matrix corresponds to the nfft, i.e. the length 
% of our timeseries block. Overlapping (if any) thus has already been taken care of,
% and the input 'overlap' is only used to compute the corresponding equivalent dof.
%
% Inputs:
%   x       - Detrended free surface elevation signal [m], organised in a matrix of dimensions nfft x nb_block
%   fs      - sampling frequency [Hz]
%   overlap - percentage overlap, used only to compute edof
%   wind    - Type of window for tappering ('rectangular', 'hann' or 'kaiser')
%
% Outputs: 
%   data    - a self-explanatory data structure containing spectra products
%             For more details, see through the code, where the data is stored.
%
% Last updated in January 2024
% Kévin Martins - kevin.martins@cnrs.fr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % --------------------- Various parameters -------------------------

  if nargin < 3, error('check inputs'); end
  if or(isempty(wind) == 1,nargin < 4)
    wind = 'rectangular';
  else
    wind = lower(wind);
  end
  nfft   = size(x,1);
  nblock = size(x,2);


  % ---------------------- Initialization ------------------------

  if (rem(nfft,2) ~= 0)
    x    = x(1:end-1,:);
    nfft = nfft-1;
  end
  freqs = [-nfft/2:nfft/2]'/nfft*fs;

  % Output fields
  data.info    = 'Energy density spectra';
  data.f       = freqs;
  data.f_info  = 'Frequency [Hz] (one-sided)';
  data.df      = abs(freqs(2)-freqs(1));
  data.df_info = 'Frequency resolution [Hz]';
  data.E       = zeros(nfft+1,1);
  data.E_info  = 'Energy density spectrum [m^2/Hz]';
    
  % Local stuff
  A = zeros(nfft+1,nblock); % Fourier coefficients, for each block
  nmid = (nfft)/2 + 1;      % Middle frequency (f = 0)


  % ---------------------- Compute FFT ----------------------------
  
  % Computing FFT (loop over blocks)
  for kk = 1:nblock
    % Preparing block (window type)
    xseg = x(:,kk);
    xseg = (xseg(:) - mean(xseg)); % De-mean
    switch wind
      case 'hann'
        ww = window(@hann,nfft); normFactor = mean(ww.^2);
        xseg = xseg.*ww / sqrt(normFactor);
      case 'kaiser'
        ww = window(@kaiser,nfft,3.5); normFactor = mean(ww.^2);
        xseg = xseg.*ww / sqrt(normFactor);
      case 'rectangular'
        ww = window(@rectwin,nfft); normFactor = mean(ww.^2);
        xseg = xseg.*ww / sqrt(normFactor);
    end
    
    % FFT
    A_loc   = fft( xseg , nfft )/nfft;
    A(:,kk) = [ A_loc(nmid:nfft,:) ; A_loc(1:nmid,:) ]; % FFTshift
    A(nmid,kk) = 0;
  end

  % Accumulating products (loop over blocks)
  for kk = 1:nblock
    % Block kk FFT
    A_loc  = A(:,kk);
    
    % Compute PSD
    data.E = data.E + abs(A_loc.^2);
  end

  % Expected values
  data.E = data.E / nblock;
  
  
  % ---------------------- Finalisation ----------------------------

  % Keeping one-sided spectrum
  data.f = data.f(nmid:end);
  data.E = 2*data.E(nmid:end)/data.df;
  
  % Number of blocks used to compute PSD
  data.nblocks      = nblock;
  data.nblocks_info = 'Number of blocks used to compute the PSD';

  % Equivalent number of degrees of freedom, based on the estimator given in Welch (1967), following Percival and Walden
  % NB: length of original timeseries might not be exact in the case data is not thrown away, but associated error is negligible
  data.edof      = fun_compute_edof( ww , nfft , nfft + fix((nblock-1)*nfft*(100-overlap)/100) , overlap );
  data.edof_info = 'Effective number of degrees of freedom';
  alpha          = 1-0.95;
  data.CI        = fix(data.edof)./chi2inv([1-alpha/2 alpha/2],fix(data.edof));
  data.CI_info   = '95% confidence interval';

  return
end
