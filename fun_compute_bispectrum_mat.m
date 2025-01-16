function [ data ] = fun_compute_bispectrum_mat( x , fs , overlap , wind )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Direct FFT-based estimation of surface elevation bispectrum [m^3].
% The function is written for data organised by blocks in a matrix, 
% which can be handy to accomodate for gappy series (e.g., lidar-derived).
% The first dimension of this matrix corresponds to the nfft, i.e. the length 
% of our timeseries block. Overlapping (if any) thus has already been taken care of,
% and the input 'overlap' is only used to compute the corresponding equivalent dof.
%
% Inputs:
%   x       - Detrended free surface elevation signal, organised in a matrix of dimensions nfft x nb_block
%   fs      - sampling frequency [Hz]
%   overlap - percentage overlap, used only to compute edof
%   wind    - Type of window for tappering ('rectangular', 'hann' or 'kaiser')
%
% Outputs: 
%   data   - a self-explanatory data structure containing bispectra products
%            For more details, see through the code, where the data is stored.
%
% Notes:
%   1 - The power bispectrum is computed as:
%                      B(f1,f2) = E[ A(f1) A(f2) A*(f1 + f2) ],
%       where A are the complex Fourier coefficients, A* denotes the conjugate 
%       of A and E is the expected value.
%   2 - Compared to the other bispectrum functions, this one contains a few more
%       bispectral products, such as the correlation coefficient as defined by Hinich et al. (1965)
%       or the standard deviation of bispectral estimates. These are useful to quantify statistical 
%       uncertainties on bispectrum-derived quantities.
%  
% Last update on January 15, 2025
% Kévin Martins - kevin.martins@cnrs.fr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % --------------------- Various parameters -------------------------

  if (isempty(wind) == 1)
    wind = 'rectangular';
  else
    wind = lower(wind);
  end
  if (exist('overlap') ~= 1)      overlap = 50; end
  nfft   = size(x,1);
  nblock = size(x,2);


  % ---------------------- Initialization ------------------------
  % notes: nfft is forced to be even, this will be useful to get frequencies
  %        centered around 0.

  if (rem(nfft,2) ~= 0)
    x    = x(1:end-1,:);
    nfft = nfft-1;
  end
  freqs = [-nfft/2:nfft/2]'/nfft*fs;

  data.info     = 'Power spectra and bispectra (not densities)';
  data.f        = freqs;
  data.f_info   = 'Frequency [Hz] (two-sided)';
  data.df       = abs(freqs(2)-freqs(1));
  data.df_info  = 'Frequency resolution [Hz]';
  data.P        = zeros(nfft+1,1);
  data.P_info   = 'Power spectrum [m^2]';
  data.B        = zeros(nfft+1,nfft+1);
  data.B_info   = 'Power bispectrum [m^3]';
  data.B_std    = zeros(nfft+1,nfft+1);
  data.B_std_info = 'Bispectrum [m^3] standard deviation, computed after Hinich and Clay (1968); Kim and Powers (1979)';
  data.Bic      = zeros(nfft+1,nfft+1);
  data.Bic_info = 'Bicoherence [-] defined after Hagihira et al. (2001)';
  data.b        = zeros(nfft+1,nfft+1);
  data.b_info   = 'Correlation coefficient [-]';


  % ---------------------- Compute FFT ----------------------------

  % Initialization
  A = zeros(nfft+1,nblock); % Fourier coefficients, for each block
  nmid = (nfft)/2 + 1;      % Middle frequency (f = 0)
  
  % Computing FFT (loop over blocks)
  for kk = 1:nblock
    % Preparing block (window type)
    xseg = x(:,kk);
    xseg = detrend(xseg(:) - mean(xseg)); % De-mean
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


  % -------------- Computing bispectrum products -------------------
  
  % Initialization
  P12  = zeros(nfft+1,nfft+1); % Double product (unit m^3), needed to compute correlation coefficient (Hinich et al., 1965)
  P3   = zeros(nfft+1,nfft+1); % Power at coupling frequency (unit m^3), needed to compute correlation coefficient (Hinich et al., 1965)
  P123 = zeros(nfft+1,nfft+1); % Triple product (unit m^3), see Appendix of Hagihira et al. (2001)

  % Dealing with f1 + f2 = f3 indices
  [ ifr1 , ifr2 ] = meshgrid( 1:nfft+1 , 1:nfft+1 );
  ifr3 = nmid + (ifr1-nmid) + (ifr2-nmid);
  ifm3val = ((ifr3 >= 1) & (ifr3 <= nfft+1)); ifr3(~ifm3val) = 1;

  % Accumulating triple products (loop over blocks)
  for kk = 1:nblock
    % Block kk FFT
    A_loc  = A(:,kk);
    CA_loc = conj(A(:,kk));
    
    % Compute bispectrum and PSD
    data.B = data.B + A_loc(ifr1) .* A_loc(ifr2) .* (CA_loc(ifr3));
    data.P = data.P + abs(A_loc.^2);
    P12    = P12  + abs(A_loc(ifr1).*abs(A_loc(ifr2))).^2;
    P3     = P3   + abs(A_loc(ifr3)).^2;
    P123   = P123 + sqrt(abs(A_loc(ifr1)).^2 .* abs(A_loc(ifr2)).^2 .* abs(A_loc(ifr3)).^2);
  end

  % Expected values
  data.B = data.B / nblock; data.B(~ifm3val) = 0;
  data.P = data.P / nblock;
  P12  = P12 / nblock;
  P3   = P3 / nblock;
  P123 = P123 / nblock; P123(~ifm3val) = 0;

  % Storing the bicoherence, the correlation coefficient and the standard deviation of bispectral estimates
  data.Bic   = abs(data.B) ./ P123;
  data.b     = abs(data.B) ./ sqrt(P12.*P3);
  data.B_std = 2*sqrt(data.P(ifr1).*data.P(ifr2).*data.P(ifr3).*(1-data.b.^2)/nblock);
  clear ifr1 ifr2 ifr3

  
  % ------------------- Skewness and asymmetry ---------------------
  % Notes: 
  %        Computation following Elgar and Guza (1985):
  %        Observations of bispectra of shoaling surface gravity wave. JFM 161, 425-448.
  %        In order to be consistent with the computation of the skewness and 
  %        asymmetry by statistical means, use a rectangular window (relatively large differences 
  %        between third-order parameters can otherwise be observed).

  % We work only in one frequency quadrant
  ifreq  = [nmid:nmid+(nfft)/2]';
  sumtmp = 6*data.B(nmid,nmid); % Initialisation with first diagonal term

  % Loop over frequencies 
  for ff = 2:length(ifreq)
    % Frequency index
    indf = ifreq(ff);
    
    % Diagonal
    sumtmp = sumtmp + 6*data.B(indf,indf);
    
    % Non-diagonal
    sumtmp = sumtmp + 12*sum(data.B(indf,nmid:indf-1));
  end
  Sk = real( sumtmp ) / ( 2*trapz(data.f(nmid:end),data.P(nmid:end)/data.df) )^(3/2);
  As = imag( sumtmp ) / ( 2*trapz(data.f(nmid:end),data.P(nmid:end)/data.df) )^(3/2);


  % ------------------ Saving data in structure --------------------
  data.Sk        = Sk;
  data.Sk_info   = 'Wave skewness [-], computed following Elgar and Guza (1985, JFM)';
  data.As        = As;
  data.As_info   = 'Wave asymmetry [-], computed following Elgar and Guza (1985, JFM)';


  % ---------------------- Finalisation ----------------------------

  % Number of blocks used to compute PSD
  data.nblocks      = nblock;
  data.nblocks_info = 'Number of blocks used to compute the PSD';

  % Equivalent number of degrees of freedom, based on the approximation of Percival and Walden (1993)
  % NB: length of original timeseries might not be exact in the case data is not thrown away, but associated error is negligible
  data.edof      = fun_compute_edof( ww , nfft , nblock*nfft - (nblock-1)*nfft*overlap/100 , overlap );
  data.edof_info = 'Effective number of degrees of freedom (Percival and Walden, 1993; their Eq. 292b)';

  % Computing the 95% confidence interval of the PSD
  alpha          = 1-0.95;
  data.P_CI      = fix(data.edof)./chi2inv([1-alpha/2 alpha/2],fix(data.edof));
  data.P_CI_info = '95% confidence interval of the PSD';

  % b95 (Haubrich, 1965)
  data.b95      = sqrt(6/data.edof);
  data.b95_info = '95% significance level on zero bicoherence: sqrt(6/edof) (Haubrich, 1965)';

  return
end
