function [ data ] = fun_compute_bispectrum( x , fs , nfft , overlap , wind )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Direct FFT-based estimation of surface elevation bispectrum [m^3].
% Here, the bicoherence is not computed, use e.g. fun_compute_bispectrum_H2001 if needed.
%
% Inputs:
%   x       - Detrended free surface elevation signal [m]
%   fs      - sampling frequency [Hz]
%   nfft    - fft length [default = 256]
%   overlap - percentage overlap (typical is 50%, 75% optimises edof)
%   wind    - Type of window for tappering ('rectangular', 'hann' or 'kaiser')
%
% Outputs: 
%   data    - a self-explanatory data structure containing bispectra products
%             For more details, see through the code, where the data is stored.
%
% Notes:
%   1 - The power bispectrum is computed as:
%                      B(f1,f2) = E[ A(f1) A(f2) A*(f1 + f2) ],
%       where A are the complex Fourier coefficients, A* denotes the conjugate 
%       of A and E is the expected value.
%   2 - Merging over frequencies is no longer optional, as it leads to unrealistic bicoherences (>1)
%
% January 15, 2025
% Kévin Martins - kevin.martins@cnrs.fr
%
% Acknowledgments:
%   Bits of this routine were inspired by codes from the hosa library and those
%   provided by Anouk de Bakker.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % --------------------- Various parameters -------------------------
  % notes: 
  %        nfft is forced to be even, this will be useful to get frequencies
  %        centered around 0.

  lx = length(x); x = x(:);
  if (exist('nfft') ~= 1)           nfft = 256; end
  if (exist('overlap') ~= 1)      overlap = 50; end
  if (isempty(wind) == 1)
    wind = 'rectangular';
  else
    wind = lower(wind);
  end
  overlap = min(99,max(overlap,0));

  nfft     = nfft - rem(nfft,2);
  eadvance = fix(nfft * overlap / 100);
  nadvance = nfft - eadvance;
  nblock   = fix((lx - eadvance) / nadvance) + 1; % +1 for not throwing any data out


  % ---------------------- Initialization ------------------------

  if (rem(nfft,2) == 0)
    freqs = [-nfft/2:nfft/2]'/nfft*fs;
  end

  data.info    = 'Power spectra and bispectra (not densities)';
  data.f       = freqs;
  data.f_info  = 'Frequency [Hz] (two-sided)';
  data.df      = abs(freqs(2)-freqs(1));
  data.df_info = 'Frequency resolution [Hz]';
  data.P       = zeros(nfft+1,1);
  data.P_info  = 'Power spectrum [m^2]';
  data.B       = zeros(nfft+1,nfft+1);
  data.B_info  = 'Power bispectrum [m^3]';


  % ---------------------- Compute FFT ----------------------------
  % Notes: 
  %        In order to be consistent with the computation of the skewness and 
  %        asymmetry by statistical means, use a rectangular window (large differences between
  %        third-order parameters can otherwise be observed).

  % Initialization
  A = zeros(nfft+1,nblock); % Fourier coefficients, for each block
  nmid = (nfft)/2 + 1;      % Middle frequency (f = 0)
  locseg = [1:nfft]';       % Indices for first block

  % Computing FFT (loop over blocks)
  for kk = 1:nblock-1
    % Preparing block kk timeseries
    % Treatment varies with the window
    % For the rectangular window, we force a certain continuity between blocks
    xseg = x(locseg);
    xseg = detrend(xseg);          % Detrend
    xseg = (xseg(:) - mean(xseg)); % De-mean
    switch wind
      case 'hann'
        ww = window(@hann,nfft); normFactor = mean(ww.^2);
        xseg = xseg.*ww / sqrt(normFactor);
      case 'kaiser'
        ww = window(@kaiser,nfft,3.5); normFactor = mean(ww.^2);
        xseg = xseg.*ww / sqrt(normFactor);
      case 'rectangular'
        % Trying to make it periodic
        count = 0;
        while abs(xseg(end)-xseg(1)) > 0.2*nanstd(xseg)
          % Updating locseg
          if kk == 1
            locseg = locseg + 1;
          else
            locseg = locseg - 1;
          end
          % Updating xseg
          xseg = x(locseg);
          count = count + 1;
        end
        if count > 1
          xseg = detrend(xseg);          % Detrend
          xseg = (xseg(:) - mean(xseg)); % De-mean
        end
        %     disp( ['block k, count = ',num2str(count)] )
        
        % Smoothing both the timeseries' head and tail
        beg_ts = xseg(1:2*fs);% length(beg_ts)
        end_ts = xseg(end-2*fs+1:end);% length(end_ts)
        merged_ts0 = [end_ts;beg_ts]; merged_ts = merged_ts0;
        dti = round( fs/8 );
        for tt = dti+1:length(merged_ts)-dti
          merged_ts(tt) = nanmean(merged_ts0(tt-dti:tt+dti));
        end
        xseg(1:2*fs) = merged_ts(end-2*fs+1:end);
        xseg(end-2*fs+1:end) = merged_ts(1:2*fs);
        
        % Final windowing
        ww = window(@rectwin,nfft); normFactor = mean(ww.^2);
        xseg = xseg.*ww / sqrt(normFactor);
      otherwise
        error('Wrong window method')          
    end
    
    % FFT
    A_loc   = fft( xseg , nfft )/nfft;
    A(:,kk) = [ A_loc(nmid:nfft,:) ; A_loc(1:nmid,:) ]; % FFTshift
    A(nmid,kk) = 0;
    
    % Indices for next block
    locseg = locseg + nadvance;
  end
  % Last block, we are not throwing any data out
  for kk = nblock:nblock
    % Preparing block kk timeseries
    % Treatment varies with the window
    % For the rectangular window, we force a certain continuity between blocks
    locseg = [length(x)-nfft+1:length(x)]';
    xseg = x(locseg);
    xseg = detrend(xseg);          % Detrend
    xseg = (xseg(:) - mean(xseg)); % De-mean
    switch wind
      case 'hann'
        ww = window(@hann,nfft); normFactor = mean(ww.^2);
        xseg = xseg.*ww / sqrt(normFactor);
      case 'kaiser'
        ww = window(@kaiser,nfft,3.5); normFactor = mean(ww.^2);
        xseg = xseg.*ww / sqrt(normFactor);
      case 'rectangular'
        % Trying to make it periodic
        count = 0;
        while abs(xseg(end)-xseg(1)) > 0.2*nanstd(xseg)
          % Updating locseg
          locseg = locseg - 1;
          % Updating xseg
          xseg = x(locseg);
          count = count + 1;
        end
        if count > 1
          xseg = detrend(xseg);          % Detrend
          xseg = (xseg(:) - mean(xseg)); % De-mean
        end
        %     disp( ['block k, count = ',num2str(count)] )
        
        % Smoothing both the timeseries' head and tail
        beg_ts = xseg(1:2*fs);% length(beg_ts)
        end_ts = xseg(end-2*fs+1:end);% length(end_ts)
        merged_ts0 = [end_ts;beg_ts]; merged_ts = merged_ts0;
        dti = round( fs/8 );
        for tt = dti+1:length(merged_ts)-dti
          merged_ts(tt) = nanmean(merged_ts0(tt-dti:tt+dti));
        end
        xseg(1:2*fs) = merged_ts(end-2*fs+1:end);
        xseg(end-2*fs+1:end) = merged_ts(1:2*fs);
        
        % Final windowing
        ww = window(@rectwin,nfft); normFactor = mean(ww.^2);
        xseg = xseg.*ww / sqrt(normFactor);
    end
    
    % FFT
    A_loc   = fft( xseg , nfft )/nfft;
    A(:,kk) = [ A_loc(nmid:nfft,:) ; A_loc(1:nmid,:) ]; % FFTshift
    A(nmid,kk) = 0;
    
    % Indices for next block
    locseg = locseg + nadvance;
  end


  % -------------- Computing bispectrum products -------------------

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
  end
  clear ifr1 ifr2 ifr3

  % Expected values
  data.B = data.B / nblock; data.B(~ifm3val) = 0;
  data.P = data.P / nblock;
  
  
  % ------------------- Skewness and asymmetry ---------------------
  % Notes: 
  %        Computation following Elgar and Guza (1985)
  %        Observations of bispectra of shoaling surface gravity wave
  %        Journal of Fluid Mechanics, 161, 425-448

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
  Sk = real( sumtmp ) / nanmean( (x-nanmean(x)).^2 )^(3/2);
  As = imag( sumtmp ) / nanmean( (x-nanmean(x)).^2 )^(3/2);


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
  data.edof      = fun_compute_edof( ww , nfft , lx , overlap ); % NB: ww already contains window information
  data.edof_info = 'Equivalent number of degrees of freedom in the PSD (Welch, 1967; see also report by Solomon, Jr., 1991)';

  % Computing the 95% confidence interval of the PSD
  alpha          = 1-0.95;
  data.P_CI      = fix(data.edof)./chi2inv([1-alpha/2 alpha/2],fix(data.edof));
  data.P_CI_info = '95% confidence interval of the PSD';

  return
end
