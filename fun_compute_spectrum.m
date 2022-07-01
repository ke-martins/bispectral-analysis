function [ data ] = fun_compute_spectrum ( x , fs , nfft , overlap , wind )
% Estimation of bispectrum [m^3] using a direct fft-based approach.
%
%	Inputs:
%	  x       - signal
%	  fs      - sampling frequency
%	  nfft    - fft length [default = 128]
%	  overlap - percentage overlap [default = 50]
%	  wind    - Type of window for tappering ('rectangular', 'hann' or 'kaiser')
%
%	Outputs: 
%   data    - a self-explanatory data structure containing spectra products
%             For more details, see through the code, where the data is stored.
%
% September 10, 2020
% Kévin Martins - kevin.martins@u-bordeaux.fr

  % --------------------- Various parameters -------------------------
  % Notes: 
  %        nfft is forced to be even, this will be useful to get frequencies
  %        centered around 0.

  lx = length(x); x = x(:);
  if (exist('nfft') ~= 1)            nfft = 128; end
  if (exist('overlap') ~= 1)      overlap = 50;  end
  if (isempty(wind) == 1)            wind = 'rectangular'; end
  overlap = min(99,max(overlap,0));

  nfft     = nfft - rem(nfft,2);
  eadvance = fix(nfft * overlap / 100);
  nadvance = nfft - eadvance;
  nblock   = fix((lx - eadvance) / nadvance) + 1; % +1 for not throwing any data out


  % ---------------------- Initialization ------------------------

  if (rem(nfft,2) == 0)
    freqs = [-nfft/2:nfft/2]'/nfft*fs;
  end

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
  locseg = [1:nfft]';       % Indices for first block


  % ---------------------- Compute FFT ----------------------------
  
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

  % Accumulating products (loop over blocks)
  for kk = 1:nblock
    % Block kk FFT
    A_loc  = A(:,kk);
    CA_loc = conj(A(:,kk));
    
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
  data.nblocks      = nblock-1;
  data.nblocks_info = 'Number of blocks used to compute the PSD';

  % Equivalent number of degrees of freedom, based on the estimator given in Welch (1967)
  data.edof      = fun_compute_edof( ww , nfft , lx , overlap ); % NB: ww already contains window information
  data.edof_info = 'Equivalent number of degrees of freedom (Welch, 1967; see also report by Solomon, Jr., 1991)';
  alpha          = 1-0.95;
  data.CI        = fix(data.edof)./chi2inv([1-alpha/2 alpha/2],fix(data.edof));
  data.CI_info   = '95% confidence interval';
end
