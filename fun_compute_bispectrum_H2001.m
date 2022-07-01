function [ data ] = fun_compute_bispectrum_H2001( x , fs , nfft , overlap , wind , mg )
% Estimation of bispectrum products using a direct fft-based approach.
%
%	Inputs:
%	  x       - signal
%	  fs      - sampling frequency
%	  nfft    - fft length [default = 128]
%	  overlap - percentage overlap [default = 50]
%	  wind    - Type of window for tappering ('rectangular', 'hann' or 'kaiser')
%	  mg      - length of spectral bandwith over which to merge
%
%	Outputs: 
%   data    - a self-explanatory data structure containing bispectra products
%             For more details, see through the code, where the data is stored.
%
% Notes:
%   1 - The power bispectrum is computed as:
%                      B(f1,f2) = E[ A(f1) A(f2) A*(f1 + f2) ],
%       where A are the complex Fourier coefficients, A* denotes the conjugate 
%       of A and E is the expected value.
%   2 - The normalisation used here is taken from Hagihira et al. (2001):
%          b(f1,f2) = |B(f1,f2)| / E[ sqrt(|A(f1)|^2 A(f2)|^2 |A(f1 + f2)|^2) ]
%       This is the most correct normalisation as it equals one for fully 
%       coupled triads (all in phase). Performing frequency merging on this 
%       version is weird in my opinion considering the expected value is 
%       taken at last on the sqrt.
%   3 - Merging over frequencies is optional; if desired, enter a frequency 
%       bandwidth in mg. Otherwise, leave it as []. Bicoherence and biphase are
%       recomputed from the merged arrays.       
%
% September 10, 2020
% Kévin Martins - kevin.martins@u-bordeaux.fr
%
% Acknowledgments:
%   Bits of this routine were inspired by codes from the hosa library and those
%   provided by Anouk de Bakker.

  % --------------------- Various parameters -------------------------
  % notes: 
  %        nfft is forced to be even, this will be useful to get frequencies
  %        centered around 0.

  lx = length(x); x = x(:);
  if (exist('nfft') ~= 1)            nfft = 128; end
  if (exist('overlap') ~= 1)      overlap = 50;  end
  if (isempty(wind) == 1)            wind = 'rectangular'; end
  if (isempty(mg) == 1)
    merge = 0;
  else
    merge = 1;
    mg = mg - rem(mg+1,2); 
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

  % Output fields
  data.info    = 'Bispectra products - the bicoherence is computed with Hagihira et al. (2001) normalisation';
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

  % Initialization
  P123 = zeros(nfft+1,nfft+1); % Triple product of the PSD (unit m^3), see Appendix of Hagihira et al. (2001)

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
    P123 = P123 + sqrt(abs(A_loc(ifr1)).^2 .* abs(A_loc(ifr2)).^2 .* abs(A_loc(ifr3)).^2);
  end

  % Expected values
  data.B = data.B / nblock; data.B(~ifm3val) = 0;
  data.P = data.P / nblock;
  P123 = P123 / nblock;     P123(~ifm3val) = 0;

  % Computing bicoherence and biphase
  if ~merge
    data.Bic = abs(data.B) ./ P123;
    data.Bic_info  = 'Bicoherence [-]';
    data.Bip = atan2(imag(data.B),real(data.B));
    data.Bip_info  = 'Biphase [-]';
  end
  clear ifr1 ifr2 ifr3


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


  % ------------------- Merging over frequencies -------------------

  if merge
    % Initialization
    mg    = mg - rem(mg+1,2);
    mm    = (mg-1)/2;                                              % Half-window for averaging
    nmid  = nfft/2 + 1;                                            % Middle frequency (f = 0)
    ifrm  = [ fliplr(nmid:-mg:1+mm) , nmid+mg : mg : nfft+1-mm ]'; % Frequency indices
    Bm    = zeros( length(ifrm) , length(ifrm) );                  % Merged bispectrum (unit m^3)
    Pm    = zeros( length(ifrm) , 1 );                             % Merged PSD (unit m^2)
    P123m = zeros(length(ifrm),length(ifrm));                      % Merged triple product of the PSDs (unit m^3)
    
    % Loop over frequencies
    for jfr1 = 1:length(ifrm)
      % Rows
      ifb1 = ifrm(jfr1); % mid of jfr1-merge-block
      
      % PSD
      Pm(jfr1) = sum(data.P(ifb1-mm:ifb1+mm));
      
      % Columns for bispectrum
      for jfr2 = 1:length(ifrm)
        ifb2 = ifrm(jfr2); % mid of jfr2-merge-block
        Bm(jfr1,jfr2)    = sum(sum(data.B(ifb1-mm:ifb1+mm,ifb2-mm:ifb2+mm)));
        P123m(jfr1,jfr2) = sum(sum(P123(ifb1-mm:ifb1+mm,ifb2-mm:ifb2+mm)));
      end
    end
    
    % Updating arrays
    data.f = freqs(ifrm);
    data.df = abs(data.f(2)-data.f(1));
    data.B = Bm;
    data.P = Pm;
    
    % Computing bicoherence
    data.Bic      = abs(data.B) ./ P123m;
    data.Bic_info = 'Bicoherence [-]';
    clear ifr1 ifr2 ifr3
   
    % Computing biphase
    data.Bip = atan2(imag(data.B),real(data.B));
    data.Bip_info  = 'Biphase [-]';
  end


  % ---------------------- Finalisation ----------------------------

  % Number of blocks used to compute PSD
  data.nblocks      = nblock-1;
  data.nblocks_info = 'Number of blocks used to compute the PSD';

  % Equivalent number of degrees of freedom, based on the estimator given in Welch (1967)
  data.edof      = fun_compute_edof( ww , nfft , lx , overlap ); % NB: ww already contains window information
  data.edof_info = 'Equivalent number of degrees of freedom in the PSD (Welch, 1967; see also report by Solomon, Jr., 1991)';

  % Update edof if frequency merging was performed
  % See Elgar (1987) for a discussion on the bias introduced to edof in that case
  if merge
    data.edof = data.edof*mg;
  end

  % Computing the 95% confidence interval of the PSD
  alpha          = 1-0.95;
  data.P_CI      = fix(data.edof)./chi2inv([1-alpha/2 alpha/2],fix(data.edof));
  data.P_CI_info = '95% confidence interval of the PSD';

  % b95 (Haubrich, 1965)
  data.b95      = sqrt(6/data.edof);
  data.b95_info = '95% significance level on zero bicoherence: sqrt(6/edof) (Haubrich, 1965)';
end
