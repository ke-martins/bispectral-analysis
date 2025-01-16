function [ v ] = fun_compute_edof( w , Ns , N , overlap )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computes the spectral estimate effective degrees of freedom, following Percival and Walden.
%
% Inputs:
%   w       - window (windowed timeseries of FFT length)
%   Ns      - bloc length for the FFT
%   N       - total number of points
%   overlap - overlap in % used to compute number of new points in each FFT
%
% Outputs:
%   v       - effective degrees of freedom in a PSD estimate, following Percival and Walden (1993, their Eq. 292b)
%
% September 25, 2024
% Kévin Martins - kevin.martins@cnrs.fr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % Appropriate normalisation
  w = w/sqrt(sum(w.^2));

  % Dealing with overlap
  n  = fix((100-overlap)/100*Ns);

  % Number of blocks
  Nb = floor((N-Ns)/n+1);

  % Computing effective degrees of freedom
  sumh = [];
  for m = 1:Nb-1
    if(Ns-m*n>=1)
      sumh(m) = (1-m/Nb)*abs(sum(w(1:Ns-m*n).*w(1+m*n:Ns)))^2;
    end
  end
  denom = 1+2*sum(sumh);
  v = 2*Nb/denom;
  
  return
end

