function res = nancumsum(x,varargin)
%function res = nancumsum(x,varargin)
%
% Call CUMSUM (v.) with all non-NaN values of X. Note however that Infinite
% values will still cause the cumulative sum to blow up.
%
% Last Saved Time-stamp: <Tue 2011-11-08 14:43:28  Lew.Gramer>

  x(isnan(x)) = 0;
  res = cumsum(x,varargin{:});

return;
