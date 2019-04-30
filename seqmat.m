function mat = seqmat(varargin)
%function mat = seqmat([sz1,sz2,...]|sz1,sz2,...)
% Return a matrix of size SZ1xSZ2x... containing a linear sequence of DOUBLE
% starting at 1 and ending at SZ1xSZ2x... Useful for testing, e.g., INTERPN.

  mat = repmat(nan,varargin{:});
  mat(:) = 1:numel(mat);

return;
