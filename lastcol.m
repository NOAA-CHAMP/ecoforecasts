function vec = lastcol(mtx)
%function vec = lastcol(mtx)
% Stub function for use in, e.g., window_func: returns a vector containing
% the last 1xM column of data from an NxM matrix 'mtx'.
% Last Saved Time-stamp: <Tue 2010-02-02 21:27:14 Eastern Standard Time gramer>

  vec = mtx(:, size(mtx, 2));

return;
