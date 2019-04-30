function datvecs = intersect_ts_data(ts1,ts2,varargin)
%function datvecs = intersect_ts_data(ts1,ts2[,INTERSECT_TSES args...])
%
% Return an Nx2 array DATVECS with the two numeric data vectors for all
% intersecting dates within the two time series structs TS1 and TS2.
%
% CALLS: INTERSECT_TSES.
%
% Last Saved Time-stamp: <Wed 2014-01-08 14:49:58 Eastern Standard Time gramer>

  [ts1,ts2] = intersect_tses(ts1,ts2,varargin{:});
  datvecs = [ts1.data(:),ts2.data(:)];

return;
