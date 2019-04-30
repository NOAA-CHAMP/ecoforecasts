function rmse = rmse_tses(ts1,ts2,tol)
%function rmse = rmse_tses(ts1,ts2[,tol])
%
% Return the root-mean-squared difference between coincident data points in
% the two time series STRUCTs TS1 and TS2. Calls INTERSECT_TSES (v.) Optional
% third arg TOL if present is passed through to INTERSECT_TSES. NOTE: Treats
% NaN in either time series as a missing value in both.
%
% Last Saved Time-stamp: <Wed 2018-08-01 17:17:17 Eastern Daylight Time gramer>

  if ( exist('tol','var') )
    [ts1,ts2] = intersect_tses(tol,ts1,ts2);
  else
    [ts1,ts2] = intersect_tses(ts1,ts2);
  end;
  goodix = find(~isnan(ts1.data) & ~isnan(ts2.data));
  rmse = sqrt(sum( (ts1.data(goodix) - ts2.data(goodix)).^2 ));

return;
