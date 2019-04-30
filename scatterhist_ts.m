function h = scatterhist_ts(ts1,ts2,nbins)
%function h = scatterhist_ts(ts1,ts2,nbins)
%
% Plot scatter and univariate histograms (with SCATTERHIST, qv.) for all
% matching timestamps of two "time series" (structs TS1 and TS2 each with
% numeric vector fields .DATE and .DATA). NBINS, if specified, is passed
% through to SCATTERHIST. 
%
% Last Saved Time-stamp: <Sat 2011-06-18 12:47:51  Lew.Gramer>

  [ix1,ix2] = intersect_dates(ts1.date,ts2.date);
  if ( nargin > 2 )
    h = scatterhist(ts1.data(ix1),ts2.data(ix2),nbins);
  else
    h = scatterhist(ts1.data(ix1),ts2.data(ix2));
  end;
  if ( nargout < 1 )
    h = [];
  end;

return;
