function [newdts,newdat] = gap_filter_ts(dts,dat,hlp,maxgap,tol,maskval)
%function [newdts,newdat] = gap_filter_ts(dts,dat,hlp,maxgap,tol,maskval)
%
% Apply 10th-order Butterworth filter in sequence backward and forward in
% time (using FILTFILT to the hourly time series in vectors DTS (DATENUM) and
% DAT (DOUBLE data), with a HLP-hour low-pass. For any gaps in the original
% time series longer than MAXGAP (DEFAULT: (3.0/24.0)), NaN-fill an HLP/2
% window around each gap in NEWDTS,NEWDAT to avoid "Gibbs ringing".
%
% CALLS (Toolbox):
%  FILTER_TS,FILTER_GAPS (Ecoforecasts); FILTFILT,BUTTER (Signal); INTERP1.
%
% Last Saved Time-stamp: <Wed 2010-11-17 14:23:06 Eastern Standard Time gramer>

  if ( ~exist('maxgap','var') || isempty(maxgap) )
    maxgap = 3.0/24.0;
  end;
  if ( ~exist('tol','var') || isempty(tol) )
    % Use the default defined in FILTER_GAPS (v.)
    tol = [];
  end;
  if ( ~exist('maskval','var') || isempty(maskval) )
    maskval = 0;
  end;

  goodix = find(isfinite(dat(:)));

  %mint = min(diff(dts(:)));
  % Just assume our time series "should be" hourly
  mint = (1.0/24.0);
  x.ts.date = [dts(1):mint:dts(end)]';
  x.ts.data = interp1(dts(goodix),dat(goodix),x.ts.date,'pchip');

  x.fts.date = x.ts.date;
  x.fts.data = filter_ts(x.ts.data,hlp);

  % Remove gaps in original from result
  x.ts.date = dts(goodix);
  x.ts.data = dat(goodix);
  x = filter_gaps(x,'ts','fts',maxgap,(hlp/24.0)/2,tol,maskval);
  newdts = x.fts.date(:);
  newdat = x.fts.data(:);
  x = []; clear x;

return;
