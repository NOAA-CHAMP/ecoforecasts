function fh = hist_ts_seasons(ts,varargin)
%function fh = hist_ts_seasons(ts,varargin)
%
% Make a figure with four subplots, each showing the HISTogram of time
% series TS for the corresponding season (v. GET_SEASON), 1 through 4.
%
% Last Saved Time-stamp: <Sun 2013-02-24 14:44:25 Eastern Standard Time gramer>

  if ( ~is_valid_ts(ts) )
    error('First arg must be a valid time series');
  end;

  fh = figure;
  if ( exist('maxigraph.m') == 2 )
    maxigraph;
  end;
  hold on;
  grid on;

  for seas = 1:4
    subplot(2,2,seas);
    dat = ts.data(get_season(ts.date)==seas);
    hist(dat,varargin{:});
    xlabel(['Season ',num2str(seas),'  (N=',num2str(numel(dat)),')']);
    dat=[]; clear dat
  end;
  linkaxes;

  if ( nargout < 1 )
    fh=[]; clear fh;
  end;

return;
