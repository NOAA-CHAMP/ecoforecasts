function wt_ts(ts,ix,varargin)
%function wt_ts(ts,ix[,WT args...])
%
% Calculate and plot a Wavelet spectral power estimate, i.e., call WT (qv.)
% on time series struct TS. Optional IX is an index- or logical vector - or
% handle to a function that accepts a time series struct and RETURNS such a
% vector: if specified, only those indices from TS are included in analysis.
%
% DEFAULT WT args: 'BlackandWhite', 'MakeFigure',true
%
% Last Saved Time-stamp: <Thu 2012-06-21 15:27:34  lew.gramer>

  if ( ~is_valid_ts(ts) )
    error('First arg must be a valid time series struct (v. IS_VALID_TS)');
  end;

  if ( ~exist('ix','var') || isempty(ix) )
    ix = 1:numel(ts.data);
  elseif ( isa(ix,'function_handle') )
    ix = ix(ts);
  elseif ( (~isnumeric(ix) && ~islogical(ix)) || ~isvector(ix) )
    error('Optional second arg IX must be a function handle, numeric- or logical-vector!');
  end;

  if ( ~isempty(varargin) )
    wtargs = varargin;
  else
    %wtargs = {'MaxScale',1024, 'BlackandWhite', 'MakeFigure',true);
    wtargs = {'BlackandWhite', 'MakeFigure',true};
  end;

  %x = [ ts.date(ix) , ts.data(ix) ];

  dt = median(diff(ts.date));

  dts = [ts.date(ix(1)):dt:ts.date(ix(end))]';
  dat = interp1(ts.date,ts.data,dts);
  if ( numel(dat) > numel(ts.data) )
    warning('Points added by interpolation: %d',numel(dat)-numel(ts.data));
  end;

  maxix = min(numel(dts),(24*366*2));
  if ( numel(dts) > maxix )
    warning('Limiting analysis to %s - %s',datestr(dts(1)),datestr(dts(maxix)));
  end;

  x = [ dts(1:maxix) , dat(1:maxix) ];


  wt(x, wtargs{:});

  xlim([dts(1)-1,dts(maxix)+1]);
  if ( exist('datetick3','file') )
    datetick3;
  else
    datetick('keeplimits');
  end;
  h = datacursormode(gcf);
  set(h, 'UpdateFcn', @wt_select_cb);

return;
