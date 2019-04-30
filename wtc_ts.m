function wtc_ts(ts1,ts2,fix1,fix2,varargin)
%function wtc_ts(ts1,ts2,fix1,fix2[,WTC args...])
%
% Calculate and plot a Wavelet Coherence Estimate, i.e., call WTC (qv.) on
% coincident values in two time series structs TS1 and TS2. Optional FIX1
% and FIX2 are index- or logical vectors - or handles to functions that
% accept a time series struct and RETURN such a vector: if specified, only
% those indices from each time series are included in the analysis.  Optional
% name-value pair 'maxPoints',N may be used to limit WTC analysis to only
% that number of data points (to avoid virtual memory problems): DEFAULT is
% one year of hourly values (24*366). All other args will be passed to WTC.
%
% DEFAULT WTC args: 'MonteCarloCount',20, 'MaxScale',1024, 'BlackandWhite'
%
% Last Saved Time-stamp: <Sun 2012-07-29 17:54:41  lew.gramer>

  if ( ~is_valid_ts(ts1) || ~is_valid_ts(ts2) )
    error('First two args must be valid time series structs (v. IS_VALID_TS)');
  end;

  if ( ~exist('fix1','var') || isempty(fix1) )
    fix1 = 1:numel(ts1.data);
  elseif ( isa(fix1,'function_handle') )
    fix1 = fix1(ts1);
  elseif ( (~isnumeric(fix1) && ~islogical(fix1)) || ~isvector(fix1) )
    error('Optional third arg FIX1 must be a function handle, numeric- or logical-vector!');
  end;
  if ( ~exist('fix2','var') || isempty(fix2) )
    fix2 = 1:numel(ts2.data);
  elseif ( isa(fix2,'function_handle') )
    fix2 = fix2(ts2);
  elseif ( (~isnumeric(fix2) && ~islogical(fix2)) || ~isvector(fix2) )
    error('Optional fourth arg FIX2 must be a function handle, numeric- or logical-vector!');
  end;

  maxPoints = (24*366*1);
  if ( isempty(varargin) )
    wtcargs = {'MonteCarloCount',20, 'MaxScale',1024, 'BlackandWhite'};

  else
    if ( ~strcmpi(varargin{1},'maxPoints') )
      wtcargs = varargin;
    else
      if ( numel(varargin) < 2 || ~isnumeric(varargin{2}) || ~isscalar(varargin{2}) )
        error('Optional string maxPoints must be followed by numeric scalar!');
      end;
      maxPoints = varargin{2};
      wtcargs = varargin(3:end);
    end;
  end;

  dt = min([min(diff(ts1.date(fix1))),min(diff(ts2.date(fix2)))]);

  [ix1,ix2] = intersect_dates(ts1.date(fix1),ts2.date(fix2),dt*0.499);
  ix1 = fix1(ix1);
  ix2 = fix2(ix2);

  %x = [ ts1.date(ix1) , ts1.data(ix1) ];
  %y = [ ts2.date(ix2) , ts2.data(ix2) ];

  %dt = min([min(diff(ts1.date(ix1))),min(diff(ts2.date(ix2)))]);

  dts1 = [ts1.date(ix1(1)):dt:ts1.date(ix1(end))]';
  dat1 = interp1(ts1.date,ts1.data,dts1);
  % dts2 = [ts2.date(ix2(1)):dt:ts2.date(ix2(end))]';
  % Make sure dates of X and Y match
  dts2 = dts1;
  dat2 = interp1(ts2.date,ts2.data,dts2);

  maxix = min(numel(dts1),maxPoints);
  if ( numel(dts1) > maxix )
    warning('Limiting analysis to %s - %s',datestr(dts1(1)),datestr(dts1(maxix)));
  end;

  x = [ dts1(1:maxix) , dat1(1:maxix) ];
  y = [ dts2(1:maxix) , dat2(1:maxix) ];


  %wt(x, 'BlackandWhite', 'MakeFigure',true);
  wtc(x,y, wtcargs{:});

  xlim([dts1(1)-1,dts1(maxix)+1]);
  datetick('keeplimits');

return;
