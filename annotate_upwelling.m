function varargout = annotate_upwelling(varargin)
%function varargout = annotate_upwelling([ax,][summerfun,][ts[,minT[,drawT[,plotstyle[,ARGS]]]]])
%
% Place markers of PLOTSTYLE (v. PLOT, DEFAULT: 'k^') at each date in time
% series TS where TS.data<=MINT (DEFAULT: PRCTILE(TS.data,3)), at value DRAWT
% (DEFAULT: MINT). If first arg is a HANDLE, make that the default AXES.
%
% If first or second arg. is a FUNCTION_HANDLE, only plots indices that are
% also returned by SUMMERFUN(TS): thus, only "northern summer upwelling" is
% annotated, by default (DEFAULT: summerfun=TS_BOREAL_WARM, v.).
%
% Optional additional ARGS are passed through to PLOT (v.)
%
% Last Saved Time-stamp: <Mon 2018-03-05 14:59:40 EST lew.gramer>

  varargout = {};
  args = varargin;
  if ( ~isempty(args) && ishandle(args{1}) )
    ax = args{1};
    args(1) = [];
  else
    ax = gca;
  end;
  if ( ~isempty(args) && isa(args{1},'function_handle') )
    summerfun = args{1};
    args(1) = [];
  else
    %summerfun = @ts_isfinite;
    summerfun = @ts_boreal_warm;
  end;
  if ( ~isempty(args) && is_ts(args{1}) )
    ts = args{1};
    args(1) = [];
  else
    error('No valid Time Series specified!');
  end;
  if ( ~isempty(args) && isnumeric(args{1}) )
    minT = args{1};
    args(1) = [];
  else
    minT = prctile(ts.data,3);
  end;
  if ( ~isempty(args) && isnumeric(args{1}) )
    drawT = args{1};
    args(1) = [];
  else
    drawT = minT;
  end;
  if ( ~isempty(args) && ischar(args{1}) )
    plotstyle = args{1};
    args(1) = [];
  else
    plotstyle = 'k^';
  end;

  ix = intersect(summerfun(ts),find(ts.data<=minT));
  hold on;
  plot(ts.date(ix),repmat(drawT,[1,numel(ix)]),plotstyle,args{:});

return;
