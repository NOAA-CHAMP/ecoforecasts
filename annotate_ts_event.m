function [lh,ax] = annotate_ts_event(varargin)
%function [lh,ax] = annotate_ts_event([ax,][filterfun,][ts[,minV[,drawV[,plotstyle[,ARGS]]]]])
%
% Place markers of PLOTSTYLE (v. PLOT, DEFAULT: 'k^') at each date in time
% series TS where TS.data<=MINV (DEFAULT: PRCTILE(TS.data,3)), at value DRAWV
% (DEFAULT: MINV). If first arg is a HANDLE, make that the default AXES. If
% ISINF(DRAWV) is True, draw markers at y=MIN(YLIM(AX)) (bottom of AXES).
%
% If first or second arg. is a FUNCTION_HANDLE, only plots indices that are
% also returned by FILTERFUN(TS) (DEFAULT: filterfun=TS_BOREAL_WARM, v.).
% Thus by default, only "northern summer events" are annotated. This is,
% e.g., useful for marking UPWELLING events (where TS is a near-bottom sea
% temperature time series), but has application for many numeric criteria.
%
% Optional additional ARGS are passed through to PLOT (v.)
%
% Last Saved Time-stamp: <Mon 2018-03-26 12:53:03 Eastern Daylight Time gramer>

  varargout = {};
  args = varargin;
  if ( ~isempty(args) && ishandle(args{1}) )
    ax = args{1};
    args(1) = [];
  else
    ax = gca;
  end;
  if ( ~isempty(args) && isa(args{1},'function_handle') )
    filterfun = args{1};
    args(1) = [];
  else
    %filterfun = @ts_isfinite;
    filterfun = @ts_boreal_warm;
  end;
  if ( ~isempty(args) && is_ts(args{1}) )
    ts = args{1};
    args(1) = [];
  else
    error('No valid Time Series specified!');
  end;
  if ( ~isempty(args) && isnumeric(args{1}) )
    minV = args{1};
    args(1) = [];
  else
    minV = prctile(ts.data,3);
  end;
  if ( ~isempty(args) && isnumeric(args{1}) )
    drawV = args{1};
    args(1) = [];

    if ( isinf(drawV) )
      if ( drawV > 0 )
        drawV = max(ylim(ax));
      else
        drawV = min(ylim(ax));
      end;
    end;
  else
    drawV = minV;
  end;
  if ( ~isempty(args) && ischar(args{1}) )
    plotstyle = args{1};
    args(1) = [];
  else
    plotstyle = 'k^';
  end;

  ix = intersect(filterfun(ts),find(ts.data<=minV));
  hold on;
  lh=plot(ts.date(ix),repmat(drawV,[1,numel(ix)]),plotstyle,args{:});

return;
