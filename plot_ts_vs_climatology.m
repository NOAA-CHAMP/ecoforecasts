function clm = plot_ts_vs_climatology(ts,varargin)
%function clm = plot_ts_vs_climatology(ts[,clm][,smoothing][,pct][,GRP_TS_ARGS])
%
% FOR A MORE GENERALIZED FUNCTION, try ECOFORECASTS/ANOMALIZE_TS_TO_TS.M

  args = varargin;

  % To use a previously calculated climatology
  if ( numel(args) > 0 && isfield(args{1},'clim') )
    clm = args{1};
    args(1) = [];
  else
    clm = [];
  end;

  % To specify whether climatology should be plotted smooth or stepped
  % This is safe because first arg to ANOMALIZE_TS is never logical
  if ( numel(args) > 0 && islogical(args{1}) )
    smoothing = args{1};
    args(1) = [];
  else
    smoothing = true;
  end;

  % To specify a different Percentile for upper and lower bound plots
  % This is safe because first arg to ANOMALIZE_TS after PCT is always non-numeric
  if ( numel(args) > 0 && isnumeric(args{1}) )
    pct = args{1};
    args(1) = [];
  else
    pct = [];
  end;


  % Calculate climatology if necessary
  if ( isempty(clm) )
    % This odd calling syntax lets us not repeat a long output sequence below
    res = cell(1,8);
    if ( isempty(pct) || isnan(pct) )
      [res{:}] = anomalize_ts(ts,args{:});
    else
      [res{:}] = anomalize_ts(ts,pct,args{:});
    end;
    [clm.anomts,clm.clim.data,clm.clim.date,clm.sd.data,clm.lopct.data,clm.hipct.data,clm.cumfun,clm.per] = res{:};
    clm.sd.date = clm.clim.date;
    clm.lopct.date = clm.clim.date;
    clm.hipct.date = clm.clim.date;
  end;

  % Embroider climatology and bounds into time series with same dates as TS
  if ( ~isfield(clm,'midts') )

    % Interpolate climatology to form a smooth curve
    if ( smoothing )
      dt = median(diff(ts.date));
      clm.smoothclim.date = [datenum(0,1,1):dt:(datenum(1,1,1)-(1.1*dt))]';
      clm.smoothclim.data = interp1(clm.clim.date,clm.clim.data,clm.smoothclim.date)';
      useclim = clm.smoothclim;
    else
      useclim = clm.clim;
    end;

    [ig,ix] = ismember(clm.cumfun(ts.date),useclim.date);

    smoothargs = {'loess'};

    clm.midts.date = ts.date;
    dts = useclim.date(ix);
    dat = useclim.data(ix);
    if ( smoothing )
      % disp(smoothargs);
      % clm.midts.data = smooth(dat,smoothargs{:});
      % disp('spline');
      % clm.midts.data = interp1(useclim.date,useclim.data,dts,'spline');
      % disp('interp1');
      % clm.midts.data = smth;
      clm.midts.data = dat;
    else
      clm.midts.data = dat;
    end;

    if ( isempty(pct) || ~isnan(pct) )
      clm.lopts.date = ts.date;
      dat = clm.midts.data + clm.lopct.data(ix);
      if ( smoothing )
        clm.lopts.data = smooth(dat,smoothargs{:});
      else
        clm.lopts.data = dat;
      end;

      clm.hipts.date = ts.date;
      dat = clm.midts.data + clm.hipct.data(ix);
      if ( smoothing )
        clm.hipts.data = smooth(dat,smoothargs{:});
      else
        clm.hipts.data = dat;
      end;
    end;
  end;

  fmg;
  if ( isfield(clm,'lopts') )
    plot_ts(ts,clm.midts,clm.lopts,':','Color',[.5,.5,.5],'LineW',1.5,clm.hipts,':','Color',[.5,.5,.5],'LineW',1.5);
    legend('TS','Climatology',['\pm',num2str(pct),'%']);
  else
    plot_ts(ts,clm.midts);
  end;

return;
