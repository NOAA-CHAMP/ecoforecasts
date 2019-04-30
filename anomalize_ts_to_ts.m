function [anomts,climts,loasdts,hiasdts,lopctts,hipctts,cumfun,per] = anomalize_ts_to_ts(ts,varargin)
%function [anomts,climts,loasdts,hiasdts,lopctts,hipctts,cumfun,per] = anomalize_ts_to_ts(ts[,smoothing][,pct][,GRP_TS args])
%
% Call ANOMALIZE_TS (v.), but return time series with time stamps identical
% to those of the ANOMTS, for climatology, climatology minus anomaly standard
% deviation, climatology plus anomaly standard deviation, PCTth (or LOPCTth)
% percentile, 100 - PCTth (or HIPCTth) percentile. Otherwise, arguments and
% return values are the same as ANOMALIZE_TS.
%
% If optional SMOOTHING is LOGICAL (v.) and true (DEFAULT), vary all returned
% time series smoothly between time stamps, e.g., for a daily climatology of
% an hourly TS, do not have a single value for all 24 values of each day.
%
% SAMPLE CALL, producing a 'year-hour' (diurnal-annual) climatology:
%  [anomts,climts,loasdts,hiasdts] = anomalize_ts_to_ts(stn.ndbc_sea_t,@get_jhour_no_leap);
%
% SAMPLE CALL, with smoothed time series and using 7th and 93rd percentiles:
%  [anomts,climts,loasdts,hiasdts] = anomalize_ts_to_ts(stn.ndbc_sea_t,true,7);
%
% NOTE that for an hourly time series TS, both of the sample calls above
% would return the same time series for climatology, "HI" and "LO" SD!
%
% Last Saved Time-stamp: <Tue 2018-02-13 15:40:20 EST lew.gramer>

  args = varargin;

  % Should climatology time series be "stepped" (have one value per CUMFUN
  % period with jumps between), or vary smoothly between timestamps? 
  if ( numel(args) > 0 && islogical(args{1}) )
    % This is safe because neither ANOMALIZE_TS nor GRP_TS can have a value
    % of type LOGICAL as their first (valid) optional arg
    smoothing = args{1};
    args(1) = [];
  else
    smoothing = true;
  end;

  [anomts,clim,tid,asd,lopct,hipct,cumfun,per] = anomalize_ts(ts,args{:});

  if ( smoothing )
    % Smooth time series based on MEDIAN sampling period of original TS
    dt = median(diff(ts.date));
    smoothtid = [datenum(0,1,1):dt:(datenum(1,1,1)-(1.1*dt))]';
    clim = interp1(tid,clim,smoothtid)';
    asd = interp1(tid,asd,smoothtid)';
    lopct = interp1(tid,lopct,smoothtid)';
    hipct = interp1(tid,hipct,smoothtid)';
    tid = smoothtid; clear smoothtid
  end;

  climts.date = anomts.date;
  climts.data = interp1(tid,clim,cumfun(anomts.date));
  loasdts.date = anomts.date;
  loasdts.data = interp1(tid,clim-asd,cumfun(anomts.date));
  hiasdts.date = anomts.date;
  hiasdts.data = interp1(tid,clim+asd,cumfun(anomts.date));
  lopctts.date = anomts.date;
  lopctts.data = interp1(tid,clim+lopct,cumfun(anomts.date));
  hipctts.date = anomts.date;
  hipctts.data = interp1(tid,clim+hipct,cumfun(anomts.date));

return;
