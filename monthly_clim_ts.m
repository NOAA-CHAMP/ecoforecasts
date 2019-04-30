function [monts,clim,anom] = monthly_clim_ts(ts)
%function [monts,clim,anom] = monthly_clim_ts(ts)
%
% Calculate robust monthly mean time series, climatology (12x1 vector), and
% anomaly for the time series TS (must have fields TS.date and TS.data).
%
% Last Saved Time-stamp: <Fri 2011-01-28 08:27:51  lew.gramer>

  [yr,mo,dy] = datevec(ts.date);
  dts = datenum(yr,mo,1);
  dat = real(ts.data);

  % Monthly mean
  monts.date = unique(dts);
  monts.data = grpstats(dat,dts);

  % nboots = 30;

  % mn = bootstrp(nboots,@grpstats,dat,dts);
  % % Bug in BOOTSTRP?? Returns square matrix, not NBOOTS x length(DAT)!
  % %monts.data = nanmean(mn)';
  % monts.data = nanmean(mn(1:nboots,:))';

  if ( nargout > 1 )
    % Climatology
    clim = grpstats(dat,mo);
    % clim = nanmean(bootstrp(nboots,@grpstats,dat,mo))';
  end;
  if ( nargout > 2 )
    % Monthly anomaly (mean minus climatology)
    [yr,mo,dy] = datevec(monts.date);
    moix = grp2idx(mo);
    anom.date = monts.date;
    anom.data = monts.data - clim(moix);
  end;

return;
