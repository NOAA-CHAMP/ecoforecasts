function stns = calc_clims(stns, basevar)
%function stns = calc_clims(stns, basevar)
%
% Calculate weekly/monthly climatologies for the time series BASEVAR: for
% each station struct in vector STNS; adds .WEEKLY/MONTHLY_CLIM/ANOM fields
% for climatologies (52x1 and 12x1 vectors) and anomalies (Nx1 time series).
%
% NOTE: Stupid MATLAB... If numel(STNS)>1 and all these fields do not exist
% already in structs STNS, an error occurs on return from this function. To
% avoid this, you may initialize these fields with a series of calls like:
%
%   >> [STNS.([BASEVAR '_weekly_clim'])] = deal([]);
%   >> [STNS.([BASEVAR '_weekly_anom'])] = deal([]);
%   >> [STNS.([BASEVAR '_monthly_clim'])] = deal([]);
%   >> [STNS.([BASEVAR '_monthly_anom'])] = deal([]);
%
% Last Saved Time-stamp: <Thu 2010-02-04 14:32:38 Eastern Standard Time gramer>

  wclim = [ basevar '_weekly_clim'];
  wanom = [ basevar '_weekly_anom'];
  mclim = [ basevar '_monthly_clim'];
  manom = [ basevar '_monthly_anom'];

  for ix = 1:length(stns)

    stns(ix).(wclim) = [];
    stns(ix).(wanom) = [];
    stns(ix).(mclim) = [];
    stns(ix).(manom) = [];

    if ( ~isempty(stns(ix).(basevar)) && ~all(isnan(stns(ix).(basevar).data(:))) )
      [yrs,mos,dys] = datevec(stns(ix).(basevar).date);
      jds = datenum(yrs,mos,dys) - datenum(yrs,1,1) + 1;
      wks = floor((jds-1)./7) + 1; wks(wks > 52) = 52;

      for wk = 1:52
        stns(ix).(wclim)(wk,1) = nanmean(stns(ix).(basevar).data(wks==wk));
      end;
      for mn = 1:12
        stns(ix).(mclim)(mn,1) = nanmean(stns(ix).(basevar).data(mos==mn));
      end;

      stns(ix).(wanom).date = stns(ix).(basevar).date;
      stns(ix).(wanom).data = stns(ix).(basevar).data - stns(ix).(wclim)(wks);

      stns(ix).(manom).date = stns(ix).(basevar).date;
      stns(ix).(manom).data = stns(ix).(basevar).data - stns(ix).(mclim)(mos);
    end;

  end;
