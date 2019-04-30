function stns = calc_cum_anom(stns, basevar, climvar, cumvar, npers, perdays)
%function stns = calc_cum_anom(stns, basevar, climvar, cumvar, npers, perdays)
%
% Use climatological anomaly of value BASEVAR (a string) for each struct in
% vector STNS, to calculate a cumulative anomaly of that variable. A vector
% field CLIMVAR must already exist in STNS. A new field CUMVAR is added to
% each struct, with the cumulative anomaly over NPERS periods of PERDAYS days
% each. Commonly used to calculate, e.g., Degree Heating Weeks on sea
% temperature: for this, the DEFAULTS of CLIMVAR=[BASEVAR '_monthly_anom'],
% CUMVAR=[BASEVAR '_dhw'], NPERS=12 (weeks), and PERDAYS=7 are handy.
%
% This code can be used to compute DHW based on the methodology of NOAA Coral
% Reef Watch, outlined in the following reference:
%  Liu, G., J. E. Meyer, I. C. Guch, and M. A. Toscano, 2001. NOAA's satellite
%  coral reef bleaching early warning products aimed at local reef sites around
%  the glob. Reef Encounter, 30: 10-13.
% Or at this Web site:
%   http://coralreefwatch.noaa.gov/satellite/methodology/methodology.html#dhw
%
% Last Saved Time-stamp: <Tue 2010-02-09 21:16:54 Eastern Standard Time gramer>

  if ( ~exist('climvar','var') || isempty(climvar) )
    climvar = [basevar '_monthly_clim'];
  end;

  if ( ~exist('npers','var') || isempty(npers) )
    npers = 12;
  end;
  if ( ~exist('perdays','var') || isempty(perdays) )
    perdays = 7;
  end;

  if ( ~exist('cumvar','var') || isempty(cumvar) )
    if ( perdays == 30 )
      cumvar = [basevar '_dhm'];
    elseif ( perdays == 1 )
      cumvar = [basevar '_dhd'];
    else
      cumvar = [basevar '_dhw'];
    end;
  end;

  anomvar = [cumvar '_anom'];

  for ix = 1:length(stns)
    stns(ix).(cumvar) = [];

    if ( ~isempty(stns(ix).(climvar)) && ~isempty(stns(ix).(basevar).date) )
      stns(ix).(anomvar).date = stns(ix).(basevar).date;
      stns(ix).(anomvar).data = ...
          stns(ix).(basevar).data - nanmax(stns(ix).(climvar)(:));

      valsperday = round(1/nanmin(diff(stns(ix).(basevar).date)));

      if ( valsperday <= 0 )
        warning('Invalid field %s.date - skipping station %d', basevar, ix);
        continue;
      end;

      % Round to nearest 0.1oC and blank out trivial anomalies
      % dat = roundn(stns(ix).(anomvar).data, -1);
      dat = stns(ix).(anomvar).data;
      %dat(dat < 1.0) = nan;
      % But using a 1-degree cutoff produces no results for MISST!
      dat(dat < 0.5) = nan;

      begdix = (valsperday*perdays*npers);

      stns(ix).(cumvar).date = stns(ix).(anomvar).date(begdix:end);
      for dix = 1:length(stns(ix).(cumvar).date)
        stns(ix).(cumvar).data(dix,1) = ...
            nansum(dat(dix:(dix+begdix-1))) / (valsperday*perdays);
      end;
    end;
  end;

return;
