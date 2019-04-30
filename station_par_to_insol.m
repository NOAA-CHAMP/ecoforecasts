function stn = station_par_to_insol(stn,parfld,insfld,uswfld,swfld,goodix)
%function stn = station_par_to_insol(stn,parfld,insfld,uswfld,swfld,goodix)
%
% Calculate new fields insolation STN.(INSFLD), upward (reflected) short
% wave radiation STN.(USWFLD), and net shortwave radiative flux STN.(SWFLD),
% using in situ (or other) PAR light field STN.(PARFLD) and published
% conversion constants. If GOODIDX is given, only use those indices from the
% time series STN.(PARFLD) to estimate D/U/SRF. CALLS: PAR_TO_INSOL.
%
% Last Saved Time-stamp: <Wed 2012-10-24 16:21:26 Eastern Daylight Time lew.gramer>

  if ( ~exist('goodix','var') || isempty(goodix) )
    goodix = 1:length(stn.(parfld).date);
  end;

  stn.(insfld).date = stn.(parfld).date(goodix);
  stn.(insfld).data = par_to_insol(stn.(parfld).data(goodix));

  % Initially, just assumed a constant Albedo, chosen in Godfrey et al (1991)
  % as 0.06. Then tried 0.04 based on regressing TOMS GSIP satellite product.
  % For a while, used empirical relationship from regressing NCEP USRF vs.
  % NCEP DSRF and NCEP wind speed. (Dunno why?!) Finally settled on SWHF from
  % AIR_SEA toolbox (v.), for comparability with other published results.

  % % alpha = 0.06;
  % alpha = 0.041;
  % stn.(uswfld).date = stn.(insfld).date;
  % stn.(uswfld).data = stn.(insfld).data .* alpha;

  % [iix,wix] = intersect_dates(stn.(insfld).date,stn.(wndfld).date);
  % alpha = insol_wind_to_albedo(stn.(insfld).data(iix),stn.(wndfld).data(wix));
  % stn.(uswfld).date = stn.(insfld).date(iix);
  % stn.(uswfld).data = stn.(insfld).data(iix) .* alpha;

  [yr,ig,ig] = datevec(stn.(insfld).date);
  yd = stn.(insfld).date - datenum(yr,1,1);
  % PAR reports an hourly mean: assume it is ~ mid-hourly value
  yd = yd - (30/(24*60));

  % NOTE: This funky call requires W longitudes to be > 0!
  [netsw,alpha] = swhf(yd,yr,-stn.lon,stn.lat,stn.(insfld).data);

  %DEBUG:
  % [netsw,alpha,sorad,sunalt,trans] = swhf_debug(yd,yr,-stn.lon,stn.lat,stn.(insfld).data);
  % stn.([swfld '_albedo']).date = stn.(insfld).date;
  % stn.([swfld '_albedo']).data = alpha;
  % stn.([swfld '_sorad']).date = stn.(insfld).date;
  % stn.([swfld '_sorad']).data = sorad;
  % stn.([swfld '_sunalt']).date = stn.(insfld).date;
  % stn.([swfld '_sunalt']).data = sunalt;
  % stn.([swfld '_trans']).date = stn.(insfld).date;
  % stn.([swfld '_trans']).data = trans;
  %DEBUG:

  if ( exist('uswfld','var') && ischar(uswfld) )
    stn.(uswfld).date = stn.(insfld).date;
    stn.(uswfld).data = stn.(insfld).data .* alpha;
  end;
  if ( exist('swfld','var') && ischar(swfld) )
    stn.(swfld).date = stn.(insfld).date;
    stn.(swfld).data = netsw;
  end;

return;
