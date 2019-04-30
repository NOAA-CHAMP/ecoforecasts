function stn = station_stokes_drift(stn,ssfld,sdfld,sufld,svfld,wsfld,wdfld,hsfld,tpfld,tdfld,dirMethod)
%function stn = station_stokes_drift(stn,ssfld,sdfld,sufld,svfld,wsfld,wdfld,hsfld,tpfld,tdfld,dirMethod)
%
% Estimate Stokes drift (residual wave-driven surface circulation) in new
% time series fields STN.(SSFLD) for speed [m/s] and STN.(SDFLD) for TARGET
% direction, based on wind speed [kts] STN.(WSFLD), and wind SOURCE direction
% STN.(WDFLD), and significant wave height [m] STN.(HSFLD). If provided, also
% uses peak wave period [s] STN.(TPFLD) to estimate cutoff frequency for
% incident waves at site, and peak wave SOURCE direction [degT] STN.(TDFLD)
% to estimate Stokes drift resultant vector direction. As with most STATION_
% functions, all time series dates need not match so long as they overlap.
% Also adds 'u' (SUFLD) and 'v' (SVFLD) vector components for Stokes drift.
%
% NOTE any of HSFLD, TPFLD,TDFLD may also be given as a scalar constant value.
%
% See Ardhuin et al, JPO, 2009, eqn. (7).
% NOTE: Unfortunately resultant direction does not seem well enumerated by
% Arduin et al., so code below reflects various guesses at Stokes direction.
%
% CALLS: STOKES_DRIFT, INTERSECT_ALL_DATES, INTERSECT_DATES.
%
% Last Saved Time-stamp: <Fri 2012-08-03 10:47:34  lew.gramer>

  % In Ardhuin, "cutoff" is Bragg wave frequency for their HF radar?? Assume
  % it is instead a "reasonable upper bound" for wave frequency at this site.
  if ( ~exist('tpfld','var') || isempty(tpfld) )
    tpfld = 2;
  end;
  % If wave direction not specified, assume it equals wind source direction
  if ( ~exist('tdfld','var') || isempty(tdfld) )
    tdfld = wdfld;
  end;
  % If  method not specified to determine Stokes drift direction, do a
  % significant wave-height weighted average of wind and wave directions
  if ( ~exist('dirMethod','var') || isempty(dirMethod) )
    dirMethod = 'blend';
  end;

  if ( isnumeric(hsfld) )
    hss.date = stn.(wsfld).date;
    hss.data = repmat(hsfld, size(stn.(wsfld).data));
  else
    hss.date = stn.(hsfld).date;
    hss.data = stn.(hsfld).data;
  end;
  if ( isnumeric(tpfld) )
    tps.date = stn.(wsfld).date;
    tps.data = repmat(tpfld, size(stn.(wsfld).data));
  else
    tps.date = stn.(tpfld).date;
    tps.data = stn.(tpfld).data;
  end;
  if ( isnumeric(tdfld) )
    tds.date = stn.(wsfld).date;
    tds.data = repmat(tdfld, size(stn.(wsfld).data));
  else
    tds.date = stn.(tdfld).date;
    tds.data = stn.(tdfld).data;
  end;


  [wsix,wdix,hsix,tpix,tdix] = ...
      intersect_all_dates([],stn.(wsfld).date,stn.(wdfld).date,hss.date,tps.date,tds.date);

  w = stn.(wsfld).data(wsix);
  wd = stn.(wdfld).data(wdix);
  hs = hss.data(hsix);
  tp = tps.data(tpix);
  % Assume quasi-Eulerian flow direction is always coherent with waves?!
  td = tds.data(tdix);

  % % Project wind onto peak wave source direction?!
  % w = w .* abs(cosd(wd - td));

  % Calculate Stokes drift speed from Ardhuin et al 2009 formula
  stn.(ssfld).date = stn.(wsfld).date(wsix);
  % NOTE: STOKES_DRIFT calls KTS2MPS on winds for us
  stn.(ssfld).data = stokes_drift(w,hs,tp);

  switch ( dirMethod )
   case 'wave',
    % Assume quasi-Eulerian flow direction is always coherent with waves?!
    d = td;

   case 'wind+35',
    % Assume quasi-Eulerian flow directed 35o to right of wind direct?!
    d = wd + 35;
    d(d>=360) = d(d>=360) - 360;

   case 'wind',
    % Assume quasi-Eulerian flow direction is always coherent with wind?!
    % (This is OK if the wind data is low-pass filtered, e.g., 72hlp)
    d = wd;

   case 'blend',
    % Quasi-Eulerian flow direction as a blend of wind and wave dirs,
    % secant-weighted by Hs: WAVEFRAC ranges from 0 to max ~88.6%.
    [tdu,tdv] = spddir_to_uv(w,td);
    [wdu,wdv] = spddir_to_uv(w,wd);
    wavefrac=min( (hs./3).*(asec(2-0.1143)), asec(2-0.1143) );
    du = (tdu.*wavefrac) + (wdu.*(1-wavefrac));
    dv = (tdv.*wavefrac) + (wdv.*(1-wavefrac));
    d = uv_to_dir(du,dv);

  end;

  % Convert wave/wind SOURCE direction into Stokes TARGET direction
  d = d - 180;
  d(d<0) = d(d<0) + 360;

  stn.(sdfld).date = tds.date(tdix);
  stn.(sdfld).data = d;


  [ssix,sdix] = intersect_dates(stn.(ssfld).date,stn.(sdfld).date);

  stn.(sufld).date = stn.(ssfld).date(ssix);
  stn.(svfld).date = stn.(ssfld).date(ssix);
  [stn.(sufld).data,stn.(svfld).data] = ...
      spddir_to_uv_curr(stn.(ssfld).data(ssix),stn.(sdfld).data(sdix));

return;
