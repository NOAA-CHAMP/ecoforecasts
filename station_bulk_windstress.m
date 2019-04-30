function stn = station_bulk_windstress(stn,taufld,wfld,wz,dfld,afld,hfld,pfld)
%function stn = station_bulk_windstress(stn,taufld,wfld,wz,dfld,afld,hfld,pfld)
%
% Calculate neutral wind stress [N/m^2] time series STN.(TAUFLD) from wind
% speed [KNOTS] time series STN.(WFLD), taken at height WZ [m]. Adds field
% STN.([WFLD '_mps']) with speed in [m/s]. If any of the optional arguments
% air temperature [oC] STN.(AFLD), Relative Humidity [%] STN.(HFLD), or air
% pressure STN.(PFLD) are given, also uses those in calculating air density
% [kg/m^3] STN.([AFLD '_dens'] and wind stress. CALLS: INTERSECT_ALL_DATES,
% INTERSECT_DATES, KTS2MPS (Ecoforecasts); STRESSTC, AIR_DENS (Air_Sea).
%
% Last Saved Time-stamp: <Sat 2017-06-17 22:51:08 Eastern Daylight Time gramer>

  if ( ~exist('taufld','var') || ~isvarname(taufld) )
    error('Second arg TAUFLD must be a valid fieldname!');
  end;
  if ( ~exist('wfld','var') || ~isfield(stn,wfld) || ~is_valid_ts(stn.(wfld)) )
    error('Third arg WFLD must name a valid time series field of STN!');
  end;

  if ( isfield(stn,taufld) )
    warning('Ecoforecasts:BulkWindstress:Overwrite',...
            'Field STN.(%s) already exists and will be overwritten',taufld);
    stn.(taufld) = [];
  end;


  % Height of wind sensor [m]
  if ( ~exist('wz','var') || isempty(wz) || ~isnumeric(wz) )
    try
      [wz,az,pz,stz,dtz,slz,dlz] = station_instrument_heights(stn.station_name);
    catch
      wz = 10;
    end;
  end;

  % Wind speed [m/s]
  vfld = [wfld '_mps'];
  stn.(vfld).date = stn.(wfld).date;
  stn.(vfld).data = kts2mps(stn.(wfld).data);

  if ( exist('afld','var') && isfield(stn,afld) )

    if ( exist('hfld','var') && isfield(stn,hfld) )
      % rho_a: Air density [kg/m^3]
      rfld = [afld '_dens'];

      if ( exist('pfld','var') && isfield(stn,pfld) )
        % This "magic" prevents really long delays when, e.g., one of these
        % time series has 20 minute sampling and another has 1 hour sampling
        tol = nanmax([nanmin(diff(stn.(afld).date)),...
                      nanmin(diff(stn.(hfld).date)),...
                      nanmin(diff(stn.(pfld).date))]) / 1.5;

        [aix,hix,pix] = ...
            intersect_all_dates(tol,stn.(afld).date,stn.(hfld).date,stn.(pfld).date);
        a.date = stn.(afld).date(aix); a.data = stn.(afld).data(aix);
        h.date = stn.(hfld).date(hix); h.data = stn.(hfld).data(hix);
        p.date = stn.(pfld).date(pix); p.data = stn.(pfld).data(pix);
        rhoa = air_dens(a.data,h.data,p.data);

      else
        [aix,hix] = intersect_dates(stn.(afld).date,stn.(hfld).date);
        a.date = stn.(afld).date(aix); a.data = stn.(afld).data(aix);
        h.date = stn.(hfld).date(hix); h.data = stn.(hfld).data(hix);
        p = [];
        rhoa = air_dens(a.data,h.data);
      end;

      stn.(rfld).date = a.date;
      stn.(rfld).data = rhoa;

      [vix,aix,rix] = ...
          intersect_all_dates([],stn.(vfld).date,stn.(afld).date,stn.(rfld).date);
      v.date = stn.(vfld).date(vix); v.data = stn.(vfld).data(vix);
      a.date = stn.(afld).date(aix); a.data = stn.(afld).data(aix);
      r.date = stn.(rfld).date(rix); r.data = stn.(rfld).data(rix);

      stn.(taufld).date = v.date;
      stn.(taufld).data = stresstc(v.data,wz,a.data,r.data);

    else

      [vix,aix] = intersect_dates(stn.(vfld).date,stn.(afld).date);
      v.date = stn.(vfld).date(vix); v.data = stn.(vfld).data(vix);
      a.date = stn.(afld).date(aix); a.data = stn.(afld).data(aix);

      stn.(taufld).date = v.date;
      stn.(taufld).data = stresstc(v.data,wz,a.data);

    end;

  else

    stn.(taufld).date = stn.(wfld).date;
    stn.(taufld).data = stresstc(stn.(wfld).data,wz);

  end;


  if ( exist('dfld','var') && isfield(stn,dfld) )

    [tix,dix] = intersect_dates(stn.(taufld).date,stn.(dfld).date);
    t.date = stn.(taufld).date(tix); t.data = stn.(taufld).data(tix);
    d.date = stn.(dfld).date(dix); d.data = stn.(dfld).data(dix);

    tauxfld = [taufld '_x'];
    tauyfld = [taufld '_y'];
    stn.(tauxfld).date = t.date;
    stn.(tauyfld).date = t.date;
    [stn.(tauxfld).data,stn.(tauyfld).data] = spddir_to_uv(t.data,d.data);

  end;

return;
