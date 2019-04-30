function stn = get_ncep_reanalysis(stn, yr, jds, vars, flds, BASEURL)
%function stn = get_ncep_reanalysis(stn, yr, jds, vars, flds, BASEURL)
%
% Load NCEP Reanalysis of atmosphere-ocean variables VARS (cellstr), directly
% from NOAA/ESRL NCEP Reanalysis netCDF file for YR, range of Julian days JD.
% Arg STN must be a struct with (at least) scalar fields STN.lon and STN.lat:
% new NCEP data are merged into STN for each value in cellstr FLDS (DEFAULT:
% variable names returned by netCDF THREDDS/dodsC ncep.reanalysis query on
% VARS.) Merging is done via timestamp, with MERGE_STATION_DATA (qv.).
%
% For a list of NCEP Reanalysis variables accessible by this routine, see:
%  http://www.esrl.noaa.gov/psd/thredds/catalog/Datasets/ncep.reanalysis/surface_gauss/catalog.html
% (Catalog may also be saved locally with this m-file, as "ncep_reanalysis_catalog.txt".)
%
% If optional BASEURL is a valid THREDDS dataset URL, use that dataset instead of
% 'http://www.esrl.noaa.gov/psd/thredds/catalog/Datasets/ncep.reanalysis/surface_gauss'.
%
% EXAMPLE:
%  >> % Get Specific Humidity at z=2m, Molasses Reef, all days in 1987 and 1988
%  >> s.lat =  25.01;
%  >> s.lon = -80.38;
%  >> s = get_ncep_reanalysis(s, 1987, [], 'shum.2m');
%  >> s = get_ncep_reanalysis(s, 1988, [], 'shum.2m');
%  >> mean(s.ncep_reanalysis_shum.data),
%  ans =
%      0.0154
%  >> plot(s.ncep_reanalysis_shum.date,s.ncep_reanalysis_shum.data);
%  >> datetick;
%
% Last Saved Time-stamp: <Thu 2010-02-18 14:56:31 Eastern Standard Time gramer>


  if ( ~iscell(vars) )
    vars = { vars };
  end;
  % CHEAT! To get actual netCDF var names, strip everything after initial Dot
  varstubs = regexprep(vars,'([^.]*)[.].*','$1');

  if ( ~exist('flds','var') || isempty(flds) )
    flds = strcat('ncep_reanalysis_', varstubs);
  end;

  if ( ~exist('BASEURL','var') || ~ischar(BASEURL) )
    BASEURL = 'http://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/ncep.reanalysis/surface_gauss';
  end;


  stnm = '';
  if ( isfield(stn,'station_name') )
    stnm = stn.station_name;
  elseif ( isfield(stn,'name') )
    stnm = stn.name;
  elseif ( isfield(stn,'code') )
    stnm = stn.code;
  end;

  if ( ~isfield(stn,'lon') || ~isfield(stn,'lat') )
    [stn.lon,stn.lat,stn.depth] = get_station_coords(stnm);
  end;

  lon = stn.lon;
  lon(lon < 0) = 360 + lon(lon < 0);
  lat = stn.lat;

  if ( ~isempty(jds) )
    % Reanalysis is every six hours, daily
    tidx = (( (jds(1)-1) * 4 ) + 1):(jds(end)*4);
  end;


  result = [];

  for ix = 1:length(vars)

    var = vars{ix};
    varstub = varstubs{ix};
    fld = flds{ix};
    if ( ~isfield(stn, fld) || isempty(stn.(fld)) )
      stn.(fld).date = [];
      stn.(fld).data = [];
    end;

    url = sprintf('%s/%s.gauss.%04d.nc', BASEURL, var, yr);
    nc = mDataset(url);
    if ( isempty(nc) )
      error('Error opening URL "%s"!', url);
    end;
    if ( ix == 1 ) % First time only - get coordinates
      if ( ~isfield(stn,'ncep_reanalysis_lonix') || ~isfield(stn,'ncep_reanalysis_latix') )
        latvar = nc{'lat'};
        [ig,latix] = min(abs(latvar(:) - lat));
        lonvar = nc{'lon'};
        [ig,lonix] = min(abs(lonvar(:) - lon));

        stn.ncep_reanalysis_lonix = lonix;
        stn.ncep_reanalysis_latix = latix;
      else
        lonix = stn.ncep_reanalysis_lonix;
        latix = stn.ncep_reanalysis_latix;
      end;
    end;
    if ( ~isempty(jds) )
      if ( ix == 1 ) % First time only - get coordinates
        tvar = nc{'time'};
        t = tvar(tidx);
      end;
      datvar = nc{varstub};
      dat = datvar(tidx,latix,lonix);
    else
      if ( ix == 1 ) % First time only - get coordinates
        tvar = nc{'time'};
        t = tvar(:);
      end;
      datvar = nc{varstub};
      dat = datvar(:,latix,lonix);
    end;
    close(nc); clear nc;

    if ( ix == 1 ) % First time only - get coordinates
      dts = datenum(1,1,1) + (t./24);
    end;

    % Replaced this hack with a call to MERGE_STATION_DATA below
    % stn.(fld).date(end+1:end+length(dat),1) = dts(:);
    % stn.(fld).data(end+1:end+length(dat),1) = dat(:);

    if ( isempty(dat) )
      warning('No data for variable "%s"!', var);
    elseif ( ~all(size(dts) == size(dat)) )
      warning('Skipping variable "%s": time mismatch!', var);
    else
      result.(fld).date = dts;
      result.(fld).data = dat;
    end;

    % This is in here to keep from overheating ESRL's servers!
    % pause(0.2);
    pause(1);

  end; % for ix = 1:length(vars)

  if ( isempty(result) )
    warning('No variables were successfully downloaded!');
  else
    stn = merge_station_data(stn, result);
  end;


return;
