function stn = station_tmd_tide(stn_or_stnm,dts)
%function stn = station_tmd_tide(stn_or_stnm,dts)
%
% Use TMD_TOOLBOX (HELP TIDE_PRED) to predict tides in [m] for station struct
% STN or station named STNM, using both the global TPXO7.2 solution, and best
% available regional sol'n(s). If STN struct given, assumes STN.station_name
% is a 5-char code for site; if STN.lon and STN.lat fields not found, calls
% GET_STATION_COORDS (v.) to get site longitude/latitude. TMD_TOOLBOX site:
%   http://polaris.esr.org/ptm_index.html
%
% NOTE: For Caribbean and Florida sites, Gulf of Mexico 1/45o Tidal Solution
% of Erofeeva and Egbert matches available in situ data well:
%   http://volkov.oce.orst.edu/tides/Mex.html
%
% Last Saved Time-stamp: <Mon 2018-08-20 12:58:08 Eastern Daylight Time gramer>

  set_more off;

  datapath = get_ecoforecasts_path('data');

  if ( isstruct(stn_or_stnm) )
    stn = stn_or_stnm;
    if ( ~isfield(stn,'station_name') )
      error('STN struct must have "station_name" field!');
    end;
    stnm = stn.station_name;
  elseif ( ischar(stn_or_stnm) )
    stnm = stn_or_stnm;
    stn.station_name = stnm;
  else
    error('First arg must be a station struct or station name string!');
  end;

  tide_fld = '';
  if ( ~exist('dts','var') || isempty(dts) )
    dts = [];
    tide_fld = 'ctd_deep_i_depth';
    if (~isfield(stn,tide_fld)); tide_fld = 'ctd_shallow_i_depth';	end;
    if (~isfield(stn,tide_fld)); tide_fld = 'tide';			end;
    if (~isfield(stn,tide_fld)); tide_fld = 'ndbc_tide';		end;
    if (~isfield(stn,tide_fld)); tide_fld = 'wave_tide_height';		end;
    if (~isfield(stn,tide_fld)); tide_fld = 'sea_t';			end;
    if (~isfield(stn,tide_fld)); tide_fld = 'ndbc_sea_t';		end;
    if (~isfield(stn,tide_fld)); tide_fld = 'ct_shallow_seatemp';	end;
    if (~isfield(stn,tide_fld)); tide_fld = 'wind1_speed';		end;
    if (~isfield(stn,tide_fld)); tide_fld = 'ndbc_wind1_speed';		end;
    if (~isfield(stn,tide_fld))
      warning('Ecoforecasts:TMDTide:NoDates',...
              'No DTS specified and no appropriate comparison field found!');
    else
      dts = stn.(tide_fld).date(1):(1/24):stn.(tide_fld).date(end);
    end;
  end;


  matfname = fullfile(datapath,[lower(stn.station_name) '_tmd_tide.mat']);

  if ( exist(matfname,'file') )
    disp(['Found ' matfname]);
    load(matfname,'station');

    if ( isfield(station,'tmd_tide') )
      stn.tmd_tide = station.tmd_tide;
      % [ix1,ix2] = intersect_dates(dts,station.tmd_tide.date);
      % stn.tmd_tide.date = station.tmd_tide.date(ix2);
      % stn.tmd_tide.data = station.tmd_tide.data(ix2);
    end;
    if ( isfield(station,'tmd_tide_u') );      stn.tmd_tide_u = station.tmd_tide_u;    end;
    if ( isfield(station,'tmd_tide_v') );      stn.tmd_tide_v = station.tmd_tide_v;    end;

    if ( isfield(station,'tpxo_tide') )
      stn.tpxo_tide = station.tpxo_tide;
      % [ix1,ix2] = intersect_dates(dts,station.tpxo_tide.date);
      % stn.tpxo_tide.date = station.tpxo_tide.date(ix2);
      % stn.tpxo_tide.data = station.tpxo_tide.data(ix2);
    end;
    if ( isfield(station,'tpxo_tide_u') );      stn.tpxo_tide_u = station.tpxo_tide_u;    end;
    if ( isfield(station,'tpxo_tide_v') );      stn.tpxo_tide_v = station.tpxo_tide_v;    end;

    if ( isfield(station,'ao_tide') )
      stn.ao_tide = station.ao_tide;
      % [ix1,ix2] = intersect_dates(dts,station.ao_tide.date);
      % stn.ao_tide.date = station.ao_tide.date(ix2);
      % stn.ao_tide.data = station.ao_tide.data(ix2);
    end;
    if ( isfield(station,'ao_tide_u') );      stn.ao_tide_u = station.ao_tide_u;    end;
    if ( isfield(station,'ao_tide_v') );      stn.ao_tide_v = station.ao_tide_v;    end;

    if ( ~isfield(stn,'depth') )
      [ig,ig,stn.depth] = get_station_coords(stn.station_name);
    end;
    if ( ~isfield(stn,'depth') || isempty(find(isfinite(stn.depth))) )
      error('Cannot find site depth(s) for station "%s"!', stn.station_name);
    end;

    if ( isfield(station,'tmd_tide_speed') && isfield(station,'tmd_tide_dir') )
      stn.tmd_tide_speed = station.tmd_tide_speed;
      stn.tmd_tide_dir = station.tmd_tide_dir;
    elseif ( isfield(station,'tmd_tide_u') && isfield(station,'tmd_tide_v') )
      stn.tmd_tide_speed.date = station.tmd_tide_u.date;
      stn.tmd_tide_speed.data = uv_to_spd(station.tmd_tide_u.data,station.tmd_tide_v.data);
      stn.tmd_tide_dir.date = station.tmd_tide_u.date;
      stn.tmd_tide_dir.data = uv_to_dir_curr(station.tmd_tide_u.data,station.tmd_tide_v.data);
    end;

    if ( isfield(station,'tpxo_tide_speed') && isfield(station,'tpxo_tide_dir') )
      stn.tpxo_tide_speed = station.tpxo_tide_speed;
      stn.tpxo_tide_dir = station.tpxo_tide_dir;
    elseif ( isfield(station,'tpxo_tide_u') && isfield(station,'tpxo_tide_v') )
      stn.tpxo_tide_speed.date = station.tpxo_tide_u.date;
      stn.tpxo_tide_speed.data = uv_to_spd(station.tpxo_tide_u.data,station.tpxo_tide_v.data);
      stn.tpxo_tide_dir.date = station.tpxo_tide_u.date;
      stn.tpxo_tide_dir.data = uv_to_dir_curr(station.tpxo_tide_u.data,station.tpxo_tide_v.data);
    end;

    if ( isfield(station,'ao_tide_speed') && isfield(station,'ao_tide_dir') )
      stn.ao_tide_speed = station.ao_tide_speed;
      stn.ao_tide_dir = station.ao_tide_dir;
    elseif ( isfield(station,'ao_tide_u') && isfield(station,'ao_tide_v') )
      stn.ao_tide_speed.date = station.ao_tide_u.date;
      stn.ao_tide_speed.data = uv_to_spd(station.ao_tide_u.data,station.ao_tide_v.data);
      stn.ao_tide_dir.date = station.ao_tide_u.date;
      stn.ao_tide_dir.data = uv_to_dir_curr(station.ao_tide_u.data,station.ao_tide_v.data);
    end;

    if ( isfield(station,'tmd_tide_i_depth') )
      stn.tmd_tide_i_depth = station.tmd_tide_i_depth;
    else
      stn.tmd_tide_i_depth.date = stn.tmd_tide.date;
      stn.tmd_tide_i_depth.data = stn.tmd_tide.data + stn.depth;
    end;

    if ( isfield(station,'tpxo_tide_i_depth') )
      stn.tpxo_tide_i_depth = station.tpxo_tide_i_depth;
    else
      stn.tpxo_tide_i_depth.date = station.tpxo_tide.date;
      stn.tpxo_tide_i_depth.data = station.tpxo_tide.data + stn.depth;
    end;

    if ( isfield(station,'ao_tide_i_depth') )
      stn.ao_tide_i_depth = station.ao_tide_i_depth;
    else
      stn.ao_tide_i_depth.date = station.ao_tide.date;
      stn.ao_tide_i_depth.data = station.ao_tide.data + stn.depth;
    end;

    station = []; clear station;


  else

    if ( ~isfield(stn,'lon') )
      [stn.lon,stn.lat,stn.depth] = get_station_coords(stn.station_name);
    end;
    if ( ~isfield(stn,'lon') || isempty(find(isfinite(stn.lon))) )
      error('Cannot find coordinates for station "%s"!', stn.station_name);
    end;
    if ( ~isfield(stn,'depth') || isempty(find(isfinite(stn.depth))) )
      error('Cannot find site depth(s) for station "%s"!', stn.station_name);
    end;

    if ( isempty(dts) )
      error('No dates were given or inferred to predict tides for!');
    end;

    %DEBUG:
    if (~isempty(tide_fld)); disp(['Will predict dates to match STN.' tide_fld]); end;

    stn.tmd_tide.date = dts(:);
    stn.tmd_tide_u.date = dts(:);
    stn.tmd_tide_v.date = dts(:);

    stn.tpxo_tide.date = dts(:);
    stn.tpxo_tide_u.date = dts(:);
    stn.tpxo_tide_v.date = dts(:);

    stn.ao_tide.date = dts(:);
    stn.ao_tide_u.date = dts(:);
    stn.ao_tide_v.date = dts(:);

    % TIDE_PRED only works easily if RUN IN the tmd_toolbox directory
    tide_mpath = which('tide_pred');
    if ( isempty(tide_mpath) )
      error('No path found to TIDE_PRED function!');
    end;
    [tide_mpath,ig] = fileparts(tide_mpath);
    %DEBUG:  dir(tidepath),

    mydir = pwd;
    cd(tide_mpath);
    stn.tmd_tide.data = tide_pred(fullfile('DATA','Model_Mex'),dts,stn.lat,stn.lon,'z');
    stn.tmd_tide_u.data = tide_pred(fullfile('DATA','Model_Mex'),dts,stn.lat,stn.lon,'u')./100;
    stn.tmd_tide_v.data = tide_pred(fullfile('DATA','Model_Mex'),dts,stn.lat,stn.lon,'v')./100;

    stn.tpxo_tide.data = tide_pred(fullfile('DATA','Model_tpxo7.2'),dts,stn.lat,stn.lon,'z');
    stn.tpxo_tide_u.data = tide_pred(fullfile('DATA','Model_tpxo7.2'),dts,stn.lat,stn.lon,'u')./100;
    stn.tpxo_tide_v.data = tide_pred(fullfile('DATA','Model_tpxo7.2'),dts,stn.lat,stn.lon,'v')./100;

    stn.ao_tide.data = tide_pred(fullfile('DATA','Model_AO'),dts,stn.lat,stn.lon,'z');
    stn.ao_tide_u.data = tide_pred(fullfile('DATA','Model_AO'),dts,stn.lat,stn.lon,'u')./100;
    stn.ao_tide_v.data = tide_pred(fullfile('DATA','Model_AO'),dts,stn.lat,stn.lon,'v')./100;
    cd(mydir);

    stn.tmd_tide_speed.date = stn.tmd_tide_u.date;
    stn.tmd_tide_speed.data = uv_to_spd(stn.tmd_tide_u.data,stn.tmd_tide_v.data);
    stn.tmd_tide_dir.date = stn.tmd_tide_u.date;
    stn.tmd_tide_dir.data = uv_to_dir_curr(stn.tmd_tide_u.data,stn.tmd_tide_v.data);
    stn.tpxo_tide_speed.date = stn.tpxo_tide_u.date;
    stn.tpxo_tide_speed.data = uv_to_spd(stn.tpxo_tide_u.data,stn.tpxo_tide_v.data);
    stn.tpxo_tide_dir.date = stn.tpxo_tide_u.date;
    stn.tpxo_tide_dir.data = uv_to_dir_curr(stn.tpxo_tide_u.data,stn.tpxo_tide_v.data);
    stn.ao_tide_speed.date = stn.ao_tide_u.date;
    stn.ao_tide_speed.data = uv_to_spd(stn.ao_tide_u.data,stn.ao_tide_v.data);
    stn.ao_tide_dir.date = stn.ao_tide_u.date;
    stn.ao_tide_dir.data = uv_to_dir_curr(stn.ao_tide_u.data,stn.ao_tide_v.data);

    stn.tmd_tide_i_depth.date = stn.tmd_tide.date;
    stn.tmd_tide_i_depth.data = stn.tmd_tide.data + stn.depth;
    stn.tpxo_tide_i_depth.date = stn.tpxo_tide.date;
    stn.tpxo_tide_i_depth.data = stn.tpxo_tide.data + stn.depth;
    stn.ao_tide_i_depth.date = stn.ao_tide.date;
    stn.ao_tide_i_depth.data = stn.ao_tide.data + stn.depth;

    disp(['Saving result to ' matfname]);
    station = [];

    station.station_name = stn.station_name;
    station.lon = stn.lon;
    station.lat = stn.lat;
    station.depth = stn.depth;

    station.tmd_tide = stn.tmd_tide; 
    station.tmd_tide_u = stn.tmd_tide_u; 
    station.tmd_tide_v = stn.tmd_tide_v; 
    station.tmd_tide_speed = stn.tmd_tide_speed; 
    station.tmd_tide_dir = stn.tmd_tide_dir; 
    station.tmd_tide_i_depth = stn.tmd_tide_i_depth; 

    station.tpxo_tide = stn.tpxo_tide;
    station.tpxo_tide_u = stn.tpxo_tide_u;
    station.tpxo_tide_v = stn.tpxo_tide_v;
    station.tpxo_tide_speed = stn.tpxo_tide_speed;
    station.tpxo_tide_dir = stn.tpxo_tide_dir;
    station.tpxo_tide_i_depth = stn.tpxo_tide_i_depth;

    station.ao_tide = stn.ao_tide;
    station.ao_tide_u = stn.ao_tide_u;
    station.ao_tide_v = stn.ao_tide_v;
    station.ao_tide_speed = stn.ao_tide_speed;
    station.ao_tide_dir = stn.ao_tide_dir;
    station.ao_tide_i_depth = stn.ao_tide_i_depth;

    save(matfname,'station');
    station = []; clear station;

  end;

  set_more;

return;
