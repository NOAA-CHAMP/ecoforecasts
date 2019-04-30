function station = reload_ndbc_data_year(stnm,yr,datapath)
%function station = reload_ndbc_data_year(stnm,yr,datapath)
%
% Call LOAD_ALL_NDBC_DATA to load all NDBC data for station named STNM
% (generally just reloads a saved MAT file). Then reload NDBC data file of
% the form SPRINTF('%s/%sh%04d.txt',DATAPATH,STNM,YR). These archived NDBC
% C-MAN/SEAKEYS/ICON station data files can be obtained for any given WMO
% 5-character code STNM from the Web site http://www.ndbc.noaa.gov.
% DEFAULT DATAPATH: 'MATLABHOME/ecoforecasts/data/'.
%
% Last Saved Time-stamp: <Sun 2011-10-23 15:31:53  lew.gramer>

  if ( ~exist('datapath','var') || isempty(datapath) )
    % Retrieve data from a path relative to this M-file's local directory
    datapath = get_ecoforecasts_path('data');
  end;

  stnm = lower(stnm);

  matfname = fullfile(datapath, [stnm '-ndbc.mat']);

  % EARLY RETURN. Justification: Greatly simplifies initial case
  if ( ~exist(matfname,'file') )
    warning('No MAT file yet: Loading *all* raw NDBC files from scratch...');
    station = load_all_ndbc_data([],stnm);
    % Little reality check
    if ( ~exist(matfname, 'file') )
      warning('LOAD_ALL_NDBC_DATA saved no MAT file "%s"?!',matfname);
    end;
    return;
  end;

  % Load this year's data
  fname = fullfile(datapath, sprintf('%sh%04d.txt',stnm,yr));
  if ( ~exist(fname, 'file') )
    error('No NDBC year-data file "%s"!', fname);
  else
    disp(['Loading NDBC "h" file: ' fname]);
    station = load_ndbc_data([], fname);
    % Sanity check
    if ( isempty(station) )
      error('No NDBC data for station "%s" year %04d!',stnm,yr);
    end;
  end;

  % Convert winds from [m/s] to [kts] - be consistent with raw data
  if ( isfield(station,'ndbc_wind1_speed') )
    station.ndbc_wind1_speed.data = station.ndbc_wind1_speed.data ./ 0.5144444444;
  end;
  if ( isfield(station,'ndbc_wind1_gust') )
    station.ndbc_wind1_gust.data = station.ndbc_wind1_gust.data ./ 0.5144444444;
  end;


  % Load previously saved historical data for this station.
  disp(['Loading ' matfname]);
  x = load(matfname, 'station');

  % First remove all OLD data for the given year
  flds = fieldnames(station);
  for fldix = 1:length(flds)
    fld = flds{fldix};
    if ( isfield(x.station,fld) && isfield(x.station.(fld),'date') )
      [oldyr,ig,ig] = datevec(x.station.(fld).date);
      x.station.(fld).date(oldyr == yr) = [];
      x.station.(fld).data(oldyr == yr) = [];
    end;
  end;

  % Then merge the year's new data with historical NDBC data
  station = merge_station_data(x.station,station);

  x = []; clear x;

  if ( isstruct(station) )
    if ( ~isfield(station,'station_name') )
      station.station_name = upper(stnm);
    end;
    disp(['Resaving ' matfname]);
    save(matfname, 'station');
  end;

  % Sanity check
  if ( ~isfield(station, 'ndbc_sea_t') )
    % NOT an ERROR: we may wish to analyze pure MET stations sometimes
    warning('Station "%s" data contained no ndbc_sea_t field', stnm);
  end;

return;
