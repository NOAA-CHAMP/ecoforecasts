function station = load_all_ndbc_data(station, stnm, yrs)
%function station = load_all_ndbc_data(station, stnm, yrs)
%
% Load all files of form ['data/' STNM 'hYYYY.txt'], where YYYY is a year.
% These archived NDBC C-MAN/SEAKEYS/ICON station data files can be obtained
% for the given the 5-character code STNM from http://www.ndbc.noaa.gov.
%
% ADDED 12 Jan 2011: Also load all files of form ['data/' STNM 'oYYYY.txt'],
% with NDBC-quality controlled "ocean" (FIO CT/CTD) data in them.
%
% Last Saved Time-stamp: <Tue 2016-10-11 15:20:25 Eastern Daylight Time lew.gramer>

  if ( ~exist('stnm','var') || isempty(stnm) )
    stnm = [];
    if ( isstruct(station) )
      if ( isfield(station,'station_name') )
        stnm = station.station_name;
      elseif  ( isfield(station,'name') )
        stnm = station.name;
      elseif  ( isfield(station,'station_code') )
        stnm = station.station_code;
      elseif  ( isfield(station,'code') )
        stnm = station.code;
      end;
    elseif  ( ~isempty(inputname(1)) )
      stnm = inputname(1);
    end;
    stnm = lower(stnm);
  end;
  if ( isempty(stnm) )
    error('STATION struct has no name and no STNM specified!');
  end;

  if ( ~exist('yrs','var') || isempty(yrs) )
    [this_year, ig, ig] = datevec(now);
    yrs = 1984:this_year;
  end;

  % Retrieve data from a path relative to this M-file's local directory
  datapath = get_ecoforecasts_path('data');

  matfname = fullfile(datapath, [stnm '-ndbc.mat']);
  if ( exist(matfname, 'file') )

    stn = station;
    disp(['Found ' matfname]);
    load(matfname, 'station');
    % Make a token attempt to merge structs if necessary
    if ( ~isempty(stn) )
      if ( ~isstruct(stn) )
        warning('Station arg was not a struct: has been overwritten!');
      else
        station = merge_station_data(stn, station);
      end;
    end;

  else

    filesfound = 0;

    disp('Loading data from NDBC "h" files...');
    for ix = 1:length(yrs)
      yr = yrs(ix);
      fname = fullfile(datapath, sprintf('%sh%d.txt', stnm, yr));
      if ( exist(fname, 'file') )
        filesfound = filesfound + 1;
        disp(fname);
        station = load_ndbc_data(station, fname);
      end;
    end;

    fname = fullfile(datapath, sprintf('%sRT.txt', stnm));
    if ( exist(fname, 'file') )
      disp('Loading data from NDBC "RT" file...');
      filesfound = filesfound + 1;
      disp(fname);
      station = load_ndbc_RT_data(station, fname);
    end;

    % Convert winds from [m/s] to [kts] - be consistent with raw data
    if ( isfield(station,'ndbc_wind1_speed') )
      station.ndbc_wind1_speed.data = mps2kts(station.ndbc_wind1_speed.data);
    end;
    if ( isfield(station,'ndbc_wind1_gust') )
      station.ndbc_wind1_gust.data = mps2kts(station.ndbc_wind1_gust.data);
    end;


    disp('Loading data from NDBC "o" files...');
    for ix = 1:length(yrs)
      yr = yrs(ix);
      fname = fullfile(datapath, sprintf('%so%d.txt', stnm, yr));
      if ( exist(fname, 'file') )
        filesfound = filesfound + 1;
        disp(fname);
        station = load_ndbc_data(station, fname);
      end;
    end;

    disp('Loading data from NDBC "a" files...');
    % ADCP (historical format) ocean currents data
    for ix = 1:length(yrs)
      yr = yrs(ix);
      fname = fullfile(datapath, sprintf('%sa%d.txt', stnm, yr));
      if ( exist(fname, 'file') )
        filesfound = filesfound + 1;
        disp(fname);
        station = load_ndbc_data(station, fname);
      end;
    end;

    disp('Loading data from NDBC "r" files...');
    % Solar radiation data
    for ix = 1:length(yrs)
      yr = yrs(ix);
      fname = fullfile(datapath, sprintf('%sr%d.txt', stnm, yr));
      if ( exist(fname, 'file') )
        filesfound = filesfound + 1;
        disp(fname);
        station = load_ndbc_data(station, fname);
      end;
    end;

    if ( filesfound > 0 )
      if ( isstruct(station) )
        if ( ~isfield(station,'station_name') )
          station.station_name = upper(stnm);
        end;
        %save(matfname, 'station');
        save(matfname, 'station', '-v7.3');
      end;
    end;

  end;

  % Sanity checks
  if ( isempty(station) )
    error('No NDBC data or empty MAT file for station "%s"?!?', stnm);
  end;
  if ( ~isfield(station, 'ndbc_sea_t') )
    % Not an error: we may wish to analyze pure MET stations some day?
    warning('Station "%s" data contained no sea temperatures?!?', stnm);
  end;

return;
