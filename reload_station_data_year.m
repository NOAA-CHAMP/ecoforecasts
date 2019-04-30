function stn = reload_station_data_year(stnm,yr,datapath,trimDiagsIntrahourlies)
%function stn = reload_station_data_year(stnm,yr,datapath,trimDiagsIntrahourlies)
%
% Load all station data for station named STNM (e.g., generally reloading a
% saved MAT file), then reload station 10-character columnar file(s) for just
% the year YR, that can be found in the directory DATAPATH. DEFAULT DATAPATH:
% 'MATLABHOME/ecoforecasts/data/'. TRIMDIAGSINTRAHOURLIES if NOT present or
% TRUE, removes DAPS and diagnostic fields (e.g., 'PTemp') and intrahourly
% winds (e.g., 'WS1_10_1' and 'WD1_10_1') before returning.
%
% Last Saved Time-stamp: <Tue 2015-04-28 14:52:56 Eastern Daylight Time gramer>

  % Retrieve data from a path relative to this M-file's local directory
  if ( ~exist('datapath', 'var') || isempty(datapath) )
    datapath = get_ecoforecasts_path('data');
  end;
  if ( ~exist(datapath, 'dir') )
      error('ReloadStationDataYear:NoDataPath', ...
            'The data source path "%s" does not exist!', ...
            datapath);
  end;

  if ( ~exist('trimDiagsIntrahourlies','var') || isempty(trimDiagsIntrahourlies) )
    trimDiagsIntrahourlies = true;
  end;


  stn = load_station_data(stnm);

  % RELOAD files for our target YR
  fpatt = sprintf('%s-%04d*.txt', stnm, yr);
  stafiles = dir(fullfile(datapath, fpatt));
  if ( isempty(stafiles) )
    error('No station 10-col files found matching "%s"!', fpatt);
  end;
  for fidx = 1:length(stafiles)
    fname = fullfile(datapath, [stafiles(fidx).name]);
    disp(sprintf('(Re?)loading 10-col data file %s', fname));
    stn = load_10col_data(fname, stn);
  end;

  % Remove fields not needed for EFs or research - SAVES MEMORY
  if ( trimDiagsIntrahourlies )
    stn = remove_diag_fields(stn);
    stn = remove_intrahour_fields(stn);
  end;

  fname = fullfile(datapath, [stnm '.mat']);
  station = stn;
  disp(['(Re?)saving station data for future runs: ' fname]);
  save(fname, 'station');
  clear station;

return;
