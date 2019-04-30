function stn = load_station_data(stn_or_stanm,datapath,trimDiagsIntrahourlies)
%function stn = load_station_data(stn_or_stanm,datapath,trimDiagsIntrahourlies)
%
% Load all station 10-character columnar files and all "CLIPS-like" BB data
% files for station STN_OR_STANM, that can be found in directory DATAPATH.
% DEFAULT DATAPATH: 'MATLABHOME/ecoforecasts/data/'. TRIMDIAGSINTRAHOURLIES
% if NOT present or TRUE, removes DAPS and diagnostic fields (e.g., 'PTemp')
% and intrahourly winds (e.g., 'WS1_10_1' and 'WD1_10_1') before returning.
%
% Last Saved Time-stamp: <Tue 2016-10-11 15:19:52 Eastern Daylight Time lew.gramer>

  set_more off;

  if ( ischar(stn_or_stanm) )
    stanm = stn_or_stanm;
    stn.station_name = stanm;
  elseif ( isstruct(stn_or_stanm) )
    stn = stn_or_stanm;
    stanm = stn.station_name;
  end;
  clear stn_or_stanm;

  % Retrieve data from a path relative to this M-file's local directory
  if ( ~exist('datapath', 'var') || isempty(datapath) )
    datapath = get_ecoforecasts_path('data');
  end;
  if ( ~exist(datapath, 'dir') )
      error('LoadStationData:NoDataPath', ...
            'The data source path "%s" does not exist!', ...
            datapath);
  end;

  if ( ~exist('trimDiagsIntrahourlies','var') || isempty(trimDiagsIntrahourlies) )
    trimDiagsIntrahourlies = true;
  end;


  % Remember that most stations have changed names over the years...
  orig_stanm = translate_station_name(stanm);
  orig_stanm = orig_stanm{:};
  if ( ~strcmpi(stanm,orig_stanm) )
    warning('Loading data for "canonical" station code "%s" instead of "%s"',orig_stanm,stanm);
    stanm = orig_stanm;
  end;
  clear orig_stanm;

  stanm = lower(stanm);

  % If possible, load station data from a MAT file, or...
  fname = fullfile(datapath, [stanm '.mat']);
  if ( exist(fname, 'file') )

    disp(['Found ' fname]);
    x = load(fname, 'station');
    flds = fieldnames(x.station);
    for fldix = 1:length(flds)
      fld = flds{fldix};
      stn.(fld) = x.station.(fld);
    end;
    x = []; clear x;
    clear fname;


    % Verify that we now have the data we think we have
    if ( ~isfield(stn, 'station_name') || ~strcmpi(stn.station_name, upper(stanm)) )
      warning('LoadStationData:WrongName', ...
              'Did we load the wrong MAT data, or is station name "%s" wrong??', ...
              stanm);
    end;

    % Remove fields not needed for EFs or research - SAVES MEMORY
    if ( trimDiagsIntrahourlies )
      stn = remove_diag_fields(stn);
      %stn = remove_intrahour_winds(stn);
      stn = remove_intrahour_fields(stn);
    end;


  % If we absolutely have to, load station data from scratch!
  else

    disp('Loading station data from raw files! Please wait...');
    station.station_name = stn.station_name;

    % Load data for actual station name, and for all known alias names
    all_stanms = translate_station_name(stanm,true);
    if ( isempty(all_stanms) )
      all_stanms = {stanm};
    else
      disp('(Note data will be loaded for all station name aliases.)');
    end;

    for ix = 1:length(all_stanms)
      cur_stanm = lower(all_stanms{ix});
      %DEBUG:
      disp(cur_stanm);

      stafiles = dir(fullfile(datapath, [cur_stanm '-*-clean.csv']));
      for fidx = 1:length(stafiles)
        fname = fullfile(datapath, [stafiles(fidx).name]);
        disp(sprintf('Loading "clean" CSV data file %s', fname));
        station = load_csv_data(station, fname);
      end;
      clear fidx;
      clear fname;
      clear stafiles;

      stafiles = dir(fullfile(datapath, [cur_stanm '-*.txt']));
      for fidx = 1:length(stafiles)
        fn = [stafiles(fidx).name];
        fname = fullfile(datapath, fn);
        if ( ~isempty(regexp(fn,'[-]ctd[-]')) )
          warning('LoadStationData:IgnoreCTD', ...
                  'Ignoring CTD download file "%s"!', ...
                  fname);
        else
          disp(sprintf('Loading 10-col data file %s', fname));
          station = load_10col_data(fname, station);
        end;
      end;
      clear fidx;
      clear fname;
      clear stafiles;

      stafiles = dir(fullfile(datapath, [cur_stanm '-*.xls']));
      for fidx = 1:length(stafiles)
        fn = [stafiles(fidx).name];
        fname = fullfile(datapath, fn);
        if ( ~isempty(regexp(fn,'[-]ctd[-]')) )
          warning('LoadStationData:IgnoreCTD', ...
                  'Ignoring CTD download file "%s"!', ...
                  fname);
        else
          disp(sprintf('Loading FTP XLS data file %s', fname));
          station = load_ftpxls_data(fname, station); 
        end;
      end;
      clear fidx;
      clear fname;
      clear stafiles;

      stafiles = dir(fullfile(datapath, [cur_stanm '*.bb']));
      if ( ~isempty(stafiles) )
        disp('Loading more station data from BB files! Please wait...');
        for fidx = 1:length(stafiles)
          fname = fullfile(datapath, [stafiles(fidx).name]);
        if ( strfind(fname, 'NDBC') )
          disp(sprintf('Skipping NDBC data file "%s"', fname));
        else
          disp(sprintf('Loading "fact" data file "%s"', fname));
          station = load_bb_data(fname, station);
        end;
        end;
      end;
      clear fidx;
      clear fname;
      clear stafiles;

    end; %for ix

    % Verify that we now have the data we think we have
    if ( isempty(station) )
      error('LoadStationData:NoData', ...
            'No raw or MAT data could be loaded for station name "%s"!', ...
            stanm);
    end;
    if ( ~isfield(station, 'Name') || ~isfield(station.Name, 'data') || ...
         ~iscell(station.Name.data) || any(~strcmpi(strtrim(station.Name.data), stanm)) )
      warning('LoadStationData:WrongName', ...
              'Did we load the wrong raw data, or is station name "%s" wrong??', ...
              stanm);
    end;

    % Remove fields not needed for EFs or research - SAVES MEMORY
    if ( trimDiagsIntrahourlies )
      station = remove_diag_fields(station);
      %station = remove_intrahour_winds(station);
      station = remove_intrahour_fields(station);
    end;

    station.station_name = upper(stanm);

    disp('Saving station data to MAT for future runs...');
    save(fullfile(datapath, [stanm '.mat']), 'station', '-v7.3');

    flds = fieldnames(station);
    for fldix = 1:length(flds)
      fld = flds{fldix};
      stn.(fld) = station.(fld);
    end;
    station = []; clear station;

  end;


  set_more;

return;
