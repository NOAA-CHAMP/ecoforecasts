function stn = load_bb_data(fname, stn, staname, plotAllResults)
%function stn = load_bb_data(fname, stn, staname, plotAllResults)
%
% Load and parse a "CLIPS fact-like" SENSOR-DATA file (extension ".bb").
% NOTE: M-func does NOT handle true CLIPS fact files - just sensor data.
%
% 'Stn' is a structure which is returned with one field for each unique
% instrument-sensor combo that is contained in the BB file. Each such field
% is in turn a struct, containing a field 'data' with either a double vector
% or a cell array of strings, as appropriate for that value, and a field
% 'date' containing a vector of datenums, the timestamps of those values.
%
% E.g., from a file that contained only lines like the following:
%     (smkf1 02/06/2006 06 037 2108 ct-shallow salinity    36.18 )
%
% we would receive a struct-of-structs with only these fields in it:
%     result.ct_shallow_salinity.date - double vector of datenums;
% and result.ct_shallow_salinity.data - double vector of salinities.
%
% A .bb file with "facts" for more than one station results in an error,
% unless the optional arg 'staname' is specified - in that case, only data
% from "facts" with that given (five-character) station code are loaded: all
% other facts in the .bb file are IGNORED.
%
% NOTE WELL: If 'stn' is not specified or is an empty array, a new struct is
% returned containing just the loaded data. If 'stn' is a struct that already
% contains a field matching the name of a column in file 'fname', an attempt
% is made to merge that existing data with the newly loaded data. If for any
% reason an existing field does not contain both a 'data' and date' fields,
% however, THAT FIELD IS REPLACED in the struct 'stn' with the new data.
%
% Last Saved Time-stamp: <Fri 2010-12-17 14:57:27  lew.gramer>
%

  result = [];

  if ( ~exist('stn', 'var') || isempty(stn) )
    stn = [];
  end;

  if ( ~exist('staname', 'var') || isempty(staname) )
    staname = [];
  end;

  if ( ~exist('plotAllResults', 'var') || isempty(plotAllResults) )
    plotAllResults = 0;
  end;

  % Load and parse each line of BB file
  fid = fopen(fname, 'r');
  if ( fid < 0 )
    error('Unable to open BB file %s', fname);
  end;
  flds = textscan(fid, '(%[^ ] %f/%f/%f %*f %*f %2f%2f %[^ ] %[^ ] %f)');
  fclose(fid);

  % Take care of any funky chars in station, instrument and sensor names
  flds{1} = regexprep(flds{1}, '\W+', '_');
  flds{7} = regexprep(flds{7}, '\W+', '_');
  flds{8} = regexprep(flds{8}, '\W+', '_');

  % Make sure we know what station we're processing for
  flds{1} = translate_station_name(flds{1});
  if ( length(unique(upper(flds{1}))) > 1 )
    % If file has multiple stations but user gave us no name, give up!
    if ( isempty(staname) )
      error('%s\n%s', ...
            'Facts can only be loaded for one station at a time!', ...
            'Please glean and separate .bb files by station code.');
    % Otherwise, keep all and only facts matching our designated station
    else
      keepidx = find(strcmpi(flds{1}, staname));
      for idx = 1:length(flds)
        flds{idx} = flds{idx}(keepidx);
      end;
    end;
  end;

  % Construct result structure from loaded (station-matching) data
  insts = unique(flds{7});
  sensors = unique(flds{8});
  for ii = 1:length(insts)
    for si = 1:length(sensors)
      idx = find(strcmp(flds{7}, insts{ii}) & strcmp(flds{8}, sensors(si)));
      if ( ~isempty(idx) )
        % Special case for e.g., 'tide_tide', 'barom_barom', etc.
        if ( strcmp(insts{ii}, sensors{si}) )
          component = sensors{si};
        else
          component = [insts{ii} '_' sensors{si}];
        end;
        result.(component).date(1:numel(idx),1) = ...
            datenum(flds{4}(idx), flds{2}(idx), ...
                    flds{3}(idx), flds{5}(idx), flds{6}(idx), 0);
        result.(component).data(1:numel(idx),1) = flds{9}(idx);
      end;
    end;
  end;

  % Make sure we have ONE value per timestamp, before plotting or merging
  sens = fieldnames(result);
  for si = 1:length(sens)
    [srtdts, oldi, newi] = unique(result.(sens{si}).date);
    result.(sens{si}).date(1:numel(oldi),1) = srtdts;
    result.(sens{si}).data(1:numel(oldi),1) = result.(sens{si}).data(oldi);
  end;

  % Finally, remove any obviously INVALID values from each time series
  sens = fieldnames(result);
  for si = 1:length(sens)
    if ( ~iscell(result.(sens{si}).data) )
      rng = valid_var_range(sens{si});
      goodidx = find( rng(1) <= result.(sens{si}).data & ...
                      result.(sens{si}).data <= rng(2) );
      result.(sens{si}).date(1:numel(goodidx),1) = result.(sens{si}).date(goodidx);
      result.(sens{si}).data(1:numel(goodidx),1) = result.(sens{si}).data(goodidx);
    end;
  end;

  % Plot loaded data before merging with existing, if caller wishes
  if ( plotAllResults )
    for si = 1:length(sens)
      figure;
      plot(result.(sens{si}).date, result.(sens{si}).data);
      datetick;
      hdl = title([sens{si} ' vs. time (from ' fname ')']);
      set(hdl, 'Interpreter', 'none');
    end;
  end;

  %DEBUG:  size(result.(sens{1}).data),

  % Merge new data with existing station data struct (if any)
  stn = merge_station_data(stn, result);

  %DEBUG:  size(stn.(sens{1}).data),

return;
