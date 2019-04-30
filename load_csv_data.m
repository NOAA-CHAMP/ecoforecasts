function stn = load_csv_data(stn, fname)
%function stn = load_csv_data(stn, fname)
% 
% Load an ICON sensor data file in an arbitrary CSV (comma separated) format
% (e.g., from manually QC'd historical data files on CoRIS). 'Stn' is a
% structure which is returned with one field for each column in the file:
% header (first line) of  file is assumed to contain the name of each
% variable, and a struct field is either created or interleaved in
% accordingly within 'stn'. Each field is in turn a struct, containing a
% field 'data' containing either a double vector or a cell array of strings,
% as appropriate for that column, and a field 'date' containing a vector of
% datenums, the timestamps of those values.
%
% NOTE WELL: If 'stn' is not specified or is an empty array, a new struct is
% returned containing just the loaded data. If 'stn' is a struct that already
% contains a field matching the name of a column in file 'fname', an attempt
% is made to merge that existing data with the newly loaded data. If for any
% reason an existing field by that name does not contain both a 'data' and
% 'date' fields, THE FIELD IS REPLACED in the struct 'stn' with the new data.
%
% Last Saved Time-stamp: <Wed 2009-02-18 21:18:54 Eastern Standard Time gramer>


  dat = importdata(fname);

  % Idiotic freakin' MATLAB - make sure data and textdata have same size!
  newdata = repmat(nan, size(dat.textdata));
  oldix = 1;
  for ix = 1:size(newdata,2)
    if ( isempty(dat.textdata{2,ix}) )
      newdata(2:end, ix) = dat.data(:, oldix);
      oldix = oldix + 1;
    end;
  end;
  dat.data = newdata;
  clear newdata;


  hdrs = dat.textdata(1,:);

  str='name'; strl=length(str); nameix = find(strncmpi(hdrs, str, strl), 1);
  % Avoid problems with stations being "renamed" over the years...
  if ( ~isempty(nameix) )
    dat.textdata(2:end, nameix) = translate_station_name(dat.textdata(2:end, nameix));
  end;

  str='date'; strl=length(str); dateix = find(strncmpi(hdrs, str, strl), 1);
  str='jyear'; strl=length(str); jyearix = find(strncmpi(hdrs, str, strl), 1);
  if ( isempty(jyearix) )
    str='year'; strl=length(str); jyearix = find(strncmpi(hdrs, str, strl), 1);
  end;
  str='jday'; strl=length(str); jdayix = find(strncmpi(hdrs, str, strl), 1);
  str='time'; strl=length(str); timeix = find(strncmpi(hdrs, str, strl), 1);
  str='hour'; strl=length(str); hourix = find(strncmpi(hdrs, str, strl), 1);
  str='min'; strl=length(str); minix = find(strncmpi(hdrs, str, strl), 1);

  for ix = 1:length(hdrs)
    % Try to 'translate' label into the equivalent ICON/G2 variable name:
    % 'WaterT' becomes 'sea_t', 'Sea1m' is 'ctd_shallow_seatemp', etc.
    hdrs{ix} = translate_var_name(hdrs{ix});
    % Take care of any funky chars in header label
    hdrs{ix} = regexprep(hdrs{ix}, '\W+', '_');
  end;


  % Always try to use the least ambiguous date format (Jyear/Jday)
  if ( ~isempty(jyearix) && ~isempty(jdayix) )
    yr = dat.data(2:end, jyearix);
    % Some formats have e.g., '2' in the Jyear column for 2002
    if ( length(find(yr < 100)) > (0.2 * length(yr)) )
      yr(50 <= yr & yr < 100) = yr(yr > 50) + 1900;
      yr(0 <= yr & yr < 50) = yr(yr < 50) + 2000;
    end;
    jd = dat.data(2:end, jdayix);
    [ig, mo, dy] = datevec(datenum(yr,1,1) + jd - 1);
  elseif ( ~isempty(dateix) )
    error('DATE PARSING NOT YET IMPLEMENTED!');
  else
    error('Found no header fields matching "date", or "jyear" and "jday"!');
  end;

  if ( ~isempty(minix) )
    mn = dat.data(2:end,minix);
  else
    mn = 0;
  end;

  if ( ~isempty(hourix) )
    hr = dat.data(2:end,hourix);
    if ( length(find(hr > 24)) > (0.2 * length(hr)) )
      % Some early formats had e.g., '1800' for 6pm, instead of '18'
      mn = mod(hr, 100);
      hr = fix(hr ./ 100);
    end;
  elseif ( ~isempty(timeix) )
    error('TIME PARSING NOT YET IMPLEMENTED!');
  else
    error('Found no header field matching "hour" or "time"!');
  end;

  result.datenum = datenum(yr, mo, dy, hr, mn, 0);

  for ix = 1:length(hdrs)
    hdr = hdrs{ix};
    result.(hdr) = [];
    result.(hdr).date = result.datenum;
    % Numeric data
    if ( isempty(dat.textdata{2,ix}) )
      result.(hdr).data = dat.data(2:end,ix);
    % Text data
    else
      result.(hdr).data = dat.textdata(2:end,ix);
    end;
  end;

  % Global datenum field no longer needed
  result = rmfield(result, 'datenum');

  % Finally, remove any obviously INVALID values from each time series
  sens = fieldnames(result);
  for isen = 1:length(sens)
    if ( ~iscell(result.(sens{isen}).data) )
      rng = valid_var_range(sens{isen});
      goodidx = find( rng(1) <= result.(sens{isen}).data & ...
                      result.(sens{isen}).data <= rng(2) & ...
                      ~isnan(result.(sens{isen}).data) );
      result.(sens{isen}).date = result.(sens{isen}).date(goodidx);
      result.(sens{isen}).data = result.(sens{isen}).data(goodidx);
    end;
  end;


  % Merge new data with existing station data struct (if any)
  % Note this CAN overwrite data in 'stn' - see m-file help.
  stn = merge_station_data(stn, result);


return;
