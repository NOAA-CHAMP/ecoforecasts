function stn = load_ftpxls_data(fname,stn)
%function stn = load_ftpxls_data(fname,stn)
%
% Load an ICON sensor data file in Excel "FTP XLS" format (as maintained at
% the site ftp://ftp.aoml.noaa.gov/ocd/pub/jankulak/). STN is a structure
% which is returned with one field for each column in the file: a "header
% line" (if one can be found) of the spreadsheet is assumed to contain the
% name of each variable, and a struct field is either created or interleaved
% in accordingly within 'stn'. Each field is in turn a struct, containing a
% field 'data' containing either a double vector or a cell array of strings,
% as appropriate for that column, and a field 'date' containing a vector of
% datenums, the timestamps of those values.
%
% NOTE WELL: If 'stn' is not specified or is an empty array, a new struct is
% returned containing just the loaded data. If 'stn' is a struct that already
% contains a field matching the name of a column in file 'fname', an attempt
% is made to merge that existing data with the newly loaded data. If for any
% reason an existing field does not contain both a 'data' and date' fields,
% however, THAT FIELD IS REPLACED in the struct 'stn' with the new data.
%
% Last Saved Time-stamp: <Tue 2015-04-28 14:38:42 Eastern Daylight Time gramer>


  if ( ~exist('stn','var') || isempty(stn) )
    stn = [];
  end;

  [ignn, ignt, rawcells] = xlsread(fname);

  % Now it's time to play "Where's the actual data?"
  hdrrow = find(strcmp(rawcells(:,1), 'Station Timestamp') & ...
                strcmp(rawcells(:,2), 'Name'));
  if ( isempty(hdrrow) )
    hdrrow = find(strcmp(rawcells(:,1), 'Station') & ...
                  strcmp(rawcells(:,2), 'Date'));
  end;
  if ( isempty(hdrrow) )
    hdrrow = find(strcmp(rawcells(:,1), 'Name') & ...
                  strcmp(rawcells(:,2), 'Date'));
  end;
  if ( isempty(hdrrow) )
    error('load_ftpxls_data:NoHeader','Cannot find header row in Excel file "%s"?', fname);
  end;
  begrow = hdrrow + 1;
  begcol = 1;

  endcol = length(rawcells(hdrrow,:));
  if ( 3 > endcol || endcol > 200 )
    error('load_ftpxls_data:NoLastCol','Cannot find last data column in Excel file "%s"?',fname);
  end;


  endrow = find(strncmp('Total observations', rawcells(:,1), ...
                        length('Total observations')), 1);
  if ( isempty(endrow) )
    %error('load_ftpxls_data:NoSummary','Cannot find last data row in Excel file "%s"?',fname);
    warning('load_ftpxls_data:NoSummary','No summary row in Excel file "%s"',fname);
    endrow = size(rawcells,1)+1;
  end;
  endrow = endrow - 1;

  % Process column headers
  hdr = rawcells(hdrrow, begcol:endcol);
  for col = begcol:endcol
    % LOWERCASE header label, and take care of any spaces or funky chars
    hdr{col} = regexprep(lower(hdr{col}), '\W+', '_');
    % Translate header label into ICON/G2-ish instrument_sensor name
    hdr{col} = translate_var_name(hdr{col});
  end;

  datecol = 0;
  hourcol = 0;
  for col = begcol:endcol
    result.(hdr{col}) = [];
    if ( strcmpi(hdr{col},'station_timestamp') )
      datecol = col;
      break;
    else
      % This old code is no longer valid for Jank's newer formats
      if ( strcmpi(hdr{col}, 'date') )
        datecol = col;
      end;
      if ( strcmpi(hdr{col}, 'hour') )
        hourcol = col;
      end;
    end;
  end;
  if ( datecol == 0 )
    error('load_ftpxls_data:NoDate','No "date" column found in file "%s"??', fname);
  end;
  if ( hourcol == 0 )
    warning('load_ftpxls_data:NoHour','No "hour" column found in file "%s"', fname);
  end;


  % Construct return structure from parsed data

  % Process the date and hour columns
  if ( hourcol )
    alltimestamps = datenum(rawcells(begrow:endrow,datecol), 'mm/dd/yyyy') ...
        + ([rawcells{begrow:endrow, hourcol}]' / 24.0);
  else
    % Stupid XLS import screws up midnight timestamps
    midnightix = find(cellfun(@isempty, strfind(rawcells(:,datecol),' ') ));
    rawcells(midnightix,datecol) = strcat(rawcells(midnightix,datecol),' 0:00:00 AM');

    alltimestamps = datenum(rawcells(begrow:endrow,datecol));
  end;

  % Process all other data columns
  nGoodCols = 0;
  for col = begcol:endcol
    % Only process non-blank columns
    if ( ~isempty(hdr{col}) )
      nGoodCols = nGoodCols + 1;
      result.(hdr{col}).date = alltimestamps;
      % Return a vector of numbers where ever possible
      nnums = numel(find(cellfun(@isnumeric, rawcells(begrow:endrow,col))));
      if ( nnums == length(result.(hdr{col}).date) )
        result.(hdr{col}).data = [rawcells{begrow:endrow, col}]';
      else
        result.(hdr{col}).data = rawcells(begrow:endrow, col);
      end;
    end;
  end;

  % Remove any obviously INVALID values from each time series
  sens = fieldnames(result);
  for isen = 1:length(sens)
    if ( ~iscell(result.(sens{isen}).data) )
      rng = valid_var_range(sens{isen});
      goodix = find( rng(1) <= result.(sens{isen}).data & ...
                      result.(sens{isen}).data <= rng(2) );
      %DEBUG:      if (length(goodix)<length(result.(sens{isen}).date)); disp({sens{isen},numel(goodix),length(result.(sens{isen}).date)}); if (strcmp(sens{isen},'ctd_deep_seatemp')); keyboard; end; end;
      if ( isempty(goodix) )
        warning('load_ftpxls_data:FieldAllBad','All "%s" data bad: %s',sens{isen},fname);
        result = rmfield(result,sens{isen});
      else
        result.(sens{isen}).date = result.(sens{isen}).date(goodix);
        result.(sens{isen}).data = result.(sens{isen}).data(goodix);
      end;
    end;
  end;

  % Merge new data with existing station data struct (if any)
  % Note this CAN overwrite data in STN - see M-file HELP.
  stn = merge_station_data(stn,result);

return;
