function result = load_webxls_data(fname, stn, plotAllResults)
% function stn = load_webxls_data(fname, stn, plotAllResults)
%
% Load an ICON sensor data file in Excel "web XLS" format (as produced by
% running an IMN query from the "classic" CREWS web page, and saving and
% resaving it in Excel per instructions on the query result page). 'Stn' is a
% structure which is returned with one field for each column in the file: a
% "header line" (if one can be found) of the spreadsheet is assumed to
% contain the name of each variable, and a struct field is either created or
% interleaved in accordingly within 'stn'. Each field is in turn a struct,
% containing a field 'data' containing either a double vector or a cell array
% of strings, as appropriate for that column, and a field 'date' containing a
% vector of datenums, the timestamps of those values.
%
% NOTE WELL: If 'stn' is not specified or is an empty array, a new struct is
% returned containing just the loaded data. If 'stn' is a struct that already
% contains a field matching the name of a column in file 'fname', an attempt
% is made to merge that existing data with the newly loaded data. If for any
% reason an existing field does not contain both a 'data' and date' fields,
% however, THAT FIELD IS REPLACED in the struct 'stn' with the new data.
%
% Last Saved Time-stamp: <Fri 2008-08-22 14:50:11 Eastern Daylight Time gramer>
%


    if ( ~exist('stn', 'var') || isempty(stn) )
        stn = [];
    end;

    if ( ~exist('plotAllResults', 'var') || isempty(plotAllResults) )
        plotAllResults = 0;
    end;


    [ignn, ignt, rawcells] = xlsread(fname);
    clear ignn ignt;


    % Now it's time to play "Where's the actual data?"

    hdrrow = find(strcmp(rawcells(:,1), 'Station') & ...
                  strcmp(rawcells(:,2), 'Date'));
    if ( isempty(hdrrow) )
      hdrrow = find(strcmp(rawcells(:,1), 'Name') & ...
                    strcmp(rawcells(:,2), 'Date'));
    end;
    if ( isempty(hdrrow) )
      error('Cannot find header row in Excel file "%s"?', fname);
    end;
    begrow = hdrrow + 1;
    begcol = 1;

    endcol = length(rawcells(hdrrow,:));
    if ( 3 > endcol || endcol > 200 )
      error('Cannot find last data column in Excel file "%s"?', fname);
    end;

    endrow = find(strncmp('Total observations', rawcells(:,1), ...
                          length('Total observations')), 1);
    if ( isempty(endrow) )
      error('Cannot find last data row in Excel file "%s"?', fname);
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
      if ( strcmpi(hdr{col}, 'date') )
        datecol = col;
      end;
      if ( strcmpi(hdr{col}, 'hour') )
        hourcol = col;
      end;
    end;
    if ( datecol == 0 )
      error('No "date" column found in file "%s"??', fname);
    end;
    if ( hourcol == 0 )
      error('No "hour" column found in file "%s"?', fname);
    end;



    % Construct return structure from parsed data

    % Process the date and hour columns
    alltimestamps = ...
        datenum(strtrim(rawcells(begrow:endrow, datecol)), 'mm/dd/yyyy') ...
        + ([rawcells{begrow:endrow, hourcol}]' / 24.0);

    % Process all other data columns
    nGoodCols = 0;
    for col = begcol:endcol
      % Only process non-blank columns
      if ( ~isempty(hdr{col}) )
        nGoodCols = nGoodCols + 1;
        % Return a vector of numbers where ever possible
        result.(hdr{col}).data = str2double(rawcells(begrow:endrow, col));
        total = length( result.(hdr{col}).data );
        nans = sum( isnan(result.(hdr{col}).data) );
        if ( nans > (total * 0.5) )
          % If most of the values in the column are non-numeric,
          % try to return a cell array of strings instead.
          result.(hdr{col}).data = rawcells(begrow:endrow, col);
        end;
        result.(hdr{col}).date = alltimestamps;
      end;
    end;


    % Remove any obviously INVALID values from each time series
    sens = fieldnames(result);
    for isen = 1:length(sens)
      if ( ~iscell(result.(sens{isen}).data) )
        rng = valid_var_range(sens{isen});
        goodidx = find( rng(1) <= result.(sens{isen}).data & ...
                        result.(sens{isen}).data <= rng(2) );
        result.(sens{isen}).date = result.(sens{isen}).date(goodidx);
        result.(sens{isen}).data = result.(sens{isen}).data(goodidx);
      end;
    end;


    % If caller wishes it, plot all numeric values
    if ( plotAllResults )
      disp(sprintf('Total cols %d, with data %d', length(hdr), nGoodCols));

      % No use in plotting date or time fields!
      sens(strcmpi(sens, 'date')) = [];
      sens(strcmpi(sens, 'year')) = [];
      sens(strcmpi(sens, 'jyear')) = [];
      sens(strcmpi(sens, 'jday')) = [];
      sens(strcmpi(sens, 'hour')) = [];
      sens(strcmpi(sens, 'minute')) = [];

      for isen = 1:length(sens)
        disp(['Plotting ' sens{isen}]);
        if ( isnumeric(result.(sens{isen}).data) )
          figure;
          plot(result.(sens{isen}).date, result.(sens{isen}).data);
          datetick;
          hdl = title([sens{isen} ' (from ' fname ')']);
          % Stupid TeX underscore substitution...
          set(hdl, 'Interpreter', 'none');
        end;
      end;
    end;

%     % Merge new data with existing station data struct (if any)
%     % Note this CAN overwrite data in 'stn' - see m-file help.
%     stn = merge_station_data(stn, result);

return;
