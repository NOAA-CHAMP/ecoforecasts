function stn = load_10col_data(fname, stn, plotAllResults, startpos)
%function stn = load_10col_data(fname, stn, plotAllResults)
%
% Load an ICON sensor data file in 10-character columnar format (standard
% output format of the ICON Java DataParser). STN is a structure which is
% returned with one field for each column in the file: header (first line) of
% file is assumed to contain the name of each variable, and a struct field is
% either created or interleaved in accordingly within STN. Each field is in
% turn a struct, containing a field '.data' with either a double vector or a
% cell array of strings, as appropriate for that column, and a field '.date'
% containing a vector of DATENUMs, the timestamps for those values.
%
% NOTE WELL: If STN is not specified or is an empty array, a new struct is
% returned containing just the loaded data. If STN is a struct that already
% contains a field matching the name of a column in file FNAME, an attempt
% is made to merge that existing data with the newly loaded data. If for any
% reason an existing field does not contain both a '.data' and '.date' field,
% however, THAT FIELD IS REPLACED in the struct STN with the new data.
%
% Last Saved Time-stamp: <Sun 2012-10-21 15:19:00 Eastern Daylight Time lew.gramer>


% NOTE: PRIVATE last argument STARTPOS, if specified, forces the function
% to FSEEK that many bytes into the file FNAME, before beginning parsing.
% This arg allows LOAD_10COL_DATA to call itself recursively on files having
% two or more sections with data in different formats (e.g., as produced by
% the dataparser whenever a station's configuration changes during a year).

    %DEBUG:    disp(['Entered ' upper(mfilename)]);

    if ( ~exist('stn', 'var') || isempty(stn) )
        stn = [];
    end;

    if ( ~exist('plotAllResults', 'var') || isempty(plotAllResults) )
        plotAllResults = false;
    end;

    if ( ~exist('startpos', 'var') || ~isnumeric(startpos) )
        startpos = 0;
    end;

    % Open and parse file
    fid = fopen(fname, 'r');
    if ( fid < 0 )
        error('Unable to open 10-char file %s', fname);
    end

    if ( startpos > 0 )
      if ( fseek(fid, startpos, -1) ~= 0 )
        error('Cannot seek to new-format position %d in 10-char file %s', ...
              startpos, fname);
      end;
    end;

    % Parse header line first, column by column
    line = fgets(fid);
    % Note: This value may change between recursive calls!
    length_of_line = length(line);
    line = regexprep(line,'[\n\r]','');

    idx = 0;
    stmpidx = 0;
    nameidx = 0;
    dateidx = 0; houridx = 0; minuteidx = 0;
    yearidx = 0; jdayidx = 0;
    for curs = 1:10:length(line)
        idx = idx + 1;
        if ( curs+9 > length(line) )
          error('Not a 10-character columnar-format file?? Near line %d: "%s"',curs,fname);
        end;
        hdr{idx} = strtrim(line(curs:curs+9));
        % Try to make sure we have what looks like a valid header line
        if ( ~isnan(str2double(hdr{idx})) )
            error('Invalid header line?? Field %d is the number "%s"!', ...
                  idx, hdr{idx});
        end;

        % Try to 'translate' label into the equivalent ICON/G2 variable name:
        % 'WaterT' becomes 'sea_t', 'Sea1m' is 'ctd_shallow_seatemp', etc.
        hdr{idx} = translate_var_name(hdr{idx});

        % Take care of any funky chars in header label
        hdr{idx} = regexprep(hdr{idx}, '\W+', '_');

        % To catch Mike.Jankulak@noaa.gov's weird 20-char first column
        if ( idx > 1 && strcmpi(strtrim(hdr{idx}), 'timestamp') )
          %DEBUG:          disp('Found timestamp hdr!');
          idx = idx - 1;
          stmpidx = idx;
        end;

        if ( strcmpi(hdr{idx}, 'name') ); nameidx = idx; end;

        if ( strcmpi(hdr{idx}, 'date') ); dateidx = idx; end;
        if ( strcmpi(hdr{idx}, 'yr') ); yearidx = idx; end;
        if ( strcmpi(hdr{idx}, 'year') ); yearidx = idx; end;
        if ( strcmpi(hdr{idx}, 'jyear') ); yearidx = idx; end;
        if ( strcmpi(hdr{idx}, 'jday') ); jdayidx = idx; end;
        if ( strcmpi(hdr{idx}, 'hr') ); houridx = idx; end;
        if ( strcmpi(hdr{idx}, 'hour') ); houridx = idx; end;
        if ( strcmpi(hdr{idx}, 'mn') ); minuteidx = idx; end;
        if ( strcmpi(hdr{idx}, 'min') ); minuteidx = idx; end;
        if ( strcmpi(hdr{idx}, 'minute') ); minuteidx = idx; end;
    end

    if ( stmpidx == 0 )
      if ( dateidx == 0 && (yearidx == 0 || jdayidx == 0) )
        error('Need one of "timestamp", "date", or "year" and "jday" in file "%s"??', fname);
      end;
      if ( houridx == 0 )
        error('No "timestamp" or "hour" column found in file "%s"??', fname);
      end;
    end;

    % Now read all data lines at one blow
    fmt = '';
    for idx = 1:length(hdr)
        if ( idx == stmpidx )
          fmt = [ fmt '%20s' ];
        else
          fmt = [ fmt '%10s' ];
        end;
    end

    vals = textscan(fid, fmt, 'Whitespace', '');

    fclose(fid);

    if ( length(vals) < 2 )
      error('Invalid 10-char file "%s"? (Less than 2 cols)', fname);
    end;

    frames = find(strcmp(strtrim(vals{1}), 'Date'));
    if ( isempty(frames) )
      % To catch Mike.Jankulak@noaa.gov's weird new 20-char first column
      if ( stmpidx > 0 )
        %DEBUG:        disp(['Seeking frames in column ' num2str(stmpidx)]);
        frames = find(cellfun(@isscalar,strfind(lower(vals{stmpidx}),'timestamp')));
      end;
    end;
    if ( ~isempty(frames) )

      % If we're in a recursion, don't confuse user with multiple warnings
      if ( startpos == 0 )
        warning('load_10col_data:FormatChange', ...
                'Found additional headers at line(s) [%s]!\n %s\n%s', ...
                num2str(frames'), fname, ...
                'Will try to process as new formats...');
      end;

      for idx = 1:length(vals)
        vals{idx} = vals{idx}(1:frames(1)-1);
      end;

    end;

    % Avoid problems with stations being "renamed" over the years...
    if ( nameidx ~= 0 )
      vals{nameidx} = translate_station_name(vals{nameidx});
    end;


    % Construct return structure from parsed data

    % Ensure 'datenum' is always the first field in result
    result.datenum = [];

    nGoodCols = 0;
    for idx = 1:length(hdr)
        % Only process non-blank columns
        if ( ~isempty(hdr{idx}) )
            nGoodCols = nGoodCols + 1;
            % Return a vector of numbers where ever possible
            result.(hdr{idx}) = str2double(vals{idx});
            total = length( result.(hdr{idx}) );
            nans = sum( isnan(result.(hdr{idx})) );
            emps = sum( isempty(vals{idx}) );
            if ( (nans > (emps * 1.1)) && (nans > (total * 0.5)) )
                % If most of the values in the column are non-numeric but
                % non-empty, then return a cell array of strings instead.
                result.(hdr{idx}) = vals{idx};
            end;
        end;
    end;


    % Calculate a MATLAB datenum for each line of the file

    if ( stmpidx ~= 0 )

      cstr = textscan(char(vals{stmpidx})', '%02.0f/%02.0f/%04.0f %02.0f:%02.0f:%02.0f');
      [mo,dy,yr,hr,mn,sc] = deal(cstr{:});
      result.datenum = datenum(yr,mo,dy,hr,mn,sc);

    else

      if ( minuteidx == 0 )
        minutes = '0';
      else
        minutes = vals{minuteidx};
      end;

      if ( dateidx == 0 )
        result.datenum = datenum(str2double(vals{yearidx}), 1, 1, ...
                                 str2double(vals{houridx}), ...
                                 str2double(minutes), 0) + ...
                                   str2double(vals{jdayidx}) - 1;
      else
        result.datenum = ...
            datenum(strcat(strtrim(vals{dateidx}), {' '}, ...
                           strtrim(vals{houridx}), ':', ...
                           strtrim(minutes)));
      end;

    end;

    if ( any(diff(result.datenum) <= 0) )
      % If we're in a recursion, don't confuse user with multiple warnings
      if ( startpos == 0 )
        warning('load_10col_data:DateOverlap', ...
                'Unsorted or overlapping timestamps in "%s"!', fname);
      end;
    end;

    % Now assign a separate hourly timestamp to EACH field in the struct
    % Need to do this to allow merging with other data files (that may
    % contain different subsets of the possible sensors) later on...
    sens = fieldnames(result);
    sens(strcmpi(sens, 'datenum')) = [];
    for isen = 1:length(sens)
        data = result.(sens{isen});
        result.(sens{isen}) = [];
        result.(sens{isen}).date = result.datenum;
        result.(sens{isen}).data = data;
    end;

    % Global datenum field no longer needed
    result = rmfield(result, 'datenum');


    % Finally, remove any obviously INVALID values from each time series
    sens = fieldnames(result);
    for isen = 1:length(sens)
        if ( ~iscell(result.(sens{isen}).data) )
            rng = valid_var_range(sens{isen});
            badix = find( ~(rng(1) <= result.(sens{isen}).data & ...
                            result.(sens{isen}).data <= rng(2)) );
            % %DEBUG:            if ~isempty(badix); fprintf(2,'Flagging %d bad points in %s\n',length(badix),sens{isen}); end;
            % result.(sens{isen}).data(badix) = nan;
            %DEBUG:            if ~isempty(badix); fprintf(2,'Removing %d bad points in %s\n',length(badix),sens{isen}); end;
            result.(sens{isen}).date(badix) = [];
            result.(sens{isen}).data(badix) = [];
        end;
    end;


    % If caller wishes it, plot all numeric values
    if ( plotAllResults )
        disp(sprintf('Total cols %d, with data %d', length(hdr), nGoodCols));

        % No use in plotting date or time fields!
        sens(strcmpi(sens, 'date')) = [];
        sens(strcmpi(sens, 'year')) = [];
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

    % Merge new data with existing station data struct (if any)
    % Note this CAN overwrite data in 'stn' - see m-file help.
    stn = merge_station_data(stn, result);

    % If this file has sections of two or more different formats (i.e., has
    % more than one header line, i.e., has more than one sequence of 10-char
    % columns in it), then have this function call itself recursively to
    % handle each such 'frame' of the file. NOTE this may have a side effect
    % of producing TWO or more plots per sensor!
    if ( ~isempty(frames) )
      % startpos = startpos + ( frames(1) * length(hdr) * 10 ) + frames(1);
      startpos = startpos + (frames(1) * length_of_line);
      stn = load_10col_data(fname, stn, plotAllResults, startpos);
    end;

return;
