ERROR
function result = parse_10col_data(hdr, vals, stn)

    % Parse header line first, column by column
    line = fgetl(fid);
    idx = 0;
    dateidx = 0; houridx = 0; minuteidx = 0;
    yearidx = 0; jdayidx = 0;
    for curs = 1:10:length(line)
        idx = idx + 1;
        hdr{idx} = strtrim(line(curs:curs+9));
        % Try to make sure we have what looks like a valid header line
        if ( ~isnan(str2double(hdr{idx})) )
            error('Invalid header line?? Field %d is the number "%s"!', ...
                  idx, hdr{idx});
        end;

        % Try to 'translate' label into the equivalent ICON/G2 variable name:
        % 'WaterT' becomes 'sea_t', 'Sea1m' becomes 'ctd_surf_seatemp', etc.
        hdr{idx} = translate_var_name(hdr{idx});

        % Take care of any funky chars in header label
        hdr{idx} = regexprep(hdr{idx}, '\W+', '_');
        if ( strcmpi(hdr{idx}, 'date') ); dateidx = idx; end;
        if ( strcmpi(hdr{idx}, 'yr') ); yearidx = idx; end;
        if ( strcmpi(hdr{idx}, 'year') ); yearidx = idx; end;
        if ( strcmpi(hdr{idx}, 'jday') ); jdayidx = idx; end;
        if ( strcmpi(hdr{idx}, 'hr') ); houridx = idx; end;
        if ( strcmpi(hdr{idx}, 'hour') ); houridx = idx; end;
        if ( strcmpi(hdr{idx}, 'mn') ); minuteidx = idx; end;
        if ( strcmpi(hdr{idx}, 'min') ); minuteidx = idx; end;
        if ( strcmpi(hdr{idx}, 'minute') ); minuteidx = idx; end;
    end

    if ( dateidx == 0 && (yearidx == 0 || jdayidx == 0) )
        error('Need one of "date", or "year" and "jday" in file "%s"??', fname);
    end;
    if ( houridx == 0 )
        error('No "hour" column found in file "%s"??', fname);
    end;

    % Construct return structure from parsed data

    % Ensure 'datenum' is always the first field in result
    result.datenum = [];

    nGoodCols = 0;
    for idx = 1:length(hdr)
        % Only process non-blank columns
        if ( length(hdr{idx}) > 0 )
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


    % Calculate a Matlab datenum for each line of the file
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
            goodidx = find( rng(1) <= result.(sens{isen}).data & ...
                            result.(sens{isen}).data <= rng(2) );
            result.(sens{isen}).date = result.(sens{isen}).date(goodidx);
            result.(sens{isen}).data = result.(sens{isen}).data(goodidx);
        end;
    end;

