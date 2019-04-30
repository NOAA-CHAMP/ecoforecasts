function stn = load_factories(fname, stn)
%function stn = load_factories(fname, stn)
%
% Load all fact factories for station 'stn' from the CSV file 'fname'.
% ASSUMES that 'fname' is in the format output by the most recent version
% of the ICON/G2 procedure 'imn-dump-fact-factories'.
%
% New fields are added to the returned struct 'stn', one for each factory.
% Each field contains the 'fact-sensor', the name of the original fact source
% variable in ICON/G2, and a cell array of fuzzy value ranges.
%
% Note 'verify_variable' is also called to build all requisite time-series.
%
% Last Saved Time-stamp: <Mon 2009-03-23 16:31:38 Eastern Daylight Time gramer>
%

    % If we don't already have data, warn user of probable errors later!
    stnWasEmpty = false;
    if ( ~exist('stn', 'var') || isempty(stn) )
        stn = [];
        stnWasEmpty = true;
        warning('%s\n%s\n%s', ...
                'Second arg "stn" was empty! Some variables required ', ...
                'by factories are not present: call VERIFY_FACTORY_VARIABLES ', ...
                'later, once all data is loaded for this "stn" struct.');
    end;

    % Open and parse file
    fid = fopen(fname, 'r');
    if ( fid < 0 )
        error('Unable to open factories CSV file %s', fname);
    end

    vals = textscan(fid, '%[^,]%[^,]%[^,]%[^,]%[^\n]', 'Delimiter', ',');

    fclose(fid);


    % Turn all funky G2-ish characters into matlab-friendly '_'
    sensors = lower(regexprep(vals{1}, {'''','\W+'}, {'','_'}));
    variables = lower(regexprep(vals{2}, {'''','\W+'}, {'','_'}));
    stati = lower(regexprep(vals{3}, {'''','\W+'}, {'','_'}));
    validjdays = lower(regexprep(vals{4}, {'''','\W+'}, {'','_'}));
    bounds = lower(vals{5});


    % Create a new field in struct 'stn.factories' for each factory
    for idx = 1:length(sensors)

        sensor = sensors{idx};
        variable = variables{idx};
        status = stati{idx};
        validjday = validjdays{idx};
        bound = bounds{idx};

        if ( ~strcmpi(status,'online') )
          fprintf(1, 'Skipping factory %s,%s%s,%s,%s\n', ...
                  sensor, variable, status, validjday, bound);
          continue;
        end;

        % HACK! ICON/G2 var names ALWAYS start with a station code: strip it
        % off from everywhere in the variable name - including INTERNALLY!
        stanm = strtok(variable, '_');
        if ( ~isempty(stanm) )
            variable = regexprep(variable, [stanm '_'], '');
        end;

        %%%% NOTE: If 'stn.factories' did not exist, it is created here...
        stn.factories.(sensor).variable = variable;
        stn.factories.(sensor).status = status;

        % At this point, we have two numbers surrounded by underscores!
        validjday = regexprep(validjday, '(sequence)*_*([0-9]*)[ _]*([0-9]*)_*', '[$2,$3]');
        stn.factories.(sensor).validjdays = eval(validjday);

        % This ugliness converts all bounds CSV fields into a single cell
        % array, with one cell for each value range, containing as members
        % the string 'fuzzy-value', and a 2x1 vector of the lower and upper
        % bounds of that fuzzy value range.
        bound = regexprep(bound, '''([^=]+)=\(([^,]+)\)', '{''$1'',[$2]}');
        stn.factories.(sensor).fuzzies = eval(['{' bound '}']);

        % Be sure all required 'derived' variables are (or get) calculated
        if ( ~stnWasEmpty )
          if ( strcmpi(stn.factories.(sensor).status,'online') )
            stn = verify_variable(stn, variable);
          end;
        end;

    end;

return;
