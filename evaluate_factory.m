function facts = evaluate_factory(stn, facname, varname, begdate, enddate)
%function facts = evaluate_factory(stn, facname, varname, begdate, enddate)
%
% Evaluate a time series and produce a struct containing a cell array of
% "fuzzy-value" strings for each 3-hour average in that time series.
% First arg 'stn' should be a station data struct already populated by calls
% to load_10col_data, load_factories, and similar m-funcs. Second arg
% 'facname' is a string - the name of the factory to evaluate, which should
% also be the name of a field in the stn.factories struct. Third arg
% 'varname' is an (optional) replacement for the name of the time series
% variable to evaluate, which should also be a field in 'stn'.
% If both 'begdate' and 'enddate' form a VALID daterange, only values between
% those two dates are evaluated and returned as facts. (And if arg 'varname'
% is empty [], the default variable name from the factory is still used.)
%
% Last Saved Time-stamp: <Sun 2010-09-26 16:01:31 Eastern Daylight Time gramer>

% MUCH SIMPLER way to evaluate factories in the near future:
% load hospital
% edges = 0:10:100;
% labels = strcat(num2str((0:10:90)','%d'),{'s'});
% AgeGroup = ordinal(hospital.Age,labels,[],edges);
% [c,labels] = summary(AgeGroup);
%
% Table = dataset({labels,'AgeGroup'},{c,'Count'});
% Table(3:6,:)
% ans = 
%     AgeGroup     Count
%     '20s'        15   
%     '30s'        41   
%     '40s'        42   
%     '50s'         2  

    facts = [];

    if ( ~isstruct(stn) || ~isfield(stn, 'factories') )
        error('First argument "stn" is not a valid struct!');
    end;
    if ( ~exist('facname', 'var') || ~isfield(stn.factories, facname) )
        error('Second argument is not a factory name in "stn"!');
    end;
    % If caller did not specify a replacement variable name, use default
    if ( ~exist('varname', 'var') || isempty(varname) )
        varname = stn.factories.(facname).variable;
    end;

    if ( ~isfield(stn, varname) )
        error('Variable name %s is not a field in "stn"!', varname);
    end;

    if ( ~exist('begdate', 'var') || ~exist('enddate', 'var') || ...
         isempty(begdate) || isempty(enddate) || begdate >= enddate )
        begdate = [];
        enddate = [];
    end;


    % Calculate 3-hour averages over whole variable history (or a subset of
    % the history, if 'begdate' and 'enddate' evaluate to a valid daterange)
    idx = 1:length(stn.(varname).date);
    if ( ~isempty(begdate) && ~isempty(enddate) )
        idx = find(begdate <= stn.(varname).date & stn.(varname).date <= enddate);
    end;
    if ( isempty(idx) )
        warning('No data for "%s" found in range %s to %s!', ...
                varname, datestr(begdate), datestr(enddate));
    end;

    [ facts.datenum, facts.averages ] = ...
        window_func(stn.(varname).date(idx), stn.(varname).data(idx), 'avg', 3);

    fuzzies = stn.factories.(facname).fuzzies;

    

    % Assign a fuzzy-value for each 3-hour average in the time series
    fz = fuzzies{1};
    hits = find(facts.averages < fz{2}(1));
    facts.fuzzies([hits],1) = {'unbelievably-low'};
    %fuzzidx([hits]) = 0;

    for idx = 1:length(fuzzies)

        fz = fuzzies{idx};

        % Check range consistency as we go
        if ( ( fz{2}(1) > fz{2}(2) ) || ...
             ( 1 < idx && fz{2}(1) ~= fuzzies{idx-1}{2}(2) ) || ...
             ( idx < length(fuzzies) && fz{2}(2) ~= fuzzies{idx+1}{2}(1) ) )
          warning('Factory "%s": fuzzy range %d inconsistent!', ...
                  facname, idx);
        end;

        % NOTE: "<= fz{2}(2)" finds values that really belong in next-higher
        % range. This is OK as we'll overwrite these next time through loop. 
        hits = find(fz{2}(1) <= facts.averages & facts.averages <= fz{2}(2));
        facts.fuzzies([hits],1) = fz(1);
        %fuzzidx([hits]) = idx;

    end;

    fz = fuzzies{end};
    % NOTE: Values exactly equal to upper-most range are considered VALID
    hits = find(fz{2}(2) < facts.averages);
    facts.fuzzies([hits],1) = {'unbelievably-high'};
    %fuzzidx([hits]) = length(fuzzies)+1;

    %facts.fuzzies(isempty(facts.fuzzies)) = {'unknown'};

    % Verify the reasonableness of our result
    if ( length(facts.fuzzies) ~= length(facts.averages) )
        warning('Returned fewer fuzzies than average values??');
    end;


return;
