function [annual, amj, jas, ond, jfm] = seasonalize_icon_data(filenamepats, fieldname, validrange, validhours)
%FUNCTION [annual, amj, jas, ond, jfm] = seasonalize_icon_data(filenamepats, fieldname, validrange, validhours)
%
% Read a series of XLS Excel files containing ICON/CMAN station data,
% including Julian days (REQUIRED), extract the column of values for
% a particular named property from all of those XLS files (e.g., sea
% temperature at 1m, salinity at 3m, PAR at surface, whatever), and
% then categorize that range of values into four seasonal groups...
%
% INPUTS:
%  filenamepats - cell array of XLS file name pattern strings
%  fieldname - name of property data column in each XLS file
%  validrange - [ minimum-acceptable-value maximum-acceptable-value ]
%      DEFAULT: All property values are acceptable
%  validhours - [ earliest-hour-of-day-to-consider latest-hour-to-consider ]
%      DEFAULT: All hours are acceptable
%
% OUTPUTS:
%  See 'help derive_prop_ranges' for the layout of each of the five return
%  structures - annual, and one each for Apr-May-Jun (amj), Jul-Aug-Sep (jas),
%  Oct-Nov-Dec (ond) and Jan-Feb-Mar (jfm), resp.
%
% DEPENDENCIES:
%  extract_icon_fields.m - load multiple columns from multiple Excel files.
%  derive_prop_ranges.m - CLIPS-specific statistical treatment of a dataset.

    if ( ~iscell(filenamepats) )
        filenamepats = {filenamepats};
    end
    if ( ~exist('validrange', 'var') )
        validrange = [-Inf +Inf];
    end
    if ( ~exist('validhours', 'var') )
        validhours = [0 24];
    end

    annual = [];
    amj = [];
    jas = [];
    ond = [];
    jfm = [];


    allfields = { [] ; [] ; [] ; [] };

    % Open each XLS file matching any element of cellstr 'filenamepats'...
    for idx = 1:length(filenamepats)

        filenamepat = filenamepats{idx};

        % ??? Later we need to do a 'dir()' here for pattern matching!
        filename = filenamepat;

        % Extract date, hour, plus the requested column of property values
        if ( strcmp(lower(fieldname), 'sst') )
          fields = extract_sst_fields(filename, {'date', 'jday', 'hour', fieldname});
        else

          fields = extract_icon_fields(filename, {'date', 'jday', 'hour', fieldname});

          % Try any known variations - otherwise, skip this particular file!
          if ( length(fields) ~= 4 )
            fields = extract_icon_fields(filename, {'year', 'julian day', 'hour', fieldname});
          end
        end


        if ( length(fields) == 4 )

          if ( iscell(fields{2}) || iscell(fields{3}) || iscell(fields{4}) )
            fprintf('BAD DATA ERROR! File "%s"\n', filename);
            keyboard;
          end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%keyboard;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

          % Build up a cell array of vectors containing ALL requested values
          allfields(1) = { [ allfields{1} ; fields{1} ] };
          allfields(2) = { [ allfields{2} ; fields{2} ] };
          allfields(3) = { [ allfields{3} ; fields{3} ] };
          allfields(4) = { [ allfields{4} ; fields{4} ] };

        end

    end

    % If we didn't get what we were after, then give up!
    if ( isempty(allfields{2}) || isempty(allfields{3}) || isempty(allfields{4}) )
      return;
    end


    % Extract only those rows from the times of day we want to consider...
    goodhours = (validhours(1) <= allfields{3} & allfields{3} <= validhours(2));
    % And extract only those rows which contain a "valid" property value...
    goodvals = (validrange(1) <= allfields{4} & allfields{4} <= validrange(2));

    goodidx = goodhours & goodvals;
    arr = allfields{1};
    fields(1) = { arr(goodidx) };
    arr = allfields{2};
    fields(2) = { arr(goodidx) };
    arr = allfields{3};
    fields(3) = { arr(goodidx) };
    arr = allfields{4};
    fields(4) = { arr(goodidx) };

    % Build an array of daily averages for the requested field
    % (??? Skip for now - CLIPS code doesn't check values per-time either!)


    %
    % Categorize this property's values into annual and seasonal
    % ranges for "high", "very high", "somewhat low", etc., etc.
    %

    annual = derive_prop_ranges([fields{4}(:)]);

    % NOTE: In ORIGINAL CLIPS code, seasons started on julian days as
    % follows: spring=60, summer=120, fall=240(!), winter=330... Per
    % Jim's request this is now just based on equinoces and solstices!
    % And the season names were changed to apply to both hemispheres.

    jasidx = find(172 < fields{2} & fields{2} <= 264);
    jas = derive_prop_ranges([fields{4}(jasidx)]);

    ondidx = find(264 < fields{2} & fields{2} <= 355);
    ond = derive_prop_ranges([fields{4}(ondidx)]);

    jfmidx = find(355 < fields{2} | fields{2} <= 80);
    jfm = derive_prop_ranges([fields{4}(jfmidx)]);

    amjidx = find(80 < fields{2} & fields{2} <= 172);
    amj = derive_prop_ranges([fields{4}(amjidx)]);

return;
