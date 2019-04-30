function stn = merge_station_data(stn, result)
%function stn = merge_station_data(stn, result)
%
% Merge all possible data from an existing station data struct STN, with a
% newly loaded set of station data RESULT. Merge done based on timestamps.
%
% NOTE WELL: If STN is not specified or is an empty array, a new struct is
% returned containing just the loaded data. If STN is a struct that already
% contains a field matching the name of a field in RESULT, an attempt will
% be made to merge that existing data with the newly loaded data. If for any
% reason an existing field does not contain both '.data' and '.date' fields,
% however, *that field is replaced* in the struct STN with the new data.
%
% ADDED Sep 2012: Now also attempts to merge 2D (profile) and 3D (e.g., SST)
% time series within existing STN struct.
%
% Last Saved Time-stamp: <Fri 2012-09-28 15:58:44 Eastern Daylight Time lew.gramer>


    % Make sure new data is well-formed
    if ( isempty(result) || ~isstruct(result) )
        warning('Second arg RESULT is not well-formed! No changes to STN.');

    % If input STN was empty, just return loaded data
    elseif ( isempty(stn) )
        stn = result;

    % Otherwise merge new data intelligently into STN using timestamps...
    % NOTE: Where no timestamp exists, *wipe out* the old field in STN!
    else

        % To attempt a merge, we MUST enforce uniqueness on all timestamps
        stn = unique_10col_data(stn);
        result = unique_10col_data(result);

        resflds = fieldnames(result);
        resflds = resflds(~strcmp(resflds, 'datenum'));
        for fi = 1:length(resflds)

            fld = resflds{fi};

            % If this is a new field, our job is easy
            if ( ~isfield(stn, fld) )
                stn.(fld) = result.(fld);

            % NOTE: If existing field has no timestamps, we *TRASH* it here!
            % (Also, this only trashes a field if its name happens to match a
            % loaded COLUMN HEADER. User-defined fields should be preserved.)
            elseif ( ~isfield(stn.(fld), 'date') || isempty(stn.(fld).date) || ...
                     ~isfield(stn.(fld), 'data') || isempty(stn.(fld).data) )
                if ( ~strcmp(fld,'station_name') && ~strcmp(fld,'lon') && ~strcmp(fld,'lat') && ~strcmp(fld,'depth') )
                  if ( isfield(stn.(fld), 'date') && isempty(stn.(fld).date) && ...
                       isfield(stn.(fld), 'data') && isempty(stn.(fld).data) )
                    warning('Field %s was empty: Replacing.', fld);
                  else
                    warning('Field %s cannot be merged: REPLACING old value!', fld);
                  end;
                  stn = rmfield(stn, fld);
                  stn.(fld) = result.(fld);
                end;

            % Otherwise, merge all values - overwriting old data with new
            else
                % First do the delicate merge of datenums
                % Use a ~9-second (1e-4 day) tolerance for matching datenums
                % (Narrow window ensures we only remove LIKELY duplicates.)
                newdatenum = union(roundn(stn.(fld).date,-4), roundn(result.(fld).date,-4));
                newdatenum = union(stn.(fld).date, result.(fld).date);
                si = find(ismember(newdatenum, stn.(fld).date));
                ri = find(ismember(newdatenum, result.(fld).date));
                stn.(fld).date = newdatenum;
                stn.(fld).data(si,1) = stn.(fld).data;
                % If the input field was a cell array, keep it that way
                if ( iscell(stn.(fld).data) == iscell(result.(fld).data) )
                    stn.(fld).data(ri,1) = result.(fld).data;
                elseif ( iscell(stn.(fld).data) )
                    stn.(fld).data(ri,1) = cellstr(num2str(result.(fld).data(:)));
                else
                    stn.(fld).data(ri,1) = str2double(result.(fld).data);
                end;
                % Frequently triggered by loading .date as row-, and .data as
                % column-vectors, for example, or vice versa. Stupid MATLAB. 
                if ( any(size(stn.(fld).date) ~= size(stn.(fld).data)) )
                    warning('Fields %s.date and .data do not match!', fld);
                    size(stn.(fld).date), size(stn.(fld).data),
                end;

                % For two- and three-dimensional time series
                if ( isfield(result.(fld), 'prof') )
                  n = size(result.(fld).prof,2);
                  if ( ~isfield(stn.(fld), 'prof') || n ~= size(stn.(fld).prof,2) )
                    error('Profile shape mismatch!');
                  end;
                  stn.(fld).prof(si,1:n) = stn.(fld).prof;
                  stn.(fld).prof(ri,1:n) = result.(fld).prof;
                end;
                if ( isfield(result.(fld), 'rawprof') )
                  n = size(result.(fld).rawprof,2);
                  if ( ~isfield(stn.(fld), 'rawprof') || n ~= size(stn.(fld).rawprof,2) )
                    error('Raw Profile shape mismatch!');
                  end;
                  stn.(fld).rawprof(si,1:n) = stn.(fld).rawprof;
                  stn.(fld).rawprof(ri,1:n) = result.(fld).rawprof;
                end;
                if ( isfield(result.(fld), 'field') )
                  if ( ndims(result.(fld).field) ~= 3 )
                    warning('Ecoforecasts:mergedNonTSField',...
                            'Result %s.field was not a time series field?!', fld);
                  else
                    n = size(result.(fld).field,2);
                    m = size(result.(fld).field,3);
                    if ( ~isfield(stn.(fld), 'field') || ndims(stn.(fld).field) ~= 3 ...
                         || n ~= size(stn.(fld).field,2) || m ~= size(stn.(fld).field,3) )
                      error('Field shape mismatch!');
                    end;
                    if ( (isfield(stn.(fld), 'lon') && any(result.(fld).lon ~= stn.(fld).lon)) ...
                         || (isfield(stn.(fld), 'lat') && any(result.(fld).lat ~= stn.(fld).lat)) )
                      error('Field LON or LAT mismatch!');
                    end;
                    stn.(fld).field(si,1:n,1:m) = stn.(fld).field;
                    stn.(fld).field(ri,1:n,1:m) = result.(fld).field;
                  end;
                end;

            end;

        end;

    end;


return;
