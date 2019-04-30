function newstn = unique_10col_data(stn)
%function newstn = unique_10col_data(stn)
%
% Eliminate all duplicate timestamps and associated data. As of 2017, also
% removes dupes from .prof, .rawprof, and if it is a time series, .field.
%
% Last Saved Time-stamp: <Mon 2017-03-06 16:17:47 Eastern Standard Time gramer>

  newstn = stn;

  flds = fieldnames(stn);
  for ifld = 1:length(flds)

    fld = flds{ifld};

    if ( ~isfield(stn.(fld), 'date') || ~isfield(stn.(fld), 'data') )

      if ( ~ismember(fld,{'station_name','lon','lat','depth'}) )
        warning('Ecoforecasts:mergedNonTS','Passing field "%s" unchanged: no .date or .data', fld);
      end;

    else

      % Find all identical timestamps
      repidx = find(diff(newstn.(fld).date) <= 0);
      while ( ~isempty(repidx) )
        % Eliminate all duplicate timestamps and data
        newstn.(fld).date(repidx) = [];
        newstn.(fld).data(repidx) = [];

        if ( isfield(newstn.(fld),'prof') )
          newstn.(fld).prof(repidx,:) = [];
        end;
        if ( isfield(newstn.(fld),'rawprof') )
          newstn.(fld).rawprof(repidx,:) = [];
        end;
        if ( isfield(newstn.(fld),'field') )
          if ( ndims(newstn.(fld).field) ~= 3 )
            warning('Ecoforecasts:mergedNonTS',...
                    'Result %s.field was not a time series field?!', fld);
          else
            newstn.(fld).field(repidx,:,:) = [];
          end;
        end;

        repidx = find(diff(newstn.(fld).date) <= 0);
      end;

    end;

  end;

return;
