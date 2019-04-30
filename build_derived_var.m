function [stn, built] = build_derived_var(stn, varname)
%function [stn, built] = build_derived_var(stn, varname)
%
% Make sure that the given station data struct 'stn' contains a field named
% 'varname'. If 'varname' does not exist, build it based on parsing of its
% name: we assume we will find somewhere in 'varname' the name of an actual
% sensor time series, as well as a time-period and operator to apply to that
% time series. E.g., if varname is 'wind1_u_1_day_average', we will first
% make sure there exist fields in stn named 'wind1_dir' and 'wind1_speed'; we
% then calculate 'wind1_u' from these; and 'wind1_u_1_day_average' from this.
%
% As of 2018 Feb 25, will also build derived time series for each colummn in
% the STN.(varname).prof field, if found.
%
% Last Saved Time-stamp: <Sun 2018-02-25 14:36:40 Eastern Standard Time gramer>

  built = 0;

  if ( isfield(stn, varname) )
    if ( isfield(stn.(varname), 'date') && isfield(stn.(varname), 'data') )
      return;
    elseif ( isempty(stn.(varname)) )
      warning('BuildDerivedVar:EmptyField',...
              'Variable "%s" exists but is empty! Replacing...', varname);
    else
      error('BuildDerivedVar:BadField',...
            'Variable "%s" already exists but is malformed??', varname);
    end;
  end;

  flds = fieldnames(stn);

  % "Derived vector" variables 'u' and 'v' are special cases!
  % First make sure all of these already exist, or get built:
  findidx = regexp(varname, '(ndbc_wind|ncep_wind|cfsr_wind|erai_wind|wind|wxt)[0-9]*(_[^_]*)_*([uv]$|[uv]_)');
  for idx = 1:length(findidx)
    built = built + 1;
    begidx = findidx(idx);
    endidx = regexp(varname(begidx:end), '_([uv]$|[uv]_)');
    endidx = begidx + endidx(1) - 1;

    instfld = varname(begidx:endidx-1);
    newfld = varname(begidx:endidx+1);

    % Don't rebuild if we already have everything in place!
    if ( ~ismember(newfld, flds) )
      if ( varname(endidx+1) == 'u' )
        stn = vector_func(stn, instfld, 'u');
        flds = fieldnames(stn);
      else
        stn = vector_func(stn, instfld, 'v');
        flds = fieldnames(stn);
      end;
      % NOTE: VECTOR_FUNC calls FILTER_GAPS on its own...
    end;
  end;


  % Break the variable name into underscore-separated fragments
  uscrs = [1 strfind(varname, '_') length(varname)];

  % Count down underscores from the end of the name back, and for each
  % fragment, search 'stn' struct's field names for a match
  findx = [];
  for idx = length(uscrs):-1:2

    findx = find(strcmp(flds, varname(1:uscrs(idx)-1)));
    % When we find a var that already exists, start building from there,
    % e.g., we find 'sea_t', so now we could build 'sea_t_1_day_deviation',
    % and then 'sea_t_1_day_deviation_3_day_average', etc., ad infinitum.
    if ( ~isempty(findx) )

      fld = flds{findx(1)};

      if ( ~isfield(stn.(fld), 'date') || ~isfield(stn.(fld), 'data') ...
           || ~isnumeric(stn.(fld).data) || ~all(size(stn.(fld).date) == size(stn.(fld).data)) )
        error('Cannot build "%s": Base var "%s" invalid!', varname, fld);
      end;

      doProf = false;
      if ( isfield(stn.(fld),'prof') && size(stn.(fld).prof,1) == numel(stn.(fld).date) )
        doProf = true;
      end; %if isfield prof

      if ( strncmp('anom', varname(uscrs(idx)+1:end), 4) )

        newfld = [fld '_anom'];
        %disp(['Building ' newfld]);

        if ( isempty(stn.(fld).date) )
          warning('Base var "%s" empty! "%s" will be also...', fld, newfld);
          stn.(newfld).date = [];
          stn.(newfld).data = [];
          if doProf; stn.(newfld).prof = []; end;
        else
          stn.(newfld).date = stn.(fld).date;
          stn.(newfld).data = stn.(fld).data - nanmean(stn.(fld).data);
          if doProf; stn.(newfld).prof = stn.(fld).prof - nanmean(stn.(fld).prof,1); end;
        end; %if ( isempty(stn.(fld).date) )

        built = built + 1;
        % Recurse - or do it again a different way. :)
        if ( idx > 2 )
          [stn, inbuilt] = build_derived_var(stn, varname);
          built = built + inbuilt;
        end; %if ( idx > 2 )

      elseif ( strncmp('qc', varname(uscrs(idx)+1:end), 2) )

        % ??? Will need to handle special cases of 'u_qc' and 'v_qc' here...

        if ( isempty(stn.(fld).date) )
          newfld = [fld '_qc'];
          warning('Base var "%s" empty! "%s" will be also...', fld, newfld);
          stn.(newfld).date = [];
          stn.(newfld).data = [];
          if doProf; stn.(newfld).prof = []; end;
        else
          stn = qa_ts(stn, fld);
        end; %if ( isempty(stn.(fld).date) ) else

        built = built + 1;
        % Recurse - one imprecation may just not be enough :)
        if ( idx > 2 )
          [stn, inbuilt] = build_derived_var(stn, varname);
          built = built + inbuilt;
        end; %if ( idx > 2 )


      else

        [n, per, wnd, fun, rest] = ...
            get_derived_formula(varname(uscrs(idx)+1:end));
        if ( ~isempty(wnd) )

          newfld = sprintf('%s_%d_%s_%s', fld, n, per, fun);
          %disp(['Building ' newfld]);

          if ( isempty(stn.(fld).date) )
            warning('Base var "%s" empty! "%s" will be also...', fld, newfld);
            stn.(newfld).date = [];
            stn.(newfld).data = [];
            if doProf; stn.(newfld).prof = []; end;

          else

            if ( ~doProf )
              [stn.(newfld).date, stn.(newfld).data] = ...
                  window_func(stn.(fld).date, stn.(fld).data, fun, wnd, 1);
            else
              [stn.(newfld).date, stn.(newfld).data, stn.(newfld).prof] = ...
                  window_func_prof(stn.(fld).date, stn.(fld).data, stn.(fld).prof, fun, wnd, 1);
            end;

            % Ensure that we don't derive nonsense from gaps in the raw data
            % BUT raw data may be daily, triad, pentad, weekly, or even monthly
            maxgap = max(3,min(diff(unique(stn.(fld).date))));
            maxgap = max(maxgap,(wnd/2/24));
            if ( ~doProf )
              stn = filter_gaps(stn, fld, newfld, maxgap, (wnd/24));
            else
              stn = filter_gaps_prof(stn, fld, newfld, maxgap, (wnd/24));
            end;
          end; %if ( isempty(stn.(fld).date) ) else

          built = built + 1;
          % Recurse - someone may have dispelled it the first time :)
          if ( length(rest) > 1 )
            [stn, inbuilt] = build_derived_var(stn, varname);
            built = built + inbuilt;
          end;

        end; %if ( ~isempty(wnd) )

      end; %if ( strncmp('anom', varname(uscrs(idx)+1:end), 4) ) else

      break;

    end; %if ( ~isempty(findx) )

  end; %for idx = length(uscrs):-1:2

return;


function [n, per, wnd, fun, rest] = get_derived_formula(varname_frag)
%INTERNAL function [n, per, wnd, fun, rest] = get_derived_formula(varname_frag)
%
% Scan the input string for a time-value, time-unit, time series operand, and
% a base variable name: results are meant to be passed to window_func (qv).

  flds = textscan(varname_frag, '%f_%[^_]_%[^_]%s');
  n = flds{1};
  per = char(flds{2});
  fun = char(flds{3});
  rest = [];
  if ( ~isempty(flds{4}) )
    rest = char(flds{4});
  end;
  wnd = [];
  if ( ~isempty(n) && isnumeric(n) && n > 0 )
    switch ( per ),
     % ICON/G2 does not allow time units of months, years, etc.
     case {'w', 'week', 'weeks'},
      wnd = n * 24.0 * 7.0;
     case {'d', 'day', 'days'},
      wnd = n * 24.0;
     case {'h', 'hour', 'hours'},
      wnd = n;
     % Shorter time units are meaningless for hourly data
     otherwise,
      wnd = [];
      warning('Unrecognized time unit %d %s', n, per);
    end;
  end;

return;
