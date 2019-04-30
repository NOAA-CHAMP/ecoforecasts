function lh = plot_ts(varargin)
%function lh = plot_ts([ax,]ts1,[ts1args,...],ts2,...)
%
% PLOT_TS is identical to PLOT (v.) except that values to be plotted are not
% specified as separate X and Y arguments, but are instead drawn from *time
% series structs* TS1,... which are expected to have .date and .data fields.
% First arg may specify the AXES to plot in, as for PLOT; otherwise, any
% argument which is NOT a time series struct (or optionally a *vector* of
% time series structs) is passed on as an arg to PLOT.
%
% Calls DATETICK (DATETICK3 if present) to display X labels as date strings.
%
% SAMPLE CALL #1:
%  >> % Plot Molasses Reef air temperature time series with thick blue line,
%  >> % Sombrero Key sea temperature with thin red dotted line & circle marker
%  >> plot_ts(mlrf1.ndbc_air_t,'LineWidth',2,smkf1.ndbc_sea_t,'ro:','LineWidth',.5);
%
% SAMPLE CALL #2:
%  >> % Plot MISST sea temperature in all stations in struct vector SITES;
%  >> % Rely on PLOT_TS to automatically distinguish by Color and Marker
%  >> plot_ts(sites.misst); legend({sites.station_name});
%
% Last Saved Time-stamp: <Wed 2018-07-11 17:01:12 Eastern Daylight Time gramer>

  argix = 1;
  if ( ishandle(varargin{1}) )
    ax = varargin{argix};
    argix = argix + 1;
  else
    ax = gca;
  end;

  co = get(ax,'ColorOrder');
  ms = '.ox+*sdv^<>ph';
  ncos = size(co,1);

  % Extract time series structs to plot, and plot args for each TS
  ts = {};
  pargs = {};
  nts = 0;
  while ( argix <= nargin )
    arg = varargin{argix};
    if ( isstruct(arg) && isfield(arg,'date') && isfield(arg,'data') )
      % Gracefully deal with vectors OR MATRICES of structs passed as a single argument
      for strcix = 1:numel(arg)
        strc = arg(strcix);
        if ( isempty(strc.date) || any(size(strc.date) ~= size(strc.date)) )
          error('Struct arguments must be non-empty time series!');
        end;
        nts = nts + 1;
        ts{nts} = strc;
        pargs{nts} = {};
      end;
    elseif ( nts == 0 )
      error('First arg after optional AXES handle must be a time series struct or vector of TS structs!');
    else
      pargs{nts} = {pargs{nts}{:} arg};
    end;
    argix = argix + 1;
  end;

  % Try to ensure that every plot has a distinct plot spec
  for tix = 1:nts
    specix = tix;
    if ( numel(pargs{tix}) > 1 && strcmpi(pargs{tix}(1),'specix') )
      specix = pargs{tix}{2};
      pargs{tix}(1:2) = [];
      if ( ~isnumeric(specix) || ~isscalar(specix) || specix < 1 || fix(specix) ~= specix )
        error('Invalid value for named arg SPECIX!');
      end;
    end;
    if ( isempty(pargs{tix}) )
      % Use default color order if user did not give a plot spec (see PLOT)
      % If user passed in >ncos Time Series, also distinguish them by Marker
      pargs{tix} = {pargs{tix}{:} 'Color',co(mod(specix-1,ncos)+1,:), 'Marker',ms(ceil(specix/ncos))};
    else
      % Plot spec, if given, must be the first plot arg for this TS
      arg = pargs{tix}{1};
      % The following mess is due to PLOT's ridiculously flexible calling sequence!
      if ( ~ischar(arg) || length(arg) > 4 )
        % User did not give a plot spec
        pargs{tix} = {'Color',co(mod(specix-1,ncos)+1,:) pargs{tix}{:}};
      elseif ( ~isempty(regexp(arg,'^[.ox+*sdv^<>ph:-]*[bgrcmykw][.ox+*sdv^<>ph:-]*$')) )
        % User gave a plot spec with a color in it - do nothing
      elseif ( ~isempty(regexp(arg,'^[.ox+*sdv^<>ph:-]*$')) )
        % User gave a plot spec WITHOUT a color in it
        pargs{tix} = {pargs{tix}{1} 'Color',co(mod(specix-1,ncos)+1,:) pargs{tix}{2:end}};
      else
        % User did not give a plot spec
        pargs{tix} = {'Color',co(mod(specix-1,ncos)+1,:) pargs{tix}{:}};
      end;
    end;
  end;

  % Turn on HOLD
  hold_state = ishold(ax);
  hold('on');
  % Plot each line
  %DEBUG:  disp(nts);
  for tix = 1:nts
    lh(tix) = plot(ax,ts{tix}.date,ts{tix}.data,pargs{tix}{:});
  end;
  % Restore previous HOLD state
  if ( ~hold_state )
    hold('off');
  end;

  % Add datestrings to X labels
  datetick('x',2,'keeplimits','keepticks');
  if ( exist('set_datetick_cursor','file') )
    set_datetick_cursor;
  end;

  % Work around bug in DATETICK2/DATETICK3 that shuffles current AXES
  axes(ax);

return;
