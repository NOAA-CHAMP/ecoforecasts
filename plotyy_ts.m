function [axen,lh1,lh2] = plotyy_ts(varargin)
%function [axen,lh1,lh2] = plotyy_ts([ax,]ts1_1,[lineopts1,...]ts2[,lineopts2,...][,plotfun1[,plotfun2]])
%
% PLOTYY_TS is similar to PLOTYY (v.) except that values to be plotted are
% not given as separate X and Y arguments, but are instead drawn from *time
% series structs* TS1,TS2 which are expected to have .date and .data fields.
% Also, PLOTYY_TS does NOT deal gracefully with plotting more than two TSes.
%
% First arg may optionally specify the AXES to plot in, as for PLOTYY. Any
% number of line/color/marker style strings or name,value pairs (v. PLOT) may
% be specified after each time series struct. Lastly, one or two function
% handles may be passed in to use in plotting each time series (v. PLOTYY).
%
% Calls DATETICK (DATETICK3 if present) to display X labels as date strings.
%
% SAMPLE CALL:
%  >> % Plot Molasses Reef air temperature time series on one AXES, thick blue line,
%  >> % Sombrero Key sea temperature on a separate AXES on the same figure panel.
%  >> plotyy_ts(mlrf1.ndbc_air_t,smkf1.ndbc_sea_t);
%
% Last Saved Time-stamp: <Wed 2018-02-21 12:46:30 EST lew.gramer>

  argix = 1;
  if ( ishandle(varargin{1}) )
    ax = varargin{argix};
    argix = argix + 1;
  else
    ax = gca;
  end;

  ts1 = varargin{argix};
  argix = argix + 1;

  ts1_args = {};
  while ( argix <= nargin && ~is_ts(varargin{argix}) )
    ts1_args{end+1} = varargin{argix};
    argix = argix + 1;
  end;

  ts2 = varargin{argix};
  argix = argix + 1;

  ts2_args = {};
  while ( argix <= nargin && ~isa(varargin{argix},'function_handle') )
    ts2_args{end+1} = varargin{argix};
    argix = argix + 1;
  end;

  if ( argix <= nargin )
    plotfun1 = varargin{argix};
    argix = argix + 1;
  else
    plotfun1 = @plot;
  end;
  if ( argix <= nargin )
    plotfun2 = varargin{argix};
    argix = argix + 1;
  else
    plotfun2 = @plot;
  end;
  
  [axen,lh1,lh2] = plotyy(ax,ts1.date,ts1.data,ts2.date,ts2.data,plotfun1,plotfun2);

  hold(axen(1),'on');
  hold(axen(2),'on');

  % If user specified line options, apply them to the line handles
  ts1_args = plotspec_check(ts1_args);
  ts2_args = plotspec_check(ts2_args);

  if ( ~isempty(ts1_args) )
    set(lh1,ts1_args{:});
    clrix = find(strcmpi(ts1_args,'Color'));
    % Special handling for line color - YAxis should match it
    if ( ~isempty(clrix) )
      set(axen(1),'YColor',ts1_args{clrix+1});
      set(axen(1),'GridColor',ts1_args{clrix+1});
    end;
  end;
  if ( ~isempty(ts2_args) )
    set(lh2,ts2_args{:});
    clrix = find(strcmpi(ts2_args,'Color'));
    % Special handling for line color - YAxis should match it
    if ( ~isempty(clrix) )
      set(axen(2),'YColor',ts2_args{clrix+1});
      set(axen(2),'GridColor',ts2_args{clrix+1});
    end;
  end;

  % Add datestrings to X labels
  if ( exist('datetick3','file') )
    datetick3;
  elseif ( exist('datetick2_5','file') )
    datetick2_5;
  elseif ( exist('datetick2','file') )
    datetick2;
  else
    datetick(axen(1),'x',2,'keeplimits','keepticks');
  end;
  if ( exist('set_datetick_cursor','file') )
    set_datetick_cursor;
  end;

return;
