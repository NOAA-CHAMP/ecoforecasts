function [hlines,haxes,hfig] = multiplot_ts(ttl,varargin)
%function [hlines,haxes,hfig] = multiplot_ts(ttl,[stn1,field1,field2,...,stn2,...])
%
% Plot time series from multiple stations/sources in one figure, as stacked
% plots with identical X-axis - v. functions MULTIPLOT, MULTIPLOT_DATETICK.
%
% EXAMPLE
% >> fwyf1 = load_station_data('fwyf1');
% >> mlrf1 = load_station_data('mlrf1');
% >> multiplot_ts('Fowey v. Molasses sea/air temps',fwyf1,'air_t','sea_t',mlrf1,'air_t','sea_t');
%
% Compare to the single-station multiplot function:
%  [stn,hlines,haxes,hfig] = MULTIPLOT_STATION(stn,fldnms,ttl,xlbl,ylbl,xlm,doverify,linspec)
%
% Note that MULTIPLOT_TS does *not* try to verify any station variables, but
% just produces an error if they don't already exist. (Cf. VERIFY_VARIABLE.)
%
% Last Saved Time-stamp: <Wed 2013-07-24 00:20:27 Eastern Daylight Time gramer>


  % % Get all name-value properties from the top of the stack
  % args = varargin;
  % if ( ischar(args{1}) )
  %   plotargs = args(1:n);
  %   args = args(n+1:end);
  % end;

  if ( ~ischar(ttl) )
    error('First arg should be a plot title!');
  end;

  X = {};
  Y = {};
  xlbl = 'Date';
  ylbl = {};
  xlm = [];
  ylms = [];
  datetick_fmt = [];
  linspec = {'b.-','g.-','r.-','c.-','m.-','y.-','k.-'};

  curstn = [];
  stno = 0;

  for ix = 1:length(varargin)

    arg = varargin{ix};

    if ( isstruct(arg) )
      curstn = arg;
      stno = stno + 1;

      if ( length(inputname(ix+1)) )
        curstnm = inputname(ix+1);
      elseif ( isfield(curstn,'station_name') )
        curstnm = curstn.station_name;
      elseif ( isfield(curstn,'name') )
        curstnm = curstn.name;
      elseif ( isfield(curstn,'name') )
        curstnm = curstn.name;
      elseif ( isfield(curstn,'station_code') )
        curstnm = curstn.station_code;
      elseif ( isfield(curstn,'code') )
        curstnm = curstn.code;
      else
        curstnm = sprintf('STN%02d', stno);
      end;

    elseif ( iscellstr(arg) )
      for ix = 1:length(arg)
        [X,Y,ylbl] = multiplot_ts_get_arg(X,Y,ylbl,curstn,curstnm,arg{ix});
      end;

    elseif ( ischar(arg) )
      [X,Y,ylbl] = multiplot_ts_get_arg(X,Y,ylbl,curstn,curstnm,arg);

    else
      warning('Unknown type for arg %d!', ix);

    end;

  end;

  if ( isempty(X) )
    error('No valid time series specified to plot!');
  end;

  ylbl = strrep(ylbl, '_', '\_');

  [hlines,haxes,hfig] = multiplot_datetick(X,Y,ttl,xlbl,ylbl,xlm,ylms,datetick_fmt,linspec);

return;


function [X,Y,ylbl] = multiplot_ts_get_arg(X,Y,ylbl,curstn,curstnm,arg);
  if ( ~isfield(curstn,arg) )
    warning('No field "%s" in struct "%s"!', arg, curstnm);
  else
    X = { X{:} curstn.(arg).date };
    Y = { Y{:} curstn.(arg).data };
    ylbl = { ylbl{:} sprintf('%s.%s', curstnm, arg) };
  end;
return;
