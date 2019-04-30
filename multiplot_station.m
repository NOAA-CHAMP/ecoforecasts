function [stn,hlines,haxes,hfig]=multiplot_station(stn,fldnms,ttl,xlbl,ylbl,xlm,ylms,doverify,linspec)
%function [stn,hlines,haxes,hfig]=multiplot_station(stn,fldnms,ttl,xlbl,ylbl,xlm,ylms,doverify,linspec)
%
% Use Matlab Exchange function 'MULTIPLOT' to plot multiple time series of
% data from the in situ station data structure 'stn', as named in string or
% or cellstr 'fldnms'.  Each such data field in turn must be a struct with
% its own '.date' and '.data' sub-fields - see, e.g., 'LOAD_STATION_DATA'.
% Optional arg YLMS is either a 1xN cell array of 1x2 vectors, or a single
% numeric 1x2 vector, that specifies y-axis limits for each of the N plots.
% If optional DOVERIFY is true (DEFAULT), try to use VERIFY_VARIABLE to build
% each field in FLDNMS. Optional cell or char array arg LINSPEC is passed to
% MULTIPLOT_DATETICK; DEFAULT is as specified by that function (qv.).
%
% Last Saved Time-stamp: <Fri 2012-07-13 19:45:04  lew.gramer>

  if ( ~iscell(fldnms) )
    fldnms = { fldnms };
  end;
  if (nargin < 3 || isempty(ttl))
    ttl = 'Station Data';
    if ( isfield(stn, 'station_name') )
      ttl = [ upper(stn.station_name) ' - ' ttl ];
    end;
  end;
  if (nargin < 4 || (~ischar(xlbl) && isempty(xlbl))); xlbl = 'Date-Time'; end;
  if (nargin < 5 || isempty(ylbl)); ylbl = strrep(lower(fldnms),'_','\_'); end;
  if (nargin < 6); xlm = []; end;
  if (nargin < 7); ylms = []; end;
  if (nargin < 8); doverify = true; end;
  if (nargin < 9); linspec = []; end;

  X = {};
  Y = {};
  ylbls = {};
  max_ndates = 0;
  for ix = 1:length(fldnms)
    fldnm = fldnms{ix};
    if ( doverify )
      stn = verify_variable(stn, fldnm);
    end;
    if ( ~isfield(stn, fldnm) || ~isfield(stn.(fldnm),'date') || isempty(stn.(fldnm).date) )
      warning('No valid field "%s" in station struct!', fldnm);
    else
      X{end+1} = stn.(fldnm).date;
      Y{end+1} = stn.(fldnm).data;
      ylbls(end+1) = ylbl(ix);
      max_ndates = max(max_ndates,length(stn.(fldnm).date));
    end;
  end;

  if ( isempty(X) )
    error('No valid fields specified to plot!');
  end;

  datetick_fmt = 2;
  % if ( max_ndates > 1500 )
  %   datetick_fmt = 12;
  % end;

  [hlines, haxes, hfig] = multiplot_datetick(X,Y,ttl,xlbl,ylbls,xlm,ylms,datetick_fmt,linspec);
  set(hfig, 'Name', ttl);

  % HACK: MATLAB 2007a puts the (often huge) station struct in ANS no matter what
  if ( nargout < 1 )
    stn = [];
    clear stn;
  end;

return;
