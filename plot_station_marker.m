function h = plot_station_marker(stn,faceclr,lineclr,fontsz,varargin)
%function h = plot_station_marker(stn,faceclr,lineclr,fontsz[,NAME1,VAL1,...])
%
% PLOT very obvious pentagon marker and if possible, TEXT label for station
% STN. STN may be a station name STRING (v. GET_STATION_FROM_STATION_NAME),
% or a STRUCT with fields .lon,.lat, and optional .station_name, or a cell
% array with two (lon,lat) or three elements (... station_name), or a numeric
% vector with at least two (lon,lat) elements. NOTE: If STN is a cell array
% containing vectors, multiple station markers will be plotted.
%
% If optional FACECLR False or empty, pentagram is only plotted in outline;
% otherwise, plot inner filled pentagram of color FACECLR (DEFAULT: 'w'). If 
% optional LINECLR False or empty (DEFAULT: 'r'), do not plot outline. If
% optional FONTSZ (DEFAULT: 12) evaluates to 0, do not TEXT station_name.
%
% All remaining arguments (NAME,VALUE pairs) are passed to TEXT (*not* PLOT).
%
% Optionally returns a vector of graphics handles from PLOT and TEXT.
%
% Last Saved Time-stamp: <Thu 2018-08-16 14:17:27 Eastern Daylight Time gramer>

  lon = [];
  lat = [];
  dep = 0;
  stnm = [];
  if ( ischar(stn) || isstring(stn) )
    stnm = char(stn);
    [lon,lat,dep] = get_station_coords(stnm);
  elseif ( iscell(stn) && numel(stn) >= 2 )
    lon = stn{1};
    lat = stn{2};
    if ( numel(stn) >= 3 && isnumeric(stn{3}) ); dep = stn{3}; end;
    if ( numel(stn) >= 3 && ischar(stn{end}) ); stnm = stn{end}; end;
  elseif ( isnumeric(stn) && numel(stn) >= 2 )
    lon = stn(1);
    lat = stn(2);
    if ( numel(stn) >= 3 ); dep = stn(3); end;
  elseif ( isstruct(stn) && isfield(stn,'lon') && isfield(stn,'lat') )
    lon = stn.lon;
    lat = stn.lat;
    if ( isfield(stn,'station_name') )
      stnm = stn.station_name;
    end;
    if ( isfield(stn,'depth') )
      dep = stn.depth;
    end;
  else
    error('Invalid first argument');
  end;
  % If we got vectors, make sure their number of elements match
  if ( numel(lon) ~= numel(lat) || (iscellstr(stnm) && numel(lon) ~= numel(stnm)) )
    error('If STN specifies arrays for LON, LAT, or station name, sizes must match');
  end;
  if ( isscalar(dep) )
    dep = repmat(dep,size(lon));
  end;

  if ( ~exist('lineclr','var') || (islogical(lineclr) && lineclr) )
    lineclr = 'r';
  end;
  if ( ~exist('faceclr','var') || (islogical(faceclr) && faceclr) )
    faceclr = 'w';
  end;
  if ( ~exist('fontsz','var') || (islogical(fontsz) && fontsz) )
    fontsz = 12;
  end;

  h=[];
  for stix = 1:numel(lon);
    alon = lon(stix);
    alat = lat(stix);
    adep = dep(stix);
    if ( faceclr )
      h(end+1) = plot3(alon,alat,adep,'wp','MarkerSize',10,'MarkerFaceColor',faceclr);
    end;
    if ( lineclr )
      h(end+1) = plot3(alon,alat,adep,'p','Color',lineclr,'MarkerSize',12);
    end;
    if ( fontsz > 0 && length(stnm) > 0 )
      if (iscellstr(stnm))	astnm = stnm{stix};
      else			astnm = stnm; end;
      if ( exist('textize') ); astnm = textize(astnm); end;
      h(end+1) = text(alon,alat,adep,...
                      [' \Leftarrow ',upper(astnm)],...
                      'Color','k','FontSize',fontsz,varargin{:});
    end;
  end;

return;
