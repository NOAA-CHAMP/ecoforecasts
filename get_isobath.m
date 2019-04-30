function coords = get_isobath(bath,deps,mindep)
%function coords = get_isobath(bath,deps[,mindep])
%
% Call CONTOURC (v.) to estimate coordinates (2xN matrix [LONS;LATS]) for the
% isobath(s) at depth(s) DEPS, from the bathymetry field BATH (a STRUCT with
% fields .lon,.lat,.field). 
%
% If DEPS is a matrix OR a cell array, COORDS is a cell array of matrices.
%
% If DEPS is a scalar and MINDEP is a numeric scalar < DEPS, returns coords
% of all points in BATH at depths within the inclusive range [MINDEP,DEPS].
%
% Last Saved Time-stamp: <Fri 2016-04-15 15:15:33 Eastern Daylight Time gramer>

  if ( ~isfield(bath,'lon') || ~isfield(bath,'lat') || ~isfield(bath,'field') )
    error('First arg BATH must be a valid bathymetry STRUCT!');
  end;
  if ( ~exist('deps','var') || isempty(deps) )
    error('Second arg DEPS must be a matrix or cell array of scalars!');
  end;

  if ( exist('mindep','var') && ~isempty(mindep) )
    if ( ~isscalar(mindep) || ~isnumeric(mindep) || ~isscalar(deps) || ~isnumeric(deps) || mindep >= deps )
      error('If MINDEP is given, it and DEPS must be numeric scalars with MINDEP < DEPS');
    else
      [ix,jx] = find(mindep <= bath.field & bath.field <= deps);
      lon = bath.lon(jx(:));
      lat = bath.lat(ix(:));
      coords = [ lon(:)' ; lat(:)' ];
    end;

  else

    if ( ~iscell(deps) )
      deps = num2cell(deps);
    end;

    for ix = 1:numel(deps)
      dep = deps{ix};
      coord = contourc(bath.lon,bath.lat,bath.field,[dep,dep]);
      if ( isempty(coord) && dep == 0 )
        % Some bathymetry fields may not quite reach shore
        coord = contourc(bath.lon,bath.lat,bath.field,[-1,-1]);
      end;
      if ( ~isempty(coord) )
        dep = coord(1,1);
        % Replace depth and length markers with NaN "breaks"
        brkix = 1;
        while ( brkix < size(coord,2) )
          len = coord(2,brkix);
          coord(:,brkix) = nan;
          brkix = brkix + len + 1;
        end;
        % Remove initial "break"
        if ( isnan(coord(1,1)) && isnan(coord(2,1)) )
          coord(:,1) = [];
        end;
      end;
      if ( numel(deps) == 1 )
        coords = coord;
      else
        coords{ix} = coord;
      end;
    end;
  end;

return;
