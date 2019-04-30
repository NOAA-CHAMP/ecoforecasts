function [LONS,LATS,varargout] = get_coords_from_args(varargin)
%function [LONS,LATS,varargout] = get_coords_from_args([STN|STNS(:)|FIELD|COORDS(1:2,:)],extra_fldnm)
%
% Extract vectors of LONS and LATS from first arg, which should be a struct
% with fields .lon,.lat, or .lons,.lats, a 2xN numeric matrix, or CxN cell.
%
% Last Saved Time-stamp: <Fri 2019-02-15 17:11:29 Eastern Standard Time gramer>

  args = varargin;
  if ( exist('varargout','var') && ~isempty(varargout) )
    varargout{1} = '';
  end;
  
  if ( isnumeric(args{1}) )
    LONS = args{1}(1,:);
    LATS = args{1}(2,:);
    
  elseif ( iscell(args{1}) )
    LONS = args{1}{1};
    LATS = args{1}{2};

  elseif ( isstruct(args{1}) )
    str = args{1};
    
    if ( isfield(str,'lon') )
      LONS = str.lon;
      LATS = str.lat;
    elseif ( isfield(str,'lons') )
      LONS = str.lons;
      LATS = str.lats;
    else
    end;
    
    if ( numel(args) > 1 )
      if ( isfield(str,args{2}) )
        varargout{1} = str.(args{2});
      end;
    end;
    
  end; %if ( isnumeric(args{1}) ) elseif ( isstruct(args{1}) ) else

return;
