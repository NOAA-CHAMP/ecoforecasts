function [c,h] = contour_field(fldstr,varargin)
%function [c,h] = contour_field(fldstr[,contourFun][,fld,ix_or_fun,VARARGIN])
%
% Use CONTOURF to plot field FLD with coordinates FLDSTR.lon and FLDSTR.lat.
% DEFAULT FLD is field FLDSTR.field. If IX_OR_FUN is a non-empty index list,
% plot SQUEEZE(FLD(IX,:,:)), if FUNCTION_HANDLE (DEFAULT: FUN=@NANMEAN), plot
% SQUEEZE(FUN{1}(FLD,FUN{:})). If FLD is 2-D, i.e., not a time series field,
% IX_OR_FUN is ignored. Sets DASPECT ratio based on latitude. If CONTOURFUN
% is a FUNCTION_HANDLE (cannot be a string), plot using that function instead
% of CONTOURF. All other args are passed on to CONTOURFUN as its args 4-N.
%
% NOTE: To specify IX_OR_FUN, arg FLD must be given first: FLD may be [].
% To pass additional args to CONTOURFUN, both FLD and IX_OR_FUN may be [].
%
% SAMPLE CALLS:
%   Plot the NANMEAN of time series (3-D) field SST.field vs. SST.lon, SST.lat:
%   >> fmg; contour_field(sst);
%   Plot the time minimum of time series series HYCOM.lon,HYCOM.lat,HYCOM.speed:
%   >> fmg; contour_field(hycom,'speed',@NANMIN);
%   Plot bathymetry (2D) map BATH, using default CONTOURF with default contours
%   >> fmg; contour_field(bath);
%   Plot bathymetry map BATH, with filled contours every 5 m from 0 to -80.
%   >> fmg; contour_field(bath,[],[],-[0:5:80]);
%   Plot the 93rd percentile of time series field WW3.lon,WW3.lat,WW3.hs.field, with 10 contour lines
%   >> fmg; contour_field(ww3,@contour,ww3.hs.field,{@prctile,93},10);
%   Plot bathymetry map BATH using function SURF
%   >> fmg; contour_field(bath,@surf);
%
% Last Saved Time-stamp: <Sun 2018-08-19 17:52:20 Eastern Daylight Time gramer>

  args = varargin;
  nargs = numel(args);

  if ( ~isfield(fldstr,'lon') || length(fldstr.lon) < 2 || ...
       ~isfield(fldstr,'lat') || length(fldstr.lat) < 2 )
    error('First arg must be a STRUCT with vector fields .lon and .lat');
  end;

  contourFun = @contourf;
  if ( nargs >= 1 && isa(args{1},'function_handle') )
    contourFun = args{1};
    args(1) = [];
    nargs = nargs - 1;
  end;

  fld = [];
  if ( nargs >= 1 )
    if ( ~isempty(args{1}) )
      fld = args{1};
    end;
    args(1) = [];
    nargs = nargs - 1;
  end;
  if ( isempty(fld) )
    if ( ~isfield(fldstr,'field') )
      error('If FLD is empty (DEFAULT), FLDSTR must have a field ".field"');
    end;
    fld = fldstr.field;
  elseif ( ischar(fld) )
    fld = fldstr.(fld);
  elseif ( ~isnumeric(fld) )
    error('FLD should either be empty (DEFAULT), a CHAR fieldname, or a numeric matrix');
  end;

  % At this point, FLD must be a numeric matrix
  if ( ~isnumeric(fld) || isvector(fld) )
    error('FLD must resolve to a numeric matrix or N-d array');
  end;

  % Caller convenience: Remove singleton dimension(s)
  fld = squeeze(fld);

  % NOTE: If FLD is 2D at this point, IX will be ignored below...
  ix = @nanmean;
  ixargs = {};
  % If caller specified IX_OR_FUN...
  if ( nargs >= 1 )
    if ( ~isempty(args{1}) )
      ix = args{1};
    end;
    args(1) = [];
    nargs = nargs - 1;
    if ( iscell(ix) )
      ixargs = ix(2:end);
      ix = ix{1};
    end;
  end;

  if ( ndims(fld) > 2 && ~isempty(ix) )
    if ( isa(ix,'function_handle') )
      fld = ix(fld,ixargs{:});
    elseif ( isnumeric(ix) || islogical(ix) )
      fld = fld(ix,:,:);
    else
      error('IX was neither a function_handle nor an index list??');
    end;
  end;
  fld = squeeze(fld);
  if ( ~ismatrix(fld) || ndims(fld) ~= 2 )
    error('FLD (under IX_OR_FUN, if specified) should resolve to a 2D matrix');
  end;

  [c,h] = feval(contourFun,fldstr.lon,fldstr.lat,fld,args{:});
  daspect([1,cosd(fldstr.lat(1)),1]);
  % If caller specified contour levels
  if ( nargs > 0 && isnumeric(args{1}) && numel(args{1}) > 2 )
    set(gca,'CLim',[min(args{1}),max(args{1})]);
  end;
  colorbar;

  if ( exist('set_surf_cursor','file') )
    set_surf_cursor;
  end;

return;
