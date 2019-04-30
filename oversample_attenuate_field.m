function [newfld,fldlat,fldlon] = oversample_attenuate_field(lostr,varargin)
%function [newfld,fldlat,fldlon] = oversample_attenuate_field(lostr,[N,][histr|{hilon,hilat,hidts}|[]],[dep],[method],[noNaNs])
%
% Use interpolation to arbitrarily increase resolution of a LOwer-resolution
% time series field STRuct LOSTR, to create a new time-series field: mask all
% land gridpoints (pts where DEP >= 0), and attenuate resulting field so that
% it reaches zero at the coast. Caller may specify a scalar factor N by which
% to increase resolution of original field: this option is useful if memory
% constraints prevent direct interpolation to the new HIgher resolution. 
%
% LOSTR may be a STRUCT with fields .lon,.lat,.date,.field, or a cell array
% with matrices in that order, i.e., sized {1xN, 1xM, 1xT, and TxMxN}, resp.
%
% Optional HISTR may be a STRUCT with fields .lon,.lat,.date, or a cell array
% with matrices in that order, i.e., sized {1xP, 1xQ, 1xR}, resp. If HISTR is
% 2-D, specify it as a 2-cell array: LOSTR.date is used. If HISTR is missing 
% or empty, do not oversample: just attenuate at native LOSTR resolution.
%
% Optional DEP may be a STRUCT with fields .lon,.lat,.field, or a cell array 
% with matrices in that order, i.e., sized {1xP, 1xQ, QxP}, resp., or DEP may
% just be a QxP field: in that case, its resolution must match that of HISTR.
%
% Optional arg METHOD is passed to INTERP2, INTERP3 (DEFAULT: '*linear'). If
% METHOD contains the string 'depth', does linear depth attenuation for all
% depths between -20 m and 0.
%
% If optional arg NONANS (DEFAULT: False), replace NaN with 0 in LOSTR.
%
% SEE ALSO oversample_field.m
%
% Last Saved Time-stamp: <Mon 2017-07-03 16:13:59 Eastern Daylight Time gramer>

  if ( iscell(lostr) )
    locel = lostr;
    lostr = [];
    lostr.lon = locel{1};
    lostr.lat = locel{2};
    lostr.date = locel{3};
    lostr.field = locel{4};
    clear locel
  elseif ( ~isfield(lostr,'field') )
    error('LOSTR must be STRUCT with fields .lon,.lat,.date,.field, or a cell array of matrices');
  end;

  args = varargin;
  if ( numel(args) < 1 )
    error('Too few arguments: must specify N or HISTR or both');
  elseif ( numel(args) > 5 )
    error('Too many arguments');
  end;

  N=[];
  if ( isscalar(args{1}) && isnumeric(args{1}) )
    N = args{1};
    % SIDE-EFFECT: If N is specified, take it off the ARGS list
    args(1) = [];
  end;

  histr = [];
  if ( numel(args) >= 1 )
    histr = args{1};
  end;
  if ( isempty(histr) )
    histr.lon = lostr.lon;
    histr.lat = lostr.lat;
    histr.date = lostr.date;
    histr.field = lostr.field;
  elseif ( iscell(histr) )
    hicel = histr;
    histr = [];
    histr.lon = hicel{1};
    histr.lat = hicel{2};
    if ( numel(hicel) < 3 )
      histr.date = lostr.date;
    else
      histr.date = hicel{3};
    end;
    clear hicel
  elseif ( ~isfield(histr,'lon') || ~isfield(histr,'lat') || ~isfield(histr,'date') )
    error('HISTR must be STRUCT with fields .lon,.lat,.date, or a cell array of matrices');
  end;

  dep = [];
  if ( numel(args) >= 2 )
    dep = args{2};
    if ( iscell(dep) )
      depcel = dep;
      dep=[]; clear dep
      deplon = depcel{1};
      deplat = depcel{2};
      dep = depcel{3};
      depcel=[]; clear depcel
    elseif ( isstruct(dep) )
      deplon = dep.lon;
      deplat = dep.lat;
      dep = dep.field;
    elseif ( ismatrix(dep) )
      deplon = histr.lon;
      deplat = histr.lat;
    else
      error('Depth arg DEP is neither CELL, STRUCT, nor MATRIX');
    end;
  end;

  method = [];
  if ( numel(args) >= 3 && ~isempty(args{3}) )
    method = args{3};
  end;
  depth_attenuate = false;
  if ( ~isempty(strfind(method,'depth')) )
    depth_attenuate = true;
    method = regexprep(method,'[,]*depth[a-z_]*[,]*','');
  end;
  if ( isempty(method) )
    % DEFAULT for INTERP[23N]
    method = '*linear';
  end;

  noNaNs = false;
  if ( numel(args) >= 4 && ~isempty(args{4}) )
    noNaNs = args{4};
  end;

  % Replace model landmask (NaNs in LOSTR) with fully attenuated gridpoints
  if ( noNaNs )
    lostr.field(isnan(lostr.field)) = 0;
  end;

  % Use bilinear interpolation to increase field resolution, ...
  if ( isempty(N) )
    fldlat = histr.lat;
    fldlon = histr.lon;
  else
    fldlat = linspace(lostr.lat(1),lostr.lat(end),numel(lostr.lat)*N);
    fldlon = linspace(lostr.lon(1),lostr.lon(end),(numel(lostr.lon)*N)+1);
  end;
  [flat,fdts,flon] = meshgrid(fldlat,histr.date,fldlon);
  warning('OFF','MATLAB:chckxy:IgnoreNaN');
  fldfld = interp3(lostr.lat,lostr.date,lostr.lon,lostr.field,flat,fdts,flon,method);
  warning('ON','MATLAB:chckxy:IgnoreNaN');

  if ( isempty(dep) )
    % If caller just wants an un-attenuated, oversampled field...
    newfld = fldfld;

  else

    % If caller wants an attenuated field...

    % ... First over- or undersample the depth field to match our resolution, ...
    flddep = interp2(deplon,deplat,dep,fldlon',fldlat,method);
    % ... Then mark land gridpoints in result as zero, ...
    fldfld(:,find(flddep >= 0)) = 0;

    % ... AND EITHER...
    if ( depth_attenuate )
      % ... Linearly attenuate (scale) field for depths in range 0 to -20 m ...
      depscl = flddep ./ -20;
      depscl(depscl < 0) = 0;
      depscl(depscl > 1) = 1;
      depscl = permute(repmat(depscl,[1,1,size(fldfld,1)]),[3,1,2]);
      newfld = fldfld .* depscl;
      depscl=[]; clear depscl;

    else
      % ... OR, Interpolate using METHOD, to attenuate field near the coastline.
      newfldsz = [length(histr.date),length(fldlat),length(fldlon)];
      x = memory; avail_doubles = x.MaxPossibleArrayBytes/8; clear x;
      warning('OFF','MATLAB:chckxy:IgnoreNaN');
      if ( avail_doubles < (prod(newfldsz) / 4) )
        disp('Insufficient memory for true 3D interpolation: this may take a while...');
        %newfld = repmat(nan,size(histr.field));
        newfld = repmat(nan,newfldsz);
        for lonix=1:numel(fldlon)
          newfld(:,:,lonix) = interp3(flat,fdts,flon,fldfld,...
                                      fldlat,histr.date,fldlon(lonix),method);
        end;
      else
        % We have enough memory to save the caller some time
        newfld = interp3(flat,fdts,flon,fldfld,fldlat,histr.date,fldlon,method);
      end;
      warning('ON','MATLAB:chckxy:IgnoreNaN');
    end;
  end;

  % ... And finally, land-mask gridpoints in result back to NaNs
  newfld(:,find(flddep >= 0)) = nan;

return;
