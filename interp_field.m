function vals = interp_field(lats,lons,fld,sitelat,sitelon,interpMethod,extrapVal,triangular,doTiming)
%function vals = interp_field(lats,lons,fld,sitelat,sitelon,[interpMethod[,extrapVal[,triangular]]])
%
% Return interpolation on a *time series field*: LATS and LONS are monotonic
% vectors of length N and M, resp.; FLD is a matrix of size DxNxM; SITELAT,
% SITELON are arrays with P elements each, e.g., P==length(SITELAT(:)), with
% coordinates of individual sites. Optional INTERPMETHOD (DEFAULT 'linear'
% for bilinear) and EXTRAPVAL are passed to INTERPN (v.), but if EXTRAPVAL is
% the string 'warning', NaN is passed as an extrapVal and the user receives a
% warning listing ALL the coordinates which lie outside of the domain. VALS
% is a DxP array (a time series of vectors) of values interpolated on FLD.
%
% If optional TRIANGULAR is true, then INTERPMETHOD is passed to GRIDDATA
% (v.) instead, to perform triangle-based 'linear' or 'cubic' interpolation.
% If TRIANGULAR, EXTRAPVAL may still be 'warning'; otherwise, it is ignored.
% NOTE: TRIANGULAR causes a loop to run over all dates: this may be slow!
%
% INTERPMETHOD may also be a function handle, e.g., @NANMEAN (function name
% strings not supported), to apply an accumulator to the 4 field grid-points
% nearest each site. If TRIANGULAR, 3 nearest points are used. Calls function
% on both spatial dimensions, e.g., as if, MIN(MIN(PERMUTE(FLD,[2 3 1]))).
% If INTERPMETHOD is a string not recognized as an option to INTERPN, INTERP2
% or GRIDDATA (v.), tries turning into a function_handle using STR2FUNC (v.).
%
% If INTERPMETHOD is a character string, then the TRIANGULAR option may also
% be specified by prepending 'triangular,' to it, e.g., 'triangular,linear'.
% In this case, the value of the TRIANGULAR arg (if passed) will be ignored.
%
% Finally, INTERPMETHOD may be a 2-, 3-, or 4-elt cell array: INTERPMETHOD{1}
% may be any of the options defined above; if it is evaluated as a function
% handle, INTERPMETHOD{2} is the X-radius of the box around SITELAT,SITELON
% passed to the function; INTERPMETHOD{3} if present is a separate Y-radius;
% INTERPMETHOD{4} if present is a minimum number of FINITE values per box -
% any YxX box with fewer finite (non-NaN) values results in a NaN gridpoint.
%
% NOTE: INTERP_FIELD also handles special cases where FLD is 1xNxM or NxM.
%
% Last Saved Time-stamp: <Thu 2019-01-17 15:42:24 Eastern Standard Time gramer>

  sitelat = sitelat(:);
  sitelon = sitelon(:);

  % Singleton special case: convert NxM array FLD to 1xNxM 3-d array
  if ( ndims(fld) == 2 && size(fld,1) == numel(lats) && size(fld,2) == numel(lons) )
    fld = reshape(fld,[1,size(fld,1),size(fld,2)]);
  end;

  if ( ndims(fld) ~= 3 || size(fld,2) ~= numel(lats) || size(fld,3) ~= numel(lons) )
    error('LATS,LONS,FLD must have N elements, M elements, and be DxNxM, resp.!');
  end;
  if ( ~isnumeric(sitelat) || ~isnumeric(sitelon) || numel(sitelat) ~= numel(sitelon) )
    error('SITELAT,SITELON must be numerical arrays with the same number of elements!');
  end;
  if ( ~exist('interpMethod','var') || isempty(interpMethod) )
    interpMethod = 'linear';
  end;
  interpXRadius = 1;
  interpYRadius = 1;
  interpMinFin = [];
  if ( iscell(interpMethod) && 2<=numel(interpMethod) && numel(interpMethod)<=4 ...
       && isnumeric(interpMethod{2}) && isnumeric(interpMethod{end}) )
    interpXRadius = interpMethod{2};
    if ( numel(interpMethod)<3 )
      interpYRadius = interpXRadius;
    else
      interpYRadius = interpMethod{3};
    end;
    if ( numel(interpMethod)>=4 )
      interpMinFin = interpMethod{4};
    end;
    interpMethod = interpMethod{1};
  end;
  if ( ~ischar(interpMethod) && ~isa(interpMethod,'function_handle') )
    error('INTERPMETHOD must be a string (v. INTERPN) or a function handle!');
  end;

  if ( ~exist('extrapVal','var') || isempty(extrapVal) )
    extrapVal = [];
  end;
  if ( ~exist('triangular','var') || isempty(triangular) )
    triangular = false;
  end;

  if ( ~exist('doTiming','var') || isempty(doTiming) )
    doTiming = false;
  end;

  if ( ischar(interpMethod) )
    if ( interpMethod(1) == '*' )
      interpMethod = interpMethod(2:end);
    end;

    % Let user prepend 'triangular' (or any variant) to specify TRIANGULAR
    [interpMethodT,interpMethodR]=strtok(interpMethod,',; ');
    if ( ~isempty(interpMethodT) && ~isempty(interpMethodR) )
      if ( lower(interpMethodT(1)) == 't' )
        triangular = true;
      end;
      interpMethod = strtok(interpMethodR,',; ');
    end;

    % Let user specify 3- or 4-gridpoint accumulator function by name also
    switch ( lower(interpMethod) )
     case {'nearest','linear','spline','cubic','v4'},
      % Do nothing - pass through to INTERPN/INTERP2/GRIDDATA
     otherwise,
      % If it would be a valid function handle, try to use it
      if ( ismember(exist(interpMethod),[2:6]) )
        interpMethod = str2func(interpMethod);
      else
        error('INTERPMETHOD bad function name "%s"',interpMethod);
      end;
    end;
  end;


  vals = [];

  d = 1:size(fld,1);

  if ( isa(interpMethod,'function_handle') )
    % Use caller's function handle on both spatial dimensions, at each site
    vals = repmat(nan,[size(fld,1),numel(sitelon)]);
    if doTiming; btic = tic; tic; end;
    for ix=1:numel(sitelon)
      if ( doTiming && mod(ix,1e5) == 0 )
        toc,
        disp(ix);
        tic,
      end;
      if ( min(lons) > sitelon(ix) || sitelon(ix) > max(lons) || ...
           min(lats) > sitelat(ix) || sitelat(ix) > max(lats) )
        if ( isempty(extrapVal) || ~isnumeric(extrapVal) )
          vals(:,ix) = NaN;
        else
          vals(:,ix) = extrapVal;
        end;
      else
        [ig,lonix] = min(abs(lons-sitelon(ix)));
        [ig,latix] = min(abs(lats-sitelat(ix)));
        lonerr = lons(lonix) - sitelon(ix);
        laterr = lats(latix) - sitelat(ix);
        if ( triangular )
          % Triangulation (below) ignores coordinate projection, so we will also
          % NOTE this by default calls @INTERPMETHOD with only *TWO* points!
          if (lonerr<laterr)
            lonix = unique([lonix lonix-sign(lonerr)]);
          else
            latix = unique([latix latix-sign(laterr)]);
          end;
        else
          % lonix = unique([lonix lonix-sign(lonerr)]);
          % latix = unique([latix latix-sign(laterr)]);
          interpXNearRadius = interpXRadius - 1;
          interpYNearRadius = interpYRadius - 1;
          % LJG 2015-05-14: if lonerr,laterr==0, only one point would be used! 
          sgnlon=sign(lonerr); if ~sgnlon; sgnlon=1; end;
          sgnlat=sign(laterr); if ~sgnlat; sgnlat=1; end;
          lonixen = unique([lonix+(sgnlon*interpXNearRadius) ...
                            lonix-(sgnlon*interpXRadius)]);
          lonix = min(lonixen):max(lonixen);
          lonix(1>lonix | lonix>size(fld,3)) = [];
          latixen = unique([latix+(sgnlat*interpYNearRadius) ...
                            latix-(sgnlat*interpYRadius)]);
          latix = min(latixen):max(latixen);
          latix(1>latix | latix>size(fld,2)) = [];
        end;

        %dat = fld(:,lonix,latix); % D@mn dyslexia!
        dat = fld(:,latix,lonix);

        % vals(:,ix) = interpMethod(dat(:,:)')';
        if ( ~isempty(interpMinFin) )
          goodix = find(sum(isfinite(dat(:,:)')',2)>=interpMinFin);
        else
          goodix = 1:size(vals,1);
        end;
        vals(goodix,ix) = interpMethod(dat(goodix,:)')';
      end;

    end; %for ix=1:numel(sitelon)
    if doTiming; toc, toc(btic), end;

  elseif ( triangular )
    % Use GRIDDATA for triangle-based linear (3-point) or cubic interpolation
    vals = repmat(nan,[size(fld,1),numel(sitelon)]);
    % WARNING: This loop can be VERY SLOW!
    %DEBUG:    warning('This may take a while...');
    if doTiming; btic = tic; end;
    for dtix=d(:)'
      vals(dtix,:) = griddata(lons,lats,squeeze(fld(dtix,:,:)),sitelon,sitelat,interpMethod);
    end;
    if doTiming; toc(btic), end;

  elseif ( numel(d) == 1 )
    if doTiming; btic = tic; end;
    % Use INTERP2 to handle singleton special case that INTERPN does not 
    % NOTE: Prime yields 1xN return value instead of Nx1...
    fld = squeeze(fld);
    if ( isempty(extrapVal) )
      vals = interp2(lons,lats,fld,sitelon,sitelat,interpMethod)';
    elseif ( strncmpi(extrapVal,'w',1) )
      vals = interp2(lons,lats,fld,sitelon,sitelat,interpMethod,NaN)';
    else
      vals = interp2(lons,lats,fld,sitelon,sitelat,interpMethod,extrapVal)';
    end; %if ( isempty(extrapVal) ) elseif else
    if doTiming; toc(btic), end;

  else
    if doTiming; btic = tic; end;
    % Use INTERPN for nearest, bilinear (4-point), or cubic interpolation
    [D,SITELAT] = ndgrid(d,sitelat);
    [D,SITELON] = ndgrid(d,sitelon);

    if ( isempty(extrapVal) )
      vals = interpn(d,lats,lons,fld,D,SITELAT,SITELON,interpMethod);
    elseif ( strncmpi(extrapVal,'w',1) )
      vals = interpn(d,lats,lons,fld,D,SITELAT,SITELON,interpMethod,NaN);
    else
      vals = interpn(d,lats,lons,fld,D,SITELAT,SITELON,interpMethod,extrapVal);
    end; %if ( isempty(extrapVal) ) elseif else
    if doTiming; toc(btic), end;

  end; %if (function_handle) elseif (triangular) elseif (numel(d) == 1) else


  if ( strncmpi(extrapVal,'w',1) )
    % If caller wants warnings for out-of-bounds sites
    dlat = max(diff(unique(lats)));
    dlon = max(diff(unique(lons)));
    badix = find( min(lons)>(sitelon+dlon) | (sitelon-dlon)>max(lons) | ...
                  min(lats)>(sitelat+dlat) | (sitelat-dlat)>max(lats) );
    if ( ~isempty(badix) )
      warning( 'INTERP_FIELD:OutOfBounds',...
               'The following coordinates are outside the domain!\n%s',...
               sprintf('%g,%g\n',[sitelat(badix),sitelon(badix)]') );
    end;
  end;

return;
