function [sites,bath] = find_ngdc_slope_sites(coords_or_stns,bath,npts,method,extrapval,fldnm)
%function [sites,bath] = find_ngdc_slope_sites(coords_or_stns,bath,npts,method,extrapval,fldnm)
%
% Find depth, slope (beta) and slope angle (beta_ang), and isobath angle
% (iso_ang - all angles in deg True), for one or more sites, from bathymetric
% data set BATH. First arg may be a cell array of coord vectors {LONS,LATS},
% a STRUCT with fields .lons and .lats, or be a STRUCT SITES having fields,
% e.g., SITES.stn1.lon,SITES.stn1.lat, SITES.stn2.lon,SITES.stn2.lat, etc.
%
% BATH if given must be STRUCT with field FLDNM (DEFAULT: 'ngdc_hires_bathy').
% N.B.: Unlike the similarly named "FIND_NGDC_SLOPE", this function expects
% the STRUCT arg BATH to have *a field* that is a bathymetry STRUCT.
%
% This code assumes BATH.(FLDNM).lon, .lat in degrees, .field depths in [m].
%
% CALLS: FIND_NGDC_SLOPE (v.) to do slope calculations; and INTERP_FIELD (v.)
% to interpolate from BATH, the depth for each site in COORDS_OR_STNS. If a
% valid bathymetry STRUCT BATH is not given, calls READ_HIRES_BATHYMETRY (v.)
% to extract bathymetry in range GEOGRAPHIC_RADIUS_M(LONS,LATS) (v.).
% 
% EXAMPLE CALL:
%  >> % Smooth depth and 90-m finite difference slope estimates to 90 m squares:
% >>  % for higher-resolution bathymetry, naturally more smoothing is required.
%  >> if ( bath.ngdc_hires_bathy.yres < 15 )
%  >>  % 10 m bathymetry: MEDIAN 9x9-point (90 m) squares with >=10 valid points
%  >>  [sites,bath] = find_ngdc_slope_sites({lons,lats},bath,9,{@nanmedian,9,9,10});
%  >> elseif ( bath.ngdc_hires_bathy.yres < 50 )
%  >>  % 30 m bathymetry: MEDIAN 3x3-point (90 m) squares with >=4 valid points
%  >>  [sites,bath] = find_ngdc_slope_sites({lons,lats},bath,3,{@nanmedian,3,3,4});
%  >> else
%  >>  % 90 m bathymetry: Simple finite difference and linear interpolation
%  >>  [sites,bath] = find_ngdc_slope_sites({lons,lats},bath,2);
%  >> end;
% 
% Last Saved Time-stamp: <Sat 2019-04-20 16:49:42 Eastern Daylight Time gramer>

  if ( ~exist('coords_or_stns','var') )
    error('Function requires at least one argument');
  end;
  if ( ~exist('bath','var') )
    bath = [];
  end;
  if ( ~exist('npts','var') )
    npts = 3;
  end;
  if ( ~exist('method','var') )
    method = '*linear';
  end;
  if ( ~exist('extrapval','var') )
    extrapval = nan;
  end;
  if ( ~exist('fldnm','var') || isempty(fldnm) )
    fldnm = 'ngdc_hires_bathy';
  end;


  sites.lons = [];
  sites.lats = [];
  if ( iscell(coords_or_stns) )
    sites.lons = coords_or_stns{1};
    sites.lats = coords_or_stns{2};
    if ( numel(coords_or_stns) >= 3 )
      sites.stnms = coords_or_stns{3};
    end;
  elseif ( isstruct(coords_or_stns) )
    stns = coords_or_stns;
    if ( numel(stns) > 1 && isfield(stns,'lon') && isfield(stns,'lat') )
      % STNS is a matrix of STRUCT
      for stix=1:numel(stns)
        if ( ~isfield(sites,'stnms') )
          sites.stnms = {};
        end;
        if ( isfield(stns,'station_name') )
          sites.stnms{end+1} = stns(stix).station_name;
        end;
        sites.lons(end+1) = stns(stix).lon;
        sites.lats(end+1) = stns(stix).lat;
      end; %for stix=1:numel(stns)

    elseif ( isfield(stns,'lons') && isfield(stns,'lats') )
      % STNS is a STRUCT with fields that are matrices
      if ( isfield(stns,'stnms') )
        sites.stnms = stns.stnms;
      end;
      sites.lons = stns.lons;
      sites.lats = stns.lats;

    else
      % STNS is a STRUCT with fields that are themselves station STRUCTs
      for cf=fieldnames(stns)'
        if ( isfield(stns.(cf{:}),'lon') && isfield(stns.(cf{:}),'lat') )
          if ( ~isfield(sites,'stnms') )
            sites.stnms = {};
          end;
          sites.stnms{end+1} = cf{:};
          sites.lons(end+1) = stns.(cf{:}).lon;
          sites.lats(end+1) = stns.(cf{:}).lat;
        end;
      end; %for cf=fieldnames(stns)'
    end; %if ( isfield(stns,'lons') && isfield(stns,'lats') ) else
  end;
  coords_or_stns = []; clear coords_or_stns;


  if ( ~isvector(sites.lons) || ~isnumeric(sites.lons) || ~isvector(sites.lats) || ~isnumeric(sites.lats) )
    error('First arg must be a cell array of coord vectors {LONS,LATS}, a STRUCT with fields .lons and .lats, or a STRUCT with site fields (i.e., having SITES.stn1.lon,SITES.stn1.lat, etc.)');
  end;


  if ( ~exist('bath','var') )
    [rx,ry] = geographic_radius_m(sites.lons,sites.lats); 
    bath.rx = rx + 1e3;
    bath.ry = ry + 1e3;
    bath.lon = mean(sites.lons(:));
    bath.lat = mean(sites.lats(:));
    bath = read_hires_bathymetry(bath,[bath.rx,bath.ry]);
  elseif ( ~isfield(bath,fldnm) )
    error('BATH if specified must be STRUCT with field .%s',fldnm);
  end; %if ( ~exist('bath','var') ) else

  % sites.depths = interp2(bath.(fldnm).lon,bath.(fldnm).lat,...
  %                        bath.(fldnm).field,sites.lons,sites.lats,method,extrapval);
  sites.depths = interp_field(bath.(fldnm).lat,bath.(fldnm).lon,...
                              bath.(fldnm).field,sites.lats,sites.lons,method,extrapval);
  [sites.betas, sites.beta_angs, sites.iso_angs, bath.(fldnm)] = ...
      find_ngdc_slope(bath.(fldnm),sites.lons,sites.lats,npts,method,extrapval);

return;
