function [results,stn,lh] = plot_bathy_transect(stn,rad,azs,bathynm)
%function [results,stn,lh] = plot_bathy_transect(stn,rad,azs,bathynm)
%
% Plot (and return in struct vector RESULTS) the bathymetric transect(s) of
% radius RAD and direction(s) AZS (degT), centered at the location of station
% specified in struct STN or of station named STN. If BATHYNM is given and is
% anything other than the DEFAULT 'ngdc_hires_bathy', extract transects from
% field STN.(BATHYNM) (or struct BATHYNM, if it has fields .lon,.lat,.field)
% instead.  If the arg BATHYNM is not given or is empty, and if the field 
% STN.ngdc_hires_bathy is not already present, calls READ_NGDC_BATHYMETRY (v.)
% to extract bathymetry (MAY BE SLOW). Calls STATION_FIELD_TRANSECT (v.) on
% station and field BATHYNM, to generate transect(s). RAD DEFAULT value is 8
% (kilometers). If AZS is non-scalar, then size(RESULTS) == [1 length(AZS)].
%
% AS OF Feb. 2018, if STN is a CELL, 1st arg is interpreted as above. Second
% cell elt. specifies site to use as "origin" (as station STRUCT or vector).
%
% Last Saved Time-stamp: <Sat 2018-03-31 17:26:14 Eastern Daylight Time gramer>

  org=[];
  if ( iscell(stn) )
    org = stn{2};
    if ( all(isfield(org,{'lon','lat'})) )
      org = [org.lon,org.lat];
    end;
    stn = stn{1};
  end;
  if ( ischar(stn) )
    stnm = stn;
    stn = []; clear stn;
    stn.station_name = stnm;
  end;
  if ( any(~isfield(stn,{'lon','lat'})) )
    if ( ~isfield(stn,'station_name') )
      error('First arg STN must have fields .lon and.lat, or be a name string');
    end;
    [stn.lon,stn.lat,stn.depth] = get_station_coords(stn.station_name);
  end;
  if ( ~exist('rad','var') || isempty(rad) )
    rad = 8;
  end;
  if ( ~exist('azs','var') || isempty(azs) )
    azs = [0];
  end;
  if ( ~exist('bathynm','var') || isempty(bathynm) )
    bathynm = 'ngdc_hires_bathy';
    if ( ~isfield(stn,bathynm) )
      % Note: If STN was a CHAR, it will be a STRUCT after this call
      stn = read_hires_bathymetry(stn,(rad+1)*1e3);
    end;
    %rng = [-rad:(0.092):rad];
  elseif ( isstruct(bathynm) )
    stn.temporary_bathy = bathynm;
    bathynm = 'temporary_bathy';
    %rng = [-rad:(0.092):rad];
  else
    %rng = linspace(-rad,rad);
  end;
  if ( ~isfield(stn,bathynm) )
    error('No bathymetry field %s',bathynm);
  elseif ( ~isfield(stn.(bathynm),'lon') || ~isfield(stn.(bathynm),'lat') )
    error('Unknown bathymetry format in STN.%s',bathynm);
  else
    latdif = min(diff(unique(stn.(bathynm).lat(:))));
    londif = min(diff(unique(stn.(bathynm).lon(:))))/cosd(stn.(bathynm).lat(1));
    mindist = min([...
        distance_wgs84(stn.(bathynm).lat(1),stn.(bathynm).lon(1),stn.(bathynm).lat(1),stn.(bathynm).lon(1)+londif),...
        distance_wgs84(stn.(bathynm).lat(1),stn.(bathynm).lon(1),stn.(bathynm).lat(1)+latdif,stn.(bathynm).lon(1))...
                  ]);
    rng = [-rad:(mindist):+rad];
  end;

  % % Caller should call FIGURE (or FMG, etc.) first!
  % figure;
  % maxigraph;
  % hold on;
  for ix = 1:length(azs)
    az = azs(ix);
    results(ix) = station_field_transect(stn,bathynm,rng,az,org);
    fld(ix,1:length(results(ix).field)) = results(ix).field(:);

    results(ix).range = rng;
    results(ix).depths = fld(ix,1:length(results(ix).field));
  end;
  lh = plot(rng,fld);
  xlabel('Distance from site [km]');
  ylabel('Depth [m]');
  legend(num2str(azs','Az=%g^o'), 'Location','Best');

  if ( nargout < 2 )
    stn = []; clear stn;
    if ( nargout < 1 )
      results = []; clear results;
    end;
  end;

return;
