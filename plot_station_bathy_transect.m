function results = plot_station_bathy_transect(stn,rad,azs,bathynm)
%function results = plot_station_bathy_transect(stn,rad,azs,bathynm)
%
% Plot (and return in struct vector RESULTS) the bathymetric transect(s) of
% radius RAD and direction(s) AZS (degT), centered at the location of station
% specified in struct STN, or of station named STN. If BATHYNM is given and is
% anything other than the DEFAULT ('ngdc_92m_bathy'), extract transects from
% field STN.(BATHYNM) (or struct BATHYNM, if it has fields .lon,.lat,.field)
% instead.  If the arg BATHYNM is not given or is empty, and if the field 
% STN.ngdc_92m_bathy is not already present, calls GET_NGDC_BATHY_STATION (v.)
% to extract bathymetry. THIS MAY BE SLOW! Calls STATION_FIELD_TRANSECT (v.) on
% station and field BATHYNM, to generate transect(s). RAD DEFAULT value is 8
% (kilometers). If AZS is non-scalar, then size(RESULTS) == [1 length(AZS)].
%
% Last Saved Time-stamp: <Fri 2010-11-19 11:08:17 Eastern Standard Time gramer>

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
    bathynm = 'ngdc_92m_bathy';
    if ( ~isfield(stn,bathynm) )
      % Note: If STN was a CHAR, it will be a STRUCT after this call
      stn = get_ngdc_bathy_station(stn);
    end;
    rng = [-rad:(0.092):rad];
  elseif ( isstruct(bathynm) )
    stn.temporary_bathy = bathynm;
    bathynm = 'temporary_bathy';
    rng = [-rad:(0.092):rad];
  else
    rng = linspace(-rad,rad);
  end;

  figure;
  maxigraph;
  hold on;
  for ix = 1:length(azs)
    az = azs(ix);
    results(ix) = station_field_transect(stn,bathynm,rng,az);
    fld(ix,1:length(results(ix).field)) = results(ix).field(:);
  end;
  plot(rng,fld);
  xlabel('Distance from site [km]');
  ylabel('Depth [m]');
  legend(num2str(azs','Az=%g^o'), 'Location','Best');

  if ( nargout < 1 )
    results = []; clear results;
  end;

return;
