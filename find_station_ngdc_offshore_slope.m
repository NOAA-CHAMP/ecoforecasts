function [stn,transects] = find_station_ngdc_offshore_slope(stn_or_stnm,rad)
%function [stn,transects] = find_station_ngdc_offshore_slope(stn_or_stnm,rad)
%
% Estimate cross-shore slope (rise / run) of station STN or site named STNM,
% in the direction of greatest descent (normally, offshore), using the NGDC
% 3-arcsecond resolution Coastal Relief Model. Does a simple linear search
% for each one degree of azimuth around given STN.lon,STN.lat coordinates,
% using GET_STATION_BATHY_TRANSECT (v.), which may be VERY SLOW: so if file
% FULLFILE(GET_ECOFORECASTS_PATH('data'),[STNM,'_ngdc_transects.mat']) exists
% already, just load transects from that file. Optional RAD (DEFAULT 0.65km
% i.e., 7 points) is the RADIUS over which to estimate slope in each azimuth.
% STN is returned with the new (or overwritten) fields .ngdc_offshore_slope,
% .offshore_direction, and .isobath_orientation (=offshore_direction-90).
%
% Last Saved Time-stamp: <Wed 2016-07-13 14:14:00 Eastern Daylight Time gramer>
%
%warning('Method not well-tested! Consider calling find_ngdc_slope.m instead');

warning('Method not well-tested! Consider calling find_ngdc_slope.m instead');

  stn = get_station_from_station_name(stn_or_stnm);

  file_radius = 1;
  if ( ~exist('rad','var') || isempty(rad) )
    rad = 0.65;
  end;
  file_radius = ceil(rad);


  if ( ~isfield(stn,'ngdc_92m_bathy') )
    stn = get_ngdc_bathy_station(stn);
  end;

  matfname = fullfile(get_ecoforecasts_path('data'),[lower(stn.station_name),...
                      '_ngdc_transects.',num2str(file_radius),'km.mat']);

  if ( exist(matfname,'file') )
    %DEBUG:    disp(['Loading ',matfname]);
    load(matfname,'transects');
  else
    transects = get_station_bathy_transect(stn,file_radius,0:359);
    %DEBUG:    disp(['Saving ',matfname]);
    save(matfname,'transects');
  end;

  dx = 92; %[m]

  max_bet=-inf;
  offshore_dir=0;
  radx = floor(rad/(dx/1e3));

  for az=1:numel(transects);
    % midx = ceil(length(transects(az).field)/2);
    % bet = mean(abs(diff(transects(az).field(midx:midx+radx)))/dx);
    [xerr,stnx] = min( ((transects(az).lon-stn.lon).^2) + ((transects(az).lat-stn.lat).^2) );
    if ( xerr>dx )
      error('Transect returned by GET_STATION_BATHY_TRANSECT is not near station!');
    end;

    % % bet = mean(abs(diff(transects(az).field(stnx:stnx+radx)))/dx);
    % bet = abs(diff(transects(az).field(stnx:stnx+radx))) / dx;
    bet = diff(transects(az).field(stnx:stnx+radx)) / dx;
    % bet = mean([max(bet(:)),min(bet(:))]);
    bet = max(bet(:));
    if ( bet > max_bet )
      max_bet = bet;
      offshore_dir = az;
    end;
  end;

  if ( nargout < 2 )
    transects=[]; clear transects;
  end;

  stn.ngdc_offshore_slope = bet;
  stn.offshore_direction = offshore_dir;
  stn.isobath_orientation = offshore_dir - 90;

return;
