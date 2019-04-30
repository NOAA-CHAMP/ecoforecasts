function stn = station_advect_field(stn,UdotdelQ,euleru,eulerv,field)
%function stn = station_advect_field(stn,UdotdelQ,euleru,eulerv,field)
%
% Advect the time series of 2-D fields STN.(FIELD).field across grid position
% of coordinates STN.lon,STN.lat, using time series of quasi-Eulerian current
% vector components STN.(EULERU) and STN.(EULERV). Save resultant time series
% as STN.(UDOTDELQ). All STN fields named in args must have .date (datenum)
% subfields. *Both* STN and STN.(FIELD) should have .lon and .lat subfields.
%
% Last Saved Time-stamp: <Tue 2010-05-25 21:11:50 Eastern Daylight Time gramer>

  lon = stn.(field).lon;
  lat = stn.(field).lat;
  dts = stn.(field).date;
  fld = stn.(field).field;
  nx = size(fld,3);
  ny = size(fld,2);

  ix1 = find(ismember(floor(dts),floor(stn.(euleru).date)));
  ix2 = find(ismember(floor(stn.(euleru).date),floor(dts)));

  dts = dts(ix1);
  fld = fld(ix1,:,:);
  % Goofy INTERP3 (below) assumes 3-d matrix is in Y,X,TIME order!
  fld = permute(fld,[2 3 1]);


  % For simplicity, assume for now that FLD is on a 1x1 km grid and that
  % STN.(EULERU).data and STN.(EULERV).data are [m/s]. Convert to [km/h].
  % Later, use DISTANCE (v.) with STN.(FIELD).lat and STN.(FIELD).lon.

  eudts = stn.(euleru).date(ix2);
  euu = stn.(euleru).data(ix2) .* (3600/1e3);
  euv = stn.(eulerv).data(ix2) .* (3600/1e3);


  centerx = repmat(interp1(lon,1:nx,stn.lon),size(euu));
  centery = repmat(interp1(lat,1:ny,stn.lat),size(euv));
  % U and V are in target direction: subtract for source point!
  sourcex = centerx - euu;
  %sourcey = centery - euv;
  % BUT rows are normally reversed in SST fields, so *add* V! Grumble...
  sourcey = centery + euv;

  stn.(UdotdelQ).date = eudts;
  % Note: Advection beyond our domain will result in NaNs
  % BUT... still need special handling here for land mask!
  stn.(UdotdelQ).data = ...
      interp3(1:nx,1:ny,dts,fld,sourcex,sourcey,eudts,'linear',nan) - ...
      interp3(1:nx,1:ny,dts,fld,centerx,centery,eudts,'linear',nan);

return;
