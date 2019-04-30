function [stn,dts] = get_fdep_stevens_metadata(stcd_or_stn);
%function [stn,dts] = get_fdep_stevens_metadata(stcd_or_stn);
%
% Return station STRUCT STN with fields .station_name,.station_id,.long_name,
% and date range DTS for Florida Dept. of Environmental Protection Stevens
% Water Monitoring Systems station with one-letter station code STCD (or may
% pass in a STRUCT with .station_name 'fdepX', where X is the station code).
%
% Last Saved Time-stamp: <Wed 2016-08-10 17:31:59 Eastern Daylight Time gramer>

  persistent fdep_stns

  % SAMPLE URL (for station K):
  %  http://www.fldep-stevens.com/export-8722375.php?t=c&txtFDate=02/01/2013&selFTime=00&txtTDate=02/28/2013&selTTime=23&selChannel=&rdoDateTime=1&rdoWTemp=1&rdoWSpeed=2&rdoATemp=1&rdoPressure=1&rdoRainfall=1&rdoWLevel=1

  fdep_stns = { ...
      'A', '8720757', datenum(2011,05,01), now,                [29 36 55.5],[-81 12 17.7],   'Bings Landing' ; ...
      'B', '8722213', datenum(2011,05,01), now,                [27 28 3.17],[-80 18 3.06],   'Binney Dock' ; ...
      'C', '8728744', datenum(2011,05,01), now,                [29 40 27.51],[-85 03 29.05], 'Dry Bar' ; ...
      'D', '8728603', datenum(2011,05,01), now,                [29 47 08.68],[-84 52 31.21], 'East Bay' ; ...
      'E', '8725081', datenum(2011,05,01), now,                [26 05 35.3],[-81 47 53.3],   'Gordon River Inlet' ; ...
      'F', '8728703', datenum(2011,05,01), now,                [29 45 20.48],[-85 00 12.69], 'Little St. Marks River' ; ...
      'G', '8721843', datenum(2011,05,01), now,                [28 05 1.23],[-80 35 31.10],  'Melbourne' ; ...
      'H', '8725114', datenum(2011,05,01), now,                [26 08 30.4],[-81 47 24.2],   'Naples Bay' ; ...
      'I', '8728732', datenum(2011,05,01), now,                [29 36 05.00],[-85 01 39.73], 'Pilots Cove' ; ...
      'J', '8721147', datenum(2011,05,01), now,                [29 03 48.5],[-80 54 57.9],   'Ponce de Leon South' ; ...
      'K', '8722375', datenum(2011,05,01), now,                [27 09 55.0],[-80 09 46.0],   'St. Lucie Inlet' ; ...
      'L', '8720494', datenum(2011,05,01), now,                [29 59 41.0],[-81 19 46.4],   'Tolomato River' ; ...
      'M', '8722125', datenum(2011,05,01), now,                [27 37 55.4],[-80 22 16.4],   'Vero Beach' ; ...
      'N', '8720554', datenum(2011,05,01), now,                [29 54 59.7],[-81 17 59.9],   'Vilano Beach' ; ...
              };

  stcd = [];
  if ( ischar(stcd_or_stn) && numel(stcd_or_stn) == 1 )
    stcd = stcd_or_stn;
    stn.station_name = ['fdep',lower(stcd)];

  elseif ( isfield(stcd_or_stn,'station_name') )
    stn = stcd_or_stn;
    if ( numel(stn.station_name) == 5 && strncmpi(stn.station_name,'fdep',4) )
      stcd = stn.station_name(5);
    end;
  end;
  if ( isempty(stcd) )
    error('First arg must be a station code letter or station STRUCT with station_name beginning FDEP');
  end;

  stix = find(cellfun(@(x)(~isempty(x)),strfind(fdep_stns(:,1),upper(stcd))));
  if ( isempty(stix) )
    error('Unknown station code %s',stcd);
  end;
  stn.station_id = fdep_stns{stix,2};
  stn.lat = dms2degrees(fdep_stns{stix,5});
  stn.lon = dms2degrees(fdep_stns{stix,6});
  stn.long_name = fdep_stns{stix,7};

  dts = fdep_stns{stix,3}:fdep_stns{stix,4};

return;
