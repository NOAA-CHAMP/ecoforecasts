function stn = read_plsf1_usgs(stn)
%function stn = read_plsf1_usgs(stn)
% Pulaski Shoals underwater temperature data,,,
% "June 14, 2009 to July 21, 2010",,,
% "Collected by: Ilsa B. Kuffner, U.S. Geological Survey, St. Petersburg Coastal and Marine Science Center, 600 4th Street S., St. Petersburg, FL 33701",,,
% Phone: 727-803-8747,,,
% Email: ikuffner@usgs.gov,,,
% Equipment:  Onset HOBO Water Temp Pro V2 ,,,
% Water depth: 18',,,
% Latitude: N 24 degrees 41.613 minutes,,,
% Longitude: W 82 degrees 46.368 minutes,,,
% ,,,
% "Time, GMT-04:00",Date,Time,DRTO

  [d,t] = textread(fullfile(get_ecoforecasts_path('data'),'Pulaski Shoals temperature data.csv'),'%s%*s%*s%f','delimiter',',','headerlines',11);

  %Overlapping values
  d(26169:26172)=[]; t(26169:26172)=[];

  dts = datenum(d);

  if ( ~exist('stn','var') || isempty(stn) )
    stn.lat =  24.6935;
    stn.lon = -82.7728;
    stn.depth = 5.5;
  end;

  stn.usgs_seatemp.date = dts(:);
  stn.usgs_seatemp.data = t(:);


return;
