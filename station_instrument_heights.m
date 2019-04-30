function [wz,az,pz,stz,dtz,slz,dlz] = station_instrument_heights(station_name)
%function [wz,az,pz,stz,dtz,slz,dlz] = station_instrument_heights(station_name)
%
% Get height/depth [m] of wind, air temperature, and air pressure sensors, of
% 'shallow' and 'deep' sea temperature + salinity package, and of 'shallow'
% and 'deep' light sensors (PAR and/or UV*) at station named 'station_name'.
%
% Last Saved Time-stamp: <Thu 2011-06-23 15:45:09  Lew.Gramer>

  % Stations are given below in reverse order of latitude/longshore distance

  % Source for most of the heights below are individual NDBC station pages, e.g.,
  %  http://www.ndbc.noaa.gov/station_page.php?station=mlrf1
  % ... and Measurement Descriptions at http://www.ndbc.noaa.gov/measdes.shtml

  switch (lower(station_name)),
   case 'lkwf1',
    d = [13.7,13.4,6.0,1,nan,nan,nan];
   case 'fwyf1',
    d = [43.9,11.0,29.3,2,nan,nan,nan];
    % Was 2.0: Increased to match GET_ALL_STATION_METADATA
    %d = [43.9,11.0,29.3,3,nan,nan,nan];
   case 'mlrf1',
    % d = [15.8,15.5,11.3,1,nan,nan,nan];
    % Was 1.0: Increased to match GET_ALL_STATION_METADATA
    % "stands in 9 feet (2.7 m) of water": http://en.wikipedia.org/wiki/Unmanned_reef_lights_of_the_Florida_Keys
    % d = [15.8,15.5,11.3,3.55,nan,nan,nan];
    d = [15.8,15.5,11.3,2.7,nan,nan,nan];
   case 'lonf1',
    %d = [7.0, 6.7, 5.8,1,nan,nan,nan];
    % Sea_T depth was 1.0! Reestimated from 1-D Wave i-depth...
    d = [ 7.0, 6.7, 5.8,1.28,nan,nan,nan];
   case 'tnrf1',
    % NUMBERS STOLEN FROM MLRF1: no metadata for TNRF1??
    % "49 feet high": http://www.uscg.mil/history/weblighthouses/LHFL.asp
    % "stands in 15 feet (4.6 m) of water": http://en.wikipedia.org/wiki/Unmanned_reef_lights_of_the_Florida_Keys
    d = [15.8,15.5,11.3,4.9,nan,nan,nan];
   case 'smkf1',
    d = [48.5,10.1,36.6,2,nan,nan,nan];
   case 'looe1',
    % NOTE: LOOE1 "borrows" meteorological data from SMKF1 for now
    d = [48.5,10.1,36.6,4.9,23.31,nan,nan];
    % d = [48.5,10.1,36.6,23.31,4.9,nan,nan];
   case 'sanf1',
    d = [45.4,12.8,38.1,1,nan,nan,nan];
   case 'plsf1',
    d = [17.7,17.4,15.8,1,nan,nan,nan];
   case 'dryf1',
    % NUMBERS STOLEN FROM PLSF1: no metadata for DRYF1??
    d = [17.7,17.4,15.8,1,nan,nan,nan];
   case '42003',
    % d = [5.0,4.0,0.0,0.6,nan,nan,nan];
    % NO?! See: http://www.ndbc.noaa.gov/data/stations/buoyht.txt
    d = [10.0,10.0,0.0,0.6,nan,nan,nan];
   case 'cmrc3',
    d = [5.5,5,5,2,5,2,5];
   case 'lsib4', %==cmrc3
    d = [5.5,5,5,2,5,2,5];
   case 'srvi2',
    d = [5.5,5,5,2,5,2,5];
   case 'lppr1',
    d = [5.5,5,5,2,5,2,5];
   case 'dbjm1',
    d = [5.5,5,5,2,5,2,5];
   case 'lciy2',
    d = [5.5,5,5,2,4.5,2,4.5];
   case 'llbp7',
    d = [5.5,5,5,2,5,2,5];


   %%%% HACK??? STOLEN FROM MLRF1
   case 'tavrk',
    d = [15.8,15.5,11.3,2.7,nan,nan,nan];
   case 'consh',
    d = [15.8,15.5,11.3,2.7,nan,nan,nan];
   case 'condp',
    d = [15.8,15.5,11.3,2.7,nan,nan,nan];

   %%%% HACK??? STOLEN FROM FWYF1
   case 'bnpin',
    d = [43.9,11.0,29.3,2,nan,nan,nan];
   case 'bnpmi',
    d = [43.9,11.0,29.3,2,nan,nan,nan];
   case 'bnppa',
    d = [43.9,11.0,29.3,2,nan,nan,nan];
   case 'bnpon',
    d = [43.9,11.0,29.3,2,nan,nan,nan];
   case 'bnpnn',
    d = [43.9,11.0,29.3,2,nan,nan,nan];

   otherwise,
    error('Have not yet coded sensor heights for station "%s"', station_name);
  end;

  wz=d(1); az=d(2); pz=d(3); stz=d(4); dtz=d(5); slz=d(6); dlz=d(7);

return;
