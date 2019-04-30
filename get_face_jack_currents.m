function jack = get_face_jack_currents(fname)
%function jack = get_face_jack_currents(fname)
%
% Read currents data from Jack Stamates' deployment of ADCPs and Tilt Current
% Meters (TCMs) in 2013-2015, particularly summers in 2014 and 2015. Looks
% for a MAT file first, and failing that calls READ_FACE_JACK_CURRENTS to
% reload the raw CSV file (see comments in that M-function).
%
% Last Saved Time-stamp: <Wed 2018-08-01 13:12:14 Eastern Daylight Time gramer>

  if ( ~exist('fname','var') || isempty(fname) )
    % If filename is needed, use default for READ_FACE_JACK_CURRENTS below
    fname = [];
  end;

  datapath = get_ecoforecasts_path('data');
  matfname = fullfile(datapath,'face_jack_currents.mat');

  if ( exist(matfname,'file') )
    disp(['Loading ',matfname]);
    load(matfname);

  else
    jack = read_face_jack_currents(fname);

    disp(['Saving ',matfname]);
    save(matfname,'jack','-v7.3');
  end;

  % Per Jack.Stamates@noaa.gov, 2017 Oct 6, "FACE ADCP bin depths.xlsx"
  jack.shallow.adcp_heights = jack.shallow.depth - ...
      [6.9,6.4,5.9,5.4,4.9,4.4,3.9,3.4,2.9,2.4,1.9,1.4,0.9,0.4,0];
  jack.deep.adcp_heights = jack.deep.depth - ...
      [21.27,20.27,19.27,18.27,17.27,16.27,15.27,14.27,13.27,12.27,11.27,10.27,9.27,8.27,7.27,6.27,5.27,4.27,3.27,2.27,1.27];

  % Having 'station_name' is useful for, e.g., modeling tides (STATION_TMD_TIDE)
  jack.shallow.station_name = 'face_shallow';
  jack.tcm1.station_name = 'face_tcm1';
  jack.tcm2.station_name = 'face_tcm2';
  jack.deep.station_name = 'face_deep';

return;
