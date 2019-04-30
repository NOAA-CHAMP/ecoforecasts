function jack = read_face_jack_currents(fname)
%function jack = read_face_jack_currents(fname)
%
% Read currents data from Jack Stamates' deployment of ADCPs and Tilt Current
% Meters (TCMs) in 2013-2015, particularly the summers of 2014 and 2015. 
%
% Last Saved Time-stamp: <Mon 2016-08-29 17:30:16 Eastern Daylight Time Lew.Gramer>

  % HUGE XLS file which I turned into a STILL HUGE CSV file...
  facepath = get_coral_path('FACE');
  if ( ~exist('fname','var') || isempty(fname) )
    fname = fullfile(facepath,'currents_2015','Shallow and deep and TCM  merge ver 2  for export.csv');
  end;

  %% Read file
  x = fileread(fname);

  %% Header line, extra lines, and first few lines of real data appear as follows:
  %Date,deep missing data,JD (UT),Deep temp,Index,East 1,East 2,East 3,East 4,East 5,East 6,East 7,East 8,East 9,East 10,East 11,East 12,East 13,East 14,East 15,East 16,East 17,East 18,East 19,East 20,East 21,North 1,North 2,North 3,North 4,North 5,North 6,North 7,North 8,North 9,North 10,North 11,North 12,North 13,North 14,North 15,North 16,North 17,North 18,North 19,North 20,North 21,Vert 1,Vert 2,Vert3,vert 4,vert 5,vert 6,vert 7,vert 8,vert 9,vert 10,vert 11,vert 12,vert 13,vert 14,vert 15,vert 16,vert 17,vert 18,vert 19,vert 20,vert 21,mag 1,mag 2,mag 3,mag 4,mag 5,mag 6,mag 7,mag 8,mag 9,mag 10,mag 11,mag 12,mag 13,mag 14,mag 15,mag 16,mag 17,mag 18,mag 19,mag 20,mag 21,dir 1,dir 2,dir 3,dir 4,dir 5,dir 6,dir 7,dir 8,dir 9,dir 10,dir 11,dir 12,dir 13,dir 14,dir 15,dir 16,dir 17,dir 18,dir 19,dir 20,dir 21,grid time,Shallow temp,Z,E1,E2,E3,E4,E5,E6,E7,E8,E9,E10,E11,E12,E13,E14,E15,N1,N2,N3,N4,N5,N6,N7,N8,N9,N10,N11,N12,N13,N14,N15,VV1,VV2,VV3,VV4,VV5,VV6,VV7,VV8,VV9,VV10,VV11,VV12,VV13,VV14,VV15,M1,M2,M3,M4,M5,M6,M7,M8,M9,M10,M11,M12,M13,M14,M15,D1,D2,D3,D4,D5,D6,D7,D8,D9,D10,D11,D12,D13,D14,D15,NewVar1,Deep East WCA,Deep North WCA,Shallow East WCA,Shallow North WCA,S has data,D has data,S and D have data,TCM has data,Shallow ADCP WCA is north,TCM 1 is north,TCM2 is north,Deep ADCP WCA  is north,NewVar15,Date,Time,T2 Speed (cm/s),T2 Bearing (degrees),T2 Velocity-N (cm/s),T2 Velocity-E (cm/s),T2   Temperature (C),T1  Speed (cm/s),T1 Bearing (degrees),T1 Vel N,T1  Velocity-E (cm/s),T1   Temperature (C)
  %11/12/13,,,,1,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,11/12/13 0:00,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,0,0,0,0,,,,,,,,,,,,,,,,,
  %11/12/13,,,,2,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,11/12/13 0:20,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,0,0,0,0,,,,,,,,,,,,,,,,,

  % (Followed by 42 similar lines...)

  %11/12/13,,11/12/13 14:20,27.08,44,-16,4,-5,-15,-17,-38,26,37,3,-17,-10,33,4,-9,14,-15,-25,-13,-65,137,-208,70,134,163,150,154,242,177,205,215,190,198,192,234,170,317,112,123,46,30,-154,-307,23,44,41,64,48,32,48,41,35,33,30,26,28,24,-46,12,17,55,49,-20,31,72,134,163,151,155,245,179,208,215,191,198,195,234,170,317,113,126,48,72,206,371,347.1,1.7,358.2,354.3,353.7,351.1,8.4,10.2,0.8,354.9,357.1,9.8,1,357,2.5,352.4,348.5,344.2,294.8,138.3,214.1,11/12/13 14:20,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,-3.277777778,171.7777778,,,0,1,0,0,,,,1,,,,,,,,,,,,,
  %11/12/13,,11/12/13 14:40,27.2,45,-17,-4,-12,-11,-5,0,-4,7,34,31,3,20,42,45,12,15,-9,1,-119,-21,-420,144,145,153,167,159,163,164,142,144,151,155,157,151,223,134,160,133,134,40,-264,-445,15,10,16,8,2,10,16,10,9,8,9,1,3,17,4,8,8,-5,-5,0,48,145,145,153,167,159,163,164,142,148,154,155,158,157,227,135,161,133,134,126,265,612,353.3,358.4,355.5,356.2,358.2,0,358.6,2.8,13.3,11.6,1.1,7.3,15.5,11.4,5.1,5.4,356.1,0.4,288.6,184.5,223.3,11/12/13 14:40,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,8.222222222,154.3888889,,,0,1,0,0,,,,1,,,,,,,,,,,,,
  %11/12/13,,11/12/13 15:00,27.19,46,14,11,31,-8,8,18,8,15,44,62,26,32,62,52,52,31,17,9,-72,-123,-85,160,196,188,238,211,215,238,248,243,253,305,298,316,302,316,313,342,313,194,-32,-234,5,8,3,2,-3,6,0,-4,0,1,3,2,-1,-3,-8,-8,-10,-13,-36,19,-13,161,196,191,238,211,216,238,248,247,260,306,300,322,306,320,315,342,313,207,127,249,5,3.2,9.4,358.1,2.2,4.8,1.9,3.5,10.3,13.8,4.9,6.1,11.1,9.8,9.3,5.7,2.8,1.6,339.6,255.4,200,11/12/13 15:00,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,26.88888889,260.8333333,,,0,1,0,0,,,,1,,,,,,,,,,,,,
  %11/12/13,,11/12/13 15:20,27.24,47,17,37,31,30,15,30,0,24,7,-12,31,16,-8,1,-17,-14,18,34,-65,-155,-118,163,174,196,162,167,174,198,204,232,237,228,250,224,266,273,233,201,143,-27,-167,-152,1,5,-2,-4,6,-3,4,9,2,10,-7,5,7,8,10,1,-12,-25,-19,37,-11,164,178,198,165,168,177,198,205,232,237,230,251,224,266,274,233,202,147,70,228,192,6,12,9,10.5,5.1,9.8,0,6.7,1.7,357.1,7.7,3.7,358,0.2,356.4,356.6,5.1,13.4,247.4,222.9,217.8,11/12/13 15:20,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,13.33333333,206.9444444,,,0,1,0,0,,,,1,,,,,,,,,,,,,
  %11/12/13,,11/12/13 15:40,27.29,48,26,-1,12,-4,-7,27,3,6,36,30,5,4,-46,-6,5,31,-77,-103,-182,-213,-517,163,165,155,186,208,196,196,211,219,232,233,252,272,268,275,268,203,182,71,-238,-472,-4,7,4,4,8,1,4,4,0,-11,7,-9,-11,-3,-1,-13,4,0,4,59,69,165,165,155,186,208,198,196,211,222,234,233,252,276,268,275,270,217,209,195,319,700,9.1,359.7,4.4,358.8,358.1,7.8,0.9,1.6,9.3,7.4,1.2,0.9,350.4,358.7,1,6.6,339.2,330.5,291.3,221.8,227.6,11/12/13 15:40,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,-3.277777778,215.7777778,,,0,1,0,0,,,,1,,,,,,,,,,,,,

  % (Followed by 52521 similar lines...)


  %% Parse file

  deep_abins = 21;
  shallow_abins = 15;

  %fmt = repmat('%s',[1,214]);
  fmt = make_format_string({2,'',1,'%s',1,'%f',1,'',5*deep_abins,'%f',1,'',(2+(5*shallow_abins)),'%f',16,'',2*5,'%f'},'','');
  C = textscan(x,fmt,'HeaderLines',1,'EndOfLine','\r\n','Delimiter',',');

  x=[]; clear x


  %% Process file

  badrow = find(cellfun(@isempty,C{1}));
  for ix = 1:numel(C)
    C{ix}(badrow) = [];
  end;

  dts = datenum(C{1},'mm/dd/yy HH:MM');

  deep_tix = 2;
  deep_aix = deep_tix + 1;
  shallow_tix = deep_aix + (5*deep_abins);
  shallow_zix = shallow_tix + 1;
  shallow_aix = shallow_zix + 1;
  tcm2_aix = shallow_aix + (5*shallow_abins);
  tcm1_aix = tcm2_aix + 5;


  % AOML  Hollywood shallow: 
  %  26 deg ,   01.149 min N      
  %  080 deg,  06.518 W           
  %  Z = 7.3m  
  % AOML   Hollywood Deep
  %  26 1.1075'N, 80 5.1716'W  
  %  Z=26.2 m   
  % TCM 1
  %  26 deg 1.103 min N   , 80 deg 5.878 min w
  %  Z=11m
  % TCM 2 
  %  26 deg 1.132 min N, 80 deg 5.545 min W
  %  Z=21 m 

  jack=[]; clear jack
  % AOML  Hollywood shallow: 
  jack.shallow.lon = dm2degrees([-080 06.518]);
  jack.shallow.lat = dm2degrees([  26 01.149]);
  jack.shallow.depth = 7.3;
  % AOML   Hollywood Deep
  jack.deep.lon = dm2degrees([-080 05.1716]);
  jack.deep.lat = dm2degrees([  26 01.1075]);
  jack.deep.depth = 26.2;
  % TCM 1
  jack.tcm1.lon = dm2degrees([-080 05.878]);
  jack.tcm1.lat = dm2degrees([  26 01.103]);
  jack.tcm1.depth = 11;
  % TCM 2
  jack.tcm2.lon = dm2degrees([-080 05.545]);
  jack.tcm2.lat = dm2degrees([  26 01.132]);
  jack.tcm2.depth = 21;


  %% Sea temperatures - limit each time series to periods when that sensor
  %  may actually have been deployed 
  cix = shallow_tix;
  goodix = find(~isnan(C{cix}),1):find(~isnan(C{cix}),1,'last'); n = numel(goodix);
  jack.shallow.seatemp.date = dts(goodix);
  jack.shallow.seatemp.data = C{cix}(goodix);

  cix = shallow_zix;
  goodix = find(~isnan(C{cix}),1):find(~isnan(C{cix}),1,'last'); n = numel(goodix);
  jack.shallow.depths.date = dts(goodix);
  jack.shallow.depths.data = C{cix}(goodix);

  cix = tcm1_aix+4;
  goodix = find(~isnan(C{cix}),1):find(~isnan(C{cix}),1,'last'); n = numel(goodix);
  jack.tcm1.seatemp.date = dts(goodix);
  jack.tcm1.seatemp.data = C{cix}(goodix);

  cix = tcm2_aix+4;
  goodix = find(~isnan(C{cix}),1):find(~isnan(C{cix}),1,'last'); n = numel(goodix);
  jack.tcm2.seatemp.date = dts(goodix);
  jack.tcm2.seatemp.data = C{cix}(goodix);

  cix = deep_tix;
  goodix = find(~isnan(C{cix}),1):find(~isnan(C{cix}),1,'last'); n = numel(goodix);
  jack.deep.seatemp.date = dts(goodix);
  jack.deep.seatemp.data = C{cix}(goodix);


  % ADCP velocity and magnitude data are in mm/sec. Direction is in degrees.   The ADCP data is compensated for a    -6deg ( 6 deg west)  magnetic offset.   The TCM are in cm/sec and are not compensated for mag var. 

  %% Tilt Current Meter 1 - only periods this instrument may actually have been deployed
  cix = tcm1_aix;
  goodix = find(~isnan(C{cix}),1):find(~isnan(C{cix}),1,'last'); n = numel(goodix);
  jack.tcm1.spd.date = dts(goodix);
  jack.tcm1.spd.data = C{cix}(goodix) ./ 100;
  jack.tcm1.raw_dir.date = dts(goodix);
  jack.tcm1.raw_dir.data = C{cix+1}(goodix);
  jack.tcm1.raw_u.date = dts(goodix);
  jack.tcm1.raw_u.data = C{cix+3}(goodix) ./ 100;
  jack.tcm1.raw_v.date = dts(goodix);
  jack.tcm1.raw_v.data = C{cix+2}(goodix) ./ 100;
  % Compensate for local magnetic variation
  jack.tcm1.dir.date = dts(goodix);
  %jack.tcm1.dir.data = jack.tcm1.raw_dir.data + 6;
  jack.tcm1.dir.data = jack.tcm1.raw_dir.data - 6;
  jack.tcm1.u.date = dts(goodix);
  jack.tcm1.v.date = dts(goodix);
  [jack.tcm1.u.data,jack.tcm1.v.data] = spddir_to_uv_curr(jack.tcm1.spd.data,jack.tcm1.dir.data);

  %% Tilt Current Meter 2 - only periods this instrument may actually have been deployed
  cix = tcm2_aix;
  goodix = find(~isnan(C{cix}),1):find(~isnan(C{cix}),1,'last'); n = numel(goodix);
  jack.tcm2.spd.date = dts(goodix);
  jack.tcm2.spd.data = C{cix}(goodix) ./ 100;
  jack.tcm2.raw_dir.date = dts(goodix);
  jack.tcm2.raw_dir.data = C{cix+1}(goodix);
  jack.tcm2.raw_u.date = dts(goodix);
  jack.tcm2.raw_u.data = C{cix+3}(goodix) ./ 100;
  jack.tcm2.raw_v.date = dts(goodix);
  jack.tcm2.raw_v.data = C{cix+2}(goodix) ./ 100;
  % Compensate for local magnetic variation
  jack.tcm2.dir.date = dts(goodix);
  %jack.tcm2.dir.data = jack.tcm2.raw_dir.data + 6;
  jack.tcm2.dir.data = jack.tcm2.raw_dir.data - 6;
  jack.tcm2.u.date = dts(goodix);
  jack.tcm2.v.date = dts(goodix);
  [jack.tcm2.u.data,jack.tcm2.v.data] = spddir_to_uv_curr(jack.tcm2.spd.data,jack.tcm2.dir.data);


  % IMPORTANT note.    ADCP data becomes biased near the surface due to acoustic side lobe interference.  
  % DO NOT use  shallow ADCP data above bin 10  and deep ADCP data above bin 18.
  shallow_top_good_bin = 10;
  deep_top_good_bin = 18;

  %% Shallow ADCP - only periods this instrument may actually have been deployed
  cix = shallow_aix;
  goodix = find(~isnan(C{cix}),1):find(~isnan(C{cix}),1,'last'); n = numel(goodix);
  shallow_sfc_bins = shallow_top_good_bin-2:shallow_top_good_bin;
  shallow_mid_bins = 5:7;
  shallow_btm_bins = 1:4;

  % Shallow ADCP - current speed
  jack.shallow.adcp_spd.date = dts(goodix);
  jack.shallow.adcp_spd.data = [];
  jack.shallow.adcp_spd.prof = repmat(nan,[n,shallow_abins]);
  for ix = 1:shallow_top_good_bin
    jack.shallow.adcp_spd.prof(:,ix) = C{cix+(shallow_abins*3)+ix-1}(goodix) ./ 1000;
  end;
  jack.shallow.adcp_spd.data = nanmean(jack.shallow.adcp_spd.prof,2);

  % Shallow ADCP - current direction
  jack.shallow.adcp_dir.date = dts(goodix);
  jack.shallow.adcp_dir.data = [];
  jack.shallow.adcp_dir.prof = repmat(nan,[n,shallow_abins]);
  for ix = 1:shallow_top_good_bin
    jack.shallow.adcp_dir.prof(:,ix) = C{cix+(shallow_abins*4)+ix-1}(goodix);
  end;

  % Shallow ADCP - current East velocity
  jack.shallow.adcp_u.date = dts(goodix);
  jack.shallow.adcp_u.data = [];
  jack.shallow.adcp_u.prof = repmat(nan,[n,shallow_abins]);
  for ix = 1:shallow_top_good_bin
    jack.shallow.adcp_u.prof(:,ix) = C{cix+ix-1}(goodix) ./ 1000;
  end;
  jack.shallow.adcp_u.data = nanmean(jack.shallow.adcp_u.prof,2);

  % Shallow ADCP - current North velocity
  jack.shallow.adcp_v.date = dts(goodix);
  jack.shallow.adcp_v.data = [];
  jack.shallow.adcp_v.prof = repmat(nan,[n,shallow_abins]);
  for ix = 1:shallow_top_good_bin
    jack.shallow.adcp_v.prof(:,ix) = C{cix+(shallow_abins*1)+ix-1}(goodix) ./ 1000;
  end;
  jack.shallow.adcp_v.data = nanmean(jack.shallow.adcp_v.prof,2);

  % Shallow ADCP - current vertical velocity
  jack.shallow.adcp_w.date = dts(goodix);
  jack.shallow.adcp_w.data = [];
  jack.shallow.adcp_w.prof = repmat(nan,[n,shallow_abins]);
  for ix = 1:shallow_top_good_bin
    jack.shallow.adcp_w.prof(:,ix) = C{cix+(shallow_abins*2)+ix-1}(goodix) ./ 1000;
  end;
  jack.shallow.adcp_w.data = nanmean(jack.shallow.adcp_w.prof,2);

  jack.shallow.adcp_dir.data = uv_to_dir_curr(jack.shallow.adcp_u.data,jack.shallow.adcp_v.data);

  % Shallow ADCP - Mid-column averages
  jack.shallow.adcp_spd_sfc.date = dts(goodix);
  jack.shallow.adcp_spd_sfc.data = nanmean(jack.shallow.adcp_spd.prof(:,shallow_sfc_bins),2);
  jack.shallow.adcp_spd_mid.date = dts(goodix);
  jack.shallow.adcp_spd_mid.data = nanmean(jack.shallow.adcp_spd.prof(:,shallow_mid_bins),2);
  jack.shallow.adcp_spd_btm.date = dts(goodix);
  jack.shallow.adcp_spd_btm.data = nanmean(jack.shallow.adcp_spd.prof(:,shallow_btm_bins),2);

  jack.shallow.adcp_u_sfc.date = dts(goodix);
  jack.shallow.adcp_u_sfc.data = nanmean(jack.shallow.adcp_u.prof(:,shallow_sfc_bins),2);
  jack.shallow.adcp_u_mid.date = dts(goodix);
  jack.shallow.adcp_u_mid.data = nanmean(jack.shallow.adcp_u.prof(:,shallow_mid_bins),2);
  jack.shallow.adcp_u_btm.date = dts(goodix);
  jack.shallow.adcp_u_btm.data = nanmean(jack.shallow.adcp_u.prof(:,shallow_btm_bins),2);

  jack.shallow.adcp_v_sfc.date = dts(goodix);
  jack.shallow.adcp_v_sfc.data = nanmean(jack.shallow.adcp_v.prof(:,shallow_sfc_bins),2);
  jack.shallow.adcp_v_mid.date = dts(goodix);
  jack.shallow.adcp_v_mid.data = nanmean(jack.shallow.adcp_v.prof(:,shallow_mid_bins),2);
  jack.shallow.adcp_v_btm.date = dts(goodix);
  jack.shallow.adcp_v_btm.data = nanmean(jack.shallow.adcp_v.prof(:,shallow_btm_bins),2);

  jack.shallow.adcp_dir_sfc.date = dts(goodix);
  jack.shallow.adcp_dir_sfc.data = uv_to_dir_curr(jack.shallow.adcp_u_sfc.data,jack.shallow.adcp_v_sfc.data);
  jack.shallow.adcp_dir_mid.date = dts(goodix);
  jack.shallow.adcp_dir_mid.data = uv_to_dir_curr(jack.shallow.adcp_u_mid.data,jack.shallow.adcp_v_mid.data);
  jack.shallow.adcp_dir_btm.date = dts(goodix);
  jack.shallow.adcp_dir_btm.data = uv_to_dir_curr(jack.shallow.adcp_u_btm.data,jack.shallow.adcp_v_btm.data);



  %% Deep ADCP - only periods this instrument may actually have been deployed
  cix = deep_aix;
  goodix = find(~isnan(C{cix}),1):find(~isnan(C{cix}),1,'last'); n = numel(goodix);
  deep_sfc_bins = deep_top_good_bin-5:deep_top_good_bin;
  deep_mid_bins = 7:12;
  deep_btm_bins = 1:6;

  % Deep ADCP - current speed
  jack.deep.adcp_spd.date = dts(goodix);
  jack.deep.adcp_spd.data = [];
  jack.deep.adcp_spd.prof = repmat(nan,[n,deep_abins]);
  for ix = 1:deep_top_good_bin
    jack.deep.adcp_spd.prof(:,ix) = C{cix+(deep_abins*3)+ix-1}(goodix) ./ 1000;
  end;
  jack.deep.adcp_spd.data = nanmean(jack.deep.adcp_spd.prof,2);

  % Deep ADCP - current direction
  jack.deep.adcp_dir.date = dts(goodix);
  jack.deep.adcp_dir.data = [];
  jack.deep.adcp_dir.prof = repmat(nan,[n,deep_abins]);
  for ix = 1:deep_top_good_bin
    jack.deep.adcp_dir.prof(:,ix) = C{cix+(deep_abins*4)+ix-1}(goodix);
  end;

  % Deep ADCP - current East velocity
  jack.deep.adcp_u.date = dts(goodix);
  jack.deep.adcp_u.data = [];
  jack.deep.adcp_u.prof = repmat(nan,[n,deep_abins]);
  for ix = 1:deep_top_good_bin
    jack.deep.adcp_u.prof(:,ix) = C{cix+ix-1}(goodix) ./ 1000;
  end;
  jack.deep.adcp_u.data = nanmean(jack.deep.adcp_u.prof,2);

  % Deep ADCP - current North velocity
  jack.deep.adcp_v.date = dts(goodix);
  jack.deep.adcp_v.data = [];
  jack.deep.adcp_v.prof = repmat(nan,[n,deep_abins]);
  for ix = 1:deep_top_good_bin
    jack.deep.adcp_v.prof(:,ix) = C{cix+(deep_abins*1)+ix-1}(goodix) ./ 1000;
  end;
  jack.deep.adcp_v.data = nanmean(jack.deep.adcp_v.prof,2);

  % Deep ADCP - current vertical velocity
  jack.deep.adcp_w.date = dts(goodix);
  jack.deep.adcp_w.data = [];
  jack.deep.adcp_w.prof = repmat(nan,[n,deep_abins]);
  for ix = 1:deep_top_good_bin
    jack.deep.adcp_w.prof(:,ix) = C{cix+(deep_abins*2)+ix-1}(goodix) ./ 1000;
  end;
  jack.deep.adcp_w.data = nanmean(jack.deep.adcp_w.prof,2);

  jack.deep.adcp_dir.data = uv_to_dir_curr(jack.deep.adcp_u.data,jack.deep.adcp_v.data);

  % Deep ADCP - Mid-column averages
  jack.deep.adcp_spd_sfc.date = dts(goodix);
  jack.deep.adcp_spd_sfc.data = nanmean(jack.deep.adcp_spd.prof(:,deep_sfc_bins),2);
  jack.deep.adcp_spd_mid.date = dts(goodix);
  jack.deep.adcp_spd_mid.data = nanmean(jack.deep.adcp_spd.prof(:,deep_mid_bins),2);
  jack.deep.adcp_spd_btm.date = dts(goodix);
  jack.deep.adcp_spd_btm.data = nanmean(jack.deep.adcp_spd.prof(:,deep_btm_bins),2);

  jack.deep.adcp_u_sfc.date = dts(goodix);
  jack.deep.adcp_u_sfc.data = nanmean(jack.deep.adcp_u.prof(:,deep_sfc_bins),2);
  jack.deep.adcp_u_mid.date = dts(goodix);
  jack.deep.adcp_u_mid.data = nanmean(jack.deep.adcp_u.prof(:,deep_mid_bins),2);
  jack.deep.adcp_u_btm.date = dts(goodix);
  jack.deep.adcp_u_btm.data = nanmean(jack.deep.adcp_u.prof(:,deep_btm_bins),2);

  jack.deep.adcp_v_sfc.date = dts(goodix);
  jack.deep.adcp_v_sfc.data = nanmean(jack.deep.adcp_v.prof(:,deep_sfc_bins),2);
  jack.deep.adcp_v_mid.date = dts(goodix);
  jack.deep.adcp_v_mid.data = nanmean(jack.deep.adcp_v.prof(:,deep_mid_bins),2);
  jack.deep.adcp_v_btm.date = dts(goodix);
  jack.deep.adcp_v_btm.data = nanmean(jack.deep.adcp_v.prof(:,deep_btm_bins),2);

  jack.deep.adcp_dir_sfc.date = dts(goodix);
  jack.deep.adcp_dir_sfc.data = uv_to_dir_curr(jack.deep.adcp_u_sfc.data,jack.deep.adcp_v_sfc.data);
  jack.deep.adcp_dir_mid.date = dts(goodix);
  jack.deep.adcp_dir_mid.data = uv_to_dir_curr(jack.deep.adcp_u_mid.data,jack.deep.adcp_v_mid.data);
  jack.deep.adcp_dir_btm.date = dts(goodix);
  jack.deep.adcp_dir_btm.data = uv_to_dir_curr(jack.deep.adcp_u_btm.data,jack.deep.adcp_v_btm.data);


  C=[]; clear C

return;
