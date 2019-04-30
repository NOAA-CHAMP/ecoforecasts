function face = get_FACE_Hazen_currents(face)
%function face = get_FACE_Hazen_currents([face])
%
% Read data processed by Jack Stamates from ADCPs deployed by Hazen & Sawyer
% and "AOML" (FACE) near Hollywood and Broward waste water outfalls (around
% 27 m isobath), and directly inshore of them (around 7 m isobath).
%
% Last Saved Time-stamp: <Sat 2018-03-24 14:38:33 Eastern Daylight Time gramer>
  
  if ( ~exist('face','var') )
    face = [];
  end;
  
  datapath = get_ecoforecasts_path('data');
  
  %adcppath = 'd:/coral/FACE/ADCP/Bottom mounted ADCP data/FACE HW and BR  2011-2012 Master ADCP data set';
  %adcppath = '/Users/lew.gramer/Documents/coral/FACE/ADCP/Bottom mounted ADCP data/FACE HW and BR  2011-2012 Master ADCP data set';
  adcppath = fullfile(get_coral_path,'FACE','ADCP','Bottom mounted ADCP data','FACE HW and BR  2011-2012 Master ADCP data set');
  

  matfname = fullfile(datapath,'FACE_Hazen_currents.mat');
  
  if ( exist(matfname,'file') )
    disp(['Loading ',matfname]);
    x = load(matfname);
    
    flds = fieldnames(x);
    for fldix = 1:numel(flds)
      fld = flds{fldix};
      if ( isfield(face,fld) )
        error('Field already exists! FACE.%s',fld);
      end;
      face.(fld) = x.(fld);
    end;
    x=[];
    
  
  else
    
    adcpfname = fullfile(adcppath,'ADCP_Master_21aug13.xlsx');
    
    disp(['Extracting ',adcpfname]);
    
    face.AOML_HW.lon = -80.10863; face.AOML_HW.lat = 26.019150; face.AOML_HW.depth = 7.3;
    face.HanS_HW.lon = -80.08619; face.HanS_HW.lat = 26.018458; face.HanS_HW.depth = 26.2;
    face.AOML_BR.lon = -80.07613; face.AOML_BR.lat = 26.264783; face.AOML_BR.depth = 8.2;
    face.HanS_BR.lon = -80.06253; face.HanS_BR.lat = 26.250000; face.HanS_BR.depth = 32.4;
    
    x = importdata(adcpfname);
    
    %x.textdata.data(2,:),
    dat = x.data.data;
    dat(dat==99999) = nan;


    %DEBUG:
    disp('Dates');
    dts = datenum(dat(2:end,1));

    % The first works on MUIRGEN, the second on MANANANN??
    if ( all(dts > datenum(100,1,1) & dts < datenum(1900,1,1)) )
      error('This code cannot be run on the Mac right now!');
      dts = dts+datenum(1900,0,0, 0,0,0);
    elseif ( all(dts < datenum(100,1,1)) )
      dts = dts+datenum(2010,1,0, 0,0,0);
    else
      error('Could not figure out correct Date format!');
    end;

    %DEBUG:
    disp('Hydrography');
    face.HanS_BR.adcp_seatemp.date = dts;
    face.HanS_BR.adcp_seatemp.data = dat(2:end,9);
    
    face.HanS_BR.adcp_depth.date = dts;
    face.HanS_BR.adcp_depth.data = dat(2:end,10);
    % For some reason, Jack switched units to centibars on this mid-series??
    face.HanS_BR.adcp_depth.data(face.HanS_BR.adcp_depth.data>100) = ...
        face.HanS_BR.adcp_depth.data(face.HanS_BR.adcp_depth.data>100)./10;
    
    face.HanS_HW.adcp_seatemp.date = dts;
    face.HanS_HW.adcp_seatemp.data = dat(2:end,11);
    
    face.HanS_HW.adcp_depth.date = dts;
    face.HanS_HW.adcp_depth.data = dat(2:end,12);
    
    face.AOML_BR.adcp_seatemp.date = dts;
    face.AOML_BR.adcp_seatemp.data = dat(2:end,13);
    
    face.AOML_BR.adcp_depth.date = dts;
    face.AOML_BR.adcp_depth.data = dat(2:end,14);
    
    face.AOML_HW.adcp_seatemp.date = dts;
    face.AOML_HW.adcp_seatemp.data = dat(2:end,15);
    
    face.AOML_HW.adcp_depth.date = dts;
    face.AOML_HW.adcp_depth.data = dat(2:end,16);
    
    
    %DEBUG:
    disp('Currents');

    % FROM "Meta data for FACE AOML and Hazen and Sawyer current  data set  5 20 2011 .docx":

    %  First 15 depths (although not every instrument has 15 valid columns in spreadsheet)

    face.HanS_HW.adcp_depth_correction = 0;
    face.HanS_HW.adcp_bin_depths = -[23,21,19,17,15,13,11,9,7,5,3,1,-1,-3,-5,] + face.HanS_HW.adcp_depth_correction;

    face.HanS_BR.adcp_depth_correction = 0;
    face.HanS_BR.adcp_bin_depths = -[28.65,26.65,24.65,22.65,20.65,18.65,16.65,14.65,12.65,10.65,8.65,6.65,4.65,2.65,0.65,] + face.HanS_BR.adcp_depth_correction;

    % Surface intensification spikes look like they are actually consistent
    % with real surface height above bottom, so... allow "depths" to top 0 
    %face.AOML_HW.adcp_depth_correction = -0.5;
    face.AOML_HW.adcp_depth_correction = 0;
    face.AOML_HW.adcp_bin_depths = -[6.9,6.4,5.9,5.4,4.9,4.4,3.9,3.4,2.9,2.4,1.9,1.4,0.9,0.4,-0.1,] + face.AOML_HW.adcp_depth_correction;;

    %face.AOML_BR.adcp_depth_correction = -1;
    face.AOML_BR.adcp_depth_correction = 0;
    face.AOML_BR.adcp_bin_depths = -[6.35,5.85,5.35,4.85,4.35,3.85,3.35,2.85,2.35,1.85,1.35,0.85,0.35,-0.15,-0.65,] + face.AOML_BR.adcp_depth_correction;


    % Column ranges with (potentially) valid data in spreadsheet
    face.HanS_HW.valid_bins = 1:14;
    face.HanS_BR.valid_bins = 1:15;
    face.AOML_HW.valid_bins = 1:15;
    face.AOML_BR.valid_bins = 1:15;


    % Cols 24-38: HZ HW E B1-B15
    face.HanS_HW.adcp_u.date = dts;
    face.HanS_HW.adcp_u.prof = dat(2:end,24+face.HanS_HW.valid_bins-1) ./ 1000; %mm/s -> m/s
    face.HanS_HW.adcp_u.data = nanmean(face.HanS_HW.adcp_u.prof,2);
    face.HanS_HW.adcp_u.depths = face.HanS_HW.adcp_bin_depths(face.HanS_HW.valid_bins);

    % Cols 39-53: HZ HW N B1-B15
    face.HanS_HW.adcp_v.date = dts;
    face.HanS_HW.adcp_v.prof = dat(2:end,39+face.HanS_HW.valid_bins-1) ./ 1000; %mm/s -> m/s
    face.HanS_HW.adcp_v.data = nanmean(face.HanS_HW.adcp_v.prof,2);
    face.HanS_HW.adcp_v.depths = face.HanS_HW.adcp_bin_depths(face.HanS_HW.valid_bins);

    % Cols 54-68: HZ HW Mag B1-B15
    % Cols 69-83: HZ HW Dir B1-B15


    % Cols 84-98: HZ BR E B1-B15
    face.HanS_BR.adcp_u.date = dts;
    face.HanS_BR.adcp_u.prof = dat(2:end,84+face.HanS_BR.valid_bins-1) ./ 1000; %mm/s -> m/s
    face.HanS_BR.adcp_u.data = nanmean(face.HanS_BR.adcp_u.prof,2);
    face.HanS_BR.adcp_u.depths = face.HanS_BR.adcp_bin_depths(face.HanS_BR.valid_bins);

    % Cols 99-113: HZ BR N B1-B15
    face.HanS_BR.adcp_v.date = dts;
    face.HanS_BR.adcp_v.prof = dat(2:end,99+face.HanS_BR.valid_bins-1) ./ 1000; %mm/s -> m/s
    face.HanS_BR.adcp_v.data = nanmean(face.HanS_BR.adcp_v.prof,2);
    face.HanS_BR.adcp_v.depths = face.HanS_BR.adcp_bin_depths(face.HanS_BR.valid_bins);
    
    % Cols 114-128: HZ BR Mag B1-B15
    % Cols 129-143: HZ BR Dir B1-B15

    
    % Cols 144-158: AOML HW E B1-B15
    face.AOML_HW.adcp_u.date = dts;
    face.AOML_HW.adcp_u.prof = dat(2:end,144+face.AOML_HW.valid_bins-1) ./ 1000; %mm/s -> m/s
    face.AOML_HW.adcp_u.data = nanmean(face.AOML_HW.adcp_u.prof,2);
    face.AOML_HW.adcp_u.depths = face.AOML_HW.adcp_bin_depths(face.AOML_HW.valid_bins);
    
    % Cols 159-173: AOML HW N B1-B15
    face.AOML_HW.adcp_v.date = dts;
    face.AOML_HW.adcp_v.prof = dat(2:end,159+face.AOML_HW.valid_bins-1) ./ 1000; %mm/s -> m/s
    face.AOML_HW.adcp_v.data = nanmean(face.AOML_HW.adcp_v.prof,2);
    face.AOML_HW.adcp_v.depths = face.AOML_HW.adcp_bin_depths(face.AOML_HW.valid_bins);
    
    % Cols 174-188: AOML HW Mag B1-B15
    % Cols 189-203: AOML HW Dir B1-B15


    % Cols 204-218: AOML BR E B1-B15
    face.AOML_BR.adcp_u.date = dts;
    face.AOML_BR.adcp_u.prof = dat(2:end,204+face.AOML_BR.valid_bins-1) ./ 1000; %mm/s -> m/s
    face.AOML_BR.adcp_u.data = nanmean(face.AOML_BR.adcp_u.prof,2);
    face.AOML_BR.adcp_u.depths = face.AOML_BR.adcp_bin_depths(face.AOML_BR.valid_bins);

    % Cols 219-233: AOML BR N B1-B15
    face.AOML_BR.adcp_v.date = dts;
    face.AOML_BR.adcp_v.prof = dat(2:end,219+face.AOML_BR.valid_bins-1) ./ 1000; %mm/s -> m/s
    face.AOML_BR.adcp_v.data = nanmean(face.AOML_BR.adcp_v.prof,2);
    face.AOML_BR.adcp_v.depths = face.AOML_BR.adcp_bin_depths(face.AOML_BR.valid_bins);

    % Cols 234-248: AOML BR Mag B1-B15
    % Cols 249-263: AOML BR Dir B1-B15

    disp(['Saving ',matfname]);
    save(matfname,'-struct','face');

  end; %if ( exist(matfname,'file') )

  if ( ~isfield(face.HanS_HW,'adcp_sfc_height') )
    face.HanS_HW = postprocess_adcp_uv(face.HanS_HW,'','',true);
    face.HanS_BR = postprocess_adcp_uv(face.HanS_BR,'','',true);
    face.AOML_HW = postprocess_adcp_uv(face.AOML_HW,'','',true);
    face.AOML_BR = postprocess_adcp_uv(face.AOML_BR,'','',true);

    % % That was a GRIDDED dataset: we are missing data for many date-times
    % face.HanS_BR.adcp_seatemp = get_FACE_Hazen_currents_filter(face.HanS_BR.adcp_seatemp);
    % face.HanS_HW.adcp_seatemp = get_FACE_Hazen_currents_filter(face.HanS_HW.adcp_seatemp);
    % face.AOML_BR.adcp_seatemp = get_FACE_Hazen_currents_filter(face.AOML_BR.adcp_seatemp);
    % face.AOML_HW.adcp_seatemp = get_FACE_Hazen_currents_filter(face.AOML_HW.adcp_seatemp);
    % face.HanS_BR.adcp_depth = get_FACE_Hazen_currents_filter(face.HanS_BR.adcp_depth);
    % face.HanS_HW.adcp_depth = get_FACE_Hazen_currents_filter(face.HanS_HW.adcp_depth);
    % face.AOML_BR.adcp_depth = get_FACE_Hazen_currents_filter(face.AOML_BR.adcp_depth);
    % face.AOML_HW.adcp_depth = get_FACE_Hazen_currents_filter(face.AOML_HW.adcp_depth);
    %
    % face.HanS_HW.adcp_u = get_FACE_Hazen_currents_filter(face.HanS_HW.adcp_u);
    % face.HanS_HW.adcp_v = get_FACE_Hazen_currents_filter(face.HanS_HW.adcp_v);
    % face.HanS_BR.adcp_u = get_FACE_Hazen_currents_filter(face.HanS_BR.adcp_u);
    % face.HanS_BR.adcp_v = get_FACE_Hazen_currents_filter(face.HanS_BR.adcp_v);
    % face.AOML_HW.adcp_u = get_FACE_Hazen_currents_filter(face.AOML_HW.adcp_u);
    % face.AOML_HW.adcp_v = get_FACE_Hazen_currents_filter(face.AOML_HW.adcp_v);
    % face.AOML_BR.adcp_u = get_FACE_Hazen_currents_filter(face.AOML_BR.adcp_u);
    % face.AOML_BR.adcp_v = get_FACE_Hazen_currents_filter(face.AOML_BR.adcp_v);

    disp(['RE-Saving ',matfname]);
    save(matfname,'-struct','face');
  end;

return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INTERNAL FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ts = get_FACE_Hazen_currents_filter(ts)
%function ts = get_FACE_Hazen_currents_filter(ts)
% Remove any NaNs from time series TS (may delete rows from .date,.data, and .prof, if it exists)
  badix = find(isnan(ts.data));
  if ( ~isempty(badix) )
    ts.date(badix) = [];
    ts.data(badix) = [];
    if ( isfield(ts,'prof') )
      ts.prof(badix,:) = [];
    end;
  end;
return;
