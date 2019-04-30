1;
%% SCRIPT catenate_nwps_mat_months_vars.m
%
% Build up and save a single STRUCT per variable, with all data for *that
% variable* from the Nearshore Wave Prediction System (NWPS) Miami grid CG1
% data, for YRMOS (DEFAULT: yrmos = [201608:201612,201701:201707]). Save to
% sprintf('DATAPATH/NWPS/nwps_mat_months_%s_%d-%d.mat',fldnm,yrmos(1),yrmos(end))).
%
% Last Saved Time-stamp: <Fri 2019-03-22 12:19:23 Eastern Daylight Time gramer>

more off;

if ( ~exist('yrmos','var') )
  yrmos = [201608:201612,201701:201707];
end;

allflds = {
    'nwps_wind_dir',    	'Wind_direction_from_which_blowing_surface' ; ...
    'nwps_wind_speed',  	'Wind_speed_surface' ; ...
    'nwps_sigwavehgt',     	'Significant_height_of_combined_wind_waves_and_swell_surface' ; ...
    'nwps_swellhgt',    	'Significant_height_of_swell_waves_surface' ; ...
    'nwps_primwavedir', 	'Primary_wave_direction_surface' ; ...
    'nwps_primwaveper', 	'Primary_wave_mean_period_surface' ; ...
    'nwps_VAR10_0_193', 	'VAR10-0-193_FROM_7-0-0_surface' ; ...
    'nwps_curr_dir',    	'Current_direction_surface' ; ...
    'nwps_curr_speed',  	'Current_speed_surface' ; ...
    'nwps_sealeveldev', 	'Deviation_of_sea_level_from_mean_surface' ; ...
    'nwps_water_depth', 	'Water_depth_surface' ; ...
          };

for fldix = 1:size(allflds,1)
  fldnm = allflds{fldix,1};
  disp(fldnm);

  matfname = fullfile(get_datapath('NWPS'),sprintf('%s_%d_%d.mat',fldnm,yrmos(1),yrmos(end)));
  if ( exist(matfname,'file') )
    disp(['SKIPPING ',matfname]);

  else
    disp(['EXTRACT data for ',matfname]);
    fld=[]; clear fld;
    fld.date = []; 
    % ALL 11 FIELDS TOGETHER are about 1,288,033,000 per month or 15 Gb / year!
    fld.(fldnm) = [];
    
    for yrmo = yrmos(:)';
      disp(yrmo); 
      yr = floor(yrmo/1e2);
      mo = yrmo - (yr*1e2);
      x=load(fullfile(get_datapath('NWPS'),['NWPS_mfl_nwps_CG1_',sprintf('%04d_%02d',yr,mo),'.mat']));
      
      if ( ~isfield(fld,'lon') )
        fld.model_name = x.model_name;
        fld.lon = x.lon;
        fld.lat = x.lat;
        lonixen = 1:numel(fld.lon);
        latixen = 1:numel(fld.lat);
      end;
      dtixen = numel(fld.date)+[1:numel(x.date)];
      fld.(fldnm).date(dtixen) = x.date;
      fld.(fldnm).field(dtixen,latixen,lonixen) = x.(fldnm).field;
      x=[]; clear x;
    end; %for yrmo = yrmos(:)';
    disp(['SAVE ',matfname]);
    save(matfname,'-struct','fld');

  end; %if ( exist(matfname,'file') ) else

end; %for fldix = 1:size(allflds,1)

more on;
