1;
%% SCRIPT append_nwps_mat_months
%
% Build up and save a single STRUCT with all Nearshore Wave Prediction System
% (NWPS) Miami region CG1 grid data for Year-Months YRMOS(1) to YRMOS(end)
% (DEFAULT: yrmos = [201608:201612,201701:201707]). Save data in MATFNAME
% (DEFAULT: sprintf('DATAPATH/NWPS/nwps_mat_months_%d-%d.mat',yrmos(1),yrmos(end)))
%
% Last Saved Time-stamp: <Fri 2019-03-22 12:08:53 Eastern Daylight Time gramer>

%for mo=1:12; if ~exist('fld','var'); fld.date=[]; end; x=load(['d:/data/NWPS/NWPS_mfl_nwps_CG1_2017_',sprintf('%02d',mo),'.mat']); fld.lon=x.lon; fld.lat=x.lat; fld.date(end+1:end+numel(x.date)) = x.date; fld.(datestr(datenum(1,mo,1),'mmm')).nwps_sigwavehgt.field = x.nwps_sigwavehgt.field; x=[]; clear x; end;

if ( ~exist('yrmos','var') )
  yrmos = [201608:201612,201701:201707];
end;

if ( ~exist('matfname','var') )
  matfname = fullfile(get_datapath('NWPS'),sprintf('nwps_mat_months_%d-%d.mat',yrmos(1),yrmos(end)));
end;

if ( exist(matfname,'file') )
  disp(['LOAD ',matfname]);
  fld = load(matfname);

else
  disp(['EXTRACT data for ',matfname]);
  
  %             date: [1×279 double]
  %              lat: [201×1 double]
  %              lon: [261×1 double]
  %       model_name: 'mfl_nwps_CG1'
  % nwps_VAR10_0_193: [1×1 struct]
  %    nwps_curr_dir: [1×1 struct]
  %  nwps_curr_speed: [1×1 struct]
  % nwps_primwavedir: [1×1 struct]
  % nwps_primwaveper: [1×1 struct]
  % nwps_sealeveldev: [1×1 struct]
  %  nwps_sigwavehgt: [1×1 struct]
  %    nwps_swellhgt: [1×1 struct]
  % nwps_water_depth: [1×1 struct]
  %    nwps_wind_dir: [1×1 struct]
  %  nwps_wind_speed: [1×1 struct]
  
  fld=[]; clear fld;
  fld.date = []; 
  % ALL 11 FIELDS TOGETHER are about 1,288,033,000 per month or 15 Gb / year!
  fld.nwps_VAR10_0_193.field = [];
  fld.nwps_curr_dir.field = [];
  fld.nwps_curr_speed.field = [];
  fld.nwps_primwavedir.field = [];
  fld.nwps_primwaveper.field = [];
  fld.nwps_sealeveldev.field = [];
  fld.nwps_sigwavehgt.field = [];
  fld.nwps_swellhgt.field = [];
  fld.nwps_water_depth.field = [];
  fld.nwps_wind_dir.field = [];
  fld.nwps_wind_speed.field = [];
  
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
    fld.date(dtixen) = x.date;
    fld.nwps_VAR10_0_193.field(dtixen,latixen,lonixen) = x.nwps_VAR10_0_193.field;
    fld.nwps_curr_dir.field(dtixen,latixen,lonixen) = x.nwps_curr_dir.field;
    fld.nwps_curr_speed.field(dtixen,latixen,lonixen) = x.nwps_curr_speed.field;
    fld.nwps_primwavedir.field(dtixen,latixen,lonixen) = x.nwps_primwavedir.field;
    fld.nwps_primwaveper.field(dtixen,latixen,lonixen) = x.nwps_primwaveper.field;
    fld.nwps_sealeveldev.field(dtixen,latixen,lonixen) = x.nwps_sealeveldev.field;
    fld.nwps_sigwavehgt.field(dtixen,latixen,lonixen) = x.nwps_sigwavehgt.field;
    fld.nwps_swellhgt.field(dtixen,latixen,lonixen) = x.nwps_swellhgt.field;
    fld.nwps_water_depth.field(dtixen,latixen,lonixen) = x.nwps_water_depth.field;
    fld.nwps_wind_dir.field(dtixen,latixen,lonixen) = x.nwps_wind_dir.field;
    fld.nwps_wind_speed.field(dtixen,latixen,lonixen) = x.nwps_wind_speed.field;
    x=[]; clear x;
  end;
  
  disp(['SAVE ',matfname]);
  save(matfname,'-struct','fld');
end;

%flds=fieldnames(fld); for fldix=1:numel(flds); fmg; contourf(squeeze(nanmean(fld.(flds{fldix}).nwps_sigwavehgt.field))); daspect([1,cosd(25),1]); titlename(flds{fldix}); pause; end;
