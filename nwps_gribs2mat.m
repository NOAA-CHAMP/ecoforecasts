1;
%% SCRIPT nwps_gribs2mat.m
%
% Read all the .grib2 Nearshore Wave Prediction System (NWPS) files for a
% given date range and produce a .MAT file for each year-month of those data.
%
% Last Saved Time-stamp: <Fri 2019-03-01 23:29:19 Eastern Standard Time gramer>


set_more off;

if ( ~exist('datapath','var') )
  %datapath = get_ecoforecasts_path('data');
  % New centrally located data directory structure
  datapath = get_ecoforecasts_path('../data');
end;
if ( ~exist('nwpspath','var') )
  nwpspath = fullfile(datapath,'NWPS');
end;

if ( ~exist('mdl','var') )
  mdl = 'mfl_nwps_CG1';
  %mdl = 'key_nwps_CG2';
end;

if ( ~exist('allflds','var') )
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
end;

switch (mdl),
 case 'mfl_nwps_CG1',		nlon=261; nlat=201;
 case 'key_nwps_CG2',		nlon=315; nlat=123;
 otherwise,          		error('Unknown MODEL %s',mdl);
end;

if ( ~exist('yrs','var') )
  %%yrs = 2016:2019;
  %yrs = 2016:2016;
  yrs = 2017:2017;
  %%yrs = 2018:2018;
end;

if ( ~exist('mos','var') )
  %mos = 1:12
  %mos = 7:10
  mos = 3:3;
end;

for yr=yrs(:)'
  for mo=mos(:)'
    matfname = fullfile(nwpspath,sprintf('NWPS_%s_%04d_%02d.mat',mdl,yr,mo));
    if ( exist(matfname,'file') )
      disp(['SKIPPING ',matfname]);
      %%%% EARLY CONTINUE
      continue;
    end;

    disp(sprintf('Extracting NWPS for %s_%04d_%02d',mdl,yr,mo));
    alldts = datenum(yr,mo,1,6,0,0):0.5:datenum(yr,mo+1,1,0,0,0);
    if ( alldts(1) < datenum(2018,1,1) )
      dhr = 3;
      ntms = 35;
    else
      dhr = 1;
      ntms = 145;
    end;
    
    field=[]; clear field
    field.model_name = mdl;
    field.date = alldts(1):(dhr/24):(alldts(end)+((ntms-1)*(dhr/24)));
    ndts = numel(field.date);
    
    for dt=alldts
      dy = get_dom(dt);
      hr = get_hour(dt);
      %nc = mDataset('d:/ecoforecasts/data/NWPS/key_nwps_CG2_20160828_0600.grib2');
      ncfname = fullfile(nwpspath,sprintf('%s_%04d%02d%02d_%02d00.grib2',mdl,yr,mo,dy,hr));
      if ( ~exist(ncfname,'file') )
        warning('Missing %s',ncfname);
        % EARLY LOOP CONTINUE
        continue;
      end; %if ( exist(matfname,'file') )
      
      %DEBUG:      disp(ncfname);
      nc = mDataset(ncfname);
      try,
        if ( ~isfield(field,'lon') )
          field.lon = cast(nc{'lon'}(:,:),'double');
          field.lat = cast(nc{'lat'}(:,:),'double');
          nlon = numel(field.lon);
          nlat = numel(field.lat);
        end;
        
        dts = dt + ( cast(nc{'time'}(:,:),'double')/24 );
        %dtix = find(ismember(field.date,dts));
        [fdtix,dtix] = intersect_dates(field.date,dts);
        if ( numel(fdtix) ~= numel(dtix) )
          error('DATE/DATA SIZE MISMATCH: %s!',ncfname);
        end;
        
        for fldix=1:size(allflds,1)
          fld = allflds{fldix,1};
          try,
            dat = cast(nc{allflds{fldix,2}}(:,:,:),'double');
          catch varME,
            close(nc); clear nc
            warning('MISSING %s from %s',allflds{fldix,2},ncfname);
            %%%% EARLY CONTINUE
            continue;
          end;
          %DEBUG:      disp(fld);
          if ( ~isfield(field,fld) )
            field.(fld).field = repmat(nan,[ndts,nlat,nlon]);
          end;
          field.(fld).field(fdtix,:,:) = dat(dtix,:,:);
        end; %for fldix=1:size(allflds,1)
      catch ME,
        close(nc); clear nc
        warning('ERROR IN %s',ncfname);
      end; %try, catch ME,
      try, close(nc); clear nc; catch ME, end;
    end; %for dt=alldts
    
    if ( ~isfield(field,'nwps_sigwavehgt') || ...
         all(isnan(field.nwps_sigwavehgt.field(:))) )
      warning('ALL Significant Wave Heights WERE NaN! "%s"',matfname);
    else
      disp(['Save ',matfname]);
      save(matfname,'-struct','field');
    end;
  end; %for mo=1:12
end; %for yr=yrs(:)'


if 0;
  fmg;
  contourf(field.lon,field.lat,squeeze(nanmean(field.nwps_sigwavehgt.field)),[0:0.1:2.0]);
  set(gca,'clim',[0,3]); colorbar('Location','East');
  daspect([1,cosd(field.lat(1)),1]);
  
  fmg;
  contourf(field.lon,field.lat,squeeze(nanmax(field.nwps_sigwavehgt.field)),[0:0.1:2.0]);
  set(gca,'clim',[0,3]); colorbar('Location','East');
  daspect([1,cosd(field.lat(1)),1]);
end;

set_more;
