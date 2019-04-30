function res = extract_xaymara_sites
%function res = extract_xaymara_sites
%
% Extract locations, seafloor depths and slopes for Xaymara Serrano's
% PhD sampling sites: locations and other site descriptors provided by
% Xaymara, seafloor data courtesy of NOAA's National Geodetic Data Ctr.
% Depth and slope are calculated from bathymetry in several ways, and stored
% both in individual Station STRUCTs (RES.ms and RES.ps) and vectors, e.g.,
% for M. cav. sites depth is in RES.mdeps, RES.mndeps, and RES.vh_mhdeps.
%
% Loads from '$ECOFORECASTS/data/xaymara_sites.mat' if present.
%
% Last Saved Time-stamp: <Sun 2017-02-05 17:28:23 Eastern Standard Time gramer>

  doPlot = false;

  datpath = get_ecoforecasts_path('data');
  matfname = fullfile(datpath,'xaymara_sites.mat');

  if ( exist(matfname,'file') )
    disp(['Loading ',matfname]);
    load(matfname);

  else
  tic,

    disp('Extracting Xaymara Keys sites');
    res.lons=[];
    res.lats=[];

    res.rawfname = fullfile(datpath,'Xaymara sites for Lew_LJG.xlsx');
    rawin = importdata(res.rawfname);

    res.nms = size(rawin.data.FLKeysSitesMcav,1);
    res.nps = size(rawin.data.FLKeysSitesPast,1);

    for stix = 1:res.nms
      res.ms(stix).station_name = rawin.textdata.FLKeysSitesMcav{stix+1,4};
      res.ms(stix).station_desc = rawin.textdata.FLKeysSitesMcav{stix+1,3};
      res.ms(stix).region = rawin.textdata.FLKeysSitesMcav{stix+1,1};
      res.ms(stix).subregion = rawin.textdata.FLKeysSitesMcav{stix+1,2};
      res.ms(stix).depth_category = rawin.textdata.FLKeysSitesMcav{stix+1,5};
      res.ms(stix).lon = rawin.data.FLKeysSitesMcav(stix,3);
      res.ms(stix).lat = rawin.data.FLKeysSitesMcav(stix,2);
      res.ms(stix).depth = -rawin.data.FLKeysSitesMcav(stix,1);
      res.ms(stix).N = rawin.data.FLKeysSitesMcav(stix,4);
      res.ms(stix).Ng = rawin.data.FLKeysSitesMcav(stix,5);
      res.mlons(stix,1) = res.ms(stix).lon;
      res.mlats(stix,1) = res.ms(stix).lat;
      res.lons(end+1,1) = res.ms(stix).lon;
      res.lats(end+1,1) = res.ms(stix).lat;
      res.mN(stix,1) = res.ms(stix).N;
      res.mNg(stix,1) = res.ms(stix).Ng;
    end;

    for stix = 1:res.nps
      res.ps(stix).station_name = rawin.textdata.FLKeysSitesPast{stix+1,4};
      res.ps(stix).station_desc = rawin.textdata.FLKeysSitesPast{stix+1,3};
      res.ps(stix).region = rawin.textdata.FLKeysSitesPast{stix+1,1};
      res.ps(stix).subregion = rawin.textdata.FLKeysSitesPast{stix+1,2};
      res.ps(stix).depth_category = rawin.textdata.FLKeysSitesPast{stix+1,5};
      res.ps(stix).lon = rawin.data.FLKeysSitesPast(stix,3);
      res.ps(stix).lat = rawin.data.FLKeysSitesPast(stix,2);
      res.ps(stix).depth = -rawin.data.FLKeysSitesPast(stix,1);
      res.ps(stix).N = rawin.data.FLKeysSitesPast(stix,4);
      res.ps(stix).Ng = rawin.data.FLKeysSitesPast(stix,5);
      res.plons(stix,1) = res.ps(stix).lon;
      res.plats(stix,1) = res.ps(stix).lat;
      res.lons(end+1,1) = res.ps(stix).lon;
      res.lats(end+1,1) = res.ps(stix).lat;
      res.pN(stix,1) = res.ps(stix).N;
      res.pNg(stix,1) = res.ps(stix).Ng;
    end;

    rawin=[]; clear rawin stix

    % disp(['NOT Saving ',matfname]);
    % %save(matfname,'nmstns','npstns','mstns','mlons','mlats','pstns','plons','plats','lons','lats','-v7.3');


    % Find geographic center point and bounding box for range of all sites
    res.lon = mean([min(res.lons),max(res.lons)]);
    res.lat = mean([min(res.lats),max(res.lats)]);
    [ig,minix] = min(res.lons);
    [ig,maxix] = max(res.lons);
    res.rady = ceil(max([distance_wgs84(res.lats(minix),res.lons(minix),res.lat,res.lons(minix)),...
                        distance_wgs84(res.lats(maxix),res.lons(maxix),res.lat,res.lons(maxix))]))*1e3;
    [ig,minix] = min(res.lats);
    [ig,maxix] = max(res.lats);
    res.radx = ceil(max([distance_wgs84(res.lats(minix),res.lons(minix),res.lats(minix),res.lon),...
                        distance_wgs84(res.lats(maxix),res.lons(maxix),res.lats(maxix),res.lon)]))*1e3;

    clear ig minix maxix

    % Map wide-view bathymetry
    if ( doPlot )
      res = plot_hires_bathymetry(res,-[0:5:50,100:100:800],[res.radx*1.05,res.rady*1.10],true,[],[],[],false);
      axis(axis);
      plot(res.lons,res.lats,'ws','MarkerFaceColor','w');
      plot(res.mlons,res.mlats,'k.','MarkerSize',0.5);
      plot(res.plons,res.plats,'r.','MarkerSize',0.5);
    else
      res = read_hires_bathymetry(res,[res.radx*1.05,res.rady*1.10],[],false);
    end;

    [LON,LAT] = meshgrid(res.ngdc_hires_bathy.lon,res.ngdc_hires_bathy.lat);
    [res.aspect,res.beta_deg,res.beta_y,res.beta_x] = gradientm(LAT,LON,res.ngdc_hires_bathy.field,wgs84Ellipsoid);
    LON=[]; LAT=[]; clear LON LAT
    res.beta = sqrt((res.beta_y.^2)+(res.beta_x.^2));

    for stix = 1:res.nms
      x.lon = res.ms(stix).lon;
      x.lat = res.ms(stix).lat;
      %x = read_hires_bathymetry(x,[4e3,4e3],[],true);
      x = read_hires_bathymetry(x,[4e3,4e3],[],false);
      res.ms(stix).ngdc_hires_bathy = x.ngdc_hires_bathy;
      x=[]; clear x

      res.mdeps(stix) = res.ms(stix).depth;
      res.mndeps(stix) = interp2(res.ngdc_hires_bathy.lon,res.ngdc_hires_bathy.lat,...
                                 res.ngdc_hires_bathy.field,...
                                 res.ms(stix).lon,res.ms(stix).lat);
      res.ms(stix).ngdc_depth = res.mndeps(stix);

      res.mbetas(stix) = interp2(res.ngdc_hires_bathy.lon,res.ngdc_hires_bathy.lat,...
                                 res.beta,res.ms(stix).lon,res.ms(stix).lat);
      res.ms(stix).ngdc_beta = res.mbetas(stix);
    end;

    for stix = 1:res.nps
      x.lon = res.ps(stix).lon;
      x.lat = res.ps(stix).lat;
      %x = read_hires_bathymetry(x,[4e3,4e3],[],true);
      x = read_hires_bathymetry(x,[4e3,4e3],[],false);
      res.ps(stix).ngdc_hires_bathy = x.ngdc_hires_bathy;
      x=[]; clear x

      res.pdeps(stix) = res.ps(stix).depth;
      res.pndeps(stix) = interp2(res.ngdc_hires_bathy.lon,res.ngdc_hires_bathy.lat,...
                                 res.ngdc_hires_bathy.field,...
                                 res.ps(stix).lon,res.ps(stix).lat);
      res.ps(stix).ngdc_depth = res.pndeps(stix);

      res.pbetas(stix) = interp2(res.ngdc_hires_bathy.lon,res.ngdc_hires_bathy.lat,...
                                 res.beta,res.ps(stix).lon,res.ps(stix).lat);
      res.ps(stix).ngdc_beta = res.pbetas(stix);
    end;


    % Finally, recalculate depths and slopes using spatial averaging of HIGHEST
    % resolution bathymetry available (10, 20, or 30 m gridpoints)
    for stix = 1:res.nms
      %DEBUG:      toc, disp(['RecalcSlope M: ',res.ms(stix).station_name]); tic,
      x.lon = res.ms(stix).lon;
      x.lat = res.ms(stix).lat;
      x = read_hires_bathymetry(x,[1e3,1e3],[],true);
      % DEBUG:      disp({res.ms(stix).station_name,min(diff(x.ngdc_hires_bathy.lat))*111e3});
      [N,smoothopt] = extract_xaymara_sites_betasmooth(x);
      % DEBUG:      disp({res.ps(stix).station_name,N,smoothopt{:}});
      x.ngdc_hires_bathy = calc_ngdc_slope(x.ngdc_hires_bathy,N);
      res.vh_mndeps(stix) = interp_field(x.ngdc_hires_bathy.lat,x.ngdc_hires_bathy.lon,...
                                         x.ngdc_hires_bathy.field,x.lat,x.lon,'linear');
      res.vh_mbetas(stix) = interp_field(x.ngdc_hires_bathy.lat,x.ngdc_hires_bathy.lon,...
                                         x.ngdc_hires_bathy.beta,x.lat,x.lon,'linear');
      res.vhs_mndeps(stix) = interp_field(x.ngdc_hires_bathy.lat,x.ngdc_hires_bathy.lon,...
                                         x.ngdc_hires_bathy.field,x.lat,x.lon,smoothopt);
      res.vhs_mbetas(stix) = interp_field(x.ngdc_hires_bathy.lat,x.ngdc_hires_bathy.lon,...
                                         x.ngdc_hires_bathy.beta,x.lat,x.lon,smoothopt);
      res.ms(stix).very_hires_bathy = x.ngdc_hires_bathy;
      res.ms(stix).very_hires_depth = res.vh_mndeps(stix);
      res.ms(stix).very_hires_beta = res.vh_mbetas(stix);
      res.ms(stix).very_hires_smooth_depth = res.vhs_mndeps(stix);
      res.ms(stix).very_hires_smooth_beta = res.vhs_mbetas(stix);
      x=[]; clear x
    end;

    for stix = 1:res.nps
      %DEBUG:      toc, disp(['RecalcSlope P: ',res.ps(stix).station_name]); tic,
      x.lon = res.ps(stix).lon;
      x.lat = res.ps(stix).lat;
      x = read_hires_bathymetry(x,[1e3,1e3],[],true);
      % DEBUG:      disp({res.ps(stix).station_name,min(diff(x.ngdc_hires_bathy.lat))*111e3});
      [N,smoothopt] = extract_xaymara_sites_betasmooth(x);
      % DEBUG:      disp({res.ps(stix).station_name,N,smoothopt{:}});
      x.ngdc_hires_bathy = calc_ngdc_slope(x.ngdc_hires_bathy,N);
      res.vh_pndeps(stix) = interp_field(x.ngdc_hires_bathy.lat,x.ngdc_hires_bathy.lon,...
                                         x.ngdc_hires_bathy.field,x.lat,x.lon,'linear');
      res.vh_pbetas(stix) = interp_field(x.ngdc_hires_bathy.lat,x.ngdc_hires_bathy.lon,...
                                         x.ngdc_hires_bathy.beta,x.lat,x.lon,'linear');
      res.vhs_pndeps(stix) = interp_field(x.ngdc_hires_bathy.lat,x.ngdc_hires_bathy.lon,...
                                         x.ngdc_hires_bathy.field,x.lat,x.lon,smoothopt);
      res.vhs_pbetas(stix) = interp_field(x.ngdc_hires_bathy.lat,x.ngdc_hires_bathy.lon,...
                                         x.ngdc_hires_bathy.beta,x.lat,x.lon,smoothopt);
      res.ps(stix).very_hires_bathy = x.ngdc_hires_bathy;
      res.ps(stix).very_hires_depth = res.vh_pndeps(stix);
      res.ps(stix).very_hires_beta = res.vh_pbetas(stix);
      res.ps(stix).very_hires_smooth_depth = res.vhs_pndeps(stix);
      res.ps(stix).very_hires_smooth_beta = res.vhs_pbetas(stix);
      x=[]; clear x
    end;


    clear stix

    disp(['Saving ',matfname]);
    save(matfname,'res','-v7.3');
  toc,
  end; %if ( exist(matfname,'file') ) else (tic, toc,)

  % How often does our bathymetry reflect Xaymara's depth estimates?
  if 0;
    clear mdd mddix
    for stix = 1:res.nms
      x = res.ms(stix);
      xb = x.very_hires_bathy;
      [LON,LAT] = meshgrid(xb.lon,xb.lat);
      % All bathymetry points within ~111 m radius (swimming distance)
      midix = find(abs(x.lon - LON(:)) < 0.001 & abs(x.lat - LAT(:)) < 0.001);
      % Find the bathymetry point with the smallest difference vs. Xaymara
      [mdd(stix),mddix(stix)] = min(abs(xb.field(midix)-x.depth));
      LON=[]; LAT=[]; x=[]; xb=[]; clear LON LAT midix x xb
    end;
    mdd,
    mean(mdd),
    % % Again, these are differences between her estimate and those from the
    % % bathymetry field that match it most CLOSELY within a 222 m square!
    % >> mdd,
    %  mdd =
    %   Columns 1 through 7
    %          0    1.1300    0.0300    0.0200    0.0100    1.5300         0
    %   Columns 8 through 14
    %     0.1200    0.0500    1.0800    0.0800    8.6300   12.4000    7.4900
    %   Columns 15 through 21
    %     0.0172    0.0173    0.8469    0.6113    1.1370    0.0012         0
    %   Columns 22 through 23
    %     3.2000   17.0000
    % >> mean(mdd)
    % ans =
    %     2.4087
  end;

return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% INTERNAL FUNCTION  %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [N,smoothopt] = extract_xaymara_sites_betasmooth(x)
  latres = min(diff(x.ngdc_hires_bathy.lat))*111e3;
  switch ( roundn(latres,1) ),
   case 10,
    N = 7;
    smoothopt = {@nanmean,12,12,18}; % Minimum 10%
   case 20,
    N = 5;
    smoothopt = {@nanmean,6,6,5}; % Minimum 14%
   case 30,
    N = 3;
    smoothopt = {@nanmean,4,4,3}; % Minimum 19%
   case 60,
    N = 3;
    smoothopt = {@nanmean,3,3,2}; % Minimum 22%
   case 90,
    N = 2;
    smoothopt = 'linear';
   otherwise,
    error('Unknown bathymetry resolution %f [m]',latres);
  end;
return;
