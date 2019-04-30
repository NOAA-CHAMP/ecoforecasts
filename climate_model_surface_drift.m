1;

if ( ~exist('CLIMATE_RFP_DIR','var') )
  CLIMATE_RFP_DIR = 'c:/Users/gramer/Documents/RFPs';
end;
if ( ~exist('dpath','var') )
  %dpath = get_ecoforecasts_path('data');
  % Ran out of room on my D: drive! Access \\dana\OCD-XFERS (Y:) instead:
  dpath = 'y:/';
end;
if ( ~exist('doFigs','var') )
  doFigs = false;
end;
if ( ~exist('doPrint','var') )
  doPrint = false;
end;
if ( doPrint )
  disp(['Will PRINT Figures to "',CLIMATE_RFP_DIR,'"']);
end;

% From: 
%  https://www.earthsystemgrid.org/project/ccsm.html
% CESM CAM5 BGC with fixed 2005 aerosols, Ocean Post Processed Data, Monthly Averages, UVEL, version 1
%  https://www.earthsystemgrid.org/dataset/ucar.cgd.ccsm4.CESM_CAM5_BGC_YY.ocn.proc.monthly_ave.UVEL.html
% CESM CAM5 BGC with fixed 2005 aerosols, Ocean Post Processed Data, Monthly Averages, UVEL2, version 1
% CESM CAM5 BGC with fixed 2005 aerosols, Ocean Post Processed Data, Monthly Averages, VVEL, version 1

% CESM CAM5 BGC with fixed 2005 aerosols, Atmosphere Post Processed Data, Monthly Averages, U, version 1
% CESM CAM5 BGC with fixed 2005 aerosols, Atmosphere Post Processed Data, Monthly Averages, V, version 1


matfname = fullfile(dpath,'b.e11.BRCP85C5CNBDRD.f09_g16.001.b.pop.h.200601-210012.mat');

if ( ~exist('u','var') )

  if ( exist(matfname,'file') )
    disp(['Loading ',matfname]);
    load(matfname);

  else

    %ang=[]; ang1=[]; ang2=[]; angdif=[]; spd=[]; u=[]; u1=[]; u2=[]; v=[]; v1=[]; v2=[]; clear all; pack
    disp('Extracting CCSM data from netCDF files...');
    %DEBUG:
    keyboard;

    % double time(time=1140);
    %   :long_name = "time";
    %   :units = "days since 0000-01-01 00:00:00";
    %   :bounds = "time_bound";
    %   :calendar = "noleap";
    %
    % >> t = cast(nc{'time'}(:),'double');
    % >> dts = datenum(0,1,1,0,0,0) + t(1:200);
    % >> datestr(dts([1,end]))
    % ans =
    % 02-Oct-2004
    % 28-Apr-2021
    % >> dep(1:4)
    % ans =
    %      5
    %     15
    %     25
    %     35

    nc = mDataset('y:/b.e11.BRCP85C5CNBDRD.f09_g16.001.b.pop.h.UVEL.200601-210012.nc');
    LAT = cast(nc{'ULAT'}(:,:),'double');
    LON = cast(nc{'ULONG'}(:,:),'double');
    dep = cast(nc{'z_t'}(:),'double')./1e2;
    deptop = cast(nc{'z_w'}(:),'double')./1e2;
    t = cast(nc{'time'}(:),'double');
    dts = datenum(0,1,1,0,0,0) + t(1:200);
    %% TOO BIG to extract!
    %u = cast(nc{'UVEL'}(1:275,1:5,:,:),'double');
    u = cast(nc{'UVEL'}(1:200,1:4,:,:),'double');
    close(nc); clear nc ans

    nc = mDataset('y:/b.e11.BRCP85C5CNBDRD.f09_g16.001.b.pop.h.VVEL.200601-210012.nc');
    % vLAT = cast(nc{'ULAT'}(:,:),'double');
    % vLON = cast(nc{'ULONG'}(:,:),'double');
    % vdep = cast(nc{'z_t'}(:),'double')./1e2;
    % vdeptop = cast(nc{'z_w'}(:),'double')./1e2;
    % vt = cast(nc{'time'}(:),'double');
    % vdts = datenum(0,1,1,0,0,0) + vt(1:200);
    %v = cast(nc{'VVEL'}(1:275,1:5,:,:),'double');
    v = cast(nc{'VVEL'}(1:200,1:4,:,:),'double');
    close(nc); clear nc

    disp(['Saving ',matfname]);
    save(matfname);
  end; %if ( exist(matfname,'file') ) else

  %DEBUG DEBUG DEBUG DEBUG:Uraw=[]; Vraw=[]; U=[]; V=[]; clear U Uraw V Vraw ULAT ULON Ut Udts
  %DEBUG DEBUG DEBUG DEBUG

  if ( ~exist('U','var') )
    % double time(time=1140);
    %   :long_name = "time";
    %   :units = "days since 2006-01-01 00:00:00";
    %   :calendar = "noleap";
    %   :bounds = "time_bnds";
    %   :_CoordinateAxisType = "Time";
    % >> datestr(dts(17:20))
    % ans =
    % 30-Jan-2006
    % 01-Mar-2006
    % 01-Apr-2006
    % 02-May-2006
    % >> datestr(Udts(1:3))
    % ans =
    % 01-Feb-2006
    % 01-Mar-2006
    % 01-Apr-2006
    % >> datestr(dts(end))
    % ans =
    % 28-Apr-2021
    % >> datestr(Udts(184))
    % ans =
    % 27-Apr-2021

    nc = mDataset(fullfile(dpath,'b.e11.BRCP85C5CNBDRD.f09_g16.001.b.cam.h0.U.200601-210012.nc'));
    ULAT = cast(nc{'lat'}(:,:),'double');
    ULON = cast(nc{'lon'}(:,:),'double');
    %Ulev = cast(nc{'lev'}(:),'double');
    Ut = cast(nc{'time'}(:),'double');
    Udts = Ut(1:184) + datenum(2006,1,1,0,0,0);
    %DEBUG:  Uraw = cast(nc{'U'}(1,:,:,:),'double');
    % Level "end" (30) is (near) the surface
    Uraw = cast(nc{'U'}(1:184,end,:,:),'double');
    close(nc); clear nc

    nc = mDataset(fullfile(dpath,'b.e11.BRCP85C5CNBDRD.f09_g16.001.b.cam.h0.V.200601-210012.nc'));
    % VLAT = cast(nc{'lat'}(:,:),'double');
    % VLON = cast(nc{'lon'}(:,:),'double');
    % %Vlev = cast(nc{'lev'}(:),'double');
    % Vt = cast(nc{'time'}(:),'double');
    % Vdts = Vt(1:184) + datenum(2006,1,1,0,0,0);
    Vraw = cast(nc{'V'}(1:184,end,:,:),'double');
    close(nc); clear nc

    U = interp3(unique(ULAT),Udts,unique(ULON),Uraw,LAT(:,1),dts,LON(1,:),'linear',nan);
    V = interp3(unique(ULAT),Udts,unique(ULON),Vraw,LAT(:,1),dts,LON(1,:),'linear',nan);

    disp(['RE-Saving ',matfname]);
    save(matfname);
  end; %if ( ~exist('U','var') )

end; %if ( ~exist('u','var') )


if ( ~exist('angdif','var') )
  disp('Deriving quantities');
  u = u./100;
  v = v./100;
  ang = uv_to_dir_curr(u,v);
  spd = uv_to_spd(u,v);

  u1=squeeze(u(:,1,:,:)); u2=squeeze(u(:,2,:,:));
  v1=squeeze(v(:,1,:,:)); v2=squeeze(v(:,2,:,:));

  angdif = uv_to_dir_curr(u2-u1,v2-v1);
  if ( doFigs )
    % Only compare vector angle differences where at least one of the two
    % vector components (u and/or v) has a magnitude of at least 1 mm/s!
    nominalix = find(abs(u1)>=0.001 | abs(v1)>=0.001);
    fmg; hist(angdif(nominalix),100); xlim([0,360]);
    if doPrint; print('-dpng',fullfile(CLIMATE_RFP_DIR,'Level_1_vs_Level_2_current_angle_difference_hist.png')); end;
  end;
end;


if ( doFigs )
  %tix=41; 
  [err,tix] = min(abs(dts - datenum(2008,1,30)));
  ixen=295:320; jxen=245:295;

  fmg;
  %plot_lores_coastline(LON(jxen,ixen)-360,LAT(jxen,ixen)); %MAP toolkit gives goofy figures
  contourf(LON(jxen,ixen)-360,LAT(jxen,ixen),squeeze(spd(tix,1,jxen,ixen)),[0:0.02:1.0]);
  colorbar('AxisLocation','in', 'Location','EastOutside');
  quiver(LON(jxen,ixen)-360,LAT(jxen,ixen),squeeze(u(tix,1,jxen,ixen)),squeeze(v(tix,1,jxen,ixen)));
  quiver(LON(jxen,ixen)-360,LAT(jxen,ixen),squeeze(u(tix,2,jxen,ixen)),squeeze(v(tix,2,jxen,ixen)));
  quiver(LON(jxen,ixen)-360,LAT(jxen,ixen),squeeze(u(tix,3,jxen,ixen)),squeeze(v(tix,3,jxen,ixen)));

  [err,Utix] = min(abs(Udts-dts(tix)));
  quiver(LON(jxen,ixen)-360,LAT(jxen,ixen),squeeze(U(Utix,jxen,ixen)),squeeze(V(Utix,jxen,ixen)),'k--');

  legend(['@',num2str(dep(1)),' m'],['@ ',num2str(dep(1)),' m'],...
         ['@',num2str(dep(2)),' m'],['@',num2str(dep(3)),' m'],'Wind',...
         'Location','Best');
  titlename(datestr(dts(tix)));
  if doPrint; print('-dpng',fullfile(CLIMATE_RFP_DIR,'currents_map_3_levels.png')); end;

  disp('Hit enter to zoom to Sargasso sub-region'); pause;
  axis([-64,-50,30,37]); daspect([1,cosd(30),1]);
  if doPrint; print('-dpng',fullfile(CLIMATE_RFP_DIR,'currents_map_3_levels_Sargasso.png')); end;

  disp('Hit enter to zoom to Antilles sub-region'); pause;
  axis([-68,-60,22,27]); daspect([1,cosd(25),1]);
  if doPrint; print('-dpng',fullfile(CLIMATE_RFP_DIR,'currents_map_3_levels_Antilles.png')); end;
end;
