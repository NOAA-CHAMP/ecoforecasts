1;

clear ans d w N2 Pv S S2 spd spd2 stnm t tmp T Ts wd z Zs

stnm = 'ne_buoy';
%stnm = 'e_buoy';

if ( ~exist('sfomc','var') || ~isfield(sfomc,stnm) )
  sfomc = []; clear sfomc
  sfomc = get_sfomc_data;
end;

[tmp,spd] = intersect_tses(sfomc.(stnm).seatemp,sfomc.(stnm).adcp_speed);

%[N2,Pv] = gsw_Nsquared(repmat(35,[numel(15:30),4]),repmat([15:30]',[1,4]),[5,15,25,35],26)
t = tmp.date;
T = tmp.prof(:,end:-1:1);
S = repmat(35,size(T));
z = -tmp.depths(end:-1:1);

% Calculate Brunt-Väisälä (buoyancy) frequency squared: N^2
N2.date = tmp.date;
[n2,Pv] = gsw_Nsquared(S',T',z',sfomc.(stnm).lat);
N2.prof = n2';
N2.depths = -Pv(:,1)';

% % Hovmoeller of buoyancy frequency squared
% fmg; surf(N2.date,N2.depths',N2.prof'); colorbar; shading interp; datetick3;
% colormap(hot); axis([tmp.date(1),datenum(1999,11,1), -43,0, 0,1, 0,0.0035]);

% % % Hovmoeller of buoyancy frequency squared as an inverse fraction of Coriolis frequency squared
% % fmg; surf(t,-Pv(:,1),real(log( (sw_f(sfomc.(stnm).lat).^2)./N2 ))); colorbar; shading interp; datetick3;
% % colormap(hot); axis([tmp.date(1),datenum(1999,11,1), -43,0, 0,1, log([1e-3,0.15])]);

spd.depths = -z;
wd = warning('OFF','MATLAB:interp2:NaNstrip');
spd.interp_prof = interp2(spd.date,sfomc.(stnm).adcp_depths,spd.prof',spd.date,spd.depths,'linear')';
warning(wd);

% Calculate vertical velocity shear squared: (du/dz)^2
S2.date = spd.date;
S2.prof = ((diff(spd.interp_prof')./diff(spd.depths)')').^2;
S2.depths = N2.depths;


% Calculate Richardson number => Ri < 0.25 suggesting Kelvin-Helmholtz instability
Ri.date = N2.date;
Ri.prof = N2.prof./S2.prof;
Ri.depths = N2.depths;
Ri.prof(Ri.prof<0) = 0;

fmg; surf(Ri.date,Ri.depths',Ri.prof'); datetick3; shading interp; set(gca,'CLim',[0,0.25]); colorbar;
%fmg; contourf(Ri.date,Ri.depths,Ri.prof',[0:0.05:0.25,100]); datetick3; shading interp; set(gca,'CLim',[0,0.25]); colorbar;

if strcmp(stnm,'e_buoy'); xlim(datenum(2000,[6,7],15)); datetick3; end;
%if strcmp(stnm,'ne_buoy'); xlim(datenum(1999,[7,8],15)); datetick3; end;
if strcmp(stnm,'ne_buoy'); xlim(datenum(1999,[7,7],[18,31])); datetick3; end;
