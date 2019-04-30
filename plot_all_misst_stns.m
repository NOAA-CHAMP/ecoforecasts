1;

more off;

if ( ~exist('stns','var') )
  stns = load_all_misst('misst.cfg','world');
  % stns = load_all_misst('misst-manm1.cfg','world');
  % stns = load_all_misst('misst-br1f1.cfg','world');
  % stns = load_all_misst('misst-ewbf1.cfg','world');
  % stns = load_all_misst('misst-mlrf1.cfg','world');
end;
if ( ~exist('fname','var') || ~exist('gsst','var') )
  % fname = 'mw_ir.fusion.2006.001.v01';
  % fname = 'mw_ir.fusion.2009.360.v02';
  % fname = 'mw_ir.fusion.2009.361.v02';
  % fname = 'mw_ir.fusion.2009.362.v02';
  fname = 'mw_ir.fusion.2009.363.rt';
  % fname = 'mw_ir.fusion.2009.364.v02';
  % fname = 'mw_ir.fusion.2009.365.v02';
  gsst = []; clear gsst;
  % gsst = read_misst(fullfile('data','misst',fname),4096,2048);
  [glon,glat,gsst,dx] = read_misst_region('world',fullfile('data','misst',fname));
end;

for ix=length(stns):-1:1;
  stn=stns(ix);
  stnm=stn.station_name;
  lonix=stn.misst_lonix;
  latix=stn.misst_latix;
  g2lonix=stn.g2_misst_lonix;
  g2latix=stn.g2_misst_latix;
  radx=8;

  poslon = lonix - (lonix(1) - radx) + 1;
  poslat = latix - (latix(1) - radx) + 1;
  g2poslon = g2lonix - (lonix(1) - radx) + 1;
  g2poslat = g2latix - (latix(1) - radx) + 1;

  % if (any(isnan(stn.misst_sst.data)))
  if (any(isnan(stn.g2_misst_sst.data)))
    disp({ix,stnm,stn.lat,stn.lon,latix(1)-1,lonix(1)-1});
    hs = [];
    figure;
    maxigraph;
    t = gsst(latix(1)-radx:latix(1)+radx,lonix(1)-radx:lonix(1)+radx);
    [ig,h] = contourf(t,[min(t(:)),16:0.5:32]);
    hs(1) = h(1);
    hold on;
    h = plot(poslon,poslat,'kd');
    hs(2) = h(1);
    % h = plot(poslon-1,poslat-1,'g*');
    h = plot(g2poslon,g2poslat,'g*');
    hs(3) = h(1);
    daspect([1 1 1]);
    caxis([16,32]);
    % legend(hs, 'MISST','Site','G2 site', 'Location','BestOutside');
    colorbar;
    legend(hs, 'MISST','Site','G2 site', 'Location','NorthWest');
    titlename([strrep(fname,'_','\_') ' ' upper(stnm) ' (' num2str(ix) ')']);
  end;

end;

clear ix stn stnm lonix latix radx t ig h hs
% clear poslon poslat g2poslon g2poslat 
more on;
