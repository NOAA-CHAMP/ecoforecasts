1;
% Extract and plot figures for mooring deployments associated with Eddy /
% Internal Wave experiment on the NOAA R/V Nancy Foster NF0914 in 2009.  See
% notes in load_NF09_mooring_hobos.m regarding this experiement and results.
%
% Last Saved Time-stamp: <Thu 2018-03-01 21:25:57 Eastern Standard Time gramer>

if ( ~exist('doFigs','var') )
  doFigs = true;
end;
if ( ~exist('doPrint','var') )
  doPrint = false;
end;
if ( ~exist('fwyf1','var') )
  fwyf1 = get_station_from_station_name('fwyf1'); fwyf1 = load_all_ndbc_data(fwyf1);
end;

if ( ~exist('nf09moor','var') || ~isfield(nf09moor,'brwd20') )
  nf09moor=[]; clear nf09moor
  matfname = fullfile(get_ecoforecasts_path('data'),'NF09_moorings.mat');

  if ( exist(matfname,'file') )
    disp(['Loading ',matfname]);
    nf09moor = load(matfname);

  else
    nf09moor = read_NF09_moorings;

    % This field was incorrectly converted by READ_HOBO_CSV.m
    % Convert dbar -> psi
    mbar = nf09moor.boca40.boca40_prs.seapres.data ./ 0.68947573;
    nf09moor.baro.date = nf09moor.boca40.boca40_prs.seapres.date;
    nf09moor.baro.data = interp1(fwyf1.ndbc_barom.date,fwyf1.ndbc_barom.data,nf09moor.baro.date);

    nf09moor.boca40.boca40_prs.seapressure.date = nf09moor.boca40.boca40_prs.seapres.date;
    nf09moor.boca40.boca40_prs.seapressure.data = (mbar - nf09moor.baro.data)./100;
    nf09moor.boca40.boca40_prs.seadepth.date = nf09moor.boca40.boca40_prs.seapres.date;
    nf09moor.boca40.boca40_prs.seadepth.data = sw_dpth(nf09moor.boca40.boca40_prs.seapressure.data,+25);

    mbar = []; clear mbar

    disp(['Saving ',matfname]);
    save(matfname,'-struct','nf09moor');
  end;
end;

if ( doFigs )
  [fh,lhs,axs] = multiplot_station(nf09moor.brwd20.brwd20_lgt,{'seatemp','par'},...
                                   'Broward 20m & 40m Mooring Data',...
                                   [],{'T [^oC]','Light, Pressure [dbar]'},...
                                   [datenum(2009,10,27),datenum(2009,11,9)],...
                                   {[18,30],[0,50]},[],'k-');
  hold(axs(1),'on');
  lhs=plot_ts(axs(1),...
              nf09moor.brwd20.brwd20_top.seatemp,nf09moor.brwd20.brwd20_btm.seatemp,...
              nf09moor.brwd40.brwd40_mid.seatemp,nf09moor.brwd40.brwd40_btm.seatemp,...
              nf09moor.brwd100.brwd100_upp.seatemp,nf09moor.brwd100.brwd100_low.seatemp,...
              nf09moor.brwd100.brwd100_btm.seatemp,...
              fwyf1.ndbc_air_t,'k:','LineWidth',1.5);
  set(lhs,'Marker','none');
  legend('20 Sfc','20 Top','20 Btm','40 Mid','40 Btm','100 Upp','100 Lwr','100 Btm','Air', ...
         'Orientation','horizontal','Location','Best');
  hold(axs(2),'on');
  plot_ts(axs(2),nf09moor.brwd40.brwd40_prs.seapres,'b-',nf09moor.brwd100.brwd100_prs.seapres,'g-');
  legend('20 Lgt','40 Pres','100 Pres', 'Orientation','horizontal','Location','Best');
  if doPrint; print('-dtiff','NF09_broward_20m_40m_100m_snail.tiff');
  else; disp('Figure not printed this time...'); end;


  [fh,lhs,axs] = multiplot_station(nf09moor.brwd20.brwd20_lgt,{'seatemp','par'},...
                                   'Boca 20m & 40m Mooring Data',...
                                   [],{'T [^oC]','Light, Pressure [dbar]'},...
                                   [datenum(2009,10,27),datenum(2009,11,9)],...
                                   {[18,30],[0,50]},[],'k-');
  hold(axs(1),'on');
  lhs=plot_ts(axs(1),...
              nf09moor.boca20.boca20_top.seatemp,nf09moor.boca20.boca20_btm.seatemp,...
              nf09moor.boca40.boca40_mid.seatemp,nf09moor.boca40.boca40_btm.seatemp,...
              fwyf1.ndbc_air_t,'k:','LineWidth',1.5);
  set(lhs,'Marker','none');
  legend('(Brwd Lgt)','20 Top','20 Btm','40 Mid','40 Btm','Air', ...
         'Orientation','horizontal','Location','Best');
  hold(axs(2),'on');
  plot_ts(axs(2),nf09moor.brwd40.brwd40_prs.seapres,'b-',nf09moor.boca40.boca40_prs.seapressure,'m-');
  legend('Brwd 20 Lgt','Brwd 40 Pres','Boca 40 Pres', 'Orientation','horizontal','Location','Best');
  if doPrint; print('-dtiff','NF09_boca_20m_40m_snail.tiff');
  else; disp('Figure not printed this time...'); end;
end;
