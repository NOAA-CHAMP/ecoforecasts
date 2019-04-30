1;
%%%% SCRIPT to find and display cold-front (or storm) events together with
%%%% corresponding heat budget (sea-surface heat flux) time series


%% LOAD DATA

%stnm = 'lkwf1';
%stnm = 'fwyf1';
stnm = 'mlrf1';
%stnm = 'lonf1';
%stnm = 'smkf1';
%stnm = 'sanf1';
%stnm = 'dryf1';

if ( ~exist('stn','var') )
  mfname = fullfile(get_ecoforecasts_path('data'),[mfilename,'-',stnm,'.mat']);
  %DEBUG:  disp(mfname);
  if ( exist(mfname,'file') )
    disp(['Load ',mfname]);
    load(mfname,'stn');
  else
    if ( strcmpi(stnm,'lkwf1') )
      stn = optimize_station_heat_budget(stnm,'erai','none','ndbc','tpxo_tide','erai');
    else
      stn = optimize_station_heat_budget(stnm,'erai','avhrr_weekly','ndbc','tpxo_tide','erai');
    end;
    stn = verify_variable(stn,'erai_ndbc_arf_1_d_sum');
    stn = verify_variable(stn,'erai_ndbc_srf_1_d_sum');
    stn = verify_variable(stn,{'ndbc_erai_erai_30a_latent_flux_term_1_d_sum','ndbc_erai_erai_30a_sensible_flux_term_1_d_sum','erai_ndbc_arf_term_1_d_sum','ndbc_erai_erai_30a_net_flux_term_1_d_sum'});
    disp(['Save ',mfname]);
    save(mfname,'stn');
  end;
end;

stn = verify_variable(stn,'ndbc_sea_t_1_d_var_3_d_avg');

stn = verify_variable(stn,'ndbc_wind1_speed_1_d_median');
stn = verify_variable(stn,'ndbc_wind1_u_1_d_median');
stn = verify_variable(stn,'ndbc_wind1_v_1_d_median');
stn = verify_variable(stn,'ndbc_wind1_u_1_d_var_0_d_asof_sum_ndbc_wind1_v_1_d_var');
stn = verify_variable(stn,'ndbc_wind1_u_3_d_var_0_d_asof_sum_ndbc_wind1_v_3_d_var');

stn = verify_variable(stn,'ndbc_barom_1_d_lp');
stn = verify_variable(stn,'ndbc_barom_3_d_lp');
stn = verify_variable(stn,'ndbc_barom_5_d_min');
stn = verify_variable(stn,'ndbc_barom_5_d_max');
stn = verify_variable(stn,'ndbc_barom_5_d_var');
stn.ndbc_barom_5_d_range = ts_op(stn.ndbc_barom_5_d_max,stn.ndbc_barom_5_d_min,'-');


% % These boxplots reveal an odd feature - barometric pressure is bimodal,
% % with a primary peak in winter, but a sharp secondary peak in July! 
% for yr=1993:2013; fh=fmg; boxplot_ts(subset_ts(stn.ndbc_barom,@(x)(find(get_year(x.date)==yr))),[],'mean',true); ylim([1000,1030]); titlename(num2str(yr)); pause; close(fh); end;



%% ISOLATE COLD FRONT "EVENTS"

% Select cold fronts based on one of several criteria

%clear jdlim
%jdlim(ts_boreal_cool(stn.ndbc_barom),1) = prctile(stn.ndbc_barom.data(ts_boreal_cool(stn.ndbc_barom)),98);
%jdlim(ts_boreal_warm(stn.ndbc_barom),1) = prctile(stn.ndbc_barom.data(ts_boreal_warm(stn.ndbc_barom)),98);

%cfld='ndbc_barom';						cts=stn.(cfld); crit=@(x)(x.data>=prctile(cts.data,99));
%cfld='ndbc_barom';						cts=stn.(cfld); crit=@(x)(x.data>=prctile(cts.data,96));
%cfld='erai_barom';						cts=stn.(cfld); crit=@(x)(x.data>=prctile(cts.data,96));
%cfld='ndbc_barom_1_d_lp';					cts=stn.(cfld); crit=@(x)(x.data>=prctile(cts.data,96));
%cfld='ndbc_barom_1_d_lp';					cts=stn.(cfld); crit=@(x)(x.data>=prctile(cts.data,92));
%cfld='ndbc_barom';						cts=stn.(cfld); crit=@(x)(x.data>=jdlim);
%cfld='ndbc_barom_5_d_range';					cts=stn.(cfld); crit=@(x)(x.data>=10);
%cfld='ndbc_barom_5_d_var';					cts=stn.(cfld); crit=@(x)(x.data>=prctile(cts.data,96));
%cfld='ndbc_barom_5_d_var';					cts=stn.(cfld); crit=@(x)(x.data>=prctile(cts.data,92));
%cfld='ndbc_wind1_u_1_d_var_0_d_asof_sum_ndbc_wind1_v_1_d_var';	cts=stn.(cfld); crit=@(x)(x.data>=prctile(cts.data,96));
%cfld='ndbc_wind1_u_3_d_var_0_d_asof_sum_ndbc_wind1_v_3_d_var';	cts=stn.(cfld); crit=@(x)(x.data>=prctile(cts.data,96));
%cfld='ndbc_erai_erai_30a_avhrr_hc_dTdt';			cts=stn.(cfld); crit=@(x)(x.data>=prctile(cts.data,96));
%cfld='ndbc_erai_erai_30a_avhrr_hc_dTdt';			cts=stn.(cfld); crit=@(x)(x.data<=prctile(cts.data,4));
%cfld='ndbc_sea_t_diff'; 					cts=stn.(cfld); crit=@(x)(x.data<=prctile(cts.data,96));
%cfld='ndbc_sea_t_diff'; 					cts=stn.(cfld); crit=@(x)(x.data<=prctile(cts.data,4));
%cfld='ndbc_sea_t_diff'; 					cts=stn.(cfld); crit=@(x)(x.data<=prctile(cts.data,1));
cfld='ndbc_sea_t_1_d_var_3_d_avg';				cts=stn.(cfld); crit=@(x)(x.data>=prctile(cts.data,96));


cfdts = unique(floor(cts.date(crit(cts))));

if (isempty(cfdts)); error('No events found!'); end;

% % Works best for 1 d wind variance
% bakdys=1; fwddys=3;
% % Trying for wind variances vs. sea temperature difference
% bakdys=2; fwddys=5;
% Trying for wind variances vs. sea temperature difference
bakdys=4; fwddys=6;
% % Works best for all others
% bakdys=3; fwddys=3;

% Isolate events (from whatever criterion) into individual non-contiguous days
cfdts(find(diff(cfdts)<3)+1) = [];


% Center each event based on a local extremum of the centering time-series

%cent_fld = 'ndbc_barom';						cent_fn = @min;
%cent_fld = 'ndbc_barom';						cent_fn = @max;
%cent_fld = 'erai_barom';						cent_fn = @min;
%cent_fld = 'erai_barom';						cent_fn = @max;
%cent_fld = 'ndbc_air_t';						cent_fn = @min;
%cent_fld = 'ndbc_air_t';						cent_fn = @max;
%cent_fld = 'ndbc_sea_t_diff';						cent_fn = @min;
%cent_fld = 'ndbc_sea_t_diff';						cent_fn = @max;
%cent_fld = 'ndbc_erai_erai_30a_latent_flux_term';			cent_fn = @min;
%cent_fld = 'ndbc_erai_erai_30a_latent_flux_term';			cent_fn = @max;
%cent_fld = 'ndbc_erai_erai_30a_net_flux_term';				cent_fn = @min;
%cent_fld = 'ndbc_erai_erai_30a_net_flux_term';				cent_fn = @max;
%cent_fld = 'b_ndbc_erai_erai_30a_net_flux_term';			cent_fn = @min;
%cent_fld = 'b_ndbc_erai_erai_30a_net_flux_term';			cent_fn = @max;
%cent_fld = 'ndbc_erai_erai_30a_avhrr_hc_dTdt';				cent_fn = @min;
%cent_fld = 'ndbc_erai_erai_30a_avhrr_hc_dTdt';				cent_fn = @max;
%cent_fld = 'ndbc_erai_erai_30a_avhrr_hc_dTdt_4_d_avg';			cent_fn = @min;
%cent_fld = 'ndbc_erai_erai_30a_avhrr_hc_dTdt_4_d_avg';			cent_fn = @max;
%cent_fld = 'ndbc_barom_5_d_var';					cent_fn = @min;
%cent_fld = 'ndbc_barom_5_d_var';					cent_fn = @min;
%cent_fld = 'ndbc_wind1_u_3_d_var_0_d_asof_sum_ndbc_wind1_v_3_d_var';	cent_fn = @min;
%cent_fld = 'ndbc_wind1_u_3_d_var_0_d_asof_sum_ndbc_wind1_v_3_d_var';	cent_fn = @max;
%cent_fld = cfld;							cent_fn = @min;
cent_fld = cfld;							cent_fn = @max;

stn = verify_variable(stn,cent_fld);

cent_ts = stn.(cent_fld);

add_flds = {'ndbc_erai_erai_30a_avhrr_hc_dTdt_4_d_avg','ndbc_wind1_u_3_d_var_0_d_asof_sum_ndbc_wind1_v_3_d_var'};
add_flds={};
nvars = numel(add_flds);


rdts = -bakdys:(1/24):fwddys;
npts = length(rdts);
vecs=[]; clear vecs
% Make room for centering time series and criterion time series
vecs = repmat(nan,[numel(cfdts),(npts*(nvars+1))]);

badix = [];
for ix = 1:numel(cfdts)
  evtix = find(cfdts(ix)-bakdys<=cent_ts.date & cent_ts.date<=cfdts(ix)+fwddys);
  if ( isempty(evtix) )
    % If no centering TS data during this event, remove it from consideration
    badix(end+1) = ix;
  else
    % Glean from time series that may not be strictly hourly
    [ig,hrix] = unique(get_yearhour(cent_ts.date(evtix)));
    %DEBUG:  if (numel(hrix)<numel(evtix)); keyboard; end;
    evtix = evtix(hrix);

    % Find extremum
    [ig,peakix] = cent_fn(cent_ts.data(evtix));
    % Recenter event at extremum
    cfdts(ix) = cent_ts.date(evtix(peakix));
    evtix = find(cfdts(ix)-bakdys<=cent_ts.date & cent_ts.date<=cfdts(ix)+fwddys);

    % Append centering time-series at end of data vector for SOM analysis
    datix = unique(round(interp1(cfdts(ix)+rdts,1:numel(rdts),get_yearhour(cent_ts.date(evtix)))));
    evtix(~isfinite(datix)) = [];
    datix(~isfinite(datix)) = [];
    vecs(ix,(npts*nvars)+datix) = cent_ts.data(evtix);

    % Catenate time-series of interest into data vector
    for varix = 1:nvars
      stn = verify_variable(stn,add_flds{varix});
      ts = stn.(add_flds{varix});
      evtix = find(cfdts(ix)-bakdys<=ts.date & ts.date<=cfdts(ix)+fwddys);
      [ig,hrix] = unique(get_yearhour(ts.date(evtix)));
      datix = unique(round(interp1(cfdts(ix)+rdts,1:numel(rdts),get_yearhour(ts.date(evtix)))));
      evtix(~isfinite(datix)) = [];
      datix(~isfinite(datix)) = [];
      vecs(ix,(npts*(varix-1))+datix) = ts.data(evtix);
    end;
  end;
end;
cfdts(badix) = [];
vecs(badix,:) = [];
clear badix


%% PLOT COLD FRONT EVENTS AS SHORT TIME SERIES

if (0)

if ( ~exist('fh','var') || ~ishandle(fh) )
  fh = [];
end;
if ( isempty(fh) )
  fh=fmg;

  % % % subplot(3,1,1); plot_ts(stn.ndbc_air_t,stn.ndbc_sea_t); xlabel('Ta, Ts'); grid on;
  % % subplot(3,1,1); plot_ts(stn.ndbc_air_t,stn.ndbc_sea_t,ts_op(stn.erai_spechumid,1000,'*')); xlabel('Ta, Ts, qa'); grid on;
  % subplot(3,1,1); plot_ts(stn.ndbc_air_t,stn.ndbc_sea_t,stn.ndbc_wind1_speed_1_d_median,ts_op(stn.erai_spechumid,1000,'*')); xlabel('Ta,Ts,U_1_d,qa'); grid on;
  subplot(3,1,1); plot_ts(stn.ndbc_air_t,stn.ndbc_sea_t,stn.ndbc_wind1_speed_1_d_median,ts_op(stn.erai_spechumid,1000,'*'),ts_fun(stn.ndbc_wind1_u_1_d_var_0_d_asof_sum_ndbc_wind1_v_1_d_var,@sqrt),'k-'); xlabel('Ta,Ts,U_1_d,qa,U_v_a_r'); grid on;
  legend('Ta','Ts','med_1_dU','qa','\sigma_1_du+\sigma_1_dv');

  % % subplot(3,1,2); plot_ts(stn.erai_ndbc_arf_1_d_sum,stn.erai_ndbc_srf_1_d_sum); xlabel(',Q_S_W'); grid on;
  % % subplot(3,1,2); plot_ts(stn.ndbc_erai_erai_30a_latent_flux,stn.ndbc_erai_erai_30a_sensible_flux,stn.ndbc_erai_erai_30a_net_flux); xlabel('Q_L_H,Q_S_H,Q_0'); grid on;
  % % stn.erai_ndbc_arf_term_dly \Sigma_1_d(Q_S_W(\gamma)+Q_L_W)
  % subplot(3,1,2); plot_ts(stn.ndbc_erai_erai_30a_latent_flux_term_1_d_sum,stn.ndbc_erai_erai_30a_sensible_flux_term_1_d_sum,stn.erai_ndbc_arf_term_1_d_sum,stn.ndbc_erai_erai_30a_net_flux_term_1_d_sum); xlabel('\Sigma_1_dQ_L_H,Q_S_H,Q_S_W(\gamma)+Q_L_W,Q_0'); grid on;
  subplot(3,1,2); plot_ts(stn.ndbc_erai_erai_30a_latent_flux,stn.ndbc_erai_erai_30a_sensible_flux,stn.erai_ndbc_arf,stn.ndbc_erai_erai_30a_net_flux); xlabel('Q_L_H,Q_S_H,Q_S_W(\gamma)+Q_L_W,Q_0'); grid on;
  legend('Q_L_H','Q_S_H','Q_S_W(\gamma)+Q_L_W','Q_0');

  % subplot(3,1,3); plot_ts(stn.ndbc_barom); xlabel('Pr'); grid on;
  subplot(3,1,3); plot_ts(stn.ndbc_barom,stn.ndbc_barom_1_d_lp,stn.ndbc_barom_3_d_lp,'-',stn.erai_barom,'k-'); xlabel('Pr'); grid on;

end;

if ( ~exist('dtix','var') )
  %dtix = 1;
  dtix = ceil(numel(cfdts)/2);
end;
for dtix = dtix:numel(cfdts)
  figure(fh);
  dt = cfdts(dtix);
  % % % %subplot(3,1,1); axis([dt-7,dt+6.9,0,40]); datetick3;
  % % % subplot(3,1,1); axis([dt-7,dt+6.9,-20,25]); datetick3;
  % % % titlename([num2str(dtix),' of ',num2str(numel(cfdts))]);
  % % % subplot(3,1,2); axis([dt-7,dt+6.9,5,29]); datetick3;
  % % % subplot(3,1,3); axis([dt-7,dt+6.9,1005,1033]); datetick3;
  % % subplot(3,1,1); axis([dt-7,dt+6.9,-600,600]); datetick3;
  % % subplot(3,1,2); axis([dt-7,dt+6.9,0,5500]); datetick3;
  % subplot(3,1,1); axis([dt-7,dt+6.9,5,33]); datetick3;
  subplot(3,1,1); axis([dt-7,dt+6.9,0,33]); datetick3;
  titlename([num2str(dtix),' of ',num2str(numel(cfdts))]);
  % subplot(3,1,2); axis([dt-7,dt+6.9,-8,8]); datetick3;
  subplot(3,1,2); axis([dt-7,dt+6.9,-1000,+1000]); datetick3;
  % subplot(3,1,3); axis([dt-7,dt+6.9,1002,1033]); datetick3;
  subplot(3,1,3); axis([dt-7,dt+6.9,1005,1030]); datetick3;
  pause;
%  break;
  if ( ~ishandle(fh) )
    disp('Quit');
    break;
  end;
end;

end; %if (0)


%[sm,pc,sc,framedts,fldnms,fhs] = som_vs_pca_ts(stn,periods,ndays,mapdims,normalrng,nodisplay,fldnms,fldabbr);
%[sm,fh] = som_tses(tses_or_vecs,mapdims,npts,fdts,fldnms,ttl,doSort)
% [sm,fh] = som_tses(vecs,[5,5],npts,[],[],...
%                    strrep([upper(stnm),' ',upper(cent_fld),'(b) based on ',upper(cfld),'(g)'],'_','\_'),false,false);
% [sm,fh] = som_tses(vecs,[5,5],npts,[],[],...
%                    strrep([upper(stnm),' ',upper(cent_fld),'(b) based on ',upper(cfld),'(g)'],'_','\_'),true,false);
[sm,fh] = som_tses(vecs,[5,5],npts,[],[],...
                   strrep([upper(stnm),' ',upper(cent_fld),'(b) based on ',upper(cfld),'(g)'],'_','\_'),false,true);
[sm,fh] = som_tses(vecs,[5,5],npts,[],[],...
                   strrep([upper(stnm),' ',upper(cent_fld),'(b) based on ',upper(cfld),'(g)'],'_','\_'),true,true);
