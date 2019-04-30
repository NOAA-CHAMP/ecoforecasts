function stn = station_monavg(stn,tfld,nffld,srfld,lrfld,lffld,sffld,rffld)

  if ( ~exist('tfld','var') || isempty(tfld) )
    tfld = 'ndbc_sea_t';
  end;
  if ( ~exist('nffld','var') || isempty(nffld) )
    % nffld = 'ndbc_ncep_30a_total_heat_flux_term';
    nffld = 'netqf';
  end;
  if ( ~exist('srfld','var') || isempty(srfld) )
    srsrc = 'ncep_srf';
    srfld = [srsrc '_term'];
    stn = station_heat_flux_term(stn,srsrc,srfld,'ndbc_sea_t',[],'tmd_tide_i_depth');
  end;
  if ( ~exist('lrfld','var') || isempty(lrfld) )
    lrsrc = 'ncep_lrf';
    lrfld = [lrsrc '_term'];
    stn = station_heat_flux_term(stn,lrsrc,lrfld,'ndbc_sea_t',[],'tmd_tide_i_depth');
  end;
  if ( ~exist('lffld','var') || isempty(lffld) )
    lfsrc = 'ndbc_ncep_30a_latent_heat_flux';
    lffld = [lfsrc '_term'];
    stn = station_heat_flux_term(stn,lfsrc,lffld,'ndbc_sea_t',[],'tmd_tide_i_depth');
  end;
  if ( ~exist('sffld','var') || isempty(sffld) )
    sfsrc = 'ndbc_ncep_30a_sensible_heat_flux';
    sffld = [sfsrc '_term'];
    stn = station_heat_flux_term(stn,sfsrc,sffld,'ndbc_sea_t',[],'tmd_tide_i_depth');
  end;
  if ( ~exist('rffld','var') || isempty(rffld) )
    rfsrc = 'ndbc_ncep_30a_rain_heat_flux';
    rffld = [rfsrc '_term'];
    if ( isfield(stn,rfsrc) )
      stn = station_heat_flux_term(stn,rfsrc,rffld,'ndbc_sea_t',[],'tmd_tide_i_depth');
    end;
  end;

  for mo = 1:12
    t(mo) = station_monavg_mean(stn,tfld,mo);
    n(mo) = station_monavg_mean(stn,nffld,mo);
    rf(mo) = station_monavg_mean(stn,srfld,mo) ...
             + station_monavg_mean(stn,lrfld,mo);
    tf(mo) = station_monavg_mean(stn,lffld,mo) ...
             + station_monavg_mean(stn,sffld,mo);
             
    if ( isfield(stn,rffld) )
      tf(mo) = tf(mo) + station_monavg_mean(stn,rffld,mo);
    end;
  end;

  figure;
%   [ax] = plotyy(1:12, t, 1:12, [n ; rf ; tf]);
  [ax] = plotyy(1:12, t, 1:12, [n(12) + cumsum(n)].*(30*24));
%   [ax] = plotyy(1:12, t, 1:12, [n.*(30*24)]);
  maxigraph;
  xlim(ax(1),[1 12]);
  xlim(ax(2),[1 12]);
  lh1 = legend(ax(1),'T_s', 'Location','NorthWest');
  lh2 = legend(ax(2),'Q_N_E_T','Q_R_A_D','Q_T_U_R', 'Location','NorthEast');

return;

function t = station_monavg_mean(stn,tfld,mo)
  yrs = [1998 1999 2003 2004 2008];

  ix = find( (get_month(stn.(tfld).date) == mo) ...
             & ismember(get_year(stn.(tfld).date), yrs) );
  dat = real(stn.(tfld).data(ix));
  dat(abs(dat) > 2500) = [];
  t = nanmean(dat);
return;
