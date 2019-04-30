function [stn,t,n,rf,tf] = station_annavg(stn,tfld,nffld,srfld,lrfld,lffld,sffld,rffld)

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
    % if ( ~isfield(stn,srfld) )
      stn = station_heat_flux_term(stn,srsrc,srfld,'ndbc_sea_t',[],'tmd_tide_i_depth');
    % end;
  end;
  if ( ~exist('lrfld','var') || isempty(lrfld) )
    lrsrc = 'ncep_lrf';
    lrfld = [lrsrc '_term'];
    % if ( ~isfield(stn,lrfld) )
      stn = station_heat_flux_term(stn,lrsrc,lrfld,'ndbc_sea_t',[],'tmd_tide_i_depth');
    % end;
  end;
  if ( ~exist('lffld','var') || isempty(lffld) )
    lfsrc = 'ndbc_ncep_30a_latent_heat_flux';
    lffld = [lfsrc '_term'];
    % if ( ~isfield(stn,lffld) )
      stn = station_heat_flux_term(stn,lfsrc,lffld,'ndbc_sea_t',[],'tmd_tide_i_depth');
    % end;
  end;
  if ( ~exist('sffld','var') || isempty(sffld) )
    sfsrc = 'ndbc_ncep_30a_sensible_heat_flux';
    sffld = [sfsrc '_term'];
    % if ( ~isfield(stn,sffld) )
      stn = station_heat_flux_term(stn,sfsrc,sffld,'ndbc_sea_t',[],'tmd_tide_i_depth');
    % end;
  end;
  if ( ~exist('rffld','var') || isempty(rffld) )
    rfsrc = 'ndbc_ncep_30a_rain_heat_flux';
    rffld = [rfsrc '_term'];
    % if ( ~isfield(stn,srfld) )
      if ( isfield(stn,rfsrc) )
        stn = station_heat_flux_term(stn,rfsrc,rffld,'ndbc_sea_t',[],'tmd_tide_i_depth');
      end;
    % end;
  end;

  yrs = [1998 1999 2002 2003 2004 2008]; %MLRF1 orig
  yrs = [1998 2002 2003 2008]; % MLRF1
  yrs = [1998 1999 2000 2001 2003 2006 2007]; % SMKF1 orig
  yrs = [1998 1999 2003 2006]; % SMKF1
  yrs = [2001 2002 2003 2006 2007 2008]; % FWYF1 orig
  yrs = [2001 2002 2006 2007]; % FWYF1

  yrs = unique(get_year(stn.(nffld).date));

  for yrix = 1:length(yrs)
    yr = yrs(yrix);
    t(yrix) = station_annavg_mean(stn,tfld,yr);
    n(yrix) = station_annavg_mean(stn,nffld,yr);
    rf(yrix) = station_annavg_mean(stn,srfld,yr) ...
             + station_annavg_mean(stn,lrfld,yr);
    tf(yrix) = station_annavg_mean(stn,lffld,yr) ...
             + station_annavg_mean(stn,sffld,yr);
             
    if ( isfield(stn,rffld) )
      tf(yrix) = tf(yrix) + station_annavg_mean(stn,rffld,yr);
    end;
  end;

  figure;
  [ax] = plotyy(yrs(2:end-1), t(2:end-1), yrs(3:end-1), [cumsum(n(2:end-2))].*(365*24));
  maxigraph;
  xlim(ax(1),yrs([2 end-1]));
  xlim(ax(2),yrs([2 end-1]));
  lh1 = legend(ax(1),'T_s', 'Location','NorthWest');
  lh2 = legend(ax(2),'Q_N_E_T','Q_R_A_D','Q_T_U_R', 'Location','NorthEast');

return;

function t = station_annavg_mean(stn,tfld,yr)
  ix = find( get_year(stn.(tfld).date) == yr );
  dat = real(stn.(tfld).data(ix));
  dat(abs(dat) > 2500) = [];
  t = nanmean(dat);
return;
