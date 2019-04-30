function stn = convert_ndbc_adcp(stn)
%function stn = convert_ndbc_adcp(stn)
%
% The numbered bins in NDBC (old) format ADCP current profiles vary from year
% to year. This processes all bins in STN, returning "adcp_spd", "adcp_dir",
% "adcp_u", and "adcp_v" fields in STN with ".date", ".data" (depth-average)
% and ".prof" (DxH) profile sub-fields, D=# unique dates, H=# unique depths.
%
% NOTE: This can result in profiles which are very difficult to plot! When
% each profile contains multiple NaNs, SURF and CONTOUR, e.g., do not work.
%
% Last Saved Time-stamp: <Sat 2013-10-19 22:06:29 Eastern Daylight Time gramer>

  dt = [];
  dp = [];
  dr = [];
  sp = [];
  for bn=1:50;
    bns = num2str(bn,'%02d');
    if ( isfield(stn,['ndbc_adcp_dir_',bns]) )
      % Convert ADCP speeds from cm/s to m/s
      stn.(['ndbc_adcp_spd_',bns]).data = stn.(['ndbc_adcp_spd_',bns]).data ./ 100;
      % Convert ADCP speed/dir to U and V components
      stn = station_spddir_to_uv(stn,['ndbc_adcp_spd_',bns],['ndbc_adcp_dir_',bns],...
                                   ['ndbc_adcp_u_',bns],['ndbc_adcp_v_',bns],true);
      n = numel(stn.(['ndbc_adcp_spd_',bns]).data);
      dt(end+1:end+n,1) = stn.(['ndbc_adcp_spd_',bns]).date;
      dp(end+1:end+n,1) = stn.(['ndbc_adcp_dep_',bns]).data;
      dr(end+1:end+n,1) = stn.(['ndbc_adcp_dir_',bns]).data;
      sp(end+1:end+n,1) = stn.(['ndbc_adcp_spd_',bns]).data;
    end;
  end;

  %First line of "41012a2002.txt":
  %2002 09 01 00     6   152  19.5     9   140  25.2    13   134  19.6    16   147  14.2    20   118  15.3    23    81  25.0    27   103  26.2    30   105  18.0    34    78   9.8    37    33   2.0

  x = sortrows([dt,dp,dr,sp],[1,2]);
  dt = x(:,1);
  dp = x(:,2);
  dr = x(:,3);
  sp = x(:,4);

  uD = unique(dt(:));
  uH = unique(dp(:));

  D = numel(uD);
  H = numel(uH);

  stn.ndbc_adcp_bin_heights = uH;

  stn.ndbc_adcp_dir.date = uD;
  stn.ndbc_adcp_speed.date = uD;
  stn.ndbc_adcp_u.date = uD;
  stn.ndbc_adcp_v.date = uD;

  stn.ndbc_adcp_dir.data = repmat(nan,[D,1]);
  stn.ndbc_adcp_dir.prof = repmat(nan,[D,H]);

  stn.ndbc_adcp_speed.data = repmat(nan,[D,1]);
  stn.ndbc_adcp_speed.prof = repmat(nan,[D,H]);

  stn.ndbc_adcp_u.data = repmat(nan,[D,1]);
  stn.ndbc_adcp_u.prof = repmat(nan,[D,H]);

  stn.ndbc_adcp_v.data = repmat(nan,[D,1]);
  stn.ndbc_adcp_v.prof = repmat(nan,[D,H]);

  for hix=1:H
    dpix = find(dp == uH(hix));
    dtix = find(ismember(uD,dt(dpix)));

    stn.ndbc_adcp_dir.prof(dtix,hix) = dr(dpix);
    stn.ndbc_adcp_speed.prof(dtix,hix) = sp(dpix);
    [stn.ndbc_adcp_u.prof(dtix,hix),stn.ndbc_adcp_v.prof(dtix,hix)] = spddir_to_uv_curr(sp(dpix),dr(dpix));
  end;

  % Whole water-column averages
  stn.ndbc_adcp_u.data(:) = nanmean(stn.ndbc_adcp_u.prof,2);
  stn.ndbc_adcp_v.data(:) = nanmean(stn.ndbc_adcp_v.prof,2);

  stn.ndbc_adcp_dir.data(:) = uv_to_dir_curr(stn.ndbc_adcp_u.data,stn.ndbc_adcp_v.data);
  stn.ndbc_adcp_speed.data(:) = uv_to_spd(stn.ndbc_adcp_u.data,stn.ndbc_adcp_v.data);

  % Though bin heights (and fieldnames) change with time, water depth shouldn't... much
  wwrg = max(uH) - min(uH);

  % Near-bottom, mid-water, and near-surface averages
  nbrg = [min(uH) + (wwrg*0.67), min(uH) + (wwrg*1.00)];
  mwrg = [min(uH) + (wwrg*0.33), min(uH) + (wwrg*0.67)];
  nsrg = [min(uH) + (wwrg*0.00), min(uH) + (wwrg*0.33)];

  % Bottom 1/3 of sampled water-column
  nbix = find(nbrg(1) <= uH & uH <= nbrg(2));
  stn.ndbc_adcp_btm_bin_depths = uH(nbix);
  stn.ndbc_adcp_btm_u.date = stn.ndbc_adcp_u.date;
  stn.ndbc_adcp_btm_u.data = nanmean(stn.ndbc_adcp_u.prof(:,nbix),2);
  stn.ndbc_adcp_btm_v.date = stn.ndbc_adcp_v.date;
  stn.ndbc_adcp_btm_v.data = nanmean(stn.ndbc_adcp_v.prof(:,nbix),2);

  stn.ndbc_adcp_btm_speed.date = stn.ndbc_adcp_btm_u.date;
  stn.ndbc_adcp_btm_speed.data = uv_to_spd(stn.ndbc_adcp_btm_u.data,stn.ndbc_adcp_btm_v.data);
  stn.ndbc_adcp_btm_dir.date = stn.ndbc_adcp_btm_u.date;
  stn.ndbc_adcp_btm_dir.data = uv_to_dir_curr(stn.ndbc_adcp_btm_u.data,stn.ndbc_adcp_btm_v.data);

  % Middle 1/3 of sampled water-column
  mwix = find(mwrg(1) < uH & uH <= mwrg(2));
  stn.ndbc_adcp_mid_bin_depths = uH(mwix);
  stn.ndbc_adcp_mid_u.date = stn.ndbc_adcp_u.date;
  stn.ndbc_adcp_mid_u.data = nanmean(stn.ndbc_adcp_u.prof(:,mwix),2);
  stn.ndbc_adcp_mid_v.date = stn.ndbc_adcp_v.date;
  stn.ndbc_adcp_mid_v.data = nanmean(stn.ndbc_adcp_v.prof(:,mwix),2);

  stn.ndbc_adcp_mid_speed.date = stn.ndbc_adcp_mid_u.date;
  stn.ndbc_adcp_mid_speed.data = uv_to_spd(stn.ndbc_adcp_mid_u.data,stn.ndbc_adcp_mid_v.data);
  stn.ndbc_adcp_mid_dir.date = stn.ndbc_adcp_mid_u.date;
  stn.ndbc_adcp_mid_dir.data = uv_to_dir_curr(stn.ndbc_adcp_mid_u.data,stn.ndbc_adcp_mid_v.data);

  % Top 1/3 of sampled water-column
  nsix = find(nsrg(1) < uH & uH <= nsrg(2));
  stn.ndbc_adcp_sfc_bin_depths = uH(nsix);
  stn.ndbc_adcp_sfc_u.date = stn.ndbc_adcp_u.date;
  stn.ndbc_adcp_sfc_u.data = nanmean(stn.ndbc_adcp_u.prof(:,nsix),2);
  stn.ndbc_adcp_sfc_v.date = stn.ndbc_adcp_v.date;
  stn.ndbc_adcp_sfc_v.data = nanmean(stn.ndbc_adcp_v.prof(:,nsix),2);

  stn.ndbc_adcp_sfc_speed.date = stn.ndbc_adcp_sfc_u.date;
  stn.ndbc_adcp_sfc_speed.data = uv_to_spd(stn.ndbc_adcp_sfc_u.data,stn.ndbc_adcp_sfc_v.data);
  stn.ndbc_adcp_sfc_dir.date = stn.ndbc_adcp_sfc_u.date;
  stn.ndbc_adcp_sfc_dir.data = uv_to_dir_curr(stn.ndbc_adcp_sfc_u.data,stn.ndbc_adcp_sfc_v.data);

return;
