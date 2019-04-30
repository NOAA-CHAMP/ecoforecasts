1;
%function [stn,wres,ures,vres] = regress_kd(stn,kdfld,wfld,ufld,vfld,season,doPlots)

doPlots = 1;

season = 'summer';


switch ( lower(season(1:3)) )
 case 'ann',
  jdays = 1:366;
 case {'sum','wet'},
  jdays = 122:306;
 case {'win','dry'},
  jdays = [1:121 307:366];
 otherwise,
  error('Season %s not yet implemented!', season);
end;

  kdfld = 'kd_1_day_bic_surf_par_bic_shallow_par';
  kdabr = 'Kd_P_A_R_1_m';
%   kdfld = 'kd_1_day_bic_surf_305nm_bic_shallow_305nm';
%   kdabr = 'Kd_U_V_B_1_m';

%   kdfld = 'kd_1_day_bic_surf_par_bic_deep_par';
%   kdabr = 'Kd_P_A_R_5_m';
%   kdfld = 'kd_1_day_bic_surf_305nm_bic_deep_305nm';
%   kdabr = 'Kd_U_V_B_5_m';


%   wfld = 'wind1_speed_3_day_lowpass';
%   wabr = 'W_s_p_d 3dLP';
%   ufld = 'wind1_u_3_day_lowpass';
%   vfld = 'wind1_v_3_day_lowpass';
  wfld = 'wind1_speed_1_day_lowpass';
  wabr = 'W_s_p_d 1dLP';
  ufld = 'wind1_u_1_day_lowpass';
  vfld = 'wind1_v_1_day_lowpass';
%   wfld = 'wind1_speed_3_day_average';
%   wabr = '\mu_3_d(W_s_p_d)';
%   ufld = 'wind1_u_3_day_average';
%   vfld = 'wind1_v_3_day_average';
%   wfld = 'wind1_speed_1_day_average';
%   wabr = '\mu_1_d(W_s_p_d)';
%   ufld = 'wind1_u_1_day_average';
%   vfld = 'wind1_v_1_day_average';

  srvi2 = verify_variable(srvi2, wfld);
  srvi2 = verify_variable(srvi2, ufld);
  srvi2 = verify_variable(srvi2, vfld);

  goodix = find(0.0 < srvi2.(kdfld).data & srvi2.(kdfld).data < 1.0);
  kdvec = datevec(srvi2.(kdfld).date);
  kdjdays = fix(srvi2.(kdfld).date(:)) - datenum(kdvec(:,1),1,1) + 1;
  seasonix = find(ismember(kdjdays,jdays));
  goodix = intersect(goodix, seasonix);

  [kix,wix] = intersect_dates(srvi2.(kdfld).date(goodix), srvi2.(wfld).date);
  [wB wStats] = robustfit(srvi2.(wfld).data(wix), srvi2.(kdfld).data(goodix(kix)));

  if ( doPlots )
    X = { srvi2.(wfld).date(wix), srvi2.(kdfld).date(goodix(kix)) };
    Y = { srvi2.(wfld).data(wix), srvi2.(kdfld).data(goodix(kix)) };
    multiplot_datetick(X, Y, ['SRVI2 Wind vs. ' kdabr ' (' season ')'], [], {wabr,kdabr});
    print('-dpng',['figs/srvi2_' season '_' kdfld '_' wfld '.png']);
    clear X;
    clear Y;
  end;

  if ( doPlots )
    figure;
    hold on;
    set(gcf, 'units','normalized', 'outerposition', [0 0 1 1]);
    plot(srvi2.(wfld).data(wix), srvi2.(kdfld).data(goodix(kix)), '.');
    [ig,minix] = min(srvi2.(wfld).data(wix)); [ig,maxix] = max(srvi2.(wfld).data(wix));
    lh = line( srvi2.(wfld).data(wix([minix maxix])), ...
               (wB(1) + (wB(2)*srvi2.(wfld).data(wix([minix maxix])))) );
    set(lh, 'color', 'red');
    legend(strrep(lower(wfld),'_','\_'), sprintf('%.5g + %.5g*X, p=%.5g, N=%d',wB(1),wB(2),wStats.p(2),length(kix)), 'Location','NorthWest');
    title([wabr ' vs. ' kdabr ' (' season ')' ]);
    print('-dpng', ['figs/srvi2_' season '_' wfld '_vs_' kdfld]);
  end;


  [kix,uix] = intersect_dates(srvi2.(kdfld).date(goodix), srvi2.(ufld).date);
  [uB uStats] = robustfit(srvi2.(ufld).data(uix), srvi2.(kdfld).data(goodix(kix)));

  if ( doPlots )
    figure;
    hold on;
    set(gcf, 'units','normalized', 'outerposition', [0 0 1 1]);
    plot(srvi2.(ufld).data(uix), srvi2.(kdfld).data(goodix(kix)), '.');
    [ig,minix] = min(srvi2.(ufld).data(uix)); [ig,maxix] = max(srvi2.(ufld).data(uix));
    lh = line( srvi2.(ufld).data(uix([minix maxix])), ...
               (uB(1) + (uB(2)*srvi2.(ufld).data(uix([minix maxix])))) );
    set(lh, 'color', 'red');
    legend(strrep(lower(ufld),'_','\_'), sprintf('%.5g + %.5g*X, p=%.5g, N=%d',uB(1),uB(2),uStats.p(2),length(kix)), 'Location','NorthWest');
    title(['W_U vs. ' kdabr ' (' season ')' ]);
    print('-dpng', ['figs/srvi2_' season '_' ufld '_vs_' kdfld]);
  end;


  [kix,vix] = intersect_dates(srvi2.(kdfld).date(goodix), srvi2.(vfld).date);
  [vB vStats] = robustfit(srvi2.(vfld).data(vix), srvi2.(kdfld).data(goodix(kix)));

  if ( doPlots )
    figure;
    hold on;
    set(gcf, 'units','normalized', 'outerposition', [0 0 1 1]);
    plot(srvi2.(vfld).data(vix), srvi2.(kdfld).data(goodix(kix)), '.');
    [ig,minix] = min(srvi2.(vfld).data(vix)); [ig,maxix] = max(srvi2.(vfld).data(vix));
    lh = line( srvi2.(vfld).data(vix([minix maxix])), ...
               (vB(1) + (vB(2)*srvi2.(vfld).data(vix([minix maxix])))) );
    set(lh, 'color', 'red');
    legend(strrep(lower(vfld),'_','\_'), sprintf('%.5g + %.5g*X, p=%.5g, N=%d',vB(1),vB(2),vStats.p(2),length(kix)), 'Location','NorthWest');
    title(['W_V vs. ' kdabr ' (' season ')' ]);
    print('-dpng', ['figs/srvi2_' season '_' vfld '_vs_' kdfld]);
  end;


%   clear kdfld;
%   clear wfld;
  clear ufld;
  clear vfld;
  clear kdvec;
  clear kdjdays;
  clear seasonix;
  clear goodix;
  clear doPlots;
  clear kdabr;
  clear wabr;
  clear uabr;
  clear vabr;
  clear kix;
  clear wix;
  clear uix;
  clear vix;
  clear lh;
  clear minix;
  clear maxix;
  clear ig;

% return;
