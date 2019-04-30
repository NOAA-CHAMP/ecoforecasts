1;

set_more off

%cstnms={'lkwf1','fwyf1','mlrf1','lonf1','smkf1','sanf1','dryf1',};
%cstnms={'fwyf1','mlrf1','lonf1','smkf1','sanf1','dryf1',};
%cstnms={'lkwf1','mlrf1','lonf1',};
cstnms={'mlrf1','lonf1',};

% Response variable
cY = {'Ts 1d', 'ndbc_sea_t_1_d_avg'};

% Control variable
cX = { ...
    'Ta', 'ndbc_air_t' ;...
    'ERAI Ta', 'erai_air_t' ;...
    'ERAI qa', 'erai_spechumid' ;...
    'ERAI QswI', 'erai_dsrf' ;...
    %'ERAI QlwI', 'erai_dlrf' ;...
    'ERAI Q0', 'erai_net_heat_flux' ;...
    %'NARR Ta', 'ncep_air_t' ;...
    %'NARR qa', 'ncep_spechumid' ;...
    %'NARR QswI', 'ncep_dsrf' ;...
    %'NARR QlwI', 'ncep_dlrf' ;...
    %'NARR Q0', 'ncep_net_heat_flux' ;...
     };


if ( ~exist('stns','var') || numel(stns) ~= numel(cstnms) )
  stns=[]; clear stns
  ta=[]; taa=[]; clear ta taa taclim tatid tad
  ts=[]; tsa=[]; clear ts tsa tsclim tstid tsd

  for cix=1:numel(cstnms)
    stns{cix} = get_station_from_station_name(cstnms{cix});
    stns{cix} = load_all_ndbc_data(stns{cix});
    stns{cix} = get_ncep_station(stns{cix},'narr');
    stns{cix} = get_erai_station(stns{cix});
    stns{cix} = adjust_erai_station(stns{cix});
  end; %for cix=1:numel(cstnms)
end; %if ( ~exist('stns','var') || numel(stns) ~= numel(cstnms) )


for xix=1:size(cX,1)

  Xs = cX{xix,1}; X = cX{xix,2};
  Ys = cY{1}; Y = cY{2};

  disp(X);

  for cix=1:numel(stns)
    % For derived variables, e.g., daily averages
    stns{cix} = verify_variable(stns{cix},X);
    stns{cix} = verify_variable(stns{cix},Y);

    % Coincident timestamps for control variable (e.g., Ta) and Ts at this site
    [ta{cix},ts{cix}] = intersect_tses(30/(24*60),stns{cix}.(X),stns{cix}.(Y));
  end; %for cix=1:numel(stns)

  % Coincident dates for X and Ts across all stations
  tas = intersect_tses(30/(24*60),ta{:});
  tss = intersect_tses(30/(24*60),ts{:});

  yrs = unique(get_year(tas{1}.date));

  % Subtract hourly (or daily) climatologies to calculate anomalies
  for cix=1:numel(stns)
    %[taa{cix},taclim{cix},tatid{cix},tad{cix}] = anomalize_ts(tas{cix},@get_jhour_no_leap);
    %[tsa{cix},tsclim{cix},tstid{cix},tsd{cix}] = anomalize_ts(tss{cix},@get_jhour_no_leap);
    [taa{cix},taclim{cix},tatid{cix},tad{cix}] = anomalize_ts(tas{cix},@get_jday_no_leap);
    [tsa{cix},tsclim{cix},tstid{cix},tsd{cix}] = anomalize_ts(tss{cix},@get_jday_no_leap);
  end; %for cix=1:numel(stns)

  % Coincident dates for X and Ts anomalies across all stations
  taas = intersect_tses(30/(24*60),taa{:});
  tsas = intersect_tses(30/(24*60),tsa{:});

  % Plot annual climatologies
  for cix=1:numel(stns)
    fmg;
    plot(tatid{cix},taclim{cix},'r-',tatid{cix},taclim{cix}-tad{cix},'r:',tatid{cix},taclim{cix}+tad{cix},'r:');
    plot(tstid{cix},tsclim{cix},'b-',tstid{cix},tsclim{cix}-tsd{cix},'b:',tstid{cix},tsclim{cix}+tsd{cix},'b:');
    axis([0,366,13,34]); datetick3('x','mmm');
    titlename([upper(stns{cix}.station_name),' climatological ',Xs,' vs. ',Ys]);
  end; %for cix=1:numel(stns)

  % Plot monthly correlations
  for cix=1:numel(stns)
    for mo=1:12;
      moix = find(get_month(taas{cix}.date)==mo);
      N(mo) = numel(moix);
      [r,p] = corrcoef(taas{cix}.data(moix),tsas{cix}.data(moix));
      R(mo) = r(1,2);
      P(mo) = p(1,2);
    end; %for mo
    disp(N);
    disp(R);
    disp(P);

    fmg;
    bar(1:12,R.^2);
    text([1:12],(R.^2)+0.05,num2str(roundn(P',-2)));
    axis([0.5,12.5,0,1]);
    annotline([],0.5);
    titlename([upper(stns{cix}.station_name),' ',Xs,' vs. ',Ys,' R^2 ',num2str(min(yrs)),'-',num2str(max(yrs))]);
  end; %for cix=1:numel(stns)

end; %for xix=1:size(cX,1)


clear cstnms ans cix cX cY mo moix n p r N P R stnm xix yix X Y Xs Ys yr yrs


set_more
