1;
% SCRIPT load_NF09_mooring_hobos.m: Extract data from raw CSV (or reload MAT
% file) for small-boat weighted mooring and benthic line sensor deployments
% coincident with the Eddy/Internal Wave experiment conducted aboard NOAA R/V
% Nancy Foster's cruise NF0914 in Fall of 2009. Goal was to have time series
% of in-water point measurements at relatively high frequency, to detect the
% impact of breaking internal wave trains on the shallow (<100 m) shelf of
% southeast Florida, that might be the result of remote interaction between
% near-mesoscale eddies and topography further down in the Florida Keys. The
% moorings did detect higher-than-inertial frequency temperature variability
% coincident with passage of a "wall" of cool water visible in satellite SST.
%
% Last Saved Time-stamp: <Fri 2017-10-20 12:55:13 Eastern Daylight Time gramer>

%On laptop MANANNAN: d:/coral/FACE/Eddy09/moorings/
facepath = get_coral_path('FACE');
moorpath = fullfile(facepath,'Eddy09','moorings');
moordatapath = fullfile(moorpath,'Data','Temp & Pressure');

matfname = fullfile(facepath,'NF09_mooring_hobos.mat');

% Mooring locations
brwd20 = dm2degrees([ 26 14.83740 ; -80 03.90250 ]);
% Broward benthic thermistors were retrieved from SHALLOWEST to DEEPEST
brwdbenth = dm2degrees([ ...
    26 14.84240 ; -80 03.89030 ;...
    26 14.84030 ; -80 03.86400 ;...
    26 14.83620 ; -80 03.82670 ;...
    26 14.83150 ; -80 03.78100 ;...
    26 14.82470 ; -80 03.75850 ;...
                   ]);
brwd40 = dm2degrees([ 26 14.81940 ; -80 03.71280 ]);
brwd100= dm2degrees([ 26 14.84820 ; -80 03.14350 ]);

boca20 = dm2degrees([ 26 20.74720 ; -80 03.23160 ]);
% Boca benthic thermistors were retrieved from DEEPEST to SHALLOWEST
bocabenth = dm2degrees([ ...
    26 20.72050 ; -80 03.06560 ;...
    26 20.72620 ; -80 03.11320 ;...
    26 20.78140 ; -80 03.11830 ;...
    26 20.76760 ; -80 03.16730 ;...
    26 20.75110 ; -80 03.20810 ;...
                   ]);
boca40 = dm2degrees([ 26 20.71930 ; -80 03.04700 ]);


if ( exist(matfname,'file') )

  disp(['Loading ' matfname]);
  load(matfname,'locs');

else

  disp('Processing raw CSV files');
  locs=repmat(struct('seatemp',[],'seapres',[],'par',[]),[0 0]);

  % % Broward Line only
  % tnos=[1:2,6:14,31,41,42];

  % % Boca Line only
  % tnos=[15:16,20,22,23:24,43];

  % BOTH Broward and Boca Lines
  tnos=[1:2,6:14,15:16,20,22,23:24,31,41,42,43];

  for ix=1:numel(tnos);
    tno=sprintf('%02d',tnos(ix));
    fname = fullfile(moordatapath,['T?',tno,'*csv']);
    %DEBUG:
    disp(fname);
    %locs(ix) = read_hobo_csv(fullfile('Data','Temp & Pressure',['T?',tno,'*csv']));
    locs(ix) = read_hobo_csv(fname);
  end;
  for ix=1:numel(tnos);
    locs(ix).tno = tnos(ix);
  end;

  % Convert back to units we know!
  lightix = find([locs.tno]==31);
  locs(lightix).light.date = locs(lightix).par.date;
  locs(lightix).light.data = locs(lightix).par.data .* 13.03;

  disp(['Saving ' matfname]);
  save(matfname,'locs');

end;


if 0;
  %%%% RESIDUAL CODE FROM 'anhobos.m' SCRIPT (same directory)
  [fh,lhs,axs] = multiplot_station(locs(lightix),{'seatemp','par'},...
                                   'Broward 20m & 40m Mooring Data(?)',...
                                   [],{'T [^oC]','Light [?], Pressure [dbar]'},...
                                   [datenum(2009,10,27),datenum(2009,11,9)],...
                                   {[26,30],[0,40]},[],'k-');
  hold(axs(2),'on');
  plot_ts(axs(2),locs(13).seapres,'b-');
  hold(axs(1),'on');
  lhs=plot_ts(axs(1),locs([1,13,2:4]).seatemp); set(lhs,'Marker','none');
  datetick3('x',2,'keeplimits');
  
  print('-dtiff','broward_20m_40m_snail.tiff');
  %disp('Figure not printed this time...');
end;

if 0;
  fmg; plot_ts(locs.seatemp); legend(num2str([locs.tno]'));
  fmg; plot_ts(locs(end-1:end).seapres); legend(num2str([locs(end-1:end).tno]'));
  for ix=1:11; fmg; plot_ts(locs(ix).seatemp,locs(end).seapres); legend(num2str(locs(ix).tno),'42: Pres'); end;
  for ixes={[1:2],[3:5],[8:11]}; ixes{:}, end;
  for ixes={[1,lightix],[2:3],[4:6],[8:11]}; fmg; plot_ts(locs(ixes{:}).seatemp,locs(end).seapres); legend(num2str(locs(ixes{:}).tno)); end;
  for ixes={[1,lightix],[2:4:6,13],[8:11,14]}; fmg; plot_ts(locs(ixes{:}).seatemp,locs(end).seapres); legend(num2str([locs(ixes{:}).tno]')); end;
  for ixes={[1:4,13],[5:7,14],[8:11,lightix]}; fmg; plot_ts(locs(ixes{:}).seatemp,locs(end).seapres); legend(num2str([locs(ixes{:}).tno]')); end;
end;
