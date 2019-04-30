1;

set_more off;

if ( ~exist('ts','var') )
  load(fullfile(get_ecoforecasts_path('data'),'T_tses.mat'),'ts');
end;

% diary(fullfile(get_ecoforecasts_path('data'),'T_tses.csv'));
for tix=1:numel(ts)
  % disp(ts(tix).station_name);
  dt = median(diff(ts(tix).date));
  mindts = 0.70 * (90/dt);

  yrs{tix} = [];
  yrn{tix} = [];
  for yr = 1984:2017;
    dtix = find(get_season(ts(tix).date)==3 & get_year(ts(tix).date)==yr);
    if ( numel(dtix) > mindts )
      yrs{tix} = [yrs{tix},yr];
      yrn{tix} = [yrn{tix},numel(unique(floor(ts(tix).date(dtix))))];
    end;
  end;
  % disp(yrs{tix});
  dp = ts(tix).ngdc_depth;
  if isempty(dp); dp = -ts(tix).depth; end;
  disp(sprintf('"%s",%g,%g,%g,%g%s',ts(tix).station_name,ts(tix).lon,ts(tix).lat,...
               dp,ts(tix).ngdc_beta,num2str(yrs{tix},',%d')));
end;
% diary off;

set_more;
