1;

doPlot = true;

if ( ~exist('stnm','var') || isempty(stnm) )
  %stnm = 'fwyf1';
  stnm = 'mlrf1';
  %stnm = 'lonf1';
  %stnm = 'smkf1';
  %stnm = 'looe1';
end;

if ( ~exist('stn','var') || ~isfield(stn,'station_name') || ~strcmpi(stnm,stn.station_name) )
  stn=[]; clear stn;
  stn = get_station_from_station_name(stnm);
end;

%fname = fullfile(get_ecoforecasts_path('data'),'usf',[upper(stnm),'_modkd_time_series.txt']);
fname = fullfile(get_ecoforecasts_path('data'),'usf',[lower(stnm),'_time_series.txt']);
fid = fopen(fname,'r');
if ( fid < 1 )
  warning('Unable to open for read "%s"',fname);
else

  fgetl(fid); % Dump the header line
  x = fscanf(fid,' %f %f %f %f %f %f %f %f %f %f %f %f\n');
  fclose(fid); clear fid

  n = numel(x)/12;
  if ( n ~= floor(n) )
    error('Bad file format "%s"',fname);
  end;
  dat = reshape(x,[12,n])';
  x=[]; clear x

  % YEAR  JDAY  kd412  kd443  kd469  kd488  kd531  kd547  kd555  kd645  kd667  kd678
  flds = {'kd412','kd443','kd469','kd488','kd531','kd547','kd555','kd645','kd667','kd678'};

  dts = datenum(dat(:,1),1,0)+dat(:,2);

  tses={};
  for fldix = 1:numel(flds)
    fld = ['usf_',flds{fldix}];
    stn.(fld).date = dts(:);
    stn.(fld).data = dat(:,fldix+2);
    tses{end+1} = stn.(fld);
  end;
  stn.usf_kd.date = dts(:);
  % %% All colors
  % %stn.usf_kd.data = nanmean(dat(:,3:end),2);
  % % All colors above orange (412-555 nm)
  % stn.usf_kd.data = nanmean(dat(:,3:end-3),2);
  % All colors except 645nm
  stn.usf_kd.data = nanmean(dat(:,[3:end-3,end-1,end]),2);

  %DEBUG
  if ( doPlot )
    fmg; hs=boxplot_tses(tses); ylim([0,1.0]);
    titlename([upper(stnm),' Satellite Kd']);
    h=[]; for hix=1:numel(hs); h(end+1)=hs{hix}(1); end;
    legend(h,flds); 
  end;

  %tses=[]; dat=[]; dts=[]; ans=[]; clear tses hs h dat dts flds fldix fld fname hix n stnm ans
  %dat=[]; dts=[]; ans=[]; clear hs h dat dts flds fldix fld fname hix n stnm ans

end;


% Fit an annual sine-wave model using least squares

% %scatter_curve_fit(get_yearday(stn.usf_kd.date),stn.usf_kd.data,'sin2'),
% f = fittype('(a*cos(((b*x)+c)*2*pi/366)) + d');
% scatter_curve_fit(get_yearday(stn.usf_kd.date),stn.usf_kd.data,f),
f = fittype('(a*cos((x-c)*2*pi/366)) + d');
scatter_curve_fit(get_yearday(stn.usf_kd.date),stn.usf_kd.data,f),
