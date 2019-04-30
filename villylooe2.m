1;
%%%% Output a tab-delimited file with a contiguous hourly time series of
%%%% depth-averaged, near-surface, and near-bottom U and V currents from the
%%%% NOAA AOML South Florida Program Looe Key Reef mooring, for Dr. Villy
%%%% Kourafalou's analysis on the coastal satellite alimetry project.
%%%%
%%%% Code by: Lew.Gramer@noaa.gov, 2011 Aug 19
%%%% CALLS: GET_LOOE1_ADCP (ICON/MATLAB Ecoforecasts Toolkit)

datapath = get_ecoforecasts_path('data');
figspath = get_ecoforecasts_path('figs');

if ( ~exist('stn','var') )
  stn = get_station_from_station_name('looe1');
end;
if ( ~isfield(stn,'adcp_u') )
  stn = get_looe1_adcp;
end;

% begdt = datenum(2008, 1, 1);
% enddt = datenum(2009, 1,25) - (0.1/24);
begdt = datenum(2009, 5, 1);
enddt = datenum(2010, 3, 1) - (0.1/24);
goodix = find( ismember(stn.adcp_sfc_bin,[25:30]) & ...
               (begdt<=stn.adcp_u.date) & (stn.adcp_u.date<=enddt) );

sfc_bin = stn.adcp_sfc_bin(goodix)';
sfc_hgt = stn.adcp_sfc_height(goodix)';

if ( isfield(stn,'adcp_last_good_bin') )
  last_good_bin = stn.adcp_last_good_bin(goodix);
else
  warning('No field STN.adcp_last_good_bin??');
  for ix = 1:length(sfc_bin)
    % Get rid of side-lobe contamination (see GET_LOOE1_ADCP code comments)
    last_good_bin(ix,1) = sfc_bin(ix) - (floor(0.06*sfc_bin(ix)) + 1);
  end;
end;

u = stn.adcp_u;
u.date = u.date(goodix);
u.data = u.data(goodix);
u.prof = u.prof(goodix,:);
u.rawprof = u.rawprof(goodix,:);

v = stn.adcp_v;
v.date = v.date(goodix);
v.data = v.data(goodix);
v.prof = v.prof(goodix,:);
v.rawprof = v.rawprof(goodix,:);

% Cross-shore component
x = stn.adcp_x;
x.date = x.date(goodix);
x.data = x.data(goodix);
x.prof = x.prof(goodix,:);
x.rawprof = x.rawprof(goodix,:);


u_avg.date = u.date;
u_avg.data = u.data;

u_sfc.date = stn.adcp_sfc_u.date(goodix);
u_sfc.data = stn.adcp_sfc_u.data(goodix);

u_btm.date = stn.adcp_btm_u.date(goodix);
u_btm.data = stn.adcp_btm_u.data(goodix);


v_avg.date = v.date;
v_avg.data = v.data;

v_sfc.date = stn.adcp_sfc_v.date(goodix);
v_sfc.data = stn.adcp_sfc_v.data(goodix);

v_btm.date = stn.adcp_btm_v.date(goodix);
v_btm.data = stn.adcp_btm_v.data(goodix);


% Simple QA - assumes *cross-shore* component must generally be quite small
%badix = find(nanmax(abs(u_avg.prof),[],2)>0.5);
badix = find( ~isfinite(u_sfc.data) | (nanmax(abs(x.prof),[],2)>1.0) );
if ( ~isempty(badix) )
  disp(['Clearing ' num2str(numel(badix)) ' bad records']);
  sfc_bin(badix) = [];
  sfc_hgt(badix) = [];
  last_good_bin(badix) = [];

  u_avg.date(badix) = [];
  u_sfc.date(badix) = [];
  u_btm.date(badix) = [];
  v_avg.date(badix) = [];
  v_sfc.date(badix) = [];
  v_btm.date(badix) = [];

  u_avg.data(badix) = [];
  u_sfc.data(badix) = [];
  u_btm.data(badix) = [];
  v_avg.data(badix) = [];
  v_sfc.data(badix) = [];
  v_btm.data(badix) = [];
end;


%%%% Make some nice plots of results

fmg; plot(u_avg.date,sfc_hgt,'.'); datetick3('x',2,'keeplimits','keepticks'); titlename('Looe Key ADCP Surface Height Estimate'); xlim(u_avg.date([1 end]));
fmg; plot_ts(u,u_sfc,u_btm); legend('Depth-Averaged','Near-Surface','Near-Bottom','Location','Best'); titlename('Looe Key ADCP U'); xlim(u_avg.date([1 end]));
% print('-dpng',fullfile(figspath,'villylooe_2008_2009_u.png'));
print('-dpng',fullfile(figspath,'villylooe_2009_2010_u.png'));

fmg; plot_ts(v,v_sfc,v_btm); legend('Depth-Averaged','Near-Surface','Near-Bottom','Location','Best'); titlename('Looe Key ADCP V'); xlim(u_avg.date([1 end]));
% print('-dpng',fullfile(figspath,'villylooe_2008_2009_v.png'));
print('-dpng',fullfile(figspath,'villylooe_2009_2010_v.png'));


% fname = fullfile(datapath,'villylooe_2008_2009.dat');
fname = fullfile(datapath,'villylooe_2009_2010.dat');
fid = fopen(fname,'w+');
if ( isempty(fid) || fid < 3 )
  error('Unable to open %s for output!',fname);
end;

fprintf(fid,'%%%% NOAA SFP PROCESSED MOORED ADCP DATA\n');
fprintf(fid,'%%%% SFP Project PI is Dr. Libby Johns, NOAA AOML\n');
fprintf(fid,'%%%% Processed by Lew.Gramer@noaa.gov on %s Eastern\n',datestr(now));
fprintf(fid,'%%%% \n');
fprintf(fid,'%%%% This data set is a product of NOAA/AOML and UM/RSMAS located in Miami, Florida.\n');
fprintf(fid,'%%%% Data contained in this file are intended for use with SoFLA-HYCOM efforts and for\n'); 
fprintf(fid,'%%%% general scientific interest.  NOAA/AOML and UM/RSMAS cannot be held liable for\n');
fprintf(fid,'%%%% the use of these data in any other manner.\n');
fprintf(fid,'%%%% \n');
fprintf(fid,'%%%% U are west-to-east vector components; V are south-to-north vector components.\n');
fprintf(fid,'%%%% U/V avg values are over all bins from bottom (bin 1) to highest good bin (just\n');
fprintf(fid,'%%%% below side-lobe contamination, 6%% + 1 bin below first surface spike); U/V sfc\n');
fprintf(fid,'%%%% are uppermost six bins below side-lobe contamination; U/V btm deepest six bins.\n');
fprintf(fid,'%%%% =================================================================== \n');
fprintf(fid,'%%%% Year\tJDay\tHour\tMinute\tSfc (m)\tU avg (m/s)\tU sfc (m/s)\tU btm (m/s)\tV avg (m/s)\tV sfc (m/s)\tV btm (m/s)\n');

for ix = 1:numel(u_avg.date)
  fprintf(fid,'%d\t%d\t%d\t%d\t%.4g\t%+.5f\t%+.5f\t%+.5f\t%+.5f\t%+.5f\t%+.5f\n',...
          get_year(u_avg.date(ix)),get_jday(u_avg.date(ix)),...
          get_hour(u_avg.date(ix)),get_minute(u_avg.date(ix)),...
          sfc_hgt(ix),...
          u_avg.data(ix),u_sfc.data(ix),u_btm.data(ix),...
          v_avg.data(ix),v_sfc.data(ix),v_btm.data(ix)...
          );
end;
fclose(fid);

disp(['Output in ' fname]);
clear ans begdt enddt fid fname ix sfc_ix
