1;
%%%% Output a tab-delimited file with a contiguous hourly time series of
%%%% depth-averaged, near-surface current Speed and Direction from the
%%%% NOAA AOML South Florida Program Looe Key Reef mooring, for analysis by
%%%% Dr. Bill Wurts of U. Kentucky on an aquaculture analysis.
%%%%
%%%% Code by: Lew.Gramer@noaa.gov, 2015 Sep 11 based on code from 2011 Aug 19
%%%% CALLS: GET_LOOE1_ADCP, GET_LOOE1_MICROCAT (Gramer PhD Thesis Toolkit)

datapath = get_ecoforecasts_path('data');
figspath = get_ecoforecasts_path('figs');

if ( ~exist('stn','var') )
  stn = get_station_from_station_name('looe1');
end;
if ( ~isfield(stn,'adcp_speed') )
  stn = get_looe1_adcp;
end;
if ( ~isfield(stn,'microcat_seatemp') )
  stn = get_looe1_microcat(stn);
end;

[adcpix,mcatix] = intersect_dates(stn.adcp_speed.date,stn.microcat_seatemp.date);

sfc_bin = stn.adcp_sfc_bin(adcpix)';
sfc_hgt = stn.adcp_sfc_height(adcpix)';

if ( isfield(stn,'adcp_last_good_bin') )
  last_good_bin = stn.adcp_last_good_bin(adcpix);
else
  warning('No field STN.adcp_last_good_bin??');
  for ix = 1:length(sfc_bin)
    % Get rid of side-lobe contamination (see GET_LOOE1_ADCP code comments)
    last_good_bin(ix,1) = sfc_bin(ix) - (floor(0.06*sfc_bin(ix)) + 1);
  end;
end;

sfc_t = stn.microcat_seatemp;
sfc_t.date = sfc_t.date(mcatix);
sfc_t.data = sfc_t.data(mcatix);

sp = stn.adcp_speed;
sp.date = sp.date(adcpix);
sp.data = sp.data(adcpix);
sp.prof = sp.prof(adcpix,:);
sp.rawprof = sp.rawprof(adcpix,:);

if ( numel(sfc_t.date) ~= numel(sp.date) )
  error('Intersection failed!');
end;

dr = stn.adcp_dir;
dr.date = dr.date(adcpix);
dr.data = dr.data(adcpix);
dr.prof = dr.prof(adcpix,:);
dr.rawprof = dr.rawprof(adcpix,:);

x = stn.adcp_x;
x.date = x.date(adcpix);
x.data = x.data(adcpix);
x.prof = x.prof(adcpix,:);
x.rawprof = x.rawprof(adcpix,:);

sfc_sp.date = stn.adcp_sfc_speed.date(adcpix);
sfc_sp.data = stn.adcp_sfc_speed.data(adcpix);

sfc_dr.date = stn.adcp_sfc_dir.date(adcpix);
sfc_dr.data = stn.adcp_sfc_dir.data(adcpix);

sfc_x.date = stn.adcp_sfc_x.date(adcpix);
sfc_x.data = stn.adcp_sfc_x.data(adcpix);


% Simple QA - assumes *cross-shore* component must generally be quite small
%badix = find(nanmax(abs(u_avg.prof),[],2)>0.5);
badix = find( ~isfinite(sfc_sp.data) | (nanmax(abs(x.prof),[],2)>1.0) );
if ( ~isempty(badix) )
  disp(['Clearing ',num2str(numel(badix)),' bad ADCP and Temp. records']);

  sfc_t.date(badix) = [];
  sfc_t.data(badix) = [];

  sfc_bin(badix) = [];
  sfc_hgt(badix) = [];
  last_good_bin(badix) = [];

  sfc_sp.date(badix) = [];
  sfc_sp.data(badix) = [];

  sfc_dr.date(badix) = [];
  sfc_dr.data(badix) = [];

  sfc_x.date(badix) = [];
  sfc_x.data(badix) = [];
end;


%%%% Make some nice plots of results

fmg; plot(sfc_sp.date,sfc_hgt,'.'); datetick3('x',2,'keeplimits','keepticks'); titlename('Looe Key ADCP Surface Height Estimate'); xlim(sfc_sp.date([1 end]));
fmg; plot_ts(sfc_sp); titlename('Looe Key ADCP Near-Surface Speed'); xlim(sfc_sp.date([1 end]));
print('-dpng',fullfile(figspath,'wurts_looe_key_near_sfc_current.png'));

fmg; plot_ts(sfc_t); titlename('Looe Key Near-Surface Temperature'); xlim(sfc_t.date([1 end]));
print('-dpng',fullfile(figspath,'wurts_looe_key_near_sfc_temp.png'));


fname = fullfile(datapath,'wurts_looe_key_near_sfc_current_and_temp.csv');
fid = fopen(fname,'w+');
if ( isempty(fid) || fid < 3 )
  error('Unable to open %s for output!',fname);
end;

fprintf(fid,'NOAA SFP PROCESSED MOORED ADCP and MICROCAT DATA\n');
fprintf(fid,'SFP Project PI is Dr. Libby Johns of NOAA AOML\n');
fprintf(fid,'Processed by Lew.Gramer@noaa.gov on %s Eastern\n',datestr(now));
fprintf(fid,'\n');
fprintf(fid,'This data set is a product of NOAA/AOML and UM/RSMAS located in Miami Florida.\n');
fprintf(fid,'Data contained in this file are intended for use with SoFLA-HYCOM efforts and for\n'); 
fprintf(fid,'general scientific interest.  NOAA/AOML and UM/RSMAS cannot be held liable for\n');
fprintf(fid,'the use of these data in any other manner.\n');
fprintf(fid,'\n');
fprintf(fid,'U are west-to-east vector components; V are south-to-north vector components.\n');
fprintf(fid,'SPEED/DIR values are over a total of six bins beginning near-surface: highest\n');
fprintf(fid,'GOOD bin below side-lobe contamination + 6%% + 1 bin below first sfc. spike.\n');
fprintf(fid,'=================================================================== \n');
fprintf(fid,'Year,Month,Day,Hour,Minute,Near-sfc temp (oC),Sfc height (m),Sfc speed (m/s),Sfc dir (deg True)\n');

for ix = 1:numel(sfc_sp.date)
  fprintf(fid,'%d,%d,%d,%d,%d,%+.5f,%.4g,%+.5f,%+.5f\n',...
          get_year(sfc_sp.date(ix)),get_month(sfc_sp.date(ix)),get_dom(sfc_sp.date(ix)),...
          get_hour(sfc_sp.date(ix)),get_minute(sfc_sp.date(ix)),...
          sfc_t.data(ix),sfc_hgt(ix),...
          sfc_sp.data(ix),sfc_dr.data(ix)...
          );
end;
fclose(fid);

disp(['Output in ' fname]);
clear ans adcpix mcatix fid fname ix sfc_ix
