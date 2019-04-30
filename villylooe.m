1;
%%%% Output a tab-delimited file with a contiguous hourly time series of
%%%% depth-averaged, near-surface, and near-bottom U and V currents from the
%%%% NOAA AOML South Florida Program Looe Key Reef mooring, for Dr. Villy
%%%% Kourafalou's analysis on the coastal satellite alimetry project.
%%%%
%%%% Code by: Lew.Gramer@noaa.gov, 2011 Aug 19
%%%% CALLS: GET_LOOE1_ADCP (ICON/MATLAB Ecoforecasts Toolkit)

if ( ~isfield(stn,'adcp_u') )
  stn = get_looe1_adcp;
end;

begdt = datenum(2007,3,1);
enddt = datenum(2008,5,1) - (0.1/24);
goodix = find( ismember(stn.adcp_sfc_bin',[25:30]) & ...
               (begdt<=stn.adcp_u.date) & (stn.adcp_u.date<=enddt) );

sfc_bin = stn.adcp_sfc_bin(goodix)';
sfc_hgt = stn.adcp_sfc_height(goodix)';
for ix = 1:length(sfc_bin)
  % Get rid of side-lobe contamination (see GET_LOOE1_ADCP code comments)
  last_good_bin(ix,1) = sfc_bin(ix) - (floor(0.06*sfc_bin(ix)) + 1);
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

u_avg.date = u.date;
u_sfc.date = u.date;
u_btm.date = u.date;
for ix = 1:numel(u_sfc.date)
    sfc_ix = last_good_bin(ix);
    u_avg.data(ix,1) = nanmean(u.prof(ix,2:sfc_ix));
    % Uppermost 3m (4 bins) of the clean profile
    u_sfc.data(ix,1) = nanmean(u.prof(ix,sfc_ix-3:sfc_ix));
    % Lowest 3m (4 bins) above presumed bottom boundary layer
    u_btm.data(ix,1) = nanmean(u.prof(ix,2:5));
end;

v_avg.date = v.date;
v_sfc.date = v.date;
v_btm.date = v.date;
for ix = 1:numel(v_sfc.date)
    sfc_ix = last_good_bin(ix);
    v_avg.data(ix,1) = nanmean(v.prof(ix,2:sfc_ix));
    % Uppermost 3m (4 bins) of the clean profile
    v_sfc.data(ix,1) = nanmean(v.prof(ix,sfc_ix-3:sfc_ix));
    % Lowest 3m (4 bins) above presumed bottom boundary layer
    v_btm.data(ix,1) = nanmean(v.prof(ix,2:5));
end;


% Simple QA - ASSUMES U is actually *cross-shore* component
%badix = find(nanmax(abs(u.prof),[],2)>0.5);
badix = find( ~isfinite(u_sfc.data) | (nanmax(abs(u.prof),[],2)>0.5) );
if ( ~isempty(badix) )
  disp(['Clearing ' num2str(numel(badix)) ' bad records']);
  sfc_bin(badix) = nan;
  sfc_hgt(badix) = nan;
  last_top_bin(badix) = nan;
  u_avg.data(badix) = nan;
  u_sfc.data(badix) = nan;
  u_btm.data(badix) = nan;
  v_avg.data(badix) = nan;
  v_sfc.data(badix) = nan;
  v_btm.data(badix) = nan;
end;


%%%% Make some nice plots of results

fmg; plot(u.date,sfc_hgt,'.'); datetick3('x',2,'keeplimits','keepticks'); titlename('Looe Key ADCP Surface Height Estimate'); xlim(u.date([1 end]));
fmg; plot_ts(u,u_sfc,u_btm); legend('Depth-Averaged','Near-Surface','Near-Bottom','Location','Best'); titlename('Looe Key ADCP U'); xlim(u.date([1 end]));
fmg; plot_ts(v,v_sfc,v_btm); legend('Depth-Averaged','Near-Surface','Near-Bottom','Location','Best'); titlename('Looe Key ADCP V'); xlim(u.date([1 end]));

fname = 'villylooe.dat';
fid = fopen(fname,'w+');
if ( isempty(fid) || fid < 3 )
  error('Unable to open %s for output!',fname);
end;

fprintf(fid,'%%%% NOAA SFP PROCESSED MOORED ADCP DATA\n');
fprintf(fid,'%%%% Processed by Lew.Gramer@noaa.gov on %s Eastern\n',datestr(now));
fprintf(fid,'%%%% \n');
fprintf(fid,'%%%% This data set is a product of NOAA/AOML and UM/RSMAS located in Miami, Florida.\n');
fprintf(fid,'%%%% Data contained in this file are intended for use with SoFLA-HYCOM efforts and for\n'); 
fprintf(fid,'%%%% general scientific interest.  NOAA/AOML and UM/RSMAS cannot be held liable for\n');
fprintf(fid,'%%%% the use of these data in any other manner.\n');
fprintf(fid,'%%%% \n');
fprintf(fid,'%%%% =================================================================== \n');
fprintf(fid,'%%%% Year\tJDay\tHour\tMinute\tSfc (m)\tU avg (m/s)\tU sfc (m/s)\tU btm (m/s)\tV avg (m/s)\tV sfc (m/s)\tV btm (m/s)\n');

for ix = 1:numel(u.date)
  fprintf(fid,'%d\t%d\t%d\t%d\t%.4g\t%+.5f\t%+.5f\t%+.5f\t%+.5f\t%+.5f\t%+.5f\n',...
          get_year(u.date(ix)),get_jday(u.date(ix)),...
          get_hour(u.date(ix)),get_minute(u.date(ix)),...
          sfc_hgt(ix),...
          u.data(ix),u_sfc.data(ix),u_btm.data(ix),...
          v.data(ix),v_sfc.data(ix),v_btm.data(ix)...
          );
end;
fclose(fid);

disp(['Output in ' fname]);
clear ans begdt enddt fid fname ix sfc_ix
