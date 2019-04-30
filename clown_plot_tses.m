function clown_plot_tses(tses)
%function clown_plot_tses(tses)
% NOT READY FOR PRIME TIME YET!

error('NOT READY FOR PRIME TIME YET!');

if ( ~exist('stn','var') || ~isfield(stn,'ndbc_sea_t') )
  stn = load_all_ndbc_data([],'lonf1');
end;
stnix2010 = find(get_year(stn.ndbc_sea_t.date)==2010&get_month(stn.ndbc_sea_t.date)==1);


% zone =    'Inner Reef'

zones = { 'Inshore',    'Offshore Patch Reef',    'Mid Channel',    'Forereef',};

for czone = zones
zone = czone{:};
disp(zone);

my_subrgns = {'Biscayne','Upper Keys','Middle Keys','Lower Keys'};


srix = find(ismember(sites.subrgn(:),my_subrgns));
ix = find(strcmp(sites.zone(srix),zone));
zix = srix(ix);
% Sort alphabetically by sub-region name
[ig,sortix] = sort(sites.subrgn(zix));
zix = zix(sortix);

ix2010 = find(get_year(sites.sst(zix(1)).date)==2010);

dts = sites.sst(zix(1)).date(ix2010);
dat = repmat(nan,[length(dts) length(zix)]);
for ix=1:length(zix)
  dat(:,ix) = sites.sst(zix(ix)).data(ix2010);
end;

linspec = {'ro','go','bo','mo', 'r*','g*','b*','m*', 'rs','gs','bs','ms', 'r^','g^','b^','m^', 'rv','gv','bv','mv',};

figure;
maxigraph;
titlename(zone);
hold on;
for ix = 1:size(dat,2)
  plot(dts,dat(:,ix),linspec{ix});
end;
plot(stn.ndbc_sea_t.date(stnix2010),stn.ndbc_sea_t.data(stnix2010),'k:','Color',[.5,.5,.5]);
ylim([8 25]);
ylabel('Combined SST [\circC]');
legs = strcat(sites.subrgn(zix),' #',num2str(sites.sno(zix)));
legs{end+1} = 'Long Key SEAKEYS';
legend(legs, 'Location','Best');
datetick3;
grid on;
print('-dpng',[zone '-clown-plot.png']);

end;
