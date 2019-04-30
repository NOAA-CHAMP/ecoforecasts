1;

set_more off

if ( ~exist('doPrint','var') || isempty(doPrint) )
  doPrint = false;
end;

fld = 'ndbc_sea_t';
%fld = 'ndbc_net_flux';
%fld = 'ndbc_air_t';
%fld = 'ndbc_wind1_speed';

%stnms = { 'lonf1','smkf1', };
%stnms = { 'fwyf1','lonf1','mlrf1','smkf1','sanf1', };
stnms = { 'ccrf1','cmdf1','cryf1','hacf1', };
%stnms = { 'fwyf1','lonf1','mlrf1','smkf1','sanf1', 'ccrf1','cmdf1','cryf1','hacf1', };

%ctlprm = 'beta';
ctlprm = 'depth';
prms = [];

if ( ~exist('stns','var') || ~isfield(stns,fld) )
  stns = {}; 
  for stix=1:numel(stnms)
    stnm = stnms{stix};
    switch ( stnm )
     case {'ccrf1','cmdf1','cryf1','hacf1','snaf1',},
      if ( ~exist(stnm,'var') )
        doFigs = false;
        extract_kuffner_hudson;
      end;
      stns{stix} = eval(stnm);
     otherwise,
      stns{stix} = get_station_from_station_name(stnm); 
      stns{stix} = load_all_ndbc_data(stns{stix});
    end;

    if ( ~isempty(strfind(fld,'flux')) )
      stns{stix} = get_erai_station(stns{stix});
      stns{stix} = station_heat_flux(stns{stix},'ndbc_wind1_speed','ndbc_air_t','erai_relhumid','ndbc_barom','ndbc_sea_t',[],[],'ndbc');
    end;

    stns{stix} = station_ngdc_offshore_slope(stns{stix}); 
    switch ( ctlprm ),
     case 'beta',
      prms(stix) = stns{stix}.ngdc_offshore_slope;
     case 'depth',
      prms(stix) = stns{stix}.depth;
    end;
  end;
end;

%peracc = @get_hour; stracc='hour';
peracc = @get_jday; stracc='jday';
%peracc = @get_year; stracc='year';

%ts1 = stns{1}.(fld);
persub = @ts_isfinite; strsub='all';
%persub = @ts_jas; strsub='summer';
%persub = @ts_jfm; strsub='winter';
%persub = @ts_boreal_warm; strsub='rainy';
%persub = @ts_boreal_cool; strsub='dry';

ts1 = subset_ts(stns{1}.(fld),persub);

tses = intersect_tses(ts1,stns{2}.(fld),stns{3}.(fld),stns{4}.(fld),stns{5}.(fld));

%pers = 1:365; % Julian Days, e.g.
pers = unique(peracc(tses{stix}.date));

sz = [numel(pers),numel(stns)];
mn = repmat(nan,sz);
mx = repmat(nan,sz);
mnd = repmat(nan,sz);
mxd = repmat(nan,sz);
p07d = repmat(nan,sz);
p93d = repmat(nan,sz);

for stix=1:numel(stns)
  for ix=1:numel(pers);
    per = pers(ix);
    perix = find(peracc(tses{stix}.date)==per);
    if ( isempty(perix) )
      dat = nan;
    else
      dat = tses{stix}.data(perix);
    end;
    mn(ix,stix) = min(dat);
    mx(ix,stix) = max(dat);

    mnd(ix,stix) = min(dat) - mean(dat);
    mxd(ix,stix) = max(dat) - mean(dat);

    p07d(ix,stix)  = prctile(dat, 7) - mean(dat);
    p93d(ix,stix) = prctile(dat,93) - mean(dat);
  end;
end;

[ig,prmix] = sort(prms,2,'descend');

fmg;
lhs = []; legs = {};
%%peroff = 5e4;
%%peroff = 10e4;
%peroff = 5e3;
%peroff = numel(pers)*3e2;
peroff = numel(pers)*15;
clr = {'k','r',[0.0,0.5,0.0],'b','m'};
for ix=prmix(:)'; 
  % Left-hand (7th percentile) cluster
  lhs(end+1) = plot(p07d(:,ix),prms(ix)+([1:numel(pers)]./peroff),'.','Color',clr{ix});
  % Right-hand (93rd percentile) cluster
  plot(p93d(:,ix),prms(ix)+([1:numel(pers)]./peroff),'.','Color',clr{ix});
  legs{end+1} = upper(stns{ix}.station_name);
end;
for ix=prmix(:)'; 
  text(0,prms(ix)+(1./peroff),[upper(stracc),'=',num2str(pers(1))]);
  text(0,prms(ix)+(numel(pers)./peroff),[upper(stracc),'=',num2str(pers(end))]);
end;
legend(lhs,legs,'Location','West');
titlename([upper(strrep(fld,'_','\_')),' ',upper(strsub),' ',upper(stracc)]);
xlabel(strrep(upper(fld),'_','\_'));
ylabel(strrep(upper(ctlprm),'_','\_'));
if (doPrint); print('-dpng',get_relative_path([ctlprm,'_vs_',fld,'_',strsub,'_',stracc,'.png'])); end;

set_more
