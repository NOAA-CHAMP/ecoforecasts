1;

fmg;
legh=[];
legs={}; 
cstnms=fieldnames(stns);
for stix=1:numel(cstnms);
  stnm=cstnms{stix};
  ts = [];
  if ( isfield(stns.(stnm),'fknms_seatemp') && station_dist(bnpon,stns.(stnm)) < 100 )
    ts = subset_ts(stns.(stnm).fknms_seatemp,...
                   @(x)(find(datenum(2010,1,1)<x.date&...
                             x.date<datenum(2010,1,31))));
    if ( ~isempty(ts.data) )
      legh(end+1)=plot_ts(stns.(stnm).fknms_seatemp,'.-',...
                          'Color',color_pick(stix),...
                          'Marker',facemarker_pick(stix));
      legs{end+1}=textize(stnm);
    end;
  end;
end;
legend(legh,legs);
xlim(datenum(2010,1,[0,32]));
datetick3;
