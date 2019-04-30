1;

if ( ~exist('stns','var') || ~isfield(stns,'FKNMS_7MILE_BR') )
  anfknms;
end;

clear tses; tses(1) = struct('date',[],'data',[]); stnms = {};

for fld=fieldnames(stns)'; 
  if (isfield(stns.(fld{:}),'fknms_seatemp') && get_year(stns.(fld{:}).fknms_seatemp.date(1))<1995); 
    %DEBUG:    disp(fld{:}); 
    stnms{end+1} = fld{:}; 
  end; 
end;
clear fld

for stix=1:numel(stnms); 
  tses(stix) = stns.(stnms{stix}).fknms_seatemp; 
end;

if 0;
  warning off; X = intersect_tses(1,tses); warning on; X{1}, % :(
  for stix=1:numel(stnms); 
    [ix1,ix2] = intersect_dates(tses(stix).date,stns.MLRF1.fknms_seatemp.date);
    disp({stnms{stix},numel(ix1),});
  end;
  w=warning('off'); X = intersect_tses((4/24),tses(1:end-1)); warning(w); X{1}, % :( :(
end;


set_more off;
clear tyrs dyrs tgap
yrs = [1990:2000,2007:2011];
str = [sprintf('%4d ',yrs),' STATION NAME'];
disp(str);
disp(repmat('#',[1,length(str)]));
flgstr = '';
for stix=1:numel(stnms)
  for yrix=1:numel(yrs)
    yr = yrs(yrix);
    tyrix = find(get_year(tses(stix).date)==yr);
    tyrs(stix,yrix) = numel(tyrix);
    dyrs(stix,yrix) = numel(unique(get_jday(tses(stix).date(tyrix))));
    if ( numel(tyrix) < 2 )
      tgap(stix,yrix) = 365;
    else
      tgap(stix,yrix) = round(max(diff(tses(stix).date(tyrix))));
    end;
  end;
  %disp([sprintf('% 3.2f',(tyrs(stix,:)/(24*365))),' ',stnms{stix}]);
  %disp([sprintf('%4d ',tyrs(stix,:)),' ',stnms{stix}]);

  flgstr = '';
  wantyrs = [1991,1996:1998,2000];
  for wantyr = wantyrs(:)'
    wantyrix = find(wantyr==yrs);
    if ( all(dyrs(stix,wantyrix) > 200) && all(tgap(stix,wantyrix) < 30) )
      flgstr = [flgstr,'*'];
    else
      flgstr = [flgstr,'_'];
    end;
  end;
  if ( numel(strfind(flgstr,'*')) == numel(wantyrs) )
      flgstr = [flgstr,'^'];
  else
      flgstr = [flgstr,' '];
  end;
  disp([sprintf('%4d ',dyrs(stix,:)),' ',flgstr,stnms{stix}]);
end;
set_more;

clear ans flgstr ix1 ix2 str stix yrix 
