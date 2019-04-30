1;

%tgtdt = datenum(2010,1,12,2,30,0); %Cold snap
tgtdt = datenum(2011,1,12,2,30,0); %Warm winter
%tgtdt = datenum(1994,8,8,22,00,0); %Cool summer
%tgtdt = datenum(1998,8,14,18,00,0); %Mass bleaching

Tnms={}; Tlons=[]; Tlats=[]; Ts=[];
for stix=1:numel(sites.stnms)
  stnm = sites.stnms{stix};
  if ( isfield(stns.(stnm),'fknms_seatemp') )
    [err,ix] = min(abs(stns.(stnm).fknms_seatemp.date-tgtdt));
    if ( err <= 2/24 && ~isnan(stns.(stnm).fknms_seatemp.data(ix)) )
      Tnms{end+1} = stnm;
      Tlons(end+1) = sites.lons(stix);
      Tlats(end+1) = sites.lats(stix);
      Ts(end+1) = stns.(stnm).fknms_seatemp.data(ix);
    end;
  end;
end;
clear stix stnm err ix

T=[]; clear T
T = griddata(Tlons,Tlats',Ts,LON,LAT,'cubic');

fmg; contourf(LON,LAT,T); colorbar; plot(Tlons,Tlats,'r.');
axis([-82.25,-80.00,24.40,25.65]);
daspect([1,cosd(LAT(1)),1]);
titlename([datestr(tgtdt),' N=',num2str(numel(Ts))]);
