function yrhrs = get_yearhour(dts)
%function yrhrs = get_yearhour(dts)
% Return DATENUM of start of the hour, for each DATENUM in DTS.
% 
% Last Saved Time-stamp: <Wed 2013-04-03 13:22:48 Eastern Daylight Time gramer>

  [yrs,mos,dys,hrs,mns,scs] = datevec(dts);
  yrhrs = datenum(yrs,mos,dys,hrs,0,0);

return;
