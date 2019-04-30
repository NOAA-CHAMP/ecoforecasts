function soy = get_startyear(dts)
%function soy = get_startyear(dts)
% Return DATENUM of start of the calendar year, for each DATENUM in DTS.
% 
% Last Saved Time-stamp: <Fri 2011-04-22 08:20:31  Lew.Gramer>

  [yrs,mos,dys] = datevec(dts);
  soy = datenum(yrs,1,1);

return;
