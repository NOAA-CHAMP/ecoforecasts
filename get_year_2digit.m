function yrs = get_year_2digit(dts)
%function yrs = get_year_2digit(dts)
% Get Year number for each datenum in DTS
  [yrs,ig,ig] = datevec(dts);
  yrs(yrs>=2100) = yrs(yrs>=2100)-2100;
  yrs(yrs>=2000) = yrs(yrs>=2000)-2000;
  yrs(yrs>1900) = yrs(yrs>=1900)-1900;
  yrs(yrs>1800) = yrs(yrs>=1800)-1800;
return;
