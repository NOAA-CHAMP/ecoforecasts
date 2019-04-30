1;

if ( ~exist('SAVE_LON','var') || ~exist('SAVE_LAT','var') || ~exist('SAVE_DAT','var') )
  error('First run SAVE_LON = LON; SAVE_LAT = LAT; SAVE_DAT = DAT;');
end;

LON=[]; LAT=[]; DAT=[]; clear LON LAT DAT

LON = SAVE_LON; LAT = SAVE_LAT; DAT = SAVE_DAT;

ix = find(-81.2 <= LON & LON <= -80.8  &  24.5 <= LAT & LAT <= 24.8);

[rix,cix] = ind2sub(size(LON),ix);

rixen = min(rix):max(rix);
cixen = min(cix):max(cix);

LON = LON(rixen,cixen);
LAT = LAT(rixen,cixen);
DAT = DAT(rixen,cixen);

clear ans cix cixen ix rix rixen

process_F010
