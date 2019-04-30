1;

[n,t,r] = xlsread('C:\Users\lew.gramer\Documents\Postdoc\gcoos_2016-01-28.csv');

goodix = find(cellfun(@length,r(:,5))==20);
stnm = num2str([r{goodix,2}]');
lats = [r(goodix,3)]';
lons = [r(goodix,4)]';
tdts = datenum(r(goodix,5),'yyyy-mm-ddTHH:MM:SS');
deps = [r(goodix,6)]';
dirs = [r{goodix,7}]';
spds = [r{goodix,8}]';
vspd = [r{goodix,9}]';
