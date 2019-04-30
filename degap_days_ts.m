function newts = degap_days_ts(ts)
%function newts = degap_days_ts(ts)
%
% Return new TS with every whole day missing more than 6 h of data removed
% Created: 2011-05-20
%
% Last Saved Time-stamp: <Mon 2015-04-13 12:59:40 Eastern Daylight Time gramer>

  newts = ts;

  maxgap = 3;
  missinghrs = 6;

  dys = floor(newts.date);
  udys = unique(dys);

  %DEBUG:
  tic, disp('GRPSTATS'); 

  n = grpstats(newts.date,dys,'numel');

  %DEBUG:
  toc,


  %DEBUG:
  tic, disp('loop'); 

  % Remove every whole day missing more than 6 hours of data
  badix = find(n<(24-missinghrs));
  %DEBUG:
  disp(['Removing ' num2str(numel(badix)) ' bad days']);
  newts.data(ismember(dys,udys(badix))) = [];
  newts.date(ismember(dys,udys(badix))) = [];

  %DEBUG:
  toc,

return;
