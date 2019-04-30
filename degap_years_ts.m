function newts = degap_years_ts(ts,missingdys,maxgap)
%function newts = degap_years_ts(ts,missingdys,maxgap)
%
% Return NEWTS by removing every year missing more than MISSINGDYS days
% (DEFAULT: 120), or having a contiguous gap > MAXGAP (DEFAULT: 30 d).
%
% Last Saved Time-stamp: <Mon 2015-04-13 13:23:29 Eastern Daylight Time gramer>

  newts = ts;

  if ( ~exist('missingdys','var') || isempty(missingdys) )
    missingdys = 120;
  end;
  if ( ~exist('maxgap','var') || isempty(maxgap) )
    maxgap = 30;
  end;

  yrs = get_year(newts.date);
  uyrs = unique(yrs);

  %DEBUG:
  tic, disp('GRPSTATS'); 

  n = grpstats(newts.date,yrs,'numel');

  %DEBUG:
  toc,

  %DEBUG:
  tic, disp('filter'); 

  % Remove every whole year missing more than MISSINGDYS of data
  badix = find(n<(365-missingdys));
  %DEBUG:
  disp(['Removing ' num2str(numel(badix)) ' bad years']);
  newts.data(ismember(yrs,uyrs(badix))) = [];
  newts.date(ismember(yrs,uyrs(badix))) = [];
  uyrs(badix) = [];

  %DEBUG:
  toc,


  %DEBUG:
  tic, disp('loop'); 

  % Remove every year with a gap longer than MAXGAP days
  badyrs = [];
  for yr = uyrs(:)';
    yrix = find(get_year(newts.date)==yr);
    dd = diff(newts.date(yrix));
    if ( max(dd(:)) > maxgap )
      badyrs(end+1) = yr;
    end;
  end;
  newts.data(ismember(get_year(newts.date),badyrs)) = [];
  newts.date(ismember(get_year(newts.date),badyrs)) = [];

  %DEBUG:
  toc,

return;
