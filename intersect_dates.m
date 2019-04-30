function [ix1,ix2] = intersect_dates(dts1, dts2, tol)
%function [ix1,ix2] = intersect_dates(dts1, dts2, tol)
%
% Return indices of all the elements in two vectors of timestamps (datenums)
% DTS1 and DTS2, that match: a "matching" timestamp is defined as one where
% the difference between the elements of the two series is less than some
% tolerance TOL. DEFAULT TOL is 30 mins and a hair: (30+(0.005/60))/(24*60).
%
% Last Saved Time-stamp: <Fri 2015-03-20 15:30:50 Eastern Daylight Time gramer>

  ix1 = [];
  ix2 = [];

  if ( isempty(dts1) || isempty(dts2) )
    warning('Ecoforecasts:IntersectDates:EmptyInput','Empty inputs');
    %%%%%%%%%%%%%%%
    % EARLY RETURN
    %%%%%%%%%%%%%%%
    return;
  end;

  if ( any(diff(dts1) <= 0) || any(diff(dts2) <= 0) )
    error('Both date series must be monotonically increasing!');
  end;
  if ( ~exist('tol','var') || isempty(tol) )
    % Default tolerance for timestamp matching is 29 minutes
    tol = (30.0+(0.005/60.0))/(24.0*60.0);
  end;
  if ( tol <= 0 )
    error('Tolerance argument should be greater than zero!');
  end;

  switched = false;
  if ( numel(dts2) > numel(dts1) )
    switched = true;
    tmp = dts1;
    dts1 = dts2;
    dts2 = tmp;
    tmp = []; clear tmp;
  end;    

  preix1 = find( ((dts2(1)-tol) <= dts1) & (dts1 <= (dts2(end)+tol)) );
  preix2 = find( ((dts1(1)-tol) <= dts2) & (dts2 <= (dts1(end)+tol)) );
  if ( isempty(preix1) || isempty(preix2) )
    return;
  end;

  % First get all the low-hanging fruit - exact timestamp matches
  [ig,exactix1,exactix2] = intersect(dts1(preix1), dts2(preix2));
  if ( ~isempty(exactix1) )
    addix1 = preix1(exactix1);
    addix2 = preix2(exactix2);
    ix1 = [ix1(:)' addix1(:)'];
    ix2 = [ix2(:)' addix2(:)'];
    preix1(exactix1) = [];
    preix2(exactix2) = [];
  end
  %DEBUG:  disp({mfilename,'exact',numel(ix1),numel(ix2)});

  % Next get all the points within TOL/2 of each other - again, easy
  % (Seems like duplicated effort, but makes final loop MUCH faster!)
  d1 = round((2*dts1(preix1))/tol);
  d2 = round((2*dts2(preix2))/tol);
  [ig,exactix1,exactix2] = intersect(d1, d2);
  if ( ~isempty(exactix1) )
    addix1 = preix1(exactix1);
    addix2 = preix2(exactix2);
    ix1 = [ix1(:)' addix1(:)'];
    ix2 = [ix2(:)' addix2(:)'];
    preix1(exactix1) = [];
    preix2(exactix2) = [];
  end
  %DEBUG:  disp({mfilename,'tol/2',numel(ix1),numel(ix2)});

  % Finally get all the match remainders one at a time
  % These are points within TOL of, but >TOL/2 from each other
  % (Unlovely but somewhat fast way to do this...)
  ini1 = 1;
  ini2 = 1;
  while ( ini1 <= length(preix1) && ini2 <= length(preix2) )
    if ( abs(dts1(preix1(ini1))-dts2(preix2(ini2))) <= tol )
      % If we hit a single, go for a grand slam!
      maxix = min(length(preix1)-ini1, length(preix2)-ini2);
      x1 = preix1(ini1:ini1+maxix);
      x2 = preix2(ini2:ini2+maxix);
      goodix = find(abs(dts1(x1) - dts2(x2)) <= tol);
      % Careful with non-contiguous matches
      lastgoodix = find(diff(goodix)>1,1);
      if ( ~isempty(lastgoodix) ); goodix = goodix(1:lastgoodix); end;

      addix1 = x1(goodix);
      addix2 = x2(goodix);
      ix1 = [ix1(:)' addix1(:)'];
      ix2 = [ix2(:)' addix2(:)'];

      ini1 = ini1 + length(goodix);
      ini2 = ini2 + length(goodix);

    else
      if ( dts1(preix1(ini1)) < dts2(preix2(ini2)) )
        ini1 = ini1 + 1;
      else
        ini2 = ini2 + 1;
      end;
    end;
  end;
  %DEBUG:  disp({mfilename,'tol',numel(ix1),numel(ix2)});

  % Make sure numel(ix1) ~= numel(ix2) after unique'ing
  [ig,newix] = unique(ix1);
  ix1 = ix1(newix);
  ix2 = ix2(newix);
  [ig,newix] = unique(ix2);
  ix1 = ix1(newix);
  ix2 = ix2(newix);

  % SORT result: Assume dates came in sorted; IX2 was just SORTed by UNIQUE.
  % Fixes a bug where non-periodic time series produced out-of-order dates:
  % Sorting nearby matches can force exact matches apart! But if TS1(i) and
  % TS2(j) are within TOL of each other,  then so are all intermediate pts.
  ix1 = sort(ix1);

  if ( switched )
    tmp = ix1;
    ix1 = ix2;
    ix2 = tmp;
    tmp = []; clear tmp;
  end;    

  %DEBUG:  disp({mfilename,'final',numel(ix1),numel(ix2)});

return;
