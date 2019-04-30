function f = findiff_ts(t,n)
%function f = findiff_ts(t,n)
%
% Calculate the N-point finite difference (DEFAULT: 3) of the time series T
% (a struct with fields T.date and T.data), and return it in F.date, F.data.
% The resulting TS has the same gaps as the time series T, except that the
% CEIL(N/2) points before and after each gap are also removed.
%
% Last Saved Time-stamp: <Fri 2015-06-26 23:21:12 Eastern Daylight Time gramer>

  % Size of (centered) finite difference template
  if (~exist('n','var') || isempty(n) )
    n = 3;
  end;

  switch ( round(n) )
   case 1,
    f.data = diff(t.data);
    endoff = 1;
   case 3,
    f.data = (-t.data(1:(end-2)) + t.data(3:end)) / 2;
    endoff = 1;
   case 5,
    f.data = (t.data(1:(end-4)) - 8*t.data(2:(end-3)) + ...
              8*t.data(4:(end-1)) - t.data(5:end)) / 12;
    endoff = 2;
   case 7,
    f.data = (-t.data(1:(end-6)) + 9*t.data(2:(end-5)) - ...
              45*t.data(3:(end-4)) + 45*t.data(5:(end-2)) - ...
              9*t.data(6:(end-1)) + t.data(7:end)) / 60;
    endoff = 3;
   case 9,
    f.data = (1/280)*t.data(1:end-8) + (-4/105)*t.data(2:end-7) + ...
             (1/5)*t.data(3:end-6) + (-4/5)*t.data(4:end-5) + ...
             (4/5)*t.data(6:end-3) + (-1/5)*t.data(7:end-2) + ...
             (4/105)*t.data(8:end-1) + (-1/280)*t.data(9:end);
    endoff = 4;
   otherwise,
    error('Do not know how to finite difference with %d-sample template!', n);
  end;

  % 1: 1        3: 2            5: 3            7: 4
  begix = ceil(n/2);
  % 1: end-1    3: end-1        5: end-2        7: end-3
  f.date = t.date(begix:end-endoff);

  % Remove any gaps, and any finite-difference artifacts from gaps
  fN = 2 * min(diff(f.date));
  gapix = find(diff(f.date) > fN);
  for ix = 1:begix
    gapix = unique([gapix-ix gapix(:) gapix+ix]);
  end;
  gapix(gapix < 1) = [];
  gapix(gapix > length(f.date)) = [];

  f.date(gapix) = [];
  f.data(gapix) = [];

return;
