function station_wtc(stn1, fld1, stn2, fld2)
%function station_wtc(stn1, fld1, stn2, fld2)
%
% Calculate and plot a Wavelet Coherence Estimate between two fields of two
% stations, i.e., call WTC (qv.) on time series STN1.(FLD1) and STN2.(FLD2).
%
% Last Saved Time-stamp: <Mon 2010-06-28 20:55:18 Eastern Daylight Time gramer>

  [ix1,ix2] = intersect_dates(stn1.(fld1).date, stn2.(fld2).date);

  % x = [ stn1.(fld1).date(ix1) , stn1.(fld1).data(ix1) ];
  % y = [ stn2.(fld2).date(ix2) , stn2.(fld2).data(ix2) ];

  dts1 = stn1.(fld1).date(ix1(1)):(1/24):stn1.(fld1).date(ix1(end));
  dat1 = interp1(stn1.(fld1).date,stn1.(fld1).data,dts1);
  dts2 = stn2.(fld2).date(ix2(1)):(1/24):stn2.(fld2).date(ix2(end));
  dat2 = interp1(stn2.(fld2).date,stn2.(fld2).data,dts2);

  maxix = min(length(dts1),(24*365*2));

  x = [ dts1(1:maxix) , dat1(1:maxix) ];
  y = [ dts2(1:maxix) , dat2(1:maxix) ];

  % wt(x, 'BlackandWhite', 'MakeFigure',true);
  wtc(x,y, 'MonteCarloCount',20, 'MaxScale',1024, 'BlackandWhite');

return;
