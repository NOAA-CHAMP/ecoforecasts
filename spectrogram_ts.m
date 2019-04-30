function [dts,per,pwr,spc] = spectrogram_ts(ts,gapmethod,maxgap,method,doPlot)
%function [dts,per,pwr,spc] = spectrogram_ts(ts[,gapmethod[,maxgap[,method[,doPlot]]]])
%
% Estimate evolutionary (time-dependent) Power Spectral Density (PSD) of time
% series struct TS (with equal-length vector fields .date and .data), using a
% Short-Time Fourier Transform (STFT). Interpolate whole time series by using
% INTERP1 with interpolation METHOD (Default: 'pchip'). However, if GAPMETHOD
% is given and 'piecewise', each segment of TS that has no gaps longer than
% MAXGAP (Default: 3/24) is processed separately. Calls SPECTROGRAM (see) on
% each contiguous sub-segment of interpolated time series. If DOPLOT, then do
% plot of the time-dependent PSD, and of a time series of super-diurnal power.
%
% DTS a vector of DATENUMs for the center of each time-period analyzed, PER a
% vector of 2*pi/frequencies analyzed, PWR the PSD, and SPC the spectrogram.
% If no output args are specified and DOPLOT is NOT FALSE, then do all plots.
%
% Last Saved Time-stamp: <Fri 2012-03-09 12:35:49  Lew.Gramer>

  if ( ~exist('gapmethod','var') || isempty(gapmethod) )
    gapmethod = 'continuous';
  end;
  if ( ~exist('maxgap','var') || isempty(maxgap) )
    maxgap = (1/24);
  end;
  if ( ~exist('method','var') || isempty(method) )
    method = 'pchip';
  end;
  if ( ~exist('doPlot','var') || isempty(doPlot) )
    if ( nargout < 1 )
      doPlot = true;
    else
      doPlot = false;
    end;
  end;

  xlbl = inputname(1);
  if ( isempty(xlbl) )
    xlbl = 'Time Series';
  end;

  %%%% ???  How do we get the WHOLE RECORD?

  switch ( lower(gapmethod) )
   case 'piecewise',
    gapix = find(diff(ts.date) > maxgap);
    begixes = [1 (gapix+1)];
    endixes = [gapix length(ts.date)];
   otherwise,
    begixes = [1];
    endixes = [length(ts.date)];
  end;

  for ix = 1:length(begixes)
    begix = begixes(ix);
    endix = endixes(ix);

    dt = ts.date(begix):(1/24):ts.date(endix);
    da = interp1(ts.date,ts.data,dt,method);

    win=1024;
    noverlap=512;
    nfft=1024;
    [s,f,t,p] = spectrogram(da,win,noverlap,nfft);

    dts = dt(1)+(t/24);
    per = 2*pi./(f*24);
    pwr = p;
    spc = s;
  end;

  if ( doPlot )
    figure; maxigraph;
    surf(dts,per,10*log10(abs(pwr)+eps),'EdgeColor','none');
    set(gca,'yscale','log');
    colormap(jet); view(0,90); datetick3;
    axis([min(dts) max(dts) min(per) max(per)]);
    titlename([strrep(xlbl,'_','\_') num2str([win,noverlap,nfft])]);

    [ig,ix24]=min(abs(per-(24/24)));
    [ig,ix12]=min(abs(per-(12/24)));
    [ig,ix8]=min(abs(per-(8/24)));
    figure;
    plot(dts,sum(pwr([ix12-1:ix8],:))./pwr([ix24],:));
    datetick3;
  end;

return;
