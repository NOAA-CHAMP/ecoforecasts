function fts = sma_ts(ts,filterIndices,smafunc)
%function fts = sma_ts(ts,filterIndices,smafunc)
%
% Apply Simple Moving Average low-pass filter to time series 'TS'. Optional
% function_handle arg SMAFUNC specifies alternate accumulator function to
% use (DEFAULT: @NANMEAN), e.g., @NANSUM for summing rather than averaging.
%
% Last Saved Time-stamp: <Fri 2011-02-18 13:43:44  Lew.Gramer>

  if ( ~exist('smafunc','var') || isempty(smafunc) )
    smafunc = @nanmean;
  end;

  if ( filterIndices == 0 )
    fts = ts;
  else
    nframes = floor(length(ts)/filterIndices);
    tempts = reshape(ts(1:(nframes*filterIndices)),[nframes filterIndices])';
    temp2ts = feval(smafunc,tempts,1);
    tempIndices = 1:filterIndices:length(ts);
    temp2ts(end+1:length(tempIndices)) = nan;
    warning('off','MATLAB:chckxy:IgnoreNaN');
    fts = interp1(tempIndices,temp2ts,1:length(ts),'spline');
    warning('on','MATLAB:chckxy:IgnoreNaN');
  end;

return;
