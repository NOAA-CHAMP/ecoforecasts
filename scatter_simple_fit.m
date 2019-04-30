function [B,Stats] = scatter_simple_fit(x,y,varargin)
%function [B,Stats] = scatter_simple_fit(x,y[,args])
%
% Like SCATTER_FIT (v.), but this version simply returns the result of a
% simple (non-robust) linear regression with REGRESS.
%
% Optional ARGS are passed through to REGRESS.
%
% Last Saved Time-stamp: <Thu 2018-04-26 12:03:22 Eastern Daylight Time gramer>

  X(1:length(x),1) = 1;
  X(1:length(x),2) = x(:);
  [B,BINT,R,RINT,Stats.regress_stats] = regress(y,X,varargin{:});
  Stats.N = min(length(find(isfinite(x))),length(find(isfinite(y))));

  Stats.yhat = B(1) + (B(2).*x);
  Stats.R2_1 = corr(y,Stats.yhat).^2;
  Stats.sse = Stats.dfe .* (Stats.robust_s.^2);
  Stats.ssr = norm(Stats.yhat - nanmean(Stats.yhat)).^2;
  Stats.R2_2 = 1 - (Stats.sse ./ (Stats.sse + Stats.ssr));
  % Rank correlation
  [Stats.R2_3,Stats.p_3] = corr(y,Stats.yhat,'type','Spearman');
  Stats.R2_3 = Stats.R2_3.^2;

return;
