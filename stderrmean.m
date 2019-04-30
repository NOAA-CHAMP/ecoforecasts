function StandardError = stderrmean(X)
%function StandardError = stderrmean(X)
%
% Calculate standard error of the mean of X

  X = X(isfinite(X));
  StandardError = std(X(:)) / sqrt(length(X(:)));

return;
