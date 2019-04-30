function gt_ts(var1, var2)
%function gt_ts(var1, var2)
% 
% "Ground-truth" (statistically compare) two time series against one another.
% 
% Last Saved Time-stamp: <Wed 2009-02-04 15:50:07 Eastern Standard Time gramer>


  if ( ~isfield(var1, 'date')  || ~isfield(var1, 'data') )
    error('Struct "var1" must have both "date" and "data" fields!');
  end;
  if ( ~isfield(var2, 'date')  || ~isfield(var2, 'data') )
    error('Struct "var2" must have both "date" and "data" fields!');
  end;

  half_hour = (0.5 / 24.0);

  % No cutoff day for now...
  cutoff_day = datenum(1987, 01, 01);

  newix = 1;
  for ix = 1:length(var1.date)
    matchix = find(abs(var1.date(ix) - var2.date) < half_hour);
    if ( ~isempty(matchix) && var1.date(ix) >= cutoff_day )
      imdts(newix) = var1.date(ix);
      imdat(newix) = var1.data(ix);
      smdts(newix) = var1.date(ix);
      smdat(newix) = mean(var2.data(matchix));
      newix = newix + 1;
    end;
  end;

  ndays = length(unique(fix(imdts)));

  X = [ ones(size(imdat)) ; imdat ]';
  Y = smdat';
  [B, BINT, R, RINT, STATS] = regress(Y, X);
  RMSE = sqrt(sum(R .^ 2));

  figure;
  hold on;
  plot(imdat, smdat, 'b.');
  plot(imdat, (B(1) + (B(2) .* X(:,2))), 'r:');
  ylim([0 max(max(imdat), max(smdat))]);
  xlabel('Var1'); ylabel('Var2');
  legend('Raw data', ...
         sprintf('var2 = %g + %g*var1: R^2=%g', B(1), B(2), STATS(1)), ...
         'Location', 'Best');
  title(sprintf('Comparison of two time series var1 and var2 means: %d days', ndays));
  %print('-dpng', [station '_groundtruth_par.png']);


return;
