function [Bs,Stats] = dump_robust_fit(stn,fitsfld,fitqflds)
%function [Bs,Stats] = dump_robust_fit(stn,fitsfld,fitqflds)
%
% Dump to the Command Window (DISP(SPRINTF(...))), statistics (R2, SI, A, B,
% RMSE) from a robust linear least-squares fit between FITSFLD and FITQFLDS.
% Scatter Index (SI) in this case is the RMSE / Std. Dev. of STN.(FITSFLD); B
% is "slope error" based on regression slope beta, i.e., 100% x abs(1-beta).
% FITQFLDS can be a CHAR or a CELLSTR (to report on multiple variables).
%
% Last Saved Time-stamp: <Thu 2013-05-16 15:04:26 Eastern Daylight Time Lew.Gramer>

  Bs = {};
  Stats = {};

  if ( ~exist('fitqflds','var') || isempty(fitqflds) )
    disp(sprintf('%-50s %-4s %-6s %-8s %-8s %-7s',['VS ',fitsfld],'R2','SI%','bias','slope%','RMSE'));
  else
    if ( ~iscell(fitqflds) )
      fitqflds = {fitqflds};
    end;
    for ix = 1:numel(fitqflds)
      tses(ix) = stn.(fitqflds{ix});
    end;
    tses(end+1) = stn.(fitsfld);

    % NOTE: This expects a struct array but returns a cell array!
    tses = intersect_tses([],tses);

    disp(sprintf('%-60s %-4s %-6s %-8s %-8s %-7s %-7s %-7s %-7s',...
                 ['VS ',fitsfld],'R2','SI%','bias','slope%','RMSE','RawRMSE','N [d]'));
    for ix = 1:numel(fitqflds)
      fitqfld = fitqflds{ix};
      [B,Stat] = scatter_fit_ts(tses{end},tses{ix},[],[],[],[],'none');
      Stat.raw_RMSE = sqrt( sum((tses{end}.data-tses{ix}.data).^2) ./ numel(tses{end}.data) );
      Bs{ix} = B;
      Stats{ix} = Stat;
      if ( Stat.p(2) > 0.01 )
        pstatStr = '*';
      else
        pstatStr = ' ';
      end;
      disp(sprintf('%-60s %4.2f %5.1f%% %+6.3f %7.3f%% %7.3f %7.3f %6.0f %s',...
                   fitqfld,Stat.R2_2,abs(100*Stat.s/nanstd(stn.(fitsfld).data)),...
                   B(1),100*abs(1-B(2)),Stat.s,Stat.raw_RMSE,...
                   Stat.N,pstatStr));
      B=[]; Stat=[]; clear B Stat
    end;
    tses = []; clear tses;
  end;

  if ( nargout < 1 )
    Bs=[]; Stats=[]; clear Bs Stats;
  end;

return;
