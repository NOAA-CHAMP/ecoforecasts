function [hl, ha] = do_tvar_peak_plots(x, Y, cmin, cmax)
  [hl, ha] = plots(x, Y, 'left', 8); 
  % stylel = repmat(strvcat('--','-.','-',':'), [2 1]);
  % stylel = repmat(strvcat('--','-'), [4 1]);
  stylel = repmat(strvcat('-'), [8 1]);
  colorl = get(gcf, 'DefaultAxesColorOrder');
  set(ha(end), 'YAxisLocation', 'right');
  for ix = 1:numel(hl)
    set(hl(ix), 'LineStyle', stylel(ix,:), 'Color', colorl(ix,:));
    set(ha(ix), 'XGrid', 'on');
    set(ha(ix), 'YColor', colorl(ix,:))
    if ( exist('cmin', 'var') && exist('cmax', 'var') )
      % NOTE: If cmin and cmax are present but empty, NO YLim is set...
      if ( ~isempty(cmin) && ~isempty(cmax) )
        ylim(ha(ix), [cmin(ix) cmax(ix)]);
      end;
    else
      ylm = max(abs(ylim(ha(ix))));
      ylim(ha(ix), [-ylm ylm]);
    end;
  end;
return;
