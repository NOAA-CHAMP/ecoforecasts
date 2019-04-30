function [sm,fh] = som_tses(tses_or_vecs,mapdims,npts,fdts,fldnms,ttl,showSort,doNorm,varargin)
%function [sm,fh] = som_tses(tses_or_vecs,mapdims,npts,fdts,fldnms,ttl,showSort,doNorm[,SOM_MAKE options...])
%
% Self-Organizing Map (Artificial Neural Network) training and analysis:
% First arg may be cell array of time series structs (TSES), or array of data
% vectors (VECS); if the former, FDTS or NPTS must specify starting date or
% length of each data vector in TSES{end}. A SOM of size MAPDIMS is trained.
% For display, each SOM node is reshaped into time series of NPTS each. If an
% array VECS is given but NPTS is not given, code plots just one time series.
%
% CALLS:
%   SOM_BMUS,SOM_HITS,SOM_MAKE (http://www.cis.hut.fi/projects/somtoolbox)
%   PLOTS,MULTIPLOT,SUPLABEL (MATLAB Central - File Exchange)
%   DATETICK3 (Ecoforecasts toolbox)
%
% Last Saved Time-stamp: <Sun 2013-11-03 18:03:30 Eastern Standard Time gramer>

  if ( ~exist('mapdims', 'var') || isempty(mapdims) )
    mapdims = [3 3];
  end;
  nmaps = mapdims(1) * mapdims(2);
  if ( ~exist('npts', 'var') || isempty(npts) )
    if ( ~isvector(tses_or_vecs) && isnumeric(tses_or_vecs) )
      % DEFAULT: each vector is composed of just one time series
      npts = size(tses_or_vecs,2);
    else
      error('If first arg is time series, arg NPTS is required');
    end;
  end;

  if ( ~exist('fldnms', 'var') || isempty(fldnms) )
    fldnms = {'',''};
  end;
  independent_fldnms = '';
  for fldix = 1:numel(fldnms)-1
    independent_fldnms = [independent_fldnms,',',fldnms{fldix}];
  end;

  if ( ~exist('ttl', 'var') || isempty(ttl) )
    ttl = '';
  end;

  if ( ~exist('showSort', 'var') || isempty(showSort) )
    showSort = true;
    %showSort = false;
  end;

  if ( ~exist('doNorm', 'var') || isempty(doNorm) )
    %doNorm = false;
    doNorm = true;
  end;

  if ( ~isempty(varargin) )
    args = varargin;
  else
    %initStr = 'randinit';
    initStr = 'lininit';
    args = {initStr, 'batch', 'tracking',0, 'shape','sheet', 'training','long', 'neigh','ep'};
  end;

  % Parse first argument
  if ( iscell(tses_or_vecs) || is_valid_ts(tses_or_vecs) )
    if ( ~iscell(tses_or_vecs) )
      tses = {tses_or_vecs};
    else
      tses = tses_or_vecs;
    end;
    nvars = numel(tses);
    if ( ~exist('fdts','var') || isempty(fdts) )
      % Divvy time series up into equal time-length vectors
      fdts = tses{end}.date(1:npts:end);
      fdts(end) = [];
    end;
    % Construct vectors from NPTS and starting dates FDTS
    for dtix=1:length(fdts)
      startdt = fdts(dtix,1);
      for vix=1:nvars
        [dterr,startix] = min(abs(tses{vix}.date-startdt));
        %if ( dterr > 2*median(diff(tses{vix}.date)); error('Gappy data!'); end;
        vecs(dtix,(npts*(vix-1))+1:(npts*vix)) = tses{vix}.data(startix:startix+(npts-1));
      end;
    end;
  elseif ( ~isvector(tses_or_vecs) && isnumeric(tses_or_vecs) )
    vecs = tses_or_vecs;
    ndts = size(vecs,1);
    nvars = (size(vecs,2)/npts);
    if ( nvars ~= floor(nvars) || nvars < 1 )
      error('size(VECS,2)/NPTS is not a positive integer!');
    end;
  else
    error('First arg must be a cell array of TSes or a vector!');
  end;
  clear tses_or_vecs;


  if ( doNorm )
    for vix=1:nvars
      x = vecs(:,(npts*(vix-1))+1:(npts*vix));
      normFac(1,(npts*(vix-1))+1:(npts*vix)) = repmat(nanstd(x(:)),[1,npts]);
      normAdd(1,(npts*(vix-1))+1:(npts*vix)) = repmat(nanmean(x(:)),[1,npts]);
    end;
    vecs = (vecs - repmat(normAdd,[size(vecs,1),1])) ./ repmat(normFac,[size(vecs,1),1]);
  end;

  disp('Training Self-Organizing Map...');

  sm = som_make(vecs, 'msize',mapdims, args{:});

  % Add Best-Matching Units vector as a member of each 'sm' struct.
  % NOTE: For SOM Unit (map-node, "mode") 'n', framedts((sm.bmus==n),[1 end])
  % are the start- and end-dates for all the frames best matched by that Unit.
  [sm.bmus, sm.qerrs] = som_bmus(sm, vecs);

  % Calculate number of matched sample frames for each "mode" (Unit).
  nummatched = som_hits(sm, vecs);
  pctmatched = 100 * (nummatched / size(vecs,1));

  % Sort BMUs by number of frames, from most-matched to least.
  [numsorted, mode_order] = sort(nummatched, 'descend');

  if ( showSort )
    % Plot "modes" from most frequently matched to least
    disp_order = mode_order;
  else
    % Plot map as trained
    disp_order = 1:numel(nummatched);
  end;

  nds = sm.codebook;
  if ( doNorm )
    nds = (sm.codebook .* repmat(normFac,[size(nds,1),1])) + repmat(normAdd,[size(nds,1),1]);
  end;

  fh = fmg;
  set(fh, 'DefaultAxesFontSize', [8]);
  cmins = min(nds); cmn = min(reshape(cmins, [npts nvars]));
  cmaxs = max(nds); cmx = max(reshape(cmaxs, [npts nvars]));
  cmin = min([-abs(cmn) ; -abs(cmx)]); cmax = max([abs(cmn) ; abs(cmx)]);
  cmin = cmn; cmax = cmx;
  %%%% ???  cmin = repmat(0, size(cmn)); cmax = cmx;
  for orderix = 1:length(disp_order)
    ix = disp_order(orderix);
    subplot(mapdims(2), mapdims(1), orderix);
    hold on;
    x = [1:npts]/24;
    cdbk = reshape(nds(ix,:), [npts nvars]);
    do_plots(x, cdbk, cmin, cmax);
    hold off;
    title(sprintf( 'Node %d (%d frames - %0.1f%%) #%d', ...
                   ix,nummatched(ix),pctmatched(ix),find(nummatched(ix)==numsorted,1) ));
  end;
  suplabel(strcat(strrep(lower(independent_fldnms),'_','\_'),','), 'x');
  suplabel(sprintf( 'Self-Organizing Map %s %dH frames (N=%d): %s%s', ...
                    ttl, npts, size(vecs,1), ...
                    strrep(lower(fldnms{end}),'_','\_') ), 't');

return;


%%
% PRIVATE function [hl, he] = do_plots(x, Y, cmin, cmax)
% CALLS: 'plots' from MATLAB Central Exchange
function [hl, ha] = do_plots(x, Y, cmin, cmax)
  % ylabs = strjust(strvcat('Hour','SeaT','AirT','WindU','WindV'),'center');
  % [hl, ha] = plots(x, Y, 'left', ylabs, 5); 
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
        if ( cmin(ix) < cmax(ix) )
          ylim(ha(ix), [cmin(ix) cmax(ix)]);
        else
          % Needed when one component (e.g., latent heat flux) is always negative
          ylim(ha(ix), [min([cmin(ix) cmax(ix)]),max([cmin(ix) cmax(ix)])]);
        end;
      end;
    else
      ylm = max(abs(ylim(ha(ix))));
      ylim(ha(ix), [-ylm ylm]);
    end;
  end;
return;
