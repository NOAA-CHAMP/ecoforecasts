function [tiphdls,cursorMode] = adddatatips(plthdl,evtXes,useNearest)
%function [tiphdls,cursorMode] = adddatatips(plthdl,evtXes,useNearest)
%
% Add DATACURSORMODE data tips at X points EVTXES in plotted object with
% handle PLTHDL (e.g., the handle return value from PLOT). If USENEAREST
% (DEFAULT: True), X values in EVTEXES need not be exact matches (==) to
% those in the plot's XData property: in those cases, uses index returned by
% MIN(ABS(XDATA-EVTEXES(ix))), with a potential warning. If USENEAREST is
% instead a string beginning with the characters 'ind' (for index), then
% evtXes is interpreted as an array of indices in XDATA instead of values.
%
% Optionally returns an array of handles TIPHDLS for all datatips created,
% and the DATACURSOR manager object cursorMode.
%
% Assumes caller previously set tool tip appearances as desired (e.g., by
% calling SET_DATETICK_CURSOR). Based on code by "Hoki" originally found at:
% https://stackoverflow.com/questions/29882186/set-data-tips-programmatically
%
% CALLS: GETFIG to find the Figure handle associated with PLTHDL.
%
% Last Saved Time-stamp: <Wed 2017-05-24 13:15:54 Eastern Daylight Time gramer>


  useIndices = false;
  if ( ~exist('useNearest','var') )
    useNearest = true;
  elseif ( strcmpi(useNearest,'ind') )
    useNearest = true;
    useIndices = true;
  end;

  fh = getfig(plthdl);
  if ( isempty(fh) )
    error('No Figure handle found for plot handle!');
  end;

  % Retrieve the DATACURSOR manager
  cursorMode = datacursormode(fh);

  xdata = get(plthdl,'XData');
  ydata = get(plthdl,'YData');
  try,
    zdata = get(plthdl,'ZData');
  catch,
    zdata = [];
  end;

  dx = median(diff(xdata));

  % Add a datatip for each event
  for ix = 1:numel(evtXes)
    if ( useIndices )
      idx = evtXes(ix);
    else
      % Find the index of the corresponding X-value
      idx = find( xdata == evtXes(ix) ) ;
      if ( isempty(idx) )
        if ( ~useNearest )
          error('No X-value matches event # %d: Try USENEAREST=True...',ix);
        end;
        [err,idx] = min(abs(xdata-evtXes(ix)));
        if ( (idx == 1 || idx == numel(evtXes)) && err > (dx*2) )
          warning('AddDataTips:Outlier','Event # %d may be outside of X range!',ix);
        end;
      end;
    end;

    if ( isempty(zdata) )
      pos = [xdata(idx), ydata(idx), 0];
    else
      pos = [xdata(idx), ydata(idx), zdata(idx)];
    end;

    tiphdls(ix) = cursorMode.createDatatip(plthdl) ;
    % Move it into the right place
    set(tiphdls(ix), 'MarkerSize',5, 'MarkerFaceColor','none', ...
                     'MarkerEdgeColor','r', 'Marker','o', 'HitTest','off', ...
                     'Position',pos);
    %update(tiphdls(ix), pos);
  end;

return;
