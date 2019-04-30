function [lhs,uyrs,fh,ax] = year_overlap_plot_ts(ts)
%function [lhs,uyrs,fh,ax] = year_overlap_plot_ts(ts)
%
% PLOT3 (v.) time series struct TS by year-day (GET_YEARDAY) and year, so
% that all year-days overlap but years are color coded.  Sets datamode cursor
% to provide useful tooltips.  Returns vector of linehandles, vector of
% unique years plotted, and figure and axes handles.
%
% Last Saved Time-stamp: <Wed 2012-03-28 10:35:55  Lew.Gramer>

  fh = gcf;
  ax = gca;

  yrs = get_year(ts.date);
  uyrs = unique(yrs);

  delt = min(diff(ts.date));
  % dts = repmat(nan,[numel(uyrs),numel([1:delt:367+(delt/2)])]);
  % dat = repmat(nan,size(dts));

  % We need HOLD ON - but respect caller's environment
  holdState = ishold(ax);
  hold on;

  co = get(ax,'ColorOrder');
  ms = '.ox+*sdv^<>ph';
  ncos = size(co,1);

  for yrix=1:numel(uyrs)
    yr = uyrs(yrix);
    % dt = [datenum(yr,1,1):delt:datenum(yr,12,31,23,59,0)]';
    % [dtix,tsix] = intersect_dates(dt,ts.date);
    % dts(yrix,dtix) = get_yearday(ts.date(tsix));
    % dat(yrix,dtix) = ts.data(tsix);
    ix = find(get_year(ts.date)==yr);
    plot3( get_yearday(ts.date(ix)),ts.data(ix),yrs(ix),...
           'Color',co(mod(yrix-1,ncos)+1,:), 'Marker',ms(ceil(yrix/ncos)) );
  end;
  % plot(dts,dat,'.');
  % datetick3;

  % Set informative data cursor (better than TYZ_SELECT_CB)
  datetick;
  h = datacursormode(gcf);
  set(h, 'UpdateFcn', @day_data_year_select_cb);

  % User may choose to re- or undo legend themselves after this call
  legend(num2str(uyrs), 'Location','Best');

  % Restore caller's HOLD state
  if ( ~holdState )
    hold off;
  end;

return;



%%%% INTERNAL FUNCTION day_data_year_select_cb

function output_txt = day_data_year_select_cb(obj,event_obj)
% Display the position of the data cursor
% obj          Currently not used (empty)
% event_obj    Handle to event object
% output_txt   Data cursor text string (string or cell array of strings).

  pos = get(event_obj,'Position');
  dvec = datevec(pos(1));
  jday = num2str( fix(pos(1) - datenum(dvec(1), 1, 1) + 1) );

  pointno = '';
  % If we can find A UNIQUE POINT corresponding to our selected data, then also
  % show user the index in our XData vector corresponding to that point
  eventattrs = get(event_obj);
  if ( isfield(eventattrs,'Target') )
    tgt = eventattrs.Target;
    attrs=get(tgt);
    if ( all(isfield(attrs,{'XData','YData'})) ...
         && ~isempty(attrs.XData) && ~isempty(attrs.YData) ...
         && all(size(attrs.XData) == size(attrs.YData)) )
      pointix = find(pos(1) == attrs.XData & pos(2) == attrs.YData);
      if ( numel(pointix) == 1 )
        pointno = [' (X#' num2str(pointix) ')'];
      end;
    end;
  end;

  % If there is a YEAR (Z-coordinate) in the position, display it
  if length(pos) > 2
    dt = datenum(pos(3),1,1) + get_yearday(pos(1));
    output_txt = {['X: ',datestr(dt),' (',jday,')'],...
                  ['Y: ',sprintf('%g', pos(2)),pointno]};
  else
    output_txt = {['X: ',datestr(pos(1),'dd-mmm HH:MM:SS'),' (',jday,')'],...
                  ['Y: ',sprintf('%g', pos(2)),pointno]};
  end;

return;
