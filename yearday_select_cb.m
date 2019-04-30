function output_txt = yearday_select_cb(obj,event_obj)
% Display the position of the data cursor
% obj          Currently not used (empty)
% event_obj    Handle to event object
% output_txt   Data cursor text string (string or cell array of strings).

pos = get(event_obj,'Position');
dvec = datevec(pos(1));
yearday = num2str( (pos(1) - datenum(dvec(1),1,1) + 1),'%.2f' );

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

% output_txt = {['X: ',datestr(pos(1),'dd-mmm'),' (',yearday,')'],...
%     ['Y: ',sprintf('%g', pos(2)),pointno]};
output_txt = {['X: ',datestr(pos(1),'dd-mmm')],...
    ['Y: ',sprintf('%g', pos(2))]};

% If there is a Z-coordinate in the position, display it as well
if length(pos) > 2
    output_txt{end+1} = ['Z: ',sprintf('%g', pos(3))];
end
