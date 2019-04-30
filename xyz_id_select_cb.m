function output_txt = xyz_id_select_cb(obj,event_obj)
% Display the position of the data cursor
% obj          Currently not used (empty)
% event_obj    Handle to event object
% output_txt   Data cursor text string (string or cell array of strings).

%get(get(event_obj, 'Target'));
pos = get(event_obj,'Position');

% If there is X/YData, also display indices
hdl = get(event_obj, 'Target');
xix=[]; yix=[];
if ( ~isempty(hdl) )
    xdata = get(hdl, 'XData');
    ydata = get(hdl, 'YData');
    xix = find(pos(1)==xdata);
    yix = find(pos(2)==ydata);
end;
if ( ~isempty(xix) && ~isempty(yix) )
  output_txt = {['X: ',sprintf('%f (#%d)', pos(1),xix)],...
                ['Y: ',sprintf('%f (#%d)', pos(2),yix)]};

else
  output_txt = {['X: ',sprintf('%f', pos(1))],...
                ['Y: ',sprintf('%f', pos(2))]};
end;

% If there is a Z-coordinate in the position, display it as well
if length(pos) > 2
    output_txt{end+1} = ['Z: ',sprintf('%f', pos(3))];
end
