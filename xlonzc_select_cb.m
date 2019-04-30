function output_txt = xlonzc_select_cb(obj,event_obj)
% Display the position of the data cursor
% obj          Currently not used (empty)
% event_obj    Handle to event object
% output_txt   Data cursor text string (string or cell array of strings).

pos = get(event_obj,'Position');

lon = pos(2);
if (lon > 180); lon = lon - 360; end;

output_txt = {['X: ',num2str(pos(1))],...
    ['Y: ',num2str(lon)]};

% If there is a Z-coordinate in the position, display it as well
if length(pos) > 2
    output_txt{end+1} = ['Z: ',sprintf('%f', pos(3))];
end

% If there is CDate (e.g., this graph was produced by scatter or scatter3,
% instead of plot or plot3), then find the C data and also display it/them.
hdl = get(event_obj, 'Target');
if ( ~isempty(hdl) )
  try
    xdata = get(hdl, 'XData');
    ydata = get(hdl, 'YData');
    zdata = get(hdl, 'ZData');
    try,
      cdata = get(hdl, 'CData');
    catch,
      cdata = '';
    end;
    % If plot caller did not mesh X and Y, do so now
    if (any(size(xdata) == 1))
      [xdata,ydata] = meshgrid(xdata,ydata);
    end;
    idx = find(xdata == pos(1) & ydata == pos(2) & zdata == pos(3));
    if ( isempty(idx) || isempty(cdata) )
      output_txt{end+1} = ['C: NO MATCH?'];
    else
      % Note comma ",": may be MORE THAN ONE matching C value!
      output_txt{end+1} = ['C: ',sprintf('%f,', cdata(idx))];
    end
  catch
  end;
end
